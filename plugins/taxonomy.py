#! /usr/bin/env python3

'''
The MIT License

Copyright (c) 2017 by Anthony Westbrook, University of New Hampshire <anthony.westbrook@unh.edu>

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in 
the Software without restriction, including without limitation the rights to 
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

# Perform taxonomic categorization on a PALADIN UniProt report

import os
import re
import operator
import argparse
import shlex
import plugins.core

TAXONOMY_DOWNLOAD = 'http://www.uniprot.org/taxonomy/?query=&sort=score&format=tab'
TAXONOMY_RAW = 'lineage-raw.dat'
TAXONOMY_LINEAGE = 'lineage.dat'

moduleDefinition = None
lineageLookup = dict()

# Plugin connection definition
def pluginConnect(passDefinition):
    passDefinition.name = "taxonomy"
    passDefinition.description = "Perform taxonomic grouping and abundance reporting"
    passDefinition.versionMajor = 1
    passDefinition.versionMinor = 1
    passDefinition.versionRevision = 2

    passDefinition.callbackInit = taxonomyInit
    passDefinition.callbackMain = taxonomyMain

    plugins.taxonomy.moduleDefinition = passDefinition

# Ensure all UniProt taxonomy data is cached
def taxonomyInit():
   # Download UniProt taxonomic lineage data
    downloadLineage()
    parseLineage()
    cacheLineage()

# Plugin main
def taxonomyMain(passArguments):
    argParser = argparse.ArgumentParser(description='PALADIN Pipeline Plugins: Taxonomy', prog='taxonomy')
    argParser.add_argument('-i', dest='input', type=str, required=True, help='PALADIN UniProt report')
    argParser.add_argument('-q', dest='quality', type=int, required=True, help='Minimum mapping quality filter')
    argParser.add_argument('-c', dest='custom', type=str, default='', help='Species-parsing regex pattern for non-Uniprot entries')
    argParser.add_argument('-t', dest='type', type=str, default='species', choices=['children', 'species'], help='Type of report')
    argParser.add_argument('-f', dest='filter', type=str, choices=['unknown', 'custom', 'group'], nargs='*', help='Filter non-standard entries')
    argParser.add_argument('-s', dest='sam', type=str, help='SAM file for reporting reads contributing to each taxonomic rank')
    levelGroup = argParser.add_mutually_exclusive_group(required=True)
    levelGroup.add_argument('-l', dest='level', type=int, help='Hierarchy level')
    levelGroup.add_argument('-r', dest='rank', type=str, help='Regex pattern for named rank')

    arguments = argParser.parse_known_args(shlex.split(passArguments))
    if arguments[0].filter == None: arguments[0].filter = []

    # Get entries
    fullEntries = plugins.core.PaladinEntry.getEntries(arguments[0].input, arguments[0].quality, arguments[0].custom)
    speciesLookup = getSpeciesLookup(fullEntries)
    taxaEntries, taxaCount = aggregateTaxa(fullEntries, arguments[0].filter)
    lineageTree = treeifyLineage(taxaEntries, taxaCount)

    # Check report combinations
    if arguments[0].type == 'children':
        if arguments[0].level != None:
            # Children of a flattened level
            renderEntries, renderCount = flattenTree(lineageTree, arguments[0].level)
            header = "Count\tAbundance\tRank {0}".format(arguments[0].level)
        else:
            # Children of a rank subtree
            rankTree = findRankSubtree(lineageTree, arguments[0].rank)
            renderEntries, renderCount = getTreeChildren(rankTree)
            header = "Count\tAbundance\tRank"
    else:
        if arguments[0].level != None:
            # All species for entire tree
            renderEntries, renderCount = filterTaxa(taxaEntries, ".*")
            renderEntries = {speciesLookup[taxon]: renderEntries[taxon] for taxon in renderEntries}
            header = "Count\tAbundance\tSpecies"
        else:
            # All species of a rank subtree
            renderEntries, renderCount = filterTaxa(taxaEntries, arguments[0].rank)
            renderEntries = {speciesLookup[taxon]: renderEntries[taxon] for taxon in renderEntries}
            header = "Count\tAbundance\tSpecies"    

    if not arguments[0].sam:
        # Render report (if non-SAM mode)
        renderAbundance(renderEntries, header, renderCount)
    else:
        # In SAM mode, generate reads report if requested
        samData = groupSam(renderEntries, arguments[0].sam, speciesLookup)
        renderSam(samData)

# Download lineage information from UniProt            
def downloadLineage():
    downloadPath = os.path.join(plugins.core.getCacheDir(plugins.taxonomy.moduleDefinition), TAXONOMY_RAW)
    plugins.core.downloadURL(TAXONOMY_DOWNLOAD, downloadPath, 'UniProt taxonomic lineage data')

# Parse and filter relevant lineage information from raw data
def parseLineage():
    inputPath = os.path.join(plugins.core.getCacheDir(plugins.taxonomy.moduleDefinition), TAXONOMY_RAW)
    outputPath = os.path.join(plugins.core.getCacheDir(plugins.taxonomy.moduleDefinition), TAXONOMY_LINEAGE)

    # Check if file has already been parsed
    if os.path.exists(outputPath): return

    with open(inputPath, 'r') as inputHandle:
        with open(outputPath, 'w') as outputHandle:
            # Parse out taxonomy for entries with Mnemonics
            inputHandle.readline()

            for line in inputHandle:
                line = line.rstrip()
                fields = line.split("\t")
                
                if len(fields) < 9: continue

                if fields[1]: 
                    outputHandle.write("{0}\t{1}\n".format(fields[1], fields[8]))

# Cache lineage data for repeated taxonomic queries
def cacheLineage():
    # Clear cache
    lineageLookup.clear()

    # Open processed lineage data
    lineagePath = os.path.join(plugins.core.getCacheDir(plugins.taxonomy.moduleDefinition), TAXONOMY_LINEAGE)
    with open(lineagePath, 'r') as fileHandle:

        # Cache lineage data
        plugins.core.sendOutput("Caching taxonomic lineage data...", 'stderr')

        for line in fileHandle:
            line = line.rstrip()
            fields = line.split("\t")
          
            if len(fields) < 2: 
                lineageLookup[fields[0]] = 'Unknown'
            else:
                lineageLookup[fields[0]] = fields[1]

# Group, filter, and count taxa
def aggregateTaxa(passEntries, passFilters):
    PaladinEntry = plugins.core.PaladinEntry
    retTaxa = dict()
    retTotal = 0

    # For each entry, filter, get species ID, and summate count
    for entry in passEntries:
        # Check for removed filters
        if "unknown" in passFilters and passEntries[entry].type == PaladinEntry.TYPE_UNKNOWN: continue
        if "custom" in passFilters and passEntries[entry].type == PaladinEntry.TYPE_CUSTOM: continue
        if "group" in passFilters and passEntries[entry].type == PaladinEntry.TYPE_UNIPROT_GROUP: continue

        key = (passEntries[entry].speciesID, passEntries[entry].speciesFull)
        if not key in retTaxa: retTaxa[key] = 0

        count = passEntries[entry].count
        retTaxa[key] += count
        retTotal += count
   
    return retTaxa, retTotal 

# Filter taxa having a specific taxonomic rank (as regex)
def filterTaxa(passTaxa, passPattern):
    retTaxa = dict()
    retCount = 0

    for taxon in passTaxa:
        if taxon[0] in lineageLookup: lineage = lineageLookup[taxon[0]]
        else: lineage = 'Unknown'

        if re.search(passPattern, lineage):
            retTaxa[taxon] = passTaxa[taxon]
            retCount += passTaxa[taxon]

    return retTaxa, retCount

# Create a speciesID to species full lookup
def getSpeciesLookup(passEntries):
    retLookup = dict()

    for entry in passEntries:
        # For virtual group entries, save as full species name
        key = (passEntries[entry].speciesID, passEntries[entry].speciesFull)

        if not key in retLookup:
            retLookup[key] = passEntries[entry].speciesFull

    return retLookup
    
# Create tree structure out of taxa/lineage data
def treeifyLineage(passTaxaEntries, passTaxaCount):
    retTree = (dict(), passTaxaCount)

    for taxon in passTaxaEntries:
        if not taxon[0] in lineageLookup: continue
        rawLineage = lineageLookup[taxon[0]]
        ranks = rawLineage.split(';')
        parent = retTree[0]

        for rank in ranks:
            rank = rank.strip()

            # Start new dictionary, initialize size cache
            if not rank in parent: parent[rank] = (dict(), 0) 

            parent[rank] = (parent[rank][0], parent[rank][1] + passTaxaEntries[taxon])
            parent = parent[rank][0]

    return retTree

# Find a specific rank instance subtree within a lineage tree
def findRankSubtree(passLineage, passRank):
    # Breadth-wise search
    for rank in passLineage[0]:
        if rank == passRank: return passLineage[0][rank]

    # Recurse if necessary
    for rank in passLineage[0]:
        retTree = findRankSubtree(passLineage[0][rank], passRank)
        if retTree: return retTree

# Get immediate children for a specific lineage rank branch
def getTreeChildren(passLineage):
    retChildren = dict()
    retCount = 0

    for rank in passLineage[0]:
        count = passLineage[0][rank][1]
        retChildren[rank] = count
        retCount += count

    return retChildren, retCount

# Flatten a lineage tree to a specific level
def flattenTree(passLineage, passLevel):
    retFlatEntries = dict()
    retFlatCount = 0

    # Negative level indicates last level
    if passLevel < 0:
        for entry in passLineage[0]:
            if not passLineage[0][entry][0]:
                retFlatEntries[entry] = passLineage[0][entry][1]
                retFlatCount += passLineage[0][entry][1]

    if passLevel == 0:
        # At level 0, extract lineage
        retFlatEntries, retFlatCount = getTreeChildren(passLineage)
    else:
        # At higher levels, recurse and append
        for entry in passLineage[0]:
            subtreeEntries, subtreeCount = flattenTree(passLineage[0][entry], passLevel - 1)
            retFlatEntries.update(subtreeEntries)
            retFlatCount += subtreeCount

    return retFlatEntries, retFlatCount
    
# Render standard dictionary abundance report
def renderAbundance(passData, passHeader, passTotal):
    # Sort entries
    sortedData = sorted(passData.items(), key = operator.itemgetter(1), reverse=True)

    # Render data
    plugins.core.sendOutput(passHeader)

    for row in sortedData:
        abundance = row[1] / passTotal * 100
        outputLine = "{0}\t{1}\t{2}".format(row[1], abundance, row[0])
        plugins.core.sendOutput(outputLine)

# Group SAM reads per taxonomic rank
def groupSam(passData, passSam, passSpecies):
    retData = dict()

    # Simplify species lookup keys
    speciesLookup = {item[0]: item[1] for item in passSpecies}

    # Read SAM entries
    samEntries = plugins.core.SamEntry.getEntries(passSam, 0)

    for entry in samEntries:
        reference = samEntries[entry].reference
        
        # Currently only handles UniProt
        if not '_' in reference: continue

        mnemonic = reference.split('_')[1]

        # Skip entries that are too ambiguous (e.g. 9ZZZZ) or recently moved by UniProt
        if not mnemonic in lineageLookup: continue
        if not mnemonic in speciesLookup: continue

        # Combine lineage and species for lookup into taxonomy results
        lineage = [item.strip() for item in lineageLookup[mnemonic].split(';')]
        lineage.append(speciesLookup[mnemonic])

        for rank in lineage:
            if rank in passData:
                query = samEntries[entry].query
                if not query in retData: retData[query] = list()
                retData[query].append(rank)

    return retData

# Render SAM reads per taxonomic rank
def renderSam(passData):
    plugins.core.sendOutput("Read\tTaxonomy")

    for read in passData:
        for rank in passData[read]:
            outputLine = "{0}\t{1}".format(read, rank)
            plugins.core.sendOutput(outputLine)
