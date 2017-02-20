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

# Analyze difference between two PALADIN alignemnts

import os
import argparse
import shlex
import operator
import re
import plugins.core

# Plugin connection definition
def pluginConnect(passDefinition):
    passDefinition.name = "difference"
    passDefinition.description = "Analyze relative differences between two PALADIN taxonomy reports"
    passDefinition.versionMajor = 1
    passDefinition.versionMinor = 0
    passDefinition.versionRevision = 0
    passDefinition.dependencies=['taxonomy']

    #passDefinition.callbackInit = proteinInit
    passDefinition.callbackMain = differenceMain
   
# Plugin main 
def differenceMain(passArguments):
    # Parse arguments 
    argParser = argparse.ArgumentParser(description='PALADIN Pipeline Plugins: Difference', prog='difference')
    argParser.add_argument('-1', dest='basisFiles', metavar=('Taxonomy', 'Sam', 'Uniprot'), type=str, nargs=3, required=True, help='Basis files of comparison (Taxonomy, Sam, UniProt)')
    argParser.add_argument('-2', dest='compareFiles', metavar=('Taxonomy', 'Sam', 'Uniprot'), type=str, nargs=3, required=True, help='Compared files of comparison (Taxonomy, Sam, UniProt')
    argParser.add_argument('-c', dest='custom', type=str, default='', help='Species-parsing regex pattern for non-UniProt entries')
    argParser.add_argument('-s', dest='simple', action='store_true', help='Simple report, do not include mapped/unmapped read contributions (alters totals)')
    arguments = argParser.parse_known_args(shlex.split(passArguments))

    runningData = list()

    # Obtain taxonomy differences
    plugins.core.sendOutput('Comparing taxonomy reports...', 'stderr')
    taxonomyData, taxonomyTotal = taxonomyCompare((arguments[0].basisFiles[0], arguments[0].compareFiles[0]))
        
    # Obtain SAM differences with Uniprot info for full report
    samData = list()
    uniprotData = list()
    if not arguments[0].simple:
        plugins.core.sendOutput('Gathering SAM data from first alignment', 'stderr')
        samEntries1 = plugins.core.SamEntry.getEntries(arguments[0].basisFiles[1], -1)
        plugins.core.sendOutput('Gathering SAM data from second alignment', 'stderr')
        samEntries2 = plugins.core.SamEntry.getEntries(arguments[0].compareFiles[1], -1)
        plugins.core.sendOutput('Comparing SAM data...', 'stderr')
        samData = samCompare(samEntries1, samEntries2)

        # Get Uniprot entries
        plugins.core.sendOutput('Gathering UniProt data from first alignment', 'stderr')
        uniprotData.append(plugins.core.PaladinEntry.getEntries(arguments[0].basisFiles[2], 0, arguments[0].custom))
        plugins.core.sendOutput('Gathering UniProt data from second alignment', 'stderr')
        uniprotData.append(plugins.core.PaladinEntry.getEntries(arguments[0].compareFiles[2], 0, arguments[0].custom))

    plugins.core.sendOutput('Comparing alignments...', 'stderr')
    combinedData = combinedCompare(taxonomyData, samData, uniprotData, plugins.taxonomy.lineageLookup, taxonomyTotal)
    renderCombined(combinedData, arguments[0].simple)

# Taxonomy report comparison
def taxonomyCompare(passReports):
    retDifference = dict()
    retTotal = 0
    reportData = (dict(), dict())

    # Read files into memory
    for idx in range(2):
        with open(passReports[idx], 'r') as fileHandle:
            # Skip header
            fileHandle.readline()

            for line in fileHandle:
                line = line.rstrip()
                fields = line.split("\t")
                reportData[idx][fields[2]] = int(fields[0])
                if idx == 1: retTotal += int(fields[0])

    # Calculate differences (set 1, set 2, shared)
    set1Keys = set(reportData[0].keys())
    set2Keys = set(reportData[1].keys())
    sharedKeys = set1Keys.intersection(set2Keys)
    set1Keys = set1Keys.difference(sharedKeys)
    set2Keys = set2Keys.difference(sharedKeys)
    
    retDifference = {key: (reportData[1][key] - reportData[0][key]) for key in sharedKeys}
    retDifference.update({key: -reportData[0][key] for key in set1Keys})
    retDifference.update({key: reportData[1][key] for key in set2Keys})

    return retDifference, retTotal

# SAM entries comparison
def samCompare(passEntries1, passEntries2):
    retDiff = dict()

    # Use set 1 as the baseRead 
    for entry in passEntries1:
        # Iterate through possible supplementary hits until neither available
        readIdx = (entry[0], 0)
        while readIdx in passEntries1 or readIdx in passEntries2:
            if not readIdx in passEntries1:
                ref1 = '*'
                ref2 = passEntries2[readIdx].reference
            elif not readIdx in passEntries2:
                ref1 = passEntries1[readIdx].reference
                ref2 = '*'
            else:
                ref1 = passEntries1[readIdx].reference
                ref2 = passEntries2[readIdx].reference

            # Note references that do not match 
            if ref1 != ref2:
                diffPair = (ref1, ref2)
                if not diffPair in retDiff: retDiff[diffPair] = 0
                retDiff[diffPair] += 1

            readIdx = (readIdx[0], readIdx[1] + 1)

    return retDiff

# Combined compare breaks down taxonomy report on per read basis, then aggregates
def combinedCompare(passTaxonomyEntries, passSamEntries, passUniprotEntries, passLineageLookup, passTotal):
    retEntries = dict()
    
    # Calculate unmapped difference and add to taxonomy set
    unmappedDiff = 0 
    for samEntry in passSamEntries:
        for setIdx in range(2):
            if samEntry[setIdx] == '*': 
                if setIdx == 0: passTotal += 1
                unmappedDiff += [1, -1][setIdx]

    # Add unmapped entry to taxonomy list
    passTaxonomyEntries['Unmapped'] = unmappedDiff

    # Iterate through each taxonomy grouping 
    for taxonEntry in passTaxonomyEntries:
        # Skip unchanged entries
        if passTaxonomyEntries[taxonEntry] == 0: continue

        # Prepare entry if new
        if not taxonEntry in retEntries: 
            retEntries[taxonEntry] = (dict(), passTaxonomyEntries[taxonEntry]/passTotal)

        # Scan SAM entries for matching lineage
        for samEntry in passSamEntries:
            # Lookup species and lineage
            samRecords = list()
            lineage = list()
            
            for setIdx in range(2):
                if samEntry[setIdx] == '*': 
                    missRecord = type('missRecord', (object,), {})()
                    missRecord.speciesID = 'Unmapped'
                    missRecord.speciesFull = 'Unmapped'
                    samRecords.append(missRecord)
                    lineage.append('Unmapped')
                else:
                    samRecords.append(passUniprotEntries[setIdx][samEntry[setIdx].split('|')[-1]])
                    if samRecords[setIdx].speciesID in passLineageLookup:
                        lineage.append([x.strip() for x in passLineageLookup[samRecords[setIdx].speciesID].split(';')])
                    else:
                        lineage.append("Unknown")
                    
            for setIdx in range(2):
                # Attempt to match on species
                if samRecords[setIdx].speciesFull == taxonEntry: 
                    samSpecies = samRecords[1-setIdx].speciesFull
                    
                    samAmount = [1, -1][setIdx] * passSamEntries[samEntry]

                    if not samSpecies in retEntries[taxonEntry][0]: retEntries[taxonEntry][0][samSpecies] = (0, 0)
                    samParts = retEntries[taxonEntry][0][samSpecies]
                    
                    if samAmount > 0: retEntries[taxonEntry][0][samSpecies] = (samParts[0] + samAmount, samParts[1])
                    else: retEntries[taxonEntry][0][samSpecies] = (samParts[0], samParts[1] + samAmount) 

                    break

                # Attempt to match on lineage
                for rankIdx in range(len(lineage[setIdx])):
                    if lineage[setIdx][rankIdx] == taxonEntry:
                        compareIdx = rankIdx
                        if compareIdx >= len(lineage[1-setIdx]): compareIdx = len(lineage[1-setIdx]) - 1
                        samLineage = lineage[1-setIdx][compareIdx]
                        samAmount = [1, -1][setIdx] * passSamEntries[samEntry]

                        if not samLineage in retEntries[taxonEntry][0]: retEntries[taxonEntry][0][samLineage] = (0, 0)
                        samParts = retEntries[taxonEntry][0][samLineage]

                        if samAmount > 0: retEntries[taxonEntry][0][samLineage] = (samParts[0] + samAmount, samParts[1])
                        else: retEntries[taxonEntry][0][samLineage] = (samParts[0], samParts[1] + samAmount)
                        break

    # Calculate percentages
    for taxonomyEntry in retEntries:
        for samEntry in retEntries[taxonomyEntry][0]:
            samParts = retEntries[taxonomyEntry][0][samEntry] 
            retEntries[taxonomyEntry][0][samEntry] = (samParts[0]/passTotal, samParts[1]/passTotal) 
 
    return retEntries

# Align stats returns 3-tuple (Different, Added, Removed)
def alignStats(passEntries1, passEntries2):
    totalDiff = 0
    totalAdd = 0
    totalDel = 0 

    for entry in passEntries1:
        if passEntries1[entry].flag & 0x4 == passEntries2[entry].flag & 0x4:
            # Matching map types
            if passEntries1[entry].reference != passEntries2[entry].reference: 
                totalDiff += 1
        else:
            if passEntries1[entry].flag & 0x4 > 0: totalAdd += 1
            else: totalDel += 1

    return (totalDiff, totalAdd, totalDel)

# Render combined report
def renderCombined(passData, passSimple):
    # Render appropriate header
    if passSimple: header = "Taxon\tDifference"
    else: header = "Taxon\tDifference\tContributor\tContribution (Pos)\tContribution (Neg)"
    plugins.core.sendOutput(header)

    sortedTaxa = sorted(passData.items(), key = lambda item: abs(item[1][1]), reverse=True)

    for taxonEntry in sortedTaxa:
        if passSimple:
            plugins.core.sendOutput("{0}\t{1}".format(taxonEntry[0], taxonEntry[1][1]))
        else: 
            # Sort SAM data
            sortedSam = sorted(taxonEntry[1][0].items(), key = lambda item: abs(item[1][0] + item[1][1]), reverse=True)
            for samEntry in sortedSam:
                plugins.core.sendOutput("{0}\t{1}\t{2}\t{3}\t{4}".format(taxonEntry[0], taxonEntry[1][1], samEntry[0], samEntry[1][0], samEntry[1][1]))
