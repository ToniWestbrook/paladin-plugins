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

# Combine multiple PALADIN outputs (both SAM and UniProt) into a single output

import os
import re
import operator
import argparse
import shlex
import plugins.core

# Plugin connection definition
def pluginConnect(passDefinition):
    passDefinition.name = 'aggregation'
    passDefinition.description = 'Combine multiple PALADIN outputs (both SAM and UniProt) into a single output'
    passDefinition.versionMajor = 1
    passDefinition.versionMinor = 0
    passDefinition.verionRevision = 0

    #passDefinition.callbackInit = aggregationInit
    passDefinition.callbackMain = aggregationMain

# Plugin main
def aggregationMain(passArguments):
    # Parse arguments 
    argParser = argparse.ArgumentParser(description='PALADIN Pipeline Plugins: Aggregation', prog='aggregation')
    argParser.add_argument('-r', dest='rootPath', metavar='ROOT', type=str, required=True, help='Root path to search')
    argParser.add_argument('-s', dest='searchPattern', metavar='PATTERN', type=str, required=True, help='SAM search pattern (infers TSV)')
    argParser.add_argument('-o', dest='output', metavar='OUTPUT', type=str, required=True, help='Output path')
    arguments = argParser.parse_known_args(shlex.split(passArguments))   

    # Build base paths
    basePaths = buildPaths(arguments[0].rootPath, arguments[0].searchPattern)

    # Aggregate data
    processData(basePaths, arguments[0].output)

# Aggregate both SAM and TSV data (may be used as API)
def processData(passBasePaths, passOutput):
    #Aggregate SAM data and write
    samData = aggregateSAM(passBasePaths)
    with open("{0}.sam".format(passOutput), 'w') as fileHandle:
        fileHandle.write(samData)
 
    # Aggregate TSV data and write
    tsvHeader, tsvEntries = aggregateTSV(passBasePaths)
    tsvEntries = finalizeTSV(tsvEntries) 

    with open("{0}_uniprot.tsv".format(passOutput), 'w') as fileHandle: 
        fileHandle.write(tsvHeader + "\n")
        for entry in tsvEntries:
            fileHandle.write("\t".join([str(strEntry) for strEntry in entry]) + "\n")

# Recurse and build path information using SAM pattern
def buildPaths(passRoot, passPattern):
    retEntries = list()

    # Iterate through all subdirectories/files, match pattern
    for root, dirs, files in os.walk(passRoot):
        for fileName in files:
            if not re.match(passPattern, fileName):continue

            # Remove SAM extension, get base file path
            baseFile = os.path.join(root, fileName)[:-4]
            retEntries.append(baseFile)

    return retEntries

# Combine SAM files
def aggregateSAM(passPaths):
    entries = list()
    headerFile = True

    plugins.core.sendOutput('Gathering SAM data...', 'stderr')

    for path in passPaths:
        path = "{0}.sam".format(path)
        plugins.core.sendOutput("Processing {0}...".format(path), 'stderr')

        with open(path, 'r') as fileHandle:
            for line in fileHandle:
                if headerFile: 
                    # First file will include headers
                    entries.append(line)
                else:
                    # Subsequence files will not
                    if line.startswith('@'): continue
                    entries.append(line)

            headerFile = False

    return ''.join(entries)

# Group TSV files
def aggregateTSV(passPaths):
    retEntries = dict()
    retHeader = ''

    plugins.core.sendOutput('Gathering UniProt data...', 'stderr')

    for path in passPaths:
        path = "{0}_uniprot.tsv".format(path)
        plugins.core.sendOutput("Processing {0}...".format(path), 'stderr')

        # Open individual report, save header
        with open(path, 'r') as fileHandle:
            retHeader = fileHandle.readline().rstrip()
                        
            for line in fileHandle:
                line = line.rstrip()
                fields = line.split("\t")
                #extended = "\t".join(fields[6:])

                # Check for first instance of entry
                existingEntry = (0, 0.0, 0.0, 0)
                if fields[4] in retEntries: existingEntry = retEntries[fields[4]]

                # Update and store
                lineCount = int(fields[0])
                lineQualMax = int(fields[3])

                newCount = existingEntry[0] + lineCount
                newAbundance = existingEntry[1] + (lineCount * float(fields[1]))
                newQualAvg = existingEntry[2] + (lineCount * float(fields[2]))
                newQualMax = existingEntry[3]
                if lineQualMax > newQualMax: newQualMax = lineQualMax

                retEntries[fields[4]] = (newCount, newAbundance, newQualAvg, newQualMax) + tuple(fields[4:])

    return retHeader, retEntries

# Recalculate and sort TSV data
def finalizeTSV(passEntries):
    retEntries = list()

    # Summate counts
    totalCount = 0
    for entry in passEntries:
        totalCount += passEntries[entry][0]

    # Sort
    entriesSorted = sorted(passEntries.items(), key = operator.itemgetter(1), reverse=True)

    # Calculate means and store new entry
    for entry in entriesSorted:
        abundance = entry[1][0] / totalCount * 100
        qualityAvg = entry[1][2] / entry[1][0]
        newEntry = (entry[1][0], abundance, qualityAvg) + entry[1][3:]
        
        retEntries.append(newEntry) 

    return retEntries

