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

# Download custom UniProt reports

import os
import argparse
import shlex
import operator
import re
import requests
import plugins.core

UNIPROT_TEMPLATE_INIT = "http://www.uniprot.org/uploadlists/"
UNIPROT_TEMPLATE_JOB = "http://www.uniprot.org/jobs/{0}.stat"
UNIPROT_TEMPLATE_RESULTS = "http://www.uniprot.org/uniprot/?query=job:{0}&format=tab&columns={1}"
UNIPROT_REQUIRED = ['entry%20name', 'id']
UNIPROT_COLUMNS = ['organism', 'protein%20names', 'genes', 'pathway', 'features', 'go', 'reviewed', 'existence', 'comments', 'database(KEGG)', 'database(GeneID)', 'database(PATRIC)', 'database(EnsemblBacteria)']

# Plugin connection definition
def pluginConnect(passDefinition):
    passDefinition.name = "uniprot"
    passDefinition.description = "Download custom UniProt reports"
    passDefinition.versionMajor = 1
    passDefinition.versionMinor = 0
    passDefinition.versionRevision = 1
    passDefinition.dependencies=[]

    #passDefinition.callbackInit = uniprotInit
    passDefinition.callbackMain = uniprotMain

# Plugin main    
def uniprotMain(passArguments):
    # Parse arguments 
    argParser = argparse.ArgumentParser(description='PALADIN Pipeline Plugins: UniProt', prog='uniprot')
    argParser.add_argument('-i', dest='input', metavar='INPUT', type=str, required=True, help='Input SAM file')
    argParser.add_argument('-c', dest='columns', metavar='COLUMN', nargs='+', default=UNIPROT_COLUMNS, type=str, help='Columns to retrieve')
    argParser.add_argument('-b', dest='batch', metavar='BATCH', type=int, default=5000, help='Number of entries per submission batch')
    argParser.add_argument('-r', dest='retry', metavar='RETRY', type=int, default=10, help='Number of times to retry after http error')

    arguments = argParser.parse_known_args(shlex.split(passArguments))

    # Aggregate SAM data
    samData, samTotal = aggregateSam(arguments[0].input)

    # Retrieve UniProt info 
    columns = UNIPROT_REQUIRED
    columns.extend(arguments[0].columns)
    uniprotData = retrieveData(list(samData.keys()), columns, arguments[0].batch, arguments[0].retry)

    # Join SAM data with Uniprot
    joinedData = joinData(samData, samTotal, uniprotData)

    # Sort and render report
    headers = "Count\tAbundance\tQuality (Average)\tQuality (Max)\tUniProtKB\tID\t{0}".format("\t".join(uniprotData['Entry name'][2:]))
    joinSorted = sorted(joinedData, reverse=True)
    renderReport(joinSorted, headers)

# Group, count and average SAM data by reference hit
def aggregateSam(passSam):
    retAggregate = dict()
    retTotal = 0

    # Retrieve mapped SAM entries for report
    plugins.core.sendOutput('Gathering SAM data...', 'stderr') 
    entries = plugins.core.SamEntry.getEntries(passSam, 0)

    # Iterate through each read/hit
    for readEntry in entries:
        samEntry = entries[readEntry]
        samReference = samEntry.reference
        if '|' in samReference: samReference = samReference.split('|')[2]

        if not samReference in retAggregate:
            retAggregate[samReference] = [0, 0.0, 0]

        mapq = samEntry.mapqual
        retAggregate[samReference][0] += 1
        retAggregate[samReference][1] += mapq
        if mapq > retAggregate[samReference][2]:
            retAggregate[samReference][2] = mapq

        retTotal += 1
 
    # Finalize averages
    for entry in retAggregate:
        retAggregate[entry][1] /= retAggregate[entry][0]

    return retAggregate, retTotal

# Retrieve data on the requested reference hits from UniProt
def retrieveData(passEntries, passColumns, passBatch, passRetry):
    retData = dict()
 
    for batchIdx in range(0, len(passEntries), passBatch):
        # Ready this batch of entries
        batchEntries = ''
        for entryIdx in range(batchIdx, batchIdx + passBatch):
            if entryIdx >= len(passEntries): break
            batchEntries += "{0},".format(passEntries[entryIdx])

        # Construct POST data
        postData = dict()
        postData['uploadQuery'] = batchEntries[:-1]
        postData['format'] = 'job'
        postData['from'] = 'ACC+ID'
        postData['to'] = 'ACC'
        postData['landingPage'] = 'false'
        
        # Submit batch query, retrieve job ID
        retry = 0

        while True:
            try:
                endBatch = batchIdx + passBatch
                if endBatch > len(passEntries): endBatch = len(passEntries)
                plugins.core.sendOutput("Fetching entries {0}:{1} of {2}...".format(batchIdx, endBatch, len(passEntries)), 'stderr')

                response = requests.post(UNIPROT_TEMPLATE_INIT, postData)
                jobID = response.text
                response.close()

                # Check for valid job ID
                if any(x in jobID for x in ["<", " "]): raise

                # Monitor if job has completed
                url = UNIPROT_TEMPLATE_JOB.format(jobID)
                jobComplete = ''

                while jobComplete != 'COMPLETED':
                    response = requests.post(url)
                    jobComplete = response.text
                    response.close()

                # Fetch data and breakout text
                url = UNIPROT_TEMPLATE_RESULTS.format(jobID, ','.join(passColumns))
                response = requests.post(url)

                for line in response.text.split("\n"):
                    line = line.rstrip()
                    fields = line.split("\t")
                    retData[fields[0]] = fields
                    
                response.close()                 

                # End retry processing
                break
            except KeyboardInterrupt: raise
            except:
                if retry > passRetry: 
                    plugins.core.sendOutput('Too many HTTP errors, quitting...', 'stderr')
                    return ''
                else: plugins.core.sendOutput('HTTP error, retrying...', 'stderr')
                retry += 1
           
    return retData

# Join local data with UniProt data
def joinData(passSam, passTotal, passUniprot):
    retData = list()

    for kbid in passSam:
        samEntry = passSam[kbid]
        joinEntry = [samEntry[0], samEntry[0] / passTotal * 100, samEntry[1], samEntry[2]]

        # Check for UniProt data
        if kbid in passUniprot: joinEntry.extend(passUniprot[kbid])
        else: joinEntry.append(kbid)

        retData.append(joinEntry)

    return retData

# Render report                
def renderReport(passData, passHeader):
    plugins.core.sendOutput(passHeader)

    for entry in passData:
        plugins.core.sendOutput("\t".join([str(item) for item in entry]))
