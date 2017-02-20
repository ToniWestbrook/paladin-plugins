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

# Generate a reference containing each protein sequence used to generate the UniRef clusters detected in a PALADIN UniProt report

import os
import gzip
import argparse
import shlex
import plugins.core

DECLUSTER_DOWNLOADS = [('UniProt ID mappings', 'idmapping.dat.gz', 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz'),
                    ('UniProt Swissprot', 'uniprot_sprot.fasta.gz', 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz'),
                    ('UniProt TrEMBL', 'uniprot_trembl.fasta.gz', 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz')]
DECLUSTER_UNIREF90 = 'idmapping-uniref90.dat'

moduleDefinition = None
idLookup = dict()
clusterLookup = dict()
sequenceLookup = dict()

# Plugin connection definition
def pluginConnect(passDefinition):
    passDefinition.name = 'decluster'
    passDefinition.description = 'Generate a reference containing each protein sequence used to generate the UniRef clusters detected in a PALADIN UniProt report'
    passDefinition.versionMajor = 1
    passDefinition.versionMinor = 0
    passDefinition.versionRevision = 0

    passDefinition.callbackInit = declusterInit
    passDefinition.callbackMain = declusterMain
   
    plugins.decluster.moduleDefinition = passDefinition

# Plugin initialization
def declusterInit():
    # Download UniProt and process data
    downloadUniprot()
    extractUniref()
    cacheLookups()

# Plugin main
def declusterMain(passArguments):
    # Parse arguments 
    argParser = argparse.ArgumentParser(description='PALADIN Pipeline Plugins: Decluster', prog='decluster')                                                                                                                             
    argParser.add_argument('-i', dest='input', metavar='INPUT', type=str, required=True, help='PALADIN UniProt report')
    argParser.add_argument('-q', dest='quality', type=int, help='Minimum mapping quality filter')
    arguments = argParser.parse_known_args(shlex.split(passArguments))

    # Get entries
    entries = plugins.core.PaladinEntry.getEntries(arguments[0].input, arguments[0].quality)
    entryIDs = [entries[entry].id for entry in entries]
    plugins.core.sendOutput("Parsed {0} entries from UniProt report".format(len(entryIDs)), 'stderr')

    # Get assocaited cluster with each entry
    clusters = getClusterInfo(entryIDs, 'id')
    plugins.core.sendOutput("Parsed {0} clusters from ID mappings".format(len(clusters)), 'stderr')

    # Get members for each cluster
    members = getClusterInfo(clusters, 'cluster')
    plugins.core.sendOutput("Parsed {0} members from ID mappings".format(len(members)), 'stderr')

    # Render sequences for all cluster members
    renderSequences(members)

# Get ID maps, and swissprot/trembl sequences 
def downloadUniprot():
    # Download and save files if necessary
    for download in DECLUSTER_DOWNLOADS:
        downloadPath = os.path.join(plugins.core.getCacheDir(plugins.decluster.moduleDefinition), download[1])
        plugins.core.downloadURL(download[2], downloadPath, download[0])

# Parse ID mappings for relevant data
def extractUniref():
    # Check if Uniref extract
    idPath = os.path.join(plugins.core.getCacheDir(plugins.decluster.moduleDefinition), DECLUSTER_DOWNLOADS[0][1]) 
    unirefPath = os.path.join(plugins.core.getCacheDir(plugins.decluster.moduleDefinition), DECLUSTER_UNIREF90)

    if os.path.exists(unirefPath): return

    with gzip.open(idPath, 'rb') as inputHandle:
        with open(unirefPath, 'w') as outputHandle:
   
            plugins.core.sendOutput("Extracting UniRef90 entries...", 'stderr') 
            # Filter for Uniref90 entries 
            for line in inputHandle:
                line = line.decode('utf-8')
                if not 'UniRef90' in line: continue
                
                outputHandle.write(line)

# Cache ID lookup data
def cacheLookups():
    # Clear cache
    idLookup.clear()
    clusterLookup.clear()
    sequenceLookup.clear()

    # Open uniref extract:
    extractPath = os.path.join(plugins.core.getCacheDir(plugins.decluster.moduleDefinition), DECLUSTER_UNIREF90)
    with open(extractPath, 'r') as fileHandle:

        # Cache UniRef90 Lookups
        plugins.core.sendOutput("Caching UniRef90 mappings...", 'stderr')

        for line in fileHandle:
            line = line.rstrip()
            fields = line.split("\t")
        
            # ID lookup is one-to-one
            idLookup[fields[0]] = fields[2]

            # Cluster lookup is one-to-many
            if not fields[2] in clusterLookup: clusterLookup[fields[2]] = list()
            clusterLookup[fields[2]].append(fields[0])

        # Open each reference
        for reference in DECLUSTER_DOWNLOADS[1:]:
            plugins.core.sendOutput("Caching {0}...".format(reference[1]), 'stderr')

            refPath = os.path.join(plugins.core.getCacheDir(plugins.decluster.moduleDefinition), reference[1])
            fileHandle = gzip.open(refPath, 'rb')
        
            activeID = ''
            for line in fileHandle:
                line = line.decode('utf-8').rstrip()

                if line.startswith('>'):
                    # Record new active ID  
                    activeID = line[line.index('|')+1:]
                    activeID = activeID[:activeID.index('|')]

                    # Start new sequence record
                    sequenceLookup[activeID] =  "{0}\n".format(line)
                else:
                    sequenceLookup[activeID] += "{0}\n".format(line)

# Get cluster entries        
def getClusterInfo(passLookup, passMode):
    retResults = list()

    for item in passLookup:
        if passMode == 'id': 
            if item in idLookup: retResults.append(idLookup[item])
        if passMode == 'cluster': 
            if item in clusterLookup: retResults += clusterLookup[item]

    return retResults

# Render sequences
def renderSequences(passMembers):
    for member in passMembers:
        plugins.core.sendOutput(sequenceLookup[member], 'stdout', '')
