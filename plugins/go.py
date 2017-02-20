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

# Perform gene ontology term analysis on a PALADIN UniProt report

import os
import plugins.core
import argparse
import shlex
import operator

# Plugin connection definition
def pluginConnect(passDefinition):
    passDefinition.name = "go"
    passDefinition.description = "Perform gene ontology term grouping and abundance reporting"
    passDefinition.versionMajor = 1
    passDefinition.versionMinor = 0
    passDefinition.versionRevision = 0

    #passDefinition.callbackInit = goInit
    passDefinition.callbackMain = goMain

# Plugin main    
def goMain(passArguments):
    # Parse arguments 
    argParser = argparse.ArgumentParser(description='PALADIN Pipeline Plugins: GO', prog='go')
    argParser.add_argument('-i', dest='input', type=str, required=True, help='PALADIN UniProt report')
    argParser.add_argument('-q', dest='quality', type=int, required=True, help='Minimum mapping quality filter')
    
    arguments = argParser.parse_known_args(shlex.split(passArguments))

    # Get entries
    entries = plugins.core.PaladinEntry.getEntries(arguments[0].input, arguments[0].quality)

    # Aggregate for ontology
    ontologyData = aggregateOntology(entries)

    # Render
    renderAbundance(ontologyData, '', 0)

# Group and count GO terms
def aggregateOntology(passEntries):
    retOntology = dict()

    for entry in passEntries:
        ontology = passEntries[entry].ontology
        if not ontology or ontology[0].rstrip() == '': continue

        for term in ontology:
            if not term in retOntology: retOntology[term] = 0
            retOntology[term] += passEntries[entry].count

    return retOntology

# Render the abundance of the GO terms
def renderAbundance(passData, passHeader, passTotal):
    sortedData = sorted(passData.items(), key = operator.itemgetter(1), reverse=True)

    for row in sortedData:
        outputLine = "{0}\t{1}".format(row[0], row[1])
        plugins.core.sendOutput(outputLine)
