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

# Distribute PALADIN execution across a cluster

import os
import argparse
import shlex
import operator
import re
import subprocess
import gzip
import plugins.core

# Plugin connection definition
def pluginConnect(passDefinition):
    passDefinition.name = "hpc"
    passDefinition.description = "Distribute PALADIN execution across a cluster"
    passDefinition.versionMajor = 1
    passDefinition.versionMinor = 0
    passDefinition.versionRevision = 0
    passDefinition.dependencies=['aggregation']

    #passDefinition.callbackInit = hpcInit
    passDefinition.callbackMain = hpcMain

# Plugin main    
def hpcMain(passArguments):
    # Parse arguments 
    argParser = argparse.ArgumentParser(description='PALADIN Pipeline Plugins: HPC', prog='hpc')
    argParser.add_argument('reference', metavar='REFERENCE', type=str, help='Reference database')
    argParser.add_argument('input', metavar='INPUT', type=str, help='Input reads')
    argParser.add_argument('output', metavar='OUTPUT', type=str, help='Output name')
    argParser.add_argument('options', metavar='OPTIONS', type=str, nargs=argparse.REMAINDER, help='PALADIN options')
    arguments = argParser.parse_known_args(shlex.split(passArguments))

    # Obtain MPI process info (import MPI here to avoid crashing on non-MPI systems)
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    processCount = comm.Get_size()
    processIdx = comm.Get_rank()

    # Process 0 is responsible for splitting reads
    if processIdx == 0: splitReads(arguments[0].input, processCount)

    # Sync to ensure all batches of reads ready for alignment
    comm.Barrier()

    # Execute PALADIN alignment
    processReads = "{0}-{1}".format(arguments[0].input, processIdx)
    processOutput = "{0}-{1}".format(arguments[0].output, processIdx)
    execute(arguments[0].reference, processReads, processOutput, ' '.join(arguments[0].options))

    # Sync to ensure all alignments complete
    comm.Barrier()

    # Aggregate result
    if processIdx == 0:
        inputPaths = [ "{0}-{1}".format(arguments[0].output, idx) for idx in range(processCount)]
        plugins.aggregation.processData(inputPaths, arguments[0].output)

        # Remove process specific files
        for idx in range(processCount):
            os.remove("{0}-{1}".format(arguments[0].input, idx))
            os.remove("{0}-{1}.sam".format(arguments[0].output, idx))
            os.remove("{0}-{1}_uniprot.tsv".format(arguments[0].output, idx))

# Execute PALADIN in the current process
def execute(passReference, passInput, passOutput, passOptions):
    command = "paladin align {0} {1} -o {2} {3}".format(passReference, passInput, passOutput, passOptions)
    output = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True)

    with open("{0}.log".format(passOutput), 'wb') as fileHandle:
        fileHandle.write(output)

# Split reads across the requested node count
def splitReads(passReads, passCount):
    # Treat as gzip file with gz extension
    if passReads.endswith('gz'): inputHandle = gzip.open(passReads, 'rb') 
    else: inputHandle =  open(passReads, 'rb') 

    # Open reads and outputs
    with inputHandle:
        outputHandles = list() 
        for outIdx in range(passCount):
            outputHandles.append(open("{0}-{1}".format(passReads, outIdx), 'w'))

        mode = -1
        outIdx = -1
        parseIdx = -1
        for line in inputHandle:
            line = line.decode('utf-8')
            if line.rstrip() == '': continue

            # Initially detect input type (0 = fasta, 1 = fastq)
            if mode == -1: mode = (1, 0)[line.startswith('>')] 
           
            if mode == 0:
                # FastA mode
                if line.startswith('>'): outIdx = (outIdx + 1) % passCount
                outputHandles[outIdx].write(line)
            else:
                # FastQ mode
                parseIdx = (parseIdx + 1) % 4
                if parseIdx == 0: outIdx = (outIdx + 1) % passCount
                outputHandles[outIdx].write(line) 

        # Close output handles         
        for handle in outputHandles: handle.close()
