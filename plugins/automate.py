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

# Automate PALADIN execution across multiple sets of reads

import os
import re
import subprocess
import argparse
import shlex
import plugins.core

# Plugin connection definition
def pluginConnect(passDefinition):
    passDefinition.name = 'automate'
    passDefinition.description = 'Automate PALADIN execution across multiple sets of reads'
    passDefinition.versionMajor = 1
    passDefinition.versionMinor = 0
    passDefinition.versionRevision = 1

    #passDefinition.callbackInit = automateInit
    passDefinition.callbackMain = automateMain

# Plugin main    
def automateMain(passArguments):
    # Parse arguments 
    argParser = argparse.ArgumentParser(description='PALADIN Pipeline Plugins: Automate', prog='automate')
    argParser.add_argument('reference', metavar='REFERENCE', type=str, help='Reference database')
    argParser.add_argument('root', metavar='ROOT', type=str, help='Root path to search') 
    argParser.add_argument('pattern', metavar='PATTERN', type=str, help='Input reads search pattern')
    argParser.add_argument('options', metavar='OPTIONS', type=str, nargs=argparse.REMAINDER, help='PALADIN options')
    arguments = argParser.parse_known_args(shlex.split(passArguments))

    for root, dirs, files in os.walk(arguments[0].root):
        for fileName in files:
            if not re.search(arguments[0].pattern, fileName): continue

            # Matching input sequence, execute PALADIN
            baseName = fileName
            if '.' in baseName: baseName = baseName[:baseName.index('.')]
            baseName = os.path.join(root, baseName)
            fullFile = os.path.join(root, fileName) 

            command = "paladin align {0} {1} -o {2} {3}".format(arguments[0].reference, fullFile, baseName, ' '.join(arguments[0].options))
            output = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True)

            with open("{0}.log".format(baseName), 'wb') as fileHandle:
                fileHandle.write(output)
