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

# Core module is responsible for plugin interopability, as well as standard API

import os
import urllib.request
import pkgutil
import importlib
import plugins
import plugins.core
import pprint
import re

# Plugin definition to be populated by each plugin via its pluginInit method at load time
class PluginDef:
    def __init__(self, passModule):
        # Plugin fields
        self.module = passModule
        self.name = ""
        self.description = ""
        self.versionMajor = 0
        self.versionMinor = 0
        self.versionRevision = 0
        self.dependencies = list()

        # Callbacks
        self.callbackInit = None
        self.callbackMain = None

# Stores data per individual SAM entry
class SamEntry:
    FIELD_QNAME = 0
    FIELD_FLAG = 1
    FIELD_RNAME = 2
    FIELD_POS = 3
    FIELD_MAPQ = 4
    FIELD_CIGAR = 5
    FIELD_RNEXT = 6
    FIELD_PNEXT = 7
    FIELD_TLEN = 8
    FIELD_SEQ = 9
    FIELD_QUAL = 10

    def __init__(self):
        self.query = ''
        self.flag = 0
        self.reference = ''
        self.pos = 0
        self.mapqual = 0
        self.cigar = ''
        self.nextref = ''
        self.nextpos = 0
        self.length = 0
        self.sequence = ''
        self.readqual = 0
        self.frame = ''

    # getEntries is the public API for obtaining SAM data (will handle caching internally)
    @staticmethod
    def getEntries(passFile, passQuality):
        # Check if in cache
        if not (passFile, passQuality) in SamEntry._entryCache:
            cache = SamEntry.generateCache(passFile, passQuality)
            SamEntry._entryCache[(passFile, passQuality)] = cache

        return SamEntry._entryCache[(passFile, passQuality)]

    # Cache this SAM file for the requested quality
    @staticmethod
    def generateCache(passFile, passQuality):
        retEntries = dict()

        # Open SAM file
        with open(passFile, 'r') as fileHandle:
            for line in fileHandle:
                line = line.rstrip()
                fields = line.split("\t")

                # Skip header and malformed lines
                if line.startswith("@"): continue
                if len(fields) < 11: continue

                # Filter for minimum quality
                if fields[SamEntry.FIELD_RNAME] == '*' and passQuality != -1: continue
                if passQuality != -1 and int(fields[SamEntry.FIELD_MAPQ]) < passQuality: continue

                # Remove PALADIN frame header since best scoring frame may change between alignments
                headerMatch = re.search("(.*?:.*?:.*?:)(.*)", fields[SamEntry.FIELD_QNAME])

                entry = SamEntry()
                entry.query = headerMatch.group(2)
                entry.flag = int(fields[SamEntry.FIELD_FLAG])

                # Fill in entry information if mapped
                if entry.isMapped:
                    entry.reference = fields[SamEntry.FIELD_RNAME]
                    entry.pos = SamEntry.getSamInt(fields[SamEntry.FIELD_POS])
                    entry.mapqual = SamEntry.getSamInt(fields[SamEntry.FIELD_MAPQ])
                    entry.cigar = fields[SamEntry.FIELD_CIGAR]
                    entry.nextref = fields[SamEntry.FIELD_RNEXT]
                    entry.nextpos = fields[SamEntry.FIELD_PNEXT]
                    entry.length = SamEntry.getSamInt(fields[SamEntry.FIELD_TLEN])
                    entry.sequence = fields[SamEntry.FIELD_SEQ]
                    entry.readqual = SamEntry.getSamInt(fields[SamEntry.FIELD_QUAL])
                    entry.frame = headerMatch.group(1)

                # Each read can have multiple non-linear/chimeric hits - store as tuple for ease of processing
                readBase = headerMatch.group(2)
                hitIdx = 0
                while (readBase, hitIdx) in retEntries: hitIdx += 1
                retEntries[(readBase, hitIdx)] = entry
                    
        return retEntries
   
    @staticmethod
    def getSamInt(passInt):
        try: return int(passInt)
        except: return 0
 
    def isMapped(): return self.flag & 0x04 > 0 

    _entryCache = dict()   

# PALADIN UniProt entry
class PaladinEntry:
    FIELD_COUNT = 0
    FIELD_ABUNDANCE = 1
    FIELD_QUALAVG = 2
    FIELD_QUALMAX = 3
    FIELD_KB = 4
    FIELD_ID = 5
    FIELD_SPECIES = 6
    FIELD_PROTEIN = 7
    FIELD_ONTOLOGY = 11
   
    TYPE_UNKNOWN = 0
    TYPE_UNIPROT_EXACT = 1
    TYPE_UNIPROT_GROUP = 2
    TYPE_CUSTOM = 3
 
    def __init__(self):
        self.type = PaladinEntry.TYPE_UNKNOWN
        self.id = 'Unknown'
        self.kb = 'Unknown'
        self.count = 0
        self.abundance = 0.0
        self.qualityAvg = 0.0
        self.qualityMax = 0
        self.speciesID = 'Unknown'
        self.speciesFull = 'Unknown'
        self.ontology = list()

    # getEntries is the public API for obtaining UniProt report data (will handle caching internally)
    @staticmethod
    def getEntries(passFile, passQuality, passPattern=None):
        # Check if in cache
        if not (passFile, passQuality, passPattern) in PaladinEntry._entryCache: 
            cache = PaladinEntry.generateCache(passFile, passQuality, passPattern)
            PaladinEntry._entryCache[(passFile, passQuality, passPattern)] = cache
        
        return PaladinEntry._entryCache[(passFile, passQuality, passPattern)]

    @staticmethod 
    def generateCache(passFile, passQuality, passPattern):
        retEntries = dict()

        # Open UniProt report, skip header
        with open(passFile, 'r') as fileHandle:
            fileHandle.readline()

            for line in fileHandle:
                line = line.rstrip()
                fields = line.split("\t")

                # Filter for minimum quality
                if float(fields[PaladinEntry.FIELD_QUALMAX]) < passQuality: continue

                entry = PaladinEntry()
                entry.count = int(fields[PaladinEntry.FIELD_COUNT])
                entry.abundance = float(fields[PaladinEntry.FIELD_ABUNDANCE])
                entry.qualAvg = float(fields[PaladinEntry.FIELD_QUALAVG])
                entry.qualMax = int(fields[PaladinEntry.FIELD_QUALMAX])
                entry.kb = fields[PaladinEntry.FIELD_KB]

                if len(fields) > 10:
                    # Existence of fields indicates a successful UniProt parse by PALADIN
                    if "_9" in entry.kb: entry.type = PaladinEntry.TYPE_UNIPROT_GROUP
                    else: entry.type = PaladinEntry.TYPE_UNIPROT_EXACT 

                    entry.speciesID = entry.kb.split('_')[1]
                    entry.speciesFull = fields[PaladinEntry.FIELD_SPECIES]
                    entry.id = fields[PaladinEntry.FIELD_ID]
                    entry.ontology = [term.strip() for term in fields[PaladinEntry.FIELD_ONTOLOGY].split(';')]
                else:
                    # Check for custom match
                    if passPattern:
                        match = re.search(passPattern, entry.kb)
                        if match: 
                            entry.type = PaladinEntry.TYPE_CUSTOM
                            entry.speciesID = match.group(1)
                            entry.speciesFull = match.group(1)
      
                retEntries[fields[PaladinEntry.FIELD_KB]] = entry

        return retEntries

    _entryCache = dict()    

# Plugins internal to core, and loaded external modules
internalPlugins = dict()
pluginModules = dict()

# Standard output and error buffers
outputStdout = list()
outputStderr = list()
consoleStdout = False
consoleStderr = True

# Search for all modules in the plugin package (directory), import each and run pluginConnect method
def connectPlugins():
    # Ready cache directory
    getCacheDir()

    # Add internal core plugins
    internalPlugins['flush'] = renderOutput
    internalPlugins['write'] = renderOutput

    # Import all external plugin modules in package (using full path)
    for importer, module, package in pkgutil.iter_modules(plugins.__path__):
        if module == "core": continue
        moduleHandle = importlib.import_module("{0}.{1}".format(plugins.__name__, module))
        pluginModules[module] = PluginDef(moduleHandle)

    # Connect to all external plugins
    for plugin in pluginModules:
        pluginModules[plugin].module.pluginConnect(pluginModules[plugin])

# Initialize plugins being used in this session
def initPlugins(passPlugins):
    initQueue = set()
    initHistory = set()

    # Scan for plugins and dependencies 
    for plugin in passPlugins:
        if not plugin in pluginModules and not plugin in internalPlugins: 
            if plugin in plugins.modulesDisabled: print("Disabled plugin: {0}".format(plugin))
            else: print("Unknown plugin: {0}".format(plugin))
            return False

        # Initialize external plugins
        if plugin in pluginModules:
            # Look for dependencies
            initQueue.update(pluginModules[plugin].dependencies)
            initQueue.add(plugin)

    # Initialize
    for plugin in initQueue:
        if pluginModules[plugin].callbackInit:
            if not plugin in initHistory:
                pluginModules[plugin].callbackInit()
                initHistory.add(plugin)

    return True

# Execute requested plugin pipeline
def execPipeline(passPipeline):
    for task in passPipeline:
        if task[0] in internalPlugins:
            # Internal plugin
            internalPlugins[task[0]](task[1].strip('"'))
        else: 
            # External plugin
            if pluginModules[task[0]].callbackMain: 
                pluginModules[task[0]].callbackMain(task[1])

# The flush internal plugin handles rendering output (to stdout or file)
def renderOutput(passFile = '', passTarget = 'stdout'):
    if passTarget == 'stdout':
        renderString = ''.join(plugins.core.outputStdout)
        del plugins.core.outputStdout[:]
    if passTarget == 'stderr':
        renderString = ''.join(plugins.core.outputStderr)
        del plugins.core.outputStderr[:]

    if not passFile:
        print(renderString, flush=True)
    else:
        with open(passFile, 'w') as fileHandle: 
            fileHandle.write(renderString)

# API - Record output into the appropriate buffer
def sendOutput(passOutput, passTarget='stdout', passEnd="\n"):
    newContent = "{0}{1}".format(passOutput, passEnd)

    # Check if output should be buffered or set to console
    console = False
    if passTarget == 'stdout' and plugins.core.consoleStdout == True: console = True
    if passTarget == 'stderr' and plugins.core.consoleStderr == True: console = True

    if console:
        print(newContent, end='', flush=True)
    else:
        if passTarget == 'stdout': plugins.core.outputStdout.append(newContent)
        if passTarget == 'stderr': plugins.core.outputStderr.append(newContent)

# API - Save HTTP resonse data to disk
def downloadURL(passURL, passFile, passDescription, passForce=False, passRender=True):
    if not os.path.exists(passFile) or passForce:
        if passRender:
            plugins.core.sendOutput("Downloading {0}...".format(passDescription), 'stderr')

        response = urllib.request.urlopen(passURL)
        fileData = response.read()
        with open(passFile, 'wb') as fileHandle:
            fileHandle.write(fileData)

# API - Get cache/settings directory
def getCacheDir(passDefinition=None):
    sub = ''
    if passDefinition: sub = passDefinition.name

    fullDirectory = "{0}/.paladin-plugins/{1}".format(os.path.expanduser("~"), sub)

    # Check for existence
    if not os.path.exists(fullDirectory):
        os.makedirs(fullDirectory)

    return fullDirectory

# API - Return value if string is integer (allows negatives)
def getInteger(passString):
    try: return int(passString)
    except: return None 
        
# API - Debugging
def debugPrint(passObject):
    pp = pprint.PrettyPrinter(indent=4)
    pp.pprint(passObject)
