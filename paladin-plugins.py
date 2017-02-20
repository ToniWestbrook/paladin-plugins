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

import argparse
import operator
import sys
import textwrap
import plugins.core

# Parse and clean arguments
def parseArguments():
    argParser = argparse.ArgumentParser(description='PALADIN Pipeline Plugins')
    argParser.add_argument('-l', dest='list', action='store_true', help='List available plugins') 
    argParser.add_argument('plugins', metavar='Plugins', type=str, nargs=argparse.REMAINDER, help='Plugin with arguments (@@plugin arguments)')
    
    retArguments = argParser.parse_args()

    # Check for no arguments, or bad syntax in plugins list
    invalidArgs = False
    if not retArguments.list and not retArguments.plugins: invalidArgs = True
    if retArguments.plugins and not retArguments.plugins[0].startswith('@@'): invalidArgs = True

    if invalidArgs:
        argParser.print_help()
        sys.exit(1)
    
    # Replace quotes for strings with spaces
    for plugIdx in range(len(retArguments.plugins)):
        plugArg = retArguments.plugins[plugIdx]
        if ' ' in plugArg: plugArg = "\"{0}\"".format(plugArg)
        retArguments.plugins[plugIdx] = plugArg
    
    return retArguments

# Extract pipeline plugin/argument pairs from arguments
def getPipeline(passArguments):
    retPipeline = list()

    strArguments = ' '.join(passArguments)
    pluginGroups = strArguments.split('@@')
    if pluginGroups[0] != '': return retPipeline
             
    for pluginGroup in pluginGroups[1:]:
        pluginGroup = pluginGroup.rstrip()
        tokens = pluginGroup.split(' ')
        plugin = tokens[0]
        arguments = ' '.join(tokens[1:])
        retPipeline.append((plugin, arguments)) 

    return retPipeline

# List plugins currently available
def listPlugins(): 
    descStart = 21
    pluginList = plugins.core.pluginModules.values()
    sortedPlugins = sorted(pluginList, key=lambda x: x.name)

    print("The following plugins are available:")

    for plugin in sortedPlugins:
        versionStr = "{0}.{1}.{2}".format(plugin.versionMajor, plugin.versionMinor, plugin.versionRevision) 
        headerText ="{0} ({1}):".format(plugin.name, versionStr)
        rawText = "{0}\n{1}{2}".format(headerText, ' ' * (descStart - len(headerText)),  plugin.description)  
     
        wrapper = textwrap.TextWrapper(initial_indent='  ', subsequent_indent = ' ' * (descStart + 3), width=78)
        renderText = wrapper.wrap(rawText)

        for line in renderText:
            print(line)

        print()

    print("The following plugins are disabled (rename to enable): {0}".format(' '.join(plugins.modulesDisabled)))
    
# Parse arguments
arguments = parseArguments()

# Connect to plugins
plugins.core.connectPlugins()

# Handle non-pipeline actions
pipelinePresent = True
if arguments.list:
    pipelinePresent = False
    listPlugins()

if pipelinePresent:
    pipeline = getPipeline(arguments.plugins) 

    # Initialize plugins
    if plugins.core.initPlugins(set(list(zip(*pipeline))[0])): 
        # Execute pipeline
        plugins.core.execPipeline(pipeline)

    # Do a final flush to standard out
    plugins.core.execPipeline([('flush', '')])
