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

# Generate plots in PNG format from pipeline generated data

import os
import sys
import shlex
import argparse
import numpy
import matplotlib
matplotlib.use('Agg') # Allows matplotlib to work in non-X
import matplotlib.pyplot as plot
import plugins.core

plotGrid = None

# Plugin connection definition
def pluginConnect(passDefinition):
    passDefinition.name = 'plotting'
    passDefinition.description = 'Generate plots in PNG format from pipeline generated data'
    passDefinition.versionMajor = 1
    passDefinition.versionMinor = 0
    passDefinition.versionRevision = 0

    passDefinition.callbackInit = plottingInit
    passDefinition.callbackMain = plottingMain

# Plugin initialization   
def plottingInit():
    setupGrid(1, 1)

# Plugin main
def plottingMain(passArguments):
    # Parse arguments 
    argParser = argparse.ArgumentParser(description='PALADIN Pipeline Plugins: Plotting', prog='plotting')
    argParser.add_argument('-i', dest='input', type=str, help='Two column input file path')
    argParser.add_argument('-o', dest='output', type=str, help='Output PNG file path')
    argParser.add_argument('-l', dest='limit', type=int, help='Limit the number of rows shown on plot')
    argParser.add_argument('-g', dest='grid', type=int, nargs=2, metavar=('ROWS', 'COLUMNS'), help='Create new grid layout')
    argParser.add_argument('-c', dest='location', type=int, nargs=2, metavar=('ROW', 'COLUMN'), help='Set current grid location')
    argParser.add_argument('-L', dest='labels', type=str, nargs=3, metavar=('TITLE', 'X-AXIS', 'Y-AXIS'), help='Plot labels')
    argParser.add_argument('-t', dest='type', choices=('pie', 'bar'), type=str, default='pie', help='Plot type')
    argParser.add_argument('-p', dest='prepend', action='store_true', help='Prepend value to labels')
    argParser.add_argument('-s', dest='size', type=int, nargs=2, metavar=('WIDTH', 'HEIGHT'), help='Size of plot')
    argParser.add_argument('-C', dest='columns', type=int, nargs=2, metavar=('VALUES', 'LABELS'), default=[2, 3], help='Columns (Default: 2 3)')
    arguments = argParser.parse_known_args(shlex.split(passArguments))
   
    # Since plotting has no required arguments, check for completely empty
    if not passArguments: 
        argParser.print_help()
        sys.exit(1)

    # Size of chart
    if arguments[0].size: setSize(arguments[0].size[0], arguments[0].size[1])

    # Grid is being setup
    if arguments[0].grid: setupGrid(arguments[0].grid[0], arguments[0].grid[1])

    # Grid location is changing
    if arguments[0].location: changeLocation(arguments[0].location[0], arguments[0].location[1])

    # Plot title
    if arguments[0].labels: setLabels(arguments[0].labels[0], arguments[0].labels[1], arguments[0].labels[2])
 
    # Plot charts
    if arguments[0].input:  
        names, values = parseData(arguments[0].input, arguments[0].prepend, arguments[0].columns)
        if arguments[0].limit: names, values = limitData(names, values, arguments[0].limit)
        if arguments[0].type == 'pie': plotPie(names, values)
        if arguments[0].type == 'bar': plotBarH(names, values)

    # Save chart
    if arguments[0].output: savePlot(arguments[0].output)

# Expects a tab-delimited, 2 column text file with header row
def parseData(passFile, passPrepend, passColumns):
    retNames = list()
    retData = list()

    with open(passFile, 'r') as fileHandle:
        fileHandle.readline()

        # Get applicable column indices
        valueCol = passColumns[0] - 1
        nameCol = passColumns[1] - 1

        # Extract data from each line of file
        for line in fileHandle:
            line = line.rstrip()
            fields = line.split("\t")
            if len(fields) < 2: continue
            if passPrepend: 
                retNames.append("({0:.2f}) {1}".format(float(fields[valueCol]), fields[nameCol]))
            else:
                retNames.append(fields[nameCol])

            retData.append(fields[valueCol])

    return retNames, retData

# Limit data to the requested number of items, binning the rest into 'Other'
def limitData(passNames, passData, passLimit):
    retNames = list()
    retData = list()

    # Copy unfiltered data
    for idx in range(passLimit):
        if idx >= len(passData): break
        retNames.append(passNames[idx])
        retData.append(passData[idx])

    # Aggregate filtered data
    otherSize = 0.0
    for idx in range(passLimit, len(passData)):
        otherSize += float(passData[idx])

    if passLimit < len(passData):
        retNames.append('Other')
        retData.append(otherSize)

    return retNames, retData

# Setup new plotting grid
def setupGrid(passRows, passCols):
    plugins.plotting.plotGrid = plot.GridSpec(passRows, passCols)

# Change current grid location
def changeLocation(passRow, passCol):
    plot.subplot(plugins.plotting.plotGrid[passRow, passCol])

# Set total plot size
def setSize(passWidth, passHeight):
    plot.figure(figsize=(passWidth, passHeight), dpi=96)

# Set title, x-axis, and y-axis labels
def setLabels(passTitle='', passX = '', passY = ''):
    plot.title(passTitle, weight='bold')
    plot.xlabel(passX, weight='bold')
    plot.ylabel(passY, weight='bold')

# Save plot to PNG file
def savePlot(passFile):
    plot.tight_layout()
    plot.savefig(passFile, dpi=96)

# Plot a pie chart
def plotPie(passNames, passData):
    plugins.core.sendOutput('Generating pie chart...', 'stderr')
    patches, text = plot.pie(passData, shadow=True, startangle=90)
    plot.axis('equal')
    plot.legend(patches, passNames, loc='lower left')

# Plot a horizontal bar chart
def plotBarH(passNames, passData):
    plugins.core.sendOutput('Generating bar chart...', 'stderr')
    floatData = [float(x) for x in passData]

    xPos = numpy.arange(len(passData))
    plot.barh(xPos, floatData[::-1], align='center')
    plot.yticks(xPos, passNames[::-1])
