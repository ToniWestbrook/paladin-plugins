#! /usr/bin/env python3

"""
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
"""

# Generate plots in PNG format from pipeline generated data

import sys
import shlex
import argparse
import numpy
import matplotlib
matplotlib.use("Agg")  # Allows matplotlib to work in non-X
import matplotlib.pyplot as plot
import core.main

plot_grid = None


def plugin_connect(definition):
    definition.name = "plotting"
    definition.description = "Generate plots in PNG format from pipeline generated data"
    definition.version_major = 1
    definition.version_minor = 1
    definition.version_revision = 0

    definition.callback_args = plotting_args
    definition.callback_init = plotting_init
    definition.callback_main = plotting_main


def plotting_args(subargs):
    # Parse arguments
    arg_parser = argparse.ArgumentParser(description="PALADIN Pipeline Plugins: Plotting", prog="plotting")
    arg_parser.add_argument("-i", dest="input", type=str, help="Two column input file path")
    arg_parser.add_argument("-o", dest="output", type=str, help="Output PNG file path")
    arg_parser.add_argument("-l", dest="limit", type=int, help="Limit the number of rows shown on plot")
    arg_parser.add_argument("-g", dest="grid", type=int, nargs=2, metavar=("ROWS", "COLUMNS"), help="Create new grid layout")
    arg_parser.add_argument("-c", dest="location", type=int, nargs=2, metavar=("ROW", "COLUMN"), help="Set current grid location")
    arg_parser.add_argument("-L", dest="labels", type=str, nargs=3, metavar=("TITLE", "X-AXIS", "Y-AXIS"), help="Plot labels")
    arg_parser.add_argument("-t", dest="type", choices=("pie", "bar"), type=str, help="Plot type")
    arg_parser.add_argument("-p", dest="prepend", action="store_true", help="Prepend value to labels")
    arg_parser.add_argument("-s", dest="size", type=int, nargs=2, metavar=("WIDTH", "HEIGHT"), help="Size of plot")
    arg_parser.add_argument("-C", dest="columns", type=int, nargs=2, metavar=("VALUES", "LABELS"), help="Columns")

    args = arg_parser.parse_known_args(shlex.split(subargs))

    # Since plotting has no required arguments, check for completely empty
    if not subargs:
        arg_parser.print_help()
        sys.exit(1)

    return args


def plotting_init():
    setup_grid(1, 1)


def plotting_main(args):
    # Size of chart
    if args[0].size:
        set_size(args[0].size[0], args[0].size[1])

    # Grid is being setup
    if args[0].grid:
        setup_grid(args[0].grid[0], args[0].grid[1])

    # Grid location is changing
    if args[0].location:
        change_location(args[0].location[0], args[0].location[1])

    # Plot title
    if args[0].labels:
        set_labels(args[0].labels[0], args[0].labels[1], args[0].labels[2])

    # Plot charts
    if args[0].input:
        names, values = parse_data(args[0].input, args[0].prepend, args[0].columns)

        if args[0].limit:
            names, values = limit_data(names, values, args[0].limit)
        if args[0].type == "pie":
            plot_pie(names, values)
        if args[0].type == "bar":
            plot_bar_h(names, values)

    # Save chart
    if args[0].output: save_plot(args[0].output)


def parse_data(filename, prepend, columns):
    """ Expects a tab-delimited, 2 column text file with header row """
    ret_names = list()
    ret_data = list()

    with open(filename, "r") as handle:
        handle.readline()

        # Get applicable column indices
        value_col = columns[0] - 1
        name_col = columns[1] - 1

        # Extract data from each line of file
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 2:
                continue

            if prepend:
                ret_names.append("({0:.2f}) {1}".format(float(fields[value_col]), fields[name_col]))
            else:
                ret_names.append(fields[name_col])

            ret_data.append(fields[value_col])

    return ret_names, ret_data


def limit_data(names, data, limit):
    """ Limit data to the requested number of items, binning the rest into "Other" """
    ret_names = list()
    ret_data = list()

    # Copy unfiltered data
    for idx in range(limit):
        if idx >= len(data):
            break

        ret_names.append(names[idx])
        ret_data.append(data[idx])

    # Aggregate filtered data
    otherSize = 0.0
    for idx in range(limit, len(data)):
        otherSize += float(data[idx])

    if limit < len(data):
        ret_names.append("Other")
        ret_data.append(otherSize)

    return ret_names, ret_data


def setup_grid(rows, cols):
    """ Setup new plotting grid """
    global plot_grid
    plot_grid = plot.GridSpec(rows, cols)


def change_location(row, col):
    """ Change current grid location """
    plot.subplot(plot_grid[row, col])


def set_size(width, height):
    """ Set total plot size """
    plot.figure(figsize=(width, height), dpi=96)


def set_labels(title="", x_label="", y_label=""):
    """ Set title, x-axis, and y-axis labels """
    plot.title(title, weight="bold")
    plot.xlabel(x_label, weight="bold")
    plot.ylabel(y_label, weight="bold")


def save_plot(filename):
    """ Save plot to PNG file """
    plot.tight_layout()
    plot.savefig(filename, dpi=96)


def plot_pie(names, data):
    """ Plot a pie chart """
    core.main.send_output("Generating pie chart...", "stderr")
    patches, text = plot.pie(data, shadow=True, startangle=90)
    plot.axis("equal")
    plot.legend(patches, names, loc="lower left")


def plot_bar_h(names, data):
    """ Plot a horizontal bar chart """
    core.main.send_output("Generating bar chart...", "stderr")
    float_data = [float(x) for x in data]

    x_pos = numpy.arange(len(data))
    plot.barh(x_pos, float_data[::-1], align="center")
    plot.yticks(x_pos, names[::-1])
