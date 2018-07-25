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

# Perform gene ontology term analysis on a PALADIN UniProt report

import argparse
import shlex
import operator
import core.main


def plugin_connect(definition):
    definition.name = "go"
    definition.description = "Perform gene ontology term grouping and abundance reporting"
    definition.version_major = 1
    definition.version_minor = 1
    definition.version_revision = 0

    definition.callback_args = go_args
    # definition.callback_init = go_init
    definition.callback_main = go_main


def go_args(subargs):
    # Parse arguments
    arg_parser = argparse.ArgumentParser(description="PALADIN Pipeline Plugins: GO", prog="go")
    arg_parser.add_argument("-i", dest="input", type=str, required=True, help="PALADIN UniProt report")
    arg_parser.add_argument("-q", dest="quality", type=int, required=True, help="Minimum mapping quality filter")

    return arg_parser.parse_known_args(shlex.split(subargs))


def go_main(args):
    # Get entries
    entries = core.main.PaladinEntry.get_entries(args[0].input, args[0].quality)

    # Aggregate for ontology
    ontology_data = aggregate_ontology(entries)

    # Render
    render_abundance(ontology_data)


def aggregate_ontology(entries):
    """ Group and count GO terms """
    ret_ontology = dict()

    for entry in entries:
        ontology = entries[entry].ontology
        if not ontology or ontology[0].rstrip() == "":
            continue

        for term in ontology:
            if term not in ret_ontology:
                ret_ontology[term] = 0

            ret_ontology[term] += entries[entry].count

    return ret_ontology


def render_abundance(data):
    """ Render the abundance of the GO terms """
    sorted_data = sorted(data.items(), key=operator.itemgetter(1), reverse=True)

    for row in sorted_data:
        output_line = "{0}\t{1}".format(row[0], row[1])
        core.main.send_output(output_line)
