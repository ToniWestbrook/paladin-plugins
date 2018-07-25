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

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OF
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

# Combine multiple PALADIN outputs (both SAM and UniProt) into a single output

import os
import re
import operator
import argparse
import shlex
import core.main

def plugin_connect(definition):
    """ Plugin connection definition """
    definition.name = "aggregation"
    definition.description = "Combine multiple PALADIN outputs (both SAM and UniProt) into a single output"
    definition.version_major = 1
    definition.version_minor = 1
    definition.version_revision = 0

    definition.callback_args = aggregation_args
    # definition.callback_init = aggregation_init
    definition.callback_main = aggregation_main


def aggregation_args(sub_args):
    """ Parse arguments """
    arg_parser = argparse.ArgumentParser(description="PALADIN Pipeline Plugins: Aggregation", prog="aggregation")
    arg_parser.add_argument("-r", dest="root_path", metavar="ROOT", type=str, required=True, help="Root path to search")
    arg_parser.add_argument("-s", dest="pattern", metavar="PATTERN", type=str, required=True, help="SAM search pattern (infers TSV)")
    arg_parser.add_argument("-o", dest="output", metavar="OUTPUT", type=str, required=True, help="Output path")
    return arg_parser.parse_known_args(shlex.split(sub_args))


def aggregation_main(args):
    # Build base paths
    base_paths = build_paths(args[0].root_path, args[0].pattern)

    # Aggregate data
    process_data(base_paths, args[0].output)


def process_data(base_paths, output_base):
    """ API - Aggregate both SAM and TSV data """
    # Aggregate SAM data and write
    sam_data = aggregate_sam(base_paths)
    with open("{0}.sam".format(output_base), "w") as handle:
        handle.write(sam_data)

    # Aggregate TSV data and write
    tsv_header, tsv_entries = aggregate_tsv(base_paths)
    tsv_entries = finalize_tsv(tsv_entries)

    with open("{0}_uniprot.tsv".format(output_base), "w") as handle:
        handle.write(tsv_header + "\n")
        for entry in tsv_entries:
            handle.write("\t".join([str(val) for val in entry]) + "\n")


def build_paths(path, pattern):
    """ Recurse and build path information using SAM pattern """
    ret_entries = list()

    # Iterate through all subdirectories/files, match pattern
    for root, _, files in os.walk(path):
        for filename in files:
            if not re.match(pattern, filename):
                continue

            # Remove SAM extension, get base file path
            base_file = os.path.join(root, filename)[:-4]
            ret_entries.append(base_file)

    return ret_entries


def aggregate_sam(paths):
    """ Combine SAM files """
    entries = list()
    header_file = True

    core.main.send_output("Gathering SAM data...", "stderr")

    for path in paths:
        path = "{0}.sam".format(path)
        core.main.send_output("Processing {0}...".format(path), "stderr")

        with open(path, "r") as handle:
            for line in handle:
                if header_file:
                    # First file will include headers
                    entries.append(line)
                else:
                    # Subsequence files will not
                    if line.startswith("@"):
                        continue
                    entries.append(line)

            header_file = False

    return "".join(entries)


def aggregate_tsv(paths):
    """ Group TSV files """
    ret_entries = dict()
    ret_header = ""

    core.main.send_output("Gathering UniProt data...", "stderr")

    for path in paths:
        path = "{0}_uniprot.tsv".format(path)
        core.main.send_output("Processing {0}...".format(path), "stderr")

        # Open individual report, save header
        with open(path, "r") as handle:
            ret_header = handle.readline().rstrip()

            for line in handle:
                fields = line.rstrip().split("\t")

                # Check for first instance of entry
                existing_entry = (0, 0.0, 0.0, 0)
                if fields[4] in ret_entries:
                    existing_entry = ret_entries[fields[4]]

                # Update and store
                line_count = int(fields[0])
                line_qual_max = int(fields[3])

                new_count = existing_entry[0] + line_count
                new_abundance = existing_entry[1] + (line_count * float(fields[1]))
                new_qual_avg = existing_entry[2] + (line_count * float(fields[2]))
                new_qual_max = existing_entry[3]
                if line_qual_max > new_qual_max:
                    new_qual_max = line_qual_max

                ret_entries[fields[4]] = (new_count, new_abundance, new_qual_avg, new_qual_max) + tuple(fields[4:])

    return ret_header, ret_entries


def finalize_tsv(entries):
    """ Recalculate and sort TSV data """
    ret_entries = list()

    # Summate counts
    total_count = 0
    for entry in entries:
        total_count += entries[entry][0]

    # Sort
    entries_sorted = sorted(entries.items(), key=operator.itemgetter(1), reverse=True)

    # Calculate means and store new entry
    for entry in entries_sorted:
        abundance = entry[1][0] / total_count * 100
        quality_avg = entry[1][2] / entry[1][0]
        new_entry = (entry[1][0], abundance, quality_avg) + entry[1][3:]

        ret_entries.append(new_entry)

    return ret_entries
