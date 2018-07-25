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

# Distribute PALADIN execution across a cluster

import os
import argparse
import shlex
import subprocess
import gzip


def plugin_connect(definition):
    definition.name = "hpc"
    definition.description = "Distribute PALADIN execution across a cluster"
    definition.version_major = 1
    definition.version_minor = 1
    definition.version_revision = 0
    definition.dependencies = ["aggregation"]

    definition.callback_args = hpc_args
    # definition.callback_init = hpc_init
    definition.callback_main = hpc_main


def hpc_args(subargs):
    # Parse arguments
    arg_parser = argparse.ArgumentParser(description="PALADIN Pipeline Plugins: HPC", prog="hpc")
    arg_parser.add_argument("reference", metavar="REFERENCE", type=str, help="Reference database")
    arg_parser.add_argument("input", metavar="INPUT", type=str, help="Input reads")
    arg_parser.add_argument("output", metavar="OUTPUT", type=str, help="Output name")
    arg_parser.add_argument("options", metavar="OPTIONS", type=str, nargs=argparse.REMAINDER, help="PALADIN options")

    return arg_parser.parse_known_args(shlex.split(subargs))


def hpc_main(args):
    # Obtain MPI process info (import MPI here to avoid crashing on non-MPI systems)
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    process_count = comm.Get_size()
    process_idx = comm.Get_rank()

    # Process 0 is responsible for splitting reads
    if process_idx == 0:
        split_reads(args[0].input, process_count)

    # Sync to ensure all batches of reads ready for alignment
    comm.Barrier()

    # Execute PALADIN alignment
    process_reads = "{0}-{1}".format(args[0].input, process_idx)
    process_output = "{0}-{1}".format(args[0].output, process_idx)
    execute(args[0].reference, process_reads, process_output, " ".join(args[0].options))

    # Sync to ensure all alignments complete
    comm.Barrier()

    # Aggregate result
    if process_idx == 0:
        input_paths = ["{0}-{1}".format(args[0].output, idx) for idx in range(process_count)]
        plugins.aggregation.process_data(input_paths, args[0].output)

        # Remove process specific files
        for idx in range(process_count):
            os.remove("{0}-{1}".format(args[0].input, idx))
            os.remove("{0}-{1}.sam".format(args[0].output, idx))
            os.remove("{0}-{1}_uniprot.tsv".format(args[0].output, idx))


def execute(reference, input_name, output_name, options):
    """ Execute PALADIN in the current process """
    command = "paladin align {0} {1} -o {2} {3}".format(reference, input_name, output_name, options)
    output = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True)

    with open("{0}.log".format(output_name), "wb") as handle:
        handle.write(output)


def split_reads(reads, count):
    """ Split reads across the requested node count """
    # Treat as gzip file with gz extension
    if reads.endswith("gz"):
        input_handle = gzip.open(reads, "rb")
    else:
        input_handle = open(reads, "rb")

    # Open reads and outputs
    with input_handle:
        output_handles = list()
        for out_idx in range(count):
            output_handles.append(open("{0}-{1}".format(reads, out_idx), "w"))

        mode = -1
        out_idx = -1
        parse_idx = -1
        for line in input_handle:
            line = line.decode("utf-8")
            if line.rstrip() == "":
                continue

            # _initially detect input type (0 = fasta, 1 = fastq)
            if mode == -1:
                mode = (1, 0)[line.startswith(">")]

            if mode == 0:
                # FastA mode
                if line.startswith(">"):
                    out_idx = (out_idx + 1) % count

                output_handles[out_idx].write(line)
            else:
                # FastQ mode
                parse_idx = (parse_idx + 1) % 4
                if parse_idx == 0:
                    out_idx = (out_idx + 1) % count

                output_handles[out_idx].write(line)

        # Close output handles
        for handle in output_handles:
            handle.close()
