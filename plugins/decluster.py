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

# Generate a reference containing each protein sequence used to generate the UniRef clusters detected in a PALADIN UniProt report

import argparse
import shlex
import sys
import core.main
from core.datastore import DataStore
from core.filestore import FileStore


module_definition = None
id_lookup = dict()
cluster_lookup = dict()
sequence_lookup = dict()


def plugin_connect(definition):
    definition.name = "decluster"
    definition.description = "Generate a reference containing each protein sequence used to generate the UniRef clusters detected in a PALADIN UniProt report"
    definition.version_major = 1
    definition.version_minor = 1
    definition.version_revision = 0
    definition.dependencies = ["crossref"]

    definition.callback_args = decluster_args
    definition.callback_init = decluster_init
    definition.callback_main = decluster_main


def decluster_args(subargs):
    # Parse arguments
    arg_parser = argparse.ArgumentParser(description="PALADIN Pipeline Plugins: Decluster", prog="decluster")
    arg_parser.add_argument("-i", dest="input", metavar="INPUT", type=str, required=True, help="PALADIN UniProt report")
    arg_parser.add_argument("-q", dest="quality", type=int, required=True, help="Minimum mapping quality filter")

    return arg_parser.parse_known_args(shlex.split(subargs))


def decluster_init():
    # Setup FileStore
    FileStore("decluster-db", "decluster-db", "decluster.db", None, FileStore.FTYPE_CACHE, FileStore.FOPT_NORMAL)
    FileStore("decluster-seqs", "decluster-swissprot", "decluster_uniprot_sprot.fasta.gz", "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz", FileStore.FTYPE_CACHE, FileStore.FOPT_GZIP_DECOMPRESS)
    FileStore("decluster-seqs", "decluster-trembl", "decluster_uniprot_trembl.fasta.gz", "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz", FileStore.FTYPE_CACHE, FileStore.FOPT_GZIP_DECOMPRESS)

    # Setup DataStore
    DataStore("decluster", FileStore.get_entry("decluster-db").path)
    DataStore.get_entry("decluster").create_table("indices", [("id", "text", "PRIMARY KEY"), ("file", "text", ""), ("pos", "integer", "")])
    DataStore.get_entry("decluster").define_query("index-lookup", "SELECT file, pos FROM indices WHERE id = ?")

    populate_database()

def decluster_main(args):
    # Get entries
    entries = core.main.PaladinEntry.get_entries(args[0].input, args[0].quality)
    cluster_ids = ["UniRef90_{0}".format(entries[entry].id) for entry in entries]
    core.main.send_output("Parsed {0} entries from UniProt report".format(len(cluster_ids)), "stderr")

    # Render sequences
    render_sequences(cluster_ids)


def populate_database():
    """ Populate sequence header indices """
    if not DataStore.get_entry("decluster").get_expired("indices", 30):
        return

    core.main.send_output("Populating UniProt sequences...", "stderr")

    # Start transaction and empty any existing data
    DataStore.get_entry("decluster").process_trans()
    DataStore.get_entry("decluster").delete_rows("indices")

    # Download each sequence file
    for entry in FileStore.get_group("decluster-seqs"):
        entry.prepare()

        with entry.get_handle("rt") as handle:
            acc = ""

            while True:
                line = handle.readline()
                if not line:
                    break

                if line.startswith(">"):
                    fields = line.rstrip().split()
                    acc = fields[0].split("|")[1]
                    DataStore.get_entry("decluster").insert_rows("indices", [(acc, entry.fid, handle.tell() - len(line))])

    # Finalize transaction and current table age
    DataStore.get_entry("decluster").process_trans()
    DataStore.get_entry("decluster").update_age("indices")


def get_sequence(acc):
    """ Lookup sequence for given acc """
    result = DataStore.get_entry("decluster").exec_query("index-lookup", [acc]).fetchone()
    if not result:
        core.main.send_output("Sequence not found for UniProt accession '{0}'".format(acc))
        sys.exit(1)

    # Seek the index position in the appropriate file
    handle = FileStore.get_entry(result[0], "decluster-seqs").get_handle()
    handle.seek(result[1])

    # Append sequence data until next header
    ret_seq = ""
    for line in handle:
        if line.startswith(">") and ret_seq:
            break

        ret_seq += line

    return ret_seq


def render_sequences(cluster_ids):
    """ Retrieve all members of requested UniRef90 clusters and render fasta data """
    # Prepare all sequence files for reading
    for entry in FileStore.get_group("decluster-seqs"):
        entry.get_handle("rt")

    for cluster_id in cluster_ids:
        for result in DataStore.get_entry("crossref").exec_query("uniprot_cross_acc", ("UniRef90", cluster_id)).fetchmany():
            core.main.send_output(get_sequence(result[0]), "stdout", "")
