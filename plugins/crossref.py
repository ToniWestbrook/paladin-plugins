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

# Provide database cross-references between IDs

import argparse
import shlex
import core.main
from core.datastore import DataStore
from core.filestore import FileStore


def plugin_connect(definition):
    definition.name = "crossref"
    definition.description = "Provide database cross-references between IDs"
    definition.version_major = 1
    definition.version_minor = 1
    definition.version_revision = 0

    definition.callback_args = crossref_args
    definition.callback_init = crossref_init
    definition.callback_main = crossref_main


def crossref_args(sub_args):
    arg_parser = argparse.ArgumentParser(description="PALADIN Pipeline Plugins: Cross Reference", prog="crossref")

    return arg_parser.parse_known_args(shlex.split(sub_args))


def crossref_init():
    # Setup FileStore
    FileStore("crossref-db", "crossref-db", "crossref.db", None, FileStore.FTYPE_CACHE, FileStore.FOPT_NORMAL)
    FileStore("crossref-uniprot", "crossref-uniprot", "idmapping.dat.gz", "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz", FileStore.FTYPE_TEMP, FileStore.FOPT_GZIP)

    # Setup DataStore
    DataStore("crossref", FileStore.get_entry("crossref-db").get_path())
    DataStore.get_entry("crossref").create_table("uniprot", [("acc", "text", ""), ("db", "text", ""), ("cross", "text", "")])
    DataStore.get_entry("crossref").define_query("uniprot_acc_cross", "SELECT cross FROM uniprot WHERE acc = ? AND db = ?")
    DataStore.get_entry("crossref").define_query("uniprot_acc_all", "SELECT db, cross FROM uniprot WHERE acc = ?")
    DataStore.get_entry("crossref").define_query("uniprot_cross_acc", "SELECT acc FROM uniprot WHERE db = ? AND cross = ?")
    DataStore.get_entry("crossref").define_query("uniprot_cross_cross", "SELECT t2.cross FROM uniprot AS t1 JOIN uniprot AS t2 ON acc WHERE t1.db = ? AND t1.cross = ? AND t2.db = ?")
    DataStore.get_entry("crossref").define_index("uniprot_acc", "uniprot", ["acc"], False)
    DataStore.get_entry("crossref").define_index("uniprot_acc_db", "uniprot", ["acc", "db"], False)
    DataStore.get_entry("crossref").define_index("uniprot_db_cross", "uniprot", ["db", "cross"], False)

    # Populate database
    populate_database()


def crossref_main(args):
    """ Crossref provides API only and no direct functionality to user """
    pass


def populate_database():
    """ Generate cross-reference database """
    if not DataStore.get_entry("crossref").get_expired("uniprot", 30):
        return

    core.main.send_output("Populating UniProt database cross-references...", "stderr")

    # Download tab delimited data
    entry = FileStore.get_entry("crossref-uniprot")
    entry.prepare()

    # Start transaction and empty any existing data
    DataStore.get_entry("crossref").drop_index("uniprot_acc")
    DataStore.get_entry("crossref").drop_index("uniprot_acc_db")
    DataStore.get_entry("crossref").drop_index("uniprot_db_cross")
    DataStore.get_entry("crossref").process_trans()
    DataStore.get_entry("crossref").delete_rows("uniprot")

    # Iterate through downloaded table and add rows
    with entry.get_handle("rt") as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                continue

            # Add to database
            DataStore.get_entry("crossref").insert_rows("uniprot", [fields])

    # Finalize transaction and current table age
    DataStore.get_entry("crossref").process_trans()
    DataStore.get_entry("crossref").create_index("uniprot_acc")
    DataStore.get_entry("crossref").create_index("uniprot_acc_db")
    DataStore.get_entry("crossref").create_index("uniprot_db_cross")
    DataStore.get_entry("crossref").update_age("uniprot")
