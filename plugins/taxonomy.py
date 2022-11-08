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

# Perform taxonomic categorization on a PALADIN UniProt report

import re
import operator
import argparse
import shlex
import core.main
from core.main import SamEntry
from core.main import PaladinEntry
from core.datastore import DataStore
from core.filestore import FileStore


def plugin_connect(definition):
    definition.name = "taxonomy"
    definition.description = "Perform taxonomic grouping and abundance reporting"
    definition.version_major = 1
    definition.version_minor = 1
    definition.version_revision = 3

    definition.callback_args = taxonomy_args
    definition.callback_init = taxonomy_init
    definition.callback_main = taxonomy_main


def taxonomy_args(sub_args):
    arg_parser = argparse.ArgumentParser(description="PALADIN Pipeline Plugins: Taxonomy", prog="taxonomy")
    arg_parser.add_argument("-i", dest="input", type=str, required=True, help="PALADIN UniProt report")
    arg_parser.add_argument("-q", dest="quality", type=int, required=True, help="Minimum mapping quality filter")
    arg_parser.add_argument("-c", dest="custom", type=str, default="", help="Species-parsing regex pattern for non-Uniprot entries")
    arg_parser.add_argument("-t", dest="type", type=str, default="species", choices=["children", "species"], help="Type of report")
    arg_parser.add_argument("-f", dest="filter", type=str, choices=["unknown", "custom", "group"], nargs="*", help="Filter non-standard entries")
    arg_parser.add_argument("-s", dest="sam", type=str, help="SAM file for reporting reads contributing to each taxonomic rank")
    levelGroup = arg_parser.add_mutually_exclusive_group(required=True)
    levelGroup.add_argument("-l", dest="level", type=int, help="Hierarchy level")
    levelGroup.add_argument("-r", dest="rank", type=str, help="Regex pattern for named rank")

    return arg_parser.parse_known_args(shlex.split(sub_args))


def taxonomy_init():
    # Setup FileStore
    FileStore("taxonomy-db", "taxonomy-db", "taxonomy.db", None, FileStore.FTYPE_CACHE, FileStore.FOPT_NORMAL)
    FileStore("taxonomy-lineage", "taxonomy-lineage", "taxonomy-lineage.dat", "https://rest.uniprot.org/taxonomy/stream?fields=id%2Cmnemonic%2Clineage&format=tsv&query=%28%2A%29", FileStore.FTYPE_TEMP, FileStore.FOPT_NORMAL)

    # Setup DataStore
    DataStore("taxonomy", FileStore.get_entry("taxonomy-db").path)
    DataStore.get_entry("taxonomy").create_table("lineage", [("mnemonic", "text", "PRIMARY KEY"), ("lineage", "text", "")])
    DataStore.get_entry("taxonomy").define_query("lineage-lookup", "SELECT lineage FROM lineage WHERE mnemonic = ?")

    # Populate database
    populate_database()


def taxonomy_main(args):
    if args[0].filter is None:
        args[0].filter = []

    # Get entries
    full_entries = PaladinEntry.get_entries(args[0].input, args[0].quality, args[0].custom)
    species_lookup = get_species_lookup(full_entries)
    taxa_entries, taxa_count = aggregate_taxa(full_entries, args[0].filter)
    lineage_tree = treeify_lineage(taxa_entries, taxa_count)

    # Check report combinations
    if args[0].type == "children":
        if args[0].level is not None:
            # Children of a flattened level
            render_entries, render_count = flatten_tree(lineage_tree, args[0].level)
            header = "Count\tAbundance\tRank {0}".format(args[0].level)
        else:
            # Children of a rank subtree
            rank_tree = find_rank_subtree(lineage_tree, args[0].rank)
            render_entries, render_count = get_tree_children(rank_tree)
            header = "Count\tAbundance\tRank"
    else:
        if args[0].level is not None:
            # All species for entire tree
            render_entries, render_count = filter_taxa(taxa_entries, ".*")
            render_entries = {species_lookup[taxon]: render_entries[taxon] for taxon in render_entries}
            header = "Count\tAbundance\tSpecies"
        else:
            # All species of a rank subtree
            render_entries, render_count = filter_taxa(taxa_entries, args[0].rank)
            render_entries = {species_lookup[taxon]: render_entries[taxon] for taxon in render_entries}
            header = "Count\tAbundance\tSpecies"

    if not args[0].sam:
        # Render report (if non-SAM mode)
        render_abundance(render_entries, header, render_count)
    else:
        # In SAM mode, generate reads report if requested
        sam_data = group_sam(render_entries, args[0].sam, species_lookup)
        render_sam(sam_data)


def populate_database():
    """ Generate (if necessary) and get lineage lookup """
    if not DataStore.get_entry("taxonomy").get_expired("lineage", 30):
        return

    core.main.send_output("Populating taxonomic lineage data...", "stderr")

    # Download tab delimited data
    entry = FileStore.get_entry("taxonomy-lineage")
    entry.prepare()

    # Start transaction and empty any existing data
    DataStore.get_entry("taxonomy").process_trans()
    DataStore.get_entry("taxonomy").delete_rows("lineage")

    # Iterate through downloaded table and add rows
    with entry.get_handle("r") as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 3 or fields[1] == "":
                continue

            # Add to database
            DataStore.get_entry("taxonomy").insert_rows("lineage", [(fields[1], fields[2])])

    # Finalize transaction and current table age
    DataStore.get_entry("taxonomy").process_trans()
    DataStore.get_entry("taxonomy").update_age("lineage")


def aggregate_taxa(entries, filters):
    """ Group, filter, and count taxa """
    ret_taxa = dict()
    ret_total = 0

    # For each entry, filter, get species ID, and summate count
    for entry in entries:
        # Check for removed filters
        if "unknown" in filters and entries[entry].type == PaladinEntry.TYPE_UNKNOWN:
            continue
        if "custom" in filters and entries[entry].type == PaladinEntry.TYPE_CUSTOM:
            continue
        if "group" in filters and entries[entry].type == PaladinEntry.TYPE_UNIPROT_GROUP:
            continue

        key = (entries[entry].species_id, entries[entry].species_full)
        if key not in ret_taxa:
            ret_taxa[key] = 0

        count = entries[entry].count
        ret_taxa[key] += count
        ret_total += count

    return ret_taxa, ret_total


def filter_taxa(taxa, pattern):
    """ Filter taxa having a specific taxonomic rank (as regex) """
    ret_taxa = dict()
    ret_count = 0

    for taxon in taxa:
        result = DataStore.get_entry("taxonomy").exec_query("lineage-lookup", [taxon[0]]).fetchone()
        lineage = result[0] if result else "Unknown"

        if re.search(pattern, lineage):
            ret_taxa[taxon] = taxa[taxon]
            ret_count += taxa[taxon]

    return ret_taxa, ret_count


def get_species_lookup(entries):
    """ Create a species_id to species full lookup """
    ret_lookup = dict()

    for entry in entries:
        # For virtual group entries, save as full species name
        key = (entries[entry].species_id, entries[entry].species_full)

        if key not in ret_lookup:
            ret_lookup[key] = entries[entry].species_full

    return ret_lookup


def treeify_lineage(taxa_entries, taxa_count):
    """ Create tree structure out of taxa/lineage data """
    ret_tree = (dict(), taxa_count)

    for taxon in taxa_entries:
        result = DataStore.get_entry("taxonomy").exec_query("lineage-lookup", [taxon[0]]).fetchone()
        if not result:
            continue

        raw_lineage = result[0]
        ranks = raw_lineage.split(";")
        parent = ret_tree[0]

        for rank in ranks:
            rank = rank.strip()

            # Start new dictionary, initialize size cache
            if rank not in parent:
                parent[rank] = (dict(), 0)

            parent[rank] = (parent[rank][0], parent[rank][1] + taxa_entries[taxon])
            parent = parent[rank][0]

    return ret_tree


def find_rank_subtree(lineage, rank):
    """ Find a specific rank instance subtree within a lineage tree """
    # Breadth-wise search
    for rank in lineage[0]:
        if rank == rank:
            return lineage[0][rank]

    # Recurse if necessary
    for rank in lineage[0]:
        ret_tree = find_rank_subtree(lineage[0][rank], rank)
        if ret_tree:
            return ret_tree


def get_tree_children(lineage):
    """ Get immediate children for a specific lineage rank branch """
    ret_children = dict()
    ret_count = 0

    for rank in lineage[0]:
        count = lineage[0][rank][1]
        ret_children[rank] = count
        ret_count += count

    return ret_children, ret_count


def flatten_tree(lineage, level):
    """ Flatten a lineage tree to a specific level """
    ret_flat_entries = dict()
    ret_flat_count = 0

    # Negative level indicates last level
    if level < 0:
        for entry in lineage[0]:
            if not lineage[0][entry][0]:
                ret_flat_entries[entry] = lineage[0][entry][1]
                ret_flat_count += lineage[0][entry][1]

    if level == 0:
        # At level 0, extract lineage
        ret_flat_entries, ret_flat_count = get_tree_children(lineage)
    else:
        # At higher levels, recurse and append
        for entry in lineage[0]:
            subtree_entries, subtree_count = flatten_tree(lineage[0][entry], level - 1)
            ret_flat_entries.update(subtree_entries)
            ret_flat_count += subtree_count

    return ret_flat_entries, ret_flat_count


def render_abundance(data, header, total):
    """ Render standard dictionary abundance report """
    # Sort entries
    sorted_data = sorted(data.items(), key=operator.itemgetter(1), reverse=True)

    # Render data
    core.main.send_output(header)

    for row in sorted_data:
        abundance = row[1] / total * 100
        output_line = "{0}\t{1}\t{2}".format(row[1], abundance, row[0])
        core.main.send_output(output_line)


def group_sam(data, sam_file, species):
    """ Group SAM reads per taxonomic rank """
    retData = dict()

    # Simplify species lookup keys
    species_lookup = {item[0]: item[1] for item in species}

    # Read SAM entries
    entries = SamEntry.get_entries(sam_file, 0)

    for entry in entries:
        reference = entries[entry].reference

        # Currently only handles UniProt
        if "_" not in reference:
            continue

        mnemonic = reference.split("_")[1]

        # Skip entries that are too ambiguous (e.g. 9ZZZZ) or recently moved by UniProt
        if mnemonic not in species_lookup:
            continue

        result = DataStore.get_entry("taxonomy").exec_query("lineage-lookup", [mnemonic]).fetchone()
        if not result:
            continue

        # Combine lineage and species for lookup into taxonomy results
        lineage = [item.strip() for item in result[0].split(";")]
        lineage.append(species_lookup[mnemonic])

        for rank in lineage:
            if rank in data:
                query = entries[entry].query
                if query not in retData:
                    retData[query] = list()
                retData[query].append(rank)

    return retData


def render_sam(data):
    """ Render SAM reads per taxonomic rank """
    core.main.send_output("Read\tTaxonomy")

    for read in data:
        for rank in data[read]:
            output_line = "{0}\t{1}".format(read, rank)
            core.main.send_output(output_line)
