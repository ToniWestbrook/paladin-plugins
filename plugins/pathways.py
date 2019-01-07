#! /usr/bin/env python3

"""
The MIT License

Copyright (c) 2019 by Anthony Westbrook, University of New Hampshire <anthony.westbrook@unh.edu>

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

# Analyze metabolic pathway participation

import argparse
import re
import shlex
from collections import Counter
import requests
import core.main
from core.filestore import FileStore
from core.datastore import DataStore

KEGG_TEMPLATE_PATHWAY = "http://rest.kegg.jp/get/{kegg_id}"


def plugin_connect(definition):
    definition.name = "pathways"
    definition.description = "Analyze metabolic pathway participation"
    definition.version_major = 1
    definition.version_minor = 0
    definition.version_revision = 0
    definition.dependencies = ["taxonomy"]

    definition.callback_args = pathways_args
    definition.callback_init = pathways_init
    definition.callback_main = pathways_main


def pathways_args(subargs):
    arg_parser = argparse.ArgumentParser(description="PALADIN Pipeline Plugins: Pathways", prog="pathways")
    arg_parser.add_argument("-i", dest="input", type=str, required=True, nargs="+", help="PALADIN UniProt report(s)")
    arg_parser.add_argument("-q", dest="quality", type=int, required=True, help="Minimum mapping quality filter")
    arg_parser.add_argument("-p", dest="path_id", type=str, help="KEGG pathway ID")
    arg_parser.add_argument("-a", dest="abstract", type=int, default=0, help="Level of EC heirarchy abstraction")
    arg_parser.add_argument("-g", dest="grouping", type=str, default="children", choices=["children", "species"], help="Taxonomy grouping")
    arg_parser.add_argument("-r", dest="rounding", type=int, default=4, help="Rounding precision")
    group_group = arg_parser.add_mutually_exclusive_group()
    group_group.add_argument("-l", dest="level", type=int, help="Grouping: taxonomic level")
    group_group.add_argument("-s", dest="search", type=str, help="Grouping: regex search pattern")

    return arg_parser.parse_known_args(shlex.split(subargs))


def pathways_init():
    # Setup FileStore
    FileStore("pathways-db", "pathways-db", "pathways.db", None, FileStore.FTYPE_CACHE, FileStore.FOPT_NORMAL)

    # Setup DataStore
    DataStore("pathways", FileStore.get_entry("pathways-db").path)
    DataStore.get_entry("pathways").create_table("enzyme", [("ec", "text", "PRIMARY KEY"), ("pathway", "text", "")])
    DataStore.get_entry("pathways").create_table("pathway", [("pathway", "text", "PRIMARY KEY"), ("info", "text", "")])
    DataStore.get_entry("pathways").define_query("enzyme-lookup", "SELECT pathway FROM enzyme WHERE ec = ?")
    DataStore.get_entry("pathways").define_query("pathway-lookup", "SELECT info FROM pathway WHERE pathway = ?")

    # Check for expired database
    if DataStore.get_entry("pathways").get_expired("enzyme", 30):
        DataStore.get_entry("pathways").delete_rows("enzyme")
        DataStore.get_entry("pathways").delete_rows("pathway")
        DataStore.get_entry("pathways").update_age("enzyme")


def pathways_main(args):
    tsv_trees = list()

    if not args[0].path_id:
        core.main.send_output("Pathway unspecified, retrieving all detected pathways...", "stderr")

    # Iterate processing for each PALADIN TSV
    for tsv_input in args[0].input:
        # Process PALADIN TSV
        entries_full = core.main.PaladinEntry.get_entries(tsv_input, args[0].quality)
        entries_ec = populate_ec(entries_full)

        # Assign entries across lineages and trim/flatten
        full_tree = build_tree(entries_ec)

        if args[0].search:
            # Regex search for rank
            trim_tree = search_tree(full_tree, args[0].search)
        else:
            # Flatten tree at level
            trim_tree = flatten_tree(full_tree, args[0].level)

        # Flatten to speices if requested
        if args[0].grouping == "species":
            trim_tree = flatten_tree(trim_tree, -1)

        tsv_trees.append(trim_tree)

    # Ensure all trees have same groupings
    for tree_src in tsv_trees:
        for group in tree_src["children"]:
            for tree_dst in tsv_trees:
                tree_dst["children"].setdefault(group, {"enzymes": Counter()})

    # Retrieve pathways
    pathways_info = retrieve_pathways(args[0].path_id, entries_ec)

    # Apply indexing and abstraction to all enzymes
    for pathway_info in pathways_info.values():
        process_enzymes(pathway_info, args[0].abstract)

    # Render participation
    if len(pathways_info) == 1:
        render_single(next(iter(pathways_info.values())), tsv_trees, args[0].rounding)
    else:
        render_all(pathways_info, tsv_trees, args[0].rounding)


def retrieve_pathways(path_id, entries_ec):
    """ Retrieve all pathways for all present enzymes """
    pathways_info = dict()

    if path_id:
        # Lookup specified pathway
        pathways_info[path_id] = cache_pathway(path_id)
    else:
        # Retrieve all pathways for all present enzymes
        for species in entries_ec.values():
            for entry_ec in species:
                # Only process enzymes with associated pathways
                pathways_ec = cache_enzyme(entry_ec)
                if "PATHWAY" not in pathways_ec:
                    continue

                # Add each pathway
                for pathway in pathways_ec["PATHWAY"]:
                    pathway = pathway.split()[0]

                    if pathway not in pathways_info:
                        # Only process pathways with enzymes
                        pathway_info = cache_pathway(pathway)
                        if "ENZYME" not in pathway_info:
                            continue

                        pathways_info[pathway] = cache_pathway(pathway)

    return pathways_info


def cache_enzyme(enzyme_id):
    """ Lookup and/or cache enzyme """
    result = DataStore.get_entry("pathways").exec_query("enzyme-lookup", [enzyme_id]).fetchone()

    if result:
        # Cache hit
        return parse_kegg(result[0])

    # Cache miss
    enzyme_raw = retrieve_kegg(enzyme_id)
    DataStore.get_entry("pathways").insert_rows("enzyme", [(enzyme_id, enzyme_raw)])
    return parse_kegg(enzyme_raw)


def cache_pathway(path_id):
    """ Lookup and/or cache pathway """
    result = DataStore.get_entry("pathways").exec_query("pathway-lookup", [path_id]).fetchone()

    if result:
        # Cache hit
        return parse_kegg(result[0])

    # Cache miss
    pathway_raw = retrieve_kegg(path_id)
    DataStore.get_entry("pathways").insert_rows("pathway", [(path_id, pathway_raw)])
    return parse_kegg(pathway_raw)


def retrieve_kegg(kegg_id):
    """ Download information given a KEGG ID (pathway or enzyme) """
    # Download raw pathway information
    response = requests.post(KEGG_TEMPLATE_PATHWAY.format(kegg_id=kegg_id))
    kegg_raw = response.text
    response.close()

    return kegg_raw


def parse_kegg(kegg_raw):
    """ Parse raw KEGG output into dictionary entries """
    kegg_info = dict()

    # Check for empty entries
    if kegg_raw.strip() == "":
        return kegg_info

    # Extract fields from result
    for line in kegg_raw.split("\n"):
        fields = re.split(" +", line)
        if fields[0] != "":
            # Lines that don't start with space are start of fields
            field = fields[0]

        kegg_info.setdefault(field, list())
        kegg_info[field].append(" ".join(fields[1:]))

    return kegg_info


def process_enzymes(pathway_info, abstract):
    """ Index and abstract enzymes, if present """
    if "ENZYME" not in pathway_info:
        return

    # Index
    pathway_info["ENZYME"] = set(pathway_info["ENZYME"])

    # Abstract
    parent_enzymes = set()
    abstract = min(abstract, 4)

    for level in range(abstract):
        for enzyme in pathway_info["ENZYME"]:
            enzyme_groups = enzyme.split(".")

            # Abstract for the current level (count dashes)
            if enzyme_groups.count("-") == level:
                enzyme_groups[-(level + 1)] = "-"
                parent_enzymes.add(".".join(enzyme_groups))

        pathway_info["ENZYME"] = pathway_info["ENZYME"].union(parent_enzymes)


def populate_ec(entries):
    """ Extract EC codes per species across PALADIN entries """
    ret_ec = dict()

    for entry in entries.values():
        match = re.search("\\(EC ([0-9-]+\\.[0-9-]+\\.[0-9-]+\\.[0-9-]+)\\)", entry.protein)
        if match:
            entry_id = (entry.species_id, entry.species_full)
            ec_id = match.groups(0)[0]
            ret_ec.setdefault(entry_id, dict())
            ret_ec[entry_id].setdefault(ec_id, 0)
            ret_ec[entry_id][ec_id] += entry.count

    return ret_ec


def build_tree(entries):
    """ Create tree structure out of taxon lineage data """
    ret_tree = {"name": "all", "children": dict(), "enzymes": Counter()}

    for entry in entries:
        result = DataStore.get_entry("taxonomy").exec_query("lineage-lookup", [entry[0]]).fetchone()
        if not result:
            continue

        raw_lineage = result[0]

        # If a UniProt species-level species ID (non 9CCCC), append species to lineage
        if not entry[0].startswith("9"):
            raw_lineage += ";{0}".format(entry[1])

        ranks = raw_lineage.split(";")
        parent = ret_tree["children"]

        for rank in ranks:
            rank = "Unknown" if rank == "" else rank.strip()

            # Initialize dictionary for this rank, add child values
            parent.setdefault(rank, {"name": rank, "children": dict(), "enzymes": Counter()})
            parent[rank]["enzymes"] += Counter(entries[entry])
            parent = parent[rank]["children"]

        # Update root level
        ret_tree["enzymes"] += Counter(entries[entry])

    return ret_tree


def flatten_tree(tree, level):
    """ Flatten a lineage tree to a specific level """
    ret_tree = {"name": "flat", "children": {}, "enzymes": {}}

    # Negative level indicates species mode, treat terminal node as level 0
    if level < 0 and len(tree["children"]) == 0:
        level = 0

    if level == 0:
        # At level 0, return current tree
        ret_tree["children"][tree["name"]] = tree
        return ret_tree
    else:
        # At higher levels, recurse and merge trees
        for entry in tree["children"].values():
            temp_tree = flatten_tree(entry, level - 1)
            ret_tree["children"].update(temp_tree["children"])

    return ret_tree


def search_tree(tree, pattern):
    """ Search for ranks matching regex pattern """
    if re.search(pattern, tree["name"]):
        return tree
    else:
        for entry in tree["children"].values():
            result = search_tree(entry, pattern)
            if "name" in result:
                return result

    return {"children": {}}


def render_single(pathway_info, trees, rounding):
    """ Render the percentage of participation per entry for a single pathway """
    # Use first report for group list
    groups = trees[0]["children"].keys()

    # Render column headers
    row_headers = ["EC"]
    for group in groups:
        for tsv_tree in enumerate(trees):
            row_headers.append("{0}:{1}".format(tsv_tree[0] + 1, group))

    core.main.send_output("\t".join(row_headers))

    # Iterate across each enzyme
    for enzyme in sorted(pathway_info["ENZYME"]):
        row_data = [enzyme]
        for group in groups:
            for tsv_tree in trees:
                row_data.append(str(tsv_tree["children"][group]["enzymes"].get(enzyme, 0)))

        core.main.send_output("\t".join(row_data))

    # Render aggregate totals
    agg_data = ["Participation"]
    for group in groups:
        for tsv_tree in trees:
            shared_ec = pathway_info["ENZYME"].intersection(tsv_tree["children"][group]["enzymes"])
            participation = str(round(len(shared_ec) / len(pathway_info["ENZYME"]), rounding))
            agg_data.append(participation)

    core.main.send_output("\t".join(agg_data))


def render_all(pathways_info, trees, rounding):
    """ Render the percentage of participation per entry for all pathways """
    # Use first report for group list
    groups = trees[0]["children"].keys()

    # Render column headers
    row_headers = ["ID", "Pathway"]
    for group in groups:
        for tsv_tree in enumerate(trees):
            row_headers.append("{0}:{1}".format(tsv_tree[0] + 1, group))

    core.main.send_output("\t".join(row_headers))

    # Iterate across each pathway
    for pathway in sorted(pathways_info):
        row_data = [pathway, pathways_info[pathway]["NAME"][0]]

        # Calculate aggregate totals for each group
        for group in groups:
            for tsv_tree in trees:
                pathway_enzymes = pathways_info[pathway]["ENZYME"]
                shared_ec = pathway_enzymes.intersection(tsv_tree["children"][group]["enzymes"])
                participation = str(round(len(shared_ec) / len(pathway_enzymes), rounding))
                row_data.append(participation)

        core.main.send_output("\t".join(row_data))
