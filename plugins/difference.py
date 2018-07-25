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

# Analyze difference between two PALADIN alignemnts

import argparse
import shlex
import core.main
from core.datastore import DataStore


def plugin_connect(definition):
    definition.name = "difference"
    definition.description = "Analyze relative differences between two PALADIN taxonomy reports"
    definition.version_major = 1
    definition.version_minor = 1
    definition.version_revision = 0
    definition.dependencies = ["taxonomy"]

    definition.callback_args = difference_args
    # definition.callback_init = difference_init
    definition.callback_main = difference_main


def difference_args(subargs):
    arg_parser = argparse.ArgumentParser(description="PALADIN Pipeline Plugins: Difference", prog="difference")
    arg_parser.add_argument("-1", dest="basis_files", metavar=("Taxonomy", "Sam", "Uniprot"), type=str, nargs=3, required=True, help="Basis files of comparison (Taxonomy, Sam, UniProt)")
    arg_parser.add_argument("-2", dest="compare_files", metavar=("Taxonomy", "Sam", "Uniprot"), type=str, nargs=3, required=True, help="Compared files of comparison (Taxonomy, Sam, UniProt")
    arg_parser.add_argument("-c", dest="custom", type=str, default="", help="Species-parsing regex pattern for non-UniProt entries")
    arg_parser.add_argument("-s", dest="simple", action="store_true", help="Simple report, do not include mapped/unmapped read contributions (alters totals)")
    return arg_parser.parse_known_args(shlex.split(subargs))


def difference_main(args):
    # Obtain taxonomy differences
    core.main.send_output("Comparing taxonomy reports...", "stderr")
    taxonomy_data, taxonomy_total = taxonomy_compare((args[0].basis_files[0], args[0].compare_files[0]))

    # Obtain SAM differences with Uniprot info for full report
    uniprot_data = list()
    if not args[0].simple:
        core.main.send_output("Gathering SAM data from first alignment", "stderr")
        sam_entries1 = core.main.SamEntry.get_entries(args[0].basis_files[1], -1)
        core.main.send_output("Gathering SAM data from second alignment", "stderr")
        sam_entries2 = core.main.SamEntry.get_entries(args[0].compare_files[1], -1)
        core.main.send_output("Comparing SAM data...", "stderr")
        sam_data = sam_compare(sam_entries1, sam_entries2)

        # Get Uniprot entries
        core.main.send_output("Gathering UniProt data from first alignment", "stderr")
        uniprot_data.append(core.main.PaladinEntry.get_entries(args[0].basis_files[2], 0, args[0].custom))
        core.main.send_output("Gathering UniProt data from second alignment", "stderr")
        uniprot_data.append(core.main.PaladinEntry.get_entries(args[0].compare_files[2], 0, args[0].custom))

    core.main.send_output("Comparing alignments...", "stderr")
    combined_data = combined_compare(taxonomy_data, sam_data, uniprot_data, taxonomy_total)
    render_combined(combined_data, args[0].simple)


def taxonomy_compare(reports):
    """ Taxonomy report comparison """
    ret_difference = dict()
    ret_total = 0
    report_data = (dict(), dict())

    # Read files into memory
    for idx in range(2):
        with open(reports[idx], "r") as handle:
            # Skip header
            handle.readline()

            for line in handle:
                fields = line.rstrip().split("\t")
                report_data[idx][fields[2]] = int(fields[0])
                if idx == 1:
                    ret_total += int(fields[0])

    # Calculate differences (set 1, set 2, shared)
    set1_keys = set(report_data[0].keys())
    set2_keys = set(report_data[1].keys())
    shared_keys = set1_keys.intersection(set2_keys)
    set1_keys = set1_keys.difference(shared_keys)
    set2_keys = set2_keys.difference(shared_keys)

    ret_difference = {key: (report_data[1][key] - report_data[0][key]) for key in shared_keys}
    ret_difference.update({key: -report_data[0][key] for key in set1_keys})
    ret_difference.update({key: report_data[1][key] for key in set2_keys})

    return ret_difference, ret_total


def sam_compare(entries1, entries2):
    """ SAM entries comparison """
    ret_diff = dict()

    # Use set 1 as the baseRead
    for entry in entries1:
        # Iterate through possible supplementary hits until neither available
        read_idx = (entry[0], 0)
        while read_idx in entries1 or read_idx in entries2:
            if read_idx not in entries1:
                ref1 = "*"
                ref2 = entries2[read_idx].reference
            elif read_idx not in entries2:
                ref1 = entries1[read_idx].reference
                ref2 = "*"
            else:
                ref1 = entries1[read_idx].reference
                ref2 = entries2[read_idx].reference

            # Note references that do not match
            if ref1 != ref2:
                diff_pair = (ref1, ref2)
                if diff_pair not in ret_diff:
                    ret_diff[diff_pair] = 0
                ret_diff[diff_pair] += 1

            read_idx = (read_idx[0], read_idx[1] + 1)

    return ret_diff


def combined_compare(taxonomy_entries, sam_entries, uniprot_entries, pass_total):
    """ Combined compare breaks down taxonomy report on per read basis, then aggregates """
    ret_entries = dict()

    # Calculate unmapped difference and add to taxonomy set
    unmapped_diff = 0
    for sam_entry in sam_entries:
        for set_idx in range(2):
            if sam_entry[set_idx] == "*":
                if set_idx == 0:
                    pass_total += 1

                unmapped_diff += [1, -1][set_idx]

    # Add unmapped entry to taxonomy list
    taxonomy_entries["Unmapped"] = unmapped_diff

    # Iterate through each taxonomy grouping
    for taxon_entry in taxonomy_entries:
        # Skip unchanged entries
        if taxonomy_entries[taxon_entry] == 0:
            continue

        # Prepare entry if new
        if taxon_entry not in ret_entries:
            ret_entries[taxon_entry] = (dict(), taxonomy_entries[taxon_entry] / pass_total)

        # Scan SAM entries for matching lineage
        for sam_entry in sam_entries:
            # Lookup species and lineage
            sam_records = list()
            lineage = list()

            for set_idx in range(2):
                if sam_entry[set_idx] == "*":
                    miss_record = type("miss_record", (object,), {})()
                    miss_record.species_id = "Unmapped"
                    miss_record.species_full = "Unmapped"
                    sam_records.append(miss_record)
                    lineage.append("Unmapped")
                else:
                    sam_records.append(uniprot_entries[set_idx][sam_entry[set_idx].split("|")[-1]])
                    species_id = sam_records[set_idx].species_id
                    result = DataStore.get_entry("taxonomy").exec_query("lineage-lookup", [species_id]).fetchone()

                    if result:
                        lineage.append([x.strip() for x in result[0].split(";")])
                    else:
                        lineage.append("Unknown")

            for set_idx in range(2):
                # Attempt to match on species
                if sam_records[set_idx].species_full == taxon_entry:
                    sam_species = sam_records[1 - set_idx].species_full
                    sam_amount = [1, -1][set_idx] * sam_entries[sam_entry]

                    if sam_species not in ret_entries[taxon_entry][0]:
                        ret_entries[taxon_entry][0][sam_species] = (0, 0)

                    sam_parts = ret_entries[taxon_entry][0][sam_species]

                    if sam_amount > 0:
                        ret_entries[taxon_entry][0][sam_species] = (sam_parts[0] + sam_amount, sam_parts[1])
                    else:
                        ret_entries[taxon_entry][0][sam_species] = (sam_parts[0], sam_parts[1] + sam_amount)

                    break

                # Attempt to match on lineage
                for rank_idx in range(len(lineage[set_idx])):
                    if lineage[set_idx][rank_idx] == taxon_entry:
                        compare_idx = rank_idx

                        if compare_idx >= len(lineage[1 - set_idx]):
                            compare_idx = len(lineage[1 - set_idx]) - 1

                        sam_lineage = lineage[1 - set_idx][compare_idx]
                        sam_amount = [1, -1][set_idx] * sam_entries[sam_entry]

                        if sam_lineage not in ret_entries[taxon_entry][0]:
                            ret_entries[taxon_entry][0][sam_lineage] = (0, 0)

                        sam_parts = ret_entries[taxon_entry][0][sam_lineage]

                        if sam_amount > 0:
                            ret_entries[taxon_entry][0][sam_lineage] = (sam_parts[0] + sam_amount, sam_parts[1])
                        else:
                            ret_entries[taxon_entry][0][sam_lineage] = (sam_parts[0], sam_parts[1] + sam_amount)

                        break

    # Calculate percentages
    for taxonomy_entry in ret_entries:
        for sam_entry in ret_entries[taxonomy_entry][0]:
            sam_parts = ret_entries[taxonomy_entry][0][sam_entry]
            ret_entries[taxonomy_entry][0][sam_entry] = (sam_parts[0] / pass_total, sam_parts[1] / pass_total)

    return ret_entries


def align_stats(entries1, entries2):
    """ Align stats returns 3-tuple (Different, Added, Removed) """
    total_diff = 0
    total_add = 0
    total_del = 0

    for entry in entries1:
        if entries1[entry].flag & 0x4 == entries2[entry].flag & 0x4:
            # Matching map types
            if entries1[entry].reference != entries2[entry].reference:
                total_diff += 1
        else:
            if entries1[entry].flag & 0x4 > 0:
                total_add += 1
            else: total_del += 1

    return (total_diff, total_add, total_del)


def render_combined(data, simple):
    """ Render combined report """
    # Render appropriate header
    if simple:
        header = "Taxon\tDifference"
    else:
        header = "Taxon\tDifference\tContributor\tContribution (Pos)\tContribution (Neg)"

    core.main.send_output(header)

    sorted_taxa = sorted(data.items(), key=lambda item: abs(item[1][1]), reverse=True)

    for taxon_entry in sorted_taxa:
        if simple:
            core.main.send_output("{0}\t{1}".format(taxon_entry[0], taxon_entry[1][1]))
        else:
            # Sort SAM data
            sorted_sam = sorted(taxon_entry[1][0].items(), key=lambda item: abs(item[1][0] + item[1][1]), reverse=True)
            for sam_entry in sorted_sam:
                core.main.send_output("{0}\t{1}\t{2}\t{3}\t{4}".format(taxon_entry[0], taxon_entry[1][1], sam_entry[0], sam_entry[1][0], sam_entry[1][1]))
