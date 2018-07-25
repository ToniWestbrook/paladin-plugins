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

# Download custom UniProt reports

import argparse
import shlex
import requests
import sys
import core.main

UNIPROT_TEMPLATE_INIT = "https://www.uniprot.org/uploadlists/"
UNIPROT_TEMPLATE_JOB = "https://www.uniprot.org/jobs/{0}.stat"
UNIPROT_TEMPLATE_RESULTS = "https://www.uniprot.org/uniprot/?query=job:{0}&format=tab&columns={1}"
UNIPROT_REQUIRED = ["entry%20name", "id"]
UNIPROT_COLUMNS = ["organism", "protein%20names", "genes", "pathway", "features", "go", "reviewed", "existence", "comments", "database(KEGG)", "database(GeneID)", "database(PATRIC)", "database(EnsemblBacteria)"]


def plugin_connect(definition):
    definition.name = "uniprot"
    definition.description = "Download custom UniProt reports"
    definition.version_major = 1
    definition.version_minor = 0
    definition.version_revision = 1
    definition.dependencies = []

    definition.callback_args = uniprot_args
    # definition.callback_init = uniprot_init
    definition.callback_main = uniprot_main


def uniprot_args(subargs):
    # Parse arguments
    arg_parser = argparse.ArgumentParser(description="PALADIN Pipeline Plugins: UniProt", prog="uniprot")
    arg_parser.add_argument("-i", dest="input", metavar="INPUT", type=str, required=True, help="Input SAM file")
    arg_parser.add_argument("-c", dest="columns", metavar="COLUMN", nargs="+", default=UNIPROT_COLUMNS, type=str, help="Columns to retrieve")
    arg_parser.add_argument("-b", dest="batch", metavar="BATCH", type=int, default=5000, help="Number of entries per submission batch")
    arg_parser.add_argument("-r", dest="retry", metavar="RETRY", type=int, default=10, help="Number of times to retry after http error")

    return arg_parser.parse_known_args(shlex.split(subargs))


def uniprot_main(args):
    # Aggregate SAM data
    sam_data, sam_total = aggregate_sam(args[0].input)

    # Retrieve UniProt info
    columns = UNIPROT_REQUIRED
    columns.extend(args[0].columns)
    uniprot_data = retrieve_data(list(sam_data.keys()), columns, args[0].batch, args[0].retry)

    # Join SAM data with Uniprot
    joined_data = join_data(sam_data, sam_total, uniprot_data)

    # Sort and render report
    print(uniprot_data)
    headers = "Count\tAbundance\tQuality (Average)\tQuality (Max)\tUniProtKB\tID\t{0}".format("\t".join(uniprot_data["Entry name"][2:]))
    join_sorted = sorted(joined_data, reverse=True)
    render_report(join_sorted, headers)


def aggregate_sam(sam_file):
    """ Group, count and average SAM data by reference hit """
    ret_aggregate = dict()
    ret_total = 0

    # Retrieve mapped SAM entries for report
    core.main.send_output("Gathering SAM data...", "stderr")
    entries = core.main.SamEntry.get_entries(sam_file, 0)

    # Iterate through each read/hit
    for read_entry in entries:
        sam_entry = entries[read_entry]
        sam_reference = sam_entry.reference
        if "|" in sam_reference:
            sam_reference = sam_reference.split("|")[2]

        if sam_reference not in ret_aggregate:
            ret_aggregate[sam_reference] = [0, 0.0, 0]

        mapq = sam_entry.mapqual
        ret_aggregate[sam_reference][0] += 1
        ret_aggregate[sam_reference][1] += mapq
        if mapq > ret_aggregate[sam_reference][2]:
            ret_aggregate[sam_reference][2] = mapq

        ret_total += 1

    # Finalize averages
    for entry in ret_aggregate:
        ret_aggregate[entry][1] /= ret_aggregate[entry][0]

    return ret_aggregate, ret_total


def retrieve_data(entries, columns, batch, max_retry):
    """ Retrieve data on the requested reference hits from UniProt """
    ret_data = dict()

    for batch_idx in range(0, len(entries), batch):
        # Ready this batch of entries
        batch_entries = ""
        for entry_idx in range(batch_idx, batch_idx + batch):
            if entry_idx >= len(entries):
                break
            batch_entries += "{0},".format(entries[entry_idx])

        # Construct POST data
        post_data = dict()
        post_data["uploadQuery"] = batch_entries[:-1]
        post_data["format"] = "job"
        post_data["from"] = "ACC+ID"
        post_data["to"] = "ACC"
        post_data["landingPage"] = "false"

        # Submit batch query, retrieve job ID
        retry = 0

        while True:
            try:
                end_batch = batch_idx + batch
                if end_batch > len(entries):
                    end_batch = len(entries)

                core.main.send_output("Fetching entries {0}:{1} of {2}...".format(batch_idx, end_batch, len(entries)), "stderr")

                response = requests.post(UNIPROT_TEMPLATE_INIT, post_data)
                job_id = response.text
                response.close()

                # Check for valid job ID
                if any(x in job_id for x in ["<", " "]):
                    raise Exception

                # Monitor if job has completed
                url = UNIPROT_TEMPLATE_JOB.format(job_id)
                job_complete = ""

                while job_complete != "COMPLETED":
                    response = requests.post(url)
                    job_complete = response.text
                    response.close()

                # Fetch data and breakout text
                url = UNIPROT_TEMPLATE_RESULTS.format(job_id, ",".join(columns))
                response = requests.post(url)

                for line in response.text.split("\n"):
                    line = line.rstrip()
                    fields = line.split("\t")
                    ret_data[fields[0]] = fields

                response.close()

                # End retry processing
                break
            except KeyboardInterrupt:
                raise
            except Exception:
                if retry > max_retry:
                    core.main.send_output("Too many HTTP errors, quitting...", "stderr")
                    sys.exit(1)
                else:
                    core.main.send_output("HTTP error, retrying...", "stderr")

                retry += 1

    return ret_data


def join_data(sam_entries, total, uniprot):
    """ Join local data with UniProt data """
    ret_data = list()

    for kbid in sam_entries:
        sam_entry = sam_entries[kbid]
        join_entry = [sam_entry[0], sam_entry[0] / total * 100, sam_entry[1], sam_entry[2]]

        # Check for UniProt data
        if kbid in uniprot:
            join_entry.extend(uniprot[kbid])
        else:
            join_entry.append(kbid)

        ret_data.append(join_entry)

    return ret_data


def render_report(data, header):
    """ Render report """
    core.main.send_output(header)

    for entry in data:
        core.main.send_output("\t".join([str(item) for item in entry]))
