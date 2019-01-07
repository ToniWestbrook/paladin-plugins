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

# Core module is responsible for plugin interopability, as well as standard API

import pkgutil
import importlib
import plugins
import pprint
import re
import sys
from core.filestore import FileStore

class PluginDef:
    """ Plugin definition to be populated by each plugin via its plugin_init method at load time """
    def __init__(self, module):
        # Plugin fields
        self.module = module
        self.name = ""
        self.description = ""
        self.version_major = 0
        self.version_minor = 0
        self.version_revision = 0
        self.dependencies = list()

        # Callbacks
        self.callback_args = None
        self.callback_init = None
        self.callback_main = None


class SamEntry:
    """ Store data per individual SAM entry """
    FIELD_QNAME = 0
    FIELD_FLAG = 1
    FIELD_RNAME = 2
    FIELD_POS = 3
    FIELD_MAPQ = 4
    FIELD_CIGAR = 5
    FIELD_RNEXT = 6
    FIELD_PNEXT = 7
    FIELD_TLEN = 8
    FIELD_SEQ = 9
    FIELD_QUAL = 10

    def __init__(self):
        self.query = ""
        self.flag = 0
        self.reference = ""
        self.pos = 0
        self.mapqual = 0
        self.cigar = ""
        self.nextref = ""
        self.nextpos = 0
        self.length = 0
        self.sequence = ""
        self.readqual = 0
        self.frame = ""

    @staticmethod
    def get_entries(filename, quality):
        """ Public API for obtaining SAM data (will handle caching internally) """
        # Check if in cache
        if not (filename, quality) in SamEntry._entries:
            cache = SamEntry.populate_entries(filename, quality)
            SamEntry._entries[(filename, quality)] = cache

        return SamEntry._entries[(filename, quality)]

    @staticmethod
    def populate_entries(filename, quality):
        """ Store SAM entries filtered for the requested quality """
        ret_entries = dict()

        # Open SAM file
        with open(filename, "r") as handle:
            for line in handle:
                fields = line.rstrip().split("\t")

                # Skip header and malformed lines
                if line.startswith("@"):
                    continue
                if len(fields) < 11:
                    continue

                # Filter for minimum quality
                if fields[SamEntry.FIELD_RNAME] == "*" and quality != -1:
                    continue
                if quality != -1 and int(fields[SamEntry.FIELD_MAPQ]) < quality:
                    continue

                # Remove PALADIN frame header since best scoring frame may change between alignments
                header_match = re.search("(.*?:.*?:.*?:)(.*)", fields[SamEntry.FIELD_QNAME])

                entry = SamEntry()
                entry.query = header_match.group(2)
                entry.flag = int(fields[SamEntry.FIELD_FLAG])

                # Fill in entry information if mapped
                if entry.is_mapped:
                    entry.reference = fields[SamEntry.FIELD_RNAME]
                    entry.pos = SamEntry.get_sam_int(fields[SamEntry.FIELD_POS])
                    entry.mapqual = SamEntry.get_sam_int(fields[SamEntry.FIELD_MAPQ])
                    entry.cigar = fields[SamEntry.FIELD_CIGAR]
                    entry.nextref = fields[SamEntry.FIELD_RNEXT]
                    entry.nextpos = fields[SamEntry.FIELD_PNEXT]
                    entry.length = SamEntry.get_sam_int(fields[SamEntry.FIELD_TLEN])
                    entry.sequence = fields[SamEntry.FIELD_SEQ]
                    entry.readqual = SamEntry.get_sam_int(fields[SamEntry.FIELD_QUAL])
                    entry.frame = header_match.group(1)

                # Each read can have multiple non-linear/chimeric hits - store as tuple for ease of processing
                read_base = header_match.group(2)
                hit_idx = 0
                while (read_base, hit_idx) in ret_entries: 
                    hit_idx += 1

                ret_entries[(read_base, hit_idx)] = entry

        return ret_entries

    @staticmethod
    def get_sam_int(val):
        if val.isdigit():
            return int(val)

        return 0

    def is_mapped(self):
        return self.flag & 0x04 > 0

    _entries = dict()


class PaladinEntry:
    """ PALADIN UniProt entry """
    FIELD_COUNT = 0
    FIELD_ABUNDANCE = 1
    FIELD_QUALAVG = 2
    FIELD_QUALMAX = 3
    FIELD_KB = 4
    FIELD_ID = 5
    FIELD_SPECIES = 6
    FIELD_PROTEIN = 7
    FIELD_ONTOLOGY = 11

    TYPE_UNKNOWN = 0
    TYPE_UNIPROT_EXACT = 1
    TYPE_UNIPROT_GROUP = 2
    TYPE_CUSTOM = 3

    def __init__(self):
        self.type = PaladinEntry.TYPE_UNKNOWN
        self.id = "Unknown"
        self.kb = "Unknown"
        self.count = 0
        self.abundance = 0.0
        self.quality_avg = 0.0
        self.quality_max = 0
        self.species_id = "Unknown"
        self.species_full = "Unknown"
        self.protein = "Unknown"
        self.ontology = list()

    @staticmethod
    def get_entries(filename, quality, pattern=None):
        """ Public API for obtaining UniProt report data (will handle caching internally) """
        # Check if in cache
        if not (filename, quality, pattern) in PaladinEntry._entries:
            cache = PaladinEntry.populate_entries(filename, quality, pattern)
            PaladinEntry._entries[(filename, quality, pattern)] = cache

        return PaladinEntry._entries[(filename, quality, pattern)]

    @staticmethod
    def populate_entries(filename, quality, pattern):
        """ Cache this UniProt report data for the requested quality """
        ret_entries = dict()

        # Open UniProt report, skip header
        with open(filename, "r") as handle:
            handle.readline()

            for line in handle:
                fields = line.rstrip().split("\t")

                # Filter for minimum quality
                if float(fields[PaladinEntry.FIELD_QUALMAX]) < quality:
                    continue

                entry = PaladinEntry()
                entry.count = int(fields[PaladinEntry.FIELD_COUNT])
                entry.abundance = float(fields[PaladinEntry.FIELD_ABUNDANCE])
                entry.qual_avg = float(fields[PaladinEntry.FIELD_QUALAVG])
                entry.qual_max = int(fields[PaladinEntry.FIELD_QUALMAX])
                entry.kb = fields[PaladinEntry.FIELD_KB]

                if len(fields) > 10:
                    # Existence of fields indicates a successful UniProt parse by PALADIN
                    if "_9" in entry.kb: entry.type = PaladinEntry.TYPE_UNIPROT_GROUP
                    else: entry.type = PaladinEntry.TYPE_UNIPROT_EXACT

                    entry.species_id = entry.kb.split("_")[1]
                    entry.species_full = fields[PaladinEntry.FIELD_SPECIES]
                    entry.id = fields[PaladinEntry.FIELD_ID]
                    entry.protein = fields[PaladinEntry.FIELD_PROTEIN]
                    entry.ontology = [term.strip() for term in fields[PaladinEntry.FIELD_ONTOLOGY].split(";")]
                else:
                    # Check for custom match
                    if pattern:
                        match = re.search(pattern, entry.kb)
                        if match:
                            entry.type = PaladinEntry.TYPE_CUSTOM
                            entry.species_id = match.group(1)
                            entry.species_full = match.group(1)

                ret_entries[fields[PaladinEntry.FIELD_KB]] = entry

        return ret_entries

    _entries = dict()


# Plugins internal to core, and loaded external modules
internal_plugins = dict()
plugin_modules = dict()

# Standard output and error buffers
output_stdout = list()
output_stderr = list()
console_stdout = False
console_stderr = True


def connect_plugins(debug):
    """ Search for all modules in the plugin package (directory), import each and run plugin_connect method """
    # Initialize File Store
    FileStore.init("pp-", "~/.paladin-plugins", ".", 30)

    # Add internal core plugins
    internal_plugins["flush"] = render_output
    internal_plugins["write"] = render_output

    # Import all external plugin modules in package (using full path)
    for importer, module, package in pkgutil.iter_modules(plugins.__path__):
        try:
            module_handle = importlib.import_module("{0}.{1}".format(plugins.__name__, module))
            if "plugin_connect" in dir(module_handle):
                plugin_modules[module] = PluginDef(module_handle)

        except Exception as exception:
            if debug:
                raise exception
            else:
                send_output("Error loading \"{0}.py\", skipping...".format(module), "stderr")

    # Connect to all external plugins
    for plugin in plugin_modules:
        plugin_modules[plugin].module.plugin_connect(plugin_modules[plugin])


def args_plugins(plugins):
    """ Run argument parsing for each plugin """
    for plugin in plugins:
        plugin_modules[plugin].callback_args()


def init_plugins(plugins):
    """ _initialize plugins being used in this session """
    init_queue = set()
    init_history = set()

    # Scan for plugins and dependencies
    for plugin in plugins:
        if plugin not in plugin_modules and plugin not in internal_plugins:
            if plugin in plugins.modules_disabled:
                print("Disabled plugin: {0}".format(plugin))
            else:
                print("Unknown plugin: {0}".format(plugin))

            return False

        # _initialize external plugins
        if plugin in plugin_modules:
            # Look for dependencies
            init_queue.update(plugin_modules[plugin].dependencies)
            init_queue.add(plugin)

    # _initialize
    for plugin in init_queue:
        if plugin_modules[plugin].callback_init:
            if plugin not in init_history:
                plugin_modules[plugin].callback_init()
                init_history.add(plugin)

    return True


#def exec_pipeline(pipeline):
#    """ Execute requested plugin pipeline """
#    for task in pipeline:
#        if task[0] in internal_plugins:
#            # Internal plugin
#            internal_plugins[task[0]](task[1].strip("\""))
#        else:
#            # External plugin
#            if plugin_modules[task[0]].callback_main:
#                plugin_modules[task[0]].callback_main(task[1])

def exec_pipeline(pipeline):
    """ Execute requested plugin pipeline """
    for task in pipeline:
        if task[0] in internal_plugins:
            # Internal plugin
            internal_plugins[task[0]](task[1].strip("\""))
        elif task[0] in plugin_modules:
            # External plugin
            plugin = plugin_modules[task[0]]

            # Process arguments (this may sys.exit if help mode)
            if plugin.callback_args:
                args = plugin.callback_args(task[1])

            # Process dependencies and initialization
            for dependency in [plugin_modules[x] for x in plugin.dependencies]:
                if dependency.callback_init:
                    dependency.callback_init()

            if plugin.callback_init:
                plugin.callback_init()

            # Execute
            if plugin.callback_main:
                plugin.callback_main(args)
        else:
            # Invalid plugin
            send_output("Invalid plugin \"{0}\"".format(task[0]), "stderr")
            sys.exit(1)


def render_output(filename="", target="stdout"):
    """ The flush internal plugin handles rendering output (to stdout or file) """
    if target == "stdout":
        render_text = "".join(output_stdout)
        del output_stdout[:]
    if target == "stderr":
        render_text = "".join(output_stderr)
        del output_stderr[:]

    if not filename:
        std_target = sys.stdout if target == "stdout" else sys.stderr
        print(render_text, flush=True, file=std_target)
    else:
        with open(filename, "w") as handle:
            handle.write(render_text)


def send_output(output_text, target="stdout", suffix="\n"):
    """ API - Record output into the appropriate buffer """
    new_content = "{0}{1}".format(output_text, suffix)

    if target == "stdout":
        if console_stdout:
            print(new_content, end="", flush=True)
        else:
            output_stdout.append(new_content)
    else:
        if console_stderr:
            print(new_content, end="", flush=True, file=sys.stderr)
        else:
            output_stderr.append(new_content)


def getInteger(val):
    """ API - Return value if string is integer (allows negatives) """
    try:
        return int(val)
    except:
        return None


def debugPrint(obj):
    """ API - Debugging """
    pp = pprint.PrettyPrinter(indent=4)
    pp.pprint(obj)
