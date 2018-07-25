#! /usr/bin/env python3

"""
The MIT License

Copyright (c) 2018 by Anthony Westbrook, University of New Hampshire <anthony.westbrook@unh.edu>

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

import argparse
import sys
import textwrap
import core.main
from core.filestore import FileStore


def parse_arguments():
    """ Parse and clean arguments """
    arg_parser = argparse.ArgumentParser(description="PALADIN Pipeline Plugins")
    arg_parser.add_argument("-l", dest="list", action="store_true", help="List available plugins")
    arg_parser.add_argument("--debug", dest="debug", action="store_true", help=argparse.SUPPRESS)
    arg_parser.add_argument("plugins", metavar="Plugins", type=str, nargs=argparse.REMAINDER, help="Plugin with arguments (@@plugin arguments)")

    ret_arguments = arg_parser.parse_args()

    # Check for no arguments, or bad syntax in plugins list
    invalid_args = False
    if not ret_arguments.list and not ret_arguments.plugins:
        invalid_args = True
    if ret_arguments.plugins and not ret_arguments.plugins[0].startswith("@@"):
        invalid_args = True

    if invalid_args:
        arg_parser.print_help()
        sys.exit(1)

    # Replace quotes for strings with spaces
    for plug_idx in range(len(ret_arguments.plugins)):
        plug_arg = ret_arguments.plugins[plug_idx]
        if " " in plug_arg:
            plug_arg = "\"{0}\"".format(plug_arg)
        ret_arguments.plugins[plug_idx] = plug_arg

    return ret_arguments


def get_pipeline(args):
    """ Extract pipeline plugin/argument pairs from arguments """
    ret_pipeline = list()

    str_arguments = " ".join(args)
    plugin_groups = str_arguments.split("@@")
    if plugin_groups[0] != "":
        return ret_pipeline

    for plugin_group in plugin_groups[1:]:
        plugin_group = plugin_group.rstrip()
        tokens = plugin_group.split(" ")
        plugin = tokens[0]
        plugin_args = " ".join(tokens[1:])
        ret_pipeline.append((plugin, plugin_args))

    return ret_pipeline


def list_plugins():
    """ List available plugins """
    desc_start = 21
    plugin_list = core.main.plugin_modules.values()
    sorted_plugins = sorted(plugin_list, key=lambda x: x.name)

    print("The following plugins are available:")

    for plugin in sorted_plugins:
        version_text = "{0}.{1}.{2}".format(plugin.version_major, plugin.version_minor, plugin.version_revision)
        header_text = "{0} ({1}):".format(plugin.name, version_text)
        raw_text = "{0}\n{1}{2}".format(header_text, " " * (desc_start - len(header_text)), plugin.description)

        wrapper = textwrap.TextWrapper(initial_indent="  ", subsequent_indent=" " * (desc_start + 3), width=78)
        render_text = wrapper.wrap(raw_text)

        for line in render_text:
            print(line)

# Parse arguments
arguments = parse_arguments()

# Connect to plugins
core.main.connect_plugins(arguments.debug)

# Handle non-pipeline actions
pipeline_present = True
if arguments.list:
    pipeline_present = False
    list_plugins()

if pipeline_present:
    pipeline = get_pipeline(arguments.plugins)

    try:
        core.main.exec_pipeline(pipeline)
        # Initialize plugins
        #if core.main.init_plugins(set(list(zip(*pipeline))[0])):
            # Execute pipeline
        #    core.main.exec_pipeline(pipeline)

        # Do a final flush to standard out
        core.main.exec_pipeline([("flush", "")])

    finally:
        # Clean up FileStore
        FileStore.destroy()
