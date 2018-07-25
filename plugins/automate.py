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

# Automate PALADIN execution across multiple sets of reads

import os
import re
import subprocess
import argparse
import shlex


def plugin_connect(definition):
    definition.name = "automate"
    definition.description = "Automate PALADIN execution across multiple sets of reads"
    definition.version_major = 1
    definition.version_minor = 1
    definition.version_revision = 0

    definition.callback_args = automate_args
    # definition.callback_init = automate_init
    definition.callback_main = automate_main


def automate_args(subargs):
    # Parse arguments
    arg_parser = argparse.ArgumentParser(description="PALADIN Pipeline Plugins: Automate", prog="automate")
    arg_parser.add_argument("reference", metavar="REFERENCE", type=str, help="Reference database")
    arg_parser.add_argument("root", metavar="ROOT", type=str, help="Root path to search")
    arg_parser.add_argument("pattern", metavar="PATTERN", type=str, help="Input reads search pattern")
    arg_parser.add_argument("options", metavar="OPTIONS", type=str, nargs=argparse.REMAINDER, help="PALADIN options")
    return arg_parser.parse_known_args(shlex.split(subargs))


def automate_main(args):
    for root, dirs, files in os.walk(args[0].root):
        for filename in files:
            if not re.search(args[0].pattern, filename):
                continue

            # Matching input sequence, execute PALADIN
            base_name = filename
            if "." in base_name:
                base_name = base_name[:base_name.index(".")]

            base_name = os.path.join(root, base_name)
            full_file = os.path.join(root, filename)

            command = "paladin align {0} {1} -o {2} {3}".format(args[0].reference, full_file, base_name, " ".join(args[0].options))
            output = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True)

            with open("{0}.log".format(base_name), "wb") as handle:
                handle.write(output)
