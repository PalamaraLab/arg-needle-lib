# This file is part of the ARG-Needle genealogical inference and
# analysis software suite.
# Copyright (C) 2023 ARG-Needle Developers.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


"""Command-line utilities arg2tskit and tskit2arg.
"""

import argparse
import logging

import arg_needle_lib
import tskit  # we assume tskit, but only load tszip if asked for

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')


def add_common_arguments(parser):
    parser.add_argument("--arg_path", help="Path to .argn (ARG format) file",
                        action="store", required=True)
    parser.add_argument("--ts_path", help="Path to .trees/.tsz (tskit/tszip format) file",
                        action="store", required=True)


def log_arguments(args):
    logging.info("Command-line args:")
    args_to_print = vars(args)
    for k in sorted(args_to_print):
        logging.info(k + ": " + str(args_to_print[k]))


def check_common_arguments(args):
    # Check --arg_path and --ts_path
    if not args.arg_path.endswith(".argn"):
        raise ValueError("--arg_path must end in .argn")
    if not (args.ts_path.endswith(".trees") or args.ts_path.endswith(".tsz")):
        raise ValueError("--ts_path must end in .trees or .tsz")


def arg2tskit():
    parser = argparse.ArgumentParser(description='Convert from arg_needle_lib.ARG to tskit.TreeSequence.')
    add_common_arguments(parser)
    args = parser.parse_args()
    log_arguments(args)
    check_common_arguments(args)

    # Get an ARG
    arg = arg_needle_lib.deserialize_arg(args.arg_path)
    # Convert to ts
    ts = arg_needle_lib.arg_to_tskit(arg, batch_size=1024)

    # Save the ts
    if args.ts_path.endswith(".trees"):
        ts.dump(args.ts_path)
    else:
        import tszip
        tszip.compress(ts, args.ts_path)


def tskit2arg():
    parser = argparse.ArgumentParser(description='Convert from tskit.TreeSequence to arg_needle_lib.ARG.')
    add_common_arguments(parser)
    args = parser.parse_args()
    log_arguments(args)
    check_common_arguments(args)

    # Get ts
    if args.ts_path.endswith(".trees"):
        ts = tskit.load(args.ts_path)
    else:
        import tszip
        ts = tszip.decompress(args.ts_path)
    # Convert to ARG
    arg = arg_needle_lib.tskit_to_arg(ts)
    # Save the ARG
    arg_needle_lib.serialize_arg(arg, args.arg_path)
