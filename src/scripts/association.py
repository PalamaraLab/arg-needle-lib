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


"""ARG-based association: Tests all bitsets or random bitsets by mutation rate

Steps:
- Read in the residualised phenotype
- Load the samples file, use overlap to make the residual vector and mask vector
- Run association on the residual vector with the calibration factor

To run as plain linear regression, pass in `--calibration_factor 1` and pass in
the phenotype with residualized covariates to `--residualised_pheno_path`.

For LMM association, we instead get the calibration factor from other output (e.g.
BOLT-LMM) and get the LOCO-residualized phenotype after using BOLT-LMM + PLINK or
other software. This is shown in the Snakefile.

Set at most one of `--min_mac` and `--min_maf`, both cannot be set.

`--sampling_rate` can function in 2 modes. If > 0, it uses that value as
a mutation rate. If = 0, it tests all branches / bitsets in the ARG for
association.

`--haps_threshold` is a p-value threshold, for which we write out all bitsets that
are more significant than this threshold (so variant i is written if and only if
p_i < threshold).

The output is a .tab.gz file containing all association statistics, and a .haps.gz
file containing those variants that are more significant than the threshold. The
`--arg_id` string is used in the labeling of variants in this output.
"""

import argparse
import logging
import numpy as np
import pandas as pd

import arg_needle_lib

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')


# Disable default help, credit to https://stackoverflow.com/a/57582191/
parser = argparse.ArgumentParser(description="Run association using the ARG.", add_help=False)
basic = parser.add_argument_group("Basic arguments")
optional = parser.add_argument_group("Optional arguments")

# Basic options
basic.add_argument("--arg_path", action="store", default="files/chr1.chunk1.argn",
    help="Path to .argn (ARG format) or .trees/.tsz (tskit/tszip format) file (default='files/chr1.chunk1.argn').")
basic.add_argument("--arg_sample_path", action="store", default="files/example.sample",
    help="Path to file listing diploid samples in the same order as haploid leaves of the ARG (default='files/example.sample').")
basic.add_argument("--arg_id", action="store", default="chr1.chunk1",
    help="Identifier for input ARG used for writing output." \
    " Must start with chromosome number followed by a period, e.g., chr1.* or chr22.* (default='chr1.chunk1').")
basic.add_argument("--residualised_pheno_path", action="store", default="files/example.pheno",
    help="Path to the phenotype, with any covariates residualised out." \
    " Should also residualise out the LOCO predictor for ARG-MLMA workflows (default='files/example.pheno').")
basic.add_argument("--out_path", action="store", default="files/results",
    help="Output path prefix for *.tab.gz and *.haps.gz output (default='files/results').")

# Add back help, credit to https://stackoverflow.com/a/57582191/
optional.add_argument(
    "-h", "--help", action="help", default=argparse.SUPPRESS, help="show this help message and exit")

# Advanced options
optional.add_argument("--sampling_rate", action="store", default=0, type=float,
    help="Mutation rate-based testing." \
    " If positive, generates random mutations on the ARG, which are tested for association." \
    " If 0, tests all branches (default=0).")
optional.add_argument("--min_mac", action="store", default=-1, type=int,
    help="Minimum MAC to test (default=-1 which means no minimum). Cannot set both min_mac and min_maf.")
optional.add_argument("--min_maf", action="store", default=-1, type=float,
    help="Minimum MAF to test (default=-1 which means no minimum). Cannot set both min_mac and min_maf.")
optional.add_argument("--max_maf", action="store", default=-1, type=float,
    help="Maximum MAF to test (default=-1 which means no maximum).")
optional.add_argument("--haps_threshold", action="store", default=0, type=float,
    help="p-value threshold, variants with p-value lower than this threshold are written to .haps.gz output (default=0)")
optional.add_argument("--random_seed", action="store", default=1, type=int,
    help="Seed to use for mutation-based testing. If zero, seeds by time (default=1).")
optional.add_argument("--calibration_factor", action="store", default=1, type=float,
    help="Calibration factor. Used in ARG-MLMA workflows, can be set to 1 to instead run linear regression (default=1).")
optional.add_argument("--descendant_list_threshold", action="store", default=-1, type=int,
    help="Advanced parameter determining how ARG branch variants are represented in memory." \
    " Uses a vector of sample IDs below and a bitset above the threshold." \
    " If less than 0, always uses a vector of sample IDs (default)." \
    " If memory is concern, setting to num_individuals/16 may reduce memory at the cost of greater runtime.")

def main():
    args = parser.parse_args()
    logging.info("Command-line args:")
    args_to_print = vars(args)
    for k in sorted(args_to_print):
        logging.info(k + ": " + str(args_to_print[k]))

    # Set up parameters
    chromosome_field = args.arg_id.split(".")[0]
    chrom_num = int(chromosome_field[3:])
    snp_prefix = args.arg_id

    # Read in the residualised phenotype
    # Columns should be FID / IID / phenotype
    pheno_df = pd.read_csv(args.residualised_pheno_path, sep="\s+|\t+", index_col=0, engine="python")
    pheno_df = pheno_df.drop(pheno_df.columns[0], axis=1)
    pheno_df.rename(columns={pheno_df.columns[0]: "PHENO"}, inplace=True)
    pheno_df = pheno_df.dropna()
    residual_dict = dict(zip(pheno_df.index, pheno_df["PHENO"]))

    # Load sample IDs
    sample_ids = []
    with open(args.arg_sample_path, "r") as infile:
        for i, line in enumerate(infile):
            if i >= 2:
                sample_ids.append(int(line.strip("\n").split(" ")[0].split("\t")[0]))

    # Make use_sample and residual lists
    use_sample = []
    residual = []
    for sample_id in sample_ids:
        if sample_id in residual_dict:
            residual.append(residual_dict[sample_id])
            use_sample.append(True)
        else:
            residual.append(-9) # The value here does not matter, it is ignored
            use_sample.append(False)

    logging.info(f"Was able to find phenotype for {sum(use_sample)} of {len(use_sample)} diploid samples")

    # Load ARG and possibly change DescendantList threshold
    if args.arg_path.endswith(".argn"):
        arg = arg_needle_lib.deserialize_arg(args.arg_path)
    elif args.arg_path.endswith(".trees"):
        import tskit
        ts = tskit.load(args.arg_path)
        arg = arg_needle_lib.tskit_to_arg(ts)
    elif args.arg_path.endswith(".tsz"):
        import tszip
        ts = tszip.decompress(args.arg_path)
        arg = arg_needle_lib.tskit_to_arg(ts)
    else:
        raise ValueError(f"Expected .argn / .trees / .tsz file, found {args.arg_path} instead.")
    arg.populate_children_and_roots()

    if args.descendant_list_threshold >= 0:
        arg_needle_lib.DescendantList.set_threshold(args.descendant_list_threshold)
        arg_needle_lib.DescendantList.print_threshold()

    if args.min_maf != -1 and args.min_mac != -1:
        raise ValueError("Can't set both min_maf and min_mac")
    elif args.min_mac != -1:
        min_maf = args.min_mac / (2 * np.sum(use_sample))
    else:
        min_maf = args.min_maf
    max_maf = args.max_maf

    logging.info(f"Using minimum MAF = {min_maf} (note: -1 means no minimum)")
    logging.info(f"Using maximum MAF = {max_maf} (note: -1 means no maximum)")
    logging.info("Starting to run association")
    if args.sampling_rate == 0:
        logging.info("Testing all clades of the ARG")
        max_chi2 = arg_needle_lib.association_diploid_all(
            arg, residual, use_sample, args.out_path,
            chrom_num, snp_prefix, min_maf=min_maf, max_maf=max_maf,
            write_bitset_threshold=args.haps_threshold,
            calibration_factor=args.calibration_factor,
            concise_pvalue=True)
    else:
        assert args.sampling_rate > 0
        logging.info(f"Testing sampled mutations with rate {args.sampling_rate} and seed {args.random_seed} (note: seed 0 means use system time to seed)")
        max_chi2 = arg_needle_lib.association_diploid_mutation(
            arg, residual, use_sample, args.out_path,
            [args.sampling_rate], args.random_seed,
            chrom_num, snp_prefix, min_maf=min_maf, max_maf=max_maf,
            write_bitset_threshold=args.haps_threshold,
            calibration_factor=args.calibration_factor,
            concise_pvalue=True)[0]
    logging.info("Done running association")
    logging.info("Check {}.tab.gz and {}.haps.gz for output.".format(
        args.out_path, args.out_path))

    logging.info("Finished!")

if __name__ == "__main__":
    main()
