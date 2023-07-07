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


"""Null model resampling-based testing for computing genome-wide significance
thresholds (Kanai et al. JHG 2016, Churchill and Doerge Genetics 1994).
"""

import getpass
import os
import sys

import click
import logging
import numpy as np
from scipy.stats.distributions import chi2
import time

import arg_needle_lib

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

@click.command()
@click.option("--arg_path", help="Path to .argn (ARG format) or .trees/.tsz (tskit/tszip format) file")
@click.option("--descendant_list_threshold", help="DescendantList threshold for determining how variants are represented. Uses vector of sample IDs below and bitset above the threshold. If less than 0, always uses vector of sample IDs (default).", default=-1, type=int)
@click.option("--min_mac", help="Minimum MAC to test (-1 means no minimum)", default=-1, type=int)
@click.option("--mutation_rate", default=1e-5, type=float,
    help="Mutation rate-based testing." \
    " If positive, generates random mutations on the ARG, which are tested for association." \
    " If 0, tests all branches (default=0).")
@click.option("--start_seed", help="Starting seed", default=1, type=int)
@click.option("--num_seeds", help="Number of seeds", default=1, type=int)
def null_model_test(arg_path, descendant_list_threshold, min_mac, mutation_rate, start_seed, num_seeds):
    """Run a GeWAS for num_seeds random phenotypes and print minimum P-values to stdout"""
    logging.info("Command-line args:")
    logging.info(f"ARG path: {arg_path}")
    logging.info(f"Descendant list threshold: {descendant_list_threshold}")
    logging.info(f"Min MAC: {min_mac}")
    logging.info(f"Mutation rate: {mutation_rate}")
    logging.info(f"Start seed: {start_seed}")
    logging.info(f"N seeds: {num_seeds}")

    # Load ARG and possibly change DescendantList threshold
    if arg_path.endswith(".argn"):
        arg = arg_needle_lib.deserialize_arg(arg_path)
    elif arg_path.endswith(".trees"):
        import tskit
        ts = tskit.load(arg_path)
        arg = arg_needle_lib.tskit_to_arg(ts)
    elif arg_path.endswith(".tsz"):
        import tszip
        ts = tszip.decompress(arg_path)
        arg = arg_needle_lib.tskit_to_arg(ts)
    else:
        raise ValueError(f"Expected .argn / .trees / .tsz file, found {arg_path} instead.")
    arg.populate_children_and_roots()

    if descendant_list_threshold >= 0:
        arg_needle_lib.DescendantList.set_threshold(descendant_list_threshold)
        arg_needle_lib.DescendantList.print_threshold()

    num_samples_haploid = arg.num_samples()
    num_samples = num_samples_haploid // 2
    min_maf = min_mac / num_samples_haploid
    logging.info(f"Using minimum MAF = {min_maf} (note: -1 means no minimum)")
    for seed_offset in range(num_seeds):
        seed = start_seed + seed_offset
        logging.info(f"Phenotype random seed = {seed}")
        np.random.seed(seed)
        phenotype = np.random.randn(num_samples)
        use_sample = [True] * (num_samples)
        if mutation_rate == 0:
            max_chi2 = arg_needle_lib.association_diploid_all(
                arg, phenotype, use_sample, "foo", 1, "foo",
                min_maf=min_maf, max_maf=-1,
                write_bitset_threshold=-1, calibration_factor=1,
                concise_pvalue=True, max_only=True)
            min_p = chi2.sf(max_chi2, 1)
            print(f"Max chi2 and min p for seed {seed}: {max_chi2} {min_p}")
        else:
            assert mutation_rate > 0
            logging.info(f"Testing sampled mutations with rate {mutation_rate} and seed {seed} (note: seed 0 means use system time to seed)")
            max_chi2 = arg_needle_lib.association_diploid_mutation(
                arg, phenotype, use_sample, "foo", [mutation_rate], seed, 1, "foo",
                min_maf=min_maf, max_maf=-1,
                write_bitset_threshold=-1, calibration_factor=1,
                concise_pvalue=True, max_only=True)[0]
            min_p = chi2.sf(max_chi2, 1)
            print(f"Max chi2 and min p for seed {seed}: {max_chi2} {min_p}")
    logging.info("Finished!")

if __name__ == "__main__":
  null_model_test()
