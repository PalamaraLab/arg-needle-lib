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


"""Run from this directory using

python3 example_grm.py
"""

# make numpy use 1 thread for timing purposes
import os
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import argparse
import logging
import msprime
import numpy as np

import arg_needle_lib
from utils import time_and_print

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

parser = argparse.ArgumentParser(description='Build ARG-GRMs.')
parser.add_argument("--random_seed", help="Seed to use for simulations and Monte Carlo ARG-GRMs (if zero, seeds Monte Carlo by time).", action="store", default=1234, type=int)
parser.add_argument("--haploid", help="GRM will be haploid instead of diploid.", action='store_true')
parser.add_argument("--alpha", help="GRM alpha parameter", action="store", default=-1, type=float)
parser.add_argument("--mc_mu", help="Mutation rate for Monte Carlo ARG-GRM.", action="store", default=5e-8, type=float)
parser.add_argument("--out_path", help="Path to store GRM in GCTA format.", action="store", default="temp")

if __name__ == '__main__':
    args = parser.parse_args()
    logging.info("Command-line args:")
    args_to_print = vars(args)
    for k in sorted(args_to_print):
        logging.info(k + ": " + str(args_to_print[k]))

    seed = args.random_seed
    diploid = not args.haploid
    alpha = args.alpha
    monte_carlo_mu = args.mc_mu
    out_path = args.out_path
    num_samples = 500
    length = 5e6

    logging.info("Starting simulation")
    ts = msprime.simulate(
        sample_size=num_samples, # number of haploid samples, must be even if using diploid GRM
        Ne=1e4,
        length=length,
        recombination_rate=1e-8,
        mutation_rate=1e-8,
        random_seed=seed)

    logging.info(str(ts.num_trees) + " trees, " + str(ts.num_nodes) + " nodes")
    arg = arg_needle_lib.tskit_to_arg(ts)
    arg.populate_children_and_roots()
    logging.info("Done with populating children and roots")

    with time_and_print("compute exact ARG GRM"):
        exact_arg_grm = arg_needle_lib.exact_arg_grm(arg, alpha=alpha, diploid=diploid)

    with time_and_print("compute Monte Carlo ARG-GRM"):
        mc_arg_grm = arg_needle_lib.monte_carlo_arg_grm(
            arg, monte_carlo_mu=monte_carlo_mu, seed=seed, alpha=alpha, diploid=diploid)

    logging.info(f"GRM shapes: {exact_arg_grm.shape} {mc_arg_grm.shape}")
    logging.info("First few entries...")
    for grm in [exact_arg_grm, mc_arg_grm]:
        print(grm[:3, :3])

    def matrix_rmse(mat1, mat2):
        return np.sqrt(np.mean(np.square(mat1 - mat2)))

    rmse = matrix_rmse(mc_arg_grm, exact_arg_grm)
    logging.info(f"RMSE for Monte Carlo ARG-GRM vs exact ARG-GRM: {rmse}")
    # RMSE decreases with higher --mc_mu

    # Now create MAF-stratified ARG-GRMs and check that the sum is the same
    maf_bin_boundaries = [0, 0.0025, 0.01, 0.05, 0.5]
    num_grms = len(maf_bin_boundaries) - 1
    stratified_grms = []
    with time_and_print("compute MAF-stratified ARG-GRMs"):
        logging.info("")
        logging.info("Computing MAF-stratified ARG-GRMs")
        for i in range(num_grms):
            logging.info("Computing for MAF bin boundaries ({}, {}]".format(
                maf_bin_boundaries[i], maf_bin_boundaries[i+1]))
            stratified_grms.append(arg_needle_lib.monte_carlo_arg_grm(
                arg, monte_carlo_mu=monte_carlo_mu, seed=seed, alpha=alpha, diploid=diploid,
                centering=False, # important to wait until after combining to center
                min_maf=maf_bin_boundaries[i], max_maf=maf_bin_boundaries[i+1]))
    combined_grm = np.sum(stratified_grms, axis=0)
    arg_needle_lib.row_column_center(combined_grm)
    arg_needle_lib.gower_center(combined_grm)
    rmse = matrix_rmse(mc_arg_grm, combined_grm)
    logging.info(f"RMSE for single ARG-GRM vs sum of MAF-stratified ARG-GRMs: {rmse}")
    # These two methods use identical mutations so the RMSE should be tiny

    logging.info(f"Writing Monte Carlo GRM to {args.out_path}.grm.*")
    arg_needle_lib.write_grm(mc_arg_grm, arg.num_mutations(), args.out_path)

    with time_and_print("compute exact ARG-GRM in a small region"):
        exact_arg_grm_region = arg_needle_lib.exact_arg_grm(arg, alpha=alpha, from_pos=length*0.49, to_pos=length*0.51,
            diploid=diploid)
        # arg_needle_lib.write_grm(exact_arg_grm_region, length, args.out_path)
