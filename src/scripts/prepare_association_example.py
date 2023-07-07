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


"""Prepare example for running association.
"""

import logging
import numpy as np
import os

import arg_needle_lib
import msprime

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')


NUM_DIPLOID_SAMPLES = 100
CHOICE_ALLELE_COUNT = 20
RANDOM_SEED = 12345
SEQUENCE_LENGTH = 5e6

def write_pheno(phen_to_write, path, header=False):
    """For association.py script, header should be False.
    """
    num_samples = phen_to_write.shape[0]
    with open(path, 'w') as outfile:
        if header:
            outfile.write('\t'.join(["FID", "IID", "pheno"]) + '\n')
        for i in range(num_samples):
            outfile.write('\t'.join([str(i+1), str(i+1), str(phen_to_write[i])]) + '\n')

def write_sample(num_samples, path):
    if not path.endswith(".sample"):
        raise ValueError("Sample path should end in .sample")
    with open(path, 'w') as outfile:
        outfile.write(' '.join(["ID_1", "ID_2", "missing"]) + '\n')
        outfile.write(' '.join(["0", "0", "0"]) + '\n')
        for i in range(num_samples):
            outfile.write(' '.join([str(i+1), str(i+1), "0"]) + '\n')

def main():
    logging.info("Preparing arg_needle_lib ARG association example")

    # Use old msprime API to simulate trees and mutations
    simulations = list(msprime.simulate(
        sample_size=NUM_DIPLOID_SAMPLES*2,
        Ne=15000,
        length=SEQUENCE_LENGTH,
        recombination_rate=1e-8,
        mutation_rate=1e-8, 
        random_seed=RANDOM_SEED,
        num_replicates=10
    ))
    simulation = simulations[0]

    allele_counts = np.array(
        [np.sum(variant.genotypes) for variant in simulation.variants()],
        dtype=np.int64
    )

    # Select a variant with the correct allele count
    logging.info("")
    logging.info("Generating phenotype from chr1.chunk1 ARG")
    np.random.seed(RANDOM_SEED)
    correct_variants = np.nonzero(allele_counts == CHOICE_ALLELE_COUNT)[0]
    if len(correct_variants) == 0:
        raise ValueError("Found no variants with allele count={}".format(CHOICE_ALLELE_COUNT))
    choice_id = np.random.choice(correct_variants, 1, replace=False)[0]
    logging.info("Chose 1 out of {} SNPs with allele count={}".format(
        len(correct_variants), CHOICE_ALLELE_COUNT))
    for i, variant in enumerate(simulation.variants()):
        if i == choice_id:
            choice_geno = variant.genotypes.reshape((NUM_DIPLOID_SAMPLES, 2)).sum(axis=-1)
            break

    # Add a small amount of noise to the phenotype
    phenotype = choice_geno + 0.05 * np.random.randn(NUM_DIPLOID_SAMPLES)

    # Set a non-carrier to NaN to include missing values
    for j in range(NUM_DIPLOID_SAMPLES):
        if choice_geno[j] == 0:
            phenotype[j] = np.nan
        break

    if not os.path.exists('files'):
       os.makedirs('files')
    write_pheno(phenotype, "files/example.pheno")
    write_sample(NUM_DIPLOID_SAMPLES, "files/example.sample")
    logging.info("Wrote files/example.pheno and files/example.sample")

    logging.info("Converting ARGs to arg_needle_lib format and serializing")
    converted_args = [arg_needle_lib.tskit_to_arg(x) for x in simulations]
    # arg_needle_lib.serialize_arg(converted_args[0], "files/chr1.chunk1.argn")
    # Write out ARGs for 5 chromosomes, 2 chunked ARGs per chromosome
    for i in range(5):
        for j in range(2):
            arg = converted_args[i*2+j]
            arg_needle_lib.serialize_arg(arg, f"files/chr{i+1}.chunk{j+1}.argn")
    logging.info("Finished!")

if __name__ == "__main__":
    main()
