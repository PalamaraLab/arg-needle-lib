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


"""Form a residual from phenotye and PLINK .sscore files

Read in a phenotype and PLINK prediction (.sscore and .sscore.vars).
Use these to form and output the residual. PLINK scoring is done as
https://www.cog-genomics.org/plink/2.0/score.

The phenotype should have 3 columns of FID / IID [unused] / phenotype

The .sscore should have 3 columns of FID / IID [unused] / SCORE1_AVG

The .sscore.vars are used to count the number of variants to undo a
normalisation within PLINK.
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

parser = argparse.ArgumentParser(description='Form residual from phenotype and prediction.')
parser.add_argument("--pheno_path", help="Input path", action="store")
parser.add_argument("--plink_score_path", help=".sscore path", action="store", default="")
parser.add_argument("--plink_score_variants_path", help=".sscore.vars path", action="store", default="")
parser.add_argument("--residualised_pheno_path", help="Output path", action="store")

if __name__ == "__main__":
    args = parser.parse_args()
    logging.info("Command-line args:")
    args_to_print = vars(args)
    for k in sorted(args_to_print):
        logging.info(k + ": " + str(args_to_print[k]))

    # Everything with phenotype, BLUP, and residual
    # BLUP is the prediction (best linear unbiased predictor)
    pheno_df = pd.read_csv(args.pheno_path, sep="\s+|\t+", index_col=0) # 2021/04/19
    pheno_df = pheno_df.drop(pheno_df.columns[0], axis=1)
    pheno_df.rename(columns={pheno_df.columns[0]: "PHENO"}, inplace=True)

    num_loco_sites = 0
    with open(args.plink_score_variants_path, "r") as infile:
        for line in infile:
            num_loco_sites += 1

    blup_df = pd.read_csv(args.plink_score_path, sep="\s+|\t+", index_col=0)
    blup_df = blup_df.drop(blup_df.columns[0], axis=1)
    blup = -blup_df["SCORE1_AVG"] * num_loco_sites * 2  # note the factor of 2, and the minus
    blup_df["BLUP"] = blup

    combined_df = pheno_df.merge(blup_df, how="inner", left_index=True, right_index=True)
    logging.info(f"Pheno df shape: {pheno_df.shape}")
    logging.info(f"BLUP df shape: {blup_df.shape}")
    logging.info(f"Combined df shape: {combined_df.shape}")
    nrows_before = combined_df.shape[0]

    combined_df = combined_df.dropna()
    nrows_after = combined_df.shape[0]
    logging.info(f"Dropped {nrows_before - nrows_after} rows containing NaNs")

    # Mean center now (using the common samples), then form the residual
    combined_df["PHENO"] -= np.mean(combined_df["PHENO"])
    combined_df["BLUP"] -= np.mean(combined_df["BLUP"])
    combined_df["RESIDUAL"] = combined_df["PHENO"] - combined_df["BLUP"]

    logging.info(combined_df.head())
    logging.info(f"Pheno std: {np.std(combined_df['PHENO'], ddof=1)}")
    logging.info(f"BLUP std: {np.std(combined_df['BLUP'], ddof=1)}")
    logging.info(f"Residual std: {np.std(combined_df['RESIDUAL'], ddof=1)}")
    logging.info(f"Correlation of pheno and BLUP: {np.corrcoef(combined_df['PHENO'], combined_df['BLUP'])[0][1]}")
    logging.info(f"Correlation of pheno and residual: {np.corrcoef(combined_df['PHENO'], combined_df['RESIDUAL'])[0][1]}")
    residual_list = list(zip(combined_df.index, combined_df["RESIDUAL"]))

    # Write out results
    with open(args.residualised_pheno_path, "w") as outfile:
        outfile.write("FID IID residualised_pheno\n")
        for (fid, residual_value) in residual_list:
            outfile.write(f"{fid} {fid} {residual_value}\n")
