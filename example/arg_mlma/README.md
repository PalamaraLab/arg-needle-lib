# ARG-MLMA example

This directory contains files used to run ARG-based mixed-linear-model association (ARG-MLMA), as described in https://doi.org/10.1038/s41588-023-01379-x.

Note that the association.py script expects a residualised phenotype where all covariates and the BLUP have been regressed out. This residualised phenotype is constructed using BOLT-LMM and PLINK, with the whole process wrapped in the Snakemake workflow. You will need to fill in a few paths from your own set-up and cluster environment to for this. The workflow also collects the BOLT-LMM calibration factor which is used to calibrate the statistics.

This directory also includes scripts to compute the genome-wide significance threshold based on null model permutation testing. The process is conducted in three steps: First, we perform a single GeWAS genealogy-wide on a random phenotype. Next, we run 1,000 random phenotypes on a small chunk of ARG to determine multiple-testing thresholds for that chunk. Finally, we use the genealogy-wide results to extrapolate this threshold to the whole genome. We do this for several MAF cut-offs, but in most cases one will want to focus on the "MAF=0" threshold.
