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

"""Python functions dealing with ARG-GRMs."""

import gzip
import logging
import numpy as np
import arg_needle_lib_pybind

__all__ = [
    "haploid_grm_to_diploid",
    "row_column_center",
    "gower_center",
    "exact_arg_grm",
    "monte_carlo_arg_grm",
    "monte_carlo_arg_grm_numpy",
    "write_grm",
]


def haploid_grm_to_diploid(hap_mat):
    """Turns a haploid GRM into a diploid GRM by treating neighboring IDs as haploid pairs.
    """
    assert len(hap_mat.shape) == 2
    assert hap_mat.shape[0] == hap_mat.shape[1]
    assert hap_mat.shape[0] % 2 == 0
    dip_mat = hap_mat[::2, ::2] + hap_mat[1::2, ::2] + hap_mat[::2, 1::2] + hap_mat[1::2, 1::2]
    return dip_mat


def row_column_center(mat):
    """Apply row-column centering to a numpy square matrix in-place.

    Also returns the matrix.

    Note: row column centering and Gower centering commute, i.e.,
    they can be applied in either order.
    """
    assert len(mat.shape) == 2
    assert mat.shape[0] == mat.shape[1]
    mat -= np.mean(mat, axis=0)[np.newaxis, :]
    mat -= np.mean(mat, axis=1)[:, np.newaxis]
    return mat


def gower_center(mat):
    """Apply row-column centering to a numpy square matrix in-place.

    Also returns the matrix.

    Note: row column centering and Gower centering commute, i.e.,
    they can be applied in either order.
    """
    assert len(mat.shape) == 2
    assert mat.shape[0] == mat.shape[1]
    N = mat.shape[0]
    # these are the same but denominator2 is faster
    # P = np.eye(N) - np.ones((N, N)) / N
    # denominator = np.trace(np.dot(np.dot(P, mat), P))
    denominator2 = np.trace(mat) - N*np.mean(mat)
    mat *= (N - 1) / denominator2
    return mat


def exact_arg_grm(arg, alpha=-1, from_pos = -1, to_pos = -1, diploid=True, centering=True):
    """
    Computes the exact ARG-GRM up to Gower centering and row-column centering.

    Arguments:
        arg: arg_needle_lib.ARG object
        alpha: if -1, the variants of branches are standardized before computing the GRM. If
          0, the variants of branches are treated as-is and not standardized. Values in between
          interpolate between these behaviors (default=-1).
        from_pos: start physical position for ARG-GRM computation (default=-1, use beginning of ARG)
        to_pos: end physical position for ARG-GRM computation (default=-1, use end of ARG)
        diploid: if True, leaves 2*i and 2*i+1 are treated as haploid pairs for sample i when
          computing the GRM (default=True)
        centering: if True, applies row-column centering and Gower centering to the GRM before
          returning (default=True)

    Returns:
        A numpy matrix with the exact ARG-GRM up to Gower centering and row-column centering.
    """
    def upper_diagonal_to_matrix(upper_diagonal, num_samples):
        arg_distance_matrix = np.zeros((num_samples, num_samples))
        for i, row in enumerate(upper_diagonal):
            for j, value in enumerate(row):
                arg_distance_matrix[i][i+j+1] = value
                arg_distance_matrix[i+j+1][i] = value
        return arg_distance_matrix

    logging.info("Computing exact ARG-GRM")
    upper_diagonal = arg_needle_lib_pybind.distance_matrix_v2(arg, alpha=alpha, from_pos=from_pos, to_pos=to_pos)
    exact_arg_grm = upper_diagonal_to_matrix(upper_diagonal, arg.num_samples())
    if diploid:
        exact_arg_grm = haploid_grm_to_diploid(exact_arg_grm)
    if centering:
        row_column_center(exact_arg_grm)
        gower_center(exact_arg_grm)
    return exact_arg_grm


def monte_carlo_arg_grm(arg, monte_carlo_mu, seed=1234, alpha=-1, diploid=True, centering=True, min_maf=0., max_maf=0.5, batch_size=256):
    """
    Computes the Monte Carlo ARG-GRM.

    This version wraps the C++ Eigen implementation and is more memory-efficient than
    the numpy version monte_carlo_arg_grm_numpy.

    Arguments:
        arg: arg_needle_lib.ARG object
        monte_carlo_mu: mutation rate used for sampling mutations to compute the ARG-GRM
        seed: random seed used for sampling (default=1234)
        alpha: if -1, all variants are standardized before computing the GRM. If 0, variants
          are treated as-is and not standardized. Values in between interpolate between these
          behaviors (default=-1).
        diploid: if True, leaves 2*i and 2*i+1 are treated as haploid pairs for sample i when
          computing the GRM (default=True)
        centering: if True, applies row-column centering and Gower centering to the GRM before
          returning (default=True)
        min_maf: exclusive minimum MAF (default=0.)
        max_maf: inclusive maximum MAF (default=0.5)
        batch_size: batch size when computing the ARG-GRM (default=256)

    Returns:
        A numpy matrix with the estimated Monte Carlo ARG-GRM.
    """
    logging.info(f"Sampling mutations using mu={monte_carlo_mu} and seed={seed}")
    arg_needle_lib_pybind.generate_mutations(arg, mu=monte_carlo_mu, random_seed=seed)
    logging.info("Computing Monte Carlo ARG-GRM using Eigen")
    grm = arg_needle_lib_pybind.compute_grm(
        arg, alpha=alpha, batch_size=batch_size, diploid=diploid, min_maf=min_maf, max_maf=max_maf)
    if centering:
        row_column_center(grm)
        gower_center(grm)
    return grm


def monte_carlo_arg_grm_numpy(arg, monte_carlo_mu, seed=1234, alpha=-1, diploid=True, centering=True, min_maf=0., max_maf=0.5, batch_size=256):
    """
    Computes the Monte Carlo ARG-GRM using numpy.

    This version is less memory-efficient, therefore the Eigen version monte_carlo_arg_grm
    is recommended instead.

    Arguments:
        arg: arg_needle_lib.ARG object
        monte_carlo_mu: mutation rate used for sampling mutations to compute the ARG-GRM
        seed: random seed used for sampling (default=1234)
        alpha: if -1, all variants are standardized before computing the GRM. If 0, variants
          are treated as-is and not standardized. Values in between interpolate between these
          behaviors (default=-1).
        diploid: if True, leaves 2*i and 2*i+1 are treated as haploid pairs for sample i when
          computing the GRM (default=True)
        centering: if True, applies row-column centering and Gower centering to the GRM before
          returning (default=True)
        min_maf: exclusive minimum MAF (default=0.)
        max_maf: inclusive maximum MAF (default=0.5)
        batch_size: batch size when computing the ARG-GRM (default=256)

    Returns:
        A numpy matrix with the estimated Monte Carlo ARG-GRM.
    """
    logging.info(f"Sampling mutations using mu={monte_carlo_mu} and seed={seed}")
    arg_needle_lib_pybind.generate_mutations(arg, mu=monte_carlo_mu, random_seed=seed)
    logging.info("Computing Monte Carlo ARG GRM using numpy")
    mutation_matrix = arg_needle_lib_pybind.get_mutations_matrix(arg)
    num_mutations = mutation_matrix.shape[0]

    num_samples = arg.num_samples()
    if diploid:
        assert num_samples % 2 == 0
        num_samples //= 2

    grm = np.zeros((num_samples, num_samples), dtype=np.float32)

    num_skipped = 0
    genotype_lists = []
    for i in range(num_mutations):
        geno = mutation_matrix[i, :]
        allele_count = int(np.sum(geno))
        af = allele_count / len(geno)
        maf = np.minimum(af, 1 - af)
        if maf == 0. or maf <= min_maf or maf > max_maf:
            num_skipped += 1
        else:
            genotype_lists.append(geno)
        if len(genotype_lists) == batch_size or i == num_mutations - 1:
            geno = np.array(genotype_lists, dtype=np.float32)
            if diploid:
                geno = geno.reshape((geno.shape[0], num_samples, 2)).sum(axis=-1)
            geno -= np.mean(geno, axis=1)[:, np.newaxis]
            # note: in very rare cases all variants are hets so std=0
            std = np.std(geno, axis=1, ddof=1)
            scales = std**alpha
            geno *= scales[:, np.newaxis]
            grm += np.dot(np.transpose(geno), geno)
            genotype_lists = []
    logging.info("{} out of {} mutations skipped due to MAF filtering".format(
        num_skipped, arg.num_mutations()))
    if centering:
        row_column_center(grm)
        gower_center(grm)
    return grm


def write_grm(grm, num_sites, path_prefix, write_binary=False):
    """
    Writes a GRM to file in a format that can be used by GCTA.

    Arguments:
        grm: numpy matrix containing the GRM
        num_sites: number of sites used to compute the GRM, expected but typically unused
          by downstream applications
        path_prefix: string, results will be written to path_prefix.grm.*
        write_binary: if True, write a binary file format that is observed to yield faster
          write times and can be read by GCTA. If False, write using gzip (default=False).
    """
    num_samples = grm.shape[0]
    with open(path_prefix + '.grm.id', 'w') as outfile:
        for i in range(num_samples):
            outfile.write('\t'.join([str(i+1), str(i+1)]) + '\n')

    if write_binary: # binary writing is much faster than gzip level 9 with similar storage
        grm_half_size = (num_samples * (num_samples + 1)) // 2
        out_grm = np.zeros(grm_half_size, dtype=np.dtype('f4')) # 'f4' means 'float32'
        current_start = 0
        for i in range(num_samples):
            out_grm[current_start:current_start+i+1] = grm[i][:(i+1)]
            current_start += i+1
        out_grm.tofile(path_prefix + '.grm.bin')
        out_num_sites = np.zeros(grm_half_size, dtype=np.dtype('f4')) # 'f4' means 'float32'
        out_num_sites += num_sites
        out_num_sites.tofile(path_prefix + '.grm.N.bin')
    else:
        with gzip.open(path_prefix + '.grm.gz', 'wb', compresslevel=6) as outfile:
            for i in range(num_samples):
                for j in range(i + 1):
                    outfile.write(('\t'.join([str(i+1), str(j+1), str(num_sites), str(grm[i][j])]) + '\n').encode())
