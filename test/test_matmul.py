# This file is part of the ARG-Needle genealogical inference and
# analysis software suite.
# Copyright (C) 2023-2025 ARG-Needle Developers.

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

import arg_needle_lib
import numpy as np
import math
import msprime

def test_matmul():
    seed = 35
    alpha = -1

    generation_time = 25
    T_OOA = 21.2e3 / generation_time
    r_CEU = 0.004
    N_CEU = 1000 / math.exp(-r_CEU * T_OOA)
    demography = msprime.Demography()
    demography.add_population(name="CEU", initial_size=N_CEU, growth_rate=r_CEU)
    n_samples = 500
    ts = msprime.sim_ancestry(
        demography=demography, samples=n_samples, recombination_rate=1e-8, sequence_length=1e7, random_seed=seed
    )

    arg = arg_needle_lib.tskit_to_arg(ts)
    arg.populate_children_and_roots()
    arg_needle_lib.generate_mutations(arg, 1e-8, seed)
    arg.populate_mutations_on_edges()
    arg_needle_lib.prepare_matmul(arg)

    np.random.seed(seed)
    n_side_in_mat = np.random.normal(0, 1, size=(20, n_samples))
    p_side_in_mat = np.random.normal(0, 1, size=(arg.num_mutations(), 20))

    hap_genotype = arg_needle_lib.get_mutations_matrix(arg).astype(np.float64)
    dip_genotype = (hap_genotype[:, ::2] + hap_genotype[:, 1::2]).T
    dip_genotype *= np.power((1 - dip_genotype.mean(axis=0, keepdims=True)/2) * dip_genotype.mean(axis=0, keepdims=True), 0.5*alpha)
    dip_genotype -= dip_genotype.mean(axis=0, keepdims=True)

    # Single and multi threaded variants of mul on mutations
    expected_res_n_side = n_side_in_mat @ dip_genotype
    arg_res_mt_n_side = arg_needle_lib.arg_matmul(arg, n_side_in_mat, axis="mutations", standardize=True, alpha=alpha, diploid=True)
    assert np.allclose(arg_res_mt_n_side, expected_res_n_side)
    arg_res_mt_n_side_mt = arg_needle_lib.arg_matmul(arg, n_side_in_mat, axis="mutations", standardize=True, alpha=alpha, diploid=True, n_threads=4)
    assert np.allclose(arg_res_mt_n_side_mt, expected_res_n_side)

    # Single and multi threaded variants of mul on samples
    expected_res_p_side = dip_genotype @ p_side_in_mat
    arg_res_mt_p_side = arg_needle_lib.arg_matmul(arg, p_side_in_mat, axis="samples", standardize=True, alpha=alpha, diploid=True)
    assert np.allclose(arg_res_mt_p_side, expected_res_p_side)
    arg_res_mt_p_side_mt = arg_needle_lib.arg_matmul(arg, p_side_in_mat, axis="samples", standardize=True, alpha=alpha, diploid=True, n_threads=4)
    assert np.allclose(arg_res_mt_p_side_mt, expected_res_p_side)
