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

import platform
import pytest
import msprime
import numpy as np
import arg_needle_lib

@pytest.mark.skipif(platform.system() != "Linux", reason="msprime.simulate depends on platform")
def test_kc2_tmrca_sv_stab():
    """Test that mutations are properly populated onto edges"""
    seed = 130222
    num_pairs = 6
    num_stabs = 1

    ts1 = msprime.simulate(sample_size=4, random_seed=seed)
    tree1 = ts1.first()
    # Generates the following tree
    # 1.23┊    6    ┊
    #     ┊  ┏━┻━┓  ┊
    # 0.85┊  ┃   5  ┊
    #     ┊  ┃  ┏┻┓ ┊
    # 0.17┊  4  ┃ ┃ ┊
    #     ┊ ┏┻┓ ┃ ┃ ┊
    # 0.00┊ 0 2 1 3 ┊
    #     0         1
    arg1 = arg_needle_lib.tskit_to_arg(ts1)
    arg1.populate_children_and_roots()
    kc_vector_1 = np.array([0, 1, 0, 0, 1, 0, 1, 1, 1, 1])
    split_vector_1 = np.array([4, 2, 4, 4, 2, 4])
    tmrca_vector_1 = np.array([tree1.time(tree1.mrca(i, j)) for i in range(4) for j in range(i + 1, 4)])

    ts2 = msprime.simulate(sample_size=4, random_seed=seed + 1)
    tree2 = ts2.first()
    # Generates the following tree
    # 2.65┊     6   ┊
    #     ┊   ┏━┻━┓ ┊
    # 0.71┊   5   ┃ ┊
    #     ┊  ┏┻━┓ ┃ ┊
    # 0.51┊  4  ┃ ┃ ┊
    #     ┊ ┏┻┓ ┃ ┃ ┊
    # 0.00┊ 0 2 1 3 ┊
    #     0         1
    arg2 = arg_needle_lib.tskit_to_arg(ts2)
    arg2.populate_children_and_roots()
    kc_vector_2 = np.array([1, 2, 0, 1, 0, 0, 1, 1, 1, 1])
    split_vector_2 = np.array([3, 2, 4, 3, 4, 4])
    tmrca_vector_2 = np.array([tree2.time(tree2.mrca(i, j)) for i in range(4) for j in range(i + 1, 4)])

    kc_squared = ((kc_vector_1 - kc_vector_2)**2).sum()
    sv_squared = ((split_vector_1 - split_vector_2)**2).sum()
    tmrca_squared = ((tmrca_vector_1 - tmrca_vector_2)**2).sum()

    metrics = arg_needle_lib.kc2_tmrca_sv_stab(arg1, arg2, num_stabs=num_stabs)

    assert metrics['kc2'] == kc_squared
    assert metrics['tmrca_mse'] == tmrca_squared / (num_pairs * num_stabs)
    assert metrics['sv2'] == sv_squared / (num_pairs * num_stabs)
