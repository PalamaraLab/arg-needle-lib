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


"""Computationally intensive tests.

These are not run by default but can be separately run if chosen (see
repository README.md).

Takes ~44 mins to run on a MacBook Air laptop.
"""

import msprime
import pytest

from arg_needle_lib import arg_to_tskit, tskit_to_arg, tskit_to_arg_thread
from tree_metrics import r2_matrix, compare_r2, robinson_foulds, kc_squared_distance


# Take 66 seconds
@pytest.mark.parametrize("seed", [10, 11, 12])
@pytest.mark.parametrize("converter", [tskit_to_arg, tskit_to_arg_thread])
@pytest.mark.parametrize("sample_size", [400])
@pytest.mark.parametrize("length", [5e3])
def test_medium_arg_checks(seed, sample_size, length, converter):
    """Convert to ARG and perform checks"""
    ts = msprime.simulate(
        sample_size=sample_size, Ne=4e4, length=length, recombination_rate=2e-8,
        mutation_rate=2e-8, random_seed=seed)

    arg = converter(ts)
    print(ts.num_trees, "trees,", ts.num_nodes, "->", len(arg.node_ids()), "nodes", end=" ", flush=True)

    if converter == tskit_to_arg:
        # With this converter, nodes always span the whole ARG, so we
        # need to disable stringent checking
        stringent = False
    else:
        stringent = True

    arg.check_basic(stringent)
    arg.populate_children_and_roots()
    arg.check_roots()
    arg.check_children(stringent)


# Takes 60 seconds for 3 seeds
# @pytest.mark.parametrize("seed", [10, 11, 12])
@pytest.mark.parametrize("seed", [11, 12])
@pytest.mark.parametrize("converter", [tskit_to_arg])
@pytest.mark.parametrize("sample_size", [1000, 10000])
@pytest.mark.parametrize("length", [5e6])
def test_big_arg_checks(seed, sample_size, length, converter):
    """Convert to ARG and perform checks"""
    ts = msprime.simulate(
        sample_size=sample_size, Ne=4e4, length=length, recombination_rate=2e-8,
        mutation_rate=2e-8, random_seed=seed)

    arg = converter(ts)
    print(ts.num_trees, "trees,", ts.num_nodes, "->", len(arg.node_ids()), "nodes", end=" ", flush=True)

    if converter == tskit_to_arg:
        # With this converter, nodes always span the whole ARG, so we
        # need to disable stringent checking
        stringent = False
    else:
        stringent = True

    arg.check_basic(stringent)
    arg.populate_children_and_roots()
    arg.check_roots()
    arg.check_children(stringent)


# Takes a long time, especially threading
@pytest.mark.parametrize("seed", [10, 11, 12])
@pytest.mark.parametrize("converter", [tskit_to_arg, tskit_to_arg_thread])
@pytest.mark.parametrize("sample_size", [100, 400])
@pytest.mark.parametrize("length", [5e3])
def test_medium_roundtrip(seed, sample_size, length, converter):
    """Roundtrip = convert to ARG then back and check metrics"""
    ts = msprime.simulate(
        sample_size=sample_size, Ne=4e4, length=length, recombination_rate=2e-8,
        mutation_rate=2e-8, random_seed=seed)

    arg = converter(ts)
    ts2 = arg_to_tskit(arg)
    print(ts.num_trees, "trees,", ts.num_nodes, "->", len(arg.node_ids()), "nodes", end=" ", flush=True)

    assert ts.num_trees == ts2.num_trees
    assert len(arg.node_ids()) == ts2.num_nodes
    for t1, t2 in zip(ts.trees(), ts2.trees()):
        r2_mat = r2_matrix(t1, t2)
        assert compare_r2(r2_mat) == (1.0, 1.0)
        assert robinson_foulds(r2_mat) == 0
        assert kc_squared_distance(t1, t2, [0.0, 0.5, 1.0]) == [0, 0, 0]
