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

"""Tests for functions defined in src/arg_utils.cpp"""

import numpy as np
import pytest
import msprime
import math
import os
import tskit

import arg_needle_lib
from tree_metrics import kc_squared_distance, tmrca_mse

metrics_double = False


def test_arg_to_newick():
    arg = arg_needle_lib.ARG(0, 100, 3)
    arg.add_sample()
    arg.add_sample()
    arg.thread_sample([0, 40], [0, 0], [3.14, 2.718])
    arg.add_sample()
    arg.thread_sample([0, 50], [1, 0], [5, 1])
    arg.populate_children_and_roots()
    newick_strings = arg_needle_lib.arg_to_newick(arg, False).split("\n")
    assert len(newick_strings) == 4
    possible_strings = ["Tree in [0, 40): (2,(0,1)3)5", "Tree in [0, 40): (2,(1,0)3)5",
                        "Tree in [0, 40): ((0,1)3,2)5", "Tree in [0, 40): ((1,0)3,2)5"]
    assert newick_strings[0] in possible_strings
    possible_strings = ["Tree in [40, 50): (2,(0,1)4)5", "Tree in [40, 50): (2,(1,0)4)5",
                        "Tree in [40, 50): ((0,1)4,2)5", "Tree in [40, 50): ((1,0)4,2)5"]
    assert newick_strings[1] in possible_strings
    possible_strings = ["Tree in [50, 100): (1,(0,2)6)4", "Tree in [50, 100): (1,(2,0)6)4",
                        "Tree in [50, 100): ((0,2)6,1)4", "Tree in [50, 100): ((2,0)6,1)4"]
    assert newick_strings[2] in possible_strings
    assert newick_strings[3] == ""

    arg = arg_needle_lib.ARG(0, 50, [0, 0, 5], [True, True, False], [[0, 2], [1, 2]], [[0, 50], [0, 50]])
    arg.populate_children_and_roots()
    newick_strings = arg_needle_lib.arg_to_newick(arg, True).split("\n")
    assert len(newick_strings) == 2
    possible_strings = ["Tree in [0, 50): (1:5,0:5)2", "Tree in [0, 50): (0:5,1:5)2"]
    assert newick_strings[0] in possible_strings
    assert newick_strings[1] == ""


def test_metrics_small():
    """TMRCA MSE, KC topology, and mutation match
    """
    arg1 = arg_needle_lib.ARG(0, 100, 4)
    arg1.add_sample()
    arg1.add_sample()
    arg1.thread_sample([0], [0], [10])
    arg1.add_sample()
    arg1.thread_sample([0], [1], [7])
    arg1.add_sample()
    arg1.thread_sample([0, 40], [2, 2], [3, 9])
    arg1.populate_children_and_roots()

    arg2 = arg_needle_lib.ARG(0, 100, 4)
    arg2.add_sample()
    arg2.add_sample()
    arg2.thread_sample([0], [0], [10])
    arg2.add_sample()
    arg2.thread_sample([0], [1], [7])
    arg2.add_sample()
    arg2.thread_sample([0, 40], [2, 2], [8, 6])
    arg2.populate_children_and_roots()

    ts1 = arg_needle_lib.arg_to_tskit(arg1)
    ts2 = arg_needle_lib.arg_to_tskit(arg2)

    # First calculation uses ARG-based TMRCA MSE
    tmrca_mse1 = arg_needle_lib.tmrca_mse(arg1, arg2)
    # Second calculation aggregates tree-based TMRCA MSE
    x = 0.4
    y = 0.6
    tmrca_mse2 = 0
    tmrca_mse2 += x*tmrca_mse(ts1.at(0), ts2.at(0))
    tmrca_mse2 += y*tmrca_mse(ts1.at(40), ts2.at(40))
    # Third calculation does it manually
    tmrca_mse3 = 0
    # Terms from TMRCA(0, 3)
    tmrca_mse3 += x*(10-10)**2 + y*(10-10)**2
    # Terms from TMRCA(1, 3)
    tmrca_mse3 += x*(7-8)**2 + y*(9-7)**2
    # Terms from TMRCA(2, 3)
    tmrca_mse3 += x*(3-8)**2 + y*(9-6)**2
    tmrca_mse3 /= 6 # number of pairs
    if metrics_double:
        assert tmrca_mse1 == tmrca_mse2 == tmrca_mse3
    else:
        assert math.isclose(tmrca_mse1, tmrca_mse2, rel_tol=1e-6)
        assert math.isclose(tmrca_mse1, tmrca_mse3, rel_tol=1e-6)
        assert math.isclose(tmrca_mse2, tmrca_mse3, rel_tol=1e-6)

    # First calculation uses ARG-based TMRCA MSE
    kc1 = arg_needle_lib.kc_topology(arg1, arg2)
    # Second calculation aggregates tree-based TMRCA MSE
    kc2 = 0
    kc2 += x*kc_squared_distance(ts1.at(0), ts2.at(0), kc_lambdas=[0])[0]
    kc2 += y*kc_squared_distance(ts1.at(40), ts2.at(40), kc_lambdas=[0])[0]
    # Third calculation does it manually
    kc3 = x*((2-1)**2 + (1-2)**2) + y*((1-2)**2 + (2-1)**2)
    assert kc1 == kc2 == kc3

    # Mutation matching on ARG 1
    assert arg_needle_lib.mutation_match(arg1, 20, [0, 0, 0, 0])
    assert arg_needle_lib.mutation_match(arg1, 20, [0, 0, 0, 1])
    assert arg_needle_lib.mutation_match(arg1, 20, [0, 0, 1, 1])
    assert not arg_needle_lib.mutation_match(arg1, 20, [0, 1, 0, 1])
    assert not arg_needle_lib.mutation_match(arg1, 20, [0, 1, 1, 0])
    assert arg_needle_lib.mutation_match(arg1, 20, [0, 1, 1, 1])
    assert not arg_needle_lib.mutation_match(arg1, 20, [1, 0, 0, 1])

    assert arg_needle_lib.mutation_match(arg1, 50, [0, 0, 0, 0])
    assert arg_needle_lib.mutation_match(arg1, 50, [0, 0, 0, 1])
    assert not arg_needle_lib.mutation_match(arg1, 50, [0, 0, 1, 1])
    assert not arg_needle_lib.mutation_match(arg1, 50, [0, 1, 0, 1])
    assert arg_needle_lib.mutation_match(arg1, 50, [0, 1, 1, 0])
    assert arg_needle_lib.mutation_match(arg1, 50, [0, 1, 1, 1])
    assert arg_needle_lib.mutation_match(arg1, 50, [1, 0, 0, 1])


def test_metrics():
    """TMRCA MSE and KC topology
    """
    sample_size = 20  # also works for 100
    ts1 = msprime.simulate(
        sample_size=sample_size, Ne=1e4, length=5e3, recombination_rate=2e-8,
        mutation_rate=2e-8, random_seed=1)
    ts2 = msprime.simulate(
        sample_size=sample_size, Ne=1e4, length=5e3, recombination_rate=2e-8,
        mutation_rate=2e-8, random_seed=2)
    arg1 = arg_needle_lib.tskit_to_arg(ts1)
    arg2 = arg_needle_lib.tskit_to_arg(ts2)

    bp1 = np.array(list(ts1.breakpoints()))
    bp2 = np.array(list(ts2.breakpoints()))
    bp_common = np.unique(np.concatenate((bp1, bp2)))
    num_regions = len(bp_common) - 1
    lengths = np.diff(bp_common)
    result = np.zeros((num_regions, 2))
    for k in range(num_regions):
        pos = bp_common[k]
        tree1 = ts1.at(pos)
        tree2 = ts2.at(pos)

        result[k, 0] = tmrca_mse(tree1, tree2)
        result[k, 1] = kc_squared_distance(tree1, tree2, [0])[0]

    tmrca_trees = np.sum(result[:, 0] * lengths) / np.sum(lengths)
    kc_trees = np.sum(result[:, 1] * lengths) / np.sum(lengths)
    tmrca_needle = arg_needle_lib.tmrca_mse(arg1, arg2)
    kc_needle = arg_needle_lib.kc_topology(arg1, arg2)
    kc_stab, tmrca_stab = arg_needle_lib.metrics_stab(arg1, arg2, num_stabs=5000)
    arg1.populate_children_and_roots()
    arg2.populate_children_and_roots()
    kc_stab_efficient, tmrca_stab_efficient, _ = arg_needle_lib.metrics_stab_efficient(
        arg1, arg2, num_stabs=5000)
    kc_stab_efficient2 = arg_needle_lib.kc2_length_stab_efficient(
        arg1, arg2, num_stabs=5000, lambdas=[0])[0]

    print(tmrca_trees, tmrca_needle, kc_trees, kc_needle)
    print(tmrca_stab, tmrca_stab_efficient, kc_stab, kc_stab_efficient, kc_stab_efficient2)
    if metrics_double:
        assert abs(tmrca_trees - tmrca_needle) < 1e-3
        assert abs(kc_trees - kc_needle) < 1e-5
        assert abs(tmrca_stab - tmrca_stab_efficient) < 1e-2
        assert abs(kc_stab - kc_stab_efficient) < 1e-5
        assert abs(kc_stab_efficient - kc_stab_efficient2) < 1e-5
    else:
        assert math.isclose(tmrca_trees, tmrca_needle, rel_tol=1e-5)
        assert math.isclose(kc_trees, kc_needle, rel_tol=1e-6)
        assert math.isclose(tmrca_stab, tmrca_stab_efficient, rel_tol=1e-3)
        assert math.isclose(kc_stab, kc_stab_efficient, rel_tol=1e-6)
        assert math.isclose(kc_stab_efficient, kc_stab_efficient2, rel_tol=1e-6)


def test_kc_merge_play():
    sample_size = 20  # also works for 100
    ts1 = msprime.simulate(
        sample_size=sample_size, Ne=1e4, length=5e3, recombination_rate=2e-8,
        mutation_rate=2e-8, random_seed=5)
    ts2 = msprime.simulate(
        sample_size=sample_size, Ne=1e4, length=5e3, recombination_rate=2e-8,
        mutation_rate=2e-8, random_seed=6)
    arg1 = arg_needle_lib.tskit_to_arg(ts1)
    arg2 = arg_needle_lib.tskit_to_arg(ts2)

    arg1.populate_children_and_roots()
    arg2.populate_children_and_roots()

    fractions = [0, 0.1, 0.2, 0.5, 1]
    results = []
    for f in fractions:
        kc_stab_efficient, _, _ = arg_needle_lib.metrics_stab_efficient(
            arg1, arg2, num_stabs=1000, merge_type=2, merge_fraction=f)
        results.append(kc_stab_efficient)

    kc_stab_efficient, _, _ = arg_needle_lib.metrics_stab_efficient(
        arg1, arg2, num_stabs=1000)

    for f, result in zip(fractions, results):
        print("Fraction {}, KC = {}".format(f, result))

    print("Regular, KC = {}".format(kc_stab_efficient))


@pytest.mark.parametrize("seed", [10, 11, 12])
@pytest.mark.parametrize("sample_size", [100])
@pytest.mark.parametrize("length", [5e3])
def test_medium_visit_identical(seed, sample_size, length):
    """Test that fast and slow visits give the same volume map"""
    ts = msprime.simulate(
        sample_size=sample_size, Ne=4e4, length=length, recombination_rate=2e-8,
        mutation_rate=2e-8, random_seed=seed)

    arg = arg_needle_lib.tskit_to_arg(ts)
    arg.populate_children_and_roots()
    print(ts.num_trees, "trees,", ts.num_nodes, "nodes", end=" ", flush=True)
    if metrics_double:
        assert arg_needle_lib.visit_identical(arg, abs_tol=1e-9, timing=False, verbose=False)
    else:
        assert arg_needle_lib.visit_identical(arg, rel_tol=1e-6, timing=False, verbose=False)


@pytest.mark.parametrize("seed", [10, 11, 12])
@pytest.mark.parametrize("sample_size", [100])
@pytest.mark.parametrize("length", [5e3])
def test_medium_precision_recall(seed, sample_size, length):
    """Test that fast and slow visits give the same volume map"""
    ts = msprime.simulate(
        sample_size=sample_size, Ne=4e4, length=length, recombination_rate=2e-8,
        mutation_rate=2e-8, random_seed=seed)

    arg = arg_needle_lib.tskit_to_arg(ts)
    arg.populate_children_and_roots()
    result = arg_needle_lib.bitset_overlap_stab(arg, arg, 5000)
    assert result[0] == result[1] == result[2];
    assert result[3] == result[4] == result[5];


def test_distance_matrix():
    arg1 = arg_needle_lib.ARG(0, 3, 4)
    arg1.add_sample()
    arg1.add_sample()
    arg1.thread_sample([0, 2], [0, 0], [2, 3])
    arg1.add_sample()
    arg1.thread_sample([0, 1, 2], [1, 1, 1], [4, 3, 1])
    arg1.add_sample()
    arg1.thread_sample([0, 1], [2, 2], [2, 4])

    arg_needle_lib.DescendantList.set_threshold(int(5e4/64))

    seed=10
    length=5e3
    sample_size=500
    ts = msprime.simulate(
        sample_size=sample_size, Ne=4e4, length=length, recombination_rate=2e-8,
        mutation_rate=0., random_seed=seed)
    print(ts.num_trees, "trees,", ts.num_nodes, "nodes", flush=True)
    arg2 = arg_needle_lib.tskit_to_arg(ts)

    def upper_diagonal_to_matrix(upper_diagonal, num_samples):
        arg_distance_matrix = np.zeros((num_samples, num_samples))
        for i, row in enumerate(upper_diagonal):
            for j, value in enumerate(row):
                arg_distance_matrix[i][i+j+1] = value
                arg_distance_matrix[i+j+1][i] = value
        return arg_distance_matrix

    for arg in [arg1, arg2]:
        arg.populate_children_and_roots()
        foo = arg_needle_lib.distance_matrix(arg)
        bar1 = upper_diagonal_to_matrix(foo, len(foo) + 1)
        print("Done with distance matrix v1")

        foo = arg_needle_lib.distance_matrix_v2(arg, alpha=0)
        bar2 = upper_diagonal_to_matrix(foo, len(foo) + 1)
        print("Done with distance matrix v2")

        assert np.allclose(bar1, bar2)

def test_tsinfer_visit():
    """Test visit_branches and visit_clades on a tsinfer ARG.
    """
    current_dir = os.path.dirname(os.path.abspath(__file__))
    input_path = os.path.join(current_dir, "../testdata/tsinfer_snp_length_5e6_samples_3e2/1.trees")
    ts = tskit.load(input_path)
    arg = arg_needle_lib.tskit_to_arg(ts)
    arg.populate_children_and_roots()
    num1, num2, num_common = arg_needle_lib.bitset_overlap_full(arg, arg, -1., -1.)
    assert num1 == num2 == num_common
    predicted_volume = arg_needle_lib.total_volume(arg)
    mu=5e-8
    ts_new = msprime.mutate(ts, rate=mu, random_seed=2)
    volume_map = arg_needle_lib.bitset_volume_map(arg)
    volume = 0
    for key in volume_map:
        volume += volume_map[key]

    mutation_map = arg_needle_lib.generate_mutations_map(arg, mu)

    print(predicted_volume*mu, volume*mu, ts_new.num_mutations, len(mutation_map))

    assert arg_needle_lib.visit_identical(arg, abs_tol=1e-9, timing=False, verbose=False)
    # number of clades has the root removed, so does number of branches
    # visit_clades does visit the root, but our function removes it
    assert num1 == len(volume_map)


def test_regional_visit():
    """Test that we obtain regional bitsets that are the same.
    """
    current_dir = os.path.dirname(os.path.abspath(__file__))
    input_path = os.path.join(current_dir, "../testdata/tsinfer_snp_length_5e6_samples_3e2/1.trees")
    ts = tskit.load(input_path)
    arg = arg_needle_lib.tskit_to_arg(ts)
    arg.populate_children_and_roots()


    test_regions = [(-1,-1), (-1, 1e5), (2e4,1e5)]
    for (i,j) in test_regions:
        num1, num2, num_common = arg_needle_lib.bitset_overlap_full(arg, arg, i,j)
        assert num1 == num2 == num_common

    # Test purposely bad conditions (and catch the error)
    bad_regions = [(10,1)]
    for (i,j) in bad_regions:
        with pytest.raises(RuntimeError):
            num1, num2, num_common = arg_needle_lib.bitset_overlap_full(arg, arg, i,j)

def test_regional_volume():
    """Test that the sum of regional volumes is equal to the total volume
    """
    seed = 42
    ts = msprime.simulate(Ne=1e4, sample_size=100, recombination_rate=5e-8, length=5e6, random_seed=seed)
    arg = arg_needle_lib.tskit_to_arg(ts)
    arg.populate_children_and_roots()

    num_regions = 3
    total_volume = arg_needle_lib.total_volume(arg)
    regional_volume = 0.
    # NOTE: we add one to make sure we have num_regions chunks
    test_regions = np.linspace(0, ts.sequence_length, num_regions+1)
    for i in range(num_regions):
        regional_volume += arg_needle_lib.local_volume(arg, test_regions[i], test_regions[i+1])
    # Check that the sum across regions are very close to one another.
    assert math.isclose(total_volume, regional_volume, rel_tol=1e-6)
