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


import pytest
import shutil

import tskit
import tszip
import msprime
import numpy as np

import arg_needle_lib
from tree_metrics import r2_matrix, compare_r2, robinson_foulds, kc_squared_distance, tmrca_mse

metrics_double = False


@pytest.mark.parametrize("seed", [2, 3, 5, 7, 12345678])
@pytest.mark.parametrize("converter", [arg_needle_lib.tskit_to_arg, arg_needle_lib.tskit_to_arg_thread])
def test_arg_checks(seed, converter):
    """Convert to ARG and perform checks"""
    ts = msprime.simulate(
        sample_size=20, Ne=1e4, length=5e3, recombination_rate=2e-8,
        mutation_rate=2e-8, random_seed=seed)

    arg = converter(ts)
    print(ts.num_trees, "trees,", ts.num_nodes, "->", len(arg.node_ids()), "nodes", end=" ", flush=True)

    if converter == arg_needle_lib.tskit_to_arg:
        # With this converter, nodes always span the whole ARG, so we
        # need to disable stringent checking
        stringent = False
    else:
        stringent = True

    arg.check_basic(stringent)
    arg.populate_children_and_roots()
    arg.check_roots()
    arg.check_children(stringent)


@pytest.mark.parametrize("seed", [2, 3, 5, 7, 12345678])
@pytest.mark.parametrize("converter", [arg_needle_lib.tskit_to_arg, arg_needle_lib.tskit_to_arg_thread])
def test_roundtrip(seed, converter):
    """Roundtrip = convert to ARG then back and check metrics"""
    ts = msprime.simulate(
        sample_size=20, Ne=1e4, length=5e3, recombination_rate=2e-8,
        mutation_rate=2e-8, random_seed=seed)

    # https://stackoverflow.com/a/3071
    arg = converter(ts)
    ts2 = arg_needle_lib.arg_to_tskit(arg)
    print(ts.num_trees, "trees,", ts.num_nodes, "->", len(arg.node_ids()), "nodes", end=" ", flush=True)

    assert ts.num_trees == ts2.num_trees
    assert len(arg.node_ids()) == ts2.num_nodes
    for t1, t2 in zip(ts.trees(), ts2.trees()):
        r2_mat = r2_matrix(t1, t2)
        assert compare_r2(r2_mat) == (1.0, 1.0)
        assert robinson_foulds(r2_mat) == 0

        if metrics_double:
            assert kc_squared_distance(t1, t2, [0.0, 0.5, 1.0]) == [0, 0, 0]
            assert tmrca_mse(t1, t2) == 0
        else:
            # The breakpoints may shift in this case
            assert kc_squared_distance(t1, t2, [0.0])[0] == 0
            assert kc_squared_distance(t1, t2, [0.5])[0] < 1e-3
            assert kc_squared_distance(t1, t2, [1.0])[0] < 1e-3
            assert tmrca_mse(t1, t2) < 1e-4


@pytest.mark.parametrize("seed", [2, 3, 5, 7, 12345678])
@pytest.mark.parametrize("permutation_seed", [314159, 271828])
def test_arg_to_tskit_permutation(seed, permutation_seed):
    """Start with ts, convert to ARG then back with permutation, and check pairwise TMRCA"""
    sample_size = 20
    ts = msprime.simulate(
        sample_size=sample_size, Ne=1e4, length=5e3, recombination_rate=2e-8,
        mutation_rate=2e-8, random_seed=seed)

    permutation_rng = np.random.default_rng(permutation_seed)
    permutation = permutation_rng.permutation(sample_size)
    arg = arg_needle_lib.tskit_to_arg(ts)
    ts2 = arg_needle_lib.arg_to_tskit(arg, sample_permutation=permutation)

    # Calculate TMRCAs of a few pairs and compare
    # Because we swapped the samples, we should get different results, but we can fix this
    # by applying the permutation to the pairs when looking up the original trees
    tmrcas1 = []
    tmrcas2 = []
    tmrcas3 = []
    positions = [5e3 // 4, 5e3 // 2, 5e3 * 3 // 4]
    pairs = [[2*i, 2*i+1] for i in range(sample_size // 2)]
    permuted_pairs = [[permutation[2*i], permutation[2*i+1]] for i in range(sample_size // 2)]

    for position in positions:
        for pair, permuted_pair in zip(pairs, permuted_pairs):
            tree1 = ts.at(position)
            tree2 = ts2.at(position)
            tmrcas1.append(tree2.tmrca(pair[0], pair[1])) # pairs in converted tree
            tmrcas2.append(tree1.tmrca(pair[0], pair[1])) # pairs in original tree
            tmrcas3.append(tree1.tmrca(permuted_pair[0], permuted_pair[1])) # permuted pairs in original tree

    assert tmrcas1 != tmrcas2
    assert tmrcas1 == tmrcas3
    assert tmrcas2 != tmrcas3


@pytest.mark.parametrize("seed", [2, 3, 5])
def test_arg_to_tskit_invalid_permutation(seed):
    """Check that invalid permutations yield an error"""
    sample_size = 20
    ts = msprime.simulate(
        sample_size=sample_size, Ne=1e4, length=5e3, recombination_rate=2e-8,
        mutation_rate=2e-8, random_seed=seed)
    arg = arg_needle_lib.tskit_to_arg(ts)

    identity_permutation = list(range(sample_size))

    permutation = identity_permutation.copy()
    permutation[0] = 1
    with pytest.raises(ValueError):
        ts2 = arg_needle_lib.arg_to_tskit(arg, sample_permutation=permutation)

    permutation = [x + 1 for x in identity_permutation]
    with pytest.raises(ValueError):
        ts2 = arg_needle_lib.arg_to_tskit(arg, sample_permutation=permutation)


def test_roundtrip_offset(tmp_path):
    """Roundtrip = convert to ARG then back. Check if offset is preserved."""
    ts = msprime.simulate(
        sample_size=20, Ne=1e4, length=5e3, recombination_rate=2e-8,
        mutation_rate=2e-8, random_seed=4242)

    arg = arg_needle_lib.tskit_to_arg(ts)
    assert arg.offset == 0

    for bad_offset in [-1, -100]:
        with pytest.raises(RuntimeError):
            # set_offset expects a nonnegative integer
            arg.set_offset(bad_offset)

    for bad_offset in [0.5, 3.14159]:
        with pytest.raises(TypeError):
            # set_offset expects a nonnegative integer
            arg.set_offset(bad_offset)

    for offset in [0, 5, 123456789]:
        arg.set_offset(offset)
        assert arg.offset == offset

        ts2 = arg_needle_lib.arg_to_tskit(arg)
        arg2 = arg_needle_lib.tskit_to_arg(ts2)
        assert arg2.offset == offset

    # Check offset is preserved through ts.dump()
    dump_path = tmp_path / "dump.trees"
    arg3 = arg_needle_lib.tskit_to_arg(ts)
    arg3.set_offset(3)
    ts3 = arg_needle_lib.arg_to_tskit(arg3)
    ts3.dump(dump_path)
    ts4 = tskit.load(dump_path)
    arg4 = arg_needle_lib.tskit_to_arg(ts4)
    assert arg4.offset == 3

    # Check offset is preserved through tszip.compress()
    compress_path = tmp_path / "compress.tsz"
    arg5 = arg_needle_lib.tskit_to_arg(ts)
    arg5.set_offset(4)
    ts5 = arg_needle_lib.arg_to_tskit(arg5)
    tszip.compress(ts5, compress_path)
    ts6 = tszip.decompress(compress_path)
    arg6 = arg_needle_lib.tskit_to_arg(ts6)
    assert arg6.offset == 4

    # Pytest keeps these around longer than we need so we manually delete
    shutil.rmtree(tmp_path)


def test_roundtrip_chromosome(tmp_path):
    """Roundtrip = convert to ARG then back. Check if chromosome number is preserved."""
    ts = msprime.simulate(
        sample_size=20, Ne=1e4, length=5e3, recombination_rate=2e-8,
        mutation_rate=2e-8, random_seed=4242)

    # Test default chromosome is set to 1
    arg = arg_needle_lib.tskit_to_arg(ts)
    assert arg.chromosome == 1

    for bad_chromosome in [0, -100]:
        with pytest.raises(RuntimeError):
            # set_chromosome expects a positive integer
            arg.set_chromosome(bad_chromosome)

    for bad_chromosome in [0.5, 3.14159]:
        with pytest.raises(TypeError):
            # set_chromosome expects a positive integer
            arg.set_chromosome(bad_chromosome)

    for chromosome in [1, 5, 123456789]:
        arg.set_chromosome(chromosome)
        assert arg.chromosome == chromosome

        ts2 = arg_needle_lib.arg_to_tskit(arg)
        arg2 = arg_needle_lib.tskit_to_arg(ts2)
        assert arg2.chromosome == chromosome

    # Check chromosome is preserved through ts.dump()
    dump_path = tmp_path / "dump.trees"
    arg3 = arg_needle_lib.tskit_to_arg(ts)
    arg3.set_chromosome(3)
    ts3 = arg_needle_lib.arg_to_tskit(arg3)
    ts3.dump(dump_path)
    ts4 = tskit.load(dump_path)
    arg4 = arg_needle_lib.tskit_to_arg(ts4)
    assert arg4.chromosome == 3

    # Check chromosome is preserved through tszip.compress()
    compress_path = tmp_path / "compress.tsz"
    arg5 = arg_needle_lib.tskit_to_arg(ts)
    arg5.set_chromosome(4)
    ts5 = arg_needle_lib.arg_to_tskit(arg5)
    tszip.compress(ts5, compress_path)
    ts6 = tszip.decompress(compress_path)
    arg6 = arg_needle_lib.tskit_to_arg(ts6)
    assert arg6.chromosome == 4

    # Pytest keeps these around longer than we need so we manually delete
    shutil.rmtree(tmp_path)

@pytest.mark.parametrize("seed", [2, 123456])
@pytest.mark.parametrize("m", [1, 3, 10, 1234])
def test_exact_mutations(seed, m):
    """Simulate exactly m mutations on the arg
       and check that it passes conversion.
    """
    ts = msprime.simulate(
        sample_size=20, Ne=1e4, length=5e3, recombination_rate=2e-8, random_seed=seed)
    arg = arg_needle_lib.tskit_to_arg(ts)
    # Simulate exactly M mutations on the ARG 
    arg_needle_lib.generate_m_mutations(arg, m, random_seed=seed)
    assert arg.num_mutations() == m
    
    # Convert back to tskit (with mutations)
    cur_ts = arg_needle_lib.arg_to_tskit(arg, mutations=True)
    assert cur_ts.num_mutations == m

def test_convert_fail():
    """Conversion from ARG to tskit should fail if:
    1. Start is nonzero
    2. Nodes are not all consecutive
    """
    arg = arg_needle_lib.ARG(20, 100, 3)
    arg.add_sample()
    arg.add_sample()
    arg.thread_sample([20, 40], [0, 0], [3.14, 2.718])
    arg.add_sample()
    arg.thread_sample([20, 50], [1, 0], [5, 1])

    with pytest.raises(RuntimeError):
        # Can only call set_offset when ARG start is zero
        arg.set_offset(5)

    with pytest.raises(ValueError):
        # Can only convert to tskit when ARG start is zero
        ts = arg_needle_lib.arg_to_tskit(arg)

    arg = arg_needle_lib.ARG(0, 100, 3)
    arg.add_sample()
    arg.add_sample()
    arg.thread_sample([0, 40], [0, 0], [3.14, 2.718])
    with pytest.raises(ValueError):
        # At this point the nodes taken should be 0, 1, 3, and 4
        # There is a gap in node IDs which causes an error
        ts = arg_needle_lib.arg_to_tskit(arg)
