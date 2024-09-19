# This file is part of the ARG-Needle genealogical inference and
# analysis software suite.
# Copyright (C) 2023-2024 ARG-Needle Developers.

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
import msprime
import arg_needle_lib


def test_map_genotypes_small():
    """Test that input genotypes are correctly converted to mutations on edges"""
    seed = 130222
    ts = msprime.simulate(sample_size=4, random_seed=seed)
    # Generates the following tree
    # 1.23┊    6    ┊
    #     ┊  ┏━┻━┓  ┊
    # 0.85┊  ┃   5  ┊
    #     ┊  ┃  ┏┻┓ ┊
    # 0.17┊  4  ┃ ┃ ┊
    #     ┊ ┏┻┓ ┃ ┃ ┊
    # 0.00┊ 0 2 1 3 ┊
    #     0         1
    arg = arg_needle_lib.tskit_to_arg(ts)
    positions = [.1, .2, .3, .4]
    genotypes = [
        [0, 1, 1, 0],
        [1, 0, 1, 0],
        [1, 1, 0, 1],
        [0, 0, 0, 0]]

    # Check we halt if no children
    with pytest.raises(RuntimeError):
        arg_needle_lib.map_genotype_to_ARG(arg, [0, 1, 1, 0], 0)

    arg.populate_children_and_roots()

    # Check mutations are correctly inserted
    for pos, g in zip(positions, genotypes):
        arg_needle_lib.map_genotype_to_ARG(arg, g, pos)

    me_0 = [m.edge for m in arg.mutations() if m.position == .1]
    expected_nodes_0 = {(2, 4), (1, 5)}
    assert set([(m.child.ID, m.parent.ID) for m in me_0]) == expected_nodes_0

    me_1 = [m.edge for m in arg.mutations() if m.position == .2]
    expected_nodes_1 = {(4, 6)}
    assert set([(m.child.ID, m.parent.ID) for m in me_1]) == expected_nodes_1

    me_2 = [m.edge for m in arg.mutations() if m.position == .3]
    expected_nodes_2 = {(0, 4), (5, 6)}
    assert set([(m.child.ID, m.parent.ID) for m in me_2]) == expected_nodes_2

    me_3 = [m.edge for m in arg.mutations() if m.position == .4]
    assert len(me_3) == 0

    # Check we halt on mutation carried by all which would segfault otherwise
    # See issue #140
    with pytest.raises(ValueError):
        arg_needle_lib.map_genotype_to_ARG(arg, [1, 1, 1, 1], 4)


def test_map_genotypes_small_async():
    """Test that input genotypes are correctly converted to mutations on edges"""
    seed = 130222
    ts = msprime.simulate(sample_size=4, random_seed=seed)
    # Generates the following tree
    # 1.23┊    6    ┊
    #     ┊  ┏━┻━┓  ┊
    # 0.85┊  ┃   5  ┊
    #     ┊  ┃  ┏┻┓ ┊
    # 0.17┊  4  ┃ ┃ ┊
    #     ┊ ┏┻┓ ┃ ┃ ┊
    # 0.00┊ 0 2 1 3 ┊
    #     0         1
    arg = arg_needle_lib.tskit_to_arg(ts)
    positions = [.1, .2, .3, .4]
    genotypes = [
        [0, 1, 1, 0],
        [1, 0, 1, 0],
        [1, 1, 0, 1],
        [0, 0, 0, 0]]

    # Check we halt if no children
    with pytest.raises(RuntimeError):
        arg_needle_lib.map_genotype_to_ARG(arg, [0, 1, 1, 0], 0)

    arg.populate_children_and_roots()

    # Check mutations are correctly inserted
    arg_needle_lib.map_genotypes_to_ARG(arg, genotypes, positions)

    me_0 = [m.edge for m in arg.mutations() if m.position == .1]
    expected_nodes_0 = {(2, 4), (1, 5)}
    assert set([(m.child.ID, m.parent.ID) for m in me_0]) == expected_nodes_0

    me_1 = [m.edge for m in arg.mutations() if m.position == .2]
    expected_nodes_1 = {(4, 6)}
    assert set([(m.child.ID, m.parent.ID) for m in me_1]) == expected_nodes_1

    me_2 = [m.edge for m in arg.mutations() if m.position == .3]
    expected_nodes_2 = {(0, 4), (5, 6)}
    assert set([(m.child.ID, m.parent.ID) for m in me_2]) == expected_nodes_2

    me_3 = [m.edge for m in arg.mutations() if m.position == .4]
    assert len(me_3) == 0

    # Check we halt on mutation carried by all which would segfault otherwise
    # See issue #140
    with pytest.raises(ValueError):
        arg_needle_lib.map_genotype_to_ARG(arg, [1, 1, 1, 1], 4)


def test_map_genotypes_big():
    """Test that input genotypes are correctly converted to mutations on edges"""
    seed = 130222
    ts = msprime.simulate(
        sample_size=400, Ne=2e4, length=1e5, recombination_rate=2e-8,
        mutation_rate=2e-8, random_seed=seed)
    G = ts.genotype_matrix()
    sites = [s.position for s in ts.sites()]

    arg = arg_needle_lib.tskit_to_arg(ts)
    arg.populate_children_and_roots()

    # Add genotypes
    for pos, g in zip(sites, G):
        arg_needle_lib.map_genotype_to_ARG(arg, g, pos)

    # For now just make sure we've added the right number of mutations.
    # Can test for accuracy once we have merged PR #137
    assert len(arg.mutations()) == len(G)


def test_map_diploid_genotypes_small():
    """Test that input genotypes are correctly converted to mutations on edges"""
    seed = 130222
    ts = msprime.simulate(sample_size=6, random_seed=seed)
    # Generates the following tree
    # 5.62┊      10     ┊
    #     ┊    ┏━━┻━━┓  ┊
    # 2.10┊    ┃     9  ┊
    #     ┊    ┃    ┏┻┓ ┊
    # 0.34┊    8    ┃ ┃ ┊
    #     ┊  ┏━┻━┓  ┃ ┃ ┊
    # 0.27┊  7   ┃  ┃ ┃ ┊
    #     ┊ ┏┻┓  ┃  ┃ ┃ ┊
    # 0.07┊ ┃ ┃  6  ┃ ┃ ┊
    #     ┊ ┃ ┃ ┏┻┓ ┃ ┃ ┊
    # 0.00┊ 0 4 1 3 2 5 ┊
    #     0             1
    arg = arg_needle_lib.tskit_to_arg(ts)
    positions = [.1, .2, .3, .4, .5]
    genotypes = [
        [1, 2, 2],
        [2, 0, 0],
        [1, 1, 1],
        [0, 2, 2],
        [2, 1, 1]]
    arg.populate_children_and_roots()

    # Check mutations are correctly inserted
    for pos, g in zip(positions, genotypes):
        arg_needle_lib.map_genotype_to_ARG_diploid(arg, g, pos)

    me_0 = [m.edge for m in arg.mutations() if m.position == positions[0]]
    expected_nodes_0 = [{(9, 10), (7, 8), (3, 6)}, {(9, 10), (6, 8), (4, 7)}]
    assert set([(m.child.ID, m.parent.ID) for m in me_0]) in expected_nodes_0

    me_1 = [m.edge for m in arg.mutations() if m.position == positions[1]]
    expected_nodes_1 = {(1, 6), (0, 7)}
    assert set([(m.child.ID, m.parent.ID) for m in me_1]) == expected_nodes_1

    me_2 = [m.edge for m in arg.mutations() if m.position == positions[2]]
    expected_nodes_2 = [{(6, 8), (5, 9)}, {(0, 7), (9, 10)}, {(7, 8), (2, 9)}]
    assert set([(m.child.ID, m.parent.ID) for m in me_2]) in expected_nodes_2

    me_3 = [m.edge for m in arg.mutations() if m.position == positions[3]]
    expected_nodes_3 = {(3, 6), (4, 7), (9, 10)}
    assert set([(m.child.ID, m.parent.ID) for m in me_3]) == expected_nodes_3

    me_4 = [m.edge for m in arg.mutations() if m.position == positions[4]]
    expected_nodes_4 = {(8, 10)}
    assert set([(m.child.ID, m.parent.ID) for m in me_4]) == expected_nodes_4

    # Check we halt on mutation carried by all which would segfault otherwise
    # See issue #140
    with pytest.raises(ValueError):
        arg_needle_lib.map_genotype_to_ARG(arg, [2, 2, 2], 4)


def test_map_genotypes_approximate():
    """Test that input genotypes are correctly converted to mutations on edges"""
    seed = 130222
    ts = msprime.simulate(sample_size=10, random_seed=seed)
    # 2.76┊           18        ┊
    #     ┊       ┏━━━━┻━━━━┓   ┊
    # 1.22┊      17         ┃   ┊
    #     ┊     ┏━┻━━┓      ┃   ┊
    # 0.93┊     ┃   16      ┃   ┊
    #     ┊     ┃   ┏┻┓     ┃   ┊
    # 0.59┊     ┃   ┃ ┃    15   ┊
    #     ┊     ┃   ┃ ┃   ┏━┻━┓ ┊
    # 0.58┊     ┃   ┃ ┃  14   ┃ ┊
    #     ┊     ┃   ┃ ┃  ┏┻━┓ ┃ ┊
    # 0.34┊    13   ┃ ┃  ┃  ┃ ┃ ┊
    #     ┊   ┏━┻━┓ ┃ ┃  ┃  ┃ ┃ ┊
    # 0.09┊  12   ┃ ┃ ┃  ┃  ┃ ┃ ┊
    #     ┊  ┏┻━┓ ┃ ┃ ┃  ┃  ┃ ┃ ┊
    # 0.08┊ 11  ┃ ┃ ┃ ┃  ┃  ┃ ┃ ┊
    #     ┊ ┏┻┓ ┃ ┃ ┃ ┃  ┃  ┃ ┃ ┊
    # 0.02┊ ┃ ┃ ┃ ┃ ┃ ┃ 10  ┃ ┃ ┊
    #     ┊ ┃ ┃ ┃ ┃ ┃ ┃ ┏┻┓ ┃ ┃ ┊
    # 0.00┊ 0 4 8 7 2 3 1 5 6 9 ┊
    #     0                     1
    arg = arg_needle_lib.tskit_to_arg(ts)
    positions = [.1, .2, .3, .4, .5]
    arg.populate_children_and_roots()
    # import pdb
    # pdb.set_trace()
    genotype1 = [1, 0, 1, 0, 1, 0, 0, 1, 1, 0]
    map1, _ = arg_needle_lib.map_genotype_to_ARG_approximate(arg, genotype1, positions[0])
    assert len(map1) == 0
    genotype2 = [1, 0, 1, 1, 1, 0, 0, 1, 1, 0]
    map2, _ = arg_needle_lib.map_genotype_to_ARG_approximate(arg, genotype2, positions[1])
    assert map2[0].child.ID == 15
    genotype3 = [1, 0, 0, 0, 1, 0, 0, 1, 1, 0]
    map3, _ = arg_needle_lib.map_genotype_to_ARG_approximate(arg, genotype3, positions[2])
    assert map3[0].child.ID == 13
    genotype4 = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0]
    map4, _ = arg_needle_lib.map_genotype_to_ARG_approximate(arg, genotype4, positions[3])
    assert len(map4) == 2


def test_mrca():
    seed = 130222
    ts = msprime.simulate(sample_size=10, random_seed=seed)
    # 2.76┊           18        ┊
    #     ┊       ┏━━━━┻━━━━┓   ┊
    # 1.22┊      17         ┃   ┊
    #     ┊     ┏━┻━━┓      ┃   ┊
    # 0.93┊     ┃   16      ┃   ┊
    #     ┊     ┃   ┏┻┓     ┃   ┊
    # 0.59┊     ┃   ┃ ┃    15   ┊
    #     ┊     ┃   ┃ ┃   ┏━┻━┓ ┊
    # 0.58┊     ┃   ┃ ┃  14   ┃ ┊
    #     ┊     ┃   ┃ ┃  ┏┻━┓ ┃ ┊
    # 0.34┊    13   ┃ ┃  ┃  ┃ ┃ ┊
    #     ┊   ┏━┻━┓ ┃ ┃  ┃  ┃ ┃ ┊
    # 0.09┊  12   ┃ ┃ ┃  ┃  ┃ ┃ ┊
    #     ┊  ┏┻━┓ ┃ ┃ ┃  ┃  ┃ ┃ ┊
    # 0.08┊ 11  ┃ ┃ ┃ ┃  ┃  ┃ ┃ ┊
    #     ┊ ┏┻┓ ┃ ┃ ┃ ┃  ┃  ┃ ┃ ┊
    # 0.02┊ ┃ ┃ ┃ ┃ ┃ ┃ 10  ┃ ┃ ┊
    #     ┊ ┃ ┃ ┃ ┃ ┃ ┃ ┏┻┓ ┃ ┃ ┊
    # 0.00┊ 0 4 8 7 2 3 1 5 6 9 ┊
    #     0                     1
    arg = arg_needle_lib.tskit_to_arg(ts)
    arg.populate_children_and_roots()
    # import pdb
    # pdb.set_trace()
    mrca1 = arg_needle_lib.most_recent_common_ancestor(arg, [0, 8, 7], 0.1)
    assert mrca1.ID == 13
    mrca2 = arg_needle_lib.most_recent_common_ancestor(arg, [0, 8], 0.1)
    assert mrca2.ID == 12
    mrca3 = arg_needle_lib.most_recent_common_ancestor(arg, [0], 0.1)
    assert mrca3.ID == 0
    mrca4 = arg_needle_lib.most_recent_common_ancestor(arg, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9], 0.1)
    assert mrca4.ID == 18
