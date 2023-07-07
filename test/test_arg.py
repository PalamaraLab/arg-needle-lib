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

import pytest
import msprime
import arg_needle_lib

def test_populate_mutations_on_edges():
    """Test that mutations are properly populated onto edges"""
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
    # tskit_to_arg does not yet support mutations
    arg = arg_needle_lib.tskit_to_arg(ts)
    # we generate mutations using arg_needle_lib and get
    # [0, 0, 1, 0]
    # [0, 0, 0, 1]
    # [0, 0, 0, 1]
    # [0, 1, 0, 0]
    # [0, 1, 0, 0]
    # (though not necessarily in that order)
    arg_needle_lib.generate_m_mutations(arg, 5, random_seed=seed)

    # Verify there are no mutations on the branches
    for node_id in range(6):
        assert len(arg.node(node_id).parent_edge_at(0.5).mutations()) == 0

    # Add mutations to edges
    arg.populate_mutations_on_edges()
    # Verify there are now mutations on the branches
    assert len(arg.node(0).parent_edge_at(0.5).mutations()) == 0
    assert len(arg.node(1).parent_edge_at(0.5).mutations()) == 2
    assert len(arg.node(2).parent_edge_at(0.5).mutations()) == 1
    assert len(arg.node(3).parent_edge_at(0.5).mutations()) == 2
    assert len(arg.node(4).parent_edge_at(0.5).mutations()) == 0
    assert len(arg.node(5).parent_edge_at(0.5).mutations()) == 0


def test_lowest_mutated_edge():
    """
    Test that the lowest parent edge carrying a mutation at
    a given site is correctly identified.
    """
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
    # tskit_to_arg does not yet support mutations
    arg = arg_needle_lib.tskit_to_arg(ts)
    # we generate mutations using arg_needle_lib and get
    # [0, 1, 0, 0] at 0.5320779165643559
    # [1, 0, 1, 0] at 0.039759496502559354
    arg_needle_lib.generate_m_mutations(arg, 2, random_seed=12345)
    arg.populate_mutations_on_edges()
    # Error tolerance is 1e-3
    assert arg.lowest_mutated_edge(2, 0.0397).child.ID == 4
    assert arg.lowest_mutated_edge(0, 0.0397).child.ID == 4
    assert arg.lowest_mutated_edge(0, 0.0387594) is None
    assert arg.lowest_mutated_edge(2, 0.5320) is None
    assert arg.lowest_mutated_edge(1, 0.5320).child.ID == 1
    assert arg.lowest_mutated_edge(3, 0.5320) is None

def test_set_and_get_sites():
    """Test that sites are properly set"""
    # Generate an ARG with 4 leaves [0, 1]
    seed = 130222
    ts = msprime.simulate(sample_size=4, random_seed=seed)
    arg = arg_needle_lib.tskit_to_arg(ts)

    # Set sites
    arg.set_sites([0.1, 0.5, 0.6])
    assert arg.num_sites() == 3

    # Overwrite sites
    arg.set_sites([0.1, 0.5, 0.6, 0.9])
    assert arg.num_sites() == 4

    # Test getter
    assert arg.get_site(0) == 0.1
    assert arg.get_site(3) == 0.9
    with pytest.raises(Exception):
        arg.get_site(4)

    # Get the closest site ID for various positions
    query_positions = [0, 0.4, 0.51, 0.7, 0.95]
    closest_site_ids = [0, 1, 1, 2, 3]
    results = []
    for position in query_positions:
        results.append(arg.get_id_of_closest_site(position))
    assert results == closest_site_ids

    with pytest.raises(Exception):
        arg.get_id_of_closest_site(-0.5)

    # Try to set unsorted sites
    with pytest.raises(Exception):
        arg.set_sites([0.5, 0.1, 0.6, 0.9])

    # Try to set out-of-bounds sites
    with pytest.raises(Exception):
        arg.set_sites([0.1, 0.5, 1.1])


def test_finding_mutations():
    """Test methods for finding mutations left/right/nearest"""

    arg = arg_needle_lib.tskit_to_arg(msprime.simulate(sample_size=40))
    arg_needle_lib.generate_m_mutations(arg, 100)
    assert arg.num_mutations() == 100

    half_way = 0.5 * (arg.start + arg.end)

    idx_left_of = arg.get_idx_of_first_mutation_left_of(half_way)
    assert idx_left_of < 99
    assert arg.mutations()[idx_left_of].position < half_way
    assert arg.mutations()[idx_left_of + 1].position >= half_way

    idx_right_of = arg.get_idx_of_first_mutation_right_of(half_way, include_equal=True)
    assert idx_right_of > 0
    assert arg.mutations()[idx_right_of].position >= half_way
    assert arg.mutations()[idx_right_of - 1].position < half_way

    idx_closest = arg.get_idx_of_mutation_closest_to(half_way)
    assert 0 < idx_closest < 99
    dist_closest = abs(arg.mutations()[idx_closest].position - half_way)
    dist_left = abs(arg.mutations()[idx_closest - 1].position - half_way)
    dist_right = abs(arg.mutations()[idx_closest + 1].position - half_way)
    assert dist_closest <= dist_left and dist_closest <= dist_right


def test_mutation_sites():
    """Test mutation_sites map"""

    arg = arg_needle_lib.tskit_to_arg(msprime.simulate(sample_size=40))
    arg_needle_lib.generate_m_mutations(arg, 10)

    # It is technically possible that two mutations may have been generated with the same
    # position, but this is essentially impossible.
    mutation_sites = arg.get_mutation_sites()
    assert len(mutation_sites) == 10

    for key, val in mutation_sites.items():
        assert key == val.get_mutations()[0].position
