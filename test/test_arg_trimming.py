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


import math

import arg_needle_lib
import msprime
import numpy as np
import pytest

seed = 42
ts = msprime.simulate(
    Ne=1e3, sample_size=100, recombination_rate=1e-7, length=5e5, random_seed=seed
)
arg = arg_needle_lib.tskit_to_arg(ts)
arg.populate_children_and_roots()

num_regions = 4
test_regions = np.linspace(0, ts.sequence_length, num_regions + 1)

paired_chunk_list = []
for i in range(num_regions):
    trimmed_arg = arg_needle_lib.trim_arg(
        arg, test_regions[i], test_regions[i + 1])
    trimmed_arg.populate_children_and_roots()
    paired_chunk_list.append(
        [trimmed_arg, arg, test_regions[i], test_regions[i + 1]]
    )


@pytest.mark.parametrize(
    "trimmed_arg, original_arg, trim_start, trim_end", paired_chunk_list
)
def test_trimmed_tree_statistics(
    trimmed_arg, original_arg, trim_start, trim_end
):
    """
    Test that the trimmed tree have the correct number of nodes, edges and volume.
    Also test that other metadata such as sequence length and start/end are correct.
    """
    # trimmed parts should have the same local volume in the original arg
    assert math.isclose(
        arg_needle_lib.local_volume(original_arg, trim_start, trim_end),
        arg_needle_lib.total_volume(trimmed_arg),
        rel_tol=1e-6
    )

    # stringent = False since tskit converter is used
    stringent = False
    trimmed_arg.check_basic(stringent)
    trimmed_arg.check_roots()
    trimmed_arg.check_children(stringent)

    # check metadata
    assert trimmed_arg.start == 0
    assert trimmed_arg.end == trim_end - trim_start
    assert trimmed_arg.offset == trim_start


def test_trimmed_tree_with_mutations():
    """Test with a manually built ARG that nodes, edges and mutations are correctly populated"""
    # Define a simple ARG
    arg = arg_needle_lib.ARG(0, 100, 3)
    arg.add_sample()
    arg.add_sample()
    arg.thread_sample([0, 60], [0, 0], [3.14, 2.718])
    arg.add_sample()
    arg.thread_sample([0, 30], [1, 0], [5, 3])

    arg.populate_children_and_roots()

    # The trees look like
    # 0-30                  30-60                   60-100
    # 5    ----      /\     3.14 ----      /\       3     ----      /\
    # 3.14 ----   /\   \    3    ----   /\   \      2.718 ----   /\   \
    # 0         0   1   2   0         0   2   1     0          0   1   2

    # Trimming without mutations and check marginal tree
    trimmed_arg = arg_needle_lib.trim_arg(arg, 30, 100)
    trimmed_arg.populate_children_and_roots()
    assert trimmed_arg.num_edges() == 7
    assert trimmed_arg.num_nodes() == 6
    assert trimmed_arg.root_at(0).node.height == 3.14
    assert trimmed_arg.root_at(40).node.height == 3
    assert trimmed_arg.mrca(0, 2, 0).height == 3
    assert trimmed_arg.mrca(1, 2, 30).height == 3
    assert trimmed_arg.mrca(0, 1, 30).height == 2.718

    # Trim with mutations
    arg_needle_lib.generate_m_mutations(arg, 30, random_seed=42)
    trimmed_arg = arg_needle_lib.trim_arg(arg, 20, 65)
    assert trimmed_arg.num_nodes() == 7
    assert trimmed_arg.num_edges() == 10

    # mutations should be sorted
    positions = []
    node_ids = []
    heights = []

    # manually trim the mutations
    for i, m in enumerate(arg.mutations()):
        if (20 <= m.position) & (m.position < 65):
            positions.append(m.position - 20)
            node_ids.append((m.edge.child.ID, m.edge.parent.ID))
            heights.append(m.height)

    for i, m in enumerate(trimmed_arg.mutations()):
        assert m.position == positions[i]
        assert trimmed_arg.get_site_positions()[i] == positions[i]
        assert m.height == heights[i]
        assert (m.edge.child.ID, m.edge.parent.ID) == node_ids[i]
