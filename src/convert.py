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

"""Converting between ARG and other formats."""

import numpy as np
import tskit

import arg_needle_lib

__all__ = [
    "arg_to_tskit",
    "tskit_to_arg_thread",
    "tskit_to_arg"
]


def arg_to_tskit(arg, batch_size=None, mutations=True, sample_permutation=None):
    """Convert arg_needle_lib ARG object to tskit.TreeSequence using tables

    For this to work properly, the node IDs in the ARG object must all be consecutive.

    Args:
        arg: arg_needle_lib.ARG object to convert
        batch_size (int): Maximum allowed size for edge list information
          mainly used to save memory (default is no batching)
        mutations (bool): If True, converts mutations as well (default=True)
        sample_permutation: If not None, then should be a permutation of values from 0
          to num_samples - 1, either as a numpy array or Python list. Sample i in the
          resulting tskit.TreeSequence will correspond to sample sample_permutation[i]
          in the arg_needle_lib.ARG.

    Returns:
        Converted tskit.TreeSequence object.

    Raises:
        ValueError: Error in ARG format, batch_size parameter, or sample_permutation.
    """
    if arg.start != 0:
        raise ValueError("Can only convert when ARG start is 0")
    if batch_size is not None:
        if batch_size <= 0 or not isinstance(batch_size, int):
            raise ValueError("Batch size must either be None or a positive integer.")
    if -1 in arg.node_ids():
        node_ids = [n for n in arg.node_ids() if n >= 0]
        num_nodes = arg.num_nodes() - 1
    else:
        node_ids = arg.node_ids()
        num_nodes = arg.num_nodes()
    if num_nodes != max(node_ids) + 1:
        raise ValueError("All nodes should be taken.")
    num_samples = len(arg.sample_names)
    if sample_permutation is not None:
        # Check that sample_permutation is a permutation of elements 0 to arg.num_samples - 1
        if len(sample_permutation) != num_samples:
            raise ValueError("Permutation should be of length equal to number of samples.")
        if sorted(sample_permutation) != np.arange(num_samples).tolist():
            raise ValueError("Not a valid permutation.")

    # Setup the full table
    tables = tskit.TableCollection(sequence_length=arg.end)
    # Lists for the node table
    flags = []
    times = []
    # Lists for the edge table
    lefts = []
    rights = []
    parents = []
    children = []
    # Counter used for batching
    edge_counter = 0
    for node_id in range(num_nodes):
        # We allow permuting the samples to facilitate different threading orders. To recover
        # the original order, we pass in a permutation over the sample IDs.
        if node_id < num_samples and sample_permutation is not None:
            permuted_node_id = sample_permutation[node_id]
        else:
            permuted_node_id = node_id

        node = arg.node(permuted_node_id)
        times.append(node.height)
        flags.append(tskit.NODE_IS_SAMPLE if arg.is_leaf(permuted_node_id) else ~tskit.NODE_IS_SAMPLE)
        for edge in node.parent_edges():
            if edge.parent.ID != -1:
                lefts.append(edge.start)
                rights.append(edge.end)
                parents.append(edge.parent.ID)
                children.append(node_id) # Since this is in tskit space, the node ID needs to be not permuted
                edge_counter += 1
        if batch_size is not None and edge_counter >= batch_size:
            # Append all of the nodes
            tables.nodes.append_columns(flags=flags, time=times)
            tables.edges.append_columns(left=lefts, right=rights, parent=parents, child=children)
            # Reset the node info lists
            flags = []
            times = []
            # Reset the edge info lists
            lefts = []
            rights = []
            parents = []
            children = []
            # Reset the counter
            edge_counter = 0
    # Append any residual columns
    tables.nodes.append_columns(flags=flags, time=times)
    tables.edges.append_columns(left=lefts, right=rights, parent=parents, child=children)
    if mutations:
        for m in arg.mutations():
            # Note: currently supports infinite-sites mutation
            height = m.height if m.height >= 0 else tskit.UNKNOWN_TIME
            site_id = tables.sites.add_row(position=m.position, ancestral_state=b'0')
            tables.mutations.add_row(site=site_id, node=m.edge.child.ID, time=height, derived_state=b'1')
    # Sort the edges in order of increasing parent height
    tables.sort()

    # Set the offset and chromosome
    offset = arg.offset if hasattr(arg, "offset") else 0
    chromosome = arg.chromosome if hasattr(arg, "chromosome") else 1
    assert isinstance(offset, int) and offset >= 0
    assert isinstance(chromosome, int) and chromosome >= 1
    meta_schema = {
        'codec': 'json',
        'properties': {
            'offset': {'type': 'integer'},
            'chromosome': {'type': 'integer'}
        }
    }
    tables.metadata_schema = tskit.MetadataSchema(meta_schema)
    tables.metadata = {"offset": offset, "chromosome": chromosome}

    return(tables.tree_sequence())


def single_threading_sequence(ts, new_id, sample_ids):
    """Generate a single threading sequence.

    There can be multiple possible threading sequences. Roughly, out of
    multiple possible nearest cousins to thread to, we choose the one with
    lowest Python hash.

    Faster implementations are probably possible; for now this is just used in
    tskit_to_arg_thread within some tests as a sanity check.

    Args:
        ts: tskit.TreeSequence object
        new_id (int): An ID for a sample in the TreeSequence
        sample_ids (list of ints): The sample IDs "already in the ARG" to compare to
          when querying this threading sequence.

    Returns:
        Tuple of positions, closest cousin IDs, and times that make up the
          threading instruction.
    """
    if isinstance(sample_ids, (list,)):
        sample_ids = set(sample_ids)
    starts = []
    cousin_ids = []
    times = []
    cousin = None
    tmrca = None
    for tree in ts.trees():
        cousins = set()
        node_id = new_id
        while cousins == set():
            node_id = tree.parent(node_id)
            descendants = list(tree.leaves(node_id))
            cousins = set(descendants) & sample_ids
        new_tmrca = tree.time(node_id)
        if cousin is None or cousin not in cousins or tmrca != new_tmrca:
            cousin = cousins.pop()
            tmrca = new_tmrca
            starts.append(tree.interval[0])
            cousin_ids.append(cousin)
            times.append(new_tmrca)

    return starts, cousin_ids, times


def tskit_to_arg_thread(ts):
    """Construct the ARG using threading operations.

    A good sanity check for the correctness of the arg_needle_lib
    threading operation.

    Args:
        ts: tskit.TreeSequence object

    Returns:
        arg_needle_lib.ARG object that was constructed to be identical using threading.
    """
    offset = 0
    chromosome = 1
    if isinstance(ts.metadata, dict):
        if "offset" in ts.metadata:
            offset = ts.metadata["offset"]
        if "chromosome" in ts.metadata:
            chromosome = ts.metadata["chromosome"]

    assert isinstance(offset, int) and offset >= 0
    assert isinstance(chromosome, int) and chromosome >= 1

    arg = arg_needle_lib.ARG(0, ts.sequence_length, ts.num_samples)
    arg.set_offset(offset)
    arg.set_chromosome(chromosome)

    arg.add_sample(str(0))
    sample_ids = set()
    sample_ids.add(0)
    for i in range(1, ts.num_samples):
        starts, cousin_ids, times = single_threading_sequence(ts, i, sample_ids)
        arg.add_sample(str(i))
        arg.thread_sample(starts, cousin_ids, times)
        sample_ids.add(i)
    return arg


def tskit_to_arg(ts, reserved_samples=-1):
    """Convert tskit.TreeSequence to arg_needle_lib.ARG using a custom constructor

    Args:
        ts: tskit.TreeSequence object
        reserved_samples (int): Parameter from ARG constructor allowing us to reserve additional
          samples for threading (default is to not reserve additional samples).

    Returns:
        An arg_needle_lib.ARG object.
    """
    num_nodes = ts.num_nodes
    num_edges = ts.num_edges

    offset = 0
    chromosome = 1
    if isinstance(ts.metadata, dict):
        if "offset" in ts.metadata:
            offset = ts.metadata["offset"]
        if "chromosome" in ts.metadata:
            chromosome = ts.metadata["chromosome"]

    assert isinstance(offset, int) and offset >= 0
    assert isinstance(chromosome, int) and chromosome >= 1

    node_heights = np.zeros(num_nodes, dtype=np.float64)
    is_sample = np.zeros(num_nodes, dtype=bool)
    edge_ids = np.zeros((num_edges, 2), dtype=int)
    edge_ranges = np.zeros((num_edges, 2), dtype=np.float64)

    node_index = 0
    edge_index = 0
    for i, node in enumerate(ts.nodes()):
        assert i == node.id
        node_heights[node_index] = node.time
        if node.flags == tskit.NODE_IS_SAMPLE:
            is_sample[node_index] = True
        else:
            is_sample[node_index] = False
        node_index += 1

    for edge in ts.edges():
        edge_ids[edge_index][0] = edge.child
        edge_ids[edge_index][1] = edge.parent
        edge_ranges[edge_index][0] = edge.left
        edge_ranges[edge_index][1] = edge.right
        edge_index += 1
    arg_sequence_length = ts.sequence_length

    arg = arg_needle_lib.ARG(
        0, arg_sequence_length, node_heights, is_sample,
        edge_ids, edge_ranges, reserved_samples)
    arg.set_offset(offset)
    arg.set_chromosome(chromosome)
    return arg
