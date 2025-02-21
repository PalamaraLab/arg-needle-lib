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

"""
ARG serialization and deserialization functions
"""

from datetime import datetime

import h5py

import numpy as np

import arg_needle_lib

__all__ = [
    "serialize_arg",
]


def serialize_arg(arg_to_serialize: arg_needle_lib.ARG, file_name: str = "arg.arg", chunk_size: int = 1000,
                  compression='gzip', compression_opts=None, node_bounds=True, mutations=None):
    """
    Serialize ARG to disk using HDF5.

    By default, gzip compression with level 9 will be used, but available options are listed below. Compression
    can be turned off by passing `compression=None`.

    Args:
        arg_to_serialize: arg_needle_lib.ARG object to serialize
        file_name: the file name to serialize the ARG to
        chunk_size: the chunk size for writing; small numbers will dramatically increase serialization time
        compression: whether to use compression; can save up to ~50% on disk at a runtime cost
            Refer to h5py documentation for details
            `gzip` provides good compression with compression_opts=9, with ~+50% runtime for serialization
            `lzf` provides poor compression, with little runtime overhead
            `szip` provides good compression, with little runtime overhead, but is not available on all systems
        compression_opts: the level of compression requested for gzip: 9 for highest compression. Should be None if
            `compression != 'gzip'`
        node_bounds: whether to serialize node start/end positions. Default: True
        mutations: whether to serialize mutations. Default: True if num_mutations > 0, otherwise False
    """
    num_nodes = arg_to_serialize.num_nodes()
    num_edges = arg_to_serialize.num_edges()
    num_mutations = arg_to_serialize.num_mutations()

    # By default, serialize mutations if they exist, but this will not happen if set to False
    if mutations is None:
        mutations = num_mutations > 0

    if mutations and num_mutations == 0:
        print(
            f'WARNING: serialization parameter `mutations` is {mutations} but the ARG contains no mutations. Setting `mutations` to False'
        )
        mutations = False

    # Create the hdf5 file, write attributes and create datasets for nodes and edges
    f = h5py.File(file_name, "w")
    f.attrs['num_nodes'] = num_nodes
    f.attrs['num_edges'] = num_edges
    f.attrs['node_bounds'] = node_bounds
    f.attrs['num_mutations'] = num_mutations
    f.attrs['mutations'] = mutations
    f.attrs['offset'] = arg_to_serialize.offset
    f.attrs['chromosome'] = arg_to_serialize.chromosome
    f.attrs['start'] = arg_to_serialize.start
    f.attrs['end'] = arg_to_serialize.end
    f.attrs['threaded_samples'] = 0
    f.attrs['datetime_created'] = datetime.now().isoformat()
    f.attrs['arg_file_version'] = 2

    if compression is None:
        dset_flags = f.create_dataset("flags", (num_nodes,), dtype=bool)
        dset_times = f.create_dataset("times", (num_nodes,), dtype=np.double)
        dset_edge_ranges = f.create_dataset("edge_ranges", (num_edges, 2), dtype=np.double)
        dset_edge_ids = f.create_dataset("edge_ids", (num_edges, 2), dtype=np.int32)

        if node_bounds:
            dset_node_bounds = f.create_dataset("node_bounds", (num_nodes, 2), dtype=np.double)

    else:
        if compression == 'gzip' and compression_opts is None:
            compression_opts = 9

        dset_flags = f.create_dataset("flags", (num_nodes,), dtype=bool, compression=compression,
                                      compression_opts=compression_opts)
        dset_times = f.create_dataset("times", (num_nodes,), dtype=np.double, compression=compression,
                                      compression_opts=compression_opts)
        dset_edge_ranges = f.create_dataset("edge_ranges", (num_edges, 2), dtype=np.double, compression=compression,
                                            compression_opts=compression_opts)
        dset_edge_ids = f.create_dataset("edge_ids", (num_edges, 2), dtype=np.int32, compression=compression,
                                         compression_opts=compression_opts)

        if node_bounds:
            dset_node_bounds = f.create_dataset("node_bounds", (num_nodes, 2), dtype=np.double, compression=compression,
                                                compression_opts=compression_opts)

    num_nodes_written = 0
    num_edges_written = 0

    # Process {chunk_size} nodes at a time, writing each chunk to file as we go
    while num_nodes_written < num_nodes:

        range_lo = num_nodes_written
        range_hi = min(num_nodes_written + chunk_size, num_nodes)
        range_len = range_hi - range_lo

        # Count the number of edges associated with this chunk
        edges_this_range = 0
        for node_idx in range(range_lo, range_hi):
            edges_this_range += len(arg_to_serialize.node(node_idx).parent_edges())

        # Create numpy arrays to organise the data
        flags = np.zeros((range_len,), dtype=bool)
        times = np.zeros((range_len,), dtype=np.double)
        node_bounds_data = np.zeros((range_len, 2), dtype=np.double)
        edge_ranges = np.zeros((edges_this_range, 2), dtype=np.double)
        edge_ids = np.zeros((edges_this_range, 2), dtype=np.int32)

        # Process all nodes and edges in this current chunk
        chunk_node_idx = 0
        chunk_edge_idx = 0
        for node_idx in range(range_lo, range_hi):
            node = arg_to_serialize.node(node_idx)
            flags[chunk_node_idx] = arg_to_serialize.is_leaf(node_idx)
            times[chunk_node_idx] = node.height
            node_bounds_data[chunk_node_idx, 0] = node.start
            node_bounds_data[chunk_node_idx, 1] = node.end
            chunk_node_idx += 1
            for edge in node.parent_edges():
                edge_ranges[chunk_edge_idx, 0] = edge.start
                edge_ranges[chunk_edge_idx, 1] = edge.end
                edge_ids[chunk_edge_idx, 0] = node_idx
                edge_ids[chunk_edge_idx, 1] = edge.parent.ID
                chunk_edge_idx += 1

        # Update threaded_samples
        f.attrs['threaded_samples'] += np.count_nonzero(flags)

        # Write data to hdf5 file
        dset_flags[num_nodes_written:num_nodes_written + range_len] = flags
        dset_times[num_nodes_written:num_nodes_written + range_len] = times
        dset_edge_ranges[num_edges_written:num_edges_written + edges_this_range, :] = edge_ranges
        dset_edge_ids[num_edges_written:num_edges_written + edges_this_range, :] = edge_ids

        if node_bounds:
            dset_node_bounds[num_nodes_written:num_nodes_written + range_len, :] = node_bounds_data

        # Increment counters
        num_nodes_written += range_len
        num_edges_written += edges_this_range

    assert num_nodes_written == num_nodes
    assert num_edges_written == num_edges

    # Write mutations if required
    if mutations:
        if compression is None:
            mut = f.create_group("mutations")
            mut_pos = mut.create_dataset("positions", (num_mutations,), dtype=np.double)
            mut_hts = mut.create_dataset("heights", (num_mutations,), dtype=np.double)
            mut_sid = mut.create_dataset("site_ids", (num_mutations,), dtype=np.int32)
            mut_edge_ids = mut.create_dataset("edge_ids", (num_mutations, 2), dtype=np.int32)

        else:
            if compression == 'gzip' and compression_opts is None:
                compression_opts = 9

            mut = f.create_group("mutations")
            mut_pos = mut.create_dataset("positions", (num_mutations,), dtype=np.double, compression=compression,
                                         compression_opts=compression_opts)
            mut_hts = mut.create_dataset("heights", (num_mutations,), dtype=np.double, compression=compression,
                                         compression_opts=compression_opts)
            mut_sid = mut.create_dataset("site_ids", (num_mutations,), dtype=np.int32, compression=compression,
                                         compression_opts=compression_opts)
            mut_edge_ids = mut.create_dataset("edge_ids", (num_mutations, 2), dtype=np.int32, compression=compression,
                                              compression_opts=compression_opts)

        # Create numpy arrays to organise the data
        positions = np.empty((num_mutations,), dtype=np.double)
        heights = np.empty((num_mutations,), dtype=np.double)
        site_ids = np.empty((num_mutations,), dtype=np.int32)
        edge_ids = np.empty((num_mutations, 2), dtype=np.int32)

        for idx, mutation in enumerate(arg_to_serialize.mutations()):
            positions[idx] = mutation.position
            heights[idx] = mutation.height
            site_ids[idx] = mutation.site_id
            edge_ids[idx, 0] = mutation.edge.child.ID
            edge_ids[idx, 1] = mutation.edge.parent.ID

        mut_pos[...] = positions
        mut_hts[...] = heights
        mut_sid[...] = site_ids
        mut_edge_ids[...] = edge_ids

    f.close()
