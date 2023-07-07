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

"""
ARG serialization and deserialization functions
"""

from datetime import datetime

import h5py
import pathlib

import numpy as np

import arg_needle_lib

__all__ = [
    "serialize_arg",
    "deserialize_arg",
    "validate_serialized_arg",
]


def serialize_arg(arg_to_serialize: arg_needle_lib.ARG, file_name: str = "arg.arg", chunk_size: int = 1000,
                  compression='gzip', compression_opts=None, node_bounds=True, mutations=None, sites=None):
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
        sites: whether to serialize sites. Default: True if sites present, otherwise False
    """
    num_nodes = arg_to_serialize.num_nodes()
    num_edges = arg_to_serialize.num_edges()
    num_sites = arg_to_serialize.num_sites()
    num_mutations = arg_to_serialize.num_mutations()

    # By default, serialize mutations if they exist, but this will not happen if set to False
    if mutations is None:
        mutations = num_mutations > 0

    if mutations and num_mutations == 0:
        print(
            f'WARNING: serialization parameter `mutations` is {mutations} but the ARG contains no mutations. Setting `mutations` to False'
        )
        mutations = False

    # By default, serialize sites if they exist, but this will not happen if set to False
    if sites is None:
        sites = num_sites > 0

    if sites and num_sites == 0:
        print(
            f'WARNING: serialization parameter `sites` is {sites} but the ARG contains no sites. Setting `sites` to False'
        )
        sites = False

    # Create the hdf5 file, write attributes and create datasets for nodes and edges
    f = h5py.File(file_name, "w")
    f.attrs['num_nodes'] = num_nodes
    f.attrs['num_edges'] = num_edges
    f.attrs['num_sites'] = num_sites
    f.attrs['node_bounds'] = node_bounds
    f.attrs['sites'] = sites
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

    # Write sites if required
    if sites:
        if compression is None:
            dset_sites = f.create_dataset("sites", (num_sites,), dtype=np.double)
        else:
            if compression == 'gzip' and compression_opts is None:
                compression_opts = 9

            dset_sites = f.create_dataset("sites", (num_sites,), dtype=np.double, compression=compression,
                                          compression_opts=compression_opts)
        
        dset_sites[...] = arg_to_serialize.get_sites()

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


def deserialize_arg(file_name: str, chunk_size: int = 1000, reserved_samples: int = -1):
    """
    Deserialize ARG from disk and create an arg_needle_lib.ARG object to return

    Args:
        file_name: the file name of the serialized ARG file
        chunk_size: the chunk size for reading; small numbers will dramatically increase deserialization time
        reserved_samples: parameter for ARG constructor allowing us to reserve additional samples for threading
            (default is to not reserve additional samples)

    Return:
        an arg_needle_lib.ARG object
    """
    validate_serialized_arg(file_name)

    with h5py.File(file_name, "r") as f:
        arg_file_version = f.attrs['arg_file_version']
    
    if arg_file_version == 1:
        return _deserialize_arg_v1(file_name, reserved_samples)
    elif arg_file_version == 2:
        return _deserialize_arg_v2(file_name, chunk_size, reserved_samples)


def _deserialize_arg_v1(file_name: str, reserved_samples: int = -1):
    """
    """
    f = h5py.File(file_name, "r")

    offset = f.attrs['offset']
    chromosome = f.attrs['chromosome']
    sequence_length = f.attrs['sequence_length']

    is_sample = f['flags'][...]
    node_heights = f['times'][...]
    edge_ids = f['edge_ids'][...]
    edge_ranges = f['edge_ranges'][...]

    arg = arg_needle_lib.ARG(0, sequence_length, node_heights, is_sample, edge_ids, edge_ranges, reserved_samples)
    arg.set_offset(offset)
    arg.set_chromosome(chromosome)
    return arg


def _deserialize_arg_v2(file_name: str, chunk_size: int = 1000, reserved_samples: int = -1):
    """
    """
    f = h5py.File(file_name, "r")

    dp = arg_needle_lib.DeserializationParams()
    dp.start = f.attrs['start']
    dp.end = f.attrs['end']
    dp.num_nodes = f.attrs['num_nodes']
    dp.offset = f.attrs['offset']
    dp.chromosome = f.attrs['chromosome']
    dp.threaded_samples = f.attrs['threaded_samples']
    dp.reserved_samples = reserved_samples

    arg = arg_needle_lib.ARG(deserialization_params=dp)

    # Process {chunk_size} nodes at a time, adding each chunk to the ARG as we go
    num_nodes = f.attrs['num_nodes']
    num_nodes_written = 0

    while num_nodes_written < num_nodes:

        range_lo = num_nodes_written
        range_hi = min(num_nodes_written + chunk_size, num_nodes)
        range_len = range_hi - range_lo

        node_heights = f['times'][num_nodes_written:num_nodes_written + range_len]
        is_sample = f['flags'][num_nodes_written:num_nodes_written + range_len]

        if f.attrs['node_bounds']:
            node_bounds_data = f['node_bounds'][num_nodes_written:num_nodes_written + range_len, :]
            arg.deserialize_add_nodes(node_heights, is_sample, node_bounds_data)
        else:
            arg.deserialize_add_nodes(node_heights, is_sample)

        num_nodes_written += range_len

    # Process {chunk_size} edges at a time, adding each chunk to the ARG as we go
    num_edges = f.attrs['num_edges']
    num_edges_written = 0
    while num_edges_written < num_edges:

        range_lo = num_edges_written
        range_hi = min(num_edges_written + chunk_size, num_edges)
        range_len = range_hi - range_lo

        edge_ids = f['edge_ids'][num_edges_written:num_edges_written + range_len, :]
        edge_ranges = f['edge_ranges'][num_edges_written:num_edges_written + range_len, :]

        arg.deserialize_add_edges(edge_ids, edge_ranges)

        num_edges_written += range_len
    
    # Add the sites; not necessary to do this in chunks
    if f.attrs['sites']:
        arg.set_sites(f['sites'][...])

    if f.attrs['mutations']:

        # Process {chunk_size} mutations at a time, adding each chunk to the ARG as we go
        num_mutations = f.attrs['num_mutations']
        num_mutations_written = 0
        while num_mutations_written < num_mutations:

            range_lo = num_mutations_written
            range_hi = min(num_mutations_written + chunk_size, num_mutations)
            range_len = range_hi - range_lo

            mut_pos = f['mutations']['positions'][num_mutations_written:num_mutations_written + range_len]
            mut_hts = f['mutations']['heights'][num_mutations_written:num_mutations_written + range_len]
            mut_sid = f['mutations']['site_ids'][num_mutations_written:num_mutations_written + range_len]
            mut_eid = f['mutations']['edge_ids'][num_mutations_written:num_mutations_written + range_len, :]

            arg.deserialize_add_mutations(mut_pos, mut_hts, mut_sid, mut_eid)

            num_mutations_written += range_len

    return arg


def validate_serialized_arg(file_name: str):
    """
    Validate a serialized ARG on disk. This method returns True if the file is valid, and False otherwise. A message
    indicating the fault will be printed if the serialized ARG file is invalid.

    Args:
        file_name: the file name of the serialized ARG file

    Return:
        whether the serialized ARG file is valid
    """

    if not pathlib.Path(file_name).is_file():
        print(f'File: {file_name} is not a valid file')
        return False

    if not h5py.is_hdf5(file_name):
        print(f'File: {file_name} is not a valid HDF5 file')
        return False

    with h5py.File(file_name, "r") as f:
        if 'arg_file_version' not in f.attrs.keys():
            print(f'File: {file_name} is not a valid arg file because it does not contain `arg_file_version` attribute')
            return False

    valid_file_versions = [1, 2]
    arg_file_version = h5py.File(file_name, "r").attrs['arg_file_version']

    if arg_file_version not in valid_file_versions:
        str_valid = ', '.join([str(x) for x in valid_file_versions])
        print(f'Arg file version ({arg_file_version}) is not supported; valid versions are {str_valid}')
        return False
    
    if arg_file_version == 1:
        return _validate_serialized_arg_v1(file_name)
    elif arg_file_version == 2:
        return _validate_serialized_arg_v2(file_name)


def _validate_serialized_arg_v1(file_name: str):
    """
    Validate v1 serialized ARG files. This method checks that the expected attributes and datasets are present in the
    file.

    Args:
        file_name: the file name of the serialized ARG file

    Return:
        whether the serialized ARG file is valid
    """
    expected_attrs = ['num_nodes', 'num_edges', 'num_mutations', 'offset', 'chromosome', 'sequence_length',
                      'datetime_created', 'arg_file_version']
    expected_dsets = ['flags', 'times', 'edge_ranges', 'edge_ids']

    with h5py.File(file_name, "r") as f:
        for attr in expected_attrs:
            if attr not in f.attrs.keys():
                print(f'Expected file {file_name} to include attribute `{attr}`')
                return False

        for dset in expected_dsets:
            if dset not in f.keys():
                print(f'Expected file {file_name} to include dataset `{dset}`')
                return False

    return True


def _validate_serialized_arg_v2(file_name: str):
    """
    Validate v2 serialized ARG files. This method checks that the expected attributes and datasets are present in the
    file.

    Args:
        file_name: the file name of the serialized ARG file

    Return:
        whether the serialized ARG file is valid
    """

    expected_attrs = ['num_nodes', 'num_edges', 'num_sites', 'node_bounds', 'sites', 'num_mutations', 'mutations',
                      'offset', 'chromosome', 'start', 'end', 'threaded_samples', 'datetime_created',
                      'arg_file_version']
    expected_dsets = ['flags', 'times', 'edge_ranges', 'edge_ids']

    with h5py.File(file_name, "r") as f:
        for attr in expected_attrs:
            if attr not in f.attrs.keys():
                print(f'Expected file {file_name} to include attribute `{attr}`')
                return False

        for dset in expected_dsets:
            if dset not in f.keys():
                print(f'Expected file {file_name} to include dataset `{dset}`')
                return False
        
        # Optional datasets
        if f.attrs['sites']:
            if 'sites' not in f.keys():
                print(f'Expected file {file_name} to include dataset `sites`')
                return False

        if f.attrs['mutations']:
            if 'mutations' not in f.keys():
                print(f'Expected file {file_name} to include dataset `mutations`')
                return False

        if f.attrs['node_bounds']:
            if 'node_bounds' not in f.keys():
                print(f'Expected file {file_name} to include dataset `node_bounds`')
                return False

    return True
