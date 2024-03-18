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
Tests for ARG serialization and deserialization functions
"""

import h5py
import msprime
import pathlib
import sys

import numpy as np

import arg_needle_lib


def test_serialize_chunk_sizes():
    """
    Test that serialised ARGS are identical regardless of the chunk size used to create them. This essentially
    determines that there are no edge-case errors in the chunking code: we test writing one node at a time and writing
    the whole ARG as a single chunk.
    """
    ts = msprime.simulate(sample_size=20, Ne=1e4, length=5e3, recombination_rate=2e-8, mutation_rate=2e-8,
                          random_seed=1234)

    arg = arg_needle_lib.tskit_to_arg(ts)
    arg_needle_lib.serialize_arg(arg, "small_chunk.arg", 1)
    arg_needle_lib.serialize_arg(arg, "med_chunk.arg", 13)
    arg_needle_lib.serialize_arg(arg, "one_chunk.arg", 100)
    assert arg_needle_lib.validate_serialized_arg("small_chunk.arg")
    assert arg_needle_lib.validate_serialized_arg("med_chunk.arg")
    assert arg_needle_lib.validate_serialized_arg("one_chunk.arg")

    f_small_chunk = h5py.File("small_chunk.arg", "r")
    f_med_chunk = h5py.File("med_chunk.arg", "r")
    f_one_chunk = h5py.File("one_chunk.arg", "r")

    assert np.allclose(f_small_chunk['flags'], f_med_chunk['flags'])
    assert np.allclose(f_small_chunk['flags'], f_one_chunk['flags'])

    assert np.allclose(f_small_chunk['times'], f_med_chunk['times'])
    assert np.allclose(f_small_chunk['times'], f_one_chunk['times'])

    assert np.allclose(f_small_chunk['edge_ranges'], f_med_chunk['edge_ranges'])
    assert np.allclose(f_small_chunk['edge_ranges'], f_one_chunk['edge_ranges'])

    assert np.allclose(f_small_chunk['edge_ids'], f_med_chunk['edge_ids'])
    assert np.allclose(f_small_chunk['edge_ids'], f_one_chunk['edge_ids'])

    f_small_chunk.close()
    f_med_chunk.close()
    f_one_chunk.close()

    # Tidy up files on disk
    pathlib.Path("small_chunk.arg").unlink()
    pathlib.Path("med_chunk.arg").unlink()
    pathlib.Path("one_chunk.arg").unlink()


def test_compression():
    """
    Test that all compression options (and no compression) work without errors
    """
    ts = msprime.simulate(sample_size=20, Ne=1e4, length=5e3, recombination_rate=2e-8, mutation_rate=2e-8,
                          random_seed=1234)
    arg = arg_needle_lib.tskit_to_arg(ts)

    arg_needle_lib.serialize_arg(arg, 'no_compression.arg', compression=None, compression_opts=None)
    arg_needle_lib.serialize_arg(arg, 'gzip.arg', compression='gzip', compression_opts=9)
    arg_needle_lib.serialize_arg(arg, 'gzip.arg')
    arg_needle_lib.serialize_arg(arg, 'lzf.arg', compression='lzf')
    assert arg_needle_lib.validate_serialized_arg("no_compression.arg")
    assert arg_needle_lib.validate_serialized_arg("gzip.arg")
    assert arg_needle_lib.validate_serialized_arg("lzf.arg")
    
    if sys.platform == 'linux':
        arg_needle_lib.serialize_arg(arg, 'szip.arg', compression='szip')
        assert arg_needle_lib.validate_serialized_arg("szip.arg")

    # Tidy up files on disk
    pathlib.Path("no_compression.arg").unlink()
    pathlib.Path("gzip.arg").unlink()
    pathlib.Path("lzf.arg").unlink()
    if sys.platform == 'linux':
        pathlib.Path("szip.arg").unlink()


def test_round_trips():
    """
    Create an ARG in memory, then serialize -> deserialize -> serialize, and check that the two copies in memory match,
    and the two copies on disk match.

    This test demonstrates that the serialization/deserialization is a round-trip process.
    """
    ts = msprime.simulate(sample_size=20, Ne=1e4, length=5e3, recombination_rate=2e-8, mutation_rate=2e-8,
                          random_seed=1234)

    arg_original = arg_needle_lib.tskit_to_arg(ts)

    arg_needle_lib.serialize_arg(arg_original, "round_trip_1.arg")
    arg_deserialized = arg_needle_lib.deserialize_arg("round_trip_1.arg")
    arg_needle_lib.serialize_arg(arg_deserialized, "round_trip_2.arg")

    assert arg_needle_lib.validate_serialized_arg("round_trip_1.arg")
    assert arg_needle_lib.validate_serialized_arg("round_trip_2.arg")

    # Check that both ARG objects in memory are the same
    assert arg_original.num_nodes() == arg_deserialized.num_nodes()
    assert arg_original.num_edges() == arg_deserialized.num_edges()

    for node_idx in range(arg_original.num_nodes()):
        node_old = arg_original.node(node_idx)
        node_new = arg_deserialized.node(node_idx)
        assert node_old.height == node_new.height
        assert node_old.start == node_new.start
        assert node_old.end == node_new.end
        assert arg_original.is_leaf(node_idx) == arg_deserialized.is_leaf(node_idx)
        assert len(node_old.parent_edges()) == len(node_new.parent_edges())

        for e1, e2 in zip(node_old.parent_edges(), node_new.parent_edges()):
            assert e1.start == e2.start
            assert e1.end == e2.end
            assert e1.parent.ID == e2.parent.ID
            assert e1.child.ID == e2.child.ID

    # Check that both ARG objects on disk are the same
    f1 = h5py.File("round_trip_1.arg", "r")
    f2 = h5py.File("round_trip_2.arg", "r")

    for key in f1.attrs.keys():
        if key != "datetime_created":
            assert f1.attrs[key] == f2.attrs[key]

    assert np.allclose(f1['flags'][...], f2['flags'][...])
    assert np.allclose(f1['times'][...], f2['times'][...])
    assert np.allclose(f1['edge_ids'][...], f2['edge_ids'][...])
    assert np.allclose(f1['edge_ranges'][...], f2['edge_ranges'][...])

    f1.close()
    f2.close()

    # Tidy up files on disk
    pathlib.Path("round_trip_1.arg").unlink()
    pathlib.Path("round_trip_2.arg").unlink()


def test_round_trips_mutations():
    """
    Create an ARG in memory, then serialize -> deserialize -> serialize, and check that the two copies in memory match,
    and the two copies on disk match, specifically with respect to mutations.

    This test demonstrates that the serialization/deserialization is a round-trip process for mutations.
    """
    ts = msprime.simulate(sample_size=20, Ne=1e4, length=5e3, recombination_rate=2e-8, mutation_rate=2e-8,
                          random_seed=1234)

    arg_original = arg_needle_lib.tskit_to_arg(ts)
    arg_needle_lib.generate_m_mutations(arg_original, 89)

    # Check several chunk sizes to make sure the chunking works
    for chunk_size in [1, 7, 100, 1000]:

        arg_needle_lib.serialize_arg(arg_original, "round_trip_1.arg")
        arg_deserialized = arg_needle_lib.deserialize_arg("round_trip_1.arg", chunk_size=chunk_size)
        arg_needle_lib.serialize_arg(arg_deserialized, "round_trip_2.arg")

        assert arg_needle_lib.validate_serialized_arg("round_trip_1.arg")
        assert arg_needle_lib.validate_serialized_arg("round_trip_2.arg")

        # Check that both ARG objects in memory are the same
        assert arg_original.num_mutations() == arg_deserialized.num_mutations()

        for m1, m2 in zip(arg_original.mutations(), arg_deserialized.mutations()):
            assert(m1.position == m2.position)
            assert(m1.height == m2.height)
            assert(m1.site_id == m2.site_id)
            assert(m1.edge.child.ID == m2.edge.child.ID)
            assert(m1.edge.parent.ID == m2.edge.parent.ID)

        # Check that both ARG objects on disk are the same
        f1 = h5py.File("round_trip_1.arg", "r")
        f2 = h5py.File("round_trip_2.arg", "r")

        for key in f1.attrs.keys():
            if key != "datetime_created":
                assert f1.attrs[key] == f2.attrs[key]

        assert np.allclose(f1['mutations']['positions'][...], f2['mutations']['positions'][...])
        assert np.allclose(f1['mutations']['heights'][...], f2['mutations']['heights'][...])
        assert np.allclose(f1['mutations']['site_ids'][...], f2['mutations']['site_ids'][...])
        assert np.allclose(f1['mutations']['edge_ids'][...], f2['mutations']['edge_ids'][...])

        assert np.allclose(f1['flags'][...], f2['flags'][...])
        assert np.allclose(f1['times'][...], f2['times'][...])
        assert np.allclose(f1['edge_ids'][...], f2['edge_ids'][...])
        assert np.allclose(f1['edge_ranges'][...], f2['edge_ranges'][...])

        f1.close()
        f2.close()

        # Tidy up files on disk
        pathlib.Path("round_trip_1.arg").unlink()
        pathlib.Path("round_trip_2.arg").unlink()


def test_validation():
    """
    Test validation of files. This method can diagnose the following problems:
     - the provided file does not exist
     - the provided file is not a valid HDF5 file
     - the provided file is HDF5 but does not contain the expected `arg_file_version` attribute
     - the provided file is an incompatible file version
     - the provided file does not contain the expected attributes or datasets
    """

    # Create a valid ARG and serialize it
    ts = msprime.simulate(sample_size=20, Ne=1e4, length=5e3, recombination_rate=2e-8, mutation_rate=2e-8,
                          random_seed=1234)
    arg = arg_needle_lib.tskit_to_arg(ts)
    arg_needle_lib.serialize_arg(arg, 'valid.arg')

    assert arg_needle_lib.validate_serialized_arg('valid.arg')

    # Check non-existent file
    assert not arg_needle_lib.validate_serialized_arg('file_that_does_not_exist')

    # Valid file that isn't an HDF5 file
    with open('valid_file_invalid_hdf5', 'w') as f:
        f.write('something')
    assert not arg_needle_lib.validate_serialized_arg('valid_file_invalid_hdf5')

    # Check valid hdf5 file that isn't an ARG file
    with h5py.File('valid_hdf5_invalid_arg', 'w') as f:
        f.attrs['something'] = 5
    assert not arg_needle_lib.validate_serialized_arg('valid_hdf5_invalid_arg')

    # Tidy up files on disk
    # pathlib.Path('valid.arg').unlink()
    # pathlib.Path('valid_file_invalid_hdf5').unlink()
    # pathlib.Path('valid_hdf5_invalid_arg').unlink()
