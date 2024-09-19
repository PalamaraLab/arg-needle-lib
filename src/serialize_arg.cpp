/*
  This file is part of the ARG-Needle genealogical inference and
  analysis software suite.
  Copyright (C) 2023-2024 ARG-Needle Developers.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "serialize_arg.hpp"

#include "deserialization_params.hpp"

#include "H5Cpp.h"

#include <iostream>
#include <filesystem>
#include <utility>
#include <vector>

namespace
{

bool check_attribute(const H5::H5File& h5_file, const std::string& expected_attr)
{
  if (!h5_file.attrExists(expected_attr)) {
    std::cerr << "Expected file " << h5_file.getFileName() << " to include attribute `" << expected_attr << "`"
              << std::endl;
    return false;
  }
  return true;
}

bool check_dataset(const H5::H5File& h5_file, const std::string& expected_dset)
{
  try {
    auto dset = h5_file.openDataSet(expected_dset);
    return true;
  } catch (const H5::Exception&) {
    std::cerr << "Expected file " << h5_file.getFileName() << " to include dataset `" << expected_dset << "`"
              << std::endl;
    return false;
  }
}

bool check_group(const H5::H5File& h5_file, const std::string& expected_group)
{
  try {
    auto group = h5_file.openGroup(expected_group);
    return true;
  } catch (const H5::Exception&) {
    std::cerr << "Expected file " << h5_file.getFileName() << " to include group `" << expected_group << "`"
              << std::endl;
    return false;
  }
}

bool read_bool_attribute(const H5::H5File& file, const std::string& attrName)
{
  bool value = false;

  try {
    // Open the attribute from the file
    const H5::Attribute attribute = file.openAttribute(attrName);

    // Read the attribute, assuming it was stored as H5T_NATIVE_HBOOL
    uint8_t buffer{};
    attribute.read(attribute.getDataType(), &buffer);

    // Convert the integer buffer to bool
    value = (buffer != 0u);
  } catch (const H5::Exception& e) {
    // Handle exceptions: attribute not found, wrong type, etc.
    std::cerr << "Error reading attribute `" << attrName << "`: " << e.getDetailMsg() << std::endl;
  }

  return value;
}

template <class T> struct dependent_false : std::false_type {
};

template <class T> inline constexpr bool dependent_false_v = dependent_false<T>::value;

template <typename T> H5::PredType getPredType()
{
  if constexpr (std::is_same_v<T, double>) {
    return H5::PredType::NATIVE_DOUBLE;
  } else if constexpr (std::is_same_v<T, int>) {
    return H5::PredType::NATIVE_INT;
  } else if constexpr (std::is_same_v<T, uint8_t>) {
    return H5::PredType::NATIVE_UINT8;
  } else {
    static_assert(dependent_false_v<T>, "Unsupported type");
  }
}

// Templated function to read a dataset into a std::vector<T>
template<typename T>
std::vector<T> read_dataset_to_vector_1d(const H5::H5File& h5_file, const std::string& dset_name, hssize_t start = 0, hssize_t stop = -1) {
  std::vector<T> data;

  try {
    H5::DataSet dataset = h5_file.openDataSet(dset_name);
    H5::DataSpace dataspace = dataset.getSpace();

    if (dataspace.getSimpleExtentNdims() != 1) {
      throw std::runtime_error("Dataset must be 1-dimensional");
    }

    hsize_t dims[1];
    dataspace.getSimpleExtentDims(dims);
    const hsize_t num_elements = dims[0];

    // Adjust stop value if necessary
    if (stop == -1 || stop > num_elements) {
      stop = num_elements;
    }

    if (start >= stop) {
      throw std::runtime_error("Invalid range: start must be less than stop");
    }

    const hsize_t block_size = stop - start;
    data.resize(block_size);

    // Define hyperslab in the dataset
    hsize_t offset[1] = {static_cast<hsize_t>(start)};
    hsize_t count[1] = {block_size};
    dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);

    // Define the memory dataspace
    H5::DataSpace memspace(1, count);
    dataset.read(data.data(), getPredType<T>(), memspace, dataspace);
  } catch (const H5::Exception& e) {
    throw std::runtime_error("Failed to read dataset `" + dset_name + "`: " + e.getDetailMsg());
  }

  return data;
}

template <typename T>
std::vector<std::array<T, 2>> read_dataset_to_vector_2d(
    const H5::H5File& h5_file, const std::string& dset_name, hssize_t start = 0, hssize_t stop = -1)
{
  std::vector<std::array<T, 2>> data;

  try {
    H5::DataSet dataset = h5_file.openDataSet(dset_name);
    H5::DataSpace dataspace = dataset.getSpace();

    if (dataspace.getSimpleExtentNdims() != 2) {
      throw std::runtime_error("Dataset must be 2-dimensional");
    }

    hsize_t dims[2];
    dataspace.getSimpleExtentDims(dims);
    if (dims[1] != 2) {
      throw std::runtime_error("Second dimension of the dataset must be 2");
    }

    hsize_t num_elements = dims[0];

    // Adjust stop value if necessary
    if (stop == -1 || stop > num_elements) {
      stop = num_elements;
    }

    if (start >= stop) {
      throw std::runtime_error("Invalid range: start must be less than stop");
    }

    hsize_t block_size = stop - start;
    data.resize(block_size);

    // Define hyperslab in the dataset
    hsize_t offset[2] = {static_cast<hsize_t>(start), 0}; // Start from 'start', covering the entire 2nd dimension
    hsize_t count[2] = {block_size, 2}; // Number of elements to read in the 1st dim and include all of the 2nd dim
    dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);

    // Define the memory dataspace
    H5::DataSpace memspace(2, count);
    dataset.read(data.data(), getPredType<T>(), memspace, dataspace);
  } catch (const H5::Exception& e) {
    throw std::runtime_error("Failed to read dataset `" + dset_name + "`: " + e.getDetailMsg());
  }

  return data;
}

int read_int_attribute(const H5::H5File& file, const std::string& attrName)
{
  int value{};

  try {
    // Open the attribute from the file
    H5::Attribute attribute = file.openAttribute(attrName);

    // Read the attribute, assuming it was stored as H5T_NATIVE_INT
    attribute.read(H5::PredType::NATIVE_INT, &value);
  } catch (const H5::Exception& e) {
    // Handle exceptions: attribute not found, wrong type, etc.
    std::cerr << "Error reading attribute `" << attrName << "`: " << e.getDetailMsg() << std::endl;
    exit(0);
  }

  return value;
}

bool validate_serialized_arg_v1(const H5::H5File& h5_file)
{
  // Expected attributes and datasets
  std::vector<std::string> expected_attrs = {"num_nodes", "num_edges", "num_mutations", "offset", "chromosome",
      "sequence_length", "datetime_created", "arg_file_version"};
  std::vector<std::string> expected_dsets = {"flags", "times", "edge_ranges", "edge_ids"};

  bool is_valid = true;

  for (const auto& attr : expected_attrs) {
    is_valid = check_attribute(h5_file, attr);
  }

  for (const auto& dset : expected_dsets) {
    is_valid = check_dataset(h5_file, dset);
  }

  return is_valid;
}

bool validate_serialized_arg_v2(const H5::H5File& h5_file)
{
  // Expected attributes and datasets
  std::vector<std::string> expected_attrs = {"num_nodes", "num_edges", "node_bounds", "num_mutations", "mutations",
      "offset", "chromosome", "start", "end", "threaded_samples", "datetime_created", "arg_file_version"};

  std::vector<std::string> expected_dsets = {"flags", "times", "edge_ranges", "edge_ids"};
  std::vector<std::string> optional_dsets = {"node_bounds"};

  std::vector<std::string> expected_groups = {};
  std::vector<std::string> optional_groups = {"mutations"};

  bool is_valid = true;

  for (const auto& attr : expected_attrs) {
    is_valid = check_attribute(h5_file, attr);
  }

  // The existence of optional datasets is marked by bool attributes of the same name
  for (const auto& dset_name : optional_dsets) {
    if (read_bool_attribute(h5_file, dset_name)) {
      expected_dsets.emplace_back(dset_name);
    }
  }

  for (const auto& dset : expected_dsets) {
    is_valid = check_dataset(h5_file, dset);
  }

  // The existence of optional groups is also marked by bool attributes of the same name
  for (const auto& group_name : optional_groups) {
    if (read_bool_attribute(h5_file, group_name)) {
      expected_groups.emplace_back(group_name);
    }
  }

  for (const auto& group : expected_groups) {
    is_valid = check_group(h5_file, group);
  }

  return is_valid;
}

ARG deserialize_arg_v1(const H5::H5File& h5_file, const int reserved_samples)
{

  // Read attributes
  const int offset = read_int_attribute(h5_file, "offset");
  const int chromosome = read_int_attribute(h5_file, "chromosome");
  const int sequence_length = read_int_attribute(h5_file, "sequence_length");
  const int num_nodes = read_int_attribute(h5_file, "num_nodes");
  const int num_edges = read_int_attribute(h5_file, "num_edges");

  // Read datasets
  std::vector<uint8_t> raw_flags = read_dataset_to_vector_1d<uint8_t>(h5_file, "flags");
  std::vector<double> raw_times = read_dataset_to_vector_1d<double>(h5_file, "times");
  std::vector<int> raw_edge_ids = read_dataset_to_vector_1d<int>(h5_file, "edge_ids");
  std::vector<double> raw_edge_ranges = read_dataset_to_vector_1d<double>(h5_file, "edge_ranges");

  assert(raw_flags.size() == num_nodes);
  assert(raw_times.size() == num_nodes);
  assert(raw_edge_ids.size() == 2 * num_edges);
  assert(raw_edge_ranges.size() == 2 * num_edges);

  // Translate from raw data to required data structures
  std::deque<bool> is_sample;
  for (const auto flag : raw_flags) {
    is_sample.push_back(flag != 0);
  }

  std::vector<std::pair<int, int>> edge_ids;
  for (std::size_t i = 0; i < raw_edge_ids.size(); i += 2) {
    edge_ids.emplace_back(raw_edge_ids[i], raw_edge_ids[i + 1]);
  }

  std::vector<std::pair<double, double>> edge_ranges;
  for (std::size_t i = 0; i < raw_edge_ranges.size(); i += 2) {
    edge_ranges.emplace_back(raw_edge_ranges[i], raw_edge_ranges[i + 1]);
  }

  // Construct the ARG object
  ARG arg(0, sequence_length, raw_times, is_sample, edge_ids, edge_ranges, reserved_samples);
  arg.set_offset(offset);
  arg.set_chromosome(chromosome);

  return arg;
}

ARG deserialize_arg_v2(const H5::H5File& h5_file, const int chunk_size, const int reserved_samples)
{
  DeserializationParams dp;
  dp.start = read_int_attribute(h5_file, "start");
  dp.end = read_int_attribute(h5_file, "end");
  dp.num_nodes = read_int_attribute(h5_file, "num_nodes");
  dp.offset = read_int_attribute(h5_file, "offset");
  dp.chromosome = read_int_attribute(h5_file, "chromosome");
  dp.threaded_samples = read_int_attribute(h5_file, "threaded_samples");
  dp.reserved_samples = reserved_samples;

  ARG arg(dp);

  // Process {chunk_size} nodes at a time, adding each chunk to the ARG as we go
  {
    const auto num_nodes = static_cast<hssize_t>(dp.num_nodes);
    hssize_t num_nodes_written = 0;

    while (num_nodes_written < num_nodes) {
      const hssize_t range_lo = num_nodes_written;
      const hssize_t range_hi = std::min(num_nodes_written + chunk_size, num_nodes);

      const auto node_heights = read_dataset_to_vector_1d<double>(h5_file, "times", range_lo, range_hi);
      const auto is_sample = read_dataset_to_vector_1d<uint8_t>(h5_file, "flags", range_lo, range_hi);

      if (read_bool_attribute(h5_file, "node_bounds")) {
        const auto node_bounds_data = read_dataset_to_vector_2d<double>(h5_file, "node_bounds", range_lo, range_hi);
        arg.deserialize_add_nodes(node_heights, is_sample, node_bounds_data);
      } else {
        arg.deserialize_add_nodes(node_heights, is_sample);
      }

      num_nodes_written += static_cast<hssize_t>(node_heights.size());
    }
  }

  // Process {chunk_size} edges at a time, adding each chunk to the ARG as we go
  {
    const auto num_edges = static_cast<hssize_t>(read_int_attribute(h5_file, "num_edges"));
    hssize_t num_edges_written = 0;

    while (num_edges_written < num_edges) {

      const hssize_t range_lo = num_edges_written;
      const hssize_t range_hi = std::min(num_edges_written + chunk_size, num_edges);

      const auto edge_id_data = read_dataset_to_vector_2d<int>(h5_file, "edge_ids", range_lo, range_hi);
      const auto edge_range_data = read_dataset_to_vector_2d<double>(h5_file, "edge_ranges", range_lo, range_hi);

      arg.deserialize_add_edges(edge_id_data, edge_range_data);

      num_edges_written += static_cast<hssize_t>(edge_id_data.size());
    }
  }

  // Process {chunk_size} mutations at a time, adding each chunk to the ARG as we go
  if (read_bool_attribute(h5_file, "mutations")) {

    const auto num_mutations = static_cast<hssize_t>(read_int_attribute(h5_file, "num_mutations"));
    hssize_t num_mutations_written = 0;

    while (num_mutations_written < num_mutations) {

      const hssize_t range_lo = num_mutations_written;
      const hssize_t range_hi = std::min(num_mutations_written + chunk_size, num_mutations);

      const auto mut_pos = read_dataset_to_vector_1d<double>(h5_file, "mutations/positions", range_lo, range_hi);
      const auto mut_hts = read_dataset_to_vector_1d<double>(h5_file, "mutations/heights", range_lo, range_hi);
      const auto mut_sid = read_dataset_to_vector_1d<int>(h5_file, "mutations/site_ids", range_lo, range_hi);
      const auto mut_eid = read_dataset_to_vector_2d<int>(h5_file, "mutations/edge_ids", range_lo, range_hi);

      arg.deserialize_add_mutations(mut_pos, mut_hts, mut_sid, mut_eid);

      num_mutations_written += static_cast<hssize_t>(mut_pos.size());
    }
  }

  return arg;
}

} // namespace

bool arg_utils::validate_serialized_arg(const std::string& file_name)
{
  // Check if file exists
  if (!std::filesystem::exists(file_name)) {
    std::cout << "File: " << file_name << " is not a valid file" << std::endl;
    return false;
  }

  if (!H5Fis_hdf5(file_name.c_str()))
  {
    std::cout << "File: " << file_name << " is not a valid HDF5 file" << std::endl;
    return false;
  }

  try {
    H5::H5File h5_file(file_name, H5F_ACC_RDONLY);

    // Check for the 'arg_file_version' attribute
    if (!check_attribute(h5_file, "arg_file_version")) {
      std::cout << "File: " << file_name
                << " is not a valid arg file because it does not contain `arg_file_version` attribute" << std::endl;
      return false;
    }

    const int arg_file_version = read_int_attribute(h5_file, "arg_file_version");

    // Validate file version
    if (arg_file_version == 1) {
      return validate_serialized_arg_v1(h5_file);
    }
    if (arg_file_version == 2) {
      return validate_serialized_arg_v2(h5_file);
    }

    std::cout << "Arg file version (" << arg_file_version << ") is not supported; valid versions are 1, 2."
              << std::endl;
    return false;

  } catch (const H5::Exception& e) {
    std::cerr << "HDF5 error on file: " << file_name << std::endl;
    std::cerr << e.getDetailMsg() << std::endl;
    return false;
  }
}

ARG arg_utils::deserialize_arg(const std::string& file_name, const int chunk_size, const int reserved_samples)
{
  if (!validate_serialized_arg(file_name)) {
    throw std::runtime_error("Invalid ARG file: " + file_name);
  }

  try {
    H5::H5File h5_file(file_name, H5F_ACC_RDONLY);

    const int arg_file_version = read_int_attribute(h5_file, "arg_file_version");

    if (arg_file_version == 1) {
      return deserialize_arg_v1(h5_file, reserved_samples);
    }
    if (arg_file_version == 2) {
      return deserialize_arg_v2(h5_file, chunk_size, reserved_samples);
    }

    throw std::logic_error(
        "Reached an unsupported arg_file_version after validation: " + std::to_string(arg_file_version));

  } catch (const H5::Exception& e) {
    std::cerr << "HDF5 error on file: " << file_name << std::endl;
    std::cerr << e.getDetailMsg() << std::endl;
    throw std::runtime_error("Unable to deserialize arg file: " + file_name);
  }
}
