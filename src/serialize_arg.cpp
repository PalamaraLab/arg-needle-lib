/*
  This file is part of the ARG-Needle genealogical inference and
  analysis software suite.
  Copyright (C) 2023 ARG-Needle Developers.

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

extern "C" {
#include "hdf5.h"
#include "H5LTpublic.h"
}
#include "H5Cpp.h"

#include <iostream>
#include <filesystem>
#include <utility>
#include <vector>

namespace
{


bool validate_serialized_arg_v1(const std::string& file_name)
{
  // Expected attributes and datasets
  std::vector<std::string> expected_attrs = {"num_nodes", "num_edges", "num_mutations", "offset", "chromosome",
      "sequence_length", "datetime_created", "arg_file_version"};
  std::vector<std::string> expected_dsets = {"flags", "times", "edge_ranges", "edge_ids"};

  // Open the HDF5 file
  const hid_t file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) {
    std::cerr << "Could not open file: " << file_name << std::endl;
    return false;
  }

  bool is_valid = true;

  // Check for the presence of expected attributes
  for (const auto& attr : expected_attrs) {
    if (H5Aexists(file_id, attr.c_str()) <= 0) {
      std::cerr << "Expected file " << file_name << " to include attribute `" << attr << "`" << std::endl;
      is_valid = false;
    }
  }

  // Check for the presence of expected datasets only if attributes are valid
  for (const auto& dset : expected_dsets) {
    if (const hid_t dset_id = H5Dopen2(file_id, dset.c_str(), H5P_DEFAULT); dset_id < 0) {
      std::cerr << "Expected file " << file_name << " to include dataset `" << dset << "`" << std::endl;
      is_valid = false;
    } else {
      H5Dclose(dset_id);
    }
  }

  H5Fclose(file_id);
  return is_valid;
}

bool validate_serialized_arg_v2(const std::string& file_name)
{
  // Expected attributes and datasets
  std::vector<std::string> expected_attrs = {"num_nodes", "num_edges", "node_bounds", "num_mutations", "mutations",
      "offset", "chromosome", "start", "end", "threaded_samples", "datetime_created", "arg_file_version"};
  std::vector<std::string> expected_dsets = {"flags", "times", "edge_ranges", "edge_ids"};
  std::vector<std::string> optional_dsets = {"mutations", "node_bounds"};

  // Open the HDF5 file
  const hid_t file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) {
    std::cerr << "Could not open file: " << file_name << std::endl;
    return false;
  }

  bool is_valid = true;

  // Check for the presence of expected attributes
  for (const auto& attr : expected_attrs) {
    if (H5Aexists(file_id, attr.c_str()) <= 0) {
      std::cerr << "Expected file " << file_name << " to include attribute `" << attr << "`" << std::endl;
      is_valid = false;
    }
  }

  // Check for the presence of optional datasets, appending them to the list of expected datasets if appropriate
  for (const auto& dset_name : optional_dsets) {
    if (const hid_t attr_id = H5Aopen(file_id, dset_name.c_str(), H5P_DEFAULT); attr_id < 0) {
      std::cerr << "Unable to open `" << dset_name << "` attribute in file: " << file_name << std::endl;
      is_valid = false;
    } else {
      int dset_flag = 0;
      H5Aread(attr_id, H5T_NATIVE_INT, &dset_flag);
      H5Aclose(attr_id);
      if (dset_flag) {
        expected_dsets.emplace_back(dset_name);
      }
    }
  }

  // Check that all expected datasets exist
  for (const auto& dset : expected_dsets) {
    const hid_t dset_id = H5Dopen2(file_id, dset.c_str(), H5P_DEFAULT);
    if (dset_id < 0) {
      std::cerr << "Expected file " << file_name << " to include dataset `" << dset << "`" << std::endl;
      is_valid = false;
    } else {
      H5Dclose(dset_id);
    }
  }

  H5Fclose(file_id);
  return is_valid;
}

ARG deserialize_arg_v1(const std::string& file_name, const int reserved_samples)
{
  // Open the HDF5 file
  const hid_t file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Read attributes
  int offset;
  H5LTget_attribute_int(file_id, ".", "offset", &offset);

  int chromosome;
  H5LTget_attribute_int(file_id, ".", "chromosome", &chromosome);

  int sequence_length;
  H5LTget_attribute_int(file_id, ".", "sequence_length", &sequence_length);

  int num_nodes;
  H5LTget_attribute_int(file_id, ".", "num_nodes", &num_nodes);

  int num_edges;
  H5LTget_attribute_int(file_id, ".", "num_edges", &num_edges);

  // Read datasets
  std::vector<uint8_t> raw_flags;
  raw_flags.resize(num_nodes);
  H5LTread_dataset(file_id, "flags", H5T_NATIVE_UCHAR, raw_flags.data());

  std::vector<double> raw_times;
  raw_times.resize(num_nodes);
  H5LTread_dataset(file_id, "times", H5T_NATIVE_DOUBLE, raw_times.data());

  std::vector<int> raw_edge_ids;
  raw_edge_ids.resize(2 * num_edges);
  H5LTread_dataset(file_id, "edge_ids", H5T_NATIVE_INT, raw_edge_ids.data());

  std::vector<double> raw_edge_ranges;
  raw_edge_ranges.resize(2 * num_edges);
  H5LTread_dataset(file_id, "edge_ranges", H5T_NATIVE_DOUBLE, raw_edge_ranges.data());

  // Close the HDF5 file
  H5Fclose(file_id);

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

ARG deserialize_arg_v2(const std::string& file_name, const int chunk_size, const int reserved_samples)
{
  H5::H5File file(file_name, H5F_ACC_RDONLY);

  auto ReadIntAttr = [&file](const std::string& attr_name) {
    const H5::Attribute attr = file.openAttribute(attr_name);
    int value;
    attr.read(attr.getIntType(), &value);
    return value;
  };

  DeserializationParams dp;
  dp.start = ReadIntAttr("start");
  dp.end = ReadIntAttr("end");
  dp.num_nodes = ReadIntAttr("num_nodes");
  dp.offset = ReadIntAttr("offset");
  dp.chromosome = ReadIntAttr("chromosome");
  dp.threaded_samples = ReadIntAttr("threaded_samples");
  dp.reserved_samples = reserved_samples;

  ARG arg(dp);

  // Assuming ARG has a constructor that takes DeserializationParams
  const int num_nodes = dp.num_nodes;
  int num_nodes_written = 0;

  // while (num_nodes_written < num_nodes) {
  //   const int range_lo = num_nodes_written;
  //   const int range_hi = std::min(num_nodes_written + chunk_size, num_nodes);
  //   const int range_len = range_hi - range_lo;
  //
  //   std::vector<double> node_heights = file.getDataSet("times").readChunk<double>(range_lo, range_len);
  //   std::vector<bool> is_sample = file.getDataSet("flags").readChunk<bool>(range_lo, range_len);
  //
  //   if (file.getAttr("node_bounds").readBool()) {
  //     auto node_bounds_data = file.getDataSet("node_bounds").readChunk<std::pair<double, double>>(range_lo, range_len);
  //     arg.deserialize_add_nodes(node_heights, is_sample, node_bounds_data);
  //   } else {
  //     arg.deserialize_add_nodes(node_heights, is_sample);
  //   }
  //
  //   num_nodes_written += range_len;
  // }
  //
  // int num_edges = file.getAttr("num_edges").readInt();
  // int num_edges_written = 0;
  // while (num_edges_written < num_edges) {
  //   int range_lo = num_edges_written;
  //   int range_hi = std::min(num_edges_written + chunk_size, num_edges);
  //   int range_len = range_hi - range_lo;
  //
  //   auto edge_ids = file.getDataSet("edge_ids").readChunk<std::pair<int, int>>(range_lo, range_len);
  //   auto edge_ranges = file.getDataSet("edge_ranges").readChunk<std::pair<double, double>>(range_lo, range_len);
  //
  //   arg.deserialize_add_edges(edge_ids, edge_ranges);
  //
  //   num_edges_written += range_len;
  // }
  //
  // if (file.getAttr("mutations").readBool()) {
  //   int num_mutations = file.getAttr("num_mutations").readInt();
  //   int num_mutations_written = 0;
  //   while (num_mutations_written < num_mutations) {
  //     int range_lo = num_mutations_written;
  //     int range_hi = std::min(num_mutations_written + chunk_size, num_mutations);
  //     int range_len = range_hi - range_lo;
  //
  //     auto mut_pos = file.getDataSet("mutations/positions").readChunk<double>(range_lo, range_len);
  //     auto mut_hts = file.getDataSet("mutations/heights").readChunk<double>(range_lo, range_len);
  //     auto mut_sid = file.getDataSet("mutations/site_ids").readChunk<int>(range_lo, range_len);
  //     auto mut_eid = file.getDataSet("mutations/edge_ids").readChunk<std::pair<int, int>>(range_lo, range_len);
  //
  //     arg.deserialize_add_mutations(mut_pos, mut_hts, mut_sid, mut_eid);
  //
  //     num_mutations_written += range_len;
  //   }
  // }

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

  // Open the HDF5 file
  const hid_t file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) {
    std::cout << "File: " << file_name << " is not a valid HDF5 file" << std::endl;
    return false;
  }

  // Check for the 'arg_file_version' attribute
  if (H5Aexists(file_id, "arg_file_version") <= 0) {
    std::cout << "File: " << file_name
              << " is not a valid arg file because it does not contain `arg_file_version` attribute" << std::endl;
    H5Fclose(file_id);
    return false;
  }

  // Read the 'arg_file_version' attribute
  const hid_t attr_id = H5Aopen(file_id, "arg_file_version", H5P_DEFAULT);
  int arg_file_version;
  H5Aread(attr_id, H5T_NATIVE_INT, &arg_file_version);
  H5Aclose(attr_id);
  H5Fclose(file_id);

  // Validate file version
  if (arg_file_version == 1) {
    return validate_serialized_arg_v1(file_name);
  }
  if (arg_file_version == 2) {
    return validate_serialized_arg_v2(file_name);
  }

  std::cout << "Arg file version (" << arg_file_version << ") is not supported; valid versions are 1, 2." << std::endl;
  return false;
}


ARG arg_utils::deserialize_arg(const std::string& file_name, const int chunk_size, const int reserved_samples) {

  if (!validate_serialized_arg(file_name)) {
    throw std::runtime_error("Invalid ARG file: " + file_name);
  }

  // Open the HDF5 file
  const hid_t file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) {
    throw std::runtime_error("Could not open file: " + file_name);
  }

  // Read the 'arg_file_version' attribute
  int arg_file_version = 0;
  if (const hid_t attr_id = H5Aopen(file_id, "arg_file_version", H5P_DEFAULT); attr_id >= 0) {
    H5Aread(attr_id, H5T_NATIVE_INT, &arg_file_version);
    H5Aclose(attr_id);
  } else {
    H5Fclose(file_id);
    throw std::runtime_error("arg_file_version attribute missing in file: " + file_name);
  }

  H5Fclose(file_id);

  // Conditional logic based on the file version
  if (arg_file_version == 1) {
    return deserialize_arg_v1(file_name, reserved_samples);
  }
  if (arg_file_version == 2) {
    return deserialize_arg_v2(file_name, chunk_size, reserved_samples);
  }

  throw std::logic_error("Reached an unsupported arg_file_version after validation: " + std::to_string(arg_file_version));
}
