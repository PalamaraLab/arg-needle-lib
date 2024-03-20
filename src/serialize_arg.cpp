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

extern "C" {
#include "hdf5.h"
}

#include <iostream>
#include <filesystem>
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
