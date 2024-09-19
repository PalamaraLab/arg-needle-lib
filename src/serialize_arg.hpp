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

#ifndef ARG_NEEDLE_LIB_SERIALIZE_ARG_H
#define ARG_NEEDLE_LIB_SERIALIZE_ARG_H

#include "arg.hpp"

#include <string>

namespace arg_utils
{

/**
 * @brief Validates the integrity of a serialized ARG file.
 *
 * This function checks a serialized ARG file for the presence of expected attributes and datasets, ensuring it conforms
 * to the specified version's format. It supports validation for multiple versions of the serialized ARG format by
 * calling version-specific validation functions. The function first checks for the file's existence and whether it is a
 * valid HDF5 file. It then reads the `arg_file_version` attribute to determine the file's format version and proceeds
 * with the appropriate validation routine.
 *
 * @param file_name The file path of the serialized ARG file to be validated.
 *
 * @return Returns `true` if the file is valid according to the specified version's format, otherwise `false`.
 */
bool validate_serialized_arg(const std::string& file_name);


ARG deserialize_arg(const std::string& file_name, int chunk_size = 1000, int reserved_samples = -1);
} // namespace arg_utils

#endif // ARG_NEEDLE_LIB_SERIALIZE_ARG_H
