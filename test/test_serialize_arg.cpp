/*
  This file is part of the ARG-Needle genealogical inference and
  analysis software suite.
  Copyright (C) 2023-2025 ARG-Needle Developers.

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

#include "arg.hpp"
#include "serialize_arg.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>

#include <iostream>
#include <sstream>


TEST_CASE("Validate serialized arg files")
{
  SECTION("Valid arg")
  {
    const std::string valid_arg_file = ARG_NEEDLE_TESTDATA_DIR "/test_serialize_arg/valid.arg";
    CHECK(arg_utils::validate_serialized_arg(valid_arg_file) == true);
  }

  SECTION("File that does not exist")
  {
    // Redirect std::cout to a stringstream, so we can test the warnings generated
    std::stringstream buffer;
    std::streambuf* prevCoutBuffer = std::cout.rdbuf(buffer.rdbuf());

    CHECK(arg_utils::validate_serialized_arg("file_that_does_not_exist") == false);


    CHECK(buffer.str() == "File: file_that_does_not_exist is not a valid file\n");
    std::cout.rdbuf(prevCoutBuffer); // Restore original buffer
  }

  SECTION("Valid file, invalid HDF5")
  {
    // Redirect std::cout to a stringstream, so we can test the warnings generated
    std::stringstream buffer;
    std::streambuf* prevCoutBuffer = std::cout.rdbuf(buffer.rdbuf());

    const std::string valid_file_invalid_hdf5 = ARG_NEEDLE_TESTDATA_DIR "/test_serialize_arg/valid_file_invalid_hdf5";

    CHECK(arg_utils::validate_serialized_arg(valid_file_invalid_hdf5) == false);
    CHECK_THAT(buffer.str(), Catch::Matchers::ContainsSubstring("valid_file_invalid_hdf5 is not a valid HDF5 file"));

    std::cout.rdbuf(prevCoutBuffer); // Restore original buffer
  }

  SECTION("Valid HDF5 file, invalid ARG")
  {
    // Redirect std::cout to a stringstream, so we can test the warnings generated
    std::stringstream buffer;
    std::streambuf* prevCoutBuffer = std::cout.rdbuf(buffer.rdbuf());

    const std::string valid_hdf5_invalid_arg = ARG_NEEDLE_TESTDATA_DIR "/test_serialize_arg/valid_hdf5_invalid_arg";

    CHECK(arg_utils::validate_serialized_arg(valid_hdf5_invalid_arg) == false);
    CHECK_THAT(buffer.str(), Catch::Matchers::ContainsSubstring("valid_hdf5_invalid_arg is not a valid arg file"));

    std::cout.rdbuf(prevCoutBuffer); // Restore original buffer
  }

  SECTION("Valid ARG")
  {
    const std::string valid_arg = ARG_NEEDLE_TESTDATA_DIR "/test_serialize_arg/valid.arg";
    CHECK(arg_utils::validate_serialized_arg(valid_arg) == true);
  }
}

TEST_CASE("Deserialize ARG")
{
  const std::string valid_arg = ARG_NEEDLE_TESTDATA_DIR "/test_serialize_arg/valid.arg";
  const ARG from_arg_file = arg_utils::deserialize_arg(valid_arg);

  CHECK(from_arg_file.num_nodes() == 50);
}

TEST_CASE("Deserialize ARG with mutations")
{
  const std::string valid_arg = ARG_NEEDLE_TESTDATA_DIR "/test_serialize_arg/valid_with_mutations.arg";
  const ARG from_arg_file = arg_utils::deserialize_arg(valid_arg);

  CHECK(from_arg_file.num_nodes() == 50);
  CHECK(from_arg_file.num_mutations() == 89);
}
