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

#include "arg.hpp"
#include "types.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <vector>

using Catch::Matchers::ContainsSubstring;

using rvec = std::vector<arg_real_t>;
using ivec = std::vector<int>;

TEST_CASE("Threading vs. adding order", "[thread]") {
  ARG arg = ARG(0, 100);
  REQUIRE_THROWS_WITH(arg.thread_sample(rvec{0}, ivec{0}, rvec{3.14}),
                      ContainsSubstring("No samples to be threaded."));
  arg.add_sample();
  REQUIRE_THROWS_WITH(arg.thread_sample(rvec{0}, ivec{0}, rvec{3.14}),
                      ContainsSubstring("No samples to be threaded."));
  arg.add_sample();
  arg.thread_sample(rvec{0}, ivec{0}, rvec{3.14});
  arg.add_sample();
  arg.thread_sample(rvec{0}, ivec{0}, rvec{2.718});
  arg.add_sample();
  REQUIRE_THROWS_WITH(arg.add_sample(), ContainsSubstring("First thread the last added sample."));
}

TEST_CASE("Threading arguments", "[thread-args]") {
  ARG arg = ARG(0, 100);
  std::vector<int> ids;
  arg.add_sample();
  arg.add_sample();

  SECTION("Sizes") {
    std::string match = "Threading vectors must be of the same nonzero size";
    REQUIRE_THROWS_WITH(
        arg.thread_sample(rvec{0, 5}, ivec{0}, rvec{3.14}), ContainsSubstring(match));
    REQUIRE_THROWS_WITH(arg.thread_sample(rvec(), ivec{0}, rvec{3.14}), ContainsSubstring(match));
    REQUIRE_THROWS_WITH(
        arg.thread_sample(rvec{0, 5}, ivec{0, 0}, rvec{3.14}), ContainsSubstring(match));
    REQUIRE_THROWS_WITH(arg.thread_sample(rvec(), ivec(), rvec()), ContainsSubstring(match));
  }
  SECTION("Sample IDs") {
    std::string match = "Threading IDs must specify samples in ARG";
    REQUIRE_THROWS_WITH(arg.thread_sample(rvec{0}, ivec{1}, rvec{3.14}), ContainsSubstring(match));
    REQUIRE_THROWS_WITH(
        arg.thread_sample(rvec{0}, ivec{200}, rvec{3.14}), ContainsSubstring(match));
    REQUIRE_THROWS_WITH(
        arg.thread_sample(rvec{0, 5}, ivec{1, 1}, rvec{3.14, 3.14}), ContainsSubstring(match));
    REQUIRE_THROWS_WITH(
        arg.thread_sample(rvec{0, 5}, ivec{0, 1}, rvec{3.14, 3.14}), ContainsSubstring(match));
  }
  SECTION("Starts 1") {
    std::string match = "Section starts must partition the region [start, end)";
    REQUIRE_THROWS_WITH(arg.thread_sample(rvec{-5}, ivec{0}, rvec{3.14}), ContainsSubstring(match));
    REQUIRE_THROWS_WITH(
        arg.thread_sample(rvec{0, 105}, ivec{0, 0}, rvec{3.14, 3.14}), ContainsSubstring(match));
    REQUIRE_THROWS_WITH(
        arg.thread_sample(rvec{105}, ivec{0}, rvec{3.14}), ContainsSubstring(match));
  }
  SECTION("Starts 2") {
    std::string match = "Section starts must be in increasing order";
    REQUIRE_THROWS_WITH(
        arg.thread_sample(rvec{0, 0}, ivec{0, 0}, rvec{3.14, 3.14}), ContainsSubstring(match));
    REQUIRE_THROWS_WITH(arg.thread_sample(rvec{0, 8, 7}, ivec{0, 0, 0}, rvec{3.14, 3.14, 3.14}),
                        ContainsSubstring(match));
  }
  SECTION("Heights") {
    std::string match = "Threading heights (times) must be positive";
    REQUIRE_THROWS_WITH(
        arg.thread_sample(rvec{0, 50}, ivec{0, 0}, rvec{3.14, -3.14}), ContainsSubstring(match));
    REQUIRE_THROWS_WITH(arg.thread_sample(rvec{0, 30, 60}, ivec{0, 0, 0}, rvec{3.14, 0, 3.14}),
                        ContainsSubstring(match));
  }
}

TEST_CASE("Sample IDs and names", "[thread]") {
  size_t num_reserved;
  SECTION("Reserving 2 IDs") {
    num_reserved = 2;
  }
  SECTION("Reserving 3 IDs") {
    num_reserved = 3;
  }
  // Note: if we only reserved 1 ID, we would still have the second ID as 1,
  // so we shouldn't add a test case for that

  ARG arg = ARG(0, 100, num_reserved);
  std::vector<int> ids;
  std::vector<std::string> names{"aaa", "bbb", "ccc", "ddd"};
  ids.push_back(arg.add_sample(names[0]));
  ids.push_back(arg.add_sample(names[1]));
  arg.thread_sample(rvec{0}, ivec{0}, rvec{3.14});
  ids.push_back(arg.add_sample(names[2]));
  arg.thread_sample(rvec{0}, ivec{0}, rvec{2.718});
  ids.push_back(arg.add_sample(names[3]));
  arg.thread_sample(rvec{0}, ivec{0}, rvec{5});
  for (size_t i = 0; i < ids.size(); ++i) {
    REQUIRE(arg.is_leaf(ids[i]));
    REQUIRE(arg.sample_names.at(ids[i]) == names[i]);
    if (i < num_reserved) {
      REQUIRE(i == ids[i]);
    }
    else {
      REQUIRE(i != ids[i]);
    }
  }
}
