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

#include <vector>

using std::vector;

typedef vector<arg_real_t> rvec;
typedef vector<int> ivec;

TEST_CASE("Threading", "[thread]") {
  SECTION("two samples") {
    ARG arg = ARG(0, 100, 0);
    arg.add_sample();
    arg.add_sample();
    arg.thread_sample(rvec{0, 50}, ivec{0, 0}, rvec{3.14, 2.718});
    arg.populate_children_and_roots();

    REQUIRE(arg.root_starts().size() == 2);
    REQUIRE(arg.root_at(0)->node->height == static_cast<arg_real_t>(3.14));
    REQUIRE(arg.root_at(50)->node->height == static_cast<arg_real_t>(2.718));
    REQUIRE(arg.mrca(0, 1, 0)->height == static_cast<arg_real_t>(3.14));
    REQUIRE(arg.mrca(0, 1, 50)->height == static_cast<arg_real_t>(2.718));
  }
  SECTION("single root but multiple trees") {
    ARG arg = ARG(0, 100, 3);
    arg.add_sample("first");
    arg.add_sample("second");
    arg.thread_sample(rvec{0}, ivec{0}, rvec{10});
    arg.add_sample("third");
    arg.thread_sample(rvec{0, 30, 70}, ivec{1, 1, 1}, rvec{3, 5, 4});
    arg.populate_children_and_roots();

    REQUIRE(arg.root_starts().size() == 1);
    REQUIRE(arg.root_at(50)->node->height == 10);
    REQUIRE(arg.mrca(0, 1, 50)->height == 10);
    REQUIRE(arg.mrca(1, 2, 15)->height == 3);
    REQUIRE(arg.mrca(1, 2, 55)->height == 5);
    REQUIRE(arg.mrca(1, 2, 85)->height == 4);
  }
  SECTION("three samples") {
    ARG arg = ARG(0, 100, 3);
    arg.add_sample("first");
    arg.add_sample("second");
    arg.thread_sample(rvec{0, 50}, ivec{0, 0}, rvec{3.14, 2.718});
    arg.add_sample("third");
    arg.thread_sample(rvec{0, 40}, ivec{1, 0}, rvec{5, 3});
    arg.populate_children_and_roots();

    REQUIRE(arg.root_starts().size() == 3);
    REQUIRE(arg.root_at(0)->node->height == 5);
    REQUIRE(arg.root_at(40)->node->height == static_cast<arg_real_t>(3.14));
    REQUIRE(arg.root_at(50)->node->height == 3);
    REQUIRE(arg.root_at(0)->node->parents.empty());
    REQUIRE(!arg.root_at(40)->node->parents.empty()); // has parents
    REQUIRE(!arg.root_at(50)->node->parents.empty()); // has parents
    REQUIRE(arg.mrca(0, 1, 0)->height == static_cast<arg_real_t>(3.14));
    REQUIRE(arg.mrca(0, 1, 50)->height == static_cast<arg_real_t>(2.718));
  }
  SECTION("three samples same location") {
    ARG arg = ARG(0, 100, 3);
    arg.add_sample("first");
    arg.add_sample("second");
    arg.thread_sample(rvec{0, 50}, ivec{0, 0}, rvec{3.14, 2.718});
    arg.add_sample("third");
    arg.thread_sample(rvec{0, 50}, ivec{1, 1}, rvec{3, 3});
    arg.populate_children_and_roots();

    REQUIRE(arg.root_starts().size() == 2);
    REQUIRE(arg.root_at(0)->node->height == static_cast<arg_real_t>(3.14));
    REQUIRE(arg.root_at(50)->node->height == 3);
    REQUIRE(arg.root_at(0)->node->parents.empty());
    REQUIRE(arg.root_at(50)->node->parents.empty());

    REQUIRE(arg.mrca(1, 2, 0)->height == 3);
    REQUIRE(arg.mrca(0, 2, 0)->height == static_cast<arg_real_t>(3.14));
    REQUIRE(arg.mrca(0, 1, 50)->height == static_cast<arg_real_t>(2.718));
    REQUIRE(arg.mrca(0, 2, 50)->height == 3);
  }
}
