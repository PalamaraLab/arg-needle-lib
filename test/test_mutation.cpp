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

#include "arg_node.hpp"
#include "arg_edge.hpp"
#include "mutation.hpp"
#include "types.hpp"

#include <catch2/catch_test_macros.hpp>

TEST_CASE("Creating a default mutation", "[mutation]")
{
  Mutation mut(nullptr, arg_real_t{1.23});

  REQUIRE(mut.edge == nullptr);
  REQUIRE(mut.position == arg_real_t{1.23});
  REQUIRE(mut.height == arg_real_t{-1.0});
  REQUIRE(mut.site_id == -1);
}

TEST_CASE("Creating a non-default mutation", "[mutation]")
{
  const auto p_child = std::make_unique<ARGNode>(0, arg_real_t{1.23}, arg_real_t{3.45}, arg_real_t{5.67});
  const auto p_parent = std::make_unique<ARGNode>(1, arg_real_t{2.34}, arg_real_t{4.56}, arg_real_t{6.78});

  const auto p_edge = std::make_unique<ARGEdge>(arg_real_t{5.01}, arg_real_t{5.55}, p_child.get(), p_parent.get());

  const Mutation mut(p_edge.get(), arg_real_t{5.47},arg_real_t{2.02}, 4);

  REQUIRE(mut.edge != nullptr);
  REQUIRE(mut.position == arg_real_t{5.47});
  REQUIRE(mut.height == arg_real_t{2.02});
  REQUIRE(mut.site_id == 4);
  REQUIRE(mut.get_midpoint_height() == arg_real_t{1.785});
}

TEST_CASE("Get midpoint height", "[mutation]")
{
  const auto p_child = std::make_unique<ARGNode>(0, arg_real_t{1.23}, arg_real_t{3.45}, arg_real_t{5.67});
  const auto p_parent = std::make_unique<ARGNode>(1, arg_real_t{2.34}, arg_real_t{4.56}, arg_real_t{6.78});

  const auto p_edge = std::make_unique<ARGEdge>(arg_real_t{5.01}, arg_real_t{5.55}, p_child.get(), p_parent.get());

  const Mutation mut_nonnull(p_edge.get(), arg_real_t{5.47},arg_real_t{2.02}, 4);
  const Mutation mut_null(nullptr, arg_real_t{0.1});

  REQUIRE(mut_nonnull.get_midpoint_height() == arg_real_t{1.785});
  REQUIRE(mut_null.get_midpoint_height() == arg_real_t{-1.0});
}
