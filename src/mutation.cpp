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

#include "mutation.hpp"

#include <cassert>

#include "arg_node.hpp"

Mutation::Mutation(ARGEdge* _edge, const arg_real_t _position, const arg_real_t _height, const int _site_id)
    : position(_position), height(_height), edge(_edge), site_id(_site_id) {

  if (edge != nullptr) {
    assert(position < edge->end && position >= edge->start);
    if (height >= 0) {
      assert(height < edge->parent->height);
      assert(height >= edge->child->height);
    }
  }
  // can't check that the site position is correct because that's stored in the ARG
}

arg_real_t Mutation::get_midpoint_height() const
{
  if (edge == nullptr) {
    std::cout << "Warning: ARGEdge* is a nullptr: midpoint height estimate will be -1.0\n";
    return arg_real_t{-1.0};
  }

  const arg_real_t lo = this->edge->child->height;
  const arg_real_t hi = this->edge->parent->height;

  // Safe(ish) midpoint in case these numbers are really large
  return lo + (hi - lo) / 2.0;
}
