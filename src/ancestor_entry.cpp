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

#include "ancestor_entry.hpp"

#include <cassert>

#include "arg_node.hpp"

AncestorEntry::AncestorEntry(ARGNode* _node, ARGNode* _ancestor, arg_real_t _start, arg_real_t _end)
    : node(_node), ancestor(_ancestor), start(_start), end(_end) {
  assert(start < end);
  assert(start >= ancestor->start);
  assert(start >= node->start);
  assert(end <= ancestor->end);
  assert(end <= node->end);
}
