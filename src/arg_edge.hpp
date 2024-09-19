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

#ifndef __ARG_EDGE_HPP_
#define __ARG_EDGE_HPP_

#include "types.hpp"

#include <iostream>
#include <map>
#include <memory>
#include <vector>

class ARGNode;

class Mutation;

class ARGEdge {

public:
  arg_real_t start, end;
  ARGNode *child, *parent;

  // Mutations are owned by ARG class (as a unique_ptr)
  std::unique_ptr<std::vector<Mutation*>> mutations;

  ARGEdge(arg_real_t _start, arg_real_t _end, ARGNode* _child, ARGNode* _parent);
  ~ARGEdge();
  // We implement an explicit copy constructor to handle the mutations unique_ptr
  explicit ARGEdge(ARGEdge const&);
  // Assigning is not allowed, because we have smart pointers!
  ARGEdge& operator=(ARGEdge const&) = delete;

  void add_mutations(std::vector<Mutation*> new_mutations);
  // These are all linear in the number of mutations, but we
  // expect few mutations per edge.
  void update_start(arg_real_t new_start);
  void update_end(arg_real_t new_end);
  bool mutated_at_site(int site_id);
  void remove_mutations_at_site(int site_id);
  std::vector<Mutation*> mutations_in_range(arg_real_t range_start, arg_real_t range_end);
  friend std::ostream& operator<<(std::ostream& os, const ARGEdge& edge);
};

#endif // __ARG_EDGE_HPP_
