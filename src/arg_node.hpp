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

#ifndef __ARG_NODE_HPP_
#define __ARG_NODE_HPP_

#include "IntervalTree.h"
#include "types.hpp"

#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

class ARGEdge;

class ARGNode {
public:
  int ID;
  arg_real_t height, start, end;
  std::map<arg_real_t, std::unique_ptr<ARGEdge>> parents;
  // we wrap this in a unique_ptr so it doesn't take up memory unnecessarily
  std::unique_ptr<IntervalTree<arg_real_t, ARGEdge*>> children;

  ARGNode(int _ID, arg_real_t _height, arg_real_t _start, arg_real_t _end);
  ~ARGNode();
  // copying and assigning is not allowed, because we have smart pointers
  ARGNode(ARGNode const&) = delete;
  ARGNode& operator=(ARGNode const&) = delete;
  // Allow move construction: https://stackoverflow.com/a/10473009
  ARGNode(ARGNode&&) = default;
  ARGNode& operator=(ARGNode&&) = default;
  ARGEdge* add_parent(arg_real_t start, arg_real_t end, ARGNode& parent);
  void remove_parent(arg_real_t start);
  void update_parent_start(arg_real_t start_old, arg_real_t start_new);
  void update_parent_end(arg_real_t start, arg_real_t end_new);
  // Get all ARGEdges that overlap a position
  std::vector<ARGEdge*> children_at(arg_real_t position) const;
  std::vector<ARGEdge*> children_overlap(arg_real_t pos_start, arg_real_t pos_end) const;
  ARGEdge* parent_edge_at(arg_real_t position) const;
  friend std::ostream& operator<<(std::ostream& os, const ARGNode& node);
};

#endif // __ARG_NODE_HPP_
