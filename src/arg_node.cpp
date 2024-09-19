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
#include "IntervalTree.h"
#include "arg_edge.hpp"
#include "utils.hpp"

#include <cassert>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

using std::cout;
using std::endl;
using std::map;
using std::ostream;
using std::string;
using std::vector;

ARGNode::ARGNode(int _ID, arg_real_t _height, arg_real_t _start, arg_real_t _end)
    : ID(_ID), height(_height), start(_start), end(_end) {
  assert(start < end);
  assert(height >= 0);
}

ARGNode::~ARGNode() {
#ifdef _DEBUG
  cout << "Deleting node " << std::to_string(ID) << endl;
#endif // _DEBUG
}

ARGEdge* ARGNode::add_parent(arg_real_t start, arg_real_t end, ARGNode& parent) {
  assert(parents.find(start) == parents.end());
  std::unique_ptr<ARGEdge> edge = std::make_unique<ARGEdge>(start, end, this, &parent);
  // Result is a pair <iterator, bool> with iterator pointing to the
  // newly inserted element and bool indicating success/failure of insertion
  auto result = parents.insert(std::make_pair(start, std::move(edge)));
  assert(result.second); // check that the insert succeeded
  return result.first->second.get();
}

void ARGNode::remove_parent(arg_real_t start) {
  parents.erase(start);
}

void ARGNode::update_parent_start(arg_real_t start_old, arg_real_t start_new) {
  assert(parents.find(start_new) == parents.end());
  std::unique_ptr<ARGEdge> edge = std::move(parents.at(start_old));
  edge->update_start(start_new);
  remove_parent(start_old);
  auto result = parents.insert(std::make_pair(start_new, std::move(edge)));
  assert(result.second); // check that the insert succeeded
}

void ARGNode::update_parent_end(arg_real_t start, arg_real_t end_new) {
  parents.at(start)->update_end(end_new);
}

vector<ARGEdge*> ARGNode::children_at(arg_real_t position) const {
  vector<ARGEdge*> result;
  if (children == nullptr) {
    return result;
  }
  vector<Interval<arg_real_t, ARGEdge*>> intervals;
  intervals = children->findOverlapping(position, position);
  for (auto interval : intervals) {
    if (interval.stop != position) {
      assert(interval.start <= position && position < interval.stop);
      result.push_back(interval.value);
    }
  }
  return result; // return by value
}

// pos_start is inclusive, pos_end is exclusive
vector<ARGEdge*> ARGNode::children_overlap(arg_real_t pos_start, arg_real_t pos_end) const {
  vector<ARGEdge*> result;
  if (children == nullptr) {
    return result;
  }
  vector<Interval<arg_real_t, ARGEdge*>> intervals;
  intervals = children->findOverlapping(pos_start, pos_end);
  for (auto interval : intervals) {
    if (interval.stop != pos_start && interval.start != pos_end) {
      assert(interval.start < pos_end && pos_start < interval.stop);
      result.push_back(interval.value);
    }
  }
  return result; // return by value
}

// Throws exception if position is not within parent range
// Returns nullptr if no parent is found
ARGEdge* ARGNode::parent_edge_at(arg_real_t position) const {
  if (position < start || position >= end) {
    throw std::logic_error(THROW_LINE("Position out of bounds."));
  }
  auto search = parents.upper_bound(position);
  if (search == parents.begin()) {
    return nullptr;
  }
  search = std::prev(search);
  ARGEdge* potential_parent = search->second.get();
  if (potential_parent->start <= position && position < potential_parent->end) {
    return potential_parent;
  }
  return nullptr;
}

ostream& operator<<(ostream& os, const ARGNode& node) {
  os << "Node " << node.ID << ": [" << node.start;
  os << ", " << node.end << "), height: " << node.height;
  os << ", parents: {";

  string subset = "";
  for (auto const& map_entry : node.parents) {
    subset += std::to_string((*map_entry.second).parent->ID) + ", ";
  }
  os << subset.substr(0, subset.size() - 2);

  os << "}";
  return os;
}
