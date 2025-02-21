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

#include "arg_edge.hpp"
#include "arg.hpp"
#include "arg_node.hpp"

#include <cassert>
#include <iostream>
#include <vector>

using std::cout;
using std::endl;
using std::ostream;
using std::vector;

ARGEdge::ARGEdge(arg_real_t _start, arg_real_t _end, ARGNode* _child, ARGNode* _parent)
    : start(_start), end(_end), child(_child), parent(_parent), mutations(nullptr) {
  assert(start < end);
  assert(start >= parent->start);
  assert(start >= child->start);
  assert(end <= parent->end);
  assert(end <= child->end);
}

ARGEdge::ARGEdge(ARGEdge const& edge)
    : start(edge.start), end(edge.end), child(edge.child), parent(edge.parent) {
  if (edge.mutations != nullptr) {
    mutations = std::make_unique<std::vector<Mutation*>>();
    for (Mutation* m : *(edge.mutations)) {
      mutations->push_back(m);
    }
  }
  else {
    mutations = nullptr;
  }
}

ARGEdge::~ARGEdge() {
#ifdef _DEBUG
  cout << "Deleting: " << *this << " (parent node may be garbage)" << endl;
#endif // _DEBUG
}

void ARGEdge::update_start(arg_real_t new_start) {
  if (new_start > start && mutations != nullptr) {
    assert(new_start < end);
    for (int i = mutations->size() - 1; i >= 0; --i) {
      if ((*mutations)[i]->position < new_start) {
        // https://stackoverflow.com/a/3487736/
        (*mutations)[i] = mutations->back();
        mutations->pop_back();
      }
    }
  }
  start = new_start;
}

void ARGEdge::update_end(arg_real_t new_end) {
  if (new_end < end && mutations != nullptr) {
    assert(new_end > start);
    for (int i = mutations->size() - 1; i >= 0; --i) {
      if ((*mutations)[i]->position >= new_end) {
        // https://stackoverflow.com/a/3487736/
        (*mutations)[i] = mutations->back();
        mutations->pop_back();
      }
    }
  }
  end = new_end;
}

void ARGEdge::add_mutations(vector<Mutation*> new_mutations) {
  // Can consider doing a site-check here, we don't want many mutations
  // for the same site.
  if (new_mutations.size() == 0) {
    return;
  }
  else if (mutations == nullptr) {
    mutations = std::make_unique<vector<Mutation*>>();
  }
  for (Mutation* mutation : new_mutations) {
    mutations->push_back(mutation);
  }
}

// This operation is linear in the number of mutations on the edge
bool ARGEdge::mutated_at_site(int site_id) {
  if (mutations != nullptr) {
    for (Mutation* mut : *mutations) {
      if (mut->site_id == site_id) {
        return true;
      }
    }
  }
  return false;
}

// Removes any/all mutations at the given site
// This operation is linear in the number of mutations on the edge
void ARGEdge::remove_mutations_at_site(int site_id) {
  if (mutations != nullptr) {
    for (int i = mutations->size() - 1; i >= 0; --i) {
      if ((*mutations)[i]->site_id == site_id) {
        // https://stackoverflow.com/a/3487736/
        (*mutations)[i] = mutations->back();
        mutations->pop_back();
      }
    }
  }
}

// Makes extra copies of the pointers
// This operation is linear in the number of mutations on the edge
std::vector<Mutation*> ARGEdge::mutations_in_range(arg_real_t range_start, arg_real_t range_end) {
  std::vector<Mutation*> range_mutations;
  if (mutations != nullptr) {
    for (Mutation* mut : *mutations) {
      if (range_start <= mut->position && mut->position < range_end) {
        range_mutations.push_back(mut);
      }
    }
  }
  return range_mutations;
}

ostream& operator<<(ostream& os, const ARGEdge& edge) {
  os << "Edge from node " << edge.child->ID << " to node " << edge.parent->ID;
  os << " with range [" << edge.start << ", " << edge.end << ")";
  return os;
}
