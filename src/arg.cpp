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

#include "arg.hpp"
#include "IntervalTree.h"
#include "arg_edge.hpp"
#include "arg_node.hpp"
#include "descendant_list.hpp"
#include "utils.hpp"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <deque>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <queue>
#include <set>
#include <stack>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

using std::cout;
using std::deque;
using std::endl;
using std::map;
using std::ostream;
using std::queue;
using std::set;
using std::stack;
using std::string;
using std::tuple;
using std::unordered_map;
using std::vector;

// Special constructor to test memory footprint
// Builds a lot of UPGMA-like separate trees
ARG::ARG(arg_real_t _start, arg_real_t _end, int num_samples, int num_trees, int _reserved_samples)
    : start(_start), end(_end) {

  threaded_samples = num_samples;
  if (_reserved_samples == -1) {
    reserved_samples = num_samples;
  }
  else if (_reserved_samples >= num_samples) {
    reserved_samples = _reserved_samples;
  }
  else {
    throw std::logic_error(THROW_LINE("Must reserve more samples than passed."));
  }

  // this could be a good use case for C++20 std::lerp
  vector<arg_real_t> breakpoints;
  arg_real_t tree_width = (end - start) / num_trees;
  for (int j = 0; j < num_trees; ++j) {
    breakpoints.push_back(start + j * tree_width);
  }
  breakpoints.push_back(end);

  for (int i = 0; i < num_samples; ++i) {
    arg_nodes.insert(std::make_pair(i, std::make_unique<ARGNode>(i, 0, start, end)));
    leaf_ids.insert(i);
    sample_names.insert({i, ""});
  }

  int node_offset = reserved_samples;
  for (int j = 0; j < num_trees; ++j) {
    for (int i = 0; i < num_samples - 1; ++i) {
      arg_nodes.insert(std::make_pair(
          node_offset + i,
          std::make_unique<ARGNode>(node_offset + i, i + 1, breakpoints[j], breakpoints[j + 1])));
      ARGNode* parent = arg_nodes.at(node_offset + i).get();
      ARGNode* child1 = arg_nodes.at(i).get();

      int child2_index = node_offset + i;
      if (i == 0) {
        child2_index = num_samples - 1;
      }
      ARGNode* child2 = arg_nodes.at(child2_index).get();

      child1->add_parent(breakpoints[j], breakpoints[j + 1], *parent);
      child2->add_parent(breakpoints[j], breakpoints[j + 1], *parent);
      num_edges_cnt += 2;
    }
    node_offset += num_samples - 1;
  }

  if (_reserved_samples == -1) {
    next_general_id = arg_nodes.size();
  }
  else if (_reserved_samples >= num_samples) {
    next_general_id = arg_nodes.size() + _reserved_samples - num_samples;
  }
}

// Special constructor to quickly instantiate when using tskit data
// Precondition: the samples are at the start, so is_sample goes from all true
// to all false. TODO: check this precondition.
ARG::ARG(arg_real_t _start, arg_real_t _end, const vector<arg_real_t>& node_heights,
         const deque<bool>& is_sample, const vector<std::pair<int, int>>& edge_ids,
         const vector<std::pair<arg_real_t, arg_real_t>>& edge_ranges, int _reserved_samples)
    : start(_start), end(_end) {
  assert(node_heights.size() == is_sample.size());
  assert(edge_ids.size() == edge_ranges.size());

  for (size_t i = 0; i < node_heights.size(); ++i) {
    if (is_sample[i]) {
      ++threaded_samples;
    }
  }
  if (_reserved_samples == -1) {
    reserved_samples = threaded_samples;
    next_general_id = node_heights.size();
  }
  else if (_reserved_samples >= threaded_samples) {
    reserved_samples = _reserved_samples;
    next_general_id = node_heights.size() + _reserved_samples - threaded_samples;
  }
  else {
    throw std::logic_error(THROW_LINE("Must reserve more samples than passed."));
  }

  for (size_t i = 0; i < node_heights.size(); ++i) {
    // here using nodes extend for the whole ARG to read from tskit format
    int node_id = i;
    if (node_id >= threaded_samples) {
      node_id += reserved_samples - threaded_samples;
    }
    arg_nodes.insert(
        std::make_pair(node_id, std::make_unique<ARGNode>(node_id, node_heights[i], start, end)));
    if (is_sample[i]) {
      assert(node_heights[i] == 0);
      leaf_ids.insert(i);
      sample_names.insert({i, ""});
    }
  }

  for (size_t i = 0; i < edge_ids.size(); ++i) {
    // if we want to allow sample reservation, but are reading from tskit format, we should
    // bump up the IDs of nodes that are not samples
    int first_id = edge_ids[i].first;
    int second_id = edge_ids[i].second;
    if (first_id >= threaded_samples) {
      first_id += reserved_samples - threaded_samples;
    }
    if (second_id >= threaded_samples) {
      second_id += reserved_samples - threaded_samples;
    }
    ARGNode* child = arg_nodes.at(first_id).get();
    ARGNode* parent = arg_nodes.at(second_id).get();
    child->add_parent(edge_ranges[i].first, edge_ranges[i].second, *parent);
    ++num_edges_cnt;
  }
}

ARG::ARG(arg_real_t _start, arg_real_t _end, int _reserved_samples)
    : start(_start), end(_end), reserved_samples(_reserved_samples) {
  next_general_id = _reserved_samples;
}

ARG::ARG(const DeserializationParams& dp)
    : start{dp.start}, end{dp.end}, offset{dp.offset}, chromosome{dp.chromosome},
      threaded_samples{dp.threaded_samples}, reserved_samples{dp.reserved_samples} {

  if (reserved_samples == -1) {
    reserved_samples = threaded_samples;
    next_general_id = dp.num_nodes;
  }
  else if (reserved_samples >= threaded_samples) {
    next_general_id = dp.num_nodes + reserved_samples - threaded_samples;
  }
  else {
    throw std::logic_error(THROW_LINE("Must reserve more samples than passed."));
  }
}

void ARG::deserialize_add_nodes(const std::vector<double>& node_heights,
                                const std::vector<uint8_t>& is_sample,
                                const std::vector<std::array<double, 2>>& node_bounds) {

  if (node_heights.size() != is_sample.size()) {
    throw std::logic_error(THROW_LINE("Expected sane number of node heights and sample flags"));
  }
  if (!node_bounds.empty() && node_bounds.size() != node_heights.size()) {
    throw std::logic_error(THROW_LINE("Expected sane number of node heights and node bounds"));
  }

  std::size_t next_id = arg_nodes.size();

  for (std::size_t i = 0ul; i < node_heights.size(); ++i) {

    // Nodes may already exist from previous chunks, so we start at the next available ID
    std::size_t node_id = next_id + i;

    if (node_id >= threaded_samples) {
      node_id += reserved_samples - threaded_samples;
    }

    const arg_real_t node_start = node_bounds.empty() ? start : node_bounds[i].front();
    const arg_real_t node_end = node_bounds.empty() ? end : node_bounds[i].back();

    arg_nodes.insert(std::make_pair(
        node_id, std::make_unique<ARGNode>(node_id, node_heights[i], node_start, node_end)));

    if (is_sample[i]) {
      assert(node_heights[i] == 0.0);
      leaf_ids.insert(node_id);
      sample_names.insert({node_id, ""});
    }
  }
}

void ARG::deserialize_add_edges(const std::vector<std::array<int, 2>>& edge_ids,
                                const std::vector<std::array<double, 2>>& edge_ranges) {

  if (edge_ids.size() != edge_ranges.size() || edge_ids.empty()) {
    throw std::logic_error(THROW_LINE("Expected non-zero & same number of edge ids and ranges"));
  }
  if (edge_ids.front().size() != 2ul || edge_ranges.front().size() != 2ul) {
    throw std::logic_error(THROW_LINE("Expected edge ids and ranges to have shape Nx2"));
  }

  for (size_t i = 0ul; i < edge_ids.size(); ++i) {

    int first_id = edge_ids[i].front();
    int second_id = edge_ids[i].back();
    if (first_id >= threaded_samples) {
      first_id += reserved_samples - threaded_samples;
    }
    if (second_id >= threaded_samples) {
      second_id += reserved_samples - threaded_samples;
    }
    ARGNode* child = arg_nodes.at(first_id).get();
    ARGNode* parent = arg_nodes.at(second_id).get();
    child->add_parent(edge_ranges[i].front(), edge_ranges[i].back(), *parent);
    ++num_edges_cnt;
  }
}

void ARG::deserialize_add_mutations(const std::vector<arg_real_t>& positions,
                                    const std::vector<arg_real_t>& heights,
                                    const std::vector<int>& site_ids,
                                    const std::vector<std::array<int, 2>>& edge_ids) {

  if (positions.size() != heights.size() || positions.size() != site_ids.size() ||
      positions.size() != edge_ids.size()) {
    throw std::logic_error(
        THROW_LINE("Expected positions, heights, site_ids and edge_ids to have the same length"));
  }

  for (size_t i = 0ul; i < positions.size(); ++i) {

    const arg_real_t this_pos = positions.at(i);
    const arg_real_t this_height = heights.at(i);

    int first_id = edge_ids[i].front();
    int second_id = edge_ids[i].back();
    if (first_id >= threaded_samples) {
      first_id += reserved_samples - threaded_samples;
    }
    if (second_id >= threaded_samples) {
      second_id += reserved_samples - threaded_samples;
    }

    ARGNode* child = arg_nodes.at(first_id).get();
    ARGNode* parent = arg_nodes.at(second_id).get();

    auto edge_it = child->parents.upper_bound(this_pos);
    if (edge_it == child->parents.begin()) {
      throw std::logic_error(THROW_LINE("Edge nodes and mutation position do not correspond."));
    }

    ARGEdge* edge_ptr = std::prev(edge_it)->second.get();

    const bool valid_id = edge_ptr->parent->ID == parent->ID;
    const bool valid_pos = edge_ptr->start <= this_pos && edge_ptr->end > this_pos;
    const bool valid_height = child->height <= this_height && parent->height > this_height;

    if (!valid_id || !valid_pos || !valid_height) {
      throw std::logic_error(THROW_LINE("Could not find correct edge for serialized mutation."));
    }

    add_mutation(edge_ptr, positions.at(i), heights.at(i), site_ids.at(i), false);
  }
  update_site_data_structures();
}

void ARG::add_mutation(ARGEdge* edge, const arg_real_t position, const arg_real_t height,
    const int site_id, const bool update_data_structures)
{
  if (mutations.empty() || position >= mutations.back()->position) {
    mutations.emplace_back(std::make_unique<Mutation>(edge, position, height, site_id));
  } else {
    auto mut_it = std::lower_bound(
        mutations.begin(), mutations.end(), std::make_unique<Mutation>(nullptr, position),
        [](const std::unique_ptr<Mutation>& lhs, const std::unique_ptr<Mutation>& rhs) {
          return *lhs < *rhs;
        });
    mutations.insert(mut_it, std::make_unique<Mutation>(edge, position, height, site_id));
  }

  mutation_sites_up_to_date = false;
  site_positions_up_to_date = false;

  if (update_data_structures) {
    update_site_data_structures();
  }
}

void ARG::update_site_data_structures() const
{
  update_mutation_sites();
  update_site_positions();
}


void ARG::reserve_n_mutations(const std::size_t num_mutations) {
  mutations.reserve(num_mutations);
}

const std::vector<std::unique_ptr<Mutation>>& ARG::get_mutations() const {
  return mutations;
}

void ARG::clear_mutations() {
  mutations.clear();
}

std::size_t ARG::get_idx_of_first_mutation_left_of(const arg_real_t physical_pos,
                                                   const bool include_equal,
                                                   const bool warn_out_of_range) const {
  if (mutations.empty()) {
    throw std::logic_error(THROW_LINE("There are no mutations."));
  }

  decltype(mutations)::const_iterator mut_it;

  if (include_equal) {
    mut_it = std::lower_bound(mutations.begin(), mutations.end(),
                              std::make_unique<Mutation>(nullptr, physical_pos),
                              [](const std::unique_ptr<Mutation>& lhs,
                                 const std::unique_ptr<Mutation>& rhs) { return *lhs <= *rhs; });

    if (mut_it != mutations.begin()) {
      --mut_it;
    }
    if (warn_out_of_range && (*mut_it)->position > physical_pos) {
      std::cout << "Warning: no mutations with position <= " << physical_pos << '\n';
    }
  }
  else {
    mut_it = std::lower_bound(mutations.cbegin(), mutations.cend(),
                              std::make_unique<Mutation>(nullptr, physical_pos),
                              [](const std::unique_ptr<Mutation>& lhs,
                                 const std::unique_ptr<Mutation>& rhs) { return *lhs < *rhs; });

    if (mut_it != mutations.begin()) {
      --mut_it;
    }
    if (warn_out_of_range && (*mut_it)->position >= physical_pos) {
      std::cout << "Warning: no mutations with position < " << physical_pos << '\n';
    }
  }

  return static_cast<std::size_t>(std::distance(mutations.begin(), mut_it));
}

std::size_t ARG::get_idx_of_first_mutation_right_of(const arg_real_t physical_pos,
                                                    const bool include_equal,
                                                    const bool warn_out_of_range) const {
  if (mutations.empty()) {
    throw std::logic_error(THROW_LINE("There are no mutations."));
  }

  decltype(mutations)::const_iterator mut_it;

  if (include_equal) {
    mut_it = std::lower_bound(mutations.cbegin(), mutations.cend(),
                              std::make_unique<Mutation>(nullptr, physical_pos),
                              [](const std::unique_ptr<Mutation>& lhs,
                                 const std::unique_ptr<Mutation>& rhs) { return *lhs < *rhs; });

    if (mut_it == mutations.end()) {
      --mut_it;
    }
    if (warn_out_of_range && (*mut_it)->position < physical_pos) {
      std::cout << "Warning: no mutations with position >= " << physical_pos << '\n';
    }
  }
  else {
    mut_it = std::lower_bound(mutations.begin(), mutations.end(),
                              std::make_unique<Mutation>(nullptr, physical_pos),
                              [](const std::unique_ptr<Mutation>& lhs,
                                 const std::unique_ptr<Mutation>& rhs) { return *lhs <= *rhs; });

    if (mut_it == mutations.end()) {
      --mut_it;
    }
    if (warn_out_of_range && (*mut_it)->position <= physical_pos) {
      std::cout << "Warning: no mutations with position > " << physical_pos << '\n';
    }
  }

  return static_cast<std::size_t>(std::distance(mutations.begin(), mut_it));
}

std::size_t ARG::get_idx_of_mutation_closest_to(arg_real_t physical_pos) const {
  if (mutations.empty()) {
    throw std::logic_error(THROW_LINE("There are no mutations."));
  }

  const std::size_t candidate_idx_1 = get_idx_of_first_mutation_left_of(physical_pos, false, false);

  if (candidate_idx_1 == mutations.size() - 1ul) {
    return candidate_idx_1;
  }

  const std::size_t candidate_idx_2 = candidate_idx_1 + 1ul;

  const arg_real_t candidate_dist_1 =
      std::fabs(mutations.at(candidate_idx_1)->position - physical_pos);
  const arg_real_t candidate_dist_2 =
      std::fabs(mutations.at(candidate_idx_2)->position - physical_pos);

  return candidate_dist_1 < candidate_dist_2 ? candidate_idx_1 : candidate_idx_2;
}

void ARG::set_offset(int _offset) {
  if (start != 0) {
    throw std::logic_error(THROW_LINE("To use a global ARG offset, start should be 0."));
  }
  if (_offset < 0) {
    throw std::logic_error(THROW_LINE("Offset must be nonnegative."));
  }
  offset = _offset;
}

void ARG::set_chromosome(int _chromosome) {
  if (_chromosome < 1) {
    throw std::logic_error(THROW_LINE("Chromosome must be a positive integer."));
  }
  chromosome = _chromosome;
}

bool ARG::is_leaf(int node_id) const {
  return (leaf_ids.find(node_id) != leaf_ids.end());
}

void ARG::populate_children_and_roots() {
  populate_children();
  populate_roots();
  // use vector always by default, but can subsequently override the threshold
  DescendantList::set_threshold(leaf_ids.size());
  // DescendantList::print_threshold();
}

void ARG::populate_children() {
  unordered_map<int, vector<Interval<arg_real_t, ARGEdge*>>> id_to_intervals;
  for (auto const& map_entry : arg_nodes) {
    ARGNode* node = map_entry.second.get();
    for (auto const& edge_entry : node->parents) {
      ARGEdge* edge = edge_entry.second.get();
      int parent_id = edge->parent->ID;
      if (id_to_intervals.find(parent_id) == id_to_intervals.end()) {
        id_to_intervals.insert({parent_id, vector<Interval<arg_real_t, ARGEdge*>>()});
      }
      id_to_intervals.at(parent_id).push_back(
          Interval<arg_real_t, ARGEdge*>(edge->start, edge->end, edge));
    }
  }

  for (auto const& map_entry : arg_nodes) {
    ARGNode* node = map_entry.second.get();
    if (id_to_intervals.find(node->ID) != id_to_intervals.end()) {
      vector<Interval<arg_real_t, ARGEdge*>> node_intervals = id_to_intervals.at(node->ID);
      // rewrites anything present
      node->children =
          std::make_unique<IntervalTree<arg_real_t, ARGEdge*>>(std::move(node_intervals));
    }
    // otherwise, we have a leaf, which should have no children, or node->children = nullptr
    // we could do an empty IntervalTree as well but it seems wasteful
  }
}

void ARG::populate_roots() {
  roots.clear(); // frees memory of any Roots, because those are in unique_ptr
  int id;
  // try to provide some determinism by using 0, which should always be leaf
  if (leaf_ids.find(0) != leaf_ids.end()) {
    id = 0;
  }
  else {
    id = *leaf_ids.begin();
    cout << "Warning: roots are calculated from leaf " << id << " instead of 0" << endl;
  }

  ARGNode* leaf_node = arg_nodes.at(id).get();
  stack<tuple<ARGNode*, arg_real_t, arg_real_t>> entries;
  entries.push(std::make_tuple(leaf_node, leaf_node->start, leaf_node->end));

  while (!entries.empty()) {
    ARGNode* node;
    arg_real_t entry_start, entry_end;
    std::tie(node, entry_start, entry_end) = entries.top();
    entries.pop();

    // COMPUTE AFFECTED PARENTS AND GAPS
    vector<ARGEdge*> affected_parents;
    vector<std::pair<arg_real_t, arg_real_t>> gaps;

    arg_real_t curr_position = entry_start;
    auto affected_first = node->parents.upper_bound(entry_start);
    // at this point, affected_first->start > entry_start or affected_first = node->parents.end()
    if (affected_first != node->parents.begin()) {
      affected_first = std::prev(affected_first);
    }
    auto affected_after_last = node->parents.lower_bound(entry_end);
    // at this point, affected_after_last->start >= entry_end or affected_after_last =
    // node->parents.end()

    for (auto map_it = affected_first; map_it != affected_after_last; ++map_it) {
      arg_real_t parent_start = map_it->second->start;
      arg_real_t parent_end = map_it->second->end;
      if (parent_end <= curr_position) {
        // we tried to get the first parent, but it doesn't overlap
        continue;
      }
      if (parent_start > curr_position) {
        gaps.push_back(std::make_pair(curr_position, parent_start));
      }
      affected_parents.push_back(map_it->second.get());
      curr_position = parent_end;
    }
    if (curr_position < entry_end) {
      gaps.push_back(std::make_pair(curr_position, entry_end));
    }

#ifdef _DEBUG
    cout << endl;
    cout << "Processing entry: [" << entry_start << ", " << entry_end << ") with node " << node->ID
         << endl;
    cout << "Affected parents:" << endl;
    for (auto affected : affected_parents) {
      cout << *affected << endl;
    }
    cout << "Gaps:" << endl;
    for (auto gap : gaps) {
      cout << "[" << gap.first << ", " << gap.second << ")" << endl;
    }
#endif // _DEBUG

    // HANDLE GAPS
    for (std::pair<arg_real_t, arg_real_t> gap : gaps) {
      roots.insert(std::make_pair(gap.first, std::make_unique<Root>(node, gap.first, gap.second)));
    }

    // HANDLE AFFECTED PARENTS
    for (ARGEdge* edge : affected_parents) {
      arg_real_t left = std::max(entry_start, edge->start);
      arg_real_t right = std::min(entry_end, edge->end);
      entries.push(std::make_tuple(edge->parent, left, right));
    }
  }
}

void ARG::clear_mutations_from_edges() {
  for (auto const& map_entry : arg_nodes) {
    ARGNode* node = map_entry.second.get();
    for (auto const& edge_entry : node->parents) {
      ARGEdge* edge = edge_entry.second.get();
      edge->mutations = nullptr;
    }
  }
}

// Iterate over mutations and add to the vector of mutations
// belonging to the edge
void ARG::populate_mutations_on_edges() {
  // First clear all mutations from all edges. Slightly expensive traversal.
  clear_mutations_from_edges();

  // Then visit all mutations and add to edge vector
  for (const std::unique_ptr<Mutation>& mutation : mutations) {
    ARGEdge* edge = mutation->edge;
    if (edge->mutations == nullptr) {
      edge->mutations = std::make_unique<std::vector<Mutation*>>();
    }
    edge->mutations->push_back(mutation.get());
  }
}

vector<arg_real_t> ARG::root_starts() const {
  if (roots.empty()) {
    throw std::logic_error(THROW_LINE("Call populate_children_and_roots() first."));
  }
  vector<arg_real_t> starts;
  for (auto const& map_entry : roots) {
    starts.push_back(map_entry.first);
  }
  return starts; // return by value
}

const Root* ARG::root_at(arg_real_t position) const {
  auto it = roots.upper_bound(position);
  if (it == roots.begin()) {
    throw std::out_of_range(THROW_LINE("No root exists here."));
  }
  it = std::prev(it);
  Root* root = it->second.get();
  assert(root->start <= position && position < root->end);
  return root;
}

// Note that we don't actually use the ARG in this function
// end is the extent to which this MRCA is guaranteed to be valid
tuple<const ARGNode*, arg_real_t>
ARG::mrca_nodes_with_end(const ARGNode& node1, const ARGNode& node2, arg_real_t position) const {
  if (position < node1.start || position >= node1.end) {
    throw std::logic_error(THROW_LINE("Invalid first node for this position."));
  }
  if (position < node2.start || position >= node2.end) {
    throw std::logic_error(THROW_LINE("Invalid second node for this position."));
  }

  const ARGNode* lower = &node1;
  const ARGNode* upper = &node2;
  if (node1.height > node2.height) {
    lower = &node2;
    upper = &node1;
  }

  arg_real_t end = std::min(lower->end, upper->end);
  while (lower != upper) {
    const ARGEdge* edge = lower->parent_edge_at(position);
    if (edge == nullptr) {
      throw std::logic_error(THROW_LINE("Unable to find an MRCA."));
    }
    end = std::min(end, edge->end);

    // Note: edge->parent does not return const
    const ARGNode* tmp = edge->parent;
    if (tmp->height <= upper->height) {
      lower = tmp;
    }
    else {
      lower = upper;
      upper = tmp;
    }
  }
  assert(end > position);
  return std::make_tuple(lower, end);
}

// Find lowest branch that carries mutation for chromosome ID
// uses a tolerance of 1e-3, which is very reasonable for practical applications
ARGEdge* ARG::lowest_mutated_edge(int haploid_id, arg_real_t position) {
  const ARGNode* node = arg_nodes.at(haploid_id).get();
  ARGEdge* edge = node->parent_edge_at(position);
  float tol = 1e-3;
  while (edge != nullptr) {
    if (edge->mutations != nullptr) {
      for (Mutation* mut : *(edge->mutations)) {
        if (mut->position - tol <= position && position < mut->position + tol) {
          return edge;
        }
      }
    }
    edge = edge->parent->parent_edge_at(position);
  }
  return edge;
}

// Find lowest branch that carries mutation for chromosome ID searching by site ID
ARGEdge* ARG::lowest_mutated_edge_by_site(int haploid_id, int site_id) {
  update_site_positions();
  return lowest_mutated_edge(haploid_id, site_positions[site_id]);
}

const ARGNode* ARG::mrca(int ID1, int ID2, arg_real_t position) const {
  const ARGNode* node;
  arg_real_t end;
  std::tie(node, end) = mrca_nodes_with_end(*arg_nodes.at(ID1), *arg_nodes.at(ID2), position);
  return node;
}

tuple<const ARGNode*, arg_real_t> ARG::mrca_with_end(int ID1, int ID2, arg_real_t position) const {
  return mrca_nodes_with_end(*arg_nodes.at(ID1), *arg_nodes.at(ID2), position);
}

set<arg_real_t> ARG::get_breakpoints() const {
  set<arg_real_t> breakpoints;
  for (auto const& map_entry : arg_nodes) {
    breakpoints.insert(map_entry.second->start);
    for (auto const& edge_entry : map_entry.second->parents) {
      breakpoints.insert(edge_entry.second->start);
    }
  }
  return breakpoints; // return by value
}

ostream& operator<<(ostream& os, const ARG& arg) {
  // Sort the nodes by ID before printing
  vector<int> keys;
  for (auto const& map_entry : arg.arg_nodes) {
    keys.push_back(map_entry.first);
  }
  std::sort(keys.begin(), keys.end());
  for (int key : keys) {
    const ARGNode& node = *(arg.arg_nodes.at(key));
    os << node << endl;
  }
  return os;
}

// Return an iterator pointing to the first mutation at or after this position, or
// the vector's end() if no such mutation exists. We allow edge to be nullptr so
// we can query mutations by position using a dummy mutation.
vector<std::unique_ptr<Mutation>>::const_iterator ARG::next_mutation(arg_real_t pos) const {
  // lower_bound uses binary search to find first mutation after position in O(logM) time
  auto low = std::lower_bound(mutations.cbegin(), mutations.cend(),
                              std::make_unique<Mutation>(nullptr, pos, 0.0),
                              [](const std::unique_ptr<Mutation>& a,
                                 const std::unique_ptr<Mutation>& b) { return *a < *b; });
  return low;
}

void ARG::check_basic(bool stringent) const {
  // Check that samples have height 0 and other nodes have height > 0
  check_node_heights();
  // Check that samples span the whole chromosome and other nodes have good spans
  check_node_spans();
  // Check that the edges in parents all look good
  check_edges();
  // Check parent relationships, make sure the gaps (areas with no parents)
  // form local roots which partition the ARG length
  check_single_parent_except_root_gaps(stringent);
  // Check num_nodes and num_edges are correct
  check_stats();
}

void ARG::check_node_heights() const {
  for (auto const& map_entry : arg_nodes) {
    ARGNode* node = map_entry.second.get();
    bool is_leaf = leaf_ids.find(node->ID) != leaf_ids.end();
    string node_id = std::to_string(node->ID);

    if (is_leaf && node->height != 0) {
      throw std::logic_error(THROW_LINE("Incorrect height for leaf node " + node_id + "."));
    }
    else if (!is_leaf && node->height <= 0) {
      throw std::logic_error(THROW_LINE("Incorrect height for non-leaf node " + node_id + "."));
    }
  }
}

void ARG::ARG::check_node_spans() const {
  if (start >= end) {
    throw std::logic_error(THROW_LINE("Bad ARG span."));
  }
  for (auto const& map_entry : arg_nodes) {
    ARGNode* node = map_entry.second.get();
    bool is_leaf = leaf_ids.find(node->ID) != leaf_ids.end();
    string node_id = std::to_string(node->ID);

    if (is_leaf && (node->start != start || node->end != end)) {
      throw std::logic_error(THROW_LINE("Bad span for leaf node " + node_id + "."));
    }
    else if (node->start >= node->end || node->start < start || node->end > end) {
      throw std::logic_error(THROW_LINE("Bad span for node " + node_id + "."));
    }
  }
}

void ARG::check_edges() const {
  for (auto const& map_entry : arg_nodes) {
    ARGNode* node = map_entry.second.get();
    for (auto const& edge_entry : node->parents) {
      ARGEdge* edge = edge_entry.second.get();
      ARGNode* parent = edge->parent;
      if (edge_entry.first != edge->start) {
        throw std::logic_error(THROW_LINE("Inconsistent edge key."));
      }
      if (edge->child != node) {
        throw std::logic_error(THROW_LINE("Inconsistent child pointer."));
      }
      if (edge->parent->height <= node->height) {
        throw std::logic_error(THROW_LINE("Parent should be older than child."));
      }
      if ((edge->start >= edge->end) || (edge->start < node->start) || (edge->end > node->end) ||
          (edge->start < parent->start) || (edge->end > parent->end)) {
        throw std::logic_error(THROW_LINE("Bad edge span."));
      }
    }
  }
}

// Return the number of nodes in the ARG
int ARG::num_nodes() const {
  return arg_nodes.size();
}

// Return the number of edges in the ARG
int ARG::num_edges() const {
  return num_edges_cnt;
}

int ARG::num_mutations() const {
  return mutations.size();
}

std::size_t ARG::get_num_sites() const {
  return site_positions.size();
}

// Re-compute stats from the ARG
// currently just num_nodes and num_edges
void ARG::check_stats() const {
  int chk_num_nodes = 0;
  int chk_num_edges = 0;
  for (auto const& map_entry : arg_nodes) {
    ARGNode* node = map_entry.second.get();
    // increment the number of edges
    chk_num_edges += node->parents.size();
    ++chk_num_nodes;
  }
  if ((chk_num_nodes != (int) arg_nodes.size())) {
    throw std::logic_error(THROW_LINE("Inconsistent number of nodes."));
  }
  if (chk_num_edges != num_edges_cnt) {
    throw std::logic_error(THROW_LINE("Inconsistent number of edges."));
  }
}

// If not stringent, it's okay to have multiple root gaps for a stretch
void ARG::check_single_parent_except_root_gaps(bool stringent) const {
  vector<std::pair<arg_real_t, arg_real_t>> gaps; // to be filled in progressively

  for (auto const& map_entry : arg_nodes) {
    ARGNode* node = map_entry.second.get();
    arg_real_t position = node->start;
    // progress through edges in sorted order
    for (auto const& edge_entry : node->parents) {
      ARGEdge* edge = edge_entry.second.get();
      if (position > edge->start) {
        throw std::logic_error(THROW_LINE("Multiple parents for a stretch."));
      }
      if (position < edge->start) {
        gaps.push_back(std::make_pair(position, edge->start));
      }
      position = edge->end;
    }
    if (position < node->end) {
      gaps.push_back(std::make_pair(position, node->end));
    }
  }

  // sort and check gaps against ARG start and end
  sort(gaps.begin(), gaps.end());
  if (gaps.empty()) {
    throw std::logic_error(THROW_LINE("Root gaps don't span the ARG."));
  }
  if (gaps[0].first != start) {
    throw std::logic_error(THROW_LINE("Root gaps don't span the ARG."));
  }
  arg_real_t pos = gaps[0].second;
  for (size_t i = 1; i < gaps.size(); ++i) {
    if (pos < gaps[i].first) {
      throw std::logic_error(THROW_LINE("Root gaps don't span the ARG."));
    }
    else if (pos > gaps[i].first && stringent) {
      throw std::logic_error(THROW_LINE("Multiple root gaps for a stretch."));
    }
    pos = std::max(pos, gaps[i].second);
  }
  if (pos != end) {
    throw std::logic_error(THROW_LINE("Root gaps don't span the ARG."));
  }
}

// The roots should partition the ARG, and the actual nodes given by the roots
// should have gaps precisely where they are labeled as roots
void ARG::check_roots() const {
  arg_real_t position = start;
  // goes through roots in sorted order by start location
  for (auto const& map_entry : roots) {
    Root* root = map_entry.second.get();
    if (map_entry.first != root->start) {
      throw std::logic_error(THROW_LINE("Inconsistent root key."));
    }
    if (position != root->start) {
      throw std::logic_error(THROW_LINE("Roots don't span the ARG."));
    }

    // Check that the node of this root has a parent gap in the root's span
    ARGNode* node = root->node;
    auto before_start = node->parents.upper_bound(root->start);
    if (before_start != node->parents.begin()) {
      before_start = std::prev(before_start);
      if (before_start->second->end > root->start) {
        throw std::logic_error(THROW_LINE("Expected a gap in parents."));
      }
    }
    auto after_start = node->parents.lower_bound(root->start);
    if (after_start != node->parents.end() && after_start->first < root->end) {
      throw std::logic_error(THROW_LINE("Expected a gap in parents."));
    }

    position = root->end;
  }

  if (position != end) {
    throw std::logic_error(THROW_LINE("Roots don't span the ARG."));
  }
}

void ARG::check_children(bool stringent) const {
  // Check child relationships, samples should have no children and other
  // nodes should have 2+ children
  check_number_of_children(stringent);
  // Bijection between child and parent data over the whole ARG
  check_parents_have_children();
  check_children_have_parents();
}

// For non-stringent check, we allow non-leaf nodes to have no children in
// ranges, but they still must have at least two children overall
void ARG::check_number_of_children(bool stringent) const {
  for (auto const& map_entry : arg_nodes) {
    const ARGNode* node = map_entry.second.get();
    if (is_leaf(node->ID)) {
      // as an extra check, we search within the ARG start and end
      vector<ARGEdge*> child_edges = node->children_overlap(start, end);
      if (!child_edges.empty()) {
        throw std::logic_error(THROW_LINE("Leaves should not have children."));
      }
    }
    else {
      vector<ARGEdge*> child_edges = node->children_overlap(start, node->start);
      if (!child_edges.empty()) {
        throw std::logic_error(THROW_LINE("Children outside of node range."));
      }
      child_edges = node->children_overlap(node->end, end);
      if (!child_edges.empty()) {
        throw std::logic_error(THROW_LINE("Children outside of node range."));
      }
      if (stringent) {
        arg_real_t position = node->start;
        while (position != node->end) {
          child_edges = node->children_at(position);
          if (child_edges.size() < 2) {
            throw std::logic_error(THROW_LINE("Non-leaf nodes need 2+ children in all ranges."));
          }
          arg_real_t next_pos = node->end;
          for (const ARGEdge* edge : child_edges) {
            next_pos = std::min(next_pos, edge->end);
          }
          position = next_pos;
        }
      }
      else {
        // weaker check, just check for 2+ children
        child_edges = node->children_overlap(node->start, node->end);
        if (child_edges.size() < 2) {
          throw std::logic_error(THROW_LINE("Non-leaf nodes need 2+ children."));
        }
      }
    }
  }
}

void ARG::check_parents_have_children() const {
  for (auto const& map_entry : arg_nodes) {
    const ARGNode* node = map_entry.second.get();
    for (auto const& edge_entry : node->parents) {
      const ARGEdge* edge = edge_entry.second.get();
      const ARGNode* parent = edge->parent;
      vector<ARGEdge*> children_of_parent = parent->children_at(edge->start);
      bool match = false;
      for (const ARGEdge* reverse_edge : children_of_parent) {
        if (reverse_edge == edge) {
          match = true;
        }
      }
      if (!match) {
        throw std::logic_error(THROW_LINE("Unable to find matching reverse edge."));
      }
    }
  }
}

void ARG::check_children_have_parents() const {
  for (auto const& map_entry : arg_nodes) {
    const ARGNode* node = map_entry.second.get();
    // as an extra check, we search within the ARG start and end
    vector<ARGEdge*> child_edges = node->children_overlap(start, end);
    for (const ARGEdge* edge : child_edges) {
      const ARGNode* child = edge->child;
      const ARGEdge* reverse_edge = child->parent_edge_at(edge->start);
      if (reverse_edge != edge) {
        throw std::logic_error(THROW_LINE("Unable to find matching reverse edge."));
      }
    }
  }
}

void ARG::check_mutations_sorted() const {
  // Will also check that mutations exist to begin with
  if (mutations.empty()) {
    throw std::logic_error(
        THROW_LINE("Need pre-existing mutations to run! Try generate_mutations_and_keep."));
  }
  if (!std::is_sorted(mutations.begin(), mutations.end())) {
    throw std::logic_error(THROW_LINE("Mutations must be sorted by position. They are sorted by "
                                      "default, so some funny business must have occurred."));
  }
}

int ARG::add_sample(string sample_name) {
  if (next_to_thread != -1) {
    throw std::logic_error(THROW_LINE("First thread the last added sample."));
  }
#ifdef _DEBUG
  cout << "Adding sample " << sample_name << endl;
#endif // _DEBUG
  int leaf_id;
  if (threaded_samples < reserved_samples) {
    leaf_id = threaded_samples;
  }
  else {
    leaf_id = next_general_id;
    ++next_general_id;
  }
  assert(arg_nodes.find(leaf_id) == arg_nodes.end());

  // create a new node for this sample
  arg_nodes.insert(std::make_pair(leaf_id, std::make_unique<ARGNode>(leaf_id, 0, start, end)));
  sample_names.insert({leaf_id, sample_name});
  leaf_ids.insert(leaf_id);
  ++threaded_samples;
  if (threaded_samples >= 2) {
    next_to_thread = leaf_id;
  }
  return leaf_id;
}

// @pre: section_starts, sample_ids, and heights are vectors of the same size
// @pre: section_starts are in increasing order and partition [start, end)
// @pre: sample_ids are leaf IDs in the ARG
// @pre: heights are >= 0
void ARG::thread_sample(vector<arg_real_t> section_starts, vector<int> sample_ids,
                        vector<arg_real_t> heights) {
  if (next_to_thread == -1) {
    throw std::logic_error(THROW_LINE("No samples to be threaded."));
  }
  size_t num_sections = section_starts.size();
  if (num_sections <= 0 || num_sections != sample_ids.size() || num_sections != heights.size()) {
    throw std::invalid_argument(THROW_LINE("Threading vectors must be of the same nonzero size"));
  }
  for (auto id : sample_ids) {
    if (leaf_ids.find(id) == leaf_ids.end() || id == next_to_thread) {
      throw std::invalid_argument(THROW_LINE("Threading IDs must specify samples in ARG"));
    }
  }
  // check sorted order
  if (section_starts[0] != start || section_starts[num_sections - 1] >= end) {
    throw std::invalid_argument(
        THROW_LINE("Section starts must partition the region [start, end)"));
  }
  for (size_t i = 0; i < num_sections - 1; ++i) {
    if (section_starts[i] >= section_starts[i + 1]) {
      throw std::invalid_argument(THROW_LINE("Section starts must be in increasing order"));
    }
  }
  // check heights are positive
  for (auto height : heights) {
    if (height <= 0) {
      throw std::invalid_argument(THROW_LINE("Threading heights (times) must be positive"));
    }
  }

  ARGNode* leaf_to_thread = arg_nodes.at(next_to_thread).get();

  // this used to be a queue, but I think stack is more intuitive
  stack<AncestorEntry> entries;
  arg_real_t section_start, section_end;
  for (size_t i = 0; i < num_sections; ++i) {
    section_start = section_starts[i];
    if (i == num_sections - 1) {
      section_end = end;
    }
    else {
      section_end = section_starts[i + 1];
    }
    ARGNode* node = arg_nodes.at(sample_ids[i]).get();
    std::unique_ptr<ARGNode> ancestor =
        std::make_unique<ARGNode>(next_general_id, heights[i], section_start, section_end);
    // add ancestor as a direct parent of leaf_to_thread in this range
    leaf_to_thread->add_parent(section_start, section_end, *ancestor);
    ++num_edges_cnt;
    // add ancestor as an ancestor of node in this range by adding an AncestorEntry
    entries.push(AncestorEntry(node, ancestor.get(), section_start, section_end));
    arg_nodes.insert(std::make_pair(next_general_id, std::move(ancestor)));
    ++next_general_id;
  }

  while (!entries.empty()) {
    process_ancestor_entry(entries); // this does all the work
  }

  next_to_thread = -1;
}

void ARG::process_ancestor_entry(stack<AncestorEntry>& entries) {
  /* Process a single AncestorEntry, popping it off and pushing on more elements
   * if needed.
   *
   * Consider the following example:
   *                                    ----P5---
   *    -P1---            -P3--
   *      -------------ancestor-------------
   *          --P2--
   *                           --P4--
   *
   * -----s-----------node------------------e------------
   *
   * We need to add ancestor as an ancestor of node, taking care of the parents
   * that may already exist in this region. The start and end points are denoted
   * by `s` and `e` respectively. (Remember that start is inclusive while end
   * is exclusive.)
   *
   * We first get all the relevant parents of the node in this region, as well
   * as the gaps, which are the regions with no parent (i.e. where node is
   * currently a root.) In the diagram, we have 5 "affected parents", P1 to P5,
   * and we have two gaps, P2 to P3 and P4 to P5. To quickly find the parents
   * and deduce the gaps, we use the lower_bound and upper_bound operations on
   * our map of parents for this node, with the start and end points as
   * querying keys.
   *
   * Dealing with the gaps is easy: we simply add ancestor as a parent of node
   * in these regions. For the other parents, if the parent height is less than
   * or equal to the ancestor height (P2 or P4), we push a new AncestorEntry
   * that says ancestor should be an ancestor of the parent in this region. This
   * is case A.
   *
   * If the parent height is greater than the ancestor height (P1, P3, P5), then
   * more work is necessary (case B). We split this possibility into 4 subcases
   * and handle them individually. In B1, the parent relationship is contained
   * within the ancestor range, as for P3. In B2, the parent relationship
   * straddles the left boundary of the ancestor range, as for P1. In B3, the
   * parent relationship straddles the right boundary of the ancestor range, as
   * for P5. In B4 (not pictured), the parent relationship would straddle both
   * the left and right boundaries, encompassing the entire ancestor range. Each
   * of these subcases is dealt with by modifying the parents of node and
   * pushing new AncestorEntry's as needed.
   */
  // AncestorEntry ancestor_entry = entries.front();
  AncestorEntry ancestor_entry = entries.top();
  entries.pop();
  ARGNode* node = ancestor_entry.node;
  ARGNode* ancestor = ancestor_entry.ancestor;
  arg_real_t entry_start = ancestor_entry.start;
  arg_real_t entry_end = ancestor_entry.end;

  /*
   * COMPUTE AFFECTED PARENTS AND GAPS
   */
  vector<ARGEdge*> affected_parents;
  vector<std::pair<arg_real_t, arg_real_t>> gaps;
  arg_real_t curr_position = entry_start;
  auto affected_first = node->parents.upper_bound(entry_start);
  // at this point, affected_first->start > entry_start or affected_first = node->parents.end()
  if (affected_first != node->parents.begin()) {
    affected_first = std::prev(affected_first);
  }
  auto affected_after_last = node->parents.lower_bound(entry_end);
  // at this point, affected_after_last->start >= entry_end or affected_after_last =
  // node->parents.end()

  for (auto map_it = affected_first; map_it != affected_after_last; ++map_it) {
    arg_real_t parent_start = map_it->second->start;
    arg_real_t parent_end = map_it->second->end;
    if (parent_end <= curr_position) {
      // we tried to get the first parent, but it doesn't overlap
      continue;
    }
    if (parent_start > curr_position) {
      gaps.push_back(std::make_pair(curr_position, parent_start));
    }
    affected_parents.push_back(map_it->second.get());
    curr_position = parent_end;
  }
  if (curr_position < entry_end) {
    gaps.push_back(std::make_pair(curr_position, entry_end));
  }

#ifdef _DEBUG
  cout << endl;
  cout << "Processing entry: [" << entry_start << ", " << entry_end << ") from " << node->ID
       << " to " << ancestor->ID << endl;
  cout << "Affected parents:" << endl;
  for (auto affected : affected_parents) {
    cout << *affected << endl;
  }
  cout << "Gaps:" << endl;
  for (auto gap : gaps) {
    cout << "[" << gap.first << ", " << gap.second << ")" << endl;
  }
#endif // _DEBUG

  /*
   * HANDLE GAPS
   */
  for (std::pair<arg_real_t, arg_real_t> gap : gaps) {
    node->add_parent(gap.first, gap.second, *ancestor);
    ++num_edges_cnt;
  }

  /*
   * HANDLE AFFECTED PARENTS
   */
  for (ARGEdge* edge : affected_parents) {
    // We currently don't allow threading to a time if there is a node there already
    // TODO: think about how to modify the equality condition
    if (edge->parent->height == ancestor->height) {
      throw std::logic_error(
          THROW_LINE("Threading where there is an existing node is not yet implemented"));
    }

    // Case A: parent height is less than ancestor height. Easy to deal with!
    if (edge->parent->height < ancestor->height) {
      arg_real_t left = std::max(entry_start, edge->start);
      arg_real_t right = std::min(entry_end, edge->end);
      entries.push(AncestorEntry(edge->parent, ancestor, left, right));
      continue;
    }

    // Case B: parent height is greater than ancestor height
    if (edge->start >= entry_start && edge->end <= entry_end) {
      // Case B1: parent edge is contained within ancestor entry
      //                              -P3--
      //     --------currentAncestor-----------
      //
      //------------currentIndividual-1----2----------------
      //
      // Actions:
      // (A) currentAncestor is new parent of currentIndividual in range 1-2.
      // (B) P3=currentAffectedParent is parent of currentAncestor in range 1-2.
      arg_real_t point1 = edge->start;
      arg_real_t point2 = edge->end;
      ARGNode* parent3 = edge->parent;
      // Note: the next line causes edge to be deleted! Make sure to get
      // out everything we need first
      node->remove_parent(point1);
      --num_edges_cnt;
      node->add_parent(point1, point2, *ancestor);
      ++num_edges_cnt;
      entries.push(AncestorEntry(ancestor, parent3, point1, point2));
    }
    else if (edge->start < entry_start && edge->end <= entry_end) {
      // Case B2: parent edge straddles left boundary of ancestor entry
      // ---P1---
      //     --------currentAncestor-----------
      //
      //-1---2---3-----currentIndividual-----------------------
      // Actions:
      // (A) P1=firstAffectedParent now only parent of currentIndividual in 1-2, instead of 1-3.
      // (B) currentAncestor is parent of currentIndividual in 2-3.
      // (C) P1 parent of currentAncestor in 2-3.
      arg_real_t point1 = edge->start;
      arg_real_t point2 = entry_start;
      arg_real_t point3 = edge->end;
      ARGNode* parent1 = edge->parent;
      node->update_parent_end(point1, point2);
      node->add_parent(point2, point3, *ancestor);
      ++num_edges_cnt;
      entries.push(AncestorEntry(ancestor, parent1, point2, point3));
    }
    else if (edge->start >= entry_start && edge->end > entry_end) {
      // Case B3: parent edge straddles right boundary of ancestor entry
      //                                    ----P5---
      //     --------currentAncestor-----------
      //
      //---------------currentIndividual----1--2-----3---------
      // Actions:
      // (A) P5=lastAffectedParent now only parent of currentIndividual in 2-3, instead of 1-3.
      // (B) currentAncestor is parent of currentIndividual in 1-2, where P5=lastAffectedParent was
      //     parent.
      // (C) P5 parent of currentAncestor in 1-2.
      arg_real_t point1 = edge->start;
      arg_real_t point2 = entry_end;
      // real point3 = edge->end;
      ARGNode* parent5 = edge->parent;
      node->update_parent_start(point1, point2);
      node->add_parent(point1, point2, *ancestor);
      ++num_edges_cnt;
      entries.push(AncestorEntry(ancestor, parent5, point1, point2));
    }
    else {
      // Case B4: parent edge straddles both boundaries of ancestor entry
      // ---P1----------------------------------------
      //     --------currentAncestor-----------
      //
      //-1---2---------currentIndividual-------3------4----------
      // Actions:
      // (A) P1=firstAffectedParent is now only parent of currentIndividual between 1-2 AND 3-4.
      // (B) currentAncestor is parent of currentIndividual between 2-3.
      // (C) P1=firstAffectedParent is parent of currentAncestor between 2-3.
      assert(edge->start < entry_start && edge->end > entry_end);
      arg_real_t point1 = edge->start;
      arg_real_t point2 = entry_start;
      arg_real_t point3 = entry_end;
      arg_real_t point4 = edge->end;
      ARGNode* parent1 = edge->parent;
      node->update_parent_end(point1, point2);
      node->add_parent(point3, point4, *parent1);
      node->add_parent(point2, point3, *ancestor);
      num_edges_cnt += 2;
      entries.push(AncestorEntry(ancestor, parent1, point2, point3));
    }
  }
}

const map<arg_real_t, Site>& ARG::get_mutation_sites() const
{
  update_mutation_sites();
  return mutation_sites;
}

void ARG::update_mutation_sites() const
{
  // This map will be out od date if we have modified the mutation vector
  if (!mutation_sites_up_to_date) {
    mutation_sites.clear();
    // Re-create the map by iterating over all mutations
    for (const auto& mutation : mutations) {
      const arg_real_t position = mutation->position;
      const auto site_iter = mutation_sites.try_emplace(mutation_sites.end(), position);
      Site& site = site_iter->second;
      site.add_mutation(mutation.get());
    }
    // The map is now definitely up-to-date
    mutation_sites_up_to_date = true;
  }
}

const std::vector<arg_real_t>& ARG::get_site_positions() const
{
  update_site_positions();
  return site_positions;
}

void ARG::update_site_positions() const
{
  // This vector may be out od date if we have modified the mutation vector
  if (!site_positions_up_to_date) {
    // Can't use mutation_sites directly as it might be out-of-date
    const auto& mutation_sites_updated = get_mutation_sites();

    site_positions.clear();
    site_positions.reserve(mutation_sites_updated.size());

    for (const auto& [key, val] : mutation_sites_updated) {
      site_positions.push_back(key);
    }

    site_positions_up_to_date = true;
  }
}
