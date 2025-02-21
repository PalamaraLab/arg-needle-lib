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

/* Higher-level ARG utilities that use public members of ARGNode / ARG

A block of code in this file is copied and modified from the
BOLT-LMM_v2.3.2 software developed by Po-Ru Loh and released under
the GNU General Public License v3.0 (GPLv3).

The license file can be found at 3rd_party/BOLT-LMM_v2.3.2/license.txt
from the root of this repository.
 */

#include "arg_utils.hpp"
#include "arg.hpp"
#include "arg_edge.hpp"
#include "arg_node.hpp"
#include "constants.hpp"
#include "descendant_list.hpp"
#include "deserialization_params.hpp"
#include "file_utils.hpp"
#include "random_utils.hpp"
#include "utils.hpp"

#include <algorithm>
#include <boost/math/distributions/chi_squared.hpp>
#include <chrono>
#include <cmath>
#include <fstream>
#include <future>
#include <functional>
#include <iostream>
#include <list>
#include <map>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <typeinfo>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using std::cout;
using std::deque;
using std::endl;
using std::ifstream;
using std::list;
using std::map;
using std::pair;
using std::set;
using std::stack;
using std::string;
using std::tuple;
using std::unordered_map;
using std::unordered_set;
using std::vector;

using Clock = std::chrono::high_resolution_clock;

using boost::dynamic_bitset;

namespace arg_utils {

unsigned validate_parallel_tasks(const unsigned num_tasks) {
    unsigned recommended_max = anl::get_default_concurrency();

    if (num_tasks == 0u) {
        std::cerr << "Warning: can't set num_tasks to 0: setting to 1\n";
        return 1u;
    }

    if (num_tasks > recommended_max) {
        std::cerr << "Warning: recommended max num_tasks is " << recommended_max
                  << ": you are using requesting " << num_tasks << '\n';
    }

    return num_tasks;
}

arg_real_t local_volume_single(const ARG &arg, arg_real_t min_pos, arg_real_t max_pos) {
    arg_real_t volume = 0;
    arg_utils::visit_branches(
            arg,
            [&volume](const DescendantList &desc_list, const DescendantList &_ignored,
                      const ARGNode *parent, const ARGNode *child, arg_real_t start, arg_real_t end) {
                (void) desc_list;
                (void) _ignored;
                arg_real_t branch_volume = (parent->height - child->height) * (end - start);
                volume += branch_volume;
            },
            min_pos, max_pos);
    return volume;
}


// Compute the number of lineages at a position at a certain height
//
// position: physical position along the chromosome
// height: height at which open lineages are counted
int num_lineages(const ARG& arg, arg_real_t position, arg_real_t height) {
  if (arg.roots.empty()) {
    throw std::logic_error(THROW_LINE("Call populate_children_and_roots() first."));
  }

  if ((position >= arg.end) || (position < arg.start)) {
    throw std::logic_error(THROW_LINE("Position out of range of ARG"));
  }
  // Get the root that overlap the position
  auto root = arg.root_at(position);
  // If height is at or above the root, there is only one open lineage
  if (height >= root->node->height) {
    return 1;
  }
  stack<const ARGEdge*> entries;
  // Get all edges that have parent / child overlapping height
  for (auto const& edge_entry : root->node->children_at(position)) {
    entries.push(edge_entry);
  }
  int n_lineages = 0;
  while (!entries.empty()) {
    const ARGEdge* edge_entry = entries.top();
    entries.pop();
    if ((edge_entry->child->height <= height) && (edge_entry->parent->height > height)) {
      ++n_lineages;
    }
    for (auto const& edge : edge_entry->child->children_at(position)) {
      // only push if top of edge is above height
      if (edge->parent->height >= height) {
        entries.push(edge);
      }
    }
  }
  return n_lineages;
}

// Trim ARG to keep nodes and edges with range within the interval [trim_start, trim_end). Start and
// end positions do not include the offset.
//
// Borrow the deserialisation constructor for now...
ARG trim_arg(ARG& arg, arg_real_t trim_start, arg_real_t trim_end) {
  // ensure a meaningful trim
  if (trim_start < arg.start) {
    throw std::logic_error(THROW_LINE("trim start is before ARG start"));
  }
  if (trim_end > arg.end) {
    throw std::logic_error(THROW_LINE("trim end is after ARG end"));
  }
  if (trim_start > trim_end) {
    throw std::logic_error(THROW_LINE("trim start is after trim end"));
  }
  vector<int> node_is_in_range(arg.num_nodes(), 0);
  vector<std::array<int, 2>> edge_ids;
  edge_ids.reserve(arg.num_edges());
  vector<std::array<double, 2>> edge_ranges;
  edge_ranges.reserve(arg.num_edges());

  for (auto const& map_entry : arg.arg_nodes) {
    auto node = map_entry.second.get();
    // when the interval [node->start, node->end) overlaps [trim_start, trim_end)
    if (node->start < trim_end && node->end > trim_start) {
      for (auto const& node_map_entry : node->parents) {
        auto edge = node_map_entry.second.get();
        // a node is in range iff it has in range edges
        if (edge->start < trim_end && edge->end > trim_start) {
          std::array<int, 2> edge_id{edge->child->ID, edge->parent->ID};
          std::array<double, 2> edge_range{std::max(edge->start - trim_start, arg_real_t(0)),
                                        std::min(edge->end, trim_end) - trim_start};
          node_is_in_range[edge->child->ID] = 1;
          node_is_in_range[edge->parent->ID] = 1;
          edge_ids.emplace_back(edge_id);
          edge_ranges.emplace_back(edge_range);
        }
      }
    }
  }

  edge_ids.shrink_to_fit();
  edge_ranges.shrink_to_fit();
  int num_nodes_in_range = std::count(node_is_in_range.begin(), node_is_in_range.end(), 1);
  vector<int> reassigned_node_id(arg.num_nodes(), 0);
  // new node id for a node is the number of in range nodes before it -- should preserve order and
  // begin with leaf nodes
  std::partial_sum(node_is_in_range.begin(), node_is_in_range.end(), reassigned_node_id.begin());
  for (auto& elem : reassigned_node_id) {
    elem -= 1; // to make ids 0-based index
  }

  for (auto& elem : edge_ids) {
    int new_child_id = reassigned_node_id[elem[0]];
    int new_parent_id = reassigned_node_id[elem[1]];
    std::array<int, 2> new_ids{new_child_id, new_parent_id};
    elem = new_ids;
  }

  // assemble all ARG data for the deserialisation constructor
  vector<double> node_heights(num_nodes_in_range, 0);
  vector<std::array<double, 2>> node_bounds(num_nodes_in_range);
  vector<uint8_t> is_sample(num_nodes_in_range, false);

  for (int old_node_id = 0; old_node_id < node_is_in_range.size(); old_node_id++) {
    auto node = arg.arg_nodes.at(old_node_id).get();
    if (node_is_in_range[node->ID]) {
      int new_node_id = reassigned_node_id[node->ID];
      is_sample[new_node_id] = arg.is_leaf(old_node_id);
      node_heights[new_node_id] = node->height;
      std::array<double, 2> current_node_bound{
          std::max(node->start - trim_start, arg_real_t(0)), std::min(node->end, trim_end) - trim_start};
      node_bounds[new_node_id] = current_node_bound;
    }
  }

  DeserializationParams dp{};
  dp.chromosome = arg.chromosome;
  dp.start = std::max(arg.start - trim_start, arg_real_t(0));
  dp.end = std::min(trim_end, arg.end) - trim_start;
  dp.offset = arg.offset + trim_start;
  dp.threaded_samples = arg.threaded_samples;
  dp.reserved_samples = arg.reserved_samples;
  dp.num_nodes = num_nodes_in_range;
  ARG trimmed_arg = ARG(dp);

  trimmed_arg.deserialize_add_nodes(node_heights, is_sample, node_bounds);
  node_heights.clear();
  is_sample.clear();
  node_bounds.clear();
  trimmed_arg.deserialize_add_edges(edge_ids, edge_ranges);
  edge_ids.clear();
  edge_ranges.clear();

  if (arg.num_mutations() > 0) {
    int num_sites_in_range = 0;
    vector<arg_real_t> positions;
    positions.reserve(arg.get_num_sites());
    vector<arg_real_t> heights;
    heights.reserve(arg.get_num_sites());
    vector<std::array<int, 2>> mutation_edge_ids;
    mutation_edge_ids.reserve(arg.get_num_sites());

    for (auto const& entry : arg.get_mutations()) {
      auto mutation = entry.get();
      if (trim_start <= mutation->position && mutation->position < trim_end) {
        num_sites_in_range++;
        positions.emplace_back(mutation->position - trim_start);
        heights.emplace_back(mutation->height);
        std::array<int, 2> mutation_edge_id{
            reassigned_node_id[mutation->edge->child->ID], reassigned_node_id[mutation->edge->parent->ID]};
        mutation_edge_ids.emplace_back(mutation_edge_id);
      }
    }
    vector<int> site_ids(num_sites_in_range);
    // fill site_ids as default
    std::iota(site_ids.begin(), site_ids.end(), 0);
    trimmed_arg.deserialize_add_mutations(positions, heights, site_ids, mutation_edge_ids);
  }
  return trimmed_arg;
}

ARG arg_from_ts_files(string node_file_name, string edge_file_name) {
  ifstream node_file(node_file_name);
  ifstream edge_file(edge_file_name);

  if (!node_file.good()) {
    throw std::invalid_argument(THROW_LINE("Invalid node file path: " + node_file_name));
  }
  if (!edge_file.good()) {
    throw std::invalid_argument(THROW_LINE("Invalid edge file path: " + edge_file_name));
  }

  string col;
  for (int i = 0; i < 6; ++i) {
    node_file >> col;
  }
  for (int i = 0; i < 4; ++i) {
    edge_file >> col;
  }

  arg_real_t start = 0, end = 0;
  vector<arg_real_t> node_heights;
  deque<bool> is_sample;
  vector<std::pair<int, int>> edge_ids;
  vector<std::pair<arg_real_t, arg_real_t>> edge_ranges;

  int a;
  bool node_is_sample;
  arg_real_t node_height;
  while (node_file >> a >> node_is_sample >> node_height >> a >> a) {
    node_heights.push_back(node_height);
    is_sample.push_back(node_is_sample);
  }

  int parent, child;
  arg_real_t left, right;
  while (edge_file >> left >> right >> parent >> child) {
    if (start == end) {
      start = left;
      end = right;
    }
    edge_ranges.push_back(std::make_pair(left, right));
    edge_ids.push_back(std::make_pair(child, parent));
  }

  return ARG(start, end, node_heights, is_sample, edge_ids, edge_ranges);
}

// To be used for calling from Python without having to convert unordered_map
bool visit_identical(const ARG& arg, arg_real_t rel_tol, arg_real_t abs_tol, bool timing,
                     bool verbose) {
  std::chrono::time_point<Clock> last_time, curr_time;
  last_time = Clock::now();
  unordered_map<DescendantList, arg_real_t, DescendantListHash> map1 = efficient_visit(arg);
  curr_time = Clock::now();
  if (timing) {
    arg_real_t split_time =
        std::chrono::duration_cast<std::chrono::milliseconds>(curr_time - last_time).count();
    cout << "Done with efficient visit in " << (split_time / 1000) << " seconds." << endl;
  }
  last_time = curr_time;
  unordered_map<string, arg_real_t> map2 = bitset_volume_map(arg);
  curr_time = Clock::now();
  if (timing) {
    arg_real_t split_time =
        std::chrono::duration_cast<std::chrono::milliseconds>(curr_time - last_time).count();
    cout << "Done with long visit in " << (split_time / 1000) << " seconds." << endl;
  }
  last_time = curr_time;
  bool result = true;
  int n = arg.leaf_ids.size();
  if (map1.size() != map2.size()) {
    cout << "Type 1: " << map1.size() << " " << map2.size() << endl;
    result = false;
  }
  else {
    for (auto const& entry : map1) {
      string key = entry.first.to_string();
      if (map2.find(key) == map2.end()) {
        cout << "Type 2 " << key << endl;
        result = false;
        break;
      }
      if (!utils::is_close(map2.at(key), entry.second, rel_tol, abs_tol)) {
        cout << "Type 3" << endl;
        cout << key << " " << entry.second << " " << map2.at(key) << " ";
        cout << utils::relative_ratio(map2.at(key), entry.second) << " ";
        cout << std::fabs(map2.at(key) - entry.second) << endl;
        result = false;
        break;
      }
    }
  }
  curr_time = Clock::now();
  if (timing) {
    arg_real_t split_time =
        std::chrono::duration_cast<std::chrono::milliseconds>(curr_time - last_time).count();
    cout << "Done with comparison in " << (split_time / 1000) << " seconds." << endl;
  }

  if (verbose) {
    cout << "Efficient visit result" << endl;
    for (auto const& entry : map1) {
      string key = entry.first.to_string();
      cout << key << " " << entry.second << endl;
    }
    cout << endl;
    cout << "Inefficient visit result" << endl;
    for (auto const& entry : map2) {
      cout << entry.first << " " << entry.second << endl;
    }
    cout << endl;
  }

  return result;
}

// To be used for calling from Python without having to convert unordered_map
void time_efficient_visit(const ARG& arg, bool timing) {
  std::chrono::time_point<Clock> last_time, curr_time;
  last_time = Clock::now();
  vector<arg_real_t> thing = allele_frequency_spectrum_volume(arg); // this is a lot faster
  curr_time = Clock::now();
  if (timing) {
    arg_real_t split_time =
        std::chrono::duration_cast<std::chrono::milliseconds>(curr_time - last_time).count();
    cout << "Done with efficient visit in " << (split_time / 1000) << " seconds." << endl;
  }
}

// used to compare that our visit is working correctly, outputs an unordered_map which
// we compare against the slow version
unordered_map<DescendantList, arg_real_t, DescendantListHash> efficient_visit(const ARG& arg) {
  unordered_map<DescendantList, arg_real_t, DescendantListHash> total_volume_map;
  visit_branches(arg, [&total_volume_map](const DescendantList& desc_list,
                                          const DescendantList& _ignored, const ARGNode* parent,
                                          const ARGNode* child, arg_real_t start, arg_real_t end) {
    (void) _ignored;
    arg_real_t branch_volume = (parent->height - child->height) * (end - start);
    arg_real_t& total_volume_entry =
        total_volume_map[desc_list]; // creates if not present, only hashes once
    total_volume_entry += branch_volume;
  });
  return total_volume_map; // return by value
}

// Naive way to get volume - does not use DescendentLists at all
arg_real_t total_volume(const ARG& arg) {
  arg_real_t volume = 0;
  for (auto const& map_entry : arg.arg_nodes) {
    ARGNode* node = map_entry.second.get();
    for (auto const& edge_entry : node->parents) {
      // iterate through every edge
      const ARGEdge* edge = edge_entry.second.get();
      const ARGNode* child = edge->child;
      const ARGNode* parent = edge->parent;
      arg_real_t start = edge->start;
      arg_real_t end = edge->end;
      volume += (parent->height - child->height) * (end - start);
    }
  }
  return volume;
}

arg_real_t local_volume(const ARG &arg, std::optional<arg_real_t> min_pos, std::optional<arg_real_t> max_pos,
                        std::optional<unsigned> num_tasks) {

    const unsigned valid_num_tasks = validate_parallel_tasks(num_tasks.value_or(anl::get_default_concurrency()));
    const arg_real_t valid_min_pos = min_pos.value_or(arg.start);
    const arg_real_t valid_max_pos = max_pos.value_or(arg.end);

    if (valid_num_tasks == 1u) {
        return local_volume_single(arg, valid_min_pos, valid_max_pos);
    }

    std::vector<std::future<arg_real_t>> results;
    const arg_real_t step_size = (valid_max_pos - valid_min_pos) / valid_num_tasks;

    for (unsigned i = 0; i < valid_num_tasks; ++i) {
        const arg_real_t lo = valid_min_pos + i * step_size;
        const arg_real_t hi = valid_min_pos + (i + 1u) * step_size;
        results.push_back(std::async(std::launch::async, local_volume_single, std::cref(arg), lo, hi));
    }

    arg_real_t res = 0.0;
    for (auto &fut: results) {
        res += fut.get();
    }

    return res;
}

// result[i] is the total volume of branches that have i descendants
vector<arg_real_t> allele_frequency_spectrum_volume(const ARG& arg) {
  size_t n = arg.leaf_ids.size();
  // 0 and n are actually impossible but may want to access those values
  vector<arg_real_t> result(n + 1, 0.);
  visit_branches(arg, [&result](const DescendantList& desc_list, const DescendantList& _ignored,
                                const ARGNode* parent, const ARGNode* child, arg_real_t start,
                                arg_real_t end) {
    (void) _ignored;
    arg_real_t branch_volume = (parent->height - child->height) * (end - start);
    size_t num_descendants = desc_list.num_values();
    result[num_descendants] += branch_volume;
  });
  return result;
}

// Since we are going to be returning everything, we do the simple memory-inefficient algorithm
//
// mu: mutation rate
// random_seed: seed for generating mutations, uses time to seed when random_seed = 0
map<arg_real_t, string> generate_mutations_map(const ARG& arg, arg_real_t mu,
                                               unsigned random_seed) {
  map<arg_real_t, string> vcf_data;
  if (random_seed == 0) {
    random_seed = std::chrono::system_clock::now().time_since_epoch().count();
  }
  std::mt19937 generator(random_seed);
  arg_utils::visit_branches(
      arg, [&mu, &vcf_data, &generator](const DescendantList& desc_list,
                                        const DescendantList& _ignored, const ARGNode* parent,
                                        const ARGNode* child, arg_real_t start, arg_real_t end) {
        (void) _ignored;
        arg_real_t branch_volume = (parent->height - child->height) * (end - start);
        unsigned num_mutations = random_utils::generate_poisson_rv(generator, mu * branch_volume);
        for (unsigned i = 0; i < num_mutations; ++i) {
          arg_real_t pos = random_utils::generate_uniform_rv(generator, start, end);
          vcf_data[pos] = desc_list.to_string();
        }
      });
  return vcf_data;
}

// Sample mutations and output them in order using a priority queue
//
// Our visit_branches function will visit branches in nondecreasing order by
// branch start. However you might visit two branches in this order and sample
// a later mutation on the first branch, then an earlier mutation on the second
// branch. Therefore, the priority queue is needed to sort the mutations. We
// flush out any mutations that are before our current branch start position.
//
// mu: mutation rate
// file_root: when empty, we output to stdout instead of writing to file
// random_seed: seed for generating mutations, uses time to seed when random_seed = 0
void generate_mutations_and_write(const ARG& arg, arg_real_t mu, string file_root,
                                  unsigned random_seed) {
  if (random_seed == 0) {
    random_seed = std::chrono::system_clock::now().time_since_epoch().count();
  }
  std::mt19937 generator(random_seed);
  size_t n = arg.leaf_ids.size();

  // open output streams if required and write samples
  file_utils::AutoGzOfstream haps_out_fstream, samples_out_fstream;
  if (file_root != "") {
    haps_out_fstream.openOrExit(file_root + ".haps.gz");
    samples_out_fstream.openOrExit(file_root + ".samples");
    if (n % 2 == 1) {
      std::cerr << "ERROR: trying to write genotype data for odd number of haploid samples."
                << endl;
      exit(1);
    }
    samples_out_fstream << "ID_1 ID_2 missing\n0 0 0" << endl;
    for (size_t i = 1; i <= n / 2; ++i) {
      samples_out_fstream << i << " " << i << " 0" << endl;
    }
  }
  else {
    cout << "position time descendants AF branch_start branch_end child_height parent_height bitset"
         << endl;
  }

  std::priority_queue<pair<int, string>, vector<pair<int, string>>, std::greater<pair<int, string>>>
      pq;
  arg_utils::visit_branches(arg, [&haps_out_fstream, &mu, &generator, n, file_root, &pq](
                                     const DescendantList& desc_list,
                                     const DescendantList& _ignored, const ARGNode* parent,
                                     const ARGNode* child, arg_real_t start, arg_real_t end) {
    (void) _ignored;
    // flush mutations that are before `start`
    while (!pq.empty() && pq.top().first < start) {
      if (file_root != "") {
        haps_out_fstream << pq.top().second << endl;
      }
      else {
        cout << pq.top().second << endl;
      }
      pq.pop();
    }
    arg_real_t branch_volume = (parent->height - child->height) * (end - start);
    size_t num_descendants = desc_list.num_values();
    unsigned num_mutations = random_utils::generate_poisson_rv(generator, mu * branch_volume);
    arg_real_t allele_frequency = ((arg_real_t) num_descendants) / ((arg_real_t) n);
    for (unsigned i = 0; i < num_mutations; i++) {
      // sample time
      arg_real_t mutation_time =
          random_utils::generate_uniform_rv(generator, child->height, parent->height);
      // sample position
      arg_real_t pos = random_utils::generate_uniform_rv(generator, start, end);
      int int_pos = (int) floor(pos);
      // put in priority_queue
      std::stringstream ss;
      string bitset_string = desc_list.to_bitset_string();
      if (file_root != "") {
        ss << "1 1:" << int_pos << "_" << mutation_time << " " << int_pos << " 1 2";
        for (unsigned i = 0; i < bitset_string.size(); i++) {
          ss << ' ' << bitset_string[i];
        }
      }
      else {
        ss << int_pos << " " << mutation_time << " " << num_descendants << " " << allele_frequency;
        ss << " " << start << " " << end << " " << child->height << " " << parent->height;
        ss << " " << bitset_string;
      }
      pq.push(std::make_pair(int_pos, ss.str()));
    }
  });
  // flush all remaining mutations
  while (!pq.empty()) {
    if (file_root != "") {
      haps_out_fstream << pq.top().second << endl;
    }
    else {
      cout << pq.top().second << endl;
    }
    pq.pop();
  }
  if (file_root != "") {
    haps_out_fstream.close();
    samples_out_fstream.close();
  }
}

// Naive sampling procedure for sampling mutations
// Key points in algorithm:
// 1. iterate through edges of the ARG randomly and keep them sampled in the ARG as a vector
// 2. sort mutations based on their position
//
// mu: mutation rate
// random_seed: seed for generating mutations, uses time to seed when random_seed = 0
//
// If it is known roughly how many mutations will be generated, you can provide a hint so that an appropraite amount of
// memory is reserved in the mutation vector.
//
// The mutations will always be sorted by physical position.
void generate_mutations(ARG& arg, arg_real_t mu, unsigned random_seed, std::size_t num_mutations_hint) {
  if (random_seed == 0) {
    random_seed = std::chrono::system_clock::now().time_since_epoch().count();
  }
  std::mt19937 generator(random_seed);
  // Clear all previous mutations
  arg.clear_mutations();

  std::vector<Mutation> new_mutations;
  if (num_mutations_hint > 0ul) {
    new_mutations.reserve(num_mutations_hint);
  }

  // Start the traversal (iterate through all the nodes & edges)
  for (auto const& map_entry : arg.arg_nodes) {
    const ARGNode* node = map_entry.second.get();
    for (auto const& edge_entry : node->parents) {
      // iterate through every edge
      ARGEdge* edge = edge_entry.second.get(); // can't be const ARGEdge*
      const ARGNode* child = edge->child;
      const ARGNode* parent = edge->parent;
      const arg_real_t start = edge->start;
      const arg_real_t end = edge->end;
      const arg_real_t branch_volume = (parent->height - child->height) * (end - start);
      // now sample a Poisson rv
      auto num_mutations =
          static_cast<unsigned>(random_utils::generate_poisson_rv(generator, mu * branch_volume));
      if (num_mutations > 0) {
        // Possible to have more than one mutation on an edge!
        for (unsigned i = 0; i < num_mutations; ++i) {
          // Uniformly sample a position
          const arg_real_t pos = random_utils::generate_uniform_rv(generator, start, end);
          // Uniformly sample an "age" on the branch
          const arg_real_t height =
              random_utils::generate_uniform_rv(generator, child->height, parent->height);

          new_mutations.emplace_back(edge, pos, height);
        }
      }
    }
  }

  // Sort the new mutations and add them to the ARG. Sorting is necessary to improve performance, as
  // add_mutation inserts into the vector, if necessary
  arg.reserve_n_mutations(new_mutations.size());
  std::sort(new_mutations.begin(), new_mutations.end());
  for (const auto& mut : new_mutations) {
    arg.add_mutation(mut.edge, mut.position, mut.height, -1, false);
  }
  arg.update_site_data_structures();
}

// Generate exactly M mutations with the following steps:
// * compute ARG volume
// * uniformly sample position of the M mutations from 0 to volume
// * re-traverse the ARG and drop the mutations when sampled volume is reached
//
// m: number of mutations
// random_seed: seed for generating mutations, uses time to seed when random_seed = 0
//
// The mutations will always be sorted by physical position.
void generate_m_mutations(ARG& arg, int m, unsigned random_seed) {
  if (m < 1) {
    throw std::logic_error(THROW_LINE("Need to simulate a positive number of mutations!"));
  }
  if (random_seed == 0) {
    random_seed = std::chrono::system_clock::now().time_since_epoch().count();
  }
  std::mt19937 generator(random_seed);

  // Clear all previous mutations and reserve space
  arg.clear_mutations();

  // Reserve enough space so that the vector never needs resizing.
  std::vector<Mutation> new_mutations;
  new_mutations.reserve(m);

  // Calculate the total volume of the ARG
  const arg_real_t total_volume = arg_utils::total_volume(arg);

  // Sample uniformly along the volume and sort
  vector<arg_real_t> volume_pos;
  volume_pos.reserve(m);
  for (int i = 0; i < m; i++) {
    volume_pos.emplace_back(random_utils::generate_uniform_rv(generator, 0.0, total_volume));
  }
  std::sort(volume_pos.begin(), volume_pos.end());

  // Keep track of the aggregated volume
  arg_real_t agg_volume = 0.0;
  size_t mut_idx = 0; // indicator of the mutation to simulate
  for (auto const& map_entry : arg.arg_nodes) {
    const ARGNode* node = map_entry.second.get();
    for (auto const& edge_entry : node->parents) {
      // iterate through every edge
      ARGEdge* edge = edge_entry.second.get(); // can't be const ARGEdge*
      const ARGNode* child = edge->child;
      const ARGNode* parent = edge->parent;
      const arg_real_t start = edge->start;
      const arg_real_t end = edge->end;
      agg_volume += (parent->height - child->height) * (end - start);
      while ((mut_idx < volume_pos.size()) && (agg_volume > volume_pos[mut_idx])) {
        const arg_real_t pos =
            random_utils::generate_uniform_rv(generator, start, end); // sample basepair position
        const arg_real_t height =
            random_utils::generate_uniform_rv(generator, child->height,
                                              parent->height); // sample age of mutation

        new_mutations.emplace_back(edge, pos, height);
        mut_idx++;
      }
    }
  }

  // Sort the new mutations and add them to the ARG. Sorting is necessary to improve performance, as
  // add_mutation inserts into the vector, if necessary
  arg.reserve_n_mutations(new_mutations.size());
  std::sort(new_mutations.begin(), new_mutations.end());
  for (const auto& mut : new_mutations) {
    arg.add_mutation(mut.edge, mut.position, mut.height, -1, false);
  }
  arg.update_site_data_structures();

  // Check that we have exactly m mutations!
  assert(arg.num_mutations() == m);
}

// Iterate through all the mutations in a physical range and return an Eigen matrix
// Note - this uses pre-existing mutations and returns bitsets
// from_pos (default = -inf) is the left boundary (physical position) of the range
// to_pos (default = +inf) is the right boundary (physical position) of the range
// include_left: whether to include mutations found at from_pos (>= instead of >)
// include_right: whether to include mutations found at to_pos (<= instead of <)
MatXui get_mutations_matrix(const ARG& arg, arg_real_t from_pos, arg_real_t to_pos,
                             bool include_left, bool include_right) {
  if (arg.get_mutations().empty()) {
    throw std::logic_error(
        THROW_LINE("Need pre-existing mutations to run! Try generate_mutations_and_keep."));
  }
  if (from_pos > to_pos) {
    throw std::logic_error(
        THROW_LINE("Left range position is higher than right range position."));
  }
  size_t n = arg.leaf_ids.size();
  if (arg.get_mutations().front()->position > to_pos || arg.get_mutations().back()->position < from_pos) {
      std::cerr << "Warning: no mutations overlapping required range.\n";
      MatXui mat(0, n);
      return mat;
  }
  int index_start = arg.get_idx_of_first_mutation_right_of(from_pos, include_left);
  int index_end = arg.get_idx_of_first_mutation_left_of(to_pos, include_right);
  // TODO: if max_time and min_time are used, might have less than (index_end - index_start)
  // mutations
  MatXui mat(index_end - index_start + 1, n);
  mat.setZero();
  int cnt = 0;
  // Visit mutations
  arg_utils::visit_mutations(
      arg, [&arg, &n, &mat, &cnt](DescendantList& desc_list, const Mutation* mutation) {
        size_t num_descendants = desc_list.num_values();
        arg_real_t allele_frequency = ((arg_real_t) num_descendants) / ((arg_real_t) n);
        arg_real_t maf = (allele_frequency <= 0.5) ? allele_frequency : 1 - allele_frequency;
        // time is handled within visit_mutations
        int int_pos = arg.offset + (int) mutation->position;
        for (int i : desc_list.values()) {
          mat(cnt, i) = 1;
        }
        cnt++;
        // also setting min_time = 0, max_time = std::numeric_limits<double>::infinity()
        // to access last two optional parameters
      },
      -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), index_start, index_end);
  return mat;
}

// Iterate through all the mutations and print out bitsets to files
// Note - this uses pre-existing mutations
void write_mutations_to_haps(const ARG& arg, string file_root, arg_real_t min_maf,
                             arg_real_t max_maf, arg_real_t min_time, arg_real_t max_time) {
  if (arg.get_mutations().empty()) {
    throw std::logic_error(
        THROW_LINE("Need pre-existing mutations to run! Try generate_mutations_and_keep."));
  }
  file_utils::AutoGzOfstream haps_out_fstream, samples_out_fstream;
  size_t n = arg.leaf_ids.size();
  if (file_root != "") {
    haps_out_fstream.openOrExit(file_root + ".haps.gz");
    samples_out_fstream.openOrExit(file_root + ".samples");
    if (n % 2 == 1) {
      std::cerr << "ERROR: trying to write genotype data for odd number of haploid samples."
                << endl;
      exit(1);
    }
    samples_out_fstream << "ID_1 ID_2 missing\n0 0 0" << endl;
    for (size_t i = 1; i <= n / 2; ++i) {
      samples_out_fstream << i << " " << i << " 0" << endl;
    }
  }
  else {
    cout << "position time descendants AF child_height parent_height bitset" << endl;
  }

  // Visit mutations
  arg_utils::visit_mutations(
      arg,
      [&arg, &n, min_maf, max_maf, &haps_out_fstream, &file_root](
          DescendantList& desc_list, const Mutation* mutation) {
        size_t num_descendants = desc_list.num_values();
        arg_real_t allele_frequency = ((arg_real_t) num_descendants) / ((arg_real_t) n);
        arg_real_t maf = (allele_frequency <= 0.5) ? allele_frequency : 1 - allele_frequency;
        // time is handled within visit_mutations
        if (maf >= min_maf && maf < max_maf) {
          string bitset_string = desc_list.to_bitset_string();
          int int_pos = arg.offset + (int) mutation->position;
          std::stringstream ss;
          if (file_root != "") {
            ss << "1 1:" << int_pos << "_" << mutation->edge->child->height << " " << int_pos
               << " 1 2";
            for (size_t i = 0; i < bitset_string.size(); ++i) {
              ss << ' ' << bitset_string[i];
            }
            haps_out_fstream << ss.str() << endl;
          }
          else {
            cout << int_pos << " " << mutation->height << " " << num_descendants << " "
                 << allele_frequency;
            cout << " " << mutation->edge->child->height << " " << mutation->edge->parent->height;
            cout << " " << bitset_string << endl;
          }
        }
      },
      min_time, max_time);
  if (file_root != "") {
    haps_out_fstream.close();
    samples_out_fstream.close();
  }
}

// Slow visit function which we currently use for tests
unordered_map<string, arg_real_t> bitset_volume_map(const ARG& arg, bool verbose) {
  if (arg.roots.empty()) {
    throw std::logic_error(THROW_LINE("Call populate_children_and_roots() first."));
  }

  unordered_map<string, arg_real_t> total_volume_map;

  arg_real_t position = arg.start;
  auto comp_height = [](const ARGNode* a, const ARGNode* b) { return a->height > b->height; };

  if (verbose) {
    cout << "ARG visit:\n"
         << "phys_from, phys_to, time_from, time_to, edge_volume, descendants\n"
         << "-----------------------------------------------------------------------------" << endl;
  }
  while (position < arg.end) {
    const Root* current_root = arg.root_at(position);
    arg_real_t position_new = current_root->end;
    // top down visit
    list<int> top_down_list;
    top_down_list.push_back(current_root->node->ID);
    std::priority_queue<ARGNode*, vector<ARGNode*>, decltype(comp_height)> block_queue(comp_height);
    unordered_set<int>
        active_blocks; // create list of those from previous generation that should be removed
    block_queue.push(arg.arg_nodes.at(current_root->node->ID).get());
    active_blocks.insert(current_root->node->ID);
    while (!top_down_list.empty()) {
      const ARGNode* current_node = arg.arg_nodes.at(top_down_list.front()).get();
      top_down_list.pop_front();
      vector<ARGEdge*> children = current_node->children_at(position);
      for (const auto& childEdge : children) {
        position_new = std::min(position_new, childEdge->end);
        if (active_blocks.find(childEdge->child->ID) == active_blocks.end()) {
          block_queue.push(arg.arg_nodes.at(childEdge->child->ID).get());
          top_down_list.push_back(childEdge->child->ID);
          active_blocks.insert(childEdge->child->ID);
        }
      }
    }
    arg_real_t edge_width = position_new - position;
    // build the lists of descendants, which are stored in node_descendants_map[node_ID]
    unordered_map<int, DescendantListOld> node_descendants_map;
    for (int i = 0; i < arg.threaded_samples; i++) {
      DescendantListOld list(i);
      node_descendants_map.insert({i, list});
    }
    while (!block_queue.empty()) {
      const ARGNode* current_node = block_queue.top();
      block_queue.pop();
      DescendantListOld* myList = &node_descendants_map[current_node->ID];
      for (const auto& childEdge : current_node->children_at(position)) {
        DescendantListOld* child_list = &node_descendants_map[childEdge->child->ID];
        myList->add(*child_list);
        arg_real_t edge_height = childEdge->parent->height - childEdge->child->height;
        arg_real_t edge_volume = edge_width * edge_height;
        string str = child_list->to_string(arg.leaf_ids.size());
        arg_real_t& total_volume_entry =
            total_volume_map[str]; // creates if not present, only hashes once
        total_volume_entry += edge_volume;
        if (verbose) {
          cout << position << "\t" << position_new << "\t" << childEdge->child->height << "\t"
               << childEdge->parent->height << "\t" << edge_volume << "\t" << str << endl;
        }
      }
      node_descendants_map[current_node->ID] = *myList;
    }
    position = position_new;
    if (verbose) {
      cout << "-----------------------------------------------------------------------------"
           << endl;
    }
  }

  return total_volume_map; // return by value
}

tuple<string, arg_real_t> newick_subtree(const ARGNode& node, arg_real_t position,
                                         arg_real_t dist_from_parent, bool verbose) {
  vector<ARGEdge*> edges = node.children_at(position);
  arg_real_t first_break = node.end;
  string result;
  if (!edges.empty()) {
    result += "(";
    for (const ARGEdge* edge : edges) {
      string result_recurse;
      arg_real_t first_break_recurse;
      std::tie(result_recurse, first_break_recurse) =
          newick_subtree(*edge->child, position, node.height - edge->child->height, verbose);
      result += result_recurse;
      result += ",";
      first_break = std::min(std::min(first_break, first_break_recurse), edge->end);
    }
    result = result.substr(0, result.size() - 1) + ")";
  }
  result += std::to_string(node.ID);
  if ((verbose) && (dist_from_parent > 0)) {
    result += ":" + std::to_string((int) dist_from_parent);
  }
  return make_tuple(result, first_break);
}

string arg_to_newick(const ARG& arg, bool verbose) {
  std::ostringstream oss;
  for (arg_real_t root_start : arg.root_starts()) {
    const Root* root = arg.root_at(root_start);

    arg_real_t position = root->start;
    arg_real_t position_new;
    while (position < root->end) {
      string newick_string;
      std::tie(newick_string, position_new) = newick_subtree(*root->node, position, 0, verbose);
      oss << "Tree in [" << position << ", " << position_new << "): ";
      oss << newick_string << endl;
      position = position_new;
    }
  }
  return oss.str();
}

// Weighted average of *squared* L2 norm over vectors
arg_real_t tmrca_mse(const ARG& arg1, const ARG& arg2) {
  // check that extents are the same
  assert(arg1.start == arg2.start && arg1.end == arg2.end);

  // check that sample IDs are the same
  assert(arg1.leaf_ids.size() == arg2.leaf_ids.size());
  for (int leaf1 : arg1.leaf_ids) {
    assert(arg2.leaf_ids.find(leaf1) != arg2.leaf_ids.end());
  }

  vector<int> leaf_ids;
  for (int leaf1 : arg1.leaf_ids) {
    leaf_ids.push_back(leaf1);
  }

  // iterate over all pairs and accumulate the results
  arg_real_t result = 0;
  int n = leaf_ids.size();
  for (int i = 0; i < n; ++i) {
    cout << utils::current_time_string() << " TMRCA against sample " << i << endl;
    for (int j = i + 1; j < n; ++j) {
      arg_real_t pos = arg1.start;
      arg_real_t end1 = pos;
      arg_real_t end2 = pos;
      const ARGNode* mrca1;
      const ARGNode* mrca2;
      while (pos != arg1.end) {
        if (end1 == pos) {
          // get next MRCA
          std::tie(mrca1, end1) = arg1.mrca_with_end(leaf_ids[i], leaf_ids[j], pos);
        }
        if (end2 == pos) {
          // get next MRCA
          std::tie(mrca2, end2) = arg2.mrca_with_end(leaf_ids[i], leaf_ids[j], pos);
        }
        arg_real_t min_end = std::min(end1, end2);
        arg_real_t diff = mrca1->height - mrca2->height;
        result += diff * diff * (min_end - pos);
        pos = min_end;
      }
    }
  }

  // average and return
  arg_real_t factor = n * (n - 1) * (arg1.end - arg1.start);
  factor /= 2;
  result /= factor;
  return result;
}

// Weighted average of *squared* L2 norm over vectors, no square root
arg_real_t kc_topology(const ARG& arg1, const ARG& arg2) {
  // check that extents are the same
  assert(arg1.start == arg2.start && arg1.end == arg2.end);

  // check that sample IDs are the same
  assert(arg1.leaf_ids.size() == arg2.leaf_ids.size());
  for (int leaf1 : arg1.leaf_ids) {
    assert(arg2.leaf_ids.find(leaf1) != arg2.leaf_ids.end());
  }

  vector<int> leaf_ids;
  for (int leaf1 : arg1.leaf_ids) {
    leaf_ids.push_back(leaf1);
  }

  set<arg_real_t> bp1 = arg1.get_breakpoints();
  set<arg_real_t> bp2 = arg2.get_breakpoints();
  vector<arg_real_t> bp;
  std::set_union(bp1.begin(), bp1.end(), bp2.begin(), bp2.end(), std::inserter(bp, bp.begin()));
  bp.push_back(arg1.end);
  cout << utils::current_time_string() << " KC number of regions: " << bp.size() - 1 << endl;

  // iterate over all pairs and accumulate the results
  arg_real_t result = 0;
  int n = leaf_ids.size();
  for (size_t k = 0; k < bp.size() - 1; ++k) {
    if (k % 1000 == 0) {
      cout << utils::current_time_string() << " KC region " << k << endl;
    }
    arg_real_t distance = 0;
    arg_real_t position = bp[k];
    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {
        const ARGNode* mrca1 = arg1.mrca(leaf_ids[i], leaf_ids[j], position);
        const ARGNode* mrca2 = arg2.mrca(leaf_ids[i], leaf_ids[j], position);
        int count1 = 0;
        int count2 = 0;

        const ARGEdge* edge;
        edge = mrca1->parent_edge_at(position);
        while (edge != nullptr) {
          count1 += 1;
          mrca1 = edge->parent;
          edge = mrca1->parent_edge_at(position);
        }
        edge = mrca2->parent_edge_at(position);
        while (edge != nullptr) {
          count2 += 1;
          mrca2 = edge->parent;
          edge = mrca2->parent_edge_at(position);
        }
        distance += (count1 - count2) * (count1 - count2);
      }
    }
    result += distance * (bp[k + 1] - position);
  }

  // average and return, note that we don't divide by the number of pairs
  result /= arg1.end - arg1.start;
  return result;
}

// Average of *squared* L2 norm over vectors, KC and TMRCA
pair<arg_real_t, arg_real_t> metrics_stab(const ARG& arg1, const ARG& arg2, int num_stabs) {
  // check that extents are the same
  assert(arg1.start == arg2.start && arg1.end == arg2.end);

  // check that sample IDs are the same
  assert(arg1.leaf_ids.size() == arg2.leaf_ids.size());
  for (int leaf1 : arg1.leaf_ids) {
    assert(arg2.leaf_ids.find(leaf1) != arg2.leaf_ids.end());
  }

  vector<int> leaf_ids;
  for (int leaf1 : arg1.leaf_ids) {
    leaf_ids.push_back(leaf1);
  }

  arg_real_t span = arg1.end - arg1.start;
  arg_real_t phi = (sqrt(5) - 1.0) / 2.0;
  arg_real_t proportion = 0;
  vector<arg_real_t> stabs;
  for (int i = 0; i < num_stabs; ++i) {
    proportion = fmod(proportion + phi, 1.0);
    stabs.push_back(proportion * span + arg1.start);
  }

  // iterate over all pairs and accumulate the results
  arg_real_t kc_distance = 0;
  arg_real_t tmrca_mse = 0;
  int n = leaf_ids.size();
  cout << utils::current_time_string() << " Total to stab: " << stabs.size() << endl;
  for (size_t k = 0; k < stabs.size(); ++k) {
    if (k % 1000 == 0) {
      cout << utils::current_time_string() << " Metrics stab " << k << endl;
    }
    arg_real_t distance = 0;
    arg_real_t position = stabs[k];
    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {
        const ARGNode* mrca1 = arg1.mrca(leaf_ids[i], leaf_ids[j], position);
        const ARGNode* mrca2 = arg2.mrca(leaf_ids[i], leaf_ids[j], position);
        tmrca_mse += (mrca1->height - mrca2->height) * (mrca1->height - mrca2->height);
        int count1 = 0;
        int count2 = 0;

        const ARGEdge* edge;
        edge = mrca1->parent_edge_at(position);
        while (edge != nullptr) {
          count1 += 1;
          mrca1 = edge->parent;
          edge = mrca1->parent_edge_at(position);
        }
        edge = mrca2->parent_edge_at(position);
        while (edge != nullptr) {
          count2 += 1;
          mrca2 = edge->parent;
          edge = mrca2->parent_edge_at(position);
        }
        distance += (count1 - count2) * (count1 - count2);
      }
    }
    kc_distance += distance;
  }

  // average, note that we don't divide by the number of pairs
  kc_distance /= num_stabs;
  // average
  int factor = n * (n - 1) / 2;
  tmrca_mse /= num_stabs;
  tmrca_mse /= factor;
  return std::make_pair(kc_distance, tmrca_mse);
}

// Average of *squared* L2 norm over vectors, KC and TMRCA
std::tuple<arg_real_t, arg_real_t, arg_real_t>
metrics_stab_efficient(const ARG& arg1, const ARG& arg2, int num_stabs, unsigned random_kc_seed,
                       int merge_type, arg_real_t merge_fraction, bool use_r2, bool use_log2) {
  // check that extents are the same
  assert(arg1.start == arg2.start && arg1.end == arg2.end);

  // check that sample IDs are the same
  assert(arg1.leaf_ids.size() == arg2.leaf_ids.size());
  for (int leaf1 : arg1.leaf_ids) {
    assert(arg2.leaf_ids.find(leaf1) != arg2.leaf_ids.end());
  }
  // IDs must be consecutive
  for (size_t i = 0; i < arg1.leaf_ids.size(); ++i) {
    assert(arg1.leaf_ids.find(i) != arg1.leaf_ids.end());
  }

  // checking roots is relatively quick compared to checking children / parents
  try {
    arg1.check_roots();
    arg2.check_roots();
  }
  catch (...) {
    std::string my_string = "ARGs don't have correct root information.";
    my_string += " Make sure to call populate_children_and_roots().";
    throw std::logic_error(THROW_LINE(my_string));
  }

  if (use_log2 && (random_kc_seed != 0)) {
    throw std::logic_error(THROW_LINE("Cannot use both log2 and random_kc options."));
  }

  assert(merge_fraction >= 0 && merge_fraction <= 1);

  // create vectors for TMRCA, KC
  int n = arg1.leaf_ids.size();
  int num_pairs = n * (n - 1) / 2;
  vector<arg_real_t> tmrca1(num_pairs, 0);
  vector<arg_real_t> tmrca2(num_pairs, 0);
  vector<arg_real_t> mrca_root1(num_pairs, 0);
  vector<arg_real_t> mrca_root2(num_pairs, 0);
  vector<arg_real_t> split_size1(num_pairs, 0);
  vector<arg_real_t> split_size2(num_pairs, 0);

  // compute stabbing positions using the golden ratio
  arg_real_t span = arg1.end - arg1.start;
  arg_real_t phi = (sqrt(5) - 1.0) / 2.0;
  arg_real_t proportion = 0;
  vector<arg_real_t> stabs;
  for (int i = 0; i < num_stabs; ++i) {
    proportion = fmod(proportion + phi, 1.0);
    stabs.push_back(proportion * span + arg1.start);
  }

  std::mt19937 gen(random_kc_seed);
  bool random_kc = (random_kc_seed != 0);

  arg_real_t kc_distance = 0;
  arg_real_t tmrca_mse = 0;
  arg_real_t split_size_mse = 0;
  cout << utils::current_time_string() << " Total to stab: " << stabs.size() << endl;
  for (size_t k = 0; k < stabs.size(); ++k) {
    if (k % 1000 == 0) {
      cout << utils::current_time_string() << " Metrics stab " << k << endl;
    }
    arg_real_t position = stabs[k];
    std::fill(tmrca1.begin(), tmrca1.end(), 0);
    std::fill(tmrca2.begin(), tmrca2.end(), 0);
    std::fill(mrca_root1.begin(), mrca_root1.end(), 0);
    std::fill(mrca_root2.begin(), mrca_root2.end(), 0);
    std::fill(split_size1.begin(), split_size1.end(), 0);
    std::fill(split_size2.begin(), split_size2.end(), 0);
    unordered_set<int> merge_ids = get_merge_ids(arg1, position, merge_type, merge_fraction, gen);
    fill_recurse(arg1.root_at(position)->node, n, position, 0, tmrca1, mrca_root1, split_size1,
                 merge_ids, use_log2, random_kc, gen);
    // Only use merge_ids for one of the ARGs, not the other!
    merge_ids.clear();
    fill_recurse(arg2.root_at(position)->node, n, position, 0, tmrca2, mrca_root2, split_size2,
                 merge_ids, use_log2, random_kc, gen);

    if (use_r2) {
      kc_distance += utils::r2(mrca_root1, mrca_root2);
    }
    else {
      kc_distance += utils::l2(mrca_root1, mrca_root2);
      // use for comparing with tsinfer paper
      // kc_distance += sqrt(utils::l2(mrca_root1, mrca_root2));
    }
    tmrca_mse += utils::l2(tmrca1, tmrca2);

    split_size_mse += utils::l2(split_size1, split_size2);
  }

  // average, note that we don't divide by the number of pairs
  kc_distance /= num_stabs;
  // average
  tmrca_mse /= num_stabs;
  tmrca_mse /= num_pairs;
  split_size_mse /= num_stabs;
  split_size_mse /= num_pairs;
  return std::tuple<arg_real_t, arg_real_t, arg_real_t>(kc_distance, tmrca_mse, split_size_mse);
}

// Average of *squared* L2 norm over KC length-aware vectors for multiple lambdas
// Mixing occurs using (1 - lambda) * topology_vector + lambda * length_vector
// No merging / breaking of polytomies is supported here, since it's not meant for tsinfer
vector<arg_real_t> kc2_length_stab_efficient(const ARG& arg1, const ARG& arg2, int num_stabs,
                                             vector<arg_real_t> lambdas) {
  // check that extents are the same
  assert(arg1.start == arg2.start && arg1.end == arg2.end);

  // check that sample IDs are the same
  assert(arg1.leaf_ids.size() == arg2.leaf_ids.size());
  for (int leaf1 : arg1.leaf_ids) {
    assert(arg2.leaf_ids.find(leaf1) != arg2.leaf_ids.end());
  }
  // IDs must be consecutive
  for (size_t i = 0; i < arg1.leaf_ids.size(); ++i) {
    assert(arg1.leaf_ids.find(i) != arg1.leaf_ids.end());
  }

  // checking roots is relatively quick compared to checking children / parents
  try {
    arg1.check_roots();
    arg2.check_roots();
  }
  catch (...) {
    std::string my_string = "ARGs don't have correct root information.";
    my_string += " Make sure to call populate_children_and_roots().";
    throw std::logic_error(THROW_LINE(my_string));
  }

  // check vector of lambdas
  for (auto lambda : lambdas) {
    assert((lambda >= 0) && (lambda <= 1));
  }

  // create vectors for TMRCA
  int n = arg1.leaf_ids.size();
  int num_pairs = n * (n - 1) / 2;
  vector<arg_real_t> tmrca1(num_pairs, 0);
  vector<arg_real_t> tmrca2(num_pairs, 0);
  vector<arg_real_t> mrca_root1(num_pairs, 0);
  vector<arg_real_t> mrca_root2(num_pairs, 0);
  vector<arg_real_t> split_size1(num_pairs, 0);
  vector<arg_real_t> split_size2(num_pairs, 0);

  // compute stabbing positions using the golden ratio
  arg_real_t span = arg1.end - arg1.start;
  arg_real_t phi = (sqrt(5) - 1.0) / 2.0;
  arg_real_t proportion = 0;
  vector<arg_real_t> stabs;
  for (int i = 0; i < num_stabs; ++i) {
    proportion = fmod(proportion + phi, 1.0);
    stabs.push_back(proportion * span + arg1.start);
  }

  bool use_log2 = false;
  unsigned random_kc_seed = 0;
  std::mt19937 gen(random_kc_seed);
  bool random_kc = (random_kc_seed != 0);

  vector<arg_real_t> kc_lambda_distances(lambdas.size(), 0);
  cout << utils::current_time_string() << " Total to stab: " << stabs.size() << endl;
  for (size_t k = 0; k < stabs.size(); ++k) {
    if (k % 1000 == 0) {
      cout << utils::current_time_string() << " Metrics stab " << k << endl;
    }
    arg_real_t position = stabs[k];
    std::fill(tmrca1.begin(), tmrca1.end(), 0);
    std::fill(tmrca2.begin(), tmrca2.end(), 0);
    std::fill(mrca_root1.begin(), mrca_root1.end(), 0);
    std::fill(mrca_root2.begin(), mrca_root2.end(), 0);
    std::fill(split_size1.begin(), split_size1.end(), 0);
    std::fill(split_size2.begin(), split_size2.end(), 0);
    unordered_set<int> merge_ids;
    // No merge IDs being used
    fill_recurse(arg1.root_at(position)->node, n, position, 0, tmrca1, mrca_root1, split_size1,
                 merge_ids, use_log2, random_kc, gen);
    fill_recurse(arg2.root_at(position)->node, n, position, 0, tmrca2, mrca_root2, split_size2,
                 merge_ids, use_log2, random_kc, gen);

    // prepare vectors for time distances
    vector<arg_real_t> kc_times1, kc_times2;
    arg_real_t time_root1 = arg1.root_at(position)->node->height;
    arg_real_t time_root2 = arg2.root_at(position)->node->height;
    for (size_t j = 0; j < tmrca1.size(); ++j) {
      kc_times1.push_back(time_root1 - tmrca1[j]);
      kc_times2.push_back(time_root2 - tmrca2[j]);
    }
    for (size_t j = 0; j < (size_t) n; ++j) {
      kc_times1.push_back(arg1.arg_nodes.at(j)->parent_edge_at(position)->parent->height);
      kc_times2.push_back(arg2.arg_nodes.at(j)->parent_edge_at(position)->parent->height);
    }

    // prepare vectors for edge distances, we need to copy mrca_root as we change the size
    vector<arg_real_t> kc_edges1 = mrca_root1;
    vector<arg_real_t> kc_edges2 = mrca_root2;
    for (size_t j = 0; j < (size_t) n; ++j) {
      kc_edges1.push_back(1);
      kc_edges2.push_back(1);
    }

    for (size_t k = 0; k < lambdas.size(); ++k) {
      arg_real_t lambda = lambdas[k];
      // create the linear combination using (1-lambda) for topology and lambda for length
      vector<arg_real_t> kc_values1(num_pairs + n, 0);
      vector<arg_real_t> kc_values2(num_pairs + n, 0);
      for (size_t j = 0; j < kc_times1.size(); ++j) {
        kc_values1[j] = (1 - lambda) * kc_edges1[j] + lambda * kc_times1[j];
        kc_values2[j] = (1 - lambda) * kc_edges2[j] + lambda * kc_times2[j];
      }
      kc_lambda_distances[k] += utils::l2(kc_values1, kc_values2);
    }
  }

  // average, note that we don't divide by the number of pairs
  for (arg_real_t& kc_lambda_distance : kc_lambda_distances) {
    kc_lambda_distance /= num_stabs;
  }
  return kc_lambda_distances;
}

// Get KC and TMRCA pairwise vectors of length (N choose 2) at a position
pair<vector<arg_real_t>, vector<arg_real_t>> kc_tmrca_vectors(const ARG& arg, arg_real_t position) {
  // IDs must be consecutive
  for (size_t i = 0; i < arg.leaf_ids.size(); ++i) {
    assert(arg.leaf_ids.find(i) != arg.leaf_ids.end());
  }

  // checking roots is relatively quick compared to checking children / parents
  try {
    arg.check_roots();
  }
  catch (...) {
    std::string my_string = "ARGs don't have correct root information.";
    my_string += " Make sure to call populate_children_and_roots().";
    throw std::logic_error(THROW_LINE(my_string));
  }

  // create vectors for TMRCA, KC
  int n = arg.leaf_ids.size();
  int num_pairs = n * (n - 1) / 2;
  vector<arg_real_t> tmrca(num_pairs, 0);
  vector<arg_real_t> mrca_root(num_pairs, 0);
  vector<arg_real_t> split_size(num_pairs, 0);

  bool use_log2 = false;
  unsigned random_kc_seed = 0;
  std::mt19937 gen(random_kc_seed);
  bool random_kc = (random_kc_seed != 0);
  unordered_set<int> merge_ids;
  // No merge IDs being used
  fill_recurse(arg.root_at(position)->node, n, position, 0, tmrca, mrca_root, split_size, merge_ids,
               use_log2, random_kc, gen);

  return std::make_pair(mrca_root, tmrca);
}

vector<int> fill_recurse(const ARGNode* node, int n, arg_real_t position, int depth,
                         vector<arg_real_t>& tmrca, vector<arg_real_t>& mrca_root,
                         vector<arg_real_t>& split_size, const unordered_set<int>& merge_ids,
                         bool use_log2, bool random_kc, std::mt19937& gen) {
  vector<vector<int>> samples_of_children;
  vector<ARGEdge*> children_edges = node->children_at(position);

  // terminal condition for leaf nodes
  if (children_edges.size() == 0) {
    assert(node->ID < n);
    return vector<int>({node->ID});
  }
  arg_real_t delta_depth;
  if (use_log2) {
    delta_depth = log2(children_edges.size());
  }
  else {
    delta_depth = 1;
  }

  if ((!random_kc) || (merge_ids.size() > 0)) {
    for (const ARGEdge* child_edge : children_edges) {
      // Recalculate delta_depth if we're merging
      if (merge_ids.size() > 0) {
        if (merge_ids.find(child_edge->child->ID) != merge_ids.end()) {
          delta_depth = 0;
        }
        else {
          delta_depth = 1;
        }
      }
      vector<int> child_samples =
          fill_recurse(child_edge->child, n, position, depth + delta_depth, tmrca, mrca_root,
                       split_size, merge_ids, use_log2, false, gen);
      samples_of_children.push_back(child_samples);
    }
    int num_leaves_total = 0;
    for (auto leaves : samples_of_children) {
      num_leaves_total += leaves.size();
    }

    // go over pairs
    for (size_t i = 0; i < samples_of_children.size() - 1; ++i) {
      for (size_t j = i + 1; j < samples_of_children.size(); ++j) {
        for (int ival : samples_of_children[i]) {
          for (int jval : samples_of_children[j]) {
            int a = std::min(ival, jval);
            int b = std::max(ival, jval);
            assert(b > a); // don't want equal
            // We need to compute the index. Derivation:
            // n*(n-1) / 2 - (n-a)*(n-a-1)/2 + (b-a-1)
            // For instance, take n = 5, a = 2, b = 4
            // Write out the pairs as
            // (0, 1), (0, 2), (0, 3), (0, 4)
            // (1, 2), (1, 3), (1, 4),
            // (2, 3), (2, 4),
            // (3, 4)
            // We want a formula that gives out 8, because the pair (2, 4) has index 8
            // in the list. First, we need to go to the row with a = 2. The number of values
            // with a >= 2 is (n-a choose 2) = (n-a)*(n-a-1)/2. So the number of values with
            // a < 2 is (n choose 2) - (n-a choose 2). Within this row, the first pair
            // is (a, a+1), so if we want to get to (a, b), we need to move forward
            // b - (a+1) = b-a-1 units. Putting this together, the index is
            // (final expression is also the same as what is used in tsinfer evaluation scripts)
            int index = b - a - 1 + (a * (2 * n - a - 1) / 2);
            split_size[index] = num_leaves_total;
            tmrca[index] = node->height;
            mrca_root[index] = depth;
          }
        }
      }
    }
  }
  else {
    int num_children = children_edges.size();

    // Build a tree of num_children in the style of R's rtree
    vector<vector<int>> binary_tree = random_binary_tree(num_children, gen);
#ifdef _DEBUG
    if (num_children > 2) {
      cout << "Tree:" << endl;
      for (size_t i = 0; i < binary_tree.size(); ++i) {
        cout << i << " " << binary_tree[i][0] << " " << binary_tree[i][1] << endl;
      }
      cout << endl;
    }
#endif // _DEBUG

    // Sort children_edges by child node IDs for greater determinism
    std::sort(
        children_edges.begin(), children_edges.end(),
        [](const ARGEdge* a, const ARGEdge* b) -> bool { return a->child->ID < b->child->ID; });

    for (int i = 0; i < num_children; ++i) {
      const ARGEdge* child_edge = children_edges[i];
      delta_depth = binary_tree[i][0];
      vector<int> child_samples =
          fill_recurse(child_edge->child, n, position, depth + delta_depth, tmrca, mrca_root,
                       split_size, merge_ids, use_log2, random_kc, gen);
      samples_of_children.push_back(child_samples);
    }
    int num_leaves_total = 0;
    for (auto leaves : samples_of_children) {
      num_leaves_total += leaves.size();
    }

    // go over pairs
    for (size_t i = 0; i < samples_of_children.size() - 1; ++i) {
      for (size_t j = i + 1; j < samples_of_children.size(); ++j) {
        int small = std::min(i, j);
        int big = std::max(i, j);
        while (small != big) {
          int parent = binary_tree[small][1];
          int temp = big;
          small = std::min(temp, parent);
          big = std::max(temp, parent);
        }
        delta_depth = binary_tree[small][0];

        for (int ival : samples_of_children[i]) {
          for (int jval : samples_of_children[j]) {
            int a = std::min(ival, jval);
            int b = std::max(ival, jval);
            assert(b > a); // don't want equal
            // see explanation above for this expression
            int index = b - a - 1 + (a * (2 * n - a - 1) / 2);
            split_size[index] = num_leaves_total;
            tmrca[index] = node->height;
            mrca_root[index] = depth + delta_depth;
          }
        }
      }
    }
  }

  // combine samples_of_children into all samples below this node and return
  vector<int> current = samples_of_children[0];
  vector<int> merged;
  for (size_t i = 1; i < samples_of_children.size(); ++i) {
    merged.clear();
    std::merge(samples_of_children[i].begin(), samples_of_children[i].end(), current.begin(),
               current.end(), std::back_inserter(merged));
    current = merged;
  }
  return merged;
}

// Top-down splitting algorithm used in R's rtree, which is what tsinfer authors use
vector<vector<int>> random_binary_tree(int k, std::mt19937& gen) {
  // first entry is depth, second entry is parent ID
  vector<vector<int>> tree(2 * k - 1, vector<int>(2, 0));

  stack<vector<int>> to_process;
  int internal_id = 2 * k - 2;
  int leaf_id = 0;
  vector<int> start_element = {k, 0, internal_id};
  to_process.push(start_element);

  while (!to_process.empty()) {
    vector<int> element = to_process.top();
    to_process.pop();
    if (element[0] == 1) {
      tree[leaf_id][0] = element[1];
      tree[leaf_id][1] = element[2];
      ++leaf_id;
    }
    else {
      tree[internal_id][0] = element[1];
      tree[internal_id][1] = element[2];

      std::uniform_int_distribution<int> distribution(1, element[0] - 1);
      int split = distribution(gen);
      vector<int> left = {split, element[1] + 1, internal_id};
      vector<int> right = {element[0] - split, element[1] + 1, internal_id};
      to_process.push(right);
      to_process.push(left);
      --internal_id;
    }
  }

  // shuffle the first n elements using the same RNG
  auto it = tree.begin();
  for (int i = 0; i < k; ++i) {
    ++it;
  }
  std::shuffle(tree.begin(), it, gen);

  return tree;
}

/*
merge_type = 0: no merging
merge_type = 1: random merging
merge_type = 2: branch height / parent height
*/
unordered_set<int> get_merge_ids(const ARG& arg, arg_real_t position, int merge_type,
                                 arg_real_t merge_fraction, std::mt19937& gen) {
  if (merge_type == 0) {
    return unordered_set<int>();
  }
  vector<pair<arg_real_t, int>> node_and_fraction;
  stack<const ARGNode*> dfs;
  dfs.push(arg.root_at(position)->node);
  while (!dfs.empty()) {
    const ARGNode* current = dfs.top();
    dfs.pop();

    vector<ARGEdge*> children_edges = current->children_at(position);
    for (const ARGEdge* child_edge : children_edges) {
      const ARGNode* child = child_edge->child;
      // we don't include any branches to leaf nodes
      if (child->height > 0) {
        dfs.push(child);
        arg_real_t fraction = 1 - child->height / current->height;
        node_and_fraction.push_back(std::make_pair(fraction, child->ID));
      }
    }
  }

  if (merge_type == 1) {
    std::shuffle(node_and_fraction.begin(), node_and_fraction.end(), gen);
  }
  else if (merge_type == 2) {
    std::sort(node_and_fraction.begin(), node_and_fraction.end());
  }
  else {
    return unordered_set<int>();
  }

  unordered_set<int> result;
  int k = (int) (node_and_fraction.size() * merge_fraction);
  for (size_t i = 0; i < (unsigned) k; ++i) {
    result.insert(node_and_fraction[i].second);
  }
  return result;
}

// Bitset overlap for ARG-wide branch precision and recall
// positional arguments are primarily used for testing
// Returns:
//     Number of bitsets in first ARG, across whole ARG
//     Number of bitsets in second ARG, across whole ARG
//     Number of bitsets in common, across whole ARG
tuple<int, int, int> bitset_overlap_full(const ARG& arg1, const ARG& arg2, arg_real_t min_position,
                                         arg_real_t max_position) {

  // check that extents are the same
  assert(arg1.start == arg2.start && arg1.end == arg2.end);

  // check that sample IDs are the same
  assert(arg1.leaf_ids.size() == arg2.leaf_ids.size());
  for (int leaf1 : arg1.leaf_ids) {
    assert(arg2.leaf_ids.find(leaf1) != arg2.leaf_ids.end());
  }
  // IDs must be consecutive
  for (size_t i = 0; i < arg1.leaf_ids.size(); ++i) {
    assert(arg1.leaf_ids.find(i) != arg1.leaf_ids.end());
  }

  // checking roots is relatively quick compared to checking children / parents
  try {
    arg1.check_roots();
    arg2.check_roots();
  }
  catch (...) {
    std::string my_string = "ARGs don't have correct root information.";
    my_string += " Make sure to call populate_children_and_roots().";
    throw std::logic_error(THROW_LINE(my_string));
  }

  size_t n = arg1.leaf_ids.size();
  unordered_set<DescendantList, DescendantListHash> bitsets1, bitsets2;

  // because we get out the values, which is non-const, we can't use
  // const DescendantList&
  arg_utils::visit_clades(
      arg1,
      [&bitsets1, n](
          DescendantList& desc_list, const ARGNode* node, arg_real_t start, arg_real_t end) {
        (void) node;
        (void) start;
        (void) end;
        // don't want roots
        if (desc_list.num_values() != n) {
          bitsets1.insert(desc_list);
        }
      },
      min_position, max_position);

  // because we get out the values, which is non-const, we can't use
  // const DescendantList&
  arg_utils::visit_clades(
      arg2,
      [&bitsets2, n](
          DescendantList& desc_list, const ARGNode* node, arg_real_t start, arg_real_t end) {
        (void) node;
        (void) start;
        (void) end;
        // don't want roots
        if (desc_list.num_values() != n) {
          bitsets2.insert(desc_list);
        }
      },
      min_position, max_position);

  int num1 = bitsets1.size();
  int num2 = bitsets2.size();
  int num_common = 0;

  for (const DescendantList& desc_list : bitsets1) {
    if (bitsets2.find(desc_list) != bitsets2.end()) {
      ++num_common;
    }
  }

  return std::make_tuple(num1, num2, num_common);
}

// Bitset overlap for branch and mutation precision and recall, sampled based on stabbing queries
// Returns:
//     Number of bitsets in first ARG across num_stabs positions
//     Number of bitsets in second ARG across num_stabs positions
//     Number of bitsets in common across num_stabs positions
//     Branch length of bitsets in first ARG across num_stabs positions
//     Branch length of bitsets in second ARG across num_stabs positions
//     Branch length of common bitsets across num_stabs positions
//
// Arguments:
//     arg2_factor: a scaling factor multiplied on to branch lengths of the second ARG
//     random_resolve_seed: if nonzero, then polytomies are randomly resolved for topology metric
//     min_mac: min MAC considered (inclusive, 0 means no filter)
//     max_mac: max MAC considered (exclusive, 0 means no filter)
tuple<int, int, int, arg_real_t, arg_real_t, arg_real_t>
bitset_overlap_stab(const ARG& arg1, const ARG& arg2, int num_stabs, arg_real_t arg2_factor,
                    unsigned random_resolve_seed, int min_mac, int max_mac) {

  // check that extents are the same
  assert(arg1.start == arg2.start && arg1.end == arg2.end);

  // check that sample IDs are the same
  assert(arg1.leaf_ids.size() == arg2.leaf_ids.size());
  for (int leaf1 : arg1.leaf_ids) {
    assert(arg2.leaf_ids.find(leaf1) != arg2.leaf_ids.end());
  }
  // IDs must be consecutive
  for (size_t i = 0; i < arg1.leaf_ids.size(); ++i) {
    assert(arg1.leaf_ids.find(i) != arg1.leaf_ids.end());
  }

  // checking roots is relatively quick compared to checking children / parents
  try {
    arg1.check_roots();
    arg2.check_roots();
  }
  catch (...) {
    std::string my_string = "ARGs don't have correct root information.";
    my_string += " Make sure to call populate_children_and_roots().";
    throw std::logic_error(THROW_LINE(my_string));
  }

  int n = arg1.leaf_ids.size();
  if (max_mac == 0) {
    max_mac = n;
  }

  std::mt19937 gen(random_resolve_seed);
  bool random_resolve = (random_resolve_seed != 0);

  // compute stabbing positions using the golden ratio
  arg_real_t span = arg1.end - arg1.start;
  arg_real_t phi = (sqrt(5) - 1.0) / 2.0;
  arg_real_t proportion = 0;
  vector<arg_real_t> stabs;
  for (int i = 0; i < num_stabs; ++i) {
    proportion = fmod(proportion + phi, 1.0);
    stabs.push_back(proportion * span + arg1.start);
  }

  int num1 = 0, num2 = 0, num_common = 0;
  arg_real_t length1 = 0, length2 = 0, length_common = 0;

  cout << utils::current_time_string() << " Total to stab: " << stabs.size() << endl;
  for (size_t k = 0; k < stabs.size(); ++k) {
    if (k % 1000 == 0) {
      cout << utils::current_time_string() << " Metrics stab " << k << endl;
    }
    arg_real_t position = stabs[k];
    unordered_map<DescendantList, pair<arg_real_t, arg_real_t>, DescendantListHash> bitset_lengths;
    fill_bitsets_recurse(
        bitset_lengths, arg1.root_at(position)->node, n, position, 0, random_resolve, gen);
    fill_bitsets_recurse(
        bitset_lengths, arg2.root_at(position)->node, n, position, 1, random_resolve, gen);

#ifdef _DEBUG
    for (const auto& element : bitset_lengths) {
      cout << element.first.to_bitset_string() << " ";
      cout << element.second.first << " ";
      cout << element.second.second << endl;
    }
    cout << endl;
#endif // _DEBUG

    for (const auto& element : bitset_lengths) {
      // check for MAC filter
      int mac = std::min((int) element.first.num_values(), n - (int) element.first.num_values());
      if ((mac < min_mac) || (mac >= max_mac)) {
        continue;
      }
      const auto& lengths = element.second;
      // explanation: a value of -1 means bitset not found, a value of 0 means bitset
      // found when randomly resolving, and a value > 0 means bitset found without needing
      // to randomly resolve
      if (lengths.first >= 0) {
        num1 += 1;
        length1 += lengths.first;
      }
      if (lengths.second >= 0) {
        num2 += 1;
        length2 += lengths.second * arg2_factor;
      }
      if ((lengths.first >= 0) && (lengths.second >= 0)) {
        num_common += 1;
      }
      arg_real_t min_length = std::min(lengths.first, lengths.second * arg2_factor);
      if (min_length > 0) { // needed to avoid adding -1
        length_common += min_length;
      }
    }
  }

  return std::make_tuple(num1, num2, num_common, length1, length2, length_common);
}

DescendantList fill_bitsets_recurse(
    unordered_map<DescendantList, pair<arg_real_t, arg_real_t>, DescendantListHash>& bitset_lengths,
    const ARGNode* node, int n, arg_real_t position, int index, bool random_resolve,
    std::mt19937& gen) {

  vector<ARGEdge*> children_edges = node->children_at(position);

  // terminal condition for leaf nodes
  if (children_edges.size() == 0) {
    assert(node->ID < n);
    return DescendantList(n, node->ID);
  }

  DescendantList desc_list(n);
  if ((!random_resolve) || children_edges.size() == 2) {
    for (const ARGEdge* child_edge : children_edges) {
      const ARGNode* child = child_edge->child;
      DescendantList child_desc_list =
          fill_bitsets_recurse(bitset_lengths, child, n, position, index, random_resolve, gen);

      arg_real_t height_diff = node->height - child->height;
      bool not_found = (bitset_lengths.find(child_desc_list) == bitset_lengths.end());
      pair<arg_real_t, arg_real_t>& entry =
          bitset_lengths[child_desc_list]; // creates if not present, only hashes once
      // the above will use value-initialization if not present, so both elements
      // of the pair would be 0
      if (index == 0) {
        entry.first = height_diff;
        if (not_found) {
          entry.second = -1;
        }
      }
      else {
        entry.second = height_diff;
        if (not_found) {
          entry.first = -1;
        }
      }

      desc_list.add(child_desc_list);
    }
  }
  else {
    vector<vector<int>> binary_tree = random_binary_tree(children_edges.size(), gen);

    // Sort children_edges by child node IDs for greater determinism
    std::sort(
        children_edges.begin(), children_edges.end(),
        [](const ARGEdge* a, const ARGEdge* b) -> bool { return a->child->ID < b->child->ID; });

    vector<DescendantList> resolved_descendant_lists;
    for (size_t i = 0; i < binary_tree.size(); ++i) {
      resolved_descendant_lists.emplace_back(n); // pushes an empty DescendantList of size n
    }

    for (size_t i = 0; i < children_edges.size(); ++i) {
      const ARGEdge* child_edge = children_edges[i];
      const ARGNode* child = child_edge->child;
      DescendantList child_desc_list =
          fill_bitsets_recurse(bitset_lengths, child, n, position, index, random_resolve, gen);

      arg_real_t height_diff = node->height - child->height;
      bool not_found = (bitset_lengths.find(child_desc_list) == bitset_lengths.end());
      pair<arg_real_t, arg_real_t>& entry =
          bitset_lengths[child_desc_list]; // creates if not present, only hashes once
      // the above will use value-initialization if not present, so both elements
      // of the pair would be 0
      if (index == 0) {
        entry.first = height_diff;
        if (not_found) {
          entry.second = -1;
        }
      }
      else {
        entry.second = height_diff;
        if (not_found) {
          entry.first = -1;
        }
      }

      // resolved_descendant_lists[i].add(child_desc_list); // not necessary
      resolved_descendant_lists[binary_tree[i][1]].add(child_desc_list);
    }
    for (size_t i = children_edges.size(); i < binary_tree.size(); ++i) {
      if (i == binary_tree.size() - 1) {
        desc_list.add(resolved_descendant_lists[i]);
      }
      else {
        DescendantList child_desc_list = resolved_descendant_lists[i];
        bool not_found = (bitset_lengths.find(child_desc_list) == bitset_lengths.end());
        pair<arg_real_t, arg_real_t>& entry =
            bitset_lengths[child_desc_list]; // creates if not present, only hashes once
        // the above will use value-initialization if not present, so both elements
        // of the pair would be 0
        if (index == 0) {
          entry.first = 0;
          if (not_found) {
            entry.second = -1;
          }
        }
        else {
          entry.second = 0;
          if (not_found) {
            entry.first = -1;
          }
        }
        resolved_descendant_lists[binary_tree[i][1]].add(child_desc_list);
      }
    }
  }

  return desc_list;
}

vector<tuple<int, arg_real_t, DescendantList>> stab_return_all_bitsets(const ARG& arg,
                                                                       arg_real_t position) {
  // check valid position
  assert(position >= arg.start && position < arg.end);

  // IDs must be consecutive
  for (size_t i = 0; i < arg.leaf_ids.size(); ++i) {
    assert(arg.leaf_ids.find(i) != arg.leaf_ids.end());
  }

  // checking roots is relatively quick compared to checking children / parents
  try {
    arg.check_roots();
  }
  catch (...) {
    std::string my_string = "ARG doesn't have correct root information.";
    my_string += " Make sure to call populate_children_and_roots().";
    throw std::logic_error(THROW_LINE(my_string));
  }

  int n = arg.leaf_ids.size();
  std::mt19937 gen(0); // unused, a dummy value
  bool random_resolve = false;

  // this data structure allows for two sets of length values, for computing overlap
  // we only need one set of values, so only need the first part of the pairs
  unordered_map<DescendantList, pair<arg_real_t, arg_real_t>, DescendantListHash> bitset_lengths;
  fill_bitsets_recurse(
      bitset_lengths, arg.root_at(position)->node, n, position, 0, random_resolve, gen);

  vector<tuple<int, arg_real_t, DescendantList>> result;
  for (const auto& element : bitset_lengths) {
    int allele_count = element.first.num_values();
    arg_real_t length = element.second.first;
    result.emplace_back(allele_count, length, element.first);
  }
  std::sort(result.begin(), result.end(),
            [](tuple<int, arg_real_t, DescendantList> a,
               tuple<int, arg_real_t, DescendantList> b) -> bool {
              if (std::get<0>(a) == std::get<0>(b)) {
                return std::get<1>(a) < std::get<1>(b);
              }
              else {
                return std::get<0>(a) < std::get<0>(b);
              }
            }); // increasing order by allele count, then by bitset length

  return result;
}

// Return the DescendantList subtending a node at a given position
// This currently uses recursion but can change to postorder using stacks
DescendantList get_bitset(const ARGNode* node, int n, arg_real_t position) {
  // Get the children of the particular node
  vector<ARGEdge*> children_edges = node->children_at(position);
  if (children_edges.size() == 0) {
    assert(node->ID < n);
    return DescendantList(n, node->ID);
  }
  DescendantList desc_list(n);
  for (const ARGEdge* child_edge : children_edges) {
    const ARGNode* child = child_edge->child;

    DescendantList child_desc_list = get_bitset(child, n, position);
    // Add in the descendent list of the child (recursive traversal)
    desc_list.add(child_desc_list);
  }
  return (desc_list);
}

// Return the DescendantList subtending a mutation
DescendantList get_carriers(const ARG& arg, const Mutation* mutation) {
  if (arg.roots.empty()) {
    throw std::logic_error(THROW_LINE("Call populate_children_and_roots() first."));
  }
  arg_real_t pos = mutation->position;
  ARGNode* child_node = mutation->edge->child;
  auto n = arg.leaf_ids.size();
  DescendantList desc_list = get_bitset(child_node, n, pos);
  return desc_list;
}

// Return an integer genotype vector subtending a mutation
vector<int> get_mutation_genotype(const ARG& arg, const Mutation* mutation, bool diploid) {
  DescendantList desc_list = get_carriers(arg, mutation);
  auto n = arg.leaf_ids.size();
  vector<int> genotype(diploid ? n / 2 : n, 0);
  for (auto id : desc_list.values()) {
    genotype.at(diploid ? id / 2 : id) += 1;
  }
  return genotype;
}

vector<arg_real_t> impute(const ARG& arg, arg_real_t position, const vector<int>& genotypes,
                          bool old) {
  int n = arg.leaf_ids.size();
  assert(n == (int) genotypes.size());
  int num_zeros = 0;
  int num_ones = 0;
  for (size_t i = 0; i < genotypes.size(); ++i) {
    if (genotypes[i] == 0) {
      num_zeros += 1;
    }
    else if (genotypes[i] == 1) {
      num_ones += 1;
    }
    else if (genotypes[i] != -1) {
      throw std::logic_error(THROW_LINE("Bad input: all values must be 0, 1, or -1."));
    }
  }
  if (num_zeros == 0 && num_ones == 0) {
    throw std::logic_error(THROW_LINE("Nothing to use to impute from."));
  }

  vector<arg_real_t> result(n, 0);
  for (size_t i = 0; i < genotypes.size(); ++i) {
    result[i] = genotypes[i];
  }
  if (num_ones == 0) {
    // impute to all zeros. Would this ever happen?
    for (auto& v : result) {
      v = 0;
    }
    return result;
  }
  if (num_zeros == 0) {
    // impute to all ones. Would this ever happen?
    for (auto& v : result) {
      v = 1;
    }
    return result;
  }

  const ARGNode* root_node = arg.root_at(position)->node;
  stack<pair<int, int>> dfs;
  dfs.push(std::make_pair(root_node->ID, -1)); // first is node ID, second is parent ID
  unordered_map<int, int> zeros_below, ones_below;
  int best_diff = n; // start at a bad value, can only get better
  arg_real_t best_height = 0;
  int best_ancestral = 0;
  int best_branch = 0; // ID of child node of this branch
  vector<arg_real_t> best_heights;
  vector<int> best_ancestrals;
  vector<int> best_branches;

  while (!dfs.empty()) {
    pair<int, int> element = dfs.top();
    int node_id = element.first;
    if (zeros_below.find(node_id) == zeros_below.end()) {
      zeros_below[node_id] = 0;
      ones_below[node_id] = 0;
      const vector<ARGEdge*> children_edges = arg.arg_nodes.at(node_id)->children_at(position);
      // condition for leaf nodes
      if (children_edges.size() == 0) {
        if (genotypes[node_id] == 0) {
          zeros_below[node_id] = 1;
        }
        else if (genotypes[node_id] == 1) {
          ones_below[node_id] = 1;
        }
      }
      else {
        for (const ARGEdge* child_edge : children_edges) {
          dfs.push(std::make_pair(child_edge->child->ID, node_id));
        }
      }
    }
    else {
      dfs.pop();
      // we don't need to consider the root node for finding the best branch
      if (element.second != -1) {
        // add values to the parents
        zeros_below[element.second] += zeros_below[element.first];
        ones_below[element.second] += ones_below[element.first];

        int node_id = element.first;
        // we don't need to consider clades that only have -1 underneath
        if ((zeros_below[node_id] > 0) || (ones_below[node_id] > 0)) {
          int same = ones_below[node_id] + (num_zeros - zeros_below[node_id]);
          int diff = num_zeros + num_ones - same;
          arg_real_t branch_height =
              arg.arg_nodes.at(element.second)->height - arg.arg_nodes.at(node_id)->height;
          if (old) {
            if ((diff < best_diff) || ((diff == best_diff) && (branch_height > best_height))) {
              best_diff = diff;
              best_height = branch_height;
              best_branch = node_id;
              best_ancestral = 0; // in this case, ancestral is 0 and derived is 1
            }
            if ((same < best_diff) || ((same == best_diff) && (branch_height > best_height))) {
              best_diff = same;
              best_height = branch_height;
              best_branch = node_id;
              best_ancestral = 1; // in this case, ancestral is 1 and derived is 0
            }
            // note that we want to do both checks, not an if / else-if
          }
          else {
            if (diff <= best_diff) {
              if (diff < best_diff) {
                best_diff = diff;
                best_heights.clear();
                best_ancestrals.clear();
                best_branches.clear();
              }
              best_heights.push_back(branch_height);
              best_branches.push_back(node_id);
              best_ancestrals.push_back(0);
            }
            if (same <= best_diff) {
              if (same < best_diff) {
                best_diff = same;
                best_heights.clear();
                best_ancestrals.clear();
                best_branches.clear();
              }
              best_heights.push_back(branch_height);
              best_branches.push_back(node_id);
              best_ancestrals.push_back(1);
            }
            // note that we want to do both checks, not an if / else-if
          }
        }
      }
    }
  }

  if (!old) {
    arg_real_t height_sum = 0;
    for (auto& v : best_heights) {
      height_sum += v;
    }

    vector<arg_real_t> height_fractions;
    for (auto& v : best_heights) {
      height_fractions.push_back(v / height_sum);
    }

    arg_real_t root_value = 0;
    for (size_t i = 0; i < height_fractions.size(); ++i) {
      root_value += best_ancestrals[i] * height_fractions[i];
    }

    stack<pair<int, arg_real_t>> new_dfs;
    new_dfs.push(std::make_pair(root_node->ID, root_value)); // first is node ID, second is dosage

    while (!new_dfs.empty()) {
      pair<int, arg_real_t> element = new_dfs.top();
      new_dfs.pop();
      int node_id = element.first;
      arg_real_t next_value = element.second;
      for (size_t i = 0; i < best_branches.size(); ++i) {
        if (node_id == best_branches[i]) {
          if (best_ancestrals[i] == 0) {
            next_value += height_fractions[i];
          }
          else {
            next_value -= height_fractions[i];
          }
        }
      }
      vector<ARGEdge*> children_edges = arg.arg_nodes.at(node_id)->children_at(position);
      // condition for leaf nodes
      if (children_edges.size() == 0) {
        // only write if the current genotype is -1
        if (result[node_id] == -1) {
          result[node_id] = next_value;
        }
      }
      else {
        for (const ARGEdge* child_edge : children_edges) {
          new_dfs.push(std::make_pair(child_edge->child->ID, next_value));
        }
      }
    }

    for (auto& v : result) {
      if (v < 0) {
        v = 0;
      }
      if (v > 1) {
        v = 1;
      }
    }

    return result;
  }

  // now that we have the best branch information, go through with a DFS and impute
  dfs.push(std::make_pair(root_node->ID, 0)); // first is node ID, second is whether below mutation
  while (!dfs.empty()) {
    pair<int, int> element = dfs.top();
    dfs.pop();
    int node_id = element.first;
    int next_status = (best_branch == node_id) ? 1 : element.second;

    vector<ARGEdge*> children_edges = arg.arg_nodes.at(node_id)->children_at(position);
    // condition for leaf nodes
    if (children_edges.size() == 0) {
      // only write if the current genotype is -1
      if (result[node_id] == -1) {
        if (next_status == 0) {
          result[node_id] = best_ancestral;
        }
        else {
          result[node_id] = 1 - best_ancestral;
        }
      }
    }
    else {
      for (const ARGEdge* child_edge : children_edges) {
        dfs.push(std::make_pair(child_edge->child->ID, next_status));
      }
    }
  }

  return result;
}

// Mutation match at a site
bool mutation_match(const ARG& arg, arg_real_t position, const vector<int>& genotypes) {
  size_t n = arg.leaf_ids.size();
  assert(n == genotypes.size());
  size_t num_ones = 0;
  for (size_t i = 0; i < genotypes.size(); ++i) {
    num_ones += genotypes[i];
  }
  if ((num_ones == 0) || (num_ones == n)) {
    return true;
  }

  const ARGNode* root_node = arg.root_at(position)->node;
  pair<bool, vector<int>> result =
      mutation_match_recurse(root_node, position, n - num_ones, num_ones, genotypes);
  return result.first;
}

pair<bool, vector<int>> mutation_match_recurse(const ARGNode* node, arg_real_t position,
                                               size_t num_zeros, size_t num_ones,
                                               const vector<int>& genotypes) {

  vector<vector<int>> samples_of_children;
  vector<ARGEdge*> children_edges = node->children_at(position);

  // terminal condition for leaf nodes
  if (children_edges.size() == 0) {
    return std::make_pair(false, vector<int>({node->ID}));
  }
  // consider the return values from the child edges underneath
  for (const ARGEdge* child_edge : children_edges) {
    pair<bool, vector<int>> child_result =
        mutation_match_recurse(child_edge->child, position, num_zeros, num_ones, genotypes);
    if (child_result.first) {
      return std::make_pair(true, vector<int>());
    }

    if (child_result.second.size() == num_zeros) {
      bool matches = true;
      for (auto v : child_result.second) {
        if (genotypes[v] != 0) {
          matches = false;
          break;
        }
      }
      if (matches) {
        return std::make_pair(true, vector<int>());
      }
    }
    if (child_result.second.size() == num_ones) {
      bool matches = true;
      for (auto v : child_result.second) {
        if (genotypes[v] != 1) {
          matches = false;
          break;
        }
      }
      if (matches) {
        return std::make_pair(true, vector<int>());
      }
    }

    samples_of_children.push_back(child_result.second);
  }

  // combine samples_of_children into all samples below this node and return
  vector<int> current = samples_of_children[0];
  vector<int> merged;
  for (size_t i = 1; i < samples_of_children.size(); ++i) {
    merged.clear();
    std::merge(samples_of_children[i].begin(), samples_of_children[i].end(), current.begin(),
               current.end(), std::back_inserter(merged));
    current = merged;
  }
  return std::make_pair(false, merged);
}

// Find a branch with the best Hamming distance to this mutation
int mutation_best(const ARG& arg, arg_real_t position, const vector<int>& genotypes,
                  unsigned random_seed) {
  int n = arg.leaf_ids.size();
  assert(n == (int) genotypes.size());
  int num_ones = 0;
  for (size_t i = 0; i < genotypes.size(); ++i) {
    num_ones += genotypes[i];
  }
  if ((num_ones <= 1) || (num_ones >= n - 1)) {
    return 0; // can always place singletons or null
  }

  std::mt19937 gen(random_seed);
  bool random = (random_seed != 0);

  const ARGNode* root_node = arg.root_at(position)->node;
  tuple<int, int, int> result =
      mutation_best_recurse(root_node, position, n - num_ones, genotypes, random, gen);
  assert(std::get<1>(result) + std::get<2>(result) == n);
  return std::get<0>(result);
}

// best Hamming distance, num_zeros, num_ones (under this branch)
tuple<int, int, int> mutation_best_recurse(const ARGNode* node, arg_real_t position,
                                           int genotypes_num_zeros, const vector<int>& genotypes,
                                           bool random, std::mt19937& gen) {

  if (!random) {
    int num_zeros = 0;
    int num_ones = 0;
    int best_diff = genotypes.size();
    vector<ARGEdge*> children_edges = node->children_at(position);
    // condition for leaf nodes
    if (children_edges.size() == 0) {
      if (genotypes[node->ID] == 0) {
        num_zeros += 1;
      }
      else {
        num_ones += 1;
      }
    }

    // consider children
    for (const ARGEdge* child_edge : children_edges) {
      tuple<int, int, int> child_result = mutation_best_recurse(
          child_edge->child, position, genotypes_num_zeros, genotypes, random, gen);
      int under_diff = std::get<0>(child_result);
      if (under_diff < best_diff) {
        best_diff = under_diff;
      }
      num_zeros += std::get<1>(child_result);
      num_ones += std::get<2>(child_result);
    }

    int n = genotypes.size();
    int same = num_ones + (genotypes_num_zeros - num_zeros);
    int diff = n - same;
    if (diff < best_diff) {
      best_diff = diff;
    }
    if (same < best_diff) {
      best_diff = same;
    }

    return std::make_tuple(best_diff, num_zeros, num_ones);
  }
  else {
    vector<ARGEdge*> children_edges = node->children_at(position);
    vector<vector<int>> zeros_and_ones;
    vector<vector<int>> binary_tree;
    int num_children = children_edges.size();
    int best_diff = genotypes.size();
    if (num_children == 0) {
      zeros_and_ones = vector<vector<int>>(1, vector<int>(2, 0));
      if (genotypes[node->ID] == 0) {
        zeros_and_ones[0][0] = 1;
      }
      else {
        zeros_and_ones[0][1] = 1;
      }
      binary_tree = {{0, -1}};
    }
    else {
      binary_tree = random_binary_tree(num_children, gen);
      // Sort children_edges by child node IDs for greater determinism
      std::sort(
          children_edges.begin(), children_edges.end(),
          [](const ARGEdge* a, const ARGEdge* b) -> bool { return a->child->ID < b->child->ID; });
      zeros_and_ones = vector<vector<int>>(binary_tree.size(), vector<int>(2, 0));
      for (int i = 0; i < num_children; ++i) {
        const ARGEdge* child_edge = children_edges[i];
        tuple<int, int, int> child_result = mutation_best_recurse(
            child_edge->child, position, genotypes_num_zeros, genotypes, random, gen);
        int child_diff = std::get<0>(child_result);
        if (child_diff < best_diff) {
          best_diff = child_diff;
        }
        int parent_index = binary_tree[i][1];
        if (parent_index != -1) {
          // shouldn't be necessary unless num_children = 1?
          zeros_and_ones[parent_index][0] += std::get<1>(child_result);
          zeros_and_ones[parent_index][1] += std::get<2>(child_result);
        }
      }
    }
    int num_zeros = 0;
    int num_ones = 0;
    for (int i = num_children; i < (int) binary_tree.size(); ++i) {
      num_zeros = zeros_and_ones[i][0];
      num_ones = zeros_and_ones[i][1];
      int parent_index = binary_tree[i][1];
      if (parent_index != -1) {
        zeros_and_ones[parent_index][0] += num_zeros;
        zeros_and_ones[parent_index][1] += num_ones;
      }

      int n = genotypes.size();
      int same = num_ones + (genotypes_num_zeros - num_zeros);
      int diff = n - same;
      if (diff < best_diff) {
        best_diff = diff;
      }
      if (same < best_diff) {
        best_diff = same;
      }
    }
    return std::make_tuple(best_diff, num_zeros, num_ones);
  }
}

Eigen::MatrixXd compute_grm(const ARG& arg, arg_real_t alpha, int batch_size, bool diploid,
                            double min_maf, double max_maf) {
  // need to have mutations to get GRM
  if (arg.get_mutations().empty()) {
    throw std::logic_error(
        THROW_LINE("Need pre-existing mutations to run! Try generate_mutations_and_keep."));
  }
  // Sample IDs must be consecutive
  for (size_t i = 0; i < arg.leaf_ids.size(); ++i) {
    assert(arg.leaf_ids.find(i) != arg.leaf_ids.end());
  }
  // if diploid check even
  if (diploid) {
    if (arg.leaf_ids.size() % 2 == 1) {
      throw std::logic_error(THROW_LINE("Haploid samples should be even for diploid=true."));
    }
  }
  // checking roots is relatively quick compared to checking children / parents
  try {
    arg.check_roots();
  }
  catch (...) {
    std::string my_string = "ARG doesn't have correct root information.";
    my_string += " Make sure to call populate_children_and_roots().";
    throw std::logic_error(THROW_LINE(my_string));
  }

  unsigned int n = arg.leaf_ids.size();
  unsigned int n_dip = n / 2;
  Eigen::MatrixXd GRM(diploid ? n_dip : n, diploid ? n_dip : n);
  GRM.setZero();
  Eigen::MatrixXd batch(batch_size, diploid ? n_dip : n);
  batch.setZero();
  unsigned int index = 0;
  unsigned int batch_index = 0;
  unsigned int num_mutations = arg.get_mutations().size();

  // Visit mutations and compute GRM entries
  arg_utils::visit_mutations(
      arg, [&GRM, &batch, &index, &batch_index, min_maf, max_maf, alpha, n, n_dip, batch_size,
            diploid, num_mutations](DescendantList& desc_list, const Mutation* mutation) {
        unsigned int allele_count = desc_list.values().size();
        unsigned int mac = std::min(allele_count, n - allele_count);
        double maf = (mac * 1.0) / n;
        if ((maf != 0) && (maf > min_maf) && (maf <= max_maf)) {
          for (int i : desc_list.values()) {
            unsigned int col = diploid ? (i / 2) : i;
            batch(batch_index, col) += 1;
          }
          batch_index += 1;
        }
        if (batch_index == batch_size || index == num_mutations - 1) {
          Eigen::MatrixXd batch_t = batch.transpose();
          Eigen::RowVectorXd mean = batch_t.colwise().mean();
          Eigen::RowVectorXd std =
              ((batch_t.rowwise() - mean).array().square().colwise().sum() / (batch_t.rows() - 1))
                  .sqrt();
          // monomorphic variants (e.g. rows of 0 in last batch) have mean =1 or =0 and std=0.
          // replace with 1 to avoid NaN (could instead filter MAC>0 in visit_mutations)
          std = (std.array() == 0).select(1, std);
          std = std.array().pow(-alpha);
          Eigen::MatrixXd result = (batch_t.rowwise() - mean).array().rowwise() / std.array();
          GRM.noalias() = GRM + result * result.transpose();
          batch.setZero();
          batch_index = 0;
        }
        index += 1;
      });
  return GRM;
}

vector<vector<arg_real_t>> distance_matrix(const ARG& arg) {
  // Sample IDs must be consecutive
  for (size_t i = 0; i < arg.leaf_ids.size(); ++i) {
    assert(arg.leaf_ids.find(i) != arg.leaf_ids.end());
  }
  // checking roots is relatively quick compared to checking children / parents
  try {
    arg.check_roots();
  }
  catch (...) {
    std::string my_string = "ARG doesn't have correct root information.";
    my_string += " Make sure to call populate_children_and_roots().";
    throw std::logic_error(THROW_LINE(my_string));
  }

  int n = arg.leaf_ids.size();
  vector<vector<arg_real_t>> upper_diagonal;
  for (int i = 0; i < n - 1; ++i) {
    upper_diagonal.push_back(vector<arg_real_t>(n - 1 - i, 0));
  }

  // because we get out the values, which is non-const, we can't use
  // const DescendantList&
  arg_utils::visit_branches(
      arg, [&upper_diagonal](DescendantList& desc_list, DescendantList& parent_desc_list,
                             const ARGNode* parent, const ARGNode* child, arg_real_t start,
                             arg_real_t end) {
        (void) child;
        arg_real_t volume = parent->height * (end - start);
        for (int v1 : desc_list.values()) {
          for (int v2 : parent_desc_list.values()) {
            assert(v1 != v2);
            int i = std::min(v1, v2);
            int j = std::max(v1, v2);
            upper_diagonal[i][j - i - 1] += volume;
          }
        }
      });

#ifdef _DEBUG
  for (int i = 0; i < n - 1; ++i) {
    for (auto j : upper_diagonal[i]) {
      cout << j << " ";
    }
    cout << endl;
  }
#endif // _DEBUG

  return upper_diagonal;
}

vector<vector<arg_real_t>> distance_matrix_v2(const ARG& arg, arg_real_t alpha,
                             arg_real_t from_pos, arg_real_t to_pos) {
  // Sample IDs must be consecutive
  for (size_t i = 0; i < arg.leaf_ids.size(); ++i) {
    assert(arg.leaf_ids.find(i) != arg.leaf_ids.end());
  }
  // checking roots is relatively quick compared to checking children / parents
  try {
    arg.check_roots();
  }
  catch (...) {
    std::string my_string = "ARG doesn't have correct root information.";
    my_string += " Make sure to call populate_children_and_roots().";
    throw std::logic_error(THROW_LINE(my_string));
  }

  int n = arg.leaf_ids.size();
  vector<vector<arg_real_t>> upper_diagonal;
  for (int i = 0; i < n - 1; ++i) {
    upper_diagonal.push_back(vector<arg_real_t>(n - 1 - i, 0));
  }

  // because we get out the values, which is non-const, we can't use
  // const DescendantList&
  arg_utils::visit_branches(
      arg, [&upper_diagonal, alpha](DescendantList& desc_list, DescendantList& parent_desc_list,
                                    const ARGNode* parent, const ARGNode* child, arg_real_t start,
                                    arg_real_t end) {
        (void) parent_desc_list;
        int n = upper_diagonal.size() + 1; // avoids having to capture more things

        arg_real_t af = ((arg_real_t) desc_list.values().size()) / ((arg_real_t) n);
        arg_real_t maf = std::min(af, 1 - af);
        // compute volume of this branch alone
        arg_real_t volume = (parent->height - child->height) * (end - start);
        arg_real_t factor = pow(maf * (1 - maf), alpha);
        volume *= factor * 0.5; // divide by 2 so it matches with the other distance matrix

        // compute list of complement IDs
        vector<int> complement;
        int i = 0;
        for (int v : desc_list.values()) {
          while (i < v) {
            complement.push_back(i);
            ++i;
          }
          ++i;
        }
        while (i < n) {
          complement.push_back(i);
          ++i;
        }

        for (int v1 : desc_list.values()) {
          for (int v2 : complement) {
            assert(v1 != v2);
            int i = std::min(v1, v2);
            int j = std::max(v1, v2);
            upper_diagonal[i][j - i - 1] += volume;
          }
        }
      }, from_pos, to_pos);

  return upper_diagonal;
}

vector<vector<vector<arg_real_t>>> distance_matrix_maf_bins(const ARG& arg,
                                                            vector<arg_real_t> maf_bins) {
  // Sample IDs must be consecutive
  for (size_t i = 0; i < arg.leaf_ids.size(); ++i) {
    assert(arg.leaf_ids.find(i) != arg.leaf_ids.end());
  }
  // checking roots is relatively quick compared to checking children / parents
  try {
    arg.check_roots();
  }
  catch (...) {
    std::string my_string = "ARG doesn't have correct root information.";
    my_string += " Make sure to call populate_children_and_roots().";
    throw std::logic_error(THROW_LINE(my_string));
  }

  if ((maf_bins.size() == 0) || (maf_bins[0] != 0)) {
    throw std::invalid_argument(THROW_LINE("MAF bins must start at 0"));
  }
  int num_bins = maf_bins.size();

  int n = arg.leaf_ids.size();
  vector<vector<vector<arg_real_t>>> upper_diagonals(num_bins, vector<vector<arg_real_t>>());
  for (int k = 0; k < num_bins; ++k) {
    for (int i = 0; i < n - 1; ++i) {
      upper_diagonals[k].push_back(vector<arg_real_t>(n - 1 - i, 0));
    }
  }

  // because we get out the values, which is non-const, we can't use
  // const DescendantList&
  arg_utils::visit_branches(
      arg, [&upper_diagonals, maf_bins](DescendantList& desc_list, DescendantList& parent_desc_list,
                                        const ARGNode* parent, const ARGNode* child,
                                        arg_real_t start, arg_real_t end) {
        (void) parent_desc_list;
        int n = upper_diagonals[0].size() + 1; // avoids having to capture more things
        int num_bins = maf_bins.size();

        arg_real_t af = ((arg_real_t) desc_list.values().size()) / ((arg_real_t) n);
        arg_real_t maf = std::min(af, 1 - af);
        // get index of which MAF bin this branch falls in
        // binary search is faster for many bins, probably not necessary here
        int index = 0;
        while ((index != num_bins - 1) && (maf >= maf_bins[index + 1])) {
          ++index;
        }
        // compute volume of this branch alone
        arg_real_t volume = (parent->height - child->height) * (end - start);

        // compute list of complement IDs
        vector<int> complement;
        int i = 0;
        for (int v : desc_list.values()) {
          while (i < v) {
            complement.push_back(i);
            ++i;
          }
          ++i;
        }
        while (i < n) {
          complement.push_back(i);
          ++i;
        }

        for (int v1 : desc_list.values()) {
          for (int v2 : complement) {
            assert(v1 != v2);
            int i = std::min(v1, v2);
            int j = std::max(v1, v2);
            upper_diagonals[index][i][j - i - 1] += volume;
          }
        }
      });

#ifdef _DEBUG
  for (int k = 0; k < num_bins; ++k) {
    cout << "Bin starting at MAF = " << maf_bins[k] << endl;
    for (int i = 0; i < n - 1; ++i) {
      for (auto j : upper_diagonals[k][i]) {
        cout << j << " ";
      }
      cout << endl;
    }
    cout << endl;
  }
#endif // _DEBUG

  return upper_diagonals;
}

// Diploid association inside the ARG, testing all clades.
//
// Writes results to file and additionally returns max chi2 over all tests.
//
// Arguments:
//   const ARG& arg: the ARG, with 2*n leaves (haploid samples)
//   const vector<real>& raw_phenotypes: length n containing the phenotypes
//       before any standardization,
//   const deque<bool>& use_sample: length n with true meaning we should include
//       the (diploid) sample in association, and false meaning to exclude. We
//       could have also coded this in the phenotype as a -9 value like PLINK,
//       but this is more principled
//   string file_root: association results are written to file_root + ".tab.gz",
//       while any bitsets below the write_bitset_threshold are additionally
//       written to file_root + ".haps.gz"
//   int chromosome: used in file output
//   string snp_prefix: used in file output, SNPs are labeled as snp_prefix +
//       str(i) where i goes from 1 to the number of clades tested
//   real min_maf: only bitsets with MAC >= round(sum(use_sample) * 2 * min_maf)
//       and additionally MAC > 0 are tested, in addition to any effects from
//       max_maf. If min_maf <= 0, it is automatically set to 0. Default
//       value = -1.
//   real max_maf: only bitsets with MAC <= round(sum(use_sample) * 2 * max_maf)
//       are tested, in addition to any effects from min_maf. If max_maf <= 0,
//       it is automatically set to 0.5. Default value = -1.
//   real write_bitset_threshold: any bitsets with p-value less than or equal to
//       the threshold are written to file_root + ".haps.gz". Currently defaults
//       to -1, which corresponds to not writing anything. This is useful if we
//       want to follow up the most significant bitsets. Setting to 1 means we
//       write everything.
//   real calibration_factor: uncalibrated chi-squared statistics are multiplied
//       by the calibration factor to yield final chi-squared statistics. Used
//       to replicate BOLT-LMM.
//   bool concise_pvalue: if false, then print a verbose pvalue, like 0.12345.
//       Otherwise, use BOLT's print formatting which computes two significant
//       digits, like 1.2E-1. BOLT also has the advantage of being able to
//       calculate p-values less than the minimum double of roughly 1e-300.
//       Defaults to true.
//   bool max_only: if false, then does not write anything to file, instead just
//       returning the max chi2 only. Defaults to true, in which case both file
//       writing and returning are on.
//   bool careful: if true, then computes a residual sum of squares from the
//       inferred beta as part of the standard error. If false, just uses norm
//       of phenotype as an approximation, like BOLT. The careful version will
//       give slightly more significant p-values, especially for significant
//       variants. Defaults to false.
arg_real_t association_diploid_all(const ARG& arg, const vector<arg_real_t>& raw_phenotypes,
                                   const deque<bool>& use_sample, string file_root, int chromosome,
                                   string snp_prefix, arg_real_t min_maf, arg_real_t max_maf,
                                   arg_real_t write_bitset_threshold, arg_real_t calibration_factor,
                                   bool concise_pvalue, bool max_only, bool careful) {
  // Sample IDs must be consecutive
  for (size_t i = 0; i < arg.leaf_ids.size(); ++i) {
    assert(arg.leaf_ids.find(i) != arg.leaf_ids.end());
  }
  // checking roots is relatively quick compared to checking children / parents
  try {
    arg.check_roots();
  }
  catch (...) {
    std::string my_string = "ARG doesn't have correct root information.";
    my_string += " Make sure to call populate_children_and_roots().";
    throw std::logic_error(THROW_LINE(my_string));
  }

  size_t n = raw_phenotypes.size();
  assert(arg.leaf_ids.size() % 2 == 0); // redundant but keeping it in
  assert(arg.leaf_ids.size() == 2 * n);
  assert(use_sample.size() == n);
  if (min_maf <= 0) {
    min_maf = 0;
  }
  if (max_maf <= 0) {
    max_maf = 0.5;
  }

  // standardize the raw_phenotypes and put into phenotypes
  // only deals with entries for which use_sample[i] is true
  vector<arg_real_t> phenotypes = utils::standardize_mask(raw_phenotypes, use_sample);
  arg_real_t active_n = 0;
  for (size_t i = 0; i < n; ++i) {
    active_n += use_sample[i]; // casts bool to 0 or 1 real value
  }

  file_utils::AutoGzOfstream results_out, haps_out;
  if (!max_only) {
    if (file_root == "") {
      throw std::logic_error(THROW_LINE("Expecting a non-empty file_root string"));
    }
    results_out.openOrExit(file_root + ".tab.gz");
    haps_out.openOrExit(file_root + ".haps.gz");
    results_out << "SNP\tCHR\tBP\tSTART_BP\tEND_BP\tMAC\tMAF\tCHISQ\tP\tSIGN"
                << "\n";
    // because we visit clades in this version, we can't get out the branch top and bottom
  }
  arg_real_t arg_offset = arg.offset;
  arg_real_t max_chi2 = 0;
  int count = 0;
  int significant_count = 0;

  // because we get out the values, which is non-const, we can't use
  // const DescendantList&
  arg_utils::visit_clades(
      arg, [&max_chi2, &phenotypes, &use_sample, &count, &arg_offset, &results_out, &haps_out,
            &significant_count, write_bitset_threshold, concise_pvalue, calibration_factor,
            max_only, file_root, snp_prefix, chromosome, min_maf, max_maf, careful, active_n](
               DescendantList& desc_list, const ARGNode* node, arg_real_t start, arg_real_t end) {
        (void) node;
        arg_real_t one_count = 0;   // count of descendant i for which i/2 is a sample
        arg_real_t one_sum = 0;     // sum of y[i/2] for descendant i and i/2 a sample
        arg_real_t cross_count = 0; // count of sample j for which 2*j and 2*j+1 are descendants

        int last_value = -5; // put a big gap to ensure correctness
        for (int v : desc_list.values()) {
          if (use_sample[v / 2]) {
            one_count += 1;
            one_sum += phenotypes[v / 2];
            if ((v ^ 1) == last_value) { // checks that last_value = 2*j and v = 2*j+1
              cross_count += 1;
            }
            last_value = v;
          }
        }
        int sign = 0;
        if (one_sum > 0) {
          sign = 1;
        }
        if (one_sum < 0) {
          sign = -1;
        }

        int mac = round(std::min(one_count, 2 * active_n - one_count));
        int min_mac = round(2 * active_n * min_maf);
        int max_mac = round(2 * active_n * max_maf);
        // MAF filter, also don't want roots
        if ((mac > 0) && (mac >= min_mac) && (mac <= max_mac)) {
          arg_real_t norm_y_sq = active_n - 1;
          arg_real_t norm_x_sq = one_count - one_count * one_count / active_n + 2 * cross_count;
          arg_real_t chi_sq_uncalibrated;
          if (careful) {
            arg_real_t denominator = norm_y_sq * norm_x_sq / (one_sum * one_sum) - 1;
            chi_sq_uncalibrated = (active_n - 2) / denominator;
            // the above is equivalent to the below verbose alternative
            // real beta_hat = one_sum / norm_x_sq;
            // real rss = norm_y_sq - beta_hat * beta_hat * norm_x_sq;
            // real sigma_sq_hat = rss / (active_n - 2);
            // real se_hat = sqrt(sigma_sq_hat / norm_x_sq);
            // real z_score = beta_hat / se_hat;
            // chi_sq_uncalibrated = z_score * z_score;
          }
          else {
            arg_real_t denominator = norm_y_sq * norm_x_sq / (one_sum * one_sum);
            // chi_sq_uncalibrated = (active_n - 2) / denominator; // this is more correct
            chi_sq_uncalibrated = (active_n - 1) / denominator; // this exactly reproduces BOLT
          }
          arg_real_t chi_sq = chi_sq_uncalibrated / calibration_factor;

          // turn chi_sq into p-value
          boost::math::chi_squared chi_sq_dist(1);
          arg_real_t p_value = boost::math::cdf(complement(chi_sq_dist, chi_sq));

          // Note: this block of code is adapted from the BOLT-LMM v2.3.2 source code,
          // developed by Po-Ru Loh and released under the GNU General Public License
          // v3.0 (GPLv3).
          //
          // The license file can be found at 3rd_party/BOLT-LMM_v2.3.2/license.txt
          // from the root of this repository.
          //
          // The else condition allows for printing p-values less than 1e-300 that
          // otherwise show up as 0, by computing the factor and exponent directly
          // from chi_sq
          char pValueBuf[100];
          if (concise_pvalue) {
            if (p_value != 0)
              sprintf(pValueBuf, "%.1E", p_value);
            else {
              double log10p = log10(2.0) - M_LOG10E * chi_sq / 2 - 0.5 * log10(chi_sq * 2 * M_PI);
              int exponent = floor(log10p);
              double fraction = pow(10.0, log10p - exponent);
              if (fraction >= 9.95) {
                fraction = 1;
                ++exponent;
              }
              sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
            }
          }

          if (chi_sq > max_chi2) {
            max_chi2 = chi_sq;
          }
          if (!max_only) {
            ++count; // 1-indexed
            int start_int = (int) start + arg_offset;
            int end_int = (int) end + arg_offset;
            int mid_int = (int) ((start + end) / 2 + arg_offset);
            arg_real_t maf = mac / (2 * active_n);
            results_out << snp_prefix << ":" << count << "\t" << chromosome << "\t";
            results_out << mid_int << "\t" << start_int << "\t" << end_int << "\t";
            results_out << mac << "\t" << maf << "\t";
            results_out << chi_sq << "\t";
            if (concise_pvalue) {
              results_out << string(pValueBuf) << "\t" << sign << "\n";
            }
            else {
              results_out << p_value << "\t" << sign << "\n";
            }

            if (p_value < write_bitset_threshold) {
              ++significant_count;
              haps_out << chromosome << " " << snp_prefix << ":" << count << " " << mid_int
                       << " 0 1";
              for (char c : desc_list.to_bitset_string()) {
                haps_out << " " << c;
              }
              haps_out << "\n";
            }
          }
        }
      });

  if (!max_only) {
    results_out.close();
    haps_out.close();
  }

  return max_chi2;
}

// Diploid association inside the ARG, using mutation-based sampling.
//
// Writes results to file and additionally returns max chi2 over all tests
// for each mutation rate.
//
// Arguments:
//   const ARG& arg: the ARG, with 2*n leaves (haploid samples)
//   const vector<real>& raw_phenotypes: length n containing the phenotypes
//       before any standardization,
//   const deque<bool>& use_sample: length n with true meaning we should include
//       the (diploid) sample in association, and false meaning to exclude. We
//       could have also coded this in the phenotype as a -9 value like PLINK,
//       but this is more principled
//   string file_root: association results are written to file_root + ".tab.gz",
//       while any bitsets below the write_bitset_threshold are additionally
//       written to file_root + ".haps.gz"
//   vector<real> mus: a vector of mutation rates to test for. The largest mutation
//       rate is then primarily used, but when returning results, the maximum chi2
//       is returned for each of the mutation rates in order
//   unsigned random_seed: seed for mutations! When equal to 0, uses the time
//       to seed instead of the seed value.
//   int chromosome: used in file output
//   string snp_prefix: used in file output, SNPs are labeled as snp_prefix +
//       str(i) where i goes from 1 to the number of clades tested
//   real min_maf: only bitsets with MAC >= round(sum(use_sample) * 2 * min_maf)
//       and additionally MAC > 0 are tested, in addition to any effects from
//       max_maf. If min_maf <= 0, it is automatically set to 0. Default
//       value = -1.
//   real max_maf: only bitsets with MAC <= round(sum(use_sample) * 2 * max_maf)
//       are tested, in addition to any effects from min_maf. If max_maf <= 0,
//       it is automatically set to 0.5. Default value = -1.
//   real write_bitset_threshold: any bitsets with p-value less than or equal to
//       the threshold are written to file_root + ".haps.gz". Currently defaults
//       to -1, which corresponds to not writing anything. This is useful if we
//       want to follow up the most significant bitsets. Setting to 1 means we
//       write everything.
//   real calibration_factor: uncalibrated chi-squared statistics are multiplied
//       by the calibration factor to yield final chi-squared statistics. Used
//       to replicate BOLT-LMM.
//   bool concise_pvalue: if false, then print a verbose pvalue, like 0.12345.
//       Otherwise, use BOLT's print formatting which computes two significant
//       digits, like 1.2E-1. BOLT also has the advantage of being able to
//       calculate p-values less than the minimum double of roughly 1e-300.
//       Defaults to true.
//   bool max_only: if false, then does not write anything to file, instead just
//       returning the max chi2 only. Defaults to true, in which case both file
//       writing and returning are on.
//   bool careful: if true, then computes a residual sum of squares from the
//       inferred beta as part of the standard error. If false, just uses norm
//       of phenotype as an approximation, like BOLT. The careful version will
//       give slightly more significant p-values, especially for significant
//       variants. Defaults to false.
vector<arg_real_t> association_diploid_mutation(
    const ARG& arg, const vector<arg_real_t>& raw_phenotypes, const deque<bool>& use_sample,
    string file_root, vector<arg_real_t> mus, unsigned random_seed, int chromosome,
    string snp_prefix, arg_real_t min_maf, arg_real_t max_maf, arg_real_t write_bitset_threshold,
    arg_real_t calibration_factor, bool concise_pvalue, bool max_only, bool careful) {
  // Sample IDs must be consecutive
  for (size_t i = 0; i < arg.leaf_ids.size(); ++i) {
    assert(arg.leaf_ids.find(i) != arg.leaf_ids.end());
  }
  // checking roots is relatively quick compared to checking children / parents
  try {
    arg.check_roots();
  }
  catch (...) {
    std::string my_string = "ARG doesn't have correct root information.";
    my_string += " Make sure to call populate_children_and_roots().";
    throw std::logic_error(THROW_LINE(my_string));
  }

  if (random_seed == 0) {
    random_seed = std::chrono::system_clock::now().time_since_epoch().count();
  }
  std::mt19937 generator(random_seed);

  size_t n = raw_phenotypes.size();
  assert(arg.leaf_ids.size() % 2 == 0); // redundant but keeping it in
  assert(arg.leaf_ids.size() == 2 * n);
  assert(use_sample.size() == n);
  if (min_maf <= 0) {
    min_maf = 0;
  }
  if (max_maf <= 0) {
    max_maf = 0.5;
  }

  if (mus.size() == 0) {
    throw std::logic_error(THROW_LINE("Must pass an array of mutation rates"));
  }
  for (size_t i = 0; i < mus.size(); ++i) {
    if (mus[i] < 0) {
      throw std::logic_error(THROW_LINE("Mutation rates must be positive"));
    }
  }

  // standardize the raw_phenotypes and put into phenotypes
  // only deals with entries for which use_sample[i] is true
  vector<arg_real_t> phenotypes = utils::standardize_mask(raw_phenotypes, use_sample);
  arg_real_t active_n = 0;
  for (size_t i = 0; i < n; ++i) {
    active_n += use_sample[i]; // casts bool to 0 or 1 real value
  }

  file_utils::AutoGzOfstream results_out, haps_out;
  if (!max_only) {
    if (file_root == "") {
      throw std::logic_error(THROW_LINE("Expecting a non-empty file_root string"));
    }
    results_out.openOrExit(file_root + ".tab.gz");
    haps_out.openOrExit(file_root + ".haps.gz");
    results_out << "SNP\tCHR\tBP\tSTART_BP\tEND_BP\tMAC\tMAF\tCHISQ\tP\tSIGN"
                   "\tMU_STAR\tHEIGHT_CHILD\tHEIGHT_PARENT\tBRANCH_VOL"
                << "\n";
  }
  arg_real_t arg_offset = arg.offset;
  vector<arg_real_t> max_chi2s(mus.size(), 0);
  int count = 0;
  int significant_count = 0;

  // because we get out the values, which is non-const, we can't use
  // const DescendantList&
  arg_utils::visit_branches(
      arg,
      [&max_chi2s, &phenotypes, &use_sample, &count, &arg_offset, &results_out, &haps_out,
       &significant_count, write_bitset_threshold, concise_pvalue, &generator, &mus,
       calibration_factor, max_only, file_root, snp_prefix, chromosome, min_maf, max_maf, careful,
       active_n](DescendantList& desc_list, DescendantList& _unused, const ARGNode* parent_node,
                 const ARGNode* child_node, arg_real_t start, arg_real_t end) {
        (void) _unused;
        arg_real_t branch_volume = (parent_node->height - child_node->height) * (end - start);
        // Now we sample mu_star ~ Expo(branch_volume). This comes from thinking about a Poisson
        // process with rate given by the volume and asking how much mu we need until we get one
        // event (mutation).
        arg_real_t mu_star = random_utils::generate_exponential_rv(generator, branch_volume);
        ;
        arg_real_t max_mu = *std::max_element(mus.begin(), mus.end());

        if (max_mu > mu_star) {
          arg_real_t one_count = 0;   // count of descendant i for which i/2 is a sample
          arg_real_t one_sum = 0;     // sum of y[i/2] for descendant i and i/2 a sample
          arg_real_t cross_count = 0; // count of sample j for which 2*j and 2*j+1 are descendants

          int last_value = -5; // put a big gap to ensure correctness
          for (int v : desc_list.values()) {
            if (use_sample[v / 2]) {
              one_count += 1;
              one_sum += phenotypes[v / 2];
              if ((v ^ 1) == last_value) { // checks that last_value = 2*j and v = 2*j+1
                cross_count += 1;
              }
              last_value = v;
            }
          }
          int sign = 0;
          if (one_sum > 0) {
            sign = 1;
          }
          if (one_sum < 0) {
            sign = -1;
          }

          int mac = round(std::min(one_count, 2 * active_n - one_count));
          int min_mac = round(2 * active_n * min_maf);
          int max_mac = round(2 * active_n * max_maf);
          // MAF filter, also don't want roots
          if ((mac > 0) && (mac >= min_mac) && (mac <= max_mac)) {
            arg_real_t norm_y_sq = active_n - 1;
            arg_real_t norm_x_sq = one_count - one_count * one_count / active_n + 2 * cross_count;
            arg_real_t chi_sq_uncalibrated;
            if (careful) {
              arg_real_t denominator = norm_y_sq * norm_x_sq / (one_sum * one_sum) - 1;
              chi_sq_uncalibrated = (active_n - 2) / denominator;
              // the above is equivalent to the below verbose alternative
              // real beta_hat = one_sum / norm_x_sq;
              // real rss = norm_y_sq - beta_hat * beta_hat * norm_x_sq;
              // real sigma_sq_hat = rss / (active_n - 2);
              // real se_hat = sqrt(sigma_sq_hat / norm_x_sq);
              // real z_score = beta_hat / se_hat;
              // chi_sq_uncalibrated = z_score * z_score;
            }
            else {
              arg_real_t denominator = norm_y_sq * norm_x_sq / (one_sum * one_sum);
              // chi_sq_uncalibrated = (active_n - 2) / denominator; // this is more correct
              chi_sq_uncalibrated = (active_n - 1) / denominator; // this exactly reproduces BOLT
            }
            arg_real_t chi_sq = chi_sq_uncalibrated / calibration_factor;

            // turn chi_sq into p-value
            boost::math::chi_squared chi_sq_dist(1);
            arg_real_t p_value = boost::math::cdf(complement(chi_sq_dist, chi_sq));

            // Note: this block of code is adapted from the BOLT-LMM v2.3.2 source code,
            // developed by Po-Ru Loh and released under the GNU General Public License
            // v3.0 (GPLv3).
            //
            // The license file can be found at 3rd_party/BOLT-LMM_v2.3.2/license.txt
            // from the root of this repository.
            //
            // The else condition allows for printing p-values less than 1e-300 that
            // otherwise show up as 0, by computing the factor and exponent directly
            // from chi_sq
            char pValueBuf[100];
            if (concise_pvalue) {
              if (p_value != 0)
                sprintf(pValueBuf, "%.1E", p_value);
              else {
                double log10p = log10(2.0) - M_LOG10E * chi_sq / 2 - 0.5 * log10(chi_sq * 2 * M_PI);
                int exponent = floor(log10p);
                double fraction = pow(10.0, log10p - exponent);
                if (fraction >= 9.95) {
                  fraction = 1;
                  ++exponent;
                }
                sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
              }
            }

            // update the max_chi2s for those with mu crossing mu_star
            for (size_t i = 0; i < mus.size(); ++i) {
              if (mus[i] > mu_star) {
                if (chi_sq > max_chi2s[i]) {
                  max_chi2s[i] = chi_sq;
                }
              }
            }
            if (!max_only) {
              ++count; // 1-indexed
              int start_int = (int) start + arg_offset;
              int end_int = (int) end + arg_offset;
              int mid_int = (int) ((start + end) / 2 + arg_offset);
              arg_real_t maf = mac / (2 * active_n);
              results_out << snp_prefix << ":" << count << "\t" << chromosome << "\t";
              results_out << mid_int << "\t" << start_int << "\t" << end_int << "\t";
              results_out << mac << "\t" << maf << "\t";
              results_out << chi_sq << "\t";
              if (concise_pvalue) {
                results_out << string(pValueBuf) << "\t";
              }
              else {
                results_out << p_value << "\t";
              }
              results_out << sign;
              // these two lines are output that association_diploid_all lacks
              results_out << "\t" << mu_star << "\t" << child_node->height;
              results_out << "\t" << parent_node->height << "\t" << branch_volume;
              results_out << "\n";

              if (p_value < write_bitset_threshold) {
                ++significant_count;
                haps_out << chromosome << " " << snp_prefix << ":" << count << " " << mid_int
                         << " 0 1";
                for (char c : desc_list.to_bitset_string()) {
                  haps_out << " " << c;
                }
                haps_out << "\n";
              }
            }
          }
        }
      });

  if (!max_only) {
    results_out.close();
    haps_out.close();
  }
  return max_chi2s;
}

// Haploid association that does all clades
//
// Returns top K chi2, volume, node_height, start, end, mean, num_descendants.
// Doesn't write to disk, spends some extra time computing clade volume.
vector<tuple<arg_real_t, arg_real_t, arg_real_t, arg_real_t, arg_real_t, arg_real_t, int>>
association_haploid(const ARG& arg, const vector<arg_real_t>& raw_phenotypes, int topk) {
  // Sample IDs must be consecutive
  for (size_t i = 0; i < arg.leaf_ids.size(); ++i) {
    assert(arg.leaf_ids.find(i) != arg.leaf_ids.end());
  }
  // checking roots is relatively quick compared to checking children / parents
  try {
    arg.check_roots();
  }
  catch (...) {
    std::string my_string = "ARG doesn't have correct root information.";
    my_string += " Make sure to call populate_children_and_roots().";
    throw std::logic_error(THROW_LINE(my_string));
  }

  vector<tuple<arg_real_t, arg_real_t, arg_real_t, arg_real_t, arg_real_t, arg_real_t, int>>
      results;
  assert(raw_phenotypes.size() == arg.leaf_ids.size());
  int n = arg.leaf_ids.size();

  // standardize the raw_phenotypes and put into phenotypes
  vector<arg_real_t> phenotypes = utils::standardize(raw_phenotypes);

  // because we get out the values, which is non-const, we can't use
  // const DescendantList&
  arg_utils::visit_clades(
      arg, [&results, &phenotypes, n](
               DescendantList& desc_list, const ARGNode* node, arg_real_t start, arg_real_t end) {
        arg_real_t one_sum = 0;
        int num_descendants = 0;
        for (int v : desc_list.values()) {
          one_sum += phenotypes[v];
          num_descendants += 1;
        }
        // don't want roots
        if (num_descendants != n) {
          arg_real_t zero_sum = -one_sum; // because phenotype was mean-centered
          arg_real_t p = ((arg_real_t) num_descendants) / ((arg_real_t) n);

          arg_real_t beta = (-p) * zero_sum + (1 - p) * one_sum;
          arg_real_t x_std = sqrt(p * (1 - p) * n / (n - 1));
          beta /= x_std * n;
          arg_real_t se2 = 1.0 / (n - 2);
          arg_real_t chi2 = (beta * beta) / se2;

          vector<tuple<arg_real_t, arg_real_t, arg_real_t>> segments;
          arg_real_t pos = start;
          while (pos < end) {
            const ARGEdge* edge = node->parent_edge_at(pos);
            arg_real_t edge_end = edge->end;
            arg_real_t branch_height = edge->parent->height - node->height;
            arg_real_t next_pos;
            if (edge_end > end) {
              next_pos = end;
            }
            else {
              next_pos = edge_end;
            }
            segments.push_back(std::make_tuple(pos, next_pos, branch_height));
            pos = next_pos;
          }

          arg_real_t clade_volume = 0;
          arg_real_t clade_mean = 0;
          for (auto segment : segments) {
            arg_real_t segment_start = std::get<0>(segment);
            arg_real_t segment_end = std::get<1>(segment);
            arg_real_t branch_height = std::get<2>(segment);
            arg_real_t segment_volume = (segment_end - segment_start) * branch_height;
            arg_real_t segment_midpoint = (segment_start + segment_end) / 2;
            clade_volume += segment_volume;
            clade_mean += segment_midpoint * segment_volume;
          }
          clade_mean /= clade_volume;

          results.push_back(std::make_tuple(
              chi2, clade_volume, node->height, start, end, clade_mean, num_descendants));
        }
      });

  // sort and return results
  if (topk <= 0 || topk >= (int) results.size()) {
    std::sort(
        results.begin(), results.end(),
        std::greater<
            tuple<arg_real_t, arg_real_t, arg_real_t, arg_real_t, arg_real_t, arg_real_t, int>>());
    return results;
  }
  // By using partial_sort instead of sort, we are O(N log K) instead of O(N log N)
  size_t actual_k = std::min<size_t>(topk, results.size());
  std::partial_sort(
      results.begin(), results.begin() + actual_k, results.end(),
      std::greater<
          tuple<arg_real_t, arg_real_t, arg_real_t, arg_real_t, arg_real_t, arg_real_t, int>>());
  results.resize(actual_k);
  return results;
}

// This function returns data to represent the max chi^2 at each genomic position
// The returned vector consists of segments that span the ARG and the max chi^2
// that can be achieved over that region. vector of (start, end, max chi^2)
//
// compress corresponds to whether to merge neighboring regions with the same
// chi2, in which case epsilon is used to determine whether two chi2 are the
// same (via comparing absolute difference)
vector<tuple<arg_real_t, arg_real_t, arg_real_t>>
association_haploid_max(const ARG& arg, const vector<arg_real_t>& raw_phenotypes, bool compress,
                        arg_real_t epsilon) {
  // Sample IDs must be consecutive
  for (size_t i = 0; i < arg.leaf_ids.size(); ++i) {
    assert(arg.leaf_ids.find(i) != arg.leaf_ids.end());
  }
  // checking roots is relatively quick compared to checking children / parents
  try {
    arg.check_roots();
  }
  catch (...) {
    std::string my_string = "ARG doesn't have correct root information.";
    my_string += " Make sure to call populate_children_and_roots().";
    throw std::logic_error(THROW_LINE(my_string));
  }

  assert(raw_phenotypes.size() == arg.leaf_ids.size());
  int n = arg.leaf_ids.size();

  // standardize the raw_phenotypes and put into phenotypes
  vector<arg_real_t> phenotypes = utils::standardize(raw_phenotypes);

  // results contains (start, end, chi-squared), originally all zeros, one for each region?
  vector<tuple<arg_real_t, arg_real_t, arg_real_t>> results;
  arg_real_t prev = -1;
  for (arg_real_t bp : arg.get_breakpoints()) {
    if (prev != -1) {
      results.push_back(std::make_tuple(prev, bp, 0));
    }
    prev = bp;
  }

  // because we get out the values, which is non-const, we can't use
  // const DescendantList&
  arg_utils::visit_clades(arg, [&results, &phenotypes, n](DescendantList& desc_list,
                                                          const ARGNode* node, arg_real_t start,
                                                          arg_real_t end) {
    (void) node;
    arg_real_t one_sum = 0;
    int num_descendants = 0;
    for (int v : desc_list.values()) {
      one_sum += phenotypes[v];
      num_descendants += 1;
    }
    // don't want roots
    if (num_descendants != n) {
      arg_real_t zero_sum = -one_sum; // because phenotype was mean-centered
      arg_real_t p = ((arg_real_t) num_descendants) / ((arg_real_t) n);

      arg_real_t beta = (-p) * zero_sum + (1 - p) * one_sum;
      arg_real_t x_std = sqrt(p * (1 - p) * n / (n - 1));
      beta /= x_std * n;
      arg_real_t se2 = 1.0 / (n - 2);
      arg_real_t chi2 = (beta * beta) / se2;

      // go from start to end and update results
      auto comp = [](tuple<arg_real_t, arg_real_t, arg_real_t> first,
                     tuple<arg_real_t, arg_real_t, arg_real_t> second) {
        return std::get<0>(first) < std::get<0>(second);
      };
      for (auto it =
               std::lower_bound(results.begin(), results.end(), std::make_tuple(start, 0, 0), comp);
           it != std::lower_bound(results.begin(), results.end(), std::make_tuple(end, 0, 0), comp);
           ++it) {
        if (std::get<2>(*it) < chi2) {
          std::get<2>(*it) = chi2;
        }
      }
    }
  });

  if (compress) {
    // merge neighbors with the same chi2
    vector<tuple<arg_real_t, arg_real_t, arg_real_t>> new_results;
    arg_real_t same_start = arg.start;
    arg_real_t last_chi2 = -1;
    for (auto result : results) {
      if (std::fabs(last_chi2 - std::get<2>(result)) > epsilon) {
        if (last_chi2 != -1) {
          new_results.push_back(std::make_tuple(same_start, std::get<0>(result), last_chi2));
        }
        same_start = std::get<0>(result);
        last_chi2 = std::get<2>(result);
      }
    }
    new_results.push_back(std::make_tuple(same_start, arg.end, last_chi2));
    return new_results;
  }

  return results;
}

// Write bitsets in .haps / .sample / .bitinfo format for converting to PLINK
//
// Ends up storing all bitsets / information in memory
//
// The bitset position is just a dummy index, and the SNP name is snp_prefix:index,
// or chromosome:index if snp_prefix is empty
//
// file_root: if empty, then outputs to stdout. Otherwise, writes to file
// compress: say we have an identical bitset that splits {[50, 100], [100, 200], [250, 300]}.
//     If true, then we write {[50, 200], [250, 300]} in *.bitinfo. If false,
//     we write all three in *.bitinfo.
// count_only: if true, then doesn't output anything, just returns a count
size_t write_bitsets_detailed(const ARG& arg, string file_root, bool diploid, int chromosome,
                              string snp_prefix, bool compress, bool count_only) {
  // Sample IDs must be consecutive
  for (size_t i = 0; i < arg.leaf_ids.size(); ++i) {
    assert(arg.leaf_ids.find(i) != arg.leaf_ids.end());
  }
  // checking roots is relatively quick compared to checking children / parents
  try {
    arg.check_roots();
  }
  catch (...) {
    std::string my_string = "ARG doesn't have correct root information.";
    my_string += " Make sure to call populate_children_and_roots().";
    throw std::logic_error(THROW_LINE(my_string));
  }

  size_t n = arg.leaf_ids.size();
  arg_real_t offset = arg.offset;
  if (diploid) {
    assert(n % 2 == 0);
  }
  if (snp_prefix == "") {
    snp_prefix = std::to_string(chromosome);
  }

  unordered_map<DescendantList, vector<vector<arg_real_t>>, DescendantListHash> bitset_extents;

  // because we get out the values, which is non-const, we can't use
  // const DescendantList&
  arg_utils::visit_clades(arg, [&bitset_extents, n](DescendantList& desc_list, const ARGNode* node,
                                                    arg_real_t start, arg_real_t end) {
    (void) node;
    // don't want roots
    if (desc_list.num_values() != n) {
      vector<vector<arg_real_t>>& entry =
          bitset_extents[desc_list]; // creates if not present, only hashes once
      entry.push_back(vector<arg_real_t>{start, end});
    }
  });

  if (!count_only) {
    // open files
    std::ofstream haps_out, sample_out, bitinfo_out;
    if (file_root != "") {
      haps_out.open(file_root + ".haps");
      sample_out.open(file_root + ".sample");
      bitinfo_out.open(file_root + ".bitinfo");
    }

    // write out everything we need
    int i = 0;
    for (const auto& element : bitset_extents) {
      ++i; // 1-indexed
      if (file_root == "") {
        cout << element.first.to_bitset_string() << " ";
      }
      else {
        haps_out << chromosome << " " << snp_prefix << ":" << i << " " << i << " 0 1";
        if (diploid) {
          for (char c : element.first.to_bitset_string()) {
            haps_out << " " << c;
          }
        }
        else {
          for (char c : element.first.to_bitset_string()) {
            // repeat the value for .haps file
            haps_out << " " << c << " " << c;
          }
        }
        haps_out << endl;
        bitinfo_out << chromosome << " " << snp_prefix << ":" << i;
      }

      if (!compress) {
        if (file_root == "") {
          for (auto start_end : element.second) {
            cout << start_end[0] << " " << start_end[1] << " ";
          }
          cout << endl;
        }
        else {
          for (auto start_end : element.second) {
            bitinfo_out << " " << start_end[0] << " " << start_end[1];
          }
          bitinfo_out << endl;
        }
      }
      else {
        arg_real_t last_start = element.second[0][0];
        arg_real_t last_end = element.second[0][0];
        for (auto start_end : element.second) {
          assert(start_end[0] >= last_end);
          if (start_end[0] > last_end) {
            if (file_root == "") {
              cout << last_start + offset << " " << last_end + offset << " ";
            }
            else {
              bitinfo_out << " " << last_start + offset << " " << last_end + offset;
            }
            last_start = start_end[0];
            last_end = start_end[1];
          }
          else {
            last_end = start_end[1];
          }
        }
        if (file_root == "") {
          cout << last_start + offset << " " << last_end + offset << endl;
        }
        else {
          bitinfo_out << " " << last_start + offset << " " << last_end + offset << endl;
        }
      }
    }

    if (file_root != "") {
      sample_out << "ID_1 ID_2 missing" << endl;
      sample_out << "0 0 0" << endl;
      for (size_t i = 1; i <= n; ++i) {
        sample_out << i << " " << i << " 0" << endl;
      }
    }

    // close files
    if (file_root != "") {
      haps_out.close();
      sample_out.close();
      bitinfo_out.close();
    }
  }

  return bitset_extents.size();
}

// Write bitsets in .haps / .sample / .bitinfo format for converting to PLINK
//
// Similar to above, but if a bitset stretch is not consecutive then it gets
// printed multiple times. therefore, it's more memory-efficient
//
// Since we only consider consecutive bitset stretch, we write out the mean
// position. Therefore, positions are not necessarily increasing!
// The SNP name is snp_prefix:index, or chromosome:index if snp_prefix is empty
//
// file_root: if empty, then outputs to stdout. Otherwise, writes to file
// count_only: if true, then doesn't output anything, just returns a count
// chromosome: used for writing
// snp_prefix: used for writing
// min_mac / max_mac: inclusive min / max minor allele count
size_t write_bitsets(const ARG& arg, string file_root, bool diploid, int chromosome,
                     string snp_prefix, size_t min_mac, size_t max_mac, bool write_dosage,
                     bool use_gz, bool count_only) {
  // Sample IDs must be consecutive
  for (size_t i = 0; i < arg.leaf_ids.size(); ++i) {
    assert(arg.leaf_ids.find(i) != arg.leaf_ids.end());
  }
  // checking roots is relatively quick compared to checking children / parents
  try {
    arg.check_roots();
  }
  catch (...) {
    std::string my_string = "ARG doesn't have correct root information.";
    my_string += " Make sure to call populate_children_and_roots().";
    throw std::logic_error(THROW_LINE(my_string));
  }

  size_t n = arg.leaf_ids.size();
  size_t num_samples = n;
  arg_real_t offset = arg.offset;
  if (diploid) {
    assert(n % 2 == 0);
    num_samples = n / 2;
  }
  if (max_mac == 0) {
    max_mac = n;
  }
  if (snp_prefix == "") {
    snp_prefix = std::to_string(chromosome);
  }
  size_t max_pq_size = (size_t) n;

  file_utils::AutoGzOfstream haps_out;
  std::ofstream sample_out, bitinfo_out;
  if (!count_only) {
    // open files
    if (file_root != "") {
      string haps_root;
      if (write_dosage) {
        assert(diploid);
        haps_root = file_root + ".dosage";
        // TODO: write .indivs file?
      }
      else {
        haps_root = file_root + ".haps";
      }
      if (use_gz) {
        haps_root = haps_root + ".gz";
      }
      haps_out.openOrExit(haps_root);
      sample_out.open(file_root + ".sample");
      bitinfo_out.open(file_root + ".bitinfo");

      sample_out << "ID_1 ID_2 missing" << endl;
      sample_out << "0 0 0" << endl;
      for (size_t i = 1; i <= num_samples; ++i) {
        sample_out << i << " " << i << " 0" << endl;
      }
      sample_out.close();
    }
  }

  int count = 0;
  auto comp_first = [](pair<arg_real_t, DescendantList> a, pair<arg_real_t, DescendantList> b) {
    return a.first > b.first;
  };
  std::priority_queue<pair<arg_real_t, DescendantList>, vector<pair<arg_real_t, DescendantList>>,
                      decltype(comp_first)>
      end_pq(comp_first);
  unordered_map<DescendantList, int, DescendantListHash> bitset_tally;
  unordered_map<DescendantList, arg_real_t, DescendantListHash> bitset_start;
  /*
  Pseudocode

  Data structures:
  Priority queue consisting of branch end & bitset, custom comp to just sort by branch end
  Unordered_map of bitset --> appearance count
  Unordered_map of bitset --> consecutive stretch start

  For each branch,
    position = start position
    For all branches that end before position
      If appearance counter goes to 0, write out using consecutive stretch start. Increment global
      count. Remove the branch from the priority queue Insert the branch
  */

  // because we get out the values, which is non-const, we can't use
  // const DescendantList&
  arg_utils::visit_clades(
      arg, [&end_pq, &bitset_tally, &bitset_start, &count, &haps_out, &bitinfo_out, &max_pq_size,
            file_root, n, offset, chromosome, snp_prefix, min_mac, max_mac, write_dosage, diploid,
            count_only](
               DescendantList& desc_list, const ARGNode* node, arg_real_t start, arg_real_t end) {
        (void) node;
        // don't want roots
        if (desc_list.num_values() != n) {
          while ((!end_pq.empty()) && (end_pq.top().first < start)) {
            arg_real_t existing_end = end_pq.top().first;
            DescendantList existing = end_pq.top().second;
            end_pq.pop();
            int& entry = bitset_tally[existing];
            --entry;
            if (entry == 0) {
              arg_real_t existing_start = bitset_start[existing];
              // delete the entries to save memory
              bitset_start.erase(bitset_start.find(existing));
              bitset_tally.erase(bitset_tally.find(existing));

              ++count; // 1-indexed
              max_pq_size = std::max(max_pq_size, end_pq.size());
              if (!count_only) {
                // write out everything we need
                if (file_root == "") {
                  cout << existing.to_bitset_string() << " ";
                  cout << existing_start + offset << " " << existing_end + offset << endl;
                }
                else {
                  int mean_position = (int) ((existing_start + existing_end) * 0.5 + offset);
                  if (write_dosage) {
                    haps_out << snp_prefix << ":" << count << " " << chromosome << " "
                             << mean_position << " 0 1";
                    const dynamic_bitset<>& existing_bitset = existing.bitset();
                    for (size_t i = 0; i < n / 2; ++i) {
                      int sum = ((int) existing_bitset[2 * i]) + ((int) existing_bitset[2 * i + 1]);
                      haps_out << " " << sum;
                    }
                  }
                  else {
                    haps_out << chromosome << " " << snp_prefix << ":" << count << " "
                             << mean_position << " 0 1";
                    if (diploid) {
                      for (char c : existing.to_bitset_string()) {
                        haps_out << " " << c;
                      }
                    }
                    else {
                      for (char c : existing.to_bitset_string()) {
                        // repeat the value for .haps file
                        haps_out << " " << c << " " << c;
                      }
                    }
                  }
                  haps_out << endl;
                  bitinfo_out << chromosome << " " << snp_prefix << ":" << count;
                  bitinfo_out << " " << existing_start + offset << " " << existing_end + offset;
                  if (true) {
                    // also write out number of descendants
                    size_t num_descendants = existing.num_values();
                    bitinfo_out << " " << num_descendants;
                  }
                  bitinfo_out << endl;
                }
              }
            }
          }

          // filter on minor allele count
          size_t bitset_ac = desc_list.num_values();
          size_t bitset_mac = std::min(bitset_ac, n - bitset_ac);
          if ((bitset_mac >= min_mac) && (bitset_mac <= max_mac)) {
            end_pq.push(std::make_pair(end, desc_list));
            int& entry = bitset_tally[desc_list]; // creates if not present, only hashes once
            if (entry == 0) {
              bitset_start[desc_list] = start; // overwrites if present, though it shouldn't be
            }
            ++entry;
          }
        }
      });

  while (!end_pq.empty()) {
    arg_real_t existing_end = end_pq.top().first;
    DescendantList existing = end_pq.top().second;
    end_pq.pop();
    int& entry = bitset_tally[existing];
    --entry;
    if (entry == 0) {
      arg_real_t existing_start = bitset_start[existing];
      // delete the entries to save memory
      bitset_start.erase(bitset_start.find(existing));
      bitset_tally.erase(bitset_tally.find(existing));

      ++count; // 1-indexed
      max_pq_size = std::max(max_pq_size, end_pq.size());
      if (!count_only) {
        // write out everything we need
        if (file_root == "") {
          cout << existing.to_bitset_string() << " ";
          cout << existing_start + offset << " " << existing_end + offset << endl;
        }
        else {
          int mean_position = (int) ((existing_start + existing_end) * 0.5 + offset);
          if (write_dosage) {
            if (write_dosage) { // should this be taken out or replaced with diploid?
              haps_out << snp_prefix << ":" << count << " " << chromosome << " " << mean_position
                       << " 0 1";
              const dynamic_bitset<>& existing_bitset = existing.bitset();
              for (size_t i = 0; i < num_samples; ++i) {
                int sum = ((int) existing_bitset[2 * i]) + ((int) existing_bitset[2 * i + 1]);
                haps_out << " " << sum;
              }
            }
          }
          else {
            haps_out << chromosome << " " << snp_prefix << ":" << count << " " << mean_position
                     << " 0 1";
            if (diploid) {
              for (char c : existing.to_bitset_string()) {
                haps_out << " " << c;
              }
            }
            else {
              for (char c : existing.to_bitset_string()) {
                // repeat the value for .haps file
                haps_out << " " << c << " " << c;
              }
            }
          }
          haps_out << endl;
          bitinfo_out << chromosome << " " << snp_prefix << ":" << count;
          bitinfo_out << " " << existing_start + offset << " " << existing_end + offset;
          if (true) {
            // also write out number of descendants
            size_t num_descendants = existing.num_values();
            bitinfo_out << " " << num_descendants;
          }
          bitinfo_out << endl;
        }
      }
    }
  }

  if (!count_only) {
    // close files
    if (file_root != "") {
      haps_out.close();
      bitinfo_out.close();
    }
  }

  cout << "max priority queue size: " << max_pq_size << endl;

  return count;
}

// writes all branches, along with time extent / position extent / MAF
void write_branches(const ARG& arg, string file_root) {
  // open haps.gz file stream
  file_utils::AutoGzOfstream haps_out_fstream;
  haps_out_fstream.openOrExit(file_root + ".haps.gz");
  // open samples file stream
  file_utils::AutoGzOfstream samples_out_fstream;
  samples_out_fstream.openOrExit(file_root + ".samples");
  size_t n = arg.leaf_ids.size();
  if (n % 2 == 1) {
    std::cerr << "ERROR: trying to write genotype data for odd number of haploid samples." << endl;
    exit(1);
  }
  samples_out_fstream << "ID_1 ID_2 missing\n0 0 0" << endl;
  for (size_t i = 1; i <= n / 2; ++i) {
    samples_out_fstream << i << "\t" << i << " 0" << endl;
  }
  // open info file stream
  file_utils::AutoGzOfstream info_out_fstream;
  info_out_fstream.openOrExit(file_root + ".branch_info");
  info_out_fstream << "# fromPos toPos fromTime toTime numDescendants MAF" << endl;

  arg_utils::visit_branches(
      arg, [&haps_out_fstream, &info_out_fstream, &arg](
               const DescendantList& desc_list, const DescendantList& _ignored,
               const ARGNode* parent, const ARGNode* child, arg_real_t start, arg_real_t end) {
        (void) _ignored;
        arg_real_t branch_volume = (parent->height - child->height) * (end - start);
        size_t num_descendants = desc_list.num_values();
        // MAF
        arg_real_t MAF = ((arg_real_t) num_descendants) / ((arg_real_t) arg.leaf_ids.size());
        int pos = start + (start + end) / 2;
        haps_out_fstream << "1 1:" << pos << " " << pos << " 1 2";
        string bitset_string = desc_list.to_bitset_string();
        std::stringstream ss;
        for (unsigned i = 0; i < bitset_string.size(); i++) {
          ss << ' ' << bitset_string[i];
        }
        haps_out_fstream << ss.str() << endl;
        info_out_fstream << start << "\t" << end << "\t" << child->height << "\t" << parent->height
                         << "\t" << num_descendants << "\t" << MAF << endl;
      });

  // close haps.gz file stream
  haps_out_fstream.close();
  // close samples file stream
  samples_out_fstream.close();
  // close branch_info file stream
  info_out_fstream.close();
}

} // namespace arg_utils
