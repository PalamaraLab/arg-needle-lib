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

#include "genotype_mapping.hpp"

#include "arg_utils.hpp"
#include "constants.hpp"
#include "descendant_list.hpp"
#include "utils.hpp"

#include <future>

namespace
{
std::vector<std::tuple<ARGEdge*, arg_real_t>> map_genotype_to_ARG_internal(
    const ARG& arg, const std::vector<int>& genotype, const arg_real_t pos)
{
  const std::size_t num_leaves = arg.leaf_ids.size();

  std::vector<std::tuple<ARGEdge*, arg_real_t>> mutations_to_add{};

  // Jot down all the carriers of the mutation
  DescendantList carriers(num_leaves, genotype);
  if (carriers.num_values() >= num_leaves) {
    throw std::invalid_argument(
        THROW_LINE("Mutations carried by all samples are currently not supported. See issue #140."));
  }

  ARGEdge* edge;
  DescendantList current_carriers(num_leaves);
  while (carriers.num_values() > 0) {
    const int leaf_id = carriers.peek();
    std::tie(edge, current_carriers) = arg_utils::highest_carrier_edge(arg, leaf_id, carriers, pos);
    mutations_to_add.emplace_back(edge, pos);
    carriers.erase(current_carriers);
  }

  return mutations_to_add;
}
} // namespace

void arg_utils::map_genotype_to_ARG(ARG& arg, const std::vector<int>& genotype, const arg_real_t pos)
{
  if (arg.roots.empty()) {
    throw std::runtime_error(THROW_LINE("Call populate_children_and_roots() first."));
  }

  const std::vector<std::tuple<ARGEdge*, arg_real_t>> mutations_to_add = map_genotype_to_ARG_internal(arg, genotype, pos);

  for (auto& [edge, pos] : mutations_to_add) {
    arg.add_mutation(edge, pos);
  }
}

void arg_utils::map_genotypes_to_ARG(
    ARG& arg, const MatXui& genotypes, const std::vector<arg_real_t>& positions, std::optional<unsigned> num_tasks)
{
  if (arg.roots.empty()) {
    throw std::runtime_error(THROW_LINE("Call populate_children_and_roots() first."));
  }

  if (genotypes.rows() != positions.size()) {
    std::cerr << "Error: The number of genotypes does not match the number of positions." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  unsigned valid_num_tasks = num_tasks.value_or(anl::get_default_concurrency());
  if (valid_num_tasks > genotypes.rows()) {
    valid_num_tasks = genotypes.rows();
  }

  const unsigned num_rows_per_task = (genotypes.rows() + valid_num_tasks - 1) / valid_num_tasks;

  std::vector<std::future<std::vector<std::tuple<ARGEdge*, arg_real_t>>>> futures;

  // Launch asynchronous tasks
  for (int i = 0; i < valid_num_tasks; ++i) {

    const unsigned row_start = i * num_rows_per_task;
    const unsigned row_end = row_start + num_rows_per_task;

    futures.push_back(std::async(std::launch::async, [&arg, &genotypes, &positions, row_start, row_end]() {
      std::vector<std::tuple<ARGEdge*, arg_real_t>> mutations_this_block{};

      for (unsigned j = row_start; j < row_end; ++j) {

        if (j >= genotypes.rows()) {
          break;
        }

        std::vector<int> genotype_vec;
        genotype_vec.reserve(genotypes.cols());
        for (int k = 0; k < genotypes.cols(); ++k) {
          genotype_vec.push_back(static_cast<int>(genotypes(j, k)));
        }

        const arg_real_t pos = positions[j];
        auto mutations_this_row = map_genotype_to_ARG_internal(arg, genotype_vec, pos);

        for (const auto mut : mutations_this_row) {
          mutations_this_block.push_back(mut);
        }
      }

      return mutations_this_block;
    }));
  }

  std::vector<std::vector<std::tuple<ARGEdge*, arg_real_t>>> all_results;
  for (auto& future : futures) {
    all_results.push_back(future.get()); // This waits for each task to complete and collects the results
  }

  for (const auto& result : all_results) {
    for (const auto& [edge, pos] : result) {
      arg.add_mutation(edge, pos);
    }
  }
}

std::tuple<std::vector<ARGEdge*>, double> arg_utils::map_genotype_to_ARG_approximate(
    const ARG& arg, const std::vector<int>& genotype, arg_real_t pos)
{
  const int allele_count = std::reduce(genotype.begin(), genotype.end());
  const std::size_t n = arg.leaf_ids.size();

  if (allele_count <= 4) {
    DescendantList carriers(n, genotype);
    if (carriers.num_values() == 0) {
      throw std::runtime_error(THROW_LINE("Can't date monomorphic mutations"));
    }
    std::vector<ARGEdge*> mapped_edges;

    ARGEdge* edge;
    DescendantList current_carriers(n);
    while (carriers.num_values() > 0) {
      // Start with a random carrier
      const int leaf_id = carriers.peek();
      std::tie(edge, current_carriers) = highest_carrier_edge(arg, leaf_id, carriers, pos);
      mapped_edges.push_back(edge);
      carriers.erase(current_carriers);
    }
    return std::make_tuple(mapped_edges, 0.0);
  }
  if (allele_count >= n - 4) {
    // parsimoniously map complement
    DescendantList carriers(n, genotype);
    carriers = carriers.complement();
    if (carriers.num_values() == 0) {
      throw std::runtime_error(THROW_LINE("Can't date monomorphic mutations"));
    }
    std::vector<ARGEdge*> mapped_edges;

    ARGEdge* edge;
    DescendantList current_carriers(n);
    while (carriers.num_values() > 0) {
      // Start with a random carrier
      int leaf_id = carriers.peek();
      std::tie(edge, current_carriers) = highest_carrier_edge(arg, leaf_id, carriers, pos);
      mapped_edges.push_back(edge);
      carriers.erase(current_carriers);
    }
    return std::make_tuple(mapped_edges, 0.0);
  }

  // algorithm goes like this:
  // start with a random carrier,
  // traverse the tree upwards and at each new internal node do "populate_mutation_mapping_scores" downwards
  // and find the best edge only keep track of a single, best branch continue until we reach the
  // mrca of all carriers
  DescendantList carriers(n, genotype);
  if (carriers.num_values() == 0) {
    throw std::runtime_error(THROW_LINE("Can't date monomorphic mutations"));
  }
  ARGNode* mrca = most_recent_common_ancestor(arg, carriers, pos);

  MutationMappingStruct mutmap{};
  populate_mutation_mapping_scores(arg, mrca, carriers, pos, mutmap);

  return mutmap.edge == nullptr ? std::tuple(std::vector<ARGEdge*>{}, mutmap.penalty)
                                : std::tuple(std::vector<ARGEdge*>{mutmap.edge}, mutmap.penalty);
}

void arg_utils::map_genotype_to_ARG_diploid(ARG& arg, const std::vector<int>& genotype, const arg_real_t pos)
{
  if (arg.roots.empty()) {
    throw std::runtime_error(THROW_LINE("Call populate_children_and_roots() first."));
  }
  std::size_t n = arg.leaf_ids.size();

  if (genotype.size() != n / 2) {
    throw std::invalid_argument(THROW_LINE("Size of diploid genotype must be exactly half the number of ARG leaves"));
  }

  // Jot down all the heterozygotes and homozygotes
  DescendantList homozygotes(n);
  DescendantList heterozygotes(n);
  for (int i = 0; i < genotype.size(); i++) {
    if (genotype[i] == 2) {
      homozygotes.set(2 * i, true);
      homozygotes.set(2 * i + 1, true);
    } else if (genotype[i] == 1) {
      heterozygotes.set(2 * i, true);
      heterozygotes.set(2 * i + 1, true);
    } else if (genotype[i] != 0) {
      throw std::runtime_error(THROW_LINE("Diploid genotypes should be 0/1/2."));
    }
  }
  if (homozygotes.num_values() >= n) {
    throw std::invalid_argument(
        THROW_LINE("Mutations carried by all samples are currently not supported. See issue #140."));
  }

  // Process all homozygotes
  ARGEdge* edge;
  DescendantList current_carriers(n);
  while (homozygotes.num_values() > 0) {
    // Start with a random carrier
    int leaf_id = homozygotes.peek();
    std::tie(edge, current_carriers) = highest_carrier_edge_diploid(arg, leaf_id, homozygotes, heterozygotes, pos);
    arg.add_mutation(edge, pos);
    homozygotes.erase(current_carriers);

    DescendantList diploid_carriers(n);
    for (const int l : current_carriers.values()) {
      diploid_carriers.set(2 * (l / 2), true);
      diploid_carriers.set(2 * (l / 2) + 1, true);
    }
    heterozygotes.erase(diploid_carriers);
  }

  // Process all remaining heterozygotes (this will be the main thingy)
  ARGEdge* edge_0;
  ARGEdge* edge_1;
  DescendantList candidates_0(n);
  DescendantList candidates_1(n);
  while (heterozygotes.num_values() > 0) {
    int leaf_id = heterozygotes.peek();
    std::tie(edge_0, candidates_0) =
        highest_carrier_edge_diploid(arg, 2 * (leaf_id / 2), homozygotes, heterozygotes, pos);
    std::tie(edge_1, candidates_1) =
        highest_carrier_edge_diploid(arg, 2 * (leaf_id / 2) + 1, homozygotes, heterozygotes, pos);
    if (candidates_0.num_values() > candidates_1.num_values()) {
      edge = edge_0;
      current_carriers = candidates_0;
    } else if (candidates_0.num_values() < candidates_1.num_values()) {
      edge = edge_1;
      current_carriers = candidates_1;
    } else if (edge_0->parent->height - edge_0->child->height > edge_1->parent->height - edge_1->child->height) {
      edge = edge_0;
      current_carriers = candidates_0;
    } else {
      edge = edge_1;
      current_carriers = candidates_1;
    }

    DescendantList diploid_carriers(n);
    for (const int l : current_carriers.values()) {
      diploid_carriers.set(2 * (l / 2), true);
      diploid_carriers.set(2 * (l / 2) + 1, true);
    }
    heterozygotes.erase(diploid_carriers);
    arg.add_mutation(edge, pos);
  }
}

ARGNode* arg_utils::most_recent_common_ancestor(
    const ARG& arg, const DescendantList& descendants, const arg_real_t position)
{
  if (arg.roots.empty()) {
    throw std::runtime_error(THROW_LINE("Call populate_children_and_roots() first."));
  }
  if (descendants.num_values() == 0) {
    return nullptr;
  }

  const std::size_t n = arg.leaf_ids.size();
  if (descendants.num_values() == n) {
    return arg.root_at(position)->node;
  }

  DescendantList tmp_desc(descendants);
  const int leaf = descendants.peek();
  DescendantList current_descendants(n, leaf);
  tmp_desc.erase(current_descendants);
  ARGNode* current_node = arg.arg_nodes.at(leaf).get();

  while (tmp_desc.num_values() > 0) {
    const int old_id = current_node->ID;

    ARGNode* next_node = current_node->parent_edge_at(position)->parent;

    for (const ARGEdge* child_edge : next_node->children_at(position)) {
      if (const ARGNode* child_node = child_edge->child; child_node->ID != old_id) {
        current_descendants.add(arg_utils::get_bitset(child_node, n, position));
      }
    }
    tmp_desc.erase(current_descendants);
    current_node = next_node;
  }
  return current_node;
}

std::pair<ARGEdge*, DescendantList> arg_utils::highest_carrier_edge(
    const ARG& arg, const int leaf_id, const DescendantList& carriers, const arg_real_t pos)
{
  // Just call the diploid function with an empty het-list
  const DescendantList het(arg.leaf_ids.size());
  return highest_carrier_edge_diploid(arg, leaf_id, carriers, het, pos);
}

std::pair<ARGEdge*, DescendantList> arg_utils::highest_carrier_edge_diploid(const ARG& arg, const int leaf_id,
    const DescendantList& homozygotes, const DescendantList& heterozygotes, const arg_real_t pos)
{
  const int n = arg.leaf_ids.size();
  DescendantList current_carriers(n, leaf_id);

  const ARGNode* node = nullptr;
  const ARGNode* parent = nullptr;
  try {
    node = arg.arg_nodes.at(leaf_id).get();
    parent = node;
  } catch (const std::out_of_range& e) {
    std::cerr << "Error: Node with ID " << leaf_id << " not found." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  ARGEdge* edge = node->parent_edge_at(pos);
  DescendantList leaves(n, leaf_id);
  DescendantList sib_leaves(n);
  DescendantList carriers = homozygotes + heterozygotes;

  // Traverse tree bottom-up
  while (carriers.includes(sib_leaves)) {
    if (heterozygotes.num_values() > 0) {
      // We need to check for instances of double-heterozygotes, i.e.,
      // whether sib_leaves contains both chromosomes for a heterozygote,
      // in which case we terminate
      DescendantList new_carriers = sib_leaves + current_carriers;
      for (const int l : heterozygotes.values()) {
        if (l % 2 == 1) {
          continue;
        }
        if (new_carriers.get(2 * (l / 2)) && new_carriers.get(2 * (l / 2) + 1)) {
          return {edge, current_carriers};
        }
      }
    }

    // If we've found all carriers, stop
    current_carriers.add(sib_leaves);
    node = parent;
    edge = node->parent_edge_at(pos); // if tree is correct this can't be null
    parent = edge->parent;

    if (current_carriers.num_values() == carriers.num_values()) {
      break;
    }

    // Find siblings at position
    sib_leaves = DescendantList(n);
    for (const auto e : parent->children_at(pos)) {
      if (e->child->ID != node->ID) {
        const ARGNode* sibling = e->child;
        sib_leaves.add(arg_utils::get_bitset(sibling, arg.leaf_ids.size(), pos));
      }
    }
    if (sib_leaves.num_values() == 0) {
      throw std::runtime_error(THROW_LINE("Couldn't find any siblings, maybe you have unary nodes?"));
    }
  }
  return {edge, current_carriers};
}

DescendantList arg_utils::populate_mutation_mapping_scores(const ARG& arg, const ARGNode* node,
    DescendantList& obs_carriers, const arg_real_t pos, MutationMappingStruct& mutmap)
{

  const std::size_t n = arg.leaf_ids.size();
  const int current_node_id = node->ID;

  DescendantList descendants(n);
  DescendantList child_descendants(n);

  if (arg.leaf_ids.find(current_node_id) == arg.leaf_ids.end()) {
    for (const ARGEdge* child_edge : node->children_at(pos)) {
      child_descendants =
          arg_utils::populate_mutation_mapping_scores(arg, child_edge->child, obs_carriers, pos, mutmap);
      descendants.add(child_descendants);
    }
  } else {
    // is leaf
    descendants = DescendantList(n, current_node_id);
  }

  const std::size_t n_desc = descendants.num_values();
  const std::size_t n_nondesc = n - n_desc;
  const std::size_t n_carriers = obs_carriers.num_values();
  const std::size_t n_noncarriers = n - n_carriers; // obs_noncarriers.num_values();

  // check if snp maps:
  // NB we only need set sizes
  const std::size_t intersect_desc_carriers = descendants.intersect(obs_carriers).num_values();
  const std::size_t intersect_desc_noncarriers = (descendants - obs_carriers).num_values();
  const std::size_t intersect_nondesc_carriers = n_carriers - intersect_desc_carriers;
  const std::size_t intersect_nondesc_noncarriers = n_noncarriers - intersect_desc_noncarriers;

  if (intersect_desc_carriers + intersect_desc_noncarriers + intersect_nondesc_carriers +
          intersect_nondesc_noncarriers !=
      n) {
    throw std::runtime_error(THROW_LINE("Something has gone wrong with the partition"));
  }

  const double node_penalty =
      static_cast<double>(intersect_nondesc_carriers + intersect_desc_noncarriers) / static_cast<double>(n);

  const bool condition_1 =
      static_cast<double>(intersect_desc_carriers) / static_cast<double>(std::max(n_desc, n_carriers)) > 0.7;
  const bool condition_2 =
      static_cast<double>(intersect_nondesc_noncarriers) / static_cast<double>(std::max(n_nondesc, n_noncarriers)) >
      0.7;
  const bool condition_3 = node_penalty < std::min(0.03, mutmap.penalty);

  if (condition_1 && condition_2 && condition_3) {
    mutmap.edge = node->parent_edge_at(pos);
    mutmap.penalty = node_penalty;
  }

  return descendants;
}
