/*
  This file is part of the ARG-Needle genealogical inference and
  analysis software suite.
  Copyright (C) 2023 ARG-Needle Developers.

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
#include "descendant_list.hpp"
#include "utils.hpp"


void arg_utils::map_genotype_to_ARG(ARG& arg, const std::vector<int>& genotype, int site_id)
{
  if (arg.roots.empty()) {
    throw std::runtime_error(THROW_LINE("Call populate_children_and_roots() first."));
  }

  const std::size_t num_leaves = arg.leaf_ids.size();

  // Jot down all the carriers of the mutation
  DescendantList carriers(num_leaves, genotype);
  const arg_real_t pos = arg.get_sites()[site_id];
  if (carriers.num_values() >= num_leaves) {
    throw std::invalid_argument(THROW_LINE(
        "Mutations carried by all samples are currently not supported. See issue #140."));
  }

  ARGEdge* edge;
  DescendantList current_carriers(num_leaves);
  while (carriers.num_values() > 0) {
    const int leaf_id = carriers.peek();
    std::tie(edge, current_carriers) = highest_carrier_edge(arg, leaf_id, carriers, pos);
    arg.add_mutation(edge, pos, -1, site_id);
    carriers.erase(current_carriers);
  }
}


std::pair<bool, std::vector<ARGEdge*>> map_genotype_to_ARG_approximate(
    ARG& arg, const std::vector<int>& genotype, arg_real_t pos, double maf_threshold)
{
  int allele_count = std::reduce(genotype.begin(), genotype.end());
  const std::size_t n = arg.leaf_ids.size();
  const double allele_frequency = static_cast<double>(allele_count) / static_cast<double>(n);

  // NB this does not consider flipped ultra rare variants (do those exist? probably, but not common)
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
      int leaf_id = carriers.peek();
      std::tie(edge, current_carriers) = arg_utils::highest_carrier_edge(arg, leaf_id, carriers, pos);
      mapped_edges.push_back(edge);
      carriers.erase(current_carriers);
    }
    // return non-flipped state and the mapped edges
    return {false, mapped_edges};
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
      std::tie(edge, current_carriers) = arg_utils::highest_carrier_edge(arg, leaf_id, carriers, pos);
      mapped_edges.push_back(edge);
      carriers.erase(current_carriers);
    }
    // return flipped state and the mapped edges
    return {true, mapped_edges};
  }
  if (allele_frequency <= maf_threshold || allele_frequency >= 1 - maf_threshold) {
    // algorithm goes like this:
    // start with a random carrier,
    // traverse the tree upwards and at each new internal node do "populate_relate_scores" downwards
    // and find the best edge only keep track of a single, best branch continue until we reach the
    // mrca of all carriers
    bool flipped = allele_frequency >= 1.0 - maf_threshold;
    DescendantList carriers(n, genotype);
    if (flipped) {
      carriers = carriers.complement();
    }
    if (carriers.num_values() == 0) {
      throw std::runtime_error(THROW_LINE("Can't date monomorphic mutations"));
    }
    ARGNode* mrca = arg_utils::most_recent_common_ancestor(arg, carriers, pos);

    arg_utils::RelateMutationMapping mutmap{};
    populate_relate_scores_lazy(arg, mrca, carriers, pos, mutmap);
    return mutmap.edge == nullptr
             ? std::pair<bool, std::vector<ARGEdge*>>(flipped, {})
             : std::pair<bool, std::vector<ARGEdge*>>(flipped, {mutmap.edge});
  }
  return {false, {}};
}


ARGNode* most_recent_common_ancestor(const ARG& arg, const DescendantList& descendants, const arg_real_t position)
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


std::pair<ARGEdge*, DescendantList> highest_carrier_edge(
    ARG& arg, const int leaf_id, const DescendantList& carriers, const arg_real_t pos)
{
  // Just call the diploid function with an empty het-list
  const DescendantList het(arg.leaf_ids.size());
  return arg_utils::highest_carrier_edge_diploid(arg, leaf_id, carriers, het, pos);
}

std::pair<ARGEdge*, DescendantList> highest_carrier_edge_diploid(
    ARG& arg, const int leaf_id, const DescendantList& homozygotes, const DescendantList& heterozygotes,
    const arg_real_t pos)
{
  const int n = arg.leaf_ids.size();
  // arg_real_t pos = arg.get_sites()[site_id];
  DescendantList current_carriers(n, leaf_id);
  const ARGNode* node = arg.arg_nodes[leaf_id].get();
  const ARGNode* parent = arg.arg_nodes[leaf_id].get();
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


DescendantList populate_relate_scores_lazy(const ARG& arg, const ARGNode* node,
                                           DescendantList& obs_carriers, const arg_real_t pos,
                                           arg_utils::RelateMutationMapping& mutmap) {

  const std::size_t n = arg.leaf_ids.size();
  const int current_node_id = node->ID;

  DescendantList descendants(n);
  DescendantList child_descendants(n);

  if (arg.leaf_ids.find(current_node_id) == arg.leaf_ids.end()) {
    for (const ARGEdge* child_edge : node->children_at(pos)) {
      child_descendants =
          arg_utils::populate_relate_scores_lazy(arg, child_edge->child, obs_carriers, pos, mutmap);
      descendants.add(child_descendants);
    }
  }
  else {
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

  if (intersect_desc_carriers + intersect_desc_noncarriers + intersect_nondesc_carriers + intersect_nondesc_noncarriers
      != n) {
    throw std::runtime_error(THROW_LINE("Something has gone wrong with the partition"));
  }

  const double node_penalty = static_cast<double>(intersect_nondesc_carriers + intersect_desc_noncarriers) / static_cast
                              <double>(n);

  const bool condition_1 = static_cast<double>(intersect_desc_carriers) / static_cast<double>(std::max(n_desc,
                               n_carriers)) > 0.7;
  const bool condition_2 = static_cast<double>(intersect_nondesc_noncarriers) / static_cast<double>(std::max(n_nondesc,
                               n_noncarriers)) > 0.7;
  const bool condition_3 = node_penalty < std::min(0.03, mutmap.penalty);

  if (condition_1 && condition_2 && condition_3) {
    mutmap.edge = node->parent_edge_at(pos);
    mutmap.penalty = node_penalty;
  }

  return descendants;
}
