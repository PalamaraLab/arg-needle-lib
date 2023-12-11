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

#include "descendant_list.hpp"
#include "utils.hpp"

void arg_utils::map_genotype_to_ARG(ARG& arg, const std::vector<int>& genotype, int site_id)
{
  if (arg.roots.empty()) {
    throw std::runtime_error(THROW_LINE("Call populate_children_and_roots() first."));
  }
  int n = arg.leaf_ids.size();

  // Jot down all the carriers of the mutation
  DescendantList carriers(n, genotype);
  arg_real_t pos = arg.get_sites()[site_id];
  if (carriers.num_values() >= n) {
    throw std::invalid_argument(THROW_LINE(
        "Mutations carried by all samples are currently not supported. See issue #140."));
  }

  ARGEdge* edge;
  DescendantList current_carriers(n);
  while (carriers.num_values() > 0) {
    // Start with a random carrier
    int leaf_id = carriers.peek();
    // DescendantList current_carriers(n, leaf_id);
    // ARGNode* node = arg.arg_nodes[leaf_id].get();
    // ARGEdge* edge = node->parent_edge_at(pos);
    // ARGNode* parent = arg.arg_nodes[leaf_id].get();
    // DescendantList leaves(n, leaf_id);
    // DescendantList sib_leaves(n);

    // // Traverse tree bottom-up
    // while (carriers.includes(sib_leaves)) {
    //   current_carriers.add(sib_leaves);
    //   node = parent;
    //   edge = node->parent_edge_at(pos); // if tree is correct this can't be null
    //   parent = edge->parent;
    //   if (current_carriers.num_values() == carriers.num_values()) {
    //     break;
    //   }

    //   // Find siblings at position
    //   ARGNode* sibling;
    //   sib_leaves = DescendantList(n);
    //   for (const auto e : parent->children_at(pos)) {
    //     if (e->child->ID != node->ID) {
    //       sibling = e->child;
    //       sib_leaves.add(arg_utils::get_bitset(sibling, arg.leaf_ids.size(), pos));
    //     }
    //   }
    //   if (sib_leaves.num_values() == 0) {
    //     throw std::runtime_error(
    //         THROW_LINE("Couldn't find any siblings, maybe you have unary nodes?"));
    //   }
    // }
    std::tie(edge, current_carriers) = highest_carrier_edge(arg, leaf_id, carriers, pos);
    arg.mutations.emplace_back(std::make_unique<Mutation>(edge, pos, -1, site_id));
    carriers.erase(current_carriers);
  }
}

