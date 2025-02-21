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

#ifndef ARG_NEEDLE_LIB_ARG_TRAVERSAL_HPP
#define ARG_NEEDLE_LIB_ARG_TRAVERSAL_HPP

#include "arg_edge.hpp"
#include "arg_node.hpp"
#include "descendant_list.hpp"
#include "mutation.hpp"
#include "types.hpp"
#include "utils.hpp"

#include <cassert>
#include <queue>
#include <stack>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace arg_utils
{

// Common function used by both visit_branches() and visit_clades().
//
// In practice, users would probably want to directly call one of
// visit_branches() or visit_clades(), but this common function helps to reduce
// code duplication.
//
//
// Idea of algorithm: we want to iterate through all the bitsets of the ARG
// and call functions on those bitsets. Therefore, we use a callback design
// pattern, where a user supplies a function that they want called at every
// bitset. The bitsets can be visited either by branch or by clade. We don't
// have objects for branches / clades. Instead we pass several arguments to
// either the branch_visit_op or clade_visit_op with the relevant information.
// The visit_branches() and visit_clades() function documentation describes the
// arguments.
//
// From each node and position, it is time-consuming to do a top-down search
// all the way down to the bottom of the ARG to get the bitset. Therefore, we
// cache relevant information about "active bitsets" in the ARG. We traverse
// the ARG from left to right, and for every position we are at we cache a
// map from active node IDs to the DescendantLists (bitsets) of those nodes
// (node_descendants_map). We delete items from the cache when we have moved
// to a new position and certain bitsets are no longer active; this is kept
// track of using a separate data structure (expiration_map) mapping
// positions to the node IDs that expire there.

// At every position, we start from the root and do a top-down postorder
// traversal which terminates at active node IDs (found using the
// node_descendants_map). This saves us from going all the way down to
// the leaves of the ARG. Then, we can go back up according to the
// postorder traversal. For each node, we merge the DescendantLists of its
// children to get its DescendantList. The minimum of the expiration
// positions of the children DescendantLists (found using the
// node_descendants_map) tells us the expiration position of this node's
// DescendantList. While we are looking up all the children, we call
// the branch_visit_op callback / functor. When we are done filling in the
// DescendantList for the node, we call the clade_visit_op callback /
// functor.
//
// In this context, visiting branches will yield the bitsets of child nodes,
// whereas visiting clades will yield the bitsets of parent nodes (as well as
// the bitsets of leaves). Branches gives more fine-grained bitsets as well as
// volume information, which is lacking when visiting clades.
//
// The min_position and max_position arguments specify the region of the ARG
// to traverse.
//
// ARG visit / traversal adapted from ARGON ARG visit algorithm
// https://github.com/pierpal/ARGON/blob/master/SRC/Argon.java#L786
//
template <typename VisitOp1, typename VisitOp2>
void visit_branches_and_clades(const ARG& arg, VisitOp1 branch_visit_op, VisitOp2 clade_visit_op,
                               arg_real_t min_position = -1, arg_real_t max_position = -1)
{
  if (arg.roots.empty()) {
    throw std::logic_error(THROW_LINE("Call populate_children_and_roots() first."));
  }
  if (min_position < arg.start) {
    min_position = arg.start;
  }
  if ((max_position < 0) || (max_position > arg.end)) {
    max_position = arg.end;
  }
  if (min_position >= max_position) {
    throw std::logic_error(THROW_LINE("Region is not well defined!"));
  }
  assert(max_position > min_position);

  int num_samples = arg.leaf_ids.size();
  std::unordered_map<int, std::pair<DescendantList, arg_real_t>> node_descendants_map;
  std::unordered_map<arg_real_t, std::vector<int>> expiration_map;
  std::stack<const ARGNode*> to_process;
  std::stack<const ARGNode*> postorder;

  // Fill in the base case for the leaf ARGNodes, calling clade_visit_op on leaves
  expiration_map.insert({max_position, std::vector<int>()});
  for (int i : arg.leaf_ids) {
    DescendantList leaf_set(num_samples, i);
    node_descendants_map.insert({i, std::make_pair(leaf_set, max_position)});
    expiration_map.at(max_position).push_back(i);
    clade_visit_op(leaf_set, arg.arg_nodes.find(i)->second.get(), min_position, max_position);
  }

  // Start at the leftmost position. With each pass through the while loop,
  // the position gets updated towards the right to be the position of the next
  // marginal tree
  arg_real_t position = min_position;
  const Root* current_root = arg.root_at(position);
  while (position < max_position) {
    // Update the root if necessary
    if (current_root->end <= position) {
      current_root = arg.root_at(position);
    }

    // Generate a postorder going up to the root using the two stacks method
    // https://www.geeksforgeeks.org/iterative-postorder-traversal/
    // The termination condition for this top-down search is nodes in the node_descendants_map
    if (node_descendants_map.find(current_root->node->ID) == node_descendants_map.end()) {
      to_process.push(current_root->node);
    }
    while (!to_process.empty()) {
      const ARGNode* node = to_process.top();
      to_process.pop();
      postorder.push(node);

      for (const ARGEdge* child_edge : node->children_at(position)) {
        const ARGNode* child_node = child_edge->child;
        if (node_descendants_map.find(child_node->ID) == node_descendants_map.end()) {
          to_process.push(child_node);
        }
      }
    }

    // Process each node in the postorder by visiting its children edges.
    // Each child will already have been filled into node_descendants_map
    arg_real_t expire_position = current_root->end;            // necessary in case we never go into while loop
    expire_position = std::min(expire_position, max_position); // make sure we don't go past end
    while (!postorder.empty()) {
      const ARGNode* node = postorder.top();
      postorder.pop();

      // We need to figure out when this node's DescendantList, which is being
      // filled in, will expire. The expiration position is the minimum of the
      // node end, the max_position, the edge end for all child edges, and the
      // child expiration position for all children.
      expire_position = node->end;
      expire_position = std::min(expire_position, max_position);
      std::vector<ARGEdge*> child_edges = node->children_at(position);
      for (const ARGEdge* child_edge : child_edges) {
        const ARGNode* child_node = child_edge->child;
        arg_real_t child_descendants_expire = node_descendants_map.at(child_node->ID).second;
        expire_position = std::min(expire_position, child_descendants_expire);
        expire_position = std::min(expire_position, child_edge->end);
      }
      // New code for handling tsinfer case
      // The branch doesn't only expire when the first child branch expires
      // It also expires if a new child edge enters
      // This code is not necessary for ARGs that consist only of binary trees,
      // but the Python test test_tsinfer_visit gives an example where it is necessary
      // TODO: profile this bit, and if it's significant, we can optionally switch
      // it off when working with binary trees.
      std::vector<ARGEdge*> child_edges_stretch = node->children_overlap(position, expire_position);
      for (const ARGEdge* child_edge_stretch : child_edges_stretch) {
        if (child_edge_stretch->start > position) {
          expire_position = std::min(expire_position, child_edge_stretch->start);
        }
      }

      // Create this ARGNode's node_descendants, while calling branch_visit_op
      DescendantList node_descendants(num_samples);
      for (const ARGEdge* child_edge : child_edges) {
        const ARGNode* child_node = child_edge->child;
        DescendantList child_descendants = node_descendants_map.at(child_node->ID).first;
        // Call the branch_visit_op callback before adding child descendants to parent
        branch_visit_op(child_descendants, node_descendants, node, child_node, position, expire_position);
        node_descendants.add(child_descendants);
      }
      node_descendants_map.insert({node->ID, std::make_pair(node_descendants, expire_position)});
      std::vector<int>& expiration_entry = expiration_map[expire_position]; // creates if not present, only hashes once
      expiration_entry.push_back(node->ID);
      // Call the clade_visit_op.
      // TODO: we may want to omit roots
      clade_visit_op(node_descendants, node, position, expire_position); // includes roots
    }
    expire_position = std::min(expire_position, current_root->end); // just to be safe
    // expire_position ends up being set to the expiration of the marginal tree

    // Increment the position and expire relevant entries from node_descendants_map
    position = expire_position;
    for (int ID : expiration_map.at(expire_position)) {
#ifdef _DEBUG
      cout << "Removing entry with end " << expire_position << " and ID " << ID << endl;
#endif // _DEBUG
      node_descendants_map.erase(node_descendants_map.find(ID));
    }
    expiration_map.erase(expiration_map.find(expire_position));
  }
  assert(expiration_map.empty());
}

// Visit all branches and call a user-supplied callback on each branch
//
// This is a templated function that takes in a callback which is called
// at each branch. Because it is templated, it needs to be in the header file
//
// A branch is a subset of an ARGEdge with no recombination events underneath.
//
// This will yield all branches, so not including the root. The branch_visit_op
// takes arguments
//
//     branch_visit_op(
//         DescendantList child_descendants,
//         DescendantList node_descendants,
//         const ARGNode* node,
//         const ARGNode* child_node,
//         real position,
//         real expire_position);
//
// Above, node->height - child_node->height can be used to get the height of
// the branch, and expire_position - position can be used to get the width of
// the branch. Additionally, child_descendants gives the list of descendants
// that would be affected by a mutation on this branch.
//
// So what is node_descendants? It should be ignored by most users. I (Brian)
// use it to compute the TMRCA and KC distance in a fancy way, and confusingly,
// it actually is filled with different values depending on the order that
// children are visited. When no children have been visited, node_descendants
// is empty. It then gets filled with child_descendants as each child gets
// visited. This allows one to record all pairs between child_descendants and
// node_descendants over the course of adding children and get all the pairwise
// relationships of leaves being joined at a node. For example:
//
//      node              .
//     / | \              .
//    c1 c2 c3            .
//   /|  /\  |\           .
//  3 1  2 4 0 5          .
//
// We would like to get all pairs of leaves that first join at node, for the
// purposes of TMRCA and KC distance. This consists of:
//     (3, 2), (3, 4), (1, 2), (1, 4)
//     (3, 0), (3, 5), (1, 0), (1, 5)
//     (2, 0), (2, 5), (4, 0), (4, 5)
//
// This is computed by calling all pairwise distances between:
//     node_descendants = {}, child_descendants = {3, 1} --> 0 pairs
//       (update node_descendants with child_descendants)
//     node_descendants = {3, 1}, child_descendants = {2, 4} --> 4 pairs
//       (update node_descendants with child_descendants)
//     node_descendants = {3, 1, 2, 4}, child_descendants = {0, 5} --> 8 pairs
//       (update node_descendants with child_descendants)
//
template <typename VisitOp>
void visit_branches(const ARG& arg, VisitOp branch_visit_op, arg_real_t min_position = -1, arg_real_t max_position = -1)
{
  visit_branches_and_clades(
      arg, branch_visit_op,
      [](DescendantList& desc_list, const ARGNode* node, arg_real_t start, arg_real_t end) {
        // no-op for clade_visit_op when visiting branches
        (void)desc_list;
        (void)node;
        (void)start;
        (void)end;
      },
      min_position, max_position);
}

// Visit all clades and call a user-supplied callback on each clade
//
// This is a templated function that takes in a callback which is called
// at each clade. Because it is templated, it needs to be in the header file
//
// A clade is a node with an extent over which there are no recombination events
// underneath. The difference between branches is that for branches, the parent
// of this node needs to remain constant, but for clades you can have the same
// constant node with recombination events above that change the parents. Thus,
// visit_clades will not yield as many duplicate bitsets.
//
// This will yield all bitsets, including the root. The clade_visit_op takes
// arguments
//
//     clade_visit_op(
//         DescendantList node_descendants,
//         const ARGNode* node,
//         real position,
//         real expire_position);
//
template <typename VisitOp>
void visit_clades(const ARG& arg, VisitOp clade_visit_op, arg_real_t min_position = -1, arg_real_t max_position = -1)
{
  visit_branches_and_clades(
      arg,
      [](DescendantList& desc_list, DescendantList& parent_desc_list, const ARGNode* parent, const ARGNode* child,
         arg_real_t start, arg_real_t end) {
        // no-op for branch_visit_op when visiting clades
        (void)desc_list;
        (void)parent_desc_list;
        (void)parent;
        (void)child;
        (void)start;
        (void)end;
      },
      clade_visit_op, min_position, max_position);
}

// Visit mutations and call a user-supplied callback on each mutation
//
// This is a templated function that takes in a callback which is called
// at each mutation. Because it is templated, it needs to be in the header file
//
// The mutation_visit_op takes
//
//     mutation_visit_op(DescendantList mutation_descendants, const Mutation* mutation);
//
// The logic here is a simplified version of that in visit_branches_and_clades().
// We use a priority queue instead of an unordered_map for keeping track of
// who has expired, which makes things O(log P) instead of O(1), where P is the
// number of active nodes. This is because we no longer iterate through every
// marginal tree start position. Also see comments in visit_branches_and_clades().
//
// (optional) mutations with age < min_time or >= max_time are skipped
// (optional) index_start specifies the index of the first mutation to be visited
// (optional) index_end specifies the index of the last mutation to be visited
//     index_end must be less than the number of mutations
//     index_end negative value means visit all mutations after index_start
//
template <typename VisitOp>
void visit_mutations(const ARG& arg, VisitOp mutation_visit_op,
    arg_real_t min_time = -std::numeric_limits<double>::infinity(),
    arg_real_t max_time = std::numeric_limits<double>::infinity(), int index_start = 0, int index_end = -1)
{
  if (arg.roots.empty()) {
    throw std::logic_error(THROW_LINE("Call populate_children_and_roots() first."));
  }
  if (arg.get_mutations().empty()) {
    std::cout << "No mutations found; did you forget to call generate_mutations_and_keep?" << std::endl;
  }
  // handle start/end conditions
  if (index_end < 0) {
    index_end = arg.num_mutations() - 1;
  }
  if (index_start > index_end) {
    throw std::logic_error(THROW_LINE("Start index must be <= end index."));
  }
  if (index_start < 0 || index_end >= arg.num_mutations()) {
    throw std::logic_error(
        THROW_LINE("Start index must at least 0 and end index must be less than the number of mutations."));
  }

  int num_samples = arg.leaf_ids.size();
  std::unordered_map<int, std::pair<DescendantList, arg_real_t>> node_descendants_map;
  std::priority_queue<std::pair<arg_real_t, int>, std::vector<std::pair<arg_real_t, int>>,
                      std::greater<std::pair<arg_real_t, int>>>
      expiration_pq;
  std::stack<const ARGNode*> to_process;
  std::stack<const ARGNode*> postorder;

  // Fill in the base case for the leaf ARGNodes
  for (int i : arg.leaf_ids) {
    DescendantList leaf_set(num_samples, i);
    node_descendants_map.insert({i, std::make_pair(leaf_set, arg.end)});
    expiration_pq.push(std::make_pair(arg.end, i));
  }

  // Iterate over all mutations, which are sorted in increasing order by position
  for (int i = index_start; i <= index_end; i++) {
    const std::unique_ptr<Mutation>& mutation_unique_ptr = arg.get_mutations()[i];
    // for (const std::unique_ptr<Mutation>& mutation_unique_ptr : arg.mutations) {
    const Mutation* mutation = mutation_unique_ptr.get();
    if (mutation->height < min_time || mutation->height >= max_time) {
      continue;
    }
    arg_real_t position = mutation->position;

    // Expire any DescendantLists ending at or before this position
    while (!expiration_pq.empty() && expiration_pq.top().first <= position) {
      node_descendants_map.erase(node_descendants_map.find(expiration_pq.top().second));
      expiration_pq.pop();
    }

    // Generate a postorder going up to this mutation using the two stacks method
    // https://www.geeksforgeeks.org/iterative-postorder-traversal/
    // The termination condition for this top-down search is nodes in the node_descendants_map
    if (node_descendants_map.find(mutation->edge->child->ID) == node_descendants_map.end()) {
      to_process.push(mutation->edge->child);
    }
    while (!to_process.empty()) {
      const ARGNode* node = to_process.top();
      to_process.pop();
      postorder.push(node);

      for (const ARGEdge* child_edge : node->children_at(position)) {
        const ARGNode* child_node = child_edge->child;
        if (node_descendants_map.find(child_node->ID) == node_descendants_map.end()) {
          to_process.push(child_node);
        }
      }
    }

    // Process each node in the postorder by visiting its children edges.
    // Each child will already have been filled into node_descendants_map
    arg_real_t expire_position;
    while (!postorder.empty()) {
      const ARGNode* node = postorder.top();
      postorder.pop();

      // We need to figure out when this node's DescendantList, which is being
      // filled in, will expire. The expiration position is the minimum of the
      // node end, the edge end for all child edges, and the child expiration
      // position for all children.
      expire_position = node->end;
      std::vector<ARGEdge*> child_edges = node->children_at(position);
      for (const ARGEdge* child_edge : child_edges) {
        const ARGNode* child_node = child_edge->child;
        arg_real_t child_descendants_expire = node_descendants_map.at(child_node->ID).second;
        expire_position = std::min(expire_position, child_descendants_expire);
        expire_position = std::min(expire_position, child_edge->end);
      }
      // New code for handling tsinfer case
      // The branch doesn't only expire when the first child branch expires
      // It also expires if a new child edge enters
      // This code is not necessary for ARGs that consist only of binary trees,
      // but the Python test test_tsinfer_visit gives an example where it is necessary
      // TODO: profile this bit, and if it's significant, we can optionally switch
      // it off when working with binary trees.
      std::vector<ARGEdge*> child_edges_stretch = node->children_overlap(position, expire_position);
      for (const ARGEdge* child_edge_stretch : child_edges_stretch) {
        if (child_edge_stretch->start > position) {
          expire_position = std::min(expire_position, child_edge_stretch->start);
        }
      }

      // Create this ARGNode's node_descendants
      DescendantList node_descendants(num_samples);
      for (const ARGEdge* child_edge : child_edges) {
        const ARGNode* child_node = child_edge->child;
        DescendantList child_descendants = node_descendants_map.at(child_node->ID).first;
        node_descendants.add(child_descendants);
      }
      node_descendants_map.insert({node->ID, std::make_pair(node_descendants, expire_position)});
      expiration_pq.push(std::make_pair(expire_position, node->ID));
    }

    // Get the Mutation's descendant list and call the visit op
    DescendantList mutation_descendants = node_descendants_map.at(mutation->edge->child->ID).first;
    mutation_visit_op(mutation_descendants, mutation);
  }
}

} // namespace arg_utils

#endif // ARG_NEEDLE_LIB_ARG_TRAVERSAL_HPP
