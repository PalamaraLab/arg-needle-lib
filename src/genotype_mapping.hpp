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

#ifndef ARG_NEEDLE_LIB_GENOTYPE_MAPPING_HPP
#define ARG_NEEDLE_LIB_GENOTYPE_MAPPING_HPP

#include "arg.hpp"
#include "arg_edge.hpp"
#include "descendant_list.hpp"
#include "types.hpp"

#include <limits>
#include <optional>
#include <utility>
#include <vector>

namespace arg_utils
{

struct MutationMappingStruct {
  ARGEdge* edge = nullptr;
  double penalty = std::numeric_limits<double>::max();
};

/**
 * @brief Maps a genotype to an ancestral recombination graph (ARG).
 *
 * This function parsimoniously integrates phased genotype information into a given ARG by assigning mutations to edges
 * that only subtend carriers. It requires the ARG to be pre-populated with children and roots.
 *
 * @param arg Reference to the ARG object to be modified.
 * @param genotype A vector containing genotype information. This must have length equal to the number of leaves in the
 *                 ARG. The genotype can only take values 0 and 1.
 * @param pos The site position.
 *
 * @throws std::runtime_error if the ARG roots are empty (indicating `populate_children_and_roots()` must be called
 * first).
 * @throws std::invalid_argument if the mutation is carried by all samples, as this is currently unsupported.
 *
 * @note This function does not support mutations carried by all samples (see issue #140).
 *
 * @details The function first checks if the ARG has been properly initialized with roots.
 * For each carrier, it locates the highest edge in the ARG such that all descendants of that edge carry the mutation.
 * It then updates the ARG by placing a mutation on that edge.
 */
void map_genotype_to_ARG(ARG& arg, const std::vector<int>& genotype, arg_real_t pos);

void map_genotypes_to_ARG(ARG& arg, const MatXui& genotypes, const std::vector<arg_real_t>& positions,
    std::optional<unsigned> num_tasks = std::nullopt);

/**
 * @brief Maps a diploid genotype to an ancestral recombination graph (ARG).
 *
 * This function integrates diploid genotype information into a given ARG. It processes both homozygous and heterozygous
 * genotypes by identifying and mapping carriers of mutations. The function requires the ARG to be pre-populated with
 * children and roots.
 *
 * @param arg Reference to the ARG object to be modified.
 * @param genotype A vector containing diploid genotype information, with values 0, 1, or 2.
 * @param pos The site position.
 *
 * @throws std::runtime_error if the ARG roots are empty (indicating `populate_children_and_roots()` must be called
 * first).
 * @throws std::invalid_argument if the size of the genotype vector is not exactly half the number of ARG leaves.
 * @throws std::invalid_argument if mutations are carried by all samples, as this is currently unsupported.
 * @throws std::runtime_error if genotype values are not 0, 1, or 2.
 *
 * @note This function does not support mutations carried by all samples (see issue #140). Genotype values must be 0
 * (non-carrier), 1 (heterozygous carrier), or 2 (homozygous carrier).
 *
 * @details The function first checks if the ARG has been properly initialized with roots. It then processes the diploid
 * genotype to identify homozygous and heterozygous carriers. For each homozygous carrier, it locates the highest edge
 * in the ARG associated with the mutation and updates the ARG by adding the mutation. The function then greedily
 * processes heterozygous carriers, determining the appropriate edge for mutation placement by comparing the number of
 * carriers or the height difference between parent and child nodes.
 */
void map_genotype_to_ARG_diploid(ARG& arg, const std::vector<int>& genotype, arg_real_t pos);

/**
 * @brief Maps a genotype to an ancestral recombination graph (ARG) approximately based on allele counts.
 *
 * This function maps a phased genotype to one or more edges in the ARG based on a heuristic from Speidel et al. (2019).
 * It does not add mutations to the ARG but returns pointers to edges where mutations would be placed under the
 * heuristic. This function assumes the minor allele is ancestral. It operates differently based on allele counts
 * and minor allele frequencies (MAFs), returning a pair of a boolean (indicating if alleles are flipped) and a vector
 * of ARGEdge pointers that the genotype maps to.
 *
 * @param arg Reference to the ARG object.
 * @param genotype A vector of integers representing the genotype.
 * @param pos The genomic position of interest.
 *
 * @return A tuple containing:
 *   A vector of pointers to ARGEdges that the genotype maps to.
 *   A double representing the mapping penalty
 *
 * @details The function operates as follows:
 *   - For allele counts ≤ 4, maps up to 4 edges parsimoniously.
 *   - For allele counts ≥ n - 4, considers the complement of carriers (flipped case) and maps accordingly.
 *   - The function returns an empty vector if no successful mapping is found.
 *
 * @throws std::runtime_error if the genotype is monomorphic (all alleles are the same), as it cannot map such
 * mutations.
 */
std::tuple<std::vector<ARGEdge*>, double> map_genotype_to_ARG_approximate(
    const ARG& arg, const std::vector<int>& genotype, arg_real_t pos);

/**
 * @brief Finds the highest edge above leaf `leaf_id` such that all descendants are carriers.
 *
 * This function is a simplified version of finding the highest carrier edge in an ARG. It delegates the task to the
 * 'highest_carrier_edge_diploid' function by providing an empty list for heterozygous carriers.
 *
 * @param arg Reference to the ARG object.
 * @param leaf_id The identifier of the leaf node from which the search begins.
 * @param carriers A DescendantList object representing the homozygous carriers.
 * @param pos The position in the genome associated with the mutation.
 *
 * @return A pair consisting of a pointer to the ARGEdge representing the highest carrier edge,
 * and a DescendantList of the current carriers found up to that edge.
 *
 * @details The function calls 'highest_carrier_edge_diploid' with the provided list of carriers treated as
 * homozygous and an empty list for heterozygous carriers. This is used in cases where the distinction
 * between homozygous and heterozygous carriers is not necessary or the carriers are known to be homozygous.
 */
std::pair<ARGEdge*, DescendantList> highest_carrier_edge(
    const ARG& arg, int leaf_id, const DescendantList& carriers, arg_real_t pos);

/**
 * @brief Finds the highest carrier edge for a given leaf in a diploid ARG.
 *
 * This function traverses the ARG (Ancestral Recombination Graph) tree from a specified leaf upwards
 * to find the highest edge where the leaf's mutation carriers are located. It differentiates between
 * homozygous and heterozygous carriers.
 *
 * Homozygotes are defined as both chromosomes (id 2*l and 2*l+1) carrying a mutation, whereas
 * heterozygotes have only one of the two chromosomes carrying a mutation.
 *
 * @param arg Reference to the ARG object.
 * @param leaf_id The identifier of the leaf node from which the search begins.
 * @param homozygotes A DescendantList object representing the homozygous carriers.
 * @param heterozygotes A DescendantList object representing the heterozygous carriers.
 * @param pos The position in the genome associated with the mutation.
 *
 * @throws std::runtime_error if no siblings are found for a given node, which might indicate the presence of unary
 * nodes.
 *
 * @return A pair consisting of a pointer to the ARGEdge representing the highest carrier edge,
 * and a DescendantList of the current carriers found up to that edge.
 *
 * @details The function performs a bottom-up traversal of the ARG tree, starting from the specified leaf.
 * It combines homozygous and heterozygous carriers to identify the relevant mutation carriers.
 * Special handling is done for double-heterozygotes to ensure accurate identification.
 * The traversal stops when all carriers are found or when specific conditions are met (e.g., double-heterozygotes).
 */
std::pair<ARGEdge*, DescendantList> highest_carrier_edge_diploid(
    const ARG& arg, int leaf_id, const DescendantList& homozygotes, const DescendantList& heterozygotes, arg_real_t pos);

/**
 * @brief Finds the most recent common ancestor (MRCA) of a set of leaves in an ARG at a specific position.
 *
 * This function traverses an ancestral recombination graph (ARG) to find the MRCA of given descendants.
 * It works by iteratively moving up the ARG from the descendants, aggregating sibling nodes until all descendants
 * are accounted for, at which point the current node is the MRCA.
 *
 * @param arg A constant reference to the ARG object.
 * @param descendants A DescendantList object representing the set of descendants for which the MRCA is sought.
 * @param position The genomic position at which to find the MRCA.
 *
 * @return A pointer to the ARGNode that represents the MRCA at the specified position. Returns nullptr if the
 * list of descendants is empty.
 *
 * @throws std::runtime_error if the ARG roots are empty, indicating that `populate_children_and_roots()` must be called
 * first.
 *
 * @details The function first checks if the ARG's roots are populated. If not, an exception is thrown.
 * If the descendants list is empty, it returns nullptr. If the descendants include all leaves in the ARG,
 * the root of the ARG at the specified position is returned as the MRCA. Otherwise, it computes the MRCA by
 * traversing the ARG tree upwards from the given descendants and aggregating their sibling nodes until all
 * are accounted for.
 *
 * The process involves checking and updating descendants lists as it moves up the ARG, ensuring that the
 * node where all original descendants have been encountered is identified as the MRCA.
 */
ARGNode* most_recent_common_ancestor(const ARG& arg, const DescendantList& descendants, arg_real_t position);

/**
 * @brief Populates scores related to genetic variant carriers in an ARG, updating a mutation mapping.
 *
 * This function processes a node within an ancestral recombination graph (ARG) to evaluate and update scores
 * related to the distribution of carriers and non-carriers of a specific genetic variant (SNP) at a given position.
 * It is designed for recursive traversal of the ARG to accumulate descendant information and calculate penalties.
 *
 * @param arg A constant reference to the ARG object.
 * @param node A pointer to the ARGNode being processed.
 * @param obs_carriers A DescendantList of observed carriers of the genetic variant.
 * @param pos The genomic position at which the evaluation is performed.
 * @param mutmap A reference to a MutationMappingStruct that is updated based on the node's evaluation.
 *
 * @return A DescendantList representing the descendants of the current node.
 *
 * @details The function works by:
 *   - Identifying whether the current node is a leaf and processing accordingly.
 *   - Recursively calling itself for child nodes, if the current node is not a leaf, to build up a DescendantList.
 *   - Calculating intersections between descendants and carriers/non-carriers.
 *   - Checking the consistency of the partitioning and throwing an error if inconsistencies are found.
 *   - Calculating a node penalty based on the ratio of certain intersecting groups.
 *   - Updating the mutation mapping (mutmap) if specific conditions based on intersections and penalties are met.
 *
 * @throws std::runtime_error if the sum of intersecting and non-intersecting groups does not equal the total number of
 * individuals, indicating an error in partitioning.
 */
DescendantList populate_mutation_mapping_scores(
    const ARG& arg, const ARGNode* node, DescendantList& obs_carriers, arg_real_t pos, MutationMappingStruct& mutmap);

} // namespace arg_utils

#endif // ARG_NEEDLE_LIB_GENOTYPE_MAPPING_HPP
