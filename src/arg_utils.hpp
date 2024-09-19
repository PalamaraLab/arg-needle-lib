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

/* Higher-level ARG utilities that use public members of ARGNode / ARG
 */

#ifndef ARG_NEEDLE_LIB_ARG_UTILS_HPP
#define ARG_NEEDLE_LIB_ARG_UTILS_HPP

#include "arg.hpp"
#include "arg_edge.hpp"
#include "arg_node.hpp"
#include "arg_traversal.hpp"
#include "descendant_list.hpp"
#include "types.hpp"
#include "utils.hpp"

#include <optional>
#include <queue>
#include <stack>
#include <string>
#include <tuple>
#include <unordered_map>


namespace arg_utils {

int num_lineages(const ARG& arg, arg_real_t position, arg_real_t height);
ARG trim_arg(ARG& arg, arg_real_t trim_start, arg_real_t trim_end);
ARG arg_from_ts_files(std::string node_file_name, std::string edge_file_name);
bool visit_identical(const ARG& arg, arg_real_t rel_tol = 1e-6, arg_real_t abs_tol = 0,
                     bool timing = true, bool verbose = false);
void time_efficient_visit(const ARG& arg, bool timing = true);
std::unordered_map<DescendantList, arg_real_t, DescendantListHash> efficient_visit(const ARG& arg);
arg_real_t total_volume(const ARG& arg);

/**
 * Computes the local volume of an ARG between min and max position.
 * NOTE: if the max_position - min_position is really large then this may result in a drop in performance relative
 * to total_volume. The main use case for this function is for testing the visit_branches function.
 *
 * @param arg the arg to calculate the local volume of
 * @param min_pos the optional start position
 * @param max_pos the optional end position
 * @param num_tasks the optional number of async tasks to use when chunking this computation
 */
arg_real_t local_volume(const ARG &arg, std::optional<arg_real_t> min_pos = std::nullopt,
                        std::optional<arg_real_t> max_pos = std::nullopt,
                        std::optional<unsigned> num_tasks = std::nullopt);


std::vector<arg_real_t> allele_frequency_spectrum_volume(const ARG& arg);
std::map<arg_real_t, std::string> generate_mutations_map(const ARG& arg, arg_real_t mu,
                                               unsigned random_seed = 0);
void generate_mutations_and_write(const ARG& arg, arg_real_t mu, std::string file_root,
                                  unsigned random_seed = 0);

/**
 *
 * @param arg the arg to generate mutations on
 * @param mu mutation rate
 * @param random_seed random seed, default 0
 * @param num_mutations_hint hint for number of mutations likely to be generated; used to reserve space in the vector
 *                           default 0. Giving a hint may improve performance if generating many mutations.
 */
void generate_mutations(ARG& arg, arg_real_t mu, unsigned random_seed = 0, std::size_t num_mutations_hint = 0);
void generate_m_mutations(ARG& arg, int m, unsigned random_seed = 0);
void write_mutations_to_haps(const ARG& arg, std::string file_root, arg_real_t min_maf = 0.,
                             arg_real_t max_maf = 1., arg_real_t min_time = 0.,
                             arg_real_t max_time = std::numeric_limits<double>::infinity());
MatXui get_mutations_matrix(const ARG& arg, arg_real_t from_pos = -std::numeric_limits<double>::infinity(),
                             arg_real_t to_pos = std::numeric_limits<double>::infinity(),
                             bool include_left = true, bool include_right = false);
std::unordered_map<std::string, arg_real_t> bitset_volume_map(const ARG& arg, bool verbose = false);
std::tuple<std::string, arg_real_t> newick_subtree(const ARGNode& node, arg_real_t position,
                                         arg_real_t dist_from_parent, bool verbose);
std::string arg_to_newick(const ARG& arg, bool verbose = false);
arg_real_t tmrca_mse(const ARG& arg1, const ARG& arg2);
arg_real_t kc_topology(const ARG& arg1, const ARG& arg2);
std::pair<arg_real_t, arg_real_t> metrics_stab(const ARG& arg1, const ARG& arg2, int num_stabs);
std::tuple<arg_real_t, arg_real_t, arg_real_t>
metrics_stab_efficient(const ARG& arg1, const ARG& arg2, int num_stabs, unsigned random_kc_seed = 0,
                       int merge_type = 0, arg_real_t merge_fraction = 0, bool use_r2 = false,
                       bool use_log2 = false);
std::vector<arg_real_t> kc2_length_stab_efficient(const ARG& arg1, const ARG& arg2, int num_stabs,
                                             std::vector<arg_real_t> lambdas = {1});
std::pair<std::vector<arg_real_t>, std::vector<arg_real_t>> kc_tmrca_vectors(const ARG& arg, arg_real_t position);
std::vector<int> fill_recurse(const ARGNode* node, int n, arg_real_t position, int depth,
                              std::vector<arg_real_t>& tmrca, std::vector<arg_real_t>& mrca_root,
                              std::vector<arg_real_t>& split_size, const std::unordered_set<int>& merge_ids,
                         bool use_log2, bool random_kc, std::mt19937& gen);
std::vector<std::vector<int>> random_binary_tree(int k, std::mt19937& gen);
std::unordered_set<int> get_merge_ids(const ARG& arg, arg_real_t position, int merge_type,
                                 arg_real_t merge_fraction, std::mt19937& gen);
std::tuple<int, int, int> bitset_overlap_full(const ARG& arg1, const ARG& arg2,
                                         arg_real_t min_position = -1,
                                         arg_real_t max_position = -1);
std::tuple<int, int, int, arg_real_t, arg_real_t, arg_real_t>
bitset_overlap_stab(const ARG& arg1, const ARG& arg2, int num_stabs, arg_real_t arg2_factor = 1,
                    unsigned random_resolve_seed = 0, int min_mac = 0, int max_mac = 0);
DescendantList fill_bitsets_recurse(
    std::unordered_map<DescendantList, std::pair<arg_real_t, arg_real_t>, DescendantListHash>& bitset_lengths,
    const ARGNode* node, int n, arg_real_t position, int index, bool random_resolve,
    std::mt19937& gen);
std::vector<std::tuple<int, arg_real_t, DescendantList>> stab_return_all_bitsets(const ARG& arg,
                                                                       arg_real_t position);
DescendantList get_bitset(const ARGNode* node, int n, arg_real_t position);
DescendantList get_carriers(const ARG& arg, const Mutation* mutation);
std::vector<int> get_mutation_genotype(const ARG& arg, const Mutation* mutation, bool diploid = false);
std::vector<arg_real_t> impute(const ARG& arg, arg_real_t position, const std::vector<int>& genotypes,
                          bool old = false);
bool mutation_match(const ARG& arg, arg_real_t position, const std::vector<int>& genotypes);
std::pair<bool, std::vector<int>> mutation_match_recurse(const ARGNode* node, arg_real_t position,
                                               size_t num_zeros, size_t num_ones,
                                               const std::vector<int>& genotypes);
int mutation_best(const ARG& arg, arg_real_t position, const std::vector<int>& genotypes,
                  unsigned random_seed = 0);
std::tuple<int, int, int> mutation_best_recurse(const ARGNode* node, arg_real_t position,
                                           int genotypes_num_zeros, const std::vector<int>& genotypes,
                                           bool random, std::mt19937& gen);
Eigen::MatrixXd compute_grm(const ARG& arg, arg_real_t alpha = 0, int batch_size = 256,
                            bool diploid = true, double min_maf = 0., double max_maf = 0.5);
std::vector<std::vector<arg_real_t>> distance_matrix(const ARG& arg); // should this be Eigen?
std::vector<std::vector<arg_real_t>> distance_matrix_v2(const ARG& arg, arg_real_t alpha = 0, arg_real_t from_pos = -1,
                                              arg_real_t to_pos = -1); // should this be Eigen?
std::vector<std::vector<std::vector<arg_real_t>>>
distance_matrix_maf_bins(const ARG& arg, std::vector<arg_real_t> maf_bins); // should this be Eigen?
arg_real_t association_diploid_all(const ARG& arg, const std::vector<arg_real_t>& raw_phenotypes,
                                   const std::deque<bool>& use_sample, std::string file_root,
                                   int chromosome = 1, std::string snp_prefix = "",
                                   arg_real_t min_maf = -1, arg_real_t max_maf = -1,
                                   arg_real_t write_bitset_threshold = -1,
                                   arg_real_t calibration_factor = 1, bool concise_pvalue = true,
                                   bool max_only = false, bool careful = false);
std::vector<arg_real_t> association_diploid_mutation(
    const ARG& arg, const std::vector<arg_real_t>& raw_phenotypes, const std::deque<bool>& use_sample,
    std::string file_root, std::vector<arg_real_t> mus, unsigned random_seed = 0, int chromosome = 1,
    std::string snp_prefix = "", arg_real_t min_maf = -1, arg_real_t max_maf = -1,
    arg_real_t write_bitset_threshold = -1, arg_real_t calibration_factor = 1,
    bool concise_pvalue = true, bool max_only = false, bool careful = false);
std::vector<std::tuple<arg_real_t, arg_real_t, arg_real_t, arg_real_t, arg_real_t, arg_real_t, int>>
association_haploid(const ARG& arg, const std::vector<arg_real_t>& raw_phenotypes, int topk = 0);
std::vector<std::tuple<arg_real_t, arg_real_t, arg_real_t>>
association_haploid_max(const ARG& arg, const std::vector<arg_real_t>& raw_phenotypes,
                        bool compress = true, arg_real_t epsilon = 1e-8);
// guaranteed to write all bitsets with no repeats
size_t write_bitsets_detailed(const ARG& arg, std::string file_root = "", bool diploid = false,
                              int chromosome = 1, std::string snp_prefix = "", bool compress = true,
                              bool count_only = false);
// guaranteed to write all bitsets, but with a possible small number of repeats
// (~1% in empirical tests)
size_t write_bitsets(const ARG& arg, std::string file_root = "", bool diploid = false,
                     int chromosome = 1, std::string snp_prefix = "", size_t min_mac = 0,
                     size_t max_mac = 0, bool write_dosage = false, bool use_gz = false,
                     bool count_only = false);
// writes all branches, which contain many bitset repeats
void write_branches(const ARG& arg, std::string file_root);



} // namespace arg_utils

#endif // ARG_NEEDLE_LIB_ARG_UTILS_HPP
