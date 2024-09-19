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

#ifndef ARG_NEEDLE_LIB_ARG_H
#define ARG_NEEDLE_LIB_ARG_H

#include "ancestor_entry.hpp"
#include "arg_edge.hpp"
#include "arg_node.hpp"
#include "deserialization_params.hpp"
#include "mutation.hpp"
#include "root.hpp"
#include "site.hpp"
#include "types.hpp"

#include <array>
#include <deque>
#include <memory>
#include <random>
#include <set>
#include <stack>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class ARG
{
private:
    int next_to_thread = -1;

    /**
     * @brief map from physical position to sites containing mutations
     */
    mutable std::map<arg_real_t, Site> mutation_sites;

    /**
     * @brief sorted vector of site positions
     */
    mutable std::vector<arg_real_t> site_positions;

    /**
      * @brief whether the mutation sites map is up to date or not
      */
    mutable bool mutation_sites_up_to_date = true;

    /**
      * @brief whether the site positions vector is up to date or not
      */
    mutable bool site_positions_up_to_date = true;

    void process_ancestor_entry(std::stack<AncestorEntry>& entries);
    void populate_children();
    void populate_roots();

    // called by check_basic()
    void check_node_heights() const;
    void check_node_spans() const;
    void check_edges() const;
    void check_single_parent_except_root_gaps(bool stringent) const;
    void check_stats() const;

    // called by check_children()
    void check_number_of_children(bool stringent) const;
    void check_parents_have_children() const;
    void check_children_have_parents() const;

    // Vector of mutations is private so that the public add_mutation can ensure the vector
    // is always sorted
    std::vector<std::unique_ptr<Mutation>> mutations;

    /**
     * @brief update the map between site positions and mutaitons.
     */
    void update_mutation_sites() const;

    /**
     * @brief update the vector of site positions.
     */
    void update_site_positions() const;

public:
    bool fragmented = false; // for now, keep this false
    arg_real_t start, end;
    int offset = 0; // used in cases where start = 0, should get set with set_offset()
    int chromosome = 1;
    int reserved_samples, threaded_samples = 0, next_general_id;
    std::unordered_map<int, std::unique_ptr<ARGNode>> arg_nodes{};
    std::unordered_map<int, std::string> sample_names;
    std::unordered_set<int> leaf_ids;
    std::map<arg_real_t, std::unique_ptr<Root>> roots;

    int num_edges_cnt = 0; // number of edges in the ARG

    // Constructor used by deserialization to create an empty ARG that data is added to a chunk at a
    // time
    explicit ARG(const DeserializationParams& dp);

    ARG(arg_real_t _start, arg_real_t _end, int _reserved_samples = 0);

    // Special constructor to quickly instantiate when reading from tskit format
    ARG(arg_real_t _start, arg_real_t _end, const std::vector<arg_real_t>& node_heights,
        const std::deque<bool>& is_sample, const std::vector<std::pair<int, int>>& edge_ids,
        const std::vector<std::pair<arg_real_t, arg_real_t>>& edge_ranges, int _reserved_samples = -1);

    // Special constructor to test memory footprint
    ARG(arg_real_t _start, arg_real_t _end, int num_samples, int num_trees,
        int _reserved_samples = -1);

    // copying and assigning is not allowed, because we have smart pointers
    ARG(ARG const&) = delete;
    ARG& operator=(ARG const&) = delete;
    // Allow move construction: https://stackoverflow.com/a/10473009
    ARG(ARG&&) = default;
    ARG& operator=(ARG&&) = default;

    void set_offset(int _offset);
    void set_chromosome(int _chromosome);
    std::vector<std::unique_ptr<Mutation>>::const_iterator next_mutation(arg_real_t pos) const;
    bool is_leaf(int node_id) const;
    int add_sample(std::string sample_name = "");
    void thread_sample(std::vector<arg_real_t> section_starts, std::vector<int> sample_ids,
                       std::vector<arg_real_t> heights);
    void populate_children_and_roots();
    void clear_mutations_from_edges();
    void populate_mutations_on_edges();
    std::vector<arg_real_t> root_starts() const;
    // Get the root that overlaps a position
    const Root* root_at(arg_real_t position) const;
    std::tuple<const ARGNode*, arg_real_t>
    mrca_nodes_with_end(const ARGNode& node1, const ARGNode& node2,
                        arg_real_t position) const; // overloading led to pybind11 issues
    const ARGNode* mrca(int ID1, int ID2, arg_real_t position) const;
    std::tuple<const ARGNode*, arg_real_t> mrca_with_end(int ID1, int ID2, arg_real_t position) const;
    ARGEdge* lowest_mutated_edge(int node_id, arg_real_t position);
    ARGEdge* lowest_mutated_edge_by_site(int node_id, int site_id);
    std::set<arg_real_t> get_breakpoints() const;
    // Outputs the nodes in sorted order
    friend std::ostream& operator<<(std::ostream& os, const ARG& arg);
    int num_nodes() const;
    int num_edges() const;
    int num_mutations() const;
    std::size_t get_num_sites() const;
    // Check these things at any time, throws exception if incorrect
    void check_basic(bool stringent = true) const;
    // Check these things after populating roots, throws exception if incorrect
    void check_roots() const;
    // Check these things after populating children, throws exception if incorrect
    void check_children(bool stringent = true) const;
    // check that the mutations are sorted (by position)
    void check_mutations_sorted() const;

    /**
     * @brief add a chunk of nodes during deserialization
     *
     */
    void deserialize_add_nodes(const std::vector<double>& node_heights,
                               const std::vector<uint8_t>& is_sample,
                               const std::vector<std::array<double, 2>>& node_bounds = {});

    /**
     * @brief add a chunk of edges during deserialization
     *
     */
    void deserialize_add_edges(const std::vector<std::array<int, 2>>& edge_ids,
                               const std::vector<std::array<double, 2>>& edge_ranges);

    /**
     * @brief add a chunk of mutations during deserialization
     *
     */
    void deserialize_add_mutations(const std::vector<arg_real_t>& positions,
                                   const std::vector<arg_real_t>& heights,
                                   const std::vector<int>& site_ids,
                                   const std::vector<std::array<int, 2>>& edge_ids);

    /**
     * @brief get map between positions and sites, updating that map if necessary
     * @return the map from physical position to sites with mutations
     */
    const std::map<arg_real_t, Site>& get_mutation_sites() const;

    /**
     * @brief get a reference to the sorted vector of site positions, updating that vector if necessary
     * @return the sorted vector of site positions
     */
    const std::vector<arg_real_t>& get_site_positions() const;

    /**
     * @brief Adds a new mutation to the ARG at a given position, height, and site ID.
     *
     * This function inserts a mutation into the ARG's list of mutations. It is inserted in a position
     * that maintains the sorted order of mutations by position. The function can also update related
     * data structures to reflect the addition of the new mutation.
     *
     * @param edge Pointer to the ARGEdge where the mutation occurred. This edge connects the mutation
     *     to the rest of the ARG.
     * @param position The physical position on the ARG where the mutation occurred. Mutations are
     *     ordered by this value.
     * @param height The height of the mutation on the ARG. This parameter is used to determine the
     *     mutation's temporal position.
     * @param site_id An integer representing the site ID associated with the mutation. This ID is
     *     used to link the mutation to a specific genomic site.
     * @param update_data_structures If true, the method updates internal data structures related to
     *     mutation sites and site positions. If false, it marks these data structures as out of date.
     *
     * @note If the `mutations` list is initially empty or if the new mutation's position is not less
     *     than any existing mutation's position, the new mutation is simply appended to the list.
     *     Otherwise, a binary search is performed to find the correct insertion point to maintain
     *     order.
     * @note If you are adding multiple mutations in quick succession, for efficiency you should set
     *     update_data_structures to false, and then call update_site_data_structures().
     */
    void add_mutation(ARGEdge* edge, arg_real_t position, arg_real_t height = -1.0, int site_id = -1,
                      bool update_data_structures = true);

    /**
     * @brief Upadte the data structurs related to mutation sites.
     *
     * This method updates the map between Positions and Sites (mutation_sites) and the vector of
     * sorted site positions (site_positions). This method is mostly an internal detail, but, may
     * be used if you are adding multiple mutations using add_mutation with parameter
     * update_data_structures set to false, after which a single call to this method will efficiently
     * update the structures.
     */
    void update_site_data_structures() const;

    /**
     * @brief reserve space in the mutation vector if you know you will be adding many mutations
     *
     */
    void reserve_n_mutations(std::size_t num_mutations);

    /**
     * @brief get a read-only view of the mutations vector
     *
     */
    const std::vector<std::unique_ptr<Mutation>>& get_mutations() const;

    /**
     * @brief clear the vector of mutations
     *
     */
    void clear_mutations();

    /**
     * @brief Get the index in the vector of mutations to the left of physical_pos.
     *
     * This is the closest mutation to physical_pos whose physical position is < physical_pos.
     * If include_equal is true, then the inequality is <=.
     *
     * @param physical_pos The physical position to compare with.
     * @param include_equal If true, the method includes mutations at the given physical position.
     * @param warn_out_of_range If true, the method issues a warning if there are no mutations to the
     *     left of physical_pos.
     * @return The index of the mutation to the left of the given physical position. If there are no
     *     mutations to the left of physical_pos, it returns the index of the first mutation (0).
     * @throws std::logic_error If there are no mutations on the ARG.
     */
    std::size_t get_idx_of_first_mutation_left_of(arg_real_t physical_pos, bool include_equal = false,
                                                  bool warn_out_of_range = true) const;

    /**
     * @brief Get the index in the vector of mutations to the right of physical_pos.
     *
     * This is the closest mutation to physical_pos whose physical position is > physical_pos.
     * If include_equal is true, then the inequality is >=.
     *
     * @param physical_pos The physical position to compare with.
     * @param include_equal If true, the method includes mutations at the given physical position.
     * @param warn_out_of_range If true, the method issues a warning if there are no mutations to the
     *     right of physical_pos.
     * @return The index of the mutation to the right of the given physical position. If there are no
     *     mutations to the right of physical_pos, it returns the index of the last mutation
     *     (#mutations-1).
     * @throws std::logic_error If there are no mutations on the ARG.
     */
    std::size_t get_idx_of_first_mutation_right_of(arg_real_t physical_pos,
                                                   bool include_equal = false,
                                                   bool warn_out_of_range = true) const;

    /**
     * @brief Get the index in the vector of mutations closest to the given physical_pos.
     *
     * @param physical_pos The physical position to compare with.
     * @return The index of the mutation closest to the given physical position.
     * @throws std::logic_error If there are no mutations on the ARG.
     */
    std::size_t get_idx_of_mutation_closest_to(arg_real_t physical_pos) const;
};

#endif // ARG_NEEDLE_LIB_ARG_H
