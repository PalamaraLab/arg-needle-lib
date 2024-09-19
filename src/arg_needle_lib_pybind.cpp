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
#include "arg_edge.hpp"
#include "arg_node.hpp"
#include "arg_utils.hpp"
#include "descendant_list.hpp"
#include "deserialization_params.hpp"
#include "genotype_mapping.hpp"
#include "mutation.hpp"
#include "root.hpp"
#include "serialize_arg.hpp"
#include "site.hpp"
#include "types.hpp"

#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <deque>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace py = pybind11;
using std::deque;
using std::vector;
using std::tuple;
using std::deque;
using std::endl;

PYBIND11_MODULE(arg_needle_lib_pybind, m) {
  py::class_<ARGNode>(m, "ARGNode")
      .def(py::init<int, arg_real_t, arg_real_t, arg_real_t>(), py::arg("ID"), py::arg("height"),
           py::arg("start"), py::arg("end"))
      .def_readonly("ID", &ARGNode::ID)
      .def_readonly("height", &ARGNode::height)
      .def_readonly("start", &ARGNode::start)
      .def_readonly("end", &ARGNode::end)
      .def("add_parent", &ARGNode::add_parent, py::arg("start"), py::arg("end"), py::arg("parent"))
      .def("remove_parent", &ARGNode::remove_parent, py::arg("start"))
      .def("update_parent_start", &ARGNode::update_parent_start, py::arg("start_old"),
           py::arg("start_new"))
      .def("update_parent_end", &ARGNode::update_parent_end, py::arg("start"), py::arg("end_new"))
      .def(
          "parent_edges",
          [](const ARGNode& node) {
            vector<ARGEdge> edges;
            for (auto const& map_entry : node.parents) {
              edges.emplace_back(*(map_entry.second));
            }
            return edges;
          },
          "Returns a sorted list of parent edges")
      .def("parent_edge_at", &ARGNode::parent_edge_at, py::arg("position"),
           py::return_value_policy::reference)
      .def(
          "parent_starts",
          [](const ARGNode& node) {
            vector<arg_real_t> starts;
            for (auto const& map_entry : node.parents) {
              starts.push_back(map_entry.first);
            }
            return starts;
          },
          "Returns a sorted list of parent starts")
      .def("__repr__", [](const ARGNode& node) {
        std::ostringstream oss;
        oss << node;
        return oss.str();
      });

  py::class_<ARGEdge>(m, "ARGEdge")
      .def_readonly("start", &ARGEdge::start)
      .def_readonly("end", &ARGEdge::end)
      .def(
          "mutations",
          [](const ARGEdge& edge) {
            vector<Mutation> mutations;
            if (edge.mutations != nullptr) {
              for (Mutation* mutation : *edge.mutations) {
                mutations.push_back(*mutation);
              }
            }
            return mutations;
          },
          "Return list of all mutation objects")
      .def_readwrite("child", &ARGEdge::child, py::return_value_policy::reference)
      .def_readwrite("parent", &ARGEdge::parent, py::return_value_policy::reference)
      .def("__repr__", [](const ARGEdge& edge) {
        std::ostringstream oss;
        oss << edge;
        return oss.str();
      });

  py::class_<Mutation>(m, "Mutation")
      .def_readonly("edge", &Mutation::edge)
      .def_readonly("position", &Mutation::position)
      .def_readonly("height", &Mutation::height)
      .def_readonly("site_id", &Mutation::site_id);

  py::class_<Site>(m, "Site")
      .def("get_mutations", &Site::get_mutations, py::return_value_policy::reference_internal)
      .def("get_position", &Site::get_position);

  py::class_<Root>(m, "Root")
      .def_readonly("start", &Root::start)
      .def_readonly("end", &Root::end)
      .def_readwrite("node", &Root::node, py::return_value_policy::reference)
      .def("__repr__", [](const Root& root) {
        std::ostringstream oss;
        oss << "[" << root.start << ", " << root.end << ")" << endl;
        oss << *root.node << endl;
        return oss.str();
      });

  py::class_<DeserializationParams>(m, "DeserializationParams")
      .def(py::init<>())
      .def_readwrite("start", &DeserializationParams::start)
      .def_readwrite("end", &DeserializationParams::end)
      .def_readwrite("num_nodes", &DeserializationParams::num_nodes)
      .def_readwrite("offset", &DeserializationParams::offset)
      .def_readwrite("chromosome", &DeserializationParams::chromosome)
      .def_readwrite("threaded_samples", &DeserializationParams::threaded_samples)
      .def_readwrite("reserved_samples", &DeserializationParams::reserved_samples);

  py::class_<ARG>(m, "ARG")
      .def(py::init<arg_real_t, arg_real_t, int>(), "Construct an empty ARG", py::arg("start"),
           py::arg("end"), py::arg("reserved_samples") = 0)
      .def(py::init<const DeserializationParams&>(), "Construct an ARG from serialized outpout",
           py::arg("deserialization_params"))
      .def(py::init<arg_real_t, arg_real_t, vector<arg_real_t>, deque<bool>,
                    vector<std::pair<int, int>>, vector<std::pair<arg_real_t, arg_real_t>>, int>(),
           "Construct an ARG with node / edge information", py::arg("start"), py::arg("end"),
           py::arg("node_heights"), py::arg("is_sample"), py::arg("edge_ids"),
           py::arg("edge_ranges"), py::arg("reserved_samples") = -1)
      .def(py::init<arg_real_t, arg_real_t, int, int, int>(),
           "Construct a test ARG for memory footprint", py::arg("start"), py::arg("end"),
           py::arg("num_samples"), py::arg("num_trees"), py::arg("reserved_samples") = -1)
      .def_readonly("start", &ARG::start)
      .def_readonly("end", &ARG::end)
      .def_readonly("threaded_samples", &ARG::threaded_samples)
      .def_readonly("offset", &ARG::offset)                     // set using set_offset
      .def_readonly("chromosome", &ARG::chromosome)             // set using set_chromosome
      .def_readonly("reserved_samples", &ARG::reserved_samples) // set using constructors
      .def(
          "node",
          [](const ARG& arg, int ID) {
            // TODO: custom exception if not found
            return arg.arg_nodes.at(ID).get();
          },
          py::return_value_policy::reference, "Get node at ID", py::arg("ID"))
      .def(
          "node_ids", // unsorted
          [](const ARG& arg) {
            vector<int> ids;
            for (auto const& map_entry : arg.arg_nodes) {
              ids.push_back(map_entry.first);
            }
            return ids;
          },
          "Return list of all node IDs")
      .def(
          "mutations",
          [](const ARG& arg) {
            vector<Mutation> mutations;
            for (auto&& m : arg.get_mutations()) {
              mutations.push_back(*m);
            }
            return mutations;
          },
          "Return list of all Mutation objects")
      .def("num_nodes", &ARG::num_nodes)
      .def("num_edges", &ARG::num_edges)
      .def("num_mutations", &ARG::num_mutations)
      .def("get_num_sites", &ARG::get_num_sites)
      .def(
          "num_samples", [](const ARG& arg) { return arg.leaf_ids.size(); },
          "Return number of samples")
      .def_readonly("leaf_ids", &ARG::leaf_ids)         // unsorted
      .def_readonly("sample_names", &ARG::sample_names) // unsorted
      .def("set_offset", &ARG::set_offset, py::arg("offset"))
      .def("set_chromosome", &ARG::set_chromosome, py::arg("chromosome"))
      // A map from each physical position (key) to a site object containing all mutations at that position (value)
      .def("get_mutation_sites", &ARG::get_mutation_sites, py::return_value_policy::reference)
      // A sorted list of physical positions
      .def("get_site_positions", &ARG::get_site_positions, py::return_value_policy::reference)
      .def("add_sample", &ARG::add_sample, py::arg("sample_name") = "")
      .def("is_leaf", &ARG::is_leaf, py::arg("node_id"))
      .def("thread_sample", &ARG::thread_sample, py::arg("section_starts"), py::arg("sample_ids"),
           py::arg("heights"))
      .def("populate_children_and_roots", &ARG::populate_children_and_roots)
      .def("populate_mutations_on_edges", &ARG::populate_mutations_on_edges)
      .def("lowest_mutated_edge", &ARG::lowest_mutated_edge, py::return_value_policy::reference)
      .def("lowest_mutated_edge_by_site", &ARG::lowest_mutated_edge_by_site,
           py::return_value_policy::reference)
      .def("root_starts", &ARG::root_starts)
      .def("root_at", &ARG::root_at, py::return_value_policy::reference, py::arg("position"))
      .def("mrca", &ARG::mrca, py::return_value_policy::reference, py::arg("ID1"), py::arg("ID2"),
           py::arg("position"))
      .def("check_basic", &ARG::check_basic, "Basic checks after initializing an ARG",
           py::arg("stringent") = true)
      .def("check_roots", &ARG::check_roots, "Checks after populating roots")
      .def("check_children", &ARG::check_children, "Checks after populating children",
           py::arg("stringent") = true)
      .def("__str__",
           [](const ARG& arg) {
             std::ostringstream oss;
             oss << arg;
             return oss.str();
           })
      .def("deserialize_add_nodes", &ARG::deserialize_add_nodes, py::arg("node_heights"),
           py::arg("is_sample"), py::arg("node_bounds") = std::vector<std::vector<arg_real_t>>{})
      .def("deserialize_add_edges", &ARG::deserialize_add_edges)
      .def("deserialize_add_mutations", &ARG::deserialize_add_mutations)
      .def("add_mutation", &ARG::add_mutation, py::arg("edge"), py::arg("position"), py::arg("height") = -1.0,
          py::arg("site_id") = -1, py::arg("update_data_structures") = false)
      .def("get_idx_of_first_mutation_left_of", &ARG::get_idx_of_first_mutation_left_of,
           py::arg("physical_pos"), py::arg("include_equal") = false,
           py::arg("warn_out_of_range") = true)
      .def("get_idx_of_first_mutation_right_of", &ARG::get_idx_of_first_mutation_right_of,
           py::arg("physical_pos"), py::arg("include_equal") = false,
           py::arg("warn_out_of_range") = true)
      .def("get_idx_of_mutation_closest_to", &ARG::get_idx_of_mutation_closest_to,
           py::arg("physical_pos"))
      .def("clear_mutations", &ARG::clear_mutations)
      .def("clear_mutations_from_edges", &ARG::clear_mutations_from_edges);

  py::class_<DescendantList>(m, "DescendantList")
      .def_static(
          "set_threshold", &DescendantList::set_threshold, "Set threshold", py::arg("threshold"))
      .def_static("print_threshold", &DescendantList::print_threshold, "Print threshold");

  // arg_utils: general ARG querying
  m.def("arg_to_newick", &arg_utils::arg_to_newick, py::arg("arg"), py::arg("verbose") = false,
        "Return a Newick representation of an ARG, `verbose` includes branch lengths");
  m.def("trim_arg", &arg_utils::trim_arg, py::arg("arg"), py::arg("trim_start"),
        py::arg("trim_end"),
        "Trim ARG down to the specified interval. Start and end positions do not include the "
        "offset.");
  m.def("num_lineages", &arg_utils::num_lineages, py::arg("arg"), py::arg("position"),
        py::arg("height"), "Count the number of lineages at a given position and height");
  m.def("total_volume", &arg_utils::total_volume, py::arg("arg"), "Get ARG volume");
  m.def("local_volume", &arg_utils::local_volume, py::arg("arg"), py::arg("min_pos") = std::nullopt,
        py::arg("max_pos") = std::nullopt, py::arg("num_tasks") = std::nullopt, "Get the local arg volume");
  // arg_utils: association
  m.def("association_diploid_all", &arg_utils::association_diploid_all, py::arg("arg"),
        py::arg("phenotype"), py::arg("use_sample"), py::arg("file_root"),
        py::arg("chromosome") = 1, py::arg("snp_prefix") = "", py::arg("min_maf") = -1,
        py::arg("max_maf") = -1, py::arg("write_bitset_threshold") = -1,
        py::arg("calibration_factor") = 1, py::arg("concise_pvalue") = true,
        py::arg("max_only") = false, py::arg("careful") = false,
        R"pbdoc(
          ARG diploid association testing all clades.

          Writes summary statistics to [file_root].tab.gz and returns the max chi2.
          Additionally, can write significant variants in [file_root].haps.gz
          if ``write_bitset_threshold`` is passed.

          Arguments:
              arg: arg_needle_lib.ARG object with 2*n leaves
              phenotype: a length n array of phenotypes in the same order as the ARG samples.
                Missing values are specified using the use_sample array and can be filled
                with any value in the phenotype.
              use_sample: a length n array of booleans, False means the phenotype is missing
              file_root: results are written to [file_root].tab.gz and [file_root].haps.gz
              chromosome: used for writing summary statistics
              snp_prefix: prefix used for writing summary statistics
              min_maf: minimum MAF to test (default=-1 which means no minimum)
              max_maf: maximum MAF to test (default=-1 which means no maximum)
              write_bitset_threshold: variants with p-value less than the threshold are written
                to [file_root].haps.gz (default=-1, no writing)
              calibration_factor: used in ARG-MLMA workflows, can be set to 1 to instead run
                linear regression (default=1)
              concise_pvalue: return only two significant digits (default=True)
              max_only: if True, does not write to file and only returns the max chi2 (default=False)
              careful: if True, computes a residual sum of squares from the inferred beta as
                part of the standard error. If False, uses the norm of the phenotype as an
                approximation, like BOLT. The careful version will give slightly more
                significant p-values, especially for significant variants. Defaults to False.

          Returns:
              The maximum chi2 across the tests.
        )pbdoc");
  m.def("association_diploid_mutation", &arg_utils::association_diploid_mutation, py::arg("arg"),
        py::arg("phenotype"), py::arg("use_sample"), py::arg("file_root"), py::arg("mus"),
        py::arg("random_seed") = 0, py::arg("chromosome") = 1, py::arg("snp_prefix") = "",
        py::arg("min_maf") = -1, py::arg("max_maf") = -1, py::arg("write_bitset_threshold") = -1,
        py::arg("calibration_factor") = 1, py::arg("concise_pvalue") = true,
        py::arg("max_only") = false, py::arg("careful") = false,
        R"pbdoc(
          ARG diploid association testing using mutation-based sampling.

          Writes summary statistics to [file_root].tab.gz and returns the max chi2
          for each mutation rate passed. Additionally, can write significant variants
          in [file_root].haps.gz if ``write_bitset_threshold`` is passed.

          One of the columns written with this option is MU_STAR. For each mutation
          tested, the MU_STAR value says "if you instead performed testing with a smaller
          mutation rate, you would sample this mutation if mu > MU_STAR". This allows one
          to filter the tab.gz output for more stringent mutation rates.

          Arguments:
              arg: arg_needle_lib.ARG object with 2*n leaves
              phenotype: a length n array of phenotypes in the same order as the ARG samples.
                Missing values are specified using the use_sample array and can be filled
                with any value in the phenotype.
              use_sample: a length n array of booleans, False means the phenotype is missing
              file_root: results are written to [file_root].tab.gz and [file_root].haps.gz
              mus: an array of mutation rates to use for testing. In practice, for the
                [file_root].tab.gz and [file_root].haps.gz output, only the largest mutation
                rate of the mus array is used. However, for returning the maximum chi2 over
                tests, a maximum chi2 is returned for each of the mutation rates in order.
              random_seed: seed for mutation-based sampling. When set to 0, uses the time
                to seed (default=0).
              chromosome: used for writing summary statistics
              snp_prefix: prefix used for writing summary statistics
              min_maf: minimum MAF to test (default=-1 which means no minimum)
              max_maf: maximum MAF to test (default=-1 which means no maximum)
              write_bitset_threshold: variants with p-value less than the threshold are written
                to [file_root].haps.gz (default=-1, no writing)
              calibration_factor: used in ARG-MLMA workflows, can be set to 1 to instead run
                linear regression (default=1)
              concise_pvalue: return only two significant digits (default=True)
              max_only: if True, does not write to file and only returns the max chi2 (default=False)
              careful: if True, computes a residual sum of squares from the inferred beta as
                part of the standard error. If False, uses the norm of the phenotype as an
                approximation, like BOLT. The careful version will give slightly more
                significant p-values, especially for significant variants. Defaults to False.

          Returns:
              An array of maximum chi2 across the tests, one for each mutation rate in ``mus``.
        )pbdoc");
  // arg_utils: mutations
  m.def("generate_mutations", &arg_utils::generate_mutations, py::arg("arg"), py::arg("mu"),
        py::arg("random_seed") = 0, py::arg("num_mutations_hint") = 0,
        "Generate mutations with a given rate and keep on the ARG");
  m.def("generate_m_mutations", &arg_utils::generate_m_mutations, py::arg("arg"), py::arg("M"),
        py::arg("random_seed") = 0,
        "Generate exactly M mutations and keep on the ARG");
  m.def("get_mutations_matrix", &arg_utils::get_mutations_matrix, py::arg("arg"),
        py::arg("from_pos") = -std::numeric_limits<double>::infinity(),
        py::arg("to_pos") = std::numeric_limits<double>::infinity(),
        py::arg("include_left") = true, py::arg("include_right") = false,
        "Return a mutation matrix corresponding to existing mutations on the ARG");
  m.def("get_genotype", &arg_utils::get_mutation_genotype, py::arg("arg"), py::arg("mutation"),
        py::arg("diploid") = false, "Get a vector representing the genotype of the mutation");
  m.def("generate_mutations_map", &arg_utils::generate_mutations_map, py::arg("arg"), py::arg("mu"),
        py::arg("random_seed") = 0, "Generate mutations and return a map");
  m.def("write_mutations_to_haps", &arg_utils::write_mutations_to_haps, py::arg("arg"),
        py::arg("file_root"), py::arg("min_maf") = 0, py::arg("max_maf") = 1,
        py::arg("min_time") = 0, py::arg("max_time") = std::numeric_limits<double>::infinity(),
        "Write mutations to disk in haps/samples format");
  // arg_utils: GRM helper functions
  m.def("distance_matrix", &arg_utils::distance_matrix, py::arg("arg"),
        "Between-sample distance matrix in upper diagonal form");
  m.def("distance_matrix_v2", &arg_utils::distance_matrix_v2, py::arg("arg"), py::arg("alpha") = 0,
        py::arg("from_pos") = -1, py::arg("to_pos") = -1, "Between-sample distance matrix in upper diagonal form");
  m.def("compute_grm", &arg_utils::compute_grm, py::arg("arg"), py::arg("alpha") = 0,
        py::arg("batch_size") = 256, py::arg("diploid") = true, py::arg("min_maf") = 0.,
        py::arg("max_maf") = 0.5, R"pbdoc(
          Compute GRM from existing mutations.

          Arguments:
              arg: arg_needle_lib.ARG object
              alpha: if -1, the variants of branches are standardized before computing the GRM. If
                0, the variants of branches are treated as-is and not standardized. Values in between
                interpolate between these behaviors (default=-1).
              batch_size: batch size for GRM computations (default=True)
              diploid: if True, leaves 2*i and 2*i+1 are treated as haploid pairs for sample i when
                computing the GRM (default=True)
              min_maf: minimum MAF for mutations to include in the GRM (default=0)
              max_maf: maximum MAF for mutations to include in the GRM (default=0.5)

          Returns:
              A numpy matrix with the GRM (no centering operations are performed).
        )pbdoc");
  m.def("distance_matrix_maf_bins", &arg_utils::distance_matrix_maf_bins, py::arg("arg"),
        py::arg("maf_bins"),
        "Between-sample distance matrices in upper diagonal form taking in MAF bins");
  // arg_utils: compare ARGs using metrics
  m.def("tmrca_mse", &arg_utils::tmrca_mse, py::arg("arg1"), py::arg("arg2"),
        "Weighted TMRCA MSE between two ARGs");
  m.def("kc_topology", &arg_utils::kc_topology, py::arg("arg1"), py::arg("arg2"),
        "Weighted KC topology distance between two ARGs");
  m.def("metrics_stab", &arg_utils::metrics_stab, py::arg("arg1"), py::arg("arg2"),
        py::arg("num_stabs"),
        "Weighted KC topology distance and TMRCA MSE between two ARGs, sampled based on stabbing "
        "queries");
  m.def("metrics_stab_efficient", &arg_utils::metrics_stab_efficient, py::arg("arg1"),
        py::arg("arg2"), py::arg("num_stabs"), py::arg("random_kc_seed") = 0,
        py::arg("merge_type") = 0, py::arg("merge_fraction") = 0, py::arg("use_r2") = false,
        py::arg("use_log2") = false,
        "Weighted KC topology distance and TMRCA MSE between two ARGs, sampled based on stabbing "
        "queries");
  m.def("kc2_length_stab_efficient", &arg_utils::kc2_length_stab_efficient, py::arg("arg1"),
        py::arg("arg2"), py::arg("num_stabs"), py::arg("lambdas") = vector<arg_real_t>{1},
        "Weighted KC length-aware distance between two ARGs, sampled based on stabbing queries");
  m.def("kc_tmrca_vectors", &arg_utils::kc_tmrca_vectors, py::arg("arg"), py::arg("position"),
        "KC and TMRCA pairwise vectors at a position");
  m.def("bitset_overlap_full", &arg_utils::bitset_overlap_full, py::arg("arg1"), py::arg("arg2"),
        py::arg("min_position") = -1, py::arg("max_position") = -1,
        "Bitset overlap for branch precision and recall, over the full ARG or subregion");
  m.def("bitset_overlap_stab", &arg_utils::bitset_overlap_stab, py::arg("arg1"), py::arg("arg2"),
        py::arg("num_stabs"), py::arg("arg2_factor") = 1, py::arg("random_resolve_seed") = 0,
        py::arg("min_mac") = 0, py::arg("max_mac") = 0,
        "Bitset overlap for branch and mutation precision and recall, sampled based on stabbing "
        "queries");
  // arg_utils: querying bitsets in the ARG
  m.def("impute", &arg_utils::impute, py::arg("arg"), py::arg("position"), py::arg("genotypes"),
        py::arg("old") = false, "ARG-based imputation described in the ARG-Needle paper");
  m.def("mutation_match", &arg_utils::mutation_match, py::arg("arg"), py::arg("position"),
        py::arg("genotypes"),
        "Boolean for whether genotypes can arise from a single mutation at position");
  m.def("mutation_best", &arg_utils::mutation_best, py::arg("arg"), py::arg("position"),
        py::arg("genotypes"), py::arg("random_seed") = 0,
        "Hamming distance for best mutation placement on ARG at a position");
  m.def("write_bitsets_detailed", &arg_utils::write_bitsets_detailed, py::arg("arg"),
        py::arg("file_root") = "", py::arg("diploid") = false, py::arg("chromosome") = 1,
        py::arg("snp_prefix") = "", py::arg("compress") = true, py::arg("count_only") = false,
        "Write out ARG branches as bitsets");
  m.def("write_bitsets", &arg_utils::write_bitsets, py::arg("arg"), py::arg("file_root") = "",
        py::arg("diploid") = false, py::arg("chromosome") = 1, py::arg("snp_prefix") = "",
        py::arg("min_mac") = 0, py::arg("max_mac") = 0, py::arg("write_dosage") = false,
        py::arg("use_gz") = false, py::arg("count_only") = false,
        "Write out ARG branches as bitsets");
  m.def("bitset_volume_map", &arg_utils::bitset_volume_map, py::arg("arg"),
        py::arg("verbose") = false, "Traverse the ARG and return a map of bitset string to volume");
  m.def(
      "stab_return_all_bitsets",
      [](const ARG& arg, arg_real_t position) {
        vector<tuple<int, arg_real_t, deque<bool>>> new_result;
        for (auto& entry : arg_utils::stab_return_all_bitsets(arg, position)) {
          new_result.emplace_back(
              std::get<0>(entry), std::get<1>(entry), std::get<2>(entry).to_deque_bool());
        }
        return new_result;
      },
      py::arg("arg"), py::arg("position"),
      "All bitsets at a position of the ARG, each tuple with values allele count / length / "
      "bitset");
  // arg_utils: functions for testing purposes
  m.def("visit_identical", &arg_utils::visit_identical, py::arg("arg"), py::arg("rel_tol") = 1e-6,
        py::arg("abs_tol") = 0, py::arg("timing") = true, py::arg("verbose") = false,
        "Check for identical ARG visit output between slow and fast versions");
  m.def("time_efficient_visit", &arg_utils::time_efficient_visit, py::arg("arg"),
        py::arg("timing") = false, "Time efficient visit routine");

    // Functions for genotype mapping
    m.def("map_genotype_to_ARG", &arg_utils::map_genotype_to_ARG, py::arg("arg"),
          py::arg("genotype"), py::arg("pos"), "Maps a genotype to an ARG");
    m.def("map_genotypes_to_ARG", &arg_utils::map_genotypes_to_ARG, py::arg("arg"),
          py::arg("genotypes"), py::arg("positions"), py::arg("num_tasks") = std::nullopt, "Maps many genotype to an ARG");
    m.def("map_genotype_to_ARG_diploid", &arg_utils::map_genotype_to_ARG_diploid, py::arg("arg"),
          py::arg("genotype"), py::arg("site_id"), "Maps a diploid genotype to an ARG");
    m.def("map_genotype_to_ARG_approximate",
          [](ARG &arg, const std::vector<int> &genotype, arg_real_t pos) {
              auto result = arg_utils::map_genotype_to_ARG_approximate(arg, genotype, pos);
              std::vector<ARGEdge> edges;
              for (const auto edge: std::get<0>(result)) {
                  edges.push_back(*edge);
              }
              return py::make_tuple(edges, std::get<1>(result));
          }, py::arg("arg"), py::arg("genotype"), py::arg("pos"),
          "Maps a genotype to an ARG approximately, based on allele counts and frequencies.");
    m.def("most_recent_common_ancestor",
          [](ARG &arg, std::vector<int> descendants, double position) {
              if (descendants.empty()) {
                  throw std::runtime_error(THROW_LINE("Descendants list cannot be empty"));
              }
              DescendantList desc(arg.leaf_ids.size(), descendants.at(0));
              for (int i = 1; i < descendants.size(); i++) {
                  desc.set(descendants.at(i), true);
              }
              return arg_utils::most_recent_common_ancestor(arg, desc, position);
          },
          py::return_value_policy::reference, py::arg("arg"), py::arg("descendants"), py::arg("position"),
          "Finds the most recent common ancestor of a set of descendants in an ARG at a specific position.");

    // serialize_arg: ARG serialization to HDF5
    m.def("validate_serialized_arg", &arg_utils::validate_serialized_arg, py::arg("file_name"),
        "Validates the integrity of a serialized ARG file.");
    m.def("deserialize_arg", &arg_utils::deserialize_arg, py::arg("file_name"), py::arg("chunk_size") = 1000,
        py::arg("reserved_samples") = -1, "Deserialize ARG from HDF5 file.");
}
