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

// Visit example ARGs in various ways

#include "arg.hpp"
#include "arg_edge.hpp"
#include "arg_node.hpp"
#include "arg_utils.hpp"
#include "descendant_list.hpp"
#include "random_utils.hpp"
#include "types.hpp"
#include "utils.hpp"

#include <cassert>
#include <functional>
#include <iostream>
#include <string>
#include <tuple>

using std::cout;
using std::endl;
using std::string;
using std::tuple;
using std::unordered_set;
using std::vector;

int main(int argc, char* argv[]) {
  string directory = ARG_NEEDLE_TESTDATA_DIR "/length_5e3_samples_2e3";
  if (argc > 1) {
    directory = (string) argv[1];
  }
  if (directory[directory.size() - 1] != '/') {
    directory += "/";
  }

  if (argc > 2) {
    string s = argv[2];
    int threshold = atoi(s.c_str());
    cout << "Setting threshold to " << threshold << endl << endl;
    DescendantList::set_threshold(threshold);
  }

  // bool verbose = false;
  // if (argc > 3 && (string) argv[3] != "0") {
  //   verbose = true;
  // }

  ARG arg = arg_utils::arg_from_ts_files(directory + "nodes.txt", directory + "edges.txt");
  arg.populate_children_and_roots();
  cout << arg.arg_nodes.size() << " nodes, " << arg.get_breakpoints().size() << " trees" << endl;
  arg_utils::time_efficient_visit(arg);

  arg_real_t volume = arg_utils::total_volume(arg);
  vector<arg_real_t> afs_volume = arg_utils::allele_frequency_spectrum_volume(arg);
  arg_real_t volume_from_afs = 0.;
  for (size_t i = 0; i < afs_volume.size(); ++i) {
    volume_from_afs += afs_volume[i];
  }
  cout << endl;
  cout << "Testing visit functions: volume from ARG " << volume;
  cout << ", volume from AFS: " << volume_from_afs << endl;

  int num_branches = 0;
  int num_clades = 0;
  int num_consecutive_bitsets = 0;
  int num_unique_bitsets = 0;
  unordered_set<DescendantList, DescendantListHash> bitsets_branches, bitsets_clades;

  num_unique_bitsets = arg_utils::write_bitsets_detailed(arg, "", false, 1, "", true, true);
  num_consecutive_bitsets =
      arg_utils::write_bitsets(arg, "", false, 1, "", 0, 0, false, false, true);
  arg_utils::visit_clades(
      arg, [&num_clades, &bitsets_clades](
               DescendantList& desc_list, const ARGNode* node, arg_real_t start, arg_real_t end) {
        bitsets_clades.insert(desc_list);
        (void) desc_list;
        (void) node;
        (void) start;
        (void) end;
        ++num_clades;
      });
  arg_utils::visit_branches(
      arg, [&num_branches, &bitsets_branches](
               DescendantList& desc_list, DescendantList& child_desc_list, const ARGNode* node,
               const ARGNode* child_node, arg_real_t start, arg_real_t end) {
        bitsets_branches.insert(desc_list);
        (void) desc_list;
        (void) child_desc_list;
        (void) node;
        (void) child_node;
        (void) start;
        (void) end;
        ++num_branches;
      });
  // number of clades should be one more than number of branches
  cout << "unique bitsets from clade visit (includes root): " << bitsets_clades.size() << endl;
  cout << "unique bitsets from branches visit: " << bitsets_branches.size() << endl;

  cout << "# unique bitsets < # consecutive bitsets < # clades < # branches" << endl;
  cout << num_unique_bitsets << " " << num_consecutive_bitsets << " ";
  cout << num_clades << " " << num_branches << endl;

  return 0;
}
