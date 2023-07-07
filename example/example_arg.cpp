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

// Build a few example ARGs and output their Newick trees

#include "arg.hpp"
#include "arg_edge.hpp"
#include "arg_node.hpp"
#include "arg_utils.hpp"
#include "types.hpp"

#include <cassert>
#include <functional>
#include <iostream>
#include <string>
#include <tuple>

using std::cout;
using std::deque;
using std::endl;
using std::string;
using std::tuple;
using std::vector;

int main(int argc, char* argv[]) {
  bool verbose = false;
  if (argc > 1 && (string) argv[1] != "0") {
    verbose = true;
  }

  ARG arg = ARG(0, 100, 3);
  arg.add_sample("first");
  arg.add_sample("second");
  arg.thread_sample(vector<arg_real_t>{0, 50}, vector<int>{0, 0}, vector<arg_real_t>{3.14, 2.718});
  arg.add_sample("third");
  arg.thread_sample(vector<arg_real_t>{0, 40}, vector<int>{1, 0}, vector<arg_real_t>{5, 3});

  cout << arg << endl;
  // Print sample names
  for (auto const& map_entry : arg.sample_names) {
    cout << map_entry.first << " : " << map_entry.second << endl;
  }
  arg.populate_children_and_roots();
  cout << "Newick ARG:" << endl;
  cout << arg_utils::arg_to_newick(arg, verbose);

  cout << endl;
  cout << "Making new ARG" << endl;
  arg = ARG(0, 100, 3);
  arg.add_sample("first");
  arg.add_sample("second");
  arg.thread_sample(vector<arg_real_t>{0}, vector<int>{0}, vector<arg_real_t>{10});
  arg.add_sample("third");
  arg.thread_sample(
      vector<arg_real_t>{0, 30, 70}, vector<int>{1, 1, 1}, vector<arg_real_t>{3, 5, 4});

  cout << arg << endl;
  arg.populate_children_and_roots();
  cout << "Newick ARG:" << endl;
  cout << arg_utils::arg_to_newick(arg, verbose);

  cout << endl;
  cout << "Making new ARG" << endl;
  arg = ARG(
      0, 100, vector<arg_real_t>{0, 0, 0, 5}, deque<bool>{true, true, true, false},
      vector<std::pair<int, int>>{std::make_pair(0, 3), std::make_pair(1, 3), std::make_pair(2, 3)},
      vector<std::pair<arg_real_t, arg_real_t>>{
          std::make_pair(0, 100), std::make_pair(0, 100), std::make_pair(0, 100)});

  cout << arg << endl;
  arg.populate_children_and_roots();
  cout << "Newick ARG:" << endl;
  cout << arg_utils::arg_to_newick(arg, verbose);

  return 0;
}
