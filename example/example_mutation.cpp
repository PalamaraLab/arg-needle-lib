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

// Test out performance of volume functions and mutation generation functions...

#include "arg.hpp"
#include "arg_edge.hpp"
#include "arg_node.hpp"
#include "arg_utils.hpp"
#include "descendant_list.hpp"
#include "random_utils.hpp"
#include "types.hpp"
#include "utils.hpp"

#include <cassert>
#include <chrono>
#include <functional>
#include <iostream>
#include <string>
#include <tuple>

using std::cout;
using std::endl;
using std::string;
using std::tuple;
using Clock = std::chrono::high_resolution_clock;

// Testing out various ways to generate mutations
int main(int argc, char* argv[]) {
  string directory = ARG_NEEDLE_TESTDATA_DIR "/length_1e6_samples_1e3";
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

  // Read in the ARG from a file
  ARG arg = arg_utils::arg_from_ts_files(directory + "nodes.txt", directory + "edges.txt");
  arg.populate_children_and_roots();
  cout << arg.arg_nodes.size() << " nodes, " << arg.get_breakpoints().size() << " trees" << endl;

  // 1. Comparing functions to obtain the volume
  // 1a.  volume function
  std::chrono::time_point<Clock> last_time, curr_time;
  last_time = Clock::now();
  arg_real_t tot_vol1 = arg_utils::total_volume(arg);
  curr_time = Clock::now();
  arg_real_t split_time =
      std::chrono::duration_cast<std::chrono::milliseconds>(curr_time - last_time).count();
  cout << "Volume: " << tot_vol1 << " (Finished volume computation in " << (split_time / 1000)
       << " seconds.)" << endl;

  // 0. Generating mutations in a random order (keeping them on the ARG)
  last_time = Clock::now();
  arg_utils::generate_mutations(arg, 1e-7, 4242);
  curr_time = Clock::now();
  split_time = std::chrono::duration_cast<std::chrono::milliseconds>(curr_time - last_time).count();
  cout << "(Simulated mutations and kept in " << (split_time / 1000) << " seconds.)" << endl;

  // Experiment 1: Searching for variants within 10 kb of midway point
  arg_real_t pos = arg.start + (arg.end - arg.start) / 2;
  last_time = Clock::now();
  // Use the iterator searching here
  for (auto it = arg.next_mutation(pos - 1e4); it != arg.next_mutation(pos + 1e4); ++it) {
    assert(((*it)->position >= pos - 1e4) && ((*it)->position <= pos + 1e4));
  }
  curr_time = Clock::now();
  split_time = std::chrono::duration_cast<std::chrono::milliseconds>(curr_time - last_time).count();
  cout << "Finished 10kb search in " << split_time << " milliseconds.)" << endl;

  // Experiment 2: Searching for variants within 100 kb of midway point
  last_time = Clock::now();
  for (auto it = arg.next_mutation(pos - 1e5); it != arg.next_mutation(pos + 1e5); ++it) {
    assert(((*it)->position >= pos - 1e5) && ((*it)->position <= pos + 1e5));
  }
  curr_time = Clock::now();
  split_time = std::chrono::duration_cast<std::chrono::milliseconds>(curr_time - last_time).count();
  cout << "Finished 100kb search in " << split_time << " milliseconds.)" << endl;

  // New output function to generate mutations and write onto the ARG
  cout << "Setting Mutation rate on O(1e-8)" << endl;
  arg_utils::generate_mutations(arg, 1.2e-8, 4242);
  cout << "Starting Mutation Generation v2..." << endl;
  last_time = Clock::now();
  arg_utils::write_mutations_to_haps(arg, "foo1");
  curr_time = Clock::now();
  split_time = std::chrono::duration_cast<std::chrono::milliseconds>(curr_time - last_time).count();
  cout << "Finished writing mut gen v2 (in " << (split_time / 1000) << " seconds.)" << endl << endl;

  // Previous function to generate and write out mutations as bitsets
  cout << "Starting Mutation Generation v1..." << endl;
  last_time = Clock::now();
  arg_utils::generate_mutations_and_write(arg, 1.2e-8, "foo2", 4242);
  curr_time = Clock::now();
  split_time = std::chrono::duration_cast<std::chrono::milliseconds>(curr_time - last_time).count();
  cout << "Finished writing mut gen v1 (in " << (split_time / 1000) << " seconds.)" << endl << endl;

  return 0;
}
