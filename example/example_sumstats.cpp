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

#include "arg.hpp"
#include "arg_utils.hpp"

#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
  // Read in the ARG from a file
  std::string directory = ARG_NEEDLE_TESTDATA_DIR "/length_1e6_samples_1e3/";
  ARG arg = arg_utils::arg_from_ts_files(directory + "nodes.txt", directory + "edges.txt");
  arg.populate_children_and_roots();

  std::cout << arg.arg_nodes.size() << " nodes, " << arg.get_breakpoints().size() << " trees" << std::endl;

  arg_utils::generate_mutations(arg, 1e-10, 18);

  std::cout << "generated " << arg.num_mutations() << " mutations" << std::endl;

  arg_utils::prepare_matmul(arg);

  Eigen::MatrixXd rand_pheno = Eigen::MatrixXd::Random(arg.leaf_ids.size()/2, 20);

  auto result = arg_utils::weighted_mut_squared_norm(arg, Eigen::MatrixXd::Random(arg.leaf_ids.size()/2, 7), true);

  return 0;
}
