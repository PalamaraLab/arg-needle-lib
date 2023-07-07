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

#include "arg_edge.hpp"
#include "arg_node.hpp"

#include <iostream>

using std::cout;
using std::endl;

int main(int argc, char* argv[]) {
  (void) argc;
  (void) argv;
  ARGNode a(0, 0, 0, 10);
  ARGNode b(1, 5, 0, 6);
  ARGNode c(2, 7, 3, 10);
  a.add_parent(0, 3, b);
  a.add_parent(3, 10, c);

  cout << a << endl;
  cout << b << endl;
  cout << c << endl;

  cout << "size of a (bytes): " << sizeof a << endl;
  cout << "size of b (bytes): " << sizeof b << endl;
  cout << "size of c (bytes): " << sizeof c << endl;

  cout << "Parents of a" << endl;
  for (auto const& map_entry : a.parents) {
    cout << map_entry.first << ", " << *map_entry.second << endl;
    cout << "size of parent edge (bytes): " << sizeof *map_entry.second << endl;
  }

  cout << endl << "Switching breakpoint to 5" << endl;
  a.update_parent_start(3, 5);
  a.update_parent_end(0, 5);

  cout << "Parents of a" << endl;
  for (auto const& map_entry : a.parents) {
    cout << map_entry.first << ", " << *map_entry.second << endl;
  }

  cout << endl << "Removing parent" << endl;
  a.remove_parent(0);

  cout << "Parents of a" << endl;
  for (auto const& map_entry : a.parents) {
    cout << map_entry.first << ", " << *map_entry.second << endl;
  }

  return 0;
}
