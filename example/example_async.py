# This file is part of the ARG-Needle genealogical inference and
# analysis software suite.
# Copyright (C) 2023 ARG-Needle Developers.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Run from this directory using

python3 example_mutation.py
"""

import msprime
import arg_needle_lib

import timeit

if __name__ == '__main__':

    # Defining parameters for the simulation
    sample_size = 10000
    length = 1e7
    seed = 1234

    print("Starting simulation")
    ts1 = msprime.simulate(
        sample_size=sample_size,
        Ne=4e4,
        length=length,
        recombination_rate=2e-8,
        random_seed=seed)

    print(str(ts1.num_trees) + " trees, " + str(ts1.num_nodes) + " nodes\n")

    # convert ts to ARG
    arg = arg_needle_lib.tskit_to_arg(ts1)
    arg.populate_children_and_roots()
    print("Done with populating children and roots\n")

    for i in range(20):
        times = timeit.repeat(lambda: arg_needle_lib.local_volume(arg, None, None, 1+i), number=1, repeat=5)
        print(f'{1+i}, {min(times)}')
