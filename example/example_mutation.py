# This file is part of the ARG-Needle genealogical inference and
# analysis software suite.
# Copyright (C) 2023-2024 ARG-Needle Developers.

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

if __name__ == '__main__':

    # Defining parameters for the simulation
    sample_size = 100
    length = 1e6
    seed = 4242

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

    print("Generating mutations with rate 1e-8")
    arg_needle_lib.generate_mutations(arg, mu=1e-8, random_seed=seed)
    print("# of mutations on ARG: %d" % arg.num_mutations())

    print("Generating exactly 100 mutations on the ARG")
    arg_needle_lib.generate_m_mutations(arg, M=100, random_seed=seed)
    print("# of mutations on ARG: %d" % arg.num_mutations())

    print("Generating exactly 10 mutations on the ARG and printing their properties...")
    arg_needle_lib.generate_m_mutations(arg, M=10, random_seed=seed)
    print("Parent_ID\tChild_ID\tPosition\tHeight")
    for m in arg.mutations():
        print('\t'.join([str(m.edge.parent.ID), str(m.edge.child.ID), str(m.position), str(m.height)]))

    print("Getting mutation matrix...")
    mut_mat = arg_needle_lib.get_mutations_matrix(arg)
    print("Mutations:\n", mut_mat)
    print("Shape:\n", mut_mat.shape)

    print("\nMutation bitsets")
    for i, m in enumerate(arg.mutations()):
        print(f"Haploid genotype for mutation {i}: ", *
              arg_needle_lib.get_genotype(arg, m, diploid=False))
        print(f"Diploid genotype for mutation {i}: ", *
              arg_needle_lib.get_genotype(arg, m, diploid=True))

    # Write mutations to a haps file
    arg_needle_lib.write_mutations_to_haps(arg, 'test')

    print("\nGenerating exactly 100 mutations and computing frequencies within ranges of physical position")
    num_mut = 100
    arg_needle_lib.generate_m_mutations(arg, M=num_mut, random_seed=seed)
    num_intervals = 10
    interval_length = length/num_intervals
    for i in range(num_intervals):
        from_pos = interval_length * i
        to_pos = interval_length * (i+1)
        print("Computing frequencies for mutations between positions", from_pos, "and", to_pos)
        mut_mat = arg_needle_lib.get_mutations_matrix(arg, from_pos, to_pos, include_left=True, include_right=(to_pos==length))
        print("Frequencies for this batch:", mut_mat.mean(1))

    print("\nComputing frequencies in batches of given size")
    batch_size = 12
    for i in range(0, num_mut, batch_size):
        print("Computing frequencies for the next", batch_size, "mutations, if available")
        from_pos = arg.mutations()[i].position
        to_pos = arg.mutations()[min(num_mut-1, i+batch_size-1)].position
        mut_mat = arg_needle_lib.get_mutations_matrix(arg, from_pos, to_pos, include_left=True, include_right=True)
        print("Frequencies for this batch:", mut_mat.mean(1))
