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


"""Utilities for comparing `tskit.Tree` objects, used in testing.

Robinson-Foulds, KC distance, TMRCA MSE, and our r^2 metric.
"""

import numpy as np

def r2_matrix(tree_a, tree_b):
    if set(tree_a.samples()) != set(tree_b.samples()):
        raise ValueError("Trees must match in their samples")
    n = tree_a.num_samples()
    if n <= 3:
        raise ValueError("No non-singleton branches to compare")

    tree = tree_a
    list_a = []
    for node in tree.nodes():
        if node == tree.root or node == tree.children(tree.root)[1]:
            continue
        subsamples = set(tree.samples(node))
        list_a.append(subsamples)

    tree = tree_b
    list_b = []
    for node in tree.nodes():
        if node == tree.root or node == tree.children(tree.root)[1]:
            continue
        subsamples = set(tree.samples(node))
        list_b.append(subsamples)

    r2 = np.zeros((len(list_a), len(list_b)))
    for i in range(len(list_a)):
        a = list_a[i]
        for j in range(len(list_b)):
            b = list_b[j]
            both = len(a & b)
            neither = n - len(a) - len(b) + both
            same = both + neither
            diff = n - same
            r2[i][j] = max(same - diff, diff - same) / n

    return(r2)

def compare_r2(r2):
    shape = np.shape(r2)
    if len(shape) != 2:
        raise ValueError("Must start with a rectangular matrix")
    if shape[0] % 2 == 0:
        raise ValueError("Expecting odd dimensions")
    n = (shape[0] + 3) / 2
    return((np.sum(np.max(r2, axis=0)) - n) / (shape[1] - n),
           (np.sum(np.max(r2, axis=1)) - n) / (shape[0] - n))

def robinson_foulds(r2, eps=1e-8):
    shape = np.shape(r2)
    if len(shape) != 2 or shape[0] != shape[1]:
        raise ValueError("Must start with a square matrix")
    if shape[0] % 2 == 0:
        raise ValueError("Expecting odd dimensions")
    return(np.sum(np.max(r2, axis=0) < 1 - eps) + np.sum(np.max(r2, axis=1) < 1 - eps))

def kc_vectors(tree):
    n = tree.num_samples()
    graph_vector = np.zeros(int(n*(n+1)/2))
    branch_vector = np.zeros(int(n*(n+1)/2))
    samples = list(set(tree.samples()))  # sorts the entries
    root_time = tree.time(tree.root)

    index = 0
    for i in range(n):
        for j in range(i+1, n):
            mrca = tree.mrca(samples[i], samples[j])
            branch_vector[index] = root_time - tree.time(mrca)
            while mrca != tree.root:
                graph_vector[index] += 1
                mrca = tree.parent(mrca)
            index += 1

    for i in range(n):
        graph_vector[index] = 1
        branch_vector[index] = tree.time(tree.parent(samples[i]))
        index += 1

    return(graph_vector, branch_vector)

def kc_squared_distance(tree_a, tree_b, kc_lambdas=[0.0, 1.0]):
    if set(tree_a.samples()) != set(tree_b.samples()):
        raise ValueError("Trees must match in their samples")
    n = tree_a.num_samples()

    graph_vector_a, branch_vector_a = kc_vectors(tree_a)
    graph_vector_b, branch_vector_b = kc_vectors(tree_b)

    kc_dists = [0 for i in range(len(kc_lambdas))]
    for i in range(len(kc_lambdas)):
        kc_lambda = kc_lambdas[i]
        vector_a = (1-kc_lambda)*graph_vector_a + kc_lambda*branch_vector_a
        vector_b = (1-kc_lambda)*graph_vector_b + kc_lambda*branch_vector_b
        kc_dist = np.sum(np.square(vector_a - vector_b))
        kc_dists[i] = kc_dist

    return kc_dists

def tmrca_mse(tree_a, tree_b):
    if set(tree_a.samples()) != set(tree_b.samples()):
        raise ValueError("Trees must match in their samples")
    n = tree_a.num_samples()
    samples = list(set(tree_a.samples()))  # sorts the entries

    vector_a = np.zeros(int(n*(n-1)/2))
    vector_b = np.zeros(int(n*(n-1)/2))
    index = 0
    for i in range(n):
        for j in range(i+1, n):
            vector_a[index] = tree_a.tmrca(samples[i], samples[j])
            vector_b[index] = tree_b.tmrca(samples[i], samples[j])
            index += 1

    return np.mean(np.square(vector_a - vector_b))
