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

"""Metrics for evaluating ARGs."""

import logging
import math
import numpy as np
import arg_needle_lib_pybind

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')


def kc2_tmrca_mse_stab(arg1, arg2, num_stabs):
    """
    Same as kc2_tmrca_sv_stab, this name is kept for backwards compatibility.
    """
    return kc2_tmrca_sv_stab(arg1, arg2, num_stabs)


def kc2_tmrca_sv_stab(arg1, arg2, num_stabs):
    """
    Compares two ARGs using stabbing queries, returning topology KC squared, TMRCA MSE,
    and the split-vector (SV) metric squared and N-normalised.

    The SV metric is computed as the L2 distance between the split-size vectors of arg1 and arg2.
    The split-size vector enumerates num_leaves(mrca(i, j)) for all leaf pairs {i, j}.
    This metric is due to Smith, https://doi.org/10.1093/sysbio/syab100.

    Arguments:
        arg1: the first ARG
        arg2: the second ARG
        num_stabs: how many stabbing queries to perform. For i from 0 to num_stabs - 1 inclusive,
          the proportional positions along the genome are then sampled as i*phi (mod 1), where phi
          is the golden ratio. See also Supplementary Note 2 in
          https://doi.org/10.1038/s41588-023-01379-x

    Returns:
        Dictionary with keys "kc2", "tmrca_mse", "sv2" and float values.
    """
    assert num_stabs > 0
    kc2, tmrca_mse, sv2 = arg_needle_lib_pybind.metrics_stab_efficient(arg1, arg2, num_stabs)
    return {"kc2": kc2, "tmrca_mse": tmrca_mse, "sv2": sv2}


def kc2_tmrca_mse_stab_average(args1, args2, num_stabs):
    """Compares two lists of ARGs with averaging, returning topology KC squared and TMRCA MSE.

    Averaging is useful for the case of running ARG-Needle with different orders of threading.

    Arguments:
        args1: the first list of ARGs
        args2: the second list of ARGs
        num_stabs: how many stabbing queries to perform. For i from 0 to num_stabs - 1 inclusive,
          the proportional positions along the genome are then sampled as i*phi (mod 1), where phi
          is the golden ratio. See also Supplementary Note 2 in
          https://doi.org/10.1038/s41588-023-01379-x

    Returns:
        Dictionary with keys "kc2" and "tmrca_mse" and float values.
    """
    if not isinstance(args1, list):
        args1 = [args1]
    if not isinstance(args2, list):
        args2 = [args2]
    assert num_stabs > 0

    span = args1[0].end - args1[0].start
    phi = (math.sqrt(5) - 1) / 2
    proportion = 0
    stabs = []
    for i in range(num_stabs):
        proportion = (proportion + phi) % 1
        stabs.append(proportion * span + args1[0].start)

    kc2 = 0
    tmrca_mse = 0
    n = args1[0].num_samples()
    num_pairs = n * (n - 1) // 2

    for i, position in enumerate(stabs):
        if i % 1000 == 0:
            logging.info("Metrics stab {}".format(i))

        mean_kc_vector1 = np.zeros(num_pairs)
        mean_tmrca_vector1 = np.zeros(num_pairs)
        mean_kc_vector2 = np.zeros(num_pairs)
        mean_tmrca_vector2 = np.zeros(num_pairs)

        for arg in args1:
            kc_vector, tmrca_vector = arg_needle_lib_pybind.kc_tmrca_vectors(arg, position)
            mean_kc_vector1 += kc_vector
            mean_tmrca_vector1 += tmrca_vector

        for arg in args2:
            kc_vector, tmrca_vector = arg_needle_lib_pybind.kc_tmrca_vectors(arg, position)
            mean_kc_vector2 += kc_vector
            mean_tmrca_vector2 += tmrca_vector

        mean_kc_vector1 /= len(args1)
        mean_tmrca_vector1 /= len(args1)
        mean_kc_vector2 /= len(args2)
        mean_tmrca_vector2 /= len(args2)

        kc2 += np.sum(np.square(mean_kc_vector1 - mean_kc_vector2))
        tmrca_mse += np.sum(np.square(mean_tmrca_vector1 - mean_tmrca_vector2))

    kc2 /= num_stabs
    tmrca_mse /= num_stabs
    tmrca_mse /= num_pairs
    return {"kc2": kc2, "tmrca_mse": tmrca_mse}


def kc2_length_stab(arg1, arg2, num_stabs, lambdas=[1]):
    """Compares two ARGs using stabbing queries, returning length-aware KC squared.

    Arguments:
        arg1: the first ARG
        arg2: the second ARG
        num_stabs: how many stabbing queries to perform. For i from 0 to num_stabs - 1 inclusive,
          the proportional positions along the genome are then sampled as i*phi (mod 1), where phi
          is the golden ratio. See also Supplementary Note 2 in
          https://doi.org/10.1038/s41588-023-01379-x
        lambdas: array-like containing values of lambdas

    Returns:
        List of values containing length-aware KC squared for each lambda.
    """
    assert num_stabs > 0
    kc2_lambdas = arg_needle_lib_pybind.kc2_length_stab_efficient(arg1, arg2, num_stabs, lambdas)
    return kc2_lambdas


def kc2_special_stab(arg, arg_true, num_stabs, special_behavior=None, merge_fraction=0, random_kc_seed=0):
    """Compares two ARGs using stabbing queries, returning topology KC squared under special behavior.

    Arguments:
        arg1: the first ARG, which can be perturbed by either merging nodes or randomly resolving
          polytomies.
        arg2: the second ARG, which can be perturbed by either merging nodes or randomly resolving
          polytomies.
        num_stabs: how many stabbing queries to perform. For i from 0 to num_stabs - 1 inclusive,
          the proportional positions along the genome are then sampled as i*phi (mod 1), where phi
          is the golden ratio. See also Supplementary Note 2 in
          https://doi.org/10.1038/s41588-023-01379-x.
        special_behavior: can be "random_resolve", "heuristic_merge", or "random_merge". Each of
          these is described in Supplementary Note 2 of https://doi.org/10.1038/s41588-023-01379-x
        merge_fraction: used by heuristic_merge and random_merge to figure out how many nodes
          to merge in each stabbing query tree. Default is 0.
        random_kc_seed: used by random_resolve and random_merge to perform random behavior. The
          C++ RNG is seeded with this seed, and then iterates through the different stabbing query
          positions, consuming randomness from the RNG as necessary. Therefore different stabbing
          query trees will have different random behavior, because the RNG is only seeded once.
          Default is 0, but for random_resolve and random_merge, we check that the seed is nonzero
          before continuing so that the user is aware of the importance of seeding.

    Returns:
        The topology KC squared value, under the appropriate special behavior.
    """
    assert num_stabs > 0
    if not (merge_fraction >= 0 and merge_fraction <= 1):
        raise ValueError("Merge fraction value of {} is out of bounds".format(merge_fraction))
    if special_behavior not in ["random_resolve", "heuristic_merge", "random_merge"]:
        raise ValueError("Special KC behavior not recognized")

    if special_behavior == "random_merge":
        if random_kc_seed == 0:
            raise ValueError("For clarity, please run with a nonzero seed")
        result = arg_needle_lib_pybind.metrics_stab_efficient(
            arg, arg_true, num_stabs, random_kc_seed=random_kc_seed,
            merge_type=1, merge_fraction=merge_fraction)

    elif special_behavior == "heuristic_merge":
        if random_kc_seed != 0:
            logging.info("random_kc_seed of {} will be ignored".format(random_kc_seed))
        result = arg_needle_lib_pybind.metrics_stab_efficient(
            arg, arg_true, num_stabs, merge_type=2, merge_fraction=merge_fraction)

    else:
        assert special_behavior == "random_resolve"
        if random_kc_seed == 0:
            raise ValueError("Sorry, cannot perform this behavior with a seed of 0")
        if merge_fraction > 0:
            logging.info("merge_fraction value of {} will be ignored".format(merge_fraction))
        result = arg_needle_lib_pybind.metrics_stab_efficient(
            arg, arg_true, num_stabs, random_kc_seed=random_kc_seed)

    # metrics_stab_efficient returns both KC squared and TMRCA MSE; we only want KC squared
    return result[0]


def rf_total_variation_stab(arg1, arg2, num_stabs, min_mac=0, max_mac=0, may_contain_polytomies=False):
    """Computes scaled Robinson-Foulds and total variation using stabbing queries.

    In this version, polytomies are not resolved, which may bias the metric. See also rf_resolve_stab.

    Arguments:
        arg1: the first ARG
        arg2: the second ARG
        num_stabs: how many stabbing queries to perform. For i from 0 to num_stabs - 1 inclusive,
          the proportional positions along the genome are then sampled as i*phi (mod 1), where phi
          is the golden ratio. See also Supplementary Note 2 in
          https://doi.org/10.1038/s41588-023-01379-x.
        min_mac: min MAC (inclusive)
        max_mac: max MAC (exclusive, value of 0 means set to maximum)
        may_contain_polytomies: if False, and the two ARGs do not have polytomies and there are no
          MAC filters, then a simple check is done verifying that calculations are done correctly.
          If there may be polytomies, then set to True. See also rf_resolve_stab which will randomly
          resolve the polytomies. Default = False.

    Returns:
        Dictionary with keys "scaled_robinson_foulds" and "total_variation" and float values, as
          well as the less used keys "branch_recall", "branch_precision", "mutation_recall", and
          "mutation_precision".
    """
    num1, num2, num_common, length1, length2, length_common = arg_needle_lib_pybind.bitset_overlap_stab(
        arg1, arg2, num_stabs, min_mac=min_mac, max_mac=max_mac)

    n = arg1.num_samples()
    if min_mac == 0 and max_mac == 0 and not may_contain_polytomies:
        assert num1 == num_stabs * (2*n - 2)
        assert num2 == num_stabs * (2*n - 2)
    # scale Robinson-Foulds to be between 0 and 1, similar to Relate
    rf = num1 + num2 - 2*num_common
    rf_denominator = num1 + num2 - 2*num_stabs*n  # remove singletons from the denominator

    result = {
        "scaled_robinson_foulds": rf * 1.0 / rf_denominator,
        "branch_recall": num_common * 1.0 / num1,
        "branch_precision": num_common * 1.0 / num2,
        "mutation_recall": length_common / length1,
        "mutation_precision": length_common / length2
    }

    arg2_factor = length1 / length2
    num1, num2, num_common, length1, length2, length_common = arg_needle_lib_pybind.bitset_overlap_stab(
        arg1, arg2, num_stabs, arg2_factor, min_mac=min_mac, max_mac=max_mac)
    logging.info("length1 = {}, length2 = {}".format(length1, length2))
    assert math.isclose(length1, length2, rel_tol=1e-4)
    result["total_variation"] = 1 - length_common / length1

    # Return the dict containing all the answers
    return result


def rf_resolve_stab(arg1, arg2, num_stabs, resolve_seeds, min_mac=0, max_mac=0):
    """Computes scaled Robinson-Foulds using stabbing queries, with polytomies randomly resolved.

    Arguments:
        arg1: the first ARG
        arg2: the second ARG
        num_stabs: how many stabbing queries to perform. For i from 0 to num_stabs - 1 inclusive,
          the proportional positions along the genome are then sampled as i*phi (mod 1), where phi
          is the golden ratio. See also Supplementary Note 2 in
          https://doi.org/10.1038/s41588-023-01379-x.
        resolve_seeds: a list of seeds to use for randomly resolving polytomies.
        min_mac: min MAC (inclusive)
        max_mac: max MAC (exclusive, value of 0 means set to maximum)

    Returns:
        Dictionary with key "scaled_robinson_foulds" and a float value, as well as the less
          used keys "branch_recall" and "branch_precision".
    """
    random_resolve_values = []
    for resolve_seed in resolve_seeds:
        random_resolve_values.append(
            arg_needle_lib_pybind.bitset_overlap_stab(
                arg1, arg2, num_stabs, random_resolve_seed=resolve_seed, min_mac=min_mac, max_mac=max_mac))
    sum_values = []
    for resolve_value_index in range(len(random_resolve_values[0])):
        sum_values.append(sum([x[resolve_value_index] for x in random_resolve_values]))

    num1, num2, num_common, length1, length2, length_common = sum_values

    n = arg1.num_samples()
    if min_mac == 0 and max_mac == 0:
        # This check is always true because we have resolved polytomies
        assert num1 == num_stabs * len(resolve_seeds) * (2*n - 2)
        assert num2 == num_stabs * len(resolve_seeds) * (2*n - 2)
    # scale Robinson-Foulds to be between 0 and 1, similar to Relate
    rf = num1 + num2 - 2*num_common
    rf_denominator = num1 + num2 - 2*num_stabs*len(resolve_seeds)*n  # remove singletons from the denominator

    result = {
        "scaled_robinson_foulds": rf * 1.0 / rf_denominator,
        "branch_recall": num_common * 1.0 / num1,
        "branch_precision": num_common * 1.0 / num2,
    }

    # Return the dict containing all the answers
    return result


def mutation_match_all(arg, ts_with_mutations, resolve_reps=0):
    """Tries to best match a given set of mutations on to an ARG, returning results.

    Args:
        arg: `arg_needle_lib.ARG` object (usually one that has been inferred)
        ts_with_mutations: `tskit.TreeSequence` object containing mutations, that we want to
          check for on `arg`
        resolve_reps (int): how many times to randomly resolve polytomies in the `arg` object
          when matching. This is unnecessary if `arg` contains no polytomies, but otherwise
          allows for matching against resolved binary trees, averaging over several random
          resolutions. Default 0, which corresponds to no resolving of polytomies.

    Returns:
        Two lists (errors, diffs), each of length equal to the number of mutations in
          `ts_with_mutations`. For each mutation, `errors` gives 0 if a mutation can be correctly
          placed on the ARG to give that genotype at that position, 1 if not, and a value in between
          if random resolution of polytomies sometimes gives a match and sometimes does not. `diffs`
          gives the Hamming distance between the best possible placed mutation and the genotype from
          `ts_with_mutations`, with the mean Hamming distance in the case of random resolution of
          polytomies.
    """
    errors = []
    diffs = []
    for k, variant in enumerate(ts_with_mutations.variants()):
        if resolve_reps == 0:
            diff = arg_needle_lib_pybind.mutation_best(arg, variant.site.position, variant.genotypes, 0)
            errors.append(int(diff > 0))
            diffs.append(diff)
        else:
            resolve_errors = []
            resolve_diffs = []
            for random_id in range(resolve_reps):
                resolve_seed = (k + 1) * resolve_reps + random_id  # changed from k to k + 1 on 2022/01/20
                diff = arg_needle_lib_pybind.mutation_best(arg, variant.site.position, variant.genotypes, resolve_seed)
                resolve_errors.append(int(diff > 0))
                resolve_diffs.append(diff)
            errors.append(np.mean(resolve_errors))
            diffs.append(np.mean(resolve_diffs))

    return errors, diffs


def mutation_match_binned(arg, ts_with_mutations, mac_thresholds, resolve_reps=0):
    """Tries to best match a given set of mutations on to an ARG, stratifying by allele frequency.

    Args:
        arg: `arg_needle_lib.ARG` object (usually one that has been inferred)
        ts_with_mutations: `tskit.TreeSequence` object containing mutations, that we want to
          check for on `arg`
        mac_thresholds: a list of minor allele count thresholds to divide up the mutations into bins.
          The number of bins will be len(mac_thresholds) + 1. The thresholds should be in increasing
          order.
        resolve_reps (int): how many times to randomly resolve polytomies in the `arg` object
          when matching. This is unnecessary if `arg` contains no polytomies, but otherwise
          allows for matching against resolved binary trees, averaging over several random
          resolutions. Default 0, which corresponds to no resolving of polytomies.

    Returns:
        Three lists (binned_errors, binned_diffs, binned_sites), each of length len(mac_thresholds) + 1.
          The MAC bins go from [0, mac_thresholds[0]) for the smallest bin, up to [mac_thresholds[-1], 0.5]
          for the largest bin. binned_errors gives the sum of errors (mutations that cannot be placed
          correctly) for all sites in the MAC bin. binned_diffs gives the sum of diffs, which is the
          Hamming distance for the best possible mutation placement, for all sites in the MAC bin.
          binned_sites gives the count of sites in each MAC bin. Note that binned_errors and binned_diffs
          can contain non-integer values if resolve_reps is not 0.
        Additionally, two lists for all errors and diffs in order of the sites are also returned
          (see function mutation_match_all).

    Raises:
        ValueError: mac_thresholds is not in increasing order.
    """
    if sorted(mac_thresholds) != mac_thresholds:
        raise ValueError("MAC thresholds should be in increasing order")

    errors, diffs = mutation_match_all(arg, ts_with_mutations, resolve_reps)

    allele_counts = []
    for variant in ts_with_mutations.variants():
        allele_counts.append(np.sum(variant.genotypes))
    allele_counts = np.array(allele_counts)
    minor_allele_counts = np.minimum(allele_counts, ts_with_mutations.num_samples - allele_counts)
    minor_bins = np.searchsorted(mac_thresholds, minor_allele_counts, side="right")

    binned_errors = [0 for i in range(len(mac_thresholds)+1)]
    binned_diffs = [0 for i in range(len(mac_thresholds)+1)]
    binned_sites = [0 for i in range(len(mac_thresholds)+1)]
    for minor_bin, error, diff in zip(minor_bins, errors, diffs):
        binned_errors[minor_bin] += error
        binned_diffs[minor_bin] += diff
        binned_sites[minor_bin] += 1

    return binned_errors, binned_diffs, binned_sites, errors, diffs
