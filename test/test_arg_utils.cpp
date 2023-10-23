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

#include "arg.hpp"
#include "arg_utils.hpp"
#include "types.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <vector>

using std::vector;

typedef vector<arg_real_t> rvec;
typedef vector<int> ivec;
typedef std::pair<arg_real_t, arg_real_t> rpair;
typedef std::pair<int, int> ipair;

TEST_CASE("ARG visit identical") {
  ARG arg(0, 42, 42); // dummy object

  SECTION("one sample") {
    arg = ARG(0, 100, 0);
    arg.add_sample();
  }
  SECTION("three samples") {
    arg = ARG(0, 100, 3);
    arg.add_sample("first");
    arg.add_sample("second");
    arg.thread_sample(rvec{0, 50}, ivec{0, 0}, rvec{3.14, 2.718});
    arg.add_sample("third");
    arg.thread_sample(rvec{0, 40}, ivec{1, 0}, rvec{5, 3});
  }
  SECTION("single root but multiple trees") {
    arg = ARG(0, 100, 3);
    arg.add_sample("first");
    arg.add_sample("second");
    arg.thread_sample(rvec{0}, ivec{0}, rvec{10});
    arg.add_sample("third");
    arg.thread_sample(rvec{0, 30, 70}, ivec{1, 1, 1}, rvec{3, 5, 4});
  }
  SECTION("multiple merges") {
    arg = ARG(0, 100, rvec{0, 0, 0, 5}, std::deque<bool>{true, true, true, false},
              vector<ipair>{ipair(0, 3), ipair(1, 3), ipair(2, 3)},
              vector<rpair>{rpair(0, 100), rpair(0, 100), rpair(0, 100)});
  }

  arg.populate_children_and_roots();
  REQUIRE(arg_utils::visit_identical(arg, 1e-9, 1e-9, false));
}

TEST_CASE("KC topology polytomies", "[!mayfail]") {
  ARG arg_tournament =
      ARG(0, 100, rvec{0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3},
          std::deque<bool>{true, true, true, true, true, true, true, true, false, false, false, false,
                           false, false, false},
          vector<ipair>{ipair(0, 8), ipair(1, 8), ipair(2, 9), ipair(3, 9), ipair(4, 10),
                        ipair(5, 10), ipair(6, 11), ipair(7, 11), ipair(8, 12), ipair(9, 12),
                        ipair(10, 13), ipair(11, 13), ipair(12, 14), ipair(13, 14)},
          vector<rpair>{rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100),
                        rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100),
                        rpair(0, 100), rpair(0, 100), rpair(0, 100),
                        rpair(0, 100)}); // https://xkcd.com/2269/
  ARG arg_polytomy = ARG(0, 100, rvec{0, 0, 0, 0, 0, 0, 0, 0, 1},
                         std::deque<bool>{true, true, true, true, true, true, true, true, false},
                         vector<ipair>{ipair(0, 8), ipair(1, 8), ipair(2, 8), ipair(3, 8),
                                       ipair(4, 8), ipair(5, 8), ipair(6, 8), ipair(7, 8)},
                         vector<rpair>{rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100),
                                       rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100)});

  arg_tournament.populate_children_and_roots();
  arg_polytomy.populate_children_and_roots();

  vector<int> tournament_distances = {
      2, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 1, 1, 1, 2};

  SECTION("regular") {
    vector<int> polytomy_distances(28, 0);
    arg_real_t kc_squared = 0;
    for (size_t i = 0; i < 28; ++i) {
      int diff = tournament_distances[i] - polytomy_distances[i];
      kc_squared += diff * diff;
    }
    arg_real_t regular_kc =
        std::get<0>(arg_utils::metrics_stab_efficient(arg_tournament, arg_polytomy, 1));
    arg_real_t regular_kc_old = arg_utils::kc_topology(arg_tournament, arg_polytomy);
    REQUIRE(std::fabs(kc_squared - regular_kc) < 1e-9);
    REQUIRE(std::fabs(kc_squared - regular_kc_old) < 1e-9);
  }

  SECTION("random_seed_2") {
    /* Random tree with seed 2
             14
           /   \
         10     13
        /  \    / \
       8    9  6   12
      /\   / \    /  \
     2  7  1  3  0    11
                     /  \
                    4    5

    Can be obtained by running

        std::mt19937 gen(2);
        vector<vector<int>> tree = arg_utils::random_binary_tree(8, gen);
        for (size_t i = 0; i < tree.size(); ++i) {
          cout << i << " " << tree[i][0] << " " << tree[i][1] << endl;
        }

    Unfortunately, this is not consistent across architectures, so we allow
    the test to fail.
     */
    vector<int> polytomy_distances = {
        0, 0, 0, 2, 2, 1, 0, 1, 2, 0, 0, 0, 1, 1, 0, 0, 0, 2, 0, 0, 0, 1, 3, 1, 0, 1, 0, 0};
    arg_real_t kc_squared = 0;
    for (size_t i = 0; i < 28; ++i) {
      int diff = tournament_distances[i] - polytomy_distances[i];
      kc_squared += diff * diff;
    }
    // we need to put the polytomy ARG first so we can know what the RNG result will be
    arg_real_t random_kc =
        std::get<0>(arg_utils::metrics_stab_efficient(arg_polytomy, arg_tournament, 1, 2));
    REQUIRE(std::fabs(kc_squared - random_kc) < 1e-9);
  }

  SECTION("random_seed_8") {
    /* Random tree with seed 8
       14
      /  \
     0   13
        /  \
       3   12
          /  \
         7    11
             /  \
            1    10
                /  \
               6    9
                   / \
                  5   8
                     / \
                    2   4

    Can be obtained by running

        std::mt19937 gen(8);
        vector<vector<int>> tree = arg_utils::random_binary_tree(8, gen);
        for (size_t i = 0; i < tree.size(); ++i) {
          cout << i << " " << tree[i][0] << " " << tree[i][1] << endl;
        }

    Unfortunately, this is not consistent across architectures, so we allow
    the test to fail.
     */
    vector<int> polytomy_distances = {
        0, 0, 0, 0, 0, 0, 0, 3, 1, 3, 3, 3, 2, 1, 6, 5, 4, 2, 1, 1, 1, 1, 5, 4, 2, 4, 2, 2};
    arg_real_t kc_squared = 0;
    for (size_t i = 0; i < 28; ++i) {
      int diff = tournament_distances[i] - polytomy_distances[i];
      kc_squared += diff * diff;
    }
    // we need to put the polytomy ARG first so we can know what the RNG result will be
    arg_real_t random_kc =
        std::get<0>(arg_utils::metrics_stab_efficient(arg_polytomy, arg_tournament, 1, 8));
    REQUIRE(std::fabs(kc_squared - random_kc) < 1e-9);
  }
}

TEST_CASE("Mutation best diff polytomies", "[!mayfail]") {
  ARG arg_polytomy = ARG(0, 100, rvec{0, 0, 0, 0, 0, 0, 0, 0, 1},
                         std::deque<bool>{true, true, true, true, true, true, true, true, false},
                         vector<ipair>{ipair(0, 8), ipair(1, 8), ipair(2, 8), ipair(3, 8),
                                       ipair(4, 8), ipair(5, 8), ipair(6, 8), ipair(7, 8)},
                         vector<rpair>{rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100),
                                       rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100)});

  arg_polytomy.populate_children_and_roots();
  arg_real_t position = 0;
  vector<int> genotypes = {1, 1, 1, 1, 0, 0, 0, 0};
  vector<int> genotypes2 = {1, 1, 0, 0, 1, 1, 0, 0};
  vector<int> genotypes3 = {0, 0, 0, 1, 0, 0, 0, 0};

  SECTION("regular") {
    int result = arg_utils::mutation_best(arg_polytomy, position, genotypes);
    int expected = 3;
    REQUIRE(result == expected);

    result = arg_utils::mutation_best(arg_polytomy, position, genotypes2);
    expected = 3;
    REQUIRE(result == expected);

    result = arg_utils::mutation_best(arg_polytomy, position, genotypes3);
    expected = 0;
    REQUIRE(result == expected);
  }

  SECTION("random_seed_2") {
    /* Random tree with seed 2
             14
           /   \
         10     13
        /  \    / \
       8    9  6   12
      /\   / \    /  \
     2  7  1  3  0    11
                     /  \
                    4    5

    Can be obtained by running

        std::mt19937 gen(2);
        vector<vector<int>> tree = arg_utils::random_binary_tree(8, gen);
        for (size_t i = 0; i < tree.size(); ++i) {
          cout << i << " " << tree[i][0] << " " << tree[i][1] << endl;
        }

    Unfortunately, this is not consistent across architectures, so we allow
    the test to fail.
     */
    unsigned int random_seed = 2;
    int result = arg_utils::mutation_best(arg_polytomy, position, genotypes, random_seed);
    // by placing the mutation between 10 and 14, we get 6 correct and 2 incorrect
    int expected = 2;
    REQUIRE(result == expected);

    result = arg_utils::mutation_best(arg_polytomy, position, genotypes2, random_seed);
    // by placing the mutation between 12 and 13, we get 7 correct and 1 incorrect
    expected = 1;
    REQUIRE(result == expected);

    result = arg_utils::mutation_best(arg_polytomy, position, genotypes3, random_seed);
    expected = 0;
    REQUIRE(result == expected);
  }

  SECTION("random_seed_8") {
    /* Random tree with seed 8
       14
      /  \
     0   13
        /  \
       3   12
          /  \
         7    11
             /  \
            1    10
                /  \
               6    9
                   / \
                  5   8
                     / \
                    2   4

    Can be obtained by running

        std::mt19937 gen(8);
        vector<vector<int>> tree = arg_utils::random_binary_tree(8, gen);
        for (size_t i = 0; i < tree.size(); ++i) {
          cout << i << " " << tree[i][0] << " " << tree[i][1] << endl;
        }

    Unfortunately, this is not consistent across architectures, so we allow
    the test to fail.
     */
    unsigned int random_seed = 8;
    int result = arg_utils::mutation_best(arg_polytomy, position, genotypes, random_seed);
    // by placing the mutation between 13 and 12, we get 6 correct and 2 incorrect
    int expected = 2;
    REQUIRE(result == expected);

    result = arg_utils::mutation_best(arg_polytomy, position, genotypes2, random_seed);
    // by placing the mutation between 9 and 10, we get 5 correct and 3 incorrect
    expected = 3;
    REQUIRE(result == expected);

    result = arg_utils::mutation_best(arg_polytomy, position, genotypes3, random_seed);
    expected = 0;
    REQUIRE(result == expected);
  }
}

TEST_CASE("Precision recall") {
  ARG arg1 = ARG(0, 100, 3);
  arg1.add_sample("first");
  arg1.add_sample("second");
  arg1.thread_sample(rvec{0, 50}, ivec{0, 0}, rvec{3.14, 2.718});
  arg1.add_sample("third");
  arg1.thread_sample(rvec{0, 40}, ivec{1, 0}, rvec{5, 3});

  SECTION("same") {
    ARG arg2 = ARG(0, 100, 3);
    arg2.add_sample("first");
    arg2.add_sample("second");
    arg2.thread_sample(rvec{0, 50}, ivec{0, 0}, rvec{3.14, 2.718});
    arg2.add_sample("third");
    arg2.thread_sample(rvec{0, 40}, ivec{1, 0}, rvec{5, 3});

    arg1.populate_children_and_roots();
    arg2.populate_children_and_roots();

    int num1, num2, num_common;
    arg_real_t length1, length2, length_common;
    std::tie(num1, num2, num_common, length1, length2, length_common) =
        arg_utils::bitset_overlap_stab(arg1, arg2, 1);
    REQUIRE(num1 == num2);
    REQUIRE(num1 == num_common);
    REQUIRE(length1 == length2);
    REQUIRE(length1 == length_common);

    std::tie(num1, num2, num_common, length1, length2, length_common) =
        arg_utils::bitset_overlap_stab(arg1, arg2, 50);
    REQUIRE(num1 == num2);
    REQUIRE(num1 == num_common);
    REQUIRE(length1 == length2);
    REQUIRE(length1 == length_common);

    std::tie(num1, num2, num_common) = arg_utils::bitset_overlap_full(arg1, arg2);
    REQUIRE(num1 == num2);
    REQUIRE(num1 == num_common);
  }

  SECTION("different, multiple merger") {
    ARG arg2 = ARG(0, 100, rvec{0, 0, 0, 5}, std::deque<bool>{true, true, true, false},
                   vector<ipair>{ipair(0, 3), ipair(1, 3), ipair(2, 3)},
                   vector<rpair>{rpair(0, 100), rpair(0, 100), rpair(0, 100)});

    arg1.populate_children_and_roots();
    arg2.populate_children_and_roots();

    int num1, num2, num_common;
    arg_real_t length1, length2, length_common;
    std::tie(num1, num2, num_common, length1, length2, length_common) =
        arg_utils::bitset_overlap_stab(arg1, arg2, 1);
    REQUIRE(num1 == 4);
    REQUIRE(num2 == 3);
    REQUIRE(num_common == 3); // all singletons from arg2 will be in arg1
    REQUIRE(length1 == 2 * 3 + 2.718);
    REQUIRE(length2 == 3 * 5);
    REQUIRE(length_common == 2 * 2.718 + 3); // only count length from singletons

    std::tie(num1, num2, num_common) = arg_utils::bitset_overlap_full(arg1, arg2);
    REQUIRE(num1 == 5);
    REQUIRE(num2 == 3);
    REQUIRE(num_common == 3);
  }

  SECTION("different, maximally disjoint") {
    ARG arg2 = ARG(0, 100, 3);
    arg2.add_sample("first");
    arg2.add_sample("second");
    arg2.thread_sample(rvec{0, 50}, ivec{0, 0}, rvec{3.14, 2.718});
    arg2.add_sample("third");
    arg2.thread_sample(rvec{0}, ivec{1}, rvec{2});

    arg1.populate_children_and_roots();
    arg2.populate_children_and_roots();

    int num1, num2, num_common;
    arg_real_t length1, length2, length_common;
    std::tie(num1, num2, num_common, length1, length2, length_common) =
        arg_utils::bitset_overlap_stab(arg1, arg2, 1);
    REQUIRE(num1 == 4);
    REQUIRE(num2 == 4);
    REQUIRE(num_common == 3); // all singletons from arg2 will be in arg1
    REQUIRE(length1 == 2 * 3 + 2.718);
    REQUIRE(length2 == 2 * 2.718 + 2);
    REQUIRE(length_common ==
            2 * 2 + 2.718); // singletons: min(2.718, 2.718) + min(2.718, 2) + min(3, 2)

    std::tie(num1, num2, num_common) = arg_utils::bitset_overlap_full(arg1, arg2);
    REQUIRE(num1 == 5);
    REQUIRE(num2 == 4);
    REQUIRE(num_common == 3);
  }
}

TEST_CASE("Imputation standard") {
  /* ARG in question:
              ___ 12 _
             /        11
            10       /  \
           / \       9   |
          /   8     / \  |
         /   / \    | |  |
        /   /   7   | |  |
       /   /   / \  | |  |
      0   1   2   3 4 5  6
   */
  ARG arg = ARG(0, 100, rvec{0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6},
                std::deque<bool>{true, true, true, true, true, true, true, false, false, false, false,
                                 false, false},
                vector<ipair>{ipair(2, 7), ipair(3, 7), ipair(1, 8), ipair(7, 8), ipair(4, 9),
                              ipair(5, 9), ipair(0, 10), ipair(8, 10), ipair(6, 11), ipair(9, 11),
                              ipair(10, 12), ipair(11, 12)},
                vector<rpair>{rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100),
                              rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100),
                              rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100)});
  arg.populate_children_and_roots();

  vector<int> start;
  vector<arg_real_t> expected;

  SECTION("all zeros") {
    start = {-1, -1, 0, -1, 0, -1, -1};
    expected = {0, 0, 0, 0, 0, 0, 0};
  }

  SECTION("all ones") {
    start = {-1, 1, -1, -1, -1, -1, -1};
    expected = {1, 1, 1, 1, 1, 1, 1};
  }

  SECTION("easy 1") {
    start = {-1, 0, 1, 1, 0, 0, 0};
    expected = {0, 0, 1, 1, 0, 0, 0};
  }

  SECTION("easy 2") {
    start = {-1, 0, 0, 1, 0, 0, 0};
    expected = {0, 0, 0, 1, 0, 0, 0};
  }

  SECTION("easy 3") {
    start = {0, 1, -1, 1, 0, 0, 0};
    expected = {0, 1, 1, 1, 0, 0, 0};
  }

  SECTION("real impute") {
    start = {1, -1, -1, 1, 0, 0, 0};
    expected = {1, 1, 1, 1, 0, 0, 0};
  }

  SECTION("error 1") { // we shouldn't correct the error but could perhaps alert
    start = {1, 1, 1, 0, -1, 0, 0};
    expected = {1, 1, 1, 0, 0, 0, 0};
  }

  SECTION("error 2") { // we shouldn't correct the error but could perhaps alert
    start = {1, 1, -1, 1, 0, 1, 0};
    expected = {1, 1, 1, 1, 0, 1, 0};
  }

  vector<arg_real_t> result = arg_utils::impute(arg, 50, start);
  REQUIRE(result.size() == expected.size());
  for (size_t i = 0; i < result.size(); ++i) {
    REQUIRE(std::fabs(result[i] - expected[i]) < 1e-6);
  }

  vector<arg_real_t> result_old = arg_utils::impute(arg, 50, start, true);
  REQUIRE(result_old.size() == expected.size());
  for (size_t i = 0; i < result_old.size(); ++i) {
    REQUIRE(result_old[i] == expected[i]);
  }
}

TEST_CASE("Imputation tiebreaker") {
  /* ARG in question:
              ___ 12 _
             /        11
            10       /  \
           / \       9   |
          /   8     / \  |
         /   / \    | |  |
        /   /   7   | |  |
       /   /   / \  | |  |
      0   1   2   3 4 5  6

  Node 10 is either closer to node 12 or closer to node 8, which affects how
  sample 0 is imputed. Note that we also need to take into account the
  distance between 11 and 12.
   */
  ARG arg(0, 42, 42); // dummy object
  vector<int> start = {-1, 1, 1, 1, 0, 0, 0};
  vector<arg_real_t> expected;
  bool old;

  SECTION("upper case old") {
    // in the upper case, sample 0 joins a bit higher, so it gets a 0 when imputed
    arg = ARG(0, 100, rvec{0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 5, 5, 6},
              std::deque<bool>{true, true, true, true, true, true, true, false, false, false, false,
                               false, false},
              vector<ipair>{ipair(2, 7), ipair(3, 7), ipair(1, 8), ipair(7, 8), ipair(4, 9),
                            ipair(5, 9), ipair(0, 10), ipair(8, 10), ipair(6, 11), ipair(9, 11),
                            ipair(10, 12), ipair(11, 12)},
              vector<rpair>{rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100),
                            rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100),
                            rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100)});
    old = true;
    expected = {0, 1, 1, 1, 0, 0, 0};
  }

  SECTION("lower case old") {
    // in the lower case, sample 0 joins a bit lower, so it gets a -1 when imputed
    arg = ARG(0, 100, rvec{0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 3, 5, 6},
              std::deque<bool>{true, true, true, true, true, true, true, false, false, false, false,
                               false, false},
              vector<ipair>{ipair(2, 7), ipair(3, 7), ipair(1, 8), ipair(7, 8), ipair(4, 9),
                            ipair(5, 9), ipair(0, 10), ipair(8, 10), ipair(6, 11), ipair(9, 11),
                            ipair(10, 12), ipair(11, 12)},
              vector<rpair>{rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100),
                            rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100),
                            rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100)});
    old = true;
    expected = {1, 1, 1, 1, 0, 0, 0};
  }

  SECTION("upper case new") {
    // in the upper case, sample 0 joins a bit higher, but now we give a probability
    arg = ARG(0, 100, rvec{0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 5, 5, 6},
              std::deque<bool>{true, true, true, true, true, true, true, false, false, false, false,
                               false, false},
              vector<ipair>{ipair(2, 7), ipair(3, 7), ipair(1, 8), ipair(7, 8), ipair(4, 9),
                            ipair(5, 9), ipair(0, 10), ipair(8, 10), ipair(6, 11), ipair(9, 11),
                            ipair(10, 12), ipair(11, 12)},
              vector<rpair>{rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100),
                            rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100),
                            rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100)});
    old = false;
    expected = {0.4, 1, 1, 1, 0, 0, 0};
  }

  SECTION("lower case new") {
    // in the lower case, sample 0 joins a bit lower, but now we give a probability
    arg = ARG(0, 100, rvec{0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 3, 5, 6},
              std::deque<bool>{true, true, true, true, true, true, true, false, false, false, false,
                               false, false},
              vector<ipair>{ipair(2, 7), ipair(3, 7), ipair(1, 8), ipair(7, 8), ipair(4, 9),
                            ipair(5, 9), ipair(0, 10), ipair(8, 10), ipair(6, 11), ipair(9, 11),
                            ipair(10, 12), ipair(11, 12)},
              vector<rpair>{rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100),
                            rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100),
                            rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100)});
    old = false;
    expected = {0.8, 1, 1, 1, 0, 0, 0};
  }

  arg.populate_children_and_roots();
  vector<arg_real_t> result = arg_utils::impute(arg, 50, start, old);
  REQUIRE(result.size() == expected.size());
  for (size_t i = 0; i < result.size(); ++i) {
    REQUIRE(std::fabs(result[i] - expected[i]) < 1e-6);
  }
}

TEST_CASE("Imputation polytomies") {
  /* ARG in question:
           ____ 11 ___
          /      |    \
        _8_     9     10
       / | \   / \   / | \
      0  1  2 3   4 5  6  7
   */
  ARG arg_polytomy = ARG(
      0, 100, rvec{0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2},
      std::deque<bool>{true, true, true, true, true, true, true, true, false, false, false, false},
      vector<ipair>{ipair(0, 8), ipair(1, 8), ipair(2, 8), ipair(3, 9), ipair(4, 9), ipair(5, 10),
                    ipair(6, 10), ipair(7, 10), ipair(8, 11), ipair(9, 11), ipair(10, 11)},
      vector<rpair>{rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100),
                    rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100),
                    rpair(0, 100)});
  arg_polytomy.populate_children_and_roots();

  vector<int> start;
  vector<arg_real_t> expected;

  SECTION("upper polytomy") {
    // note that for this case 1 is ancestral and 0 is derived. can change later
    start = {-1, 1, -1, 1, -1, 0, 0, -1};
    expected = {1, 1, 1, 1, 1, 0, 0, 0};
  }

  SECTION("lower polytomy") {
    start = {1, 0, -1, -1, -1, -1, 0, -1};
    expected = {1, 0, 0, 0, 0, 0, 0, 0};
  }

  vector<arg_real_t> result = arg_utils::impute(arg_polytomy, 50, start);
  REQUIRE(result.size() == expected.size());
  for (size_t i = 0; i < result.size(); ++i) {
    REQUIRE(std::fabs(result[i] - expected[i]) < 1e-6);
  }

  vector<arg_real_t> result_old = arg_utils::impute(arg_polytomy, 50, start, true);
  REQUIRE(result_old.size() == expected.size());
  for (size_t i = 0; i < result_old.size(); ++i) {
    REQUIRE(result_old[i] == expected[i]);
  }
}

TEST_CASE("Volume Tests") {
  // Test of function to quickly get volume
  ARG arg(0, 42, 42); // dummy object

  SECTION("two samples - single tree") {
    arg = ARG(0, 100, 2);
    arg.add_sample("first");
    arg.add_sample("second");
    arg.thread_sample(rvec{0}, ivec{0}, rvec{10});
    // expected volume should be 2 * 10 * 100 = 2000
    arg_real_t expected = 2000;
    arg_real_t result = arg_utils::total_volume(arg);
    REQUIRE(result == expected);
  }

  SECTION("two samples - multiple trees") {
    arg = ARG(0, 100, 2);
    arg.add_sample("first");
    arg.add_sample("second");
    arg.thread_sample(rvec{0, 25, 50, 75}, ivec{0, 0, 0, 0}, rvec{10, 5, 2, 5});
    arg_real_t expected = 2 * (10 * 25 + 5 * 25 + 2 * 25 + 5 * 25);
    arg_real_t result = arg_utils::total_volume(arg);
    REQUIRE(result == expected);
  }

  SECTION("three samples - two trees - same topology") {
    arg = ARG(0, 100, 3);
    arg.add_sample("first");
    arg.add_sample("second");
    arg.thread_sample(rvec{0}, ivec{0}, rvec{10});
    arg.add_sample("third");
    arg.thread_sample(rvec{0, 50}, ivec{0, 0}, rvec{5, 7});
    arg_real_t expected = 2 * (10 * 100) + 5 * 50 + 7 * 50;
    arg_real_t result = arg_utils::total_volume(arg);
    REQUIRE(result == expected);
  }
}

TEST_CASE("Mutation Tests") {

  SECTION("Simulate 1 mutation") {
    ARG arg = ARG(0, 100, 3);
    arg.add_sample("first");
    arg.add_sample("second");
    arg.thread_sample(rvec{0}, ivec{0}, rvec{10});
    arg.add_sample("third");
    arg.thread_sample(rvec{0, 50}, ivec{0, 0}, rvec{5, 7});
    int expected = 1;
    arg_utils::generate_m_mutations(arg, 1, 0);
    int result = arg.num_mutations();
    REQUIRE(result == expected);
  }

  SECTION("Simulate 10 mutations") {
    ARG arg = ARG(0, 100, 3);
    arg.add_sample("first");
    arg.add_sample("second");
    arg.thread_sample(rvec{0}, ivec{0}, rvec{10});
    arg.add_sample("third");
    arg.thread_sample(rvec{0, 50}, ivec{0, 0}, rvec{5, 7});
    int expected = 10;
    arg_utils::generate_m_mutations(arg, 10, 0);
    int result = arg.num_mutations();
    REQUIRE(result == expected);
  }

  SECTION("Simulate mutations") {
    ARG arg = ARG(0, 100, 3);
    arg.add_sample("first");
    arg.add_sample("second");
    arg.thread_sample(rvec{0}, ivec{0}, rvec{10});
    arg.add_sample("third");
    arg.thread_sample(rvec{0, 50}, ivec{0, 0}, rvec{5, 7});
    arg_utils::generate_mutations(arg, 1e-2, 0, 1234ul);
    int result = arg.num_mutations();
    REQUIRE(result >= 0);
  }

  SECTION("Test adding mutations results in sorted vector") {
    ARG arg = ARG(0, 100, 3);
    arg.add_mutation(nullptr, 0.1);
    arg.add_mutation(nullptr, 0.3);
    arg.add_mutation(nullptr, 0.2);

    REQUIRE(arg.get_mutations().size() == 3ul);

    // The mutations should always be sorted by position no matter what order they're added
    REQUIRE(arg.get_mutations().at(0ul)->position == 0.1);
    REQUIRE(arg.get_mutations().at(1ul)->position == 0.2);
    REQUIRE(arg.get_mutations().at(2ul)->position == 0.3);
  }
}

TEST_CASE("Number of Lineages Tests") {
  /* ARG in question:
             ____ 11 ____
            /     |      \
           /      |       10
          /       9       /|\
        _8_      / \     / | \
       / | \    /   \   /  |  \
      0  1  2  3     4 5   6   7
   */
  ARG arg_polytomy = ARG(
      0, 100, rvec{0, 0, 0, 0, 0, 0, 0, 0, 1, 1.3, 1.7, 2},
      std::deque<bool>{true, true, true, true, true, true, true, true, false, false, false, false},
      vector<ipair>{ipair(0, 8), ipair(1, 8), ipair(2, 8), ipair(3, 9), ipair(4, 9), ipair(5, 10),
                    ipair(6, 10), ipair(7, 10), ipair(8, 11), ipair(9, 11), ipair(10, 11)},
      vector<rpair>{rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100),
                    rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100), rpair(0, 100),
                    rpair(0, 100)});
  arg_polytomy.populate_children_and_roots();
  arg_real_t height;
  int expected;

  SECTION("above root") {
    height = 5;
    expected = 1;
  }

  SECTION("right on the nose") {
    height = 2;
    expected = 1;
  }

  SECTION("just below root") {
    height = 1.9;
    expected = 3;
  }

  SECTION("right on 10") {
    height = 1.7;
    expected = 3;
  }

  SECTION("between 9 and 10") {
    height = 1.5;
    expected = 5;
  }

  SECTION("right on 9") {
    height = 1.3;
    expected = 5;
  }

  SECTION("between 8 and 9") {
    height = 1.2;
    expected = 6;
  }
  SECTION("right on 8") {
    height = 1;
    expected = 6;
  }

  SECTION("between 8 and bottom") {
    height = 0.5;
    expected = 8;
  }

  SECTION("at the bottom") {
    height = 0;
    expected = 8;
  }

  int result = arg_utils::num_lineages(arg_polytomy, 50, height);
  REQUIRE(result == expected);
}

TEST_CASE("Test local volume async") {

    const std::string nodes_file_name = ARG_NEEDLE_TESTDATA_DIR "/length_1e6_samples_1e3/nodes.txt";
    const std::string edges_file_name = ARG_NEEDLE_TESTDATA_DIR "/length_1e6_samples_1e3/edges.txt";

    // Read in the ARG from a file
    ARG arg = arg_utils::arg_from_ts_files(nodes_file_name, edges_file_name);
    arg.populate_children_and_roots();

    // Test total volume calculation with varying numbers of tasks
    arg_real_t vol_v1 = arg_utils::total_volume(arg);
    arg_real_t vol_v2 = arg_utils::local_volume(arg, arg.start, arg.end);
    arg_real_t vol_v3 = arg_utils::local_volume(arg, arg.start, arg.end, 2u);
    arg_real_t vol_v4 = arg_utils::local_volume(arg, arg.start, arg.end, 3u);
    arg_real_t vol_v5 = arg_utils::local_volume(arg, arg.start, arg.end, 10u);
    arg_real_t vol_v6 = arg_utils::local_volume(arg, arg.start, arg.end, 100u);

    REQUIRE(vol_v1 == Catch::Approx(vol_v2).margin(1e-12));
    REQUIRE(vol_v1 == Catch::Approx(vol_v3).margin(1e-12));
    REQUIRE(vol_v1 == Catch::Approx(vol_v4).margin(1e-12));
    REQUIRE(vol_v1 == Catch::Approx(vol_v5).margin(1e-12));
    REQUIRE(vol_v1 == Catch::Approx(vol_v6).margin(1e-12));

    // Test local volume calculation on part of the length, with varying numbers of tasks
    arg_real_t vol_v7 = arg_utils::local_volume(arg, 2e5, 4e5);
    arg_real_t vol_v8 = arg_utils::local_volume(arg, 2e5, 4e5, 3u);
    arg_real_t vol_v9 = arg_utils::local_volume(arg, 2e5, 4e5, 7u);

    REQUIRE(vol_v7 == Catch::Approx(vol_v8).margin(1e-12));
    REQUIRE(vol_v7 == Catch::Approx(vol_v9).margin(1e-12));
}
