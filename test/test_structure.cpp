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
#include "types.hpp"

#include <catch2/catch_test_macros.hpp>

#include <vector>

using rvec = std::vector<arg_real_t>;
using ivec = std::vector<int>;
using rpair = std::pair<arg_real_t, arg_real_t>;
using ipair = std::pair<int, int>;

TEST_CASE("ARG structure", "[check][thread]") {
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
              std::vector<ipair>{ipair(0, 3), ipair(1, 3), ipair(2, 3)},
              std::vector<rpair>{rpair(0, 100), rpair(0, 100), rpair(0, 100)});
  }

  REQUIRE_NOTHROW(arg.check_basic());
  arg.populate_children_and_roots();
  REQUIRE_NOTHROW(arg.check_roots());
  REQUIRE_NOTHROW(arg.check_children());
}

TEST_CASE("Bad ARGs basic check", "[check]") {
  ARG arg(0, 42, 42); // dummy object

  // Make some bad ARGs with the special constructor and check for failure
  // Many of these require asserts to be turned off
#ifdef NDEBUG
  // there's no good way to check node spans
  SECTION("ARG span") {
    arg = ARG(100, 0);
    arg.add_sample();
  }
  SECTION("heights") {
    SECTION("leaf height") {
      arg = ARG(0, 100, rvec{3.14}, std::deque<bool>{true}, std::vector<ipair>{}, std::vector<rpair>{});
    }
    SECTION("non-leaf height") {
      arg =
          ARG(0, 100, rvec{0, -2.718, -2.718}, std::deque<bool>{true, false, false},
              std::vector<ipair>{ipair(1, 0), ipair(2, 0)}, std::vector<rpair>{rpair(0, 100), rpair(0, 100)});
    }
  }
  SECTION("edges") {
    SECTION("edges of too long span") {
      arg = ARG(0, 100, rvec{0, 0, 2.718}, std::deque<bool>{true, true, false},
                std::vector<ipair>{ipair(0, 2), ipair(1, 2)},
                std::vector<rpair>{rpair(0, 105), rpair(-5, 100)});
    }
    SECTION("upside down edges") {
      arg =
          ARG(0, 100, rvec{0, 0, 2.718}, std::deque<bool>{true, true, false},
              std::vector<ipair>{ipair(2, 0), ipair(2, 1)}, std::vector<rpair>{rpair(0, 100), rpair(0, 100)});
    }
    SECTION("edges of same height") {
      arg =
          ARG(0, 100, rvec{0, 0, 2.718}, std::deque<bool>{true, true, false},
              std::vector<ipair>{ipair(0, 2), ipair(1, 0)}, std::vector<rpair>{rpair(0, 100), rpair(0, 100)});
    }
  }
#endif

  SECTION("single parent except root gaps") {
    SECTION("multiple roots") {
      arg = ARG(0, 100, rvec{0, 0, 2.718}, std::deque<bool>{true, true, false},
                std::vector<ipair>{ipair(0, 2), ipair(1, 2), ipair(1, 2)},
                std::vector<rpair>{rpair(0, 100), rpair(0, 50), rpair(70, 100)});
      // when we turn off stringent checking, it should pass
      REQUIRE_NOTHROW(arg.check_basic(false));
    }
  }

  REQUIRE_THROWS(arg.check_basic());
}

// Make some bad ARGs with the special constructor and check for failure
TEST_CASE("Bad ARGs children check", "[check]") {
  ARG arg(0, 42, 42); // dummy object

  SECTION("0 children") {
    arg = ARG(0, 100, rvec{0, 0, 2.718, 3.14}, std::deque<bool>{true, true, false, false},
              std::vector<ipair>{ipair(0, 2), ipair(1, 2)}, std::vector<rpair>{rpair(0, 100), rpair(0, 100)});
    REQUIRE_THROWS(arg.check_basic()); // fails local root partition test
    // when we turn off stringent checking, it should pass
    REQUIRE_NOTHROW(arg.check_basic(false));
  }
  SECTION("1 child") {
    SECTION("whole extent") {
      arg = ARG(0, 100, rvec{0, 0, 2.718, 3.14}, std::deque<bool>{true, true, false, false},
                std::vector<ipair>{ipair(0, 2), ipair(1, 2), ipair(2, 3)},
                std::vector<rpair>{rpair(0, 100), rpair(0, 100), rpair(0, 100)});
      REQUIRE_NOTHROW(arg.check_basic());
    }
    SECTION("partial extent") {
      arg = ARG(0, 100, rvec{0, 0, 2.718, 3.14}, std::deque<bool>{true, true, false, false},
                std::vector<ipair>{ipair(0, 3), ipair(1, 3), ipair(0, 2), ipair(1, 2), ipair(2, 3)},
                std::vector<rpair>{
                    rpair(0, 50), rpair(0, 50), rpair(50, 100), rpair(50, 100), rpair(50, 100)});
      REQUIRE_THROWS(arg.check_basic()); // fails local root partition test
      // when we turn off stringent checking, it should pass
      REQUIRE_NOTHROW(arg.check_basic(false));
    }
  }

  arg.populate_children_and_roots();
  REQUIRE_NOTHROW(arg.check_roots());
  REQUIRE_THROWS(arg.check_children());
}
