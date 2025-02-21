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

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>

#include <sstream>
#include <vector>
#include <boost/iostreams/filter/zlib.hpp>


TEST_CASE("Finding mutations left")
{

  SECTION("Empty mutation vector") {
    ARG arg = ARG(0, 100);
    CHECK_THROWS_WITH(arg.get_idx_of_first_mutation_left_of(1.0, false), Catch::Matchers::ContainsSubstring("There are no mutations"));
    CHECK_THROWS_WITH(arg.get_idx_of_first_mutation_left_of(1.0, true), Catch::Matchers::ContainsSubstring("There are no mutations"));
  }

  SECTION("One mutation in vector, don't include equal") {

    // Redirect std::cout to a stringstream, so we can test the warnings generated
    std::stringstream buffer;
    std::streambuf* prevCoutBuffer = std::cout.rdbuf(buffer.rdbuf());

    ARG arg = ARG(0, 100);
    arg.add_mutation(nullptr, 1.23);

    CHECK(arg.get_idx_of_first_mutation_left_of(1.0) == 0ul);
    CHECK(arg.get_idx_of_first_mutation_left_of(2.0) == 0ul);

    CHECK(buffer.str() == "Warning: no mutations with position < 1\n");
    std::cout.rdbuf(prevCoutBuffer); // Restore original buffer
  }

  SECTION("One mutation in vector, include equal") {

    // Redirect std::cout to a stringstream, so we can test the warnings generated
    std::stringstream buffer;
    std::streambuf* prevCoutBuffer = std::cout.rdbuf(buffer.rdbuf());

    ARG arg = ARG(0, 100);
    arg.add_mutation(nullptr, 1.23);

    CHECK(arg.get_idx_of_first_mutation_left_of(0.8, true) == 0ul);
    CHECK(arg.get_idx_of_first_mutation_left_of(1.23, true) == 0ul);
    CHECK(arg.get_idx_of_first_mutation_left_of(2.0, true) == 0ul);

    CHECK(buffer.str() == "Warning: no mutations with position <= 0.8\n");
    std::cout.rdbuf(prevCoutBuffer); // Restore original buffer
  }

  SECTION("Multiple mutations in vector, don't include equal") {

    // Redirect std::cout to a stringstream, so we can test the warnings generated
    std::stringstream buffer;
    std::streambuf* prevCoutBuffer = std::cout.rdbuf(buffer.rdbuf());

    ARG arg = ARG(0, 100);
    arg.add_mutation(nullptr, 1.23);
    arg.add_mutation(nullptr, 2.34);
    arg.add_mutation(nullptr, 3.45);

    CHECK(arg.get_idx_of_first_mutation_left_of(0.9) == 0ul);
    CHECK(arg.get_idx_of_first_mutation_left_of(2.0) == 0ul);
    CHECK(arg.get_idx_of_first_mutation_left_of(3.44) == 1ul);
    CHECK(arg.get_idx_of_first_mutation_left_of(3.45) == 1ul);
    CHECK(arg.get_idx_of_first_mutation_left_of(3.46) == 2ul);
    CHECK(arg.get_idx_of_first_mutation_left_of(4.0) == 2ul);
    CHECK(arg.get_idx_of_first_mutation_left_of(5.0) == 2ul);

    CHECK(buffer.str() == "Warning: no mutations with position < 0.9\n");
    std::cout.rdbuf(prevCoutBuffer); // Restore original buffer
  }

  SECTION("Multiple mutations in vector, include equal") {

    // Redirect std::cout to a stringstream, so we can test the warnings generated
    std::stringstream buffer;
    std::streambuf* prevCoutBuffer = std::cout.rdbuf(buffer.rdbuf());

    ARG arg = ARG(0, 100);
    arg.add_mutation(nullptr, 1.23);
    arg.add_mutation(nullptr, 2.34);
    arg.add_mutation(nullptr, 3.45);

    CHECK(arg.get_idx_of_first_mutation_left_of(0.9, true) == 0ul);
    CHECK(arg.get_idx_of_first_mutation_left_of(2.33, true) == 0ul);
    CHECK(arg.get_idx_of_first_mutation_left_of(2.34, true) == 1ul);
    CHECK(arg.get_idx_of_first_mutation_left_of(2.35, true) == 1ul);

    CHECK(buffer.str() == "Warning: no mutations with position <= 0.9\n");
    std::cout.rdbuf(prevCoutBuffer); // Restore original buffer
  }
}

TEST_CASE("Finding mutations right") {

  SECTION("Empty mutation vector") {
    ARG arg = ARG(0, 100);
    CHECK_THROWS_WITH(arg.get_idx_of_first_mutation_right_of(1.0, false), Catch::Matchers::ContainsSubstring("There are no mutations"));
    CHECK_THROWS_WITH(arg.get_idx_of_first_mutation_right_of(1.0, true), Catch::Matchers::ContainsSubstring("There are no mutations"));
  }

  SECTION("One mutation in vector, don't include equal") {

    // Redirect std::cout to a stringstream, so we can test the warnings generated
    std::stringstream buffer;
    std::streambuf* prevCoutBuffer = std::cout.rdbuf(buffer.rdbuf());

    ARG arg = ARG(0, 100);
    arg.add_mutation(nullptr, 1.23);

    CHECK(arg.get_idx_of_first_mutation_right_of(1.0) == 0ul);
    CHECK(arg.get_idx_of_first_mutation_right_of(2.0) == 0ul);

    CHECK(buffer.str() == "Warning: no mutations with position > 2\n");
    std::cout.rdbuf(prevCoutBuffer); // Restore original buffer
  }

  SECTION("One mutation in vector, include equal") {

    // Redirect std::cout to a stringstream, so we can test the warnings generated
    std::stringstream buffer;
    std::streambuf* prevCoutBuffer = std::cout.rdbuf(buffer.rdbuf());

    ARG arg = ARG(0, 100);
    arg.add_mutation(nullptr, 1.23);

    CHECK(arg.get_idx_of_first_mutation_right_of(0.8, true) == 0ul);
    CHECK(arg.get_idx_of_first_mutation_right_of(1.23, true) == 0ul);
    CHECK(arg.get_idx_of_first_mutation_right_of(2.0, true) == 0ul);

    CHECK(buffer.str() == "Warning: no mutations with position >= 2\n");
    std::cout.rdbuf(prevCoutBuffer); // Restore original buffer
  }

  SECTION("Multiple mutations in vector, don't include equal") {

    // Redirect std::cout to a stringstream, so we can test the warnings generated
    std::stringstream buffer;
    std::streambuf* prevCoutBuffer = std::cout.rdbuf(buffer.rdbuf());

    ARG arg = ARG(0, 100);
    arg.add_mutation(nullptr, 1.23);
    arg.add_mutation(nullptr, 2.34);
    arg.add_mutation(nullptr, 3.45);

    CHECK(arg.get_idx_of_first_mutation_right_of(1.0) == 0ul);
    CHECK(arg.get_idx_of_first_mutation_right_of(2.0) == 1ul);
    CHECK(arg.get_idx_of_first_mutation_right_of(3.44) == 2ul);
    CHECK(arg.get_idx_of_first_mutation_right_of(3.45) == 2ul);
    CHECK(arg.get_idx_of_first_mutation_right_of(3.46) == 2ul);

    CHECK(buffer.str() == "Warning: no mutations with position > 3.45\nWarning: no mutations with position > 3.46\n");
    std::cout.rdbuf(prevCoutBuffer); // Restore original buffer
  }

  SECTION("Multiple mutations in vector, include equal") {

    // Redirect std::cout to a stringstream, so we can test the warnings generated
    std::stringstream buffer;
    std::streambuf* prevCoutBuffer = std::cout.rdbuf(buffer.rdbuf());

    ARG arg = ARG(0, 100);
    arg.add_mutation(nullptr, 1.23);
    arg.add_mutation(nullptr, 2.34);
    arg.add_mutation(nullptr, 3.45);

    CHECK(arg.get_idx_of_first_mutation_right_of(1.0, true) == 0ul);
    CHECK(arg.get_idx_of_first_mutation_right_of(2.33, true) == 1ul);
    CHECK(arg.get_idx_of_first_mutation_right_of(2.34, true) == 1ul);
    CHECK(arg.get_idx_of_first_mutation_right_of(2.35, true) == 2ul);
    CHECK(arg.get_idx_of_first_mutation_right_of(5.0, true) == 2ul);

    CHECK(buffer.str() == "Warning: no mutations with position >= 5\n");
    std::cout.rdbuf(prevCoutBuffer); // Restore original buffer
  }
}

TEST_CASE("Finding closest mutation") {

  SECTION("Empty mutation vector") {
    ARG arg = ARG(0, 100);
    CHECK_THROWS_WITH(arg.get_idx_of_mutation_closest_to(1.0), Catch::Matchers::ContainsSubstring("There are no mutations"));
  }

  SECTION("One mutation in vector") {
    ARG arg = ARG(0, 100);
    arg.add_mutation(nullptr, 1.23);

    CHECK(arg.get_idx_of_mutation_closest_to(-1.0) == 0ul);
    CHECK(arg.get_idx_of_mutation_closest_to(0.0) == 0ul);
    CHECK(arg.get_idx_of_mutation_closest_to(1e10) == 0ul);
  }

  SECTION("Multiple mutations in vector") {

    ARG arg = ARG(0, 100);
    arg.add_mutation(nullptr, 0.5);
    arg.add_mutation(nullptr, 1.23);
    arg.add_mutation(nullptr, 2.0);
    arg.add_mutation(nullptr, 3.0);

    CHECK(arg.get_idx_of_mutation_closest_to(0.0) == 0ul);
    CHECK(arg.get_idx_of_mutation_closest_to(1.22) == 1ul);
    CHECK(arg.get_idx_of_mutation_closest_to(1.23) == 1ul);
    CHECK(arg.get_idx_of_mutation_closest_to(1.24) == 1ul);
    CHECK(arg.get_idx_of_mutation_closest_to(1.99) == 2ul);
    CHECK(arg.get_idx_of_mutation_closest_to(2.0) == 2ul);
    CHECK(arg.get_idx_of_mutation_closest_to(2.01) == 2ul);
    CHECK(arg.get_idx_of_mutation_closest_to(2.5) == 3ul);
    CHECK(arg.get_idx_of_mutation_closest_to(2.99) == 3ul);
    CHECK(arg.get_idx_of_mutation_closest_to(3.0) == 3ul);
    CHECK(arg.get_idx_of_mutation_closest_to(3.01) == 3ul);
    CHECK(arg.get_idx_of_mutation_closest_to(1e10) == 3ul);
  }
}

TEST_CASE("Test mutation_sites map") {

  SECTION("Empty mutation vector") {
    ARG arg = ARG(0, 100);
    CHECK(arg.get_mutation_sites().empty());
    CHECK(arg.get_site_positions().empty());
  }

  SECTION("Multiple mutations in vector") {

    ARG arg = ARG(0, 100);
    arg.add_mutation(nullptr, 1.0);
    arg.add_mutation(nullptr, 2.0);
    arg.add_mutation(nullptr, 3.0);

    std::map<arg_real_t, Site> mutation_sites = arg.get_mutation_sites();
    CHECK(mutation_sites.size() == 3ul);

    std::vector<arg_real_t> site_positions = arg.get_site_positions();
    CHECK(site_positions.size() == 3ul);

    CHECK(mutation_sites.at(1.0).get_mutations().size() == 1ul);
    CHECK(mutation_sites.at(1.0).get_mutations().front()->position == 1.0);
    CHECK(mutation_sites.at(1.0).get_position() == 1.0);

    CHECK(mutation_sites.at(2.0).get_mutations().size() == 1ul);
    CHECK(mutation_sites.at(2.0).get_mutations().front()->position == 2.0);
    CHECK(mutation_sites.at(2.0).get_position() == 2.0);

    CHECK(mutation_sites.at(3.0).get_mutations().size() == 1ul);
    CHECK(mutation_sites.at(3.0).get_mutations().front()->position == 3.0);
    CHECK(mutation_sites.at(3.0).get_position() == 3.0);

    CHECK(site_positions.at(0ul) == 1.0);
    CHECK(site_positions.at(1ul) == 2.0);
    CHECK(site_positions.at(2ul) == 3.0);

    // Add another mutation and check we're still OK
    arg.add_mutation(nullptr, 1.5);
    mutation_sites = arg.get_mutation_sites();
    site_positions = arg.get_site_positions();

    CHECK(mutation_sites.at(1.5).get_mutations().size() == 1ul);

    CHECK(site_positions.at(0ul) == 1.0);
    CHECK(site_positions.at(1ul) == 1.5);
    CHECK(site_positions.at(2ul) == 2.0);
    CHECK(site_positions.at(3ul) == 3.0);
  }

  SECTION("Multiple mutations at a single site") {

    ARG arg = ARG(0, 100);
    arg.add_mutation(nullptr, 1.0);
    arg.add_mutation(nullptr, 1.0);
    arg.add_mutation(nullptr, 2.0);
    arg.add_mutation(nullptr, 2.0);
    arg.add_mutation(nullptr, 2.0);
    arg.add_mutation(nullptr, 3.0);
    arg.add_mutation(nullptr, 4.0);

    const std::map<arg_real_t, Site> mutation_sites = arg.get_mutation_sites();
    CHECK(mutation_sites.size() == 4ul);

    const std::vector<arg_real_t> site_positions = arg.get_site_positions();
    CHECK(site_positions.size() == 4ul);

    CHECK(mutation_sites.at(1.0).get_mutations().size() == 2ul);
    CHECK(mutation_sites.at(1.0).get_mutations().at(0)->position == 1.0);
    CHECK(mutation_sites.at(1.0).get_mutations().at(1)->position == 1.0);

    CHECK(mutation_sites.at(2.0).get_mutations().size() == 3ul);
    CHECK(mutation_sites.at(2.0).get_mutations().at(0)->position == 2.0);
    CHECK(mutation_sites.at(2.0).get_mutations().at(1)->position == 2.0);
    CHECK(mutation_sites.at(2.0).get_mutations().at(2)->position == 2.0);

    CHECK(mutation_sites.at(3.0).get_mutations().size() == 1ul);
    CHECK(mutation_sites.at(3.0).get_mutations().at(0)->position == 3.0);

    CHECK(mutation_sites.at(4.0).get_mutations().size() == 1ul);
    CHECK(mutation_sites.at(4.0).get_mutations().at(0)->position == 4.0);

    CHECK(site_positions.at(0ul) == 1.0);
    CHECK(site_positions.at(1ul) == 2.0);
    CHECK(site_positions.at(2ul) == 3.0);
    CHECK(site_positions.at(3ul) == 4.0);
  }
}
