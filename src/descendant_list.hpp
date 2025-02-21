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

#ifndef __DESCENDANT_LIST_HPP_
#define __DESCENDANT_LIST_HPP_

#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS // https://stackoverflow.com/a/3897217/

#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <deque>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

/**
 * @brief Represents a list of descendants.
 *
 * This class encapsulates functionality for managing a list of descendants,
 * which can be represented either as a bitset or as an ordered list of integers.
 */
class DescendantList {

private:
  size_t n;
  static size_t threshold; // default value 64 defined in descendant_list.cpp
  mutable std::vector<int> ordered_ids;
  boost::dynamic_bitset<> db;
  bool using_bitset = false;
  void switch_to_db();

public:

  /**
   * @brief Constructor for DescendantList initializing with a specified size.
   * @param _n The size of the DescendantList.
   */
  explicit DescendantList(size_t _n);

  /**
   * @brief Constructor for DescendantList initializing with size and a leaf ID.
   * @param _n The size of the DescendantList.
   * @param _leaf_id The leaf ID to be added to the DescendantList.
   */
  DescendantList(size_t _n, int _leaf_id);

  /**
   * @brief Constructor for DescendantList initializing with size and a vector representing a bitset.
   * @param _n The size of the DescendantList.
   * @param bitset A vector representing a bitset.
   */
  DescendantList(size_t _n, const std::vector<int>& bitset);

  /**
   * @brief Constructor for DescendantList initializing with size and a boost dynamic bitset.
   * @param _n The size of the DescendantList.
   * @param bitset A boost dynamic bitset.
   */
  DescendantList(size_t _n, const boost::dynamic_bitset<>& bitset);

  /**
   * @brief Retrieves the top element of the DescendantList.
   * @return The lowest ID, or index of the first set bit (if using the bitset storage).
   */
  [[nodiscard]] int peek() const;

  int get(const int i) const;
  void set(const int i, bool v);
  void add(const DescendantList& dList);
  void erase(const DescendantList& other);
  bool includes(DescendantList& other);
  DescendantList intersect(DescendantList& other);
  const std::vector<int>& values() const;
  const boost::dynamic_bitset<>& bitset();
  DescendantList complement();
  std::size_t num_values() const;
  std::string to_string() const;
  std::string to_bitset_string() const;
  std::deque<bool> to_deque_bool() const;

  /**
   * @brief Computes a hash value for the DescendantList.
   * @return The computed hash value.
   */
  std::size_t hash() const;

  static void set_threshold(std::size_t _threshold);
  static void print_threshold();

  /**
   * @brief Overloads the equality operator for DescendantList.
   * @param other Another DescendantList to compare with.
   * @return True if both DescendantLists are equal, false otherwise.
   */
  bool operator==(const DescendantList& other) const;

  /**
   * @brief Overloads the addition operator for DescendantList.
   * @param other Another DescendantList to add.
   * @return The resulting DescendantList after addition.
   */
  DescendantList operator+(const DescendantList& other) const;

  /**
   * @brief Overloads the subtraction operator for DescendantList.
   * @param other Another DescendantList to subtract from this list.
   * @return The resulting DescendantList after subtraction.
   */
  DescendantList operator-(const DescendantList& other) const;

  /**
   * @brief Overloads the stream insertion operator for DescendantList.
   * @param os The output stream.
   * @param dList The DescendantList to be inserted into the stream.
   * @return The modified output stream.
   */
  friend std::ostream& operator<<(std::ostream& os, const DescendantList& dList);
};

struct DescendantListHash {
  size_t operator()(const DescendantList& desc_list) const {
    return desc_list.hash();
  }
};

class DescendantListOld {

private:
  std::unordered_set<int> ids;

public:
  DescendantListOld();
  DescendantListOld(int leaf_id);
  void add(const DescendantListOld& dList);
  std::string to_string(int num_elements) const;
  friend std::ostream& operator<<(std::ostream& os, const DescendantListOld& dList);
};

#endif // __DESCENDANT_LIST_HPP_
