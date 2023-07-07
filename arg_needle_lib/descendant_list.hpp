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

class DescendantList {

private:
  size_t n;
  static size_t threshold; // default value 64 defined in descendant_list.cpp
  std::vector<int> ordered_ids;
  boost::dynamic_bitset<> db;
  bool using_bitset = false;

public:
  explicit DescendantList(size_t _n);
  DescendantList(size_t _n, int _leaf_id);
  void add(const DescendantList& dList);
  const std::vector<int>& values();
  const boost::dynamic_bitset<>& bitset();
  size_t num_values() const;
  std::string to_string() const;
  std::string to_bitset_string() const;
  std::deque<bool> to_deque_bool() const;
  size_t hash() const;
  static void set_threshold(size_t _threshold);
  static void print_threshold();
  bool operator==(const DescendantList& other) const;
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
