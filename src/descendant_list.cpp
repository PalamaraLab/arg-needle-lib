/*
  This file is part of the ARG-Needle genealogical inference and
  analysis software suite.
  Copyright (C) 2023-2024 ARG-Needle Developers.

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

#include "descendant_list.hpp"
#include "utils.hpp"

#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>
#include <iostream>

using boost::dynamic_bitset;
using std::cout;
using std::deque;
using std::endl;
using std::ostream;
using std::string;
using std::vector;

size_t DescendantList::threshold = 64;

DescendantList::DescendantList(size_t _n, int _leaf_id) : n(_n) {
  if (threshold <= 1) {
    using_bitset = true;
    db = dynamic_bitset<>(n);
    db.set(_leaf_id, true);
  }
  else {
    ordered_ids.reserve(1);
    ordered_ids.push_back(_leaf_id);
  }
}

DescendantList::DescendantList(size_t _n) : n(_n) {
  if (threshold <= 0) {
    using_bitset = true;
    db = dynamic_bitset<>(n);
  }
}

DescendantList::DescendantList(size_t _n, const dynamic_bitset<>& bitset) : n(_n), db(bitset) {
  using_bitset = true;
}

DescendantList::DescendantList(size_t _n, const std::vector<int>& bitset) : n(_n) {
  if (bitset.size() != n) {
    throw std::invalid_argument(THROW_LINE("Bitset has wrong size."));
  }

  // Count carriers
  int num_desc = 0;
  for (const int b : bitset) {
    if (b > 1 || b < 0) {
      throw std::invalid_argument(THROW_LINE("Only bitsets with value 0 and 1 are accepted."));
    }
    else if (b == 1) {
      ++num_desc;
    }
  }

  using_bitset = num_desc >= threshold;
  if (using_bitset) {
    // Initialize dynamic bitset
    db = dynamic_bitset<>(n);
    for (int i = 0; i < bitset.size(); i++) {
      if (bitset[i]) {
        db.set(i, true);
      }
    }
  }
  else {
    ordered_ids.reserve(num_desc);
    // Initialize set
    for (int i = 0; i < bitset.size(); i++) {
      if (bitset[i]) {
        ordered_ids.push_back(i);
      }
    }
  }
}

void DescendantList::switch_to_db() {
  if (!using_bitset) {
    db = dynamic_bitset<>(n);
    for (int i : ordered_ids) {
      db.set(i, true);
    }
    using_bitset = true;
    ordered_ids.clear();
  }
}

int DescendantList::peek() const{
  if (num_values() == 0) {
    throw std::runtime_error(THROW_LINE("No elements in descendant list."));
  }
  if (using_bitset) {
    return db.find_first();
  }
  return ordered_ids[0];
}

int DescendantList::get(const int i) const {
  if (i < 0 || i >= n) {
    throw std::invalid_argument(THROW_LINE("Index out of bounds."));
  }
  if (using_bitset) {
    return db[i];
  }
  return 1 ? std::find(ordered_ids.begin(), ordered_ids.end(), i) != ordered_ids.end() : 0;
}

void DescendantList::set(const int i, bool v) {
  if (i < 0 || i >= n) {
    throw std::invalid_argument(THROW_LINE("Index out of bounds."));
  }
  if (using_bitset) {
    db.set(i, v);
  }
  else {
    if (v) {
      ordered_ids.insert(std::upper_bound(ordered_ids.begin(), ordered_ids.end(), i), i);
    }
    else {
      ordered_ids.erase(std::remove(ordered_ids.begin(), ordered_ids.end(), i), ordered_ids.end());
    }
  }
  if (!using_bitset && num_values() >= threshold) {
    switch_to_db();
  }
}

// Add input bitset to this bitset
void DescendantList::add(const DescendantList& other) {
  assert(n == other.n);

  if (using_bitset && other.using_bitset) {
    db |= other.db;
  }
  else if (using_bitset && !other.using_bitset) {
    for (int i : other.ordered_ids) {
      db.set(i, true);
    }
  }
  else if (other.using_bitset || num_values() + other.num_values() >= threshold) {
    db = dynamic_bitset<>(n);
    for (int i : ordered_ids) {
      db.set(i, true);
    }
    if (other.using_bitset) {
      db |= other.db;
    }
    else {
      for (int i : other.ordered_ids) {
        db.set(i, true);
      }
    }
    using_bitset = true;
    ordered_ids.clear();
  }
  else {
    vector<int> tmp;
    tmp.reserve(ordered_ids.size() + other.ordered_ids.size());
    // TODO: remove duplicates if they exist, even though they shouldn't
    std::merge(ordered_ids.begin(), ordered_ids.end(), other.ordered_ids.begin(),
               other.ordered_ids.end(), std::back_inserter(tmp));
    tmp.resize(tmp.end() - tmp.begin());
    ordered_ids = std::move(tmp);
  }

  // If we've passed the threshold, convert container to db.
  if (!using_bitset && num_values() >= threshold) {
    switch_to_db();
  }
}

// Subtract input bitset from this bitset
void DescendantList::erase(const DescendantList& other) {
  assert(n == other.n);
  if (using_bitset && other.using_bitset) {
    db -= other.db;
  }
  else if (using_bitset && !other.using_bitset) {
    for (const int i : other.ordered_ids) {
      db.set(i, false);
    }
  }
  else if (!using_bitset && other.using_bitset) {
    vector<int> tmp;
    for (int i : ordered_ids) {
      if (!other.db[i]) {
        tmp.push_back(i);
      }
    }
    ordered_ids = std::move(tmp);
  }
  else {
    vector<int> tmp;
    std::set_difference(ordered_ids.begin(), ordered_ids.end(), other.ordered_ids.begin(),
                        other.ordered_ids.end(), std::back_inserter(tmp));
    ordered_ids = std::move(tmp);
  }

  // If we've dropped below threshold, convert to set-based DescList
  if (using_bitset && num_values() < threshold) {
    ordered_ids.clear();
    for (size_t i = db.find_first(); i < n; i = db.find_next(i)) {
      ordered_ids.push_back(i);
    }
    db.clear();
    using_bitset = false;
    return;
  }
}

// Returns true iff input bitset is a subset of this bitset
bool DescendantList::includes(DescendantList& other) {
  if (using_bitset && other.using_bitset) {
    return other.db.is_subset_of(db);
  }
  else if (using_bitset && !other.using_bitset) {
    for (int i : other.ordered_ids) {
      if (!db[i]) {
        return false;
      }
    }
  }
  else if (!using_bitset && other.using_bitset) {
    vector<int> intersection;
    intersection.reserve(std::min(ordered_ids.size(), other.num_values()));
    for (int i : ordered_ids) {
      if (other.db[i]) {
        intersection.push_back(i);
      }
    }
    return intersection.size() == other.num_values();
  }
  else {
    vector<int> intersection;
    intersection.reserve(std::min(ordered_ids.size(), other.ordered_ids.size()));
    std::set_intersection(ordered_ids.begin(), ordered_ids.end(), other.ordered_ids.begin(),
                          other.ordered_ids.end(), std::back_inserter(intersection));
    return intersection.size() == other.ordered_ids.size();
  }
  return true;
  // else {
  //   return intersect(other).num_values() == other.num_values();
  // }
}

DescendantList DescendantList::intersect(DescendantList& other) {
  dynamic_bitset<> db_other = other.bitset();
  dynamic_bitset<> intersection(n);
  for (const auto v : values()) {
    if (db_other[v]) {
      intersection.set(v, true);
    }
  }
  return DescendantList(n, intersection);
}

const vector<int>& DescendantList::values() const{
  if (using_bitset) {
    ordered_ids.clear();
    for (size_t index = db.find_first(); index < n; index = db.find_next(index)) {
      ordered_ids.push_back(index);
    }
  }
  return ordered_ids;
}

const dynamic_bitset<>& DescendantList::bitset() {
  if (!using_bitset) {
    db = dynamic_bitset<>(n);
    for (int i : ordered_ids) {
      db.set(i, true);
    }
  }
  return db;
}

DescendantList DescendantList::complement() {
  std::vector<int> complement(n, 1);
  for (auto v : values()) {
    complement[v] = false;
  }
  return DescendantList(n, complement);
}

size_t DescendantList::num_values() const {
  if (using_bitset) {
    return db.count();
  }
  else {
    return ordered_ids.size();
  }
}

string DescendantList::to_string() const {
  if (using_bitset) {
    string s;
    boost::to_string(db, s);
    reverse(s.begin(), s.end());
    return s;
  }
  else {
    string s(n, '0');
    for (int i : ordered_ids) {
      s[i] = '1';
    }
    return s;
  }
}

string DescendantList::to_bitset_string() const {
  if (using_bitset) {
    return this->to_string();
  }
  else {
    dynamic_bitset<> bs(n);
    for (int i : ordered_ids) {
      bs.set(i, true);
    }
    string s;
    boost::to_string(bs, s);
    reverse(s.begin(), s.end());
    return s;
  }
}

deque<bool> DescendantList::to_deque_bool() const {
  deque<bool> result(n, false);
  if (using_bitset) {
    for (size_t i = 0; i < n; ++i) {
      if (db[i]) {
        result[i] = true;
      }
    }
  }
  else {
    for (int i : ordered_ids) {
      result[i] = true;
    }
  }
  return result;
}

// https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector/
// https://stackoverflow.com/a/3897217/ (for boost::dynamic_bitset)
size_t DescendantList::hash() const {
  if (using_bitset) {
    return boost::hash_value(db.m_bits);
  }
  else {
    size_t seed = ordered_ids.size();
    for (auto& i : ordered_ids) {
      seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    // boost::hash_range(seed, ordered_ids.begin(), ordered_ids.end());
    return seed;
  }
}

void DescendantList::set_threshold(size_t _threshold) {
  threshold = _threshold;
}

void DescendantList::print_threshold() {
  cout << "DescendantList threshold is currently " << threshold << endl;
}

bool DescendantList::operator==(const DescendantList& other) const {
  if ((n != other.n) || (using_bitset != other.using_bitset)) {
    return false;
  }
  if (using_bitset) {
    return db == other.db;
  }
  else {
    return ordered_ids == other.ordered_ids;
  }
}

DescendantList DescendantList::operator+(const DescendantList& other) const {
  if (n != other.n) {
    std::invalid_argument(THROW_LINE("Bitsets have different sizes."));
  }
  DescendantList out_list(*this);
  out_list.add(other);
  return out_list;
}

DescendantList DescendantList::operator-(const DescendantList& other) const {
  if (n != other.n) {
    std::invalid_argument(THROW_LINE("Bitsets have different sizes."));
  }
  DescendantList out_list(*this);
  out_list.erase(other);
  return out_list;
}

// looks pretty slow, with many string allocations!
ostream& operator<<(ostream& os, const DescendantList& desc_list) {
  os << "{";
  string subset = "";
  for (auto const& leafID : desc_list.ordered_ids) {
    subset += std::to_string(leafID) + ", ";
  }
  os << subset.substr(0, subset.size() - 2);
  os << "}";
  return os;
}

DescendantListOld::DescendantListOld(int _leaf_id) {
  ids.insert(_leaf_id);
}

DescendantListOld::DescendantListOld() {
}

void DescendantListOld::add(const DescendantListOld& dList) {
  ids.insert(dList.ids.begin(), dList.ids.end());
}

string DescendantListOld::to_string(int num_elements) const {
  string res(num_elements, '0');
  for (int i : ids) {
    res[i] = '1';
  }
  return res;
}

ostream& operator<<(ostream& os, const DescendantListOld& dList) {
  os << "{";
  string subset = "";
  for (auto const& leafID : dList.ids) {
    subset += std::to_string(leafID) + ", ";
  }
  os << subset.substr(0, subset.size() - 2);
  os << "}";
  return os;
}
