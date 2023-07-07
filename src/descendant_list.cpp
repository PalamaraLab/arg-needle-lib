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

#include "descendant_list.hpp"

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
  if (1 >= threshold) {
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
  if (0 >= threshold) {
    using_bitset = true;
    db = dynamic_bitset<>(n);
  }
}

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
  else {
    if (other.using_bitset || num_values() + other.num_values() >= threshold) {
      using_bitset = true;
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
    }
    else {
      vector<int> tmp;
      tmp.reserve(ordered_ids.size() + other.ordered_ids.size());
      // TODO: remove duplicates if they exist, even though they shouldn't
      std::merge(ordered_ids.begin(), ordered_ids.end(), other.ordered_ids.begin(),
                 other.ordered_ids.end(), std::back_inserter(tmp));
      ordered_ids = std::move(tmp);
    }
  }
}

const vector<int>& DescendantList::values() {
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

// looks pretty slow, with many string allocations
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
