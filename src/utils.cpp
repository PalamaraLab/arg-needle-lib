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

#include "utils.hpp"

#include <cassert>
#include <chrono>
#include <cmath>
#include <ctime>
#include <deque>
#include <stdexcept>
#include <string>
#include <vector>

using std::deque;
using std::string;
using std::vector;

namespace utils {

string current_time_string() {
  auto now = std::chrono::system_clock::now();
  std::time_t curr_time = std::chrono::system_clock::to_time_t(now);

  char output[100];
  if (std::strftime(output, sizeof(output), "%Y-%m-%d %X", std::localtime(&curr_time))) {
    string my_string(output);
    return my_string;
  }
  else {
    throw std::runtime_error(THROW_LINE("Trying to get time failed."));
  }
}

int string_to_int(string s) {
  return atoi(s.c_str());
}

int arg_to_int(char* arg) {
  return string_to_int((string) arg);
}

// technically the squared L2 norm
arg_real_t l2(const vector<arg_real_t>& v1, const vector<arg_real_t>& v2) {
  assert(v1.size() == v2.size());
  arg_real_t result = 0;
  for (size_t i = 0; i < v1.size(); ++i) {
    arg_real_t diff = v2[i] - v1[i];
    result += diff * diff;
  }
  return result;
}

arg_real_t r2(const vector<arg_real_t>& v1, const vector<arg_real_t>& v2) {
  assert(v1.size() == v2.size());
  arg_real_t size = static_cast<arg_real_t>(v1.size());

  // std::accumulate requires C++20
  // real v1_mean = std::accumulate(v1.begin(), v1.end(), 0.0) / size;
  // real v2_mean = std::accumulate(v2.begin(), v2.end(), 0.0) / size;

  arg_real_t v1_mean = 0;
  arg_real_t v2_mean = 0;
  for (size_t i = 0; i < size; ++i) {
    v1_mean += v1[i];
    v2_mean += v2[i];
  }
  v1_mean /= size;
  v2_mean /= size;

  arg_real_t cov = 0;
  arg_real_t sd1 = 0;
  arg_real_t sd2 = 0;
  for (size_t i = 0; i < size; ++i) {
    arg_real_t v1_i_centered = v1[i] - v1_mean;
    arg_real_t v2_i_centered = v2[i] - v2_mean;

    sd1 += v1_i_centered * v1_i_centered;
    sd2 += v2_i_centered * v2_i_centered;
    cov += v1_i_centered * v2_i_centered;
  }
  sd1 = sqrt(sd1 / size);
  sd2 = sqrt(sd2 / size);

  // For now, we're not doing anything special if sd1 or sd2 is 0
  cov = cov / sd1 / sd2 / size;

  return cov * cov;
}

arg_real_t relative_ratio(arg_real_t a, arg_real_t b) {
  // designed to emulate Python's math.isclose
  // https://docs.python.org/3/library/math.html#math.isclose
  if (std::abs(a) >= abs(b)) {
    return std::abs(a - b) / std::abs(a);
  }
  else {
    return std::abs(a - b) / std::abs(b);
  }
}

bool is_close(arg_real_t a, arg_real_t b, arg_real_t rel_tol, arg_real_t abs_tol) {
  // designed to emulate Python's math.isclose
  // https://docs.python.org/3/library/math.html#math.isclose
  assert(rel_tol > 0);
  assert(abs_tol >= 0);
  return std::abs(a - b) <= std::max(rel_tol * std::max(std::abs(a), std::abs(b)), abs_tol);
}

vector<arg_real_t> standardize(const vector<arg_real_t>& raw) {
  assert(raw.size() >= 1);
  arg_real_t mean = 0;
  arg_real_t std = 0;
  for (arg_real_t value : raw) {
    mean += value;
  }
  mean /= ((arg_real_t) raw.size());
  for (arg_real_t value : raw) {
    arg_real_t diff = value - mean;
    std += diff * diff;
  }
  std /= ((arg_real_t) raw.size() - 1);
  std = sqrt(std);
  vector<arg_real_t> standardized;
  for (arg_real_t value : raw) {
    arg_real_t new_value = (value - mean) / std;
    standardized.push_back(new_value);
  }
  return standardized;
}

// standardizes only looking at the entries for which use_sample is true
// for entries where use_sample is false, the corresponding value gets set to -9
vector<arg_real_t> standardize_mask(const vector<arg_real_t>& raw, const deque<bool>& use_sample) {
  assert(raw.size() >= 1);
  assert(use_sample.size() == raw.size());
  arg_real_t mean = 0;
  arg_real_t std = 0;
  arg_real_t active_n = 0;
  for (size_t i = 0; i < raw.size(); ++i) {
    if (use_sample[i]) {
      active_n += 1;
      mean += raw[i];
    }
  }
  mean /= active_n;
  for (size_t i = 0; i < raw.size(); ++i) {
    if (use_sample[i]) {
      arg_real_t diff = raw[i] - mean;
      std += diff * diff;
    }
  }
  std /= (active_n - 1);
  std = sqrt(std);
  vector<arg_real_t> standardized;
  for (size_t i = 0; i < raw.size(); ++i) {
    if (use_sample[i]) {
      arg_real_t new_value = (raw[i] - mean) / std;
      standardized.push_back(new_value);
    }
    else {
      standardized.push_back(-9);
    }
  }
  return standardized;
}

} // namespace utils
