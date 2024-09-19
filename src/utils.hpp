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

#ifndef ARG_NEEDLE_LIB_UTILS_H
#define ARG_NEEDLE_LIB_UTILS_H

#include "types.hpp"

#include <deque>
#include <map>
#include <memory>
#include <string>
#include <vector>


// Utility for exceptions
#define THROW_LINE(a) (std::string(__FILE__) + ":" + std::to_string(__LINE__) + ": " + a)

namespace utils {

std::string current_time_string();

int string_to_int(std::string s);

int arg_to_int(char* arg);

// technically the squared L2 norm
arg_real_t l2(const std::vector<arg_real_t>& v1, const std::vector<arg_real_t>& v2);

arg_real_t r2(const std::vector<arg_real_t>& v1, const std::vector<arg_real_t>& v2);

arg_real_t relative_ratio(arg_real_t a, arg_real_t b);

bool is_close(arg_real_t a, arg_real_t b, arg_real_t rel_tol = 1e-9, arg_real_t abs_tol = 0);

std::vector<arg_real_t> standardize(const std::vector<arg_real_t>& raw);

std::vector<arg_real_t> standardize_mask(const std::vector<arg_real_t>& raw, const std::deque<bool>& use_sample);

// work in progress, not working yet
// template <typename T>
// const typename map<real, T>::iterator map_floor(map<real, T> container, real position) {
//   const typename map<real, T>::iterator it = container.upper_bound(position);
//   // auto map_floor(map<real, T> container, real position) {
//   //   auto it = container.upper_bound(position);
//   if (it == container.begin()) {
//     return container.end();
//   }
//   it = std::prev(it);
//   return it;
// }

} // namespace utils

#endif // ARG_NEEDLE_LIB_UTILS_H
