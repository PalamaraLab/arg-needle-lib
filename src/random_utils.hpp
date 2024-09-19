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

#ifndef __RANDOM_UTILS_HPP_
#define __RANDOM_UTILS_HPP_

#include "types.hpp"

#include <random>

namespace random_utils {

arg_real_t generate_exponential_rv(std::mt19937& generator, arg_real_t rate);
arg_real_t generate_poisson_rv(std::mt19937& generator, arg_real_t rate);
arg_real_t generate_uniform_rv(std::mt19937& generator, arg_real_t from, arg_real_t to);

} // namespace random_utils

#endif // __RANDOM_UTILS_HPP_
