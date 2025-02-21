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

#include "random_utils.hpp"

#include <random>

namespace random_utils {

arg_real_t generate_exponential_rv(std::mt19937& generator, arg_real_t rate) {
  std::exponential_distribution<arg_real_t> distribution(rate);
  return distribution(generator);
}

arg_real_t generate_poisson_rv(std::mt19937& generator, arg_real_t rate) {
  std::poisson_distribution<int> distribution(rate);
  return distribution(generator);
}

arg_real_t generate_uniform_rv(std::mt19937& generator, arg_real_t from, arg_real_t to) {
  std::uniform_real_distribution<arg_real_t> distribution(from, to);
  return distribution(generator);
}

} // namespace random_utils
