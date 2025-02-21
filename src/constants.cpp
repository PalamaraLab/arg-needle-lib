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

#include "constants.hpp"

#include <thread>

unsigned anl::get_default_concurrency() {

  static const unsigned default_concurrency = [] {
    // Compute the runtime constant. This lambda is executed only once.
    unsigned result = std::thread::hardware_concurrency();
    return result == 0u ? anl::default_max_tasks : result;
  }();

  return default_concurrency;
}
