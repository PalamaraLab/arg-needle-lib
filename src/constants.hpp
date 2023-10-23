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

#ifndef ARG_NEEDLE_LIB_CONSTANTS_H
#define ARG_NEEDLE_LIB_CONSTANTS_H

namespace anl {

    /**
     * The default max max for potentially concurrent tasks, used if the number of cores cannot be determined at
     * runtime.
     */
    constexpr unsigned default_max_tasks = 16u;

    /**
     * Calculate the default number of potentially concurrent tasks for utilities that may be run asynchronously.
     *
     * This is one less than the number of cores. If the number of cores can't be calculated, this function returns
     * one less than default_max_tasks. If this this would be 0, then 1 is returned instead.
     *
     * The result is cached to ensure the computation is performed at most once each time arg-needle-lib is run.
     *
     * @return the default number of potentially concurrent tasks
     */
    unsigned get_default_concurrency();

}

#endif // ARG_NEEDLE_LIB_CONSTANTS_H
