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

#ifndef ARG_NEEDLE_LIB_SITE_H
#define ARG_NEEDLE_LIB_SITE_H

#include "mutation.hpp"

/**
 * @class Site
 * @brief Represents a site at which there may be one or more mutations.
 */
class Site {
private:
  /**
   * @brief Vector of mutations present on this site.
   */
  std::vector<Mutation*> mutations;

public:
  /**
   * @brief Add a mutation to this site.
   *
   * @param mutation pointer to the mutation to add to this site.
   */
  void add_mutation(Mutation* mutation);

  /**
   * @brief Get the mutations present at this site.
   *
   * @return the mutations present at this site.
   */
  [[nodiscard]] const std::vector<Mutation*>& get_mutations() const;

  /**
   * @brief Get the position of this site.
   *
   * @return the position of this site.
   */
  [[nodiscard]] arg_real_t get_position() const;
};

#endif // ARG_NEEDLE_LIB_SITE_H
