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

#ifndef ARG_NEEDLE_LIB_ROOT_H
#define ARG_NEEDLE_LIB_ROOT_H

#include "arg_node.hpp"
#include "types.hpp"

/**
 * @class Root
 * @brief Represents a Root.
 */
class Root {
public:
  /**
   * @brief Pointer to an ARGNode object.
   */
  ARGNode* node;

  /**
   * @brief The start position.
   */
  arg_real_t start;

  /**
   * @brief The end position.
   */
  arg_real_t end;

  /**
   * @brief Construct a new Root object.
   *
   * @param _node Pointer to an ARGNode object.
   * @param _start The start position.
   * @param _end The end position.
   */
  Root(ARGNode* _node, arg_real_t _start, arg_real_t _end);
};

#endif // ARG_NEEDLE_LIB_ROOT_H
