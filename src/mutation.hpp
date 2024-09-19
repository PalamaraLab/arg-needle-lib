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

#ifndef ARG_NEEDLE_LIB_MUTATION_H
#define ARG_NEEDLE_LIB_MUTATION_H

#include "arg_edge.hpp"
#include "types.hpp"

/**
 * @class Mutation
 * @brief Represents a Mutation.
 */
class Mutation {
public:
  /**
   * @brief The position of the mutation.
   */
  double position;

  /**
   * @brief The height of the mutation.
   */
  double height;

  /**
   * @brief The ARGEdge associated with the mutation.
   */
  ARGEdge* edge;

  /**
   * @brief The site id of the mutation.
   */
  int site_id;

  /**
   * @brief Construct a new Mutation object.
   * @param _edge Pointer to an ARGEdge object.
   * @param _position Position of the mutation.
   * @param _height Height of the mutation, default value is -1.0.
   * @param _site_id Site id of the mutation, default value is -1.
   */
  Mutation(ARGEdge* _edge, double _position, double _height = -1.0, int _site_id = -1);

  /**
   * @brief Less than operator, compares the position of two Mutation objects.
   * @param other The other Mutation object to compare with.
   * @return true If position of this object is less than the other.
   * @return false Otherwise.
   */
  bool operator<(const Mutation& other) const {
    return position < other.position;
  }

  /**
   * @brief Less than or equal operator, compares the position of two Mutation objects.
   * @param other The other Mutation object to compare with.
   * @return true If position of this object is less than or equal to the other.
   * @return false Otherwise.
   */
  bool operator<=(const Mutation& other) const {
    return position <= other.position;
  }

  /**
   * @brief Greater than operator, compares the position of two Mutation objects.
   * @param other The other Mutation object to compare with.
   * @return true If position of this object is greater than or equal to the other.
   * @return false Otherwise.
   */
  bool operator>(const Mutation& other) const {
    return position > other.position;
  }

  /**
   * @brief Less than or equal operator, compares the position of two Mutation objects.
   * @param other The other Mutation object to compare with.
   * @return true If position of this object is less than or equal to the other.
   * @return false Otherwise.
   */
  bool operator>=(const Mutation& other) const {
    return position >= other.position;
  }
};

#endif // ARG_NEEDLE_LIB_MUTATION_H
