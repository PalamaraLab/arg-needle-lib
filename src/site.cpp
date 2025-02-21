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

#include "site.hpp"

#include "mutation.hpp"
#include "utils.hpp"

void Site::add_mutation(Mutation* mutation) {
  mutations.push_back(mutation);
}

const std::vector<Mutation*>& Site::get_mutations() const {
  return mutations;
}


arg_real_t Site::get_position() const {

  if (mutations.empty()) {
    throw std::logic_error(THROW_LINE("Site has no mutations so does not have an associated position..."));
  }

  return mutations.front()->position;
}
