# This file is part of the ARG-Needle genealogical inference and
# analysis software suite.
# Copyright (C) 2023-2025 ARG-Needle Developers.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


add_executable(example_arg example_arg.cpp)
target_link_libraries(example_arg PRIVATE arg_needle_lib)

add_executable(example_mutation example_mutation.cpp)
target_link_libraries(example_mutation PRIVATE arg_needle_lib)

add_executable(example_node example_node.cpp)
target_link_libraries(example_node PRIVATE arg_needle_lib)

add_executable(example_visit example_visit.cpp)
target_link_libraries(example_visit PRIVATE arg_needle_lib)
