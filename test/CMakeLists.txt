# This file is part of the ARG-Needle genealogical inference and
# analysis software suite.
# Copyright (C) 2023-2024 ARG-Needle Developers.

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

Include(FetchContent)
FetchContent_Declare(
        Catch2
        GIT_SHALLOW TRUE
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG 53d0d913a422d356b23dd927547febdf69ee9081 # version 3.5.0
)
FetchContent_MakeAvailable(Catch2)
message(STATUS "Obtained Catch2 ${Catch2_VERSION}")

set(
        test_src
        test_arg.cpp
        test_arg_utils.cpp
        test_mutation.cpp
        test_serialize_arg.cpp
        test_simple.cpp
        test_structure.cpp
        test_threading.cpp
)

add_executable(unit_tests ${test_src})
target_link_libraries(unit_tests PRIVATE arg_needle_lib Catch2::Catch2WithMain)

# automatically discover tests
list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)
include(CTest)
include(Catch)
catch_discover_tests(unit_tests)
