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

# Find boost
find_package(Boost COMPONENTS iostreams REQUIRED)

# Find Eigen
find_package(Eigen3 CONFIG REQUIRED)

# Find the appropriate system threads library
find_package(Threads REQUIRED)

# Find HDF5
find_package(HDF5 COMPONENTS CXX REQUIRED)

set(
        arg_needle_lib_src
        ancestor_entry.cpp
        arg.cpp
        arg_edge.cpp
        arg_node.cpp
        arg_utils.cpp
        constants.cpp
        descendant_list.cpp
        file_utils.cpp
        genotype_mapping.cpp
        mutation.cpp
        random_utils.cpp
        root.cpp
        serialize_arg.cpp
        site.cpp
        utils.cpp
)

set(
        arg_needle_lib_hdr
        ancestor_entry.hpp
        arg.hpp
        arg_edge.hpp
        arg_node.hpp
        arg_utils.hpp
        constants.hpp
        descendant_list.hpp
        deserialization_params.hpp
        file_utils.hpp
        genotype_mapping.hpp
        mutation.hpp
        IntervalTree.h
        random_utils.hpp
        root.hpp
        serialize_arg.hpp
        site.hpp
        types.hpp
        utils.hpp
        arg_traversal.hpp)

add_library(arg_needle_lib STATIC ${arg_needle_lib_src} ${arg_needle_lib_hdr})

target_include_directories(
        arg_needle_lib INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:>
)

set_target_properties(arg_needle_lib PROPERTIES PUBLIC_HEADER "${arg_needle_lib_hdr}")


target_link_libraries(arg_needle_lib PRIVATE Boost::headers Boost::iostreams Threads::Threads HDF5::HDF5)
target_link_libraries(arg_needle_lib PUBLIC Eigen3::Eigen)

# Turning off warnings for now
# target_link_libraries(arg_needle_lib PRIVATE project_warnings)

# This seems necessary for macOS...
target_include_directories(arg_needle_lib PUBLIC ${Boost_INCLUDE_DIRS})

if (PYTHON_BINDINGS)
    set_target_properties(arg_needle_lib PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
    pybind11_add_module(arg_needle_lib_pybind arg_needle_lib_pybind.cpp)
    target_link_libraries(arg_needle_lib_pybind PRIVATE arg_needle_lib)
    install(TARGETS arg_needle_lib_pybind LIBRARY DESTINATION .)
endif ()
