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

cmake_minimum_required(VERSION 3.16)
message(STATUS "Using CMake version ${CMAKE_VERSION}")
cmake_policy(VERSION 3.16)

project(arg_needle_lib LANGUAGES CXX VERSION 1.1.3)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(PYBIND11_NEWPYTHON ON)

# Project settings including default build type
include(cmake/ProjectSettings.cmake)

# Link this 'library' to use the warnings specified in CompilerWarnings.cmake
add_library(project_warnings INTERFACE)
include(cmake/CompilerWarnings.cmake)
set_project_warnings(project_warnings)

# Sanitiser options if supported by compiler
include(cmake/Sanitisers.cmake)
enable_sanitisers(project_settings)

# allow for static analysis options
include(cmake/StaticAnalysers.cmake)

option(BUILD_SHARED_LIBS "Enable compilation of shared libraries" OFF)

set(arg_needle_testdata_dir ${CMAKE_CURRENT_SOURCE_DIR}/testdata)
add_definitions(-DARG_NEEDLE_TESTDATA_DIR=\"${arg_needle_testdata_dir}\")

option(PYTHON_BINDINGS "Whether to build the python bindings" ON)
if (PYTHON_BINDINGS)
    include(FetchContent)
    FetchContent_Declare(
            pybind11
            GIT_REPOSITORY https://github.com/pybind/pybind11
            GIT_TAG 7c33cdc2d39c7b99a122579f53bc94c8eb3332ff # Version 2.13.5
    )
    FetchContent_GetProperties(pybind11)
    if (NOT pybind11_POPULATED)
        FetchContent_Populate(pybind11)
        add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
    endif ()
endif ()

add_subdirectory(src)
add_subdirectory(example)

option(ENABLE_TESTING "Enable Test Builds" ON)
if (ENABLE_TESTING)
    enable_testing()
    add_subdirectory(test)
endif ()

# Clang format target
find_program(
        CLANG_FORMAT_EXE NAMES "clang-format-16" "clang-format-15" "clang-format-14" "clang-format-13" "clang-format-12"
        "clang-format-11" "clang-format-10" "clang-format" DOC "Path to clang-format executable"
)
if (NOT CLANG_FORMAT_EXE)
    message(STATUS "clang-format not found.")
else ()
    message(STATUS "clang-format found: ${CLANG_FORMAT_EXE}")
    add_custom_target(
            arg_needle_lib_clang_format
            COMMAND ${CLANG_FORMAT_EXE} -style=file -i ${CMAKE_SOURCE_DIR}/src/*pp
            COMMAND ${CLANG_FORMAT_EXE} -style=file -i ${CMAKE_SOURCE_DIR}/example/*pp
            COMMAND ${CLANG_FORMAT_EXE} -style=file -i ${CMAKE_SOURCE_DIR}/test/*pp
    )
endif ()

# Python testing targets
find_package (Python3 3.9 COMPONENTS Interpreter)
if (NOT Python3_Interpreter_FOUND)
    message(STATUS "Python 3.9+ interpreter not found.")
else ()
    message(STATUS "Python ${Python3_VERSION} interpreter found: ${Python3_EXECUTABLE}")
    set (BIG_PYTHON_TEST ${CMAKE_SOURCE_DIR}/test/test_big_args.py)
    add_custom_target(
            pytest
            COMMAND ${Python3_EXECUTABLE} -m pytest ${CMAKE_SOURCE_DIR}/test/ -v -s --ignore=${BIG_PYTHON_TEST}
    )
    add_custom_target(
            pytest_long
            COMMAND ${Python3_EXECUTABLE} -m pytest ${BIG_PYTHON_TEST} -v -s
    )
endif ()
