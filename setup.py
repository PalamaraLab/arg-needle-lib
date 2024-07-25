# This file is part of the ARG-Needle genealogical inference and
# analysis software suite.
# Copyright (C) 2023 ARG-Needle Developers.

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

# Based on https://github.com/pybind/cmake_example

import os
import platform
import sys
import subprocess

from setuptools import setup, Extension, find_namespace_packages
from setuptools.command.build_ext import build_ext

# Convert distutils Windows platform specifiers to CMake -A arguments
PLAT_TO_CMAKE = {
    "win32": "Win32",
    "win-amd64": "x64",
    "win-arm32": "ARM",
    "win-arm64": "ARM64",
}


# A CMakeExtension needs a sourcedir instead of a file list.
# The name must be the _single_ output extension from the CMake build.
# If you need multiple extensions, see scikit-build.
class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        # required for auto-detection of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        cfg = "Debug" if self.debug else "Release"

        # CMake lets you override the generator - we need to check this.
        # Can be set with Conda-Build, for example.
        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")

        # Set Python_EXECUTABLE instead if you use PYBIND11_FINDPYTHON
        # EXAMPLE_VERSION_INFO shows you how to pass a value into the C++ code
        # from Python.
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",
            "-DWARNINGS_AS_ERRORS=OFF",
            "-DENABLE_TESTING=OFF",
            "-DMAKE_DOCS=OFF"
        ]

        if platform.system() != "Darwin":
            cmake_args.append("-DBoost_NO_BOOST_CMAKE=ON")  # from arni: o/w boost 1.74 gets confused re. mtx

        build_args = []

        if self.compiler.compiler_type != "msvc":
            # Using Ninja-build since it a) is available as a wheel and b)
            # multithreads automatically. MSVC would require all variables be
            # exported for Ninja to pick it up, which is a little tricky to do.
            # Users can override the generator with CMAKE_GENERATOR in CMake
            # 3.15+.
            if not cmake_generator:
                cmake_args += ["-GNinja"]

        else:

            # Single config generators are handled "normally"
            single_config = any(x in cmake_generator for x in {"NMake", "Ninja"})

            # CMake allows an arch-in-generator style for backward compatibility
            contains_arch = any(x in cmake_generator for x in {"ARM", "Win64"})

            # Specify the arch if using MSVC generator, but only if it doesn't
            # contain a backward-compatibility arch spec already in the
            # generator name.
            if not single_config and not contains_arch:
                cmake_args += ["-A", PLAT_TO_CMAKE[self.plat_name]]

            # Multi-config generators have a different way to specify configs
            if not single_config:
                cmake_args += [
                    "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(cfg.upper(), extdir)
                ]
                build_args += ["--config", cfg]

        # Set CMAKE_BUILD_PARALLEL_LEVEL to control the parallel build level
        # across all generators.
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            # self.parallel is a Python 3 only way to set parallel jobs by hand
            # using -j in the build_ext call, not supported by pip or PyPA-build.
            if hasattr(self, "parallel") and self.parallel:
                # CMake 3.12+ only.
                build_args += ["-j{}".format(self.parallel)]

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        subprocess.check_call(
            ["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp
        )
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args, cwd=self.build_temp
        )


with open('PyPI_README.md', encoding='utf-8') as f:
    long_description = f.read()

with open('RELEASE_NOTES.md', encoding='utf-8') as f:
    release_notes = f.read()

setup(
    name='arg-needle-lib',
    version='1.0.2',
    author='ARG-Needle Developers',
    url='https://palamaralab.github.io/software/argneedle/',
    install_requires=[
        'click',
        'h5py',
        'msprime>=1.0.0',
        'numpy<2',
        'pandas',
        'scipy',
        'tskit>=0.1.5',
    ],
    extras_require={
        'dev': [
            'pytest',
            'tszip>=0.2.1',
        ],
        'docs': [
            'sphinx',
            'sphinx-rtd-theme',
        ],
    },
    description='Ancestral recombination graph (ARG) data structure and operations.',
    packages=['arg_needle_lib', 'arg_needle_lib.scripts'],
    long_description='\n'.join([long_description, release_notes]),
    long_description_content_type='text/markdown',
    ext_modules=[CMakeExtension('arg_needle_lib')],
    cmdclass=dict(build_ext=CMakeBuild),
    entry_points = {
        'console_scripts': [
            'arg_association=arg_needle_lib.scripts.association:main',
            'arg_association_prepare_example=arg_needle_lib.scripts.prepare_association_example:main',
            'arg2tskit=arg_needle_lib.scripts.convert:arg2tskit',
            'tskit2arg=arg_needle_lib.scripts.convert:tskit2arg',
        ],
    },
    zip_safe=False,
)
