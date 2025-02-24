# arg-needle-lib

A C++ library for representing ARGs, including threading operations and ARG-based association, with a Python API exposed through [`pybind11`](https://github.com/pybind/pybind11).
See [this paper](https://www.nature.com/articles/s41588-023-01379-x) for details.

Binaries for the latest release of `arg-needle-lib` can be obtained by running:
```
pip install arg-needle-lib
```

## Documentation

Please see the [ARG-Needle manual](https://palamaralab.github.io/software/argneedle/) for all usage instructions and documentation.
This document is mainly intended to describe compiling the `arg-needle-lib` library from the source files in this repository.

## Local pip install

Locally pip installing `arg-needle-lib` will build the C++ library from source files and install the Python bindings.
Building the library requires you to have `boost` installed.
As long as boost is installed in a standard location (e.g. by your system package manager using `brew install boost`, `sudo apt-get install libboost-all-dev`, `aptitude search boost`, etc), it will be found automatically.
The python bindings rely on [pybind11](https://pybind11.readthedocs.io/en/latest/) which is obtained automatically: you do not need to have this installed.

Follow these steps to perform the pip install:
```
python3 -m venv venv
source venv/bin/activate

pip install --upgrade pip setuptools wheel
pip install .
```

If you are a developer you may want to substitute the last line with
```
pip install ".[dev]"
```
which will install some additional dependencies such as `pytest` and `sphinx`.
You should rerun the pip install for any local changes to be built.

## Python tests

To run tests, after installing the package with pip (see above), run the following from root:
```
mkdir build && cd build
cmake ..
make pytest
```
There is also a test with big ARGs that takes around an hour to run, which can be run using `make pytest_long` instead of the last command.

## Compiling and Dependencies

To directly build the C++ library (without installing the Python API) you must have `boost` installed as above.

The `arg-needle-lib` library uses [CMake](https://cmake.org/).
From the `arg_needle_lib` directory, configure and build the project:

```
mkdir build && cd build
cmake ..
cmake --build . --parallel 4
```

This will build the library, all examples and, by default, the unit tests.
You can then run the unit tests:

```
ctest
```

or an example:

```
./example/example_arg
```

## For developers: making a release

- Bump the version number in [pyproject.toml](pyproject.toml), [CMakeLists.txt](CMakeLists.txt), and [docs/conf.py](docs/conf.py)
- Update [RELEASE_NOTES.md](RELEASE_NOTES.md)
- Push changes and check that all [GitHub workflows](https://github.com/PalamaraLab/arg_needle_lib/actions) pass
- Tag the commit in Git using syntax `vX.Y.Z`
- Make a release on GitHub, which should trigger a new build that will upload Python wheels to PyPI
- If building wheels for a new Python version:
  - Update classifiers in [pyproject.toml](pyproject.toml)
  - Check whether version of `pypa/cibuildwheel` needs updating in [build-wheels.yml](.github/workflows/build-wheels.yml)

## License

`arg-needle-lib` is distributed under the GNU General Public License v3.0 (GPLv3). For any questions or comments on `arg-needle-lib`, please contact Pier Palamara using `<lastname>@stats.ox.ac.uk`.

## Acknowledgements

`arg-needle-lib` is developed by (in alphabetical order) Arjun Biddanda, Fergus Cooper, Árni Freyr Gunnarsson, Pier Francesco Palamara, Sinan Shi, Brian C. Zhang, and Jiazheng Zhu.

The file `src/IntervalTree.h` is copied from https://github.com/ekg/intervaltree/tree/e8082c74a6f5c18de99d8b4cc4a55e2e62a1150d, developed by Erik Garrison and released under the MIT License.

The files `src/file_utils.hpp` and `src/file_utils.cpp` are adapted from analogous files in the Eagle software, which can be found at https://github.com/poruloh/Eagle/tree/master/src. Eagle was developed by Po-Ru Loh and released under the GNU General Public License v3.0 (GPLv3).

The file `src/arg_utils.cpp` contains a block of code that is adapted from the BOLT-LMM_v2.3.2 software. BOLT-LMM v2.3.2 was developed by Po-Ru Loh and released under the GNU General Public License v3.0 (GPLv3).

All third-party licenses can be found under the `3rd_party/` directory.
