# arg-needle-lib Release Notes

## v1.1.0 (2024-09-19)

### Major changes

- Added methods for genotype mapping
  - map_genotype_to_ARG, taking a single genotype as a vector/list and a position
  - map_genotypes_to_ARG, taking a matrix of genotypes, and a vector/list of positions
- Method to get sorted vector/list of positions from ARG object is renamed to get_site_positions from get_sites 

### Other changes

- Improve documentation
- Deserialization now performed in C++ rather than Python
- Python infrastructure modernized to replace setup.py with pyproject.toml


## v1.0.2 (2023-09-29)

### Breaking changes

None

### Other changes

- Improve documentation.
- Build Python wheels for macOS arm64 (Apple Silicon)
- Build wheels for Python 3.12


## v1.0.1 (2023-07-14)

### Breaking changes

None

### Other changes

- Clean-up and release of the source code.
- Improved mutation class.
- Added from-to parameters to some functions.


## v1.0.0 (2023-03-07)

Initial PyPI release of arg-needle-lib.
