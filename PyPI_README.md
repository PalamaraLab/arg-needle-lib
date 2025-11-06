# arg-needle-lib

This repository contains arg-needle-lib, which implements an ARG data structure and ARG-based analyses such as genealogy-wide association.

Prebuilt CPython wheels are available for Linux (compatible with glibc ≥ 2.28) and macOS (built on macOS 15 for x86_64 and macOS 14 for arm64).

| Platform \ CPython          | ≤3.8 | 3.9 | 3.10 | 3.11 | 3.12 | 3.13 | 3.14 |
|-----------------------------| ---- | --- | ---- | ---- | ---- | ---- | ---- |
| Linux x86_64                | ❌    | ✅   | ✅    | ✅    | ✅    | ✅    | ✅    |
| Linux aarch64               | ❌    | ✅   | ✅    | ✅    | ✅    | ✅    | ✅    |
| macOS Intel (x86_64)        | ❌    | ✅   | ✅    | ✅    | ✅    | ✅    | ✅    |
| macOS Apple Silicon (arm64) | ❌    | ✅   | ✅    | ✅    | ✅    | ✅    | ✅    |

## Quickstart

### Install the Python module from PyPI

Most functionality is available through a Python module which can be installed with:

```bash
pip install arg-needle-lib
```

This Python module is currently available on Linux and macOS.

### Documentation

Please see the [ARG-Needle manual](https://palamaralab.github.io/software/argneedle/) for all usage instructions and documentation.

## License

arg-needle-lib is distributed under the GNU General Public License v3.0 (GPLv3). For any questions or comments on arg-needle-lib, please contact Pier Palamara using `<lastname>@stats.ox.ac.uk`.

## Acknowledgements

arg-needle-lib is developed by (in alphabetical order) Arjun Biddanda, Fergus Cooper, Árni Freyr Gunnarsson, Pier Francesco Palamara, Sinan Shi, Brian C. Zhang, and Jiazheng Zhu.

## Reference

If you use this software, please cite:

B. C. Zhang, A. Biddanda, Á. F. Gunnarsson, F. Cooper, P. F. Palamara, Biobank-scale inference of ancestral recombination graphs enables genealogical analysis of complex traits. [Nature Genetics, 2023](https://www.nature.com/articles/s41588-023-01379-x).
