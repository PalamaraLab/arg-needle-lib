"""
Check that all files used in building project have consistent version number
"""

import re

def _read_version_triplet(filename, regex):
    with open(filename) as file:
        for line in file:
            regex = re.compile(regex)
            match = re.search(regex, line)
            if match:
                major = int(match[1])
                minor = int(match[2])
                patch = int(match[3])
                triplet = (major, minor, patch)
                print(f" - {filename}: {triplet}")
                return triplet

    raise ValueError(f"Cannot find valid version triplet in {filename}")


if __name__ == "__main__":
    print("Checking project version files are consistent...")

    cmake_version = _read_version_triplet("CMakeLists.txt", "project\(.* (\d+)\.(\d+)\.(\d+)\)")
    docs_version = _read_version_triplet("docs/conf.py", "release.*v(\d+)\.(\d+)\.(\d+)")
    pyproject_version = _read_version_triplet("pyproject.toml", "version.*\"(\d+)\.(\d+)\.(\d+)")

    if cmake_version != docs_version:
        raise ValueError(f"CMake/docs mismatch: {cmake_version}/{docs_version}")

    if cmake_version != pyproject_version:
        raise ValueError(f"CMake/pyproject mismatch: {cmake_version}/{pyproject_version}")

    print(f"All versioned files match.")
    print("Done!")
