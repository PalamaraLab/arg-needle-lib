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


import sys

if __name__ == "__main__":
    logfile = sys.argv[1]

    calibration_lines = []
    with open(logfile, 'r') as rf:
        for line in rf:
            if "alibration" in line:
                calibration_lines.append(line.strip())

    if len(calibration_lines) == 1:
        # Here we didn't run the LMM so there's no calibration
        print("1")
    elif len(calibration_lines) > 1:
        # Here we get the last calibration line if we ran the LMM
        best_line = calibration_lines[-2].split()
        # We need to invert the value reported by BOLT (value to multiply by)
        # to get the calibration factor (value to divide by)
        print(1 / best_line[-2])
    else:
        # Wrong input
        raise RuntimeError(f"No calibration factor found in {logfile}")
