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


"""Utilities for examples.
"""

from contextlib import contextmanager
import time

# https://docs.python.org/3/library/time.html#time.perf_counter
# See https://stackoverflow.com/a/64401876/
default_timer = time.perf_counter

@contextmanager
def time_and_print(description=None):
    """Times a block of code and prints after finishing.

    Inspired by https://stackoverflow.com/a/30024601.
    """
    start = default_timer()
    yield
    end = default_timer()
    time_in_seconds = end - start
    if description is None or description == "":
        print(f"Time elapsed (seconds): {time_in_seconds}")
    else:
        print(f"Time to {description} (seconds): {time_in_seconds}")
