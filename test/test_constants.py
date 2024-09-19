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


from arg_needle_lib import ANL_NODE_IS_SAMPLE, ANL_NODE_IS_NOT_SAMPLE

def test_constant_values_and_types():

    assert isinstance(ANL_NODE_IS_SAMPLE, int)
    assert isinstance(ANL_NODE_IS_NOT_SAMPLE, int)

    assert ANL_NODE_IS_SAMPLE == 1
    assert ANL_NODE_IS_NOT_SAMPLE == 0
