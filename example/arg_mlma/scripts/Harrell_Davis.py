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
import warnings
import numpy as np
import scipy.stats as st
import scipy.stats.mstats as ms

def bootstrap_quant(data, quant=0.95, replicates=1000):
    stats = np.zeros((replicates,1))
    for i in range(replicates):
        boots = np.random.choice(data, len(data))
        stats[i] = ms.hdquantiles(boots, prob=[quant])
    CI = ms.hdquantiles(stats, prob=[0.025, 0.5, 0.975])
    return [CI[0], CI[1], CI[2], len(data)]


def compute_ci(res):
    P_sig = ms.hdquantiles(res, prob=[0.05])
    # If we want the Harrel-Davis SE we can use this instead:
    # P_se = ms.hdquantiles_sd(res, prob=[0.05])
    CI = bootstrap_quant(res, quant=0.05, replicates=1000)
    return [P_sig[0], CI[0], CI[2], CI[3]]


def compute_cis(ps, maxReplicates=1000):
    res = [p for p in ps[0: maxReplicates]]
    assert len(res) > 0
    CI = compute_ci(res)
    return CI, res


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    ps = []
    for p in sys.stdin:
        ps.append(float(p))
    CI, res = compute_cis(ps, maxReplicates=1000)
    print('\t'.join(map(str, CI)))
