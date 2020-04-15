"""
FILE: math.py
LAST MODIFIED: 24-12-2015
DESCRIPTION: Commonly used math functions

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import numpy as np
import scipy.stats as stats


def norm(v):
    """
    Normalise a vector v
    """
    return np.divide(v, mag(v))


def norms(v):
    """
    Normalise a list of vectors v
    """
    return np.divide(v, mag(v)[:, np.newaxis])


def mag(v):
    return np.sqrt((np.array(v) ** 2.0).sum(-1))


def angle(v1, v2, tol=1e-9):
    x = (v1 * v2).sum(-1) / (mag(v1) * mag(v2))
    if abs(1.0 - x) <= tol:
        return 0.0
    else:
        return np.arccos(x)


def rms(x):
    return np.sqrt((x * x).mean(-1))


def meanConfidenceInterval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = a.shape[0]
    m, se = np.mean(a), stats.sem(a)
    h = se * stats.t.ppf((1 + confidence) / 2., n - 1)
    return m, m - h, m + h


def trimAngle(x):
    return np.mod(x, 2 * np.pi * np.sign(x))
