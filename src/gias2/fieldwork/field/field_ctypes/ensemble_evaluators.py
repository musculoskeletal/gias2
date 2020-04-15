"""
FILE: ensemble_evaluators.py
LAST MODIFIED: 24-12-2015
DESCRIPTION: cython field evaluators
    
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
"""

"""
python functions wrapping c code for evaluating elements
"""
import ctypes as C

from numpy import zeros

cEval = C.CDLL('/home/yams/PhD/fieldwork/field/field_ctypes/ensemble_evaluators_c.so')
fpter = C.POINTER(C.c_float)


def eval(B, P):
    l, w = B.shape
    E = zeros(w, dtype=float)

    c_B = (fpter * B.shape[0])(*[row.ctypes.data_as(fpter) for row in B])
    c_P = P.ctypes.data_as(C.POINTER(C.c_double))
    c_E = E.ctypes.data_as(C.POINTER(C.c_double))

    cEval.eval(c_B, c_P, c_E, l, w)
    return E
