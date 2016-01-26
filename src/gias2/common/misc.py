"""
FILE: misc.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: miscellaneous functions

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import sys
import numpy as np

def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """

    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out

class ProgressOutput:
    
    def __init__(self, task, total):
        self.task = task
        self.total = total
        self.value = 0
        
    def progress(self, value, comment=''):
        self.value = value        
        percent = int(100.*value/self.total)
        outcomment = ''
        if len(comment)>0:
            outcomment = ": {}".format(comment)
        sys.stdout.write("Progress: {} : {:2d}\% {}  \r".format(self.task, percent, outcomment))
        sys.stdout.flush()
        
    def output(self, comment=''):
        outcomment = ''
        if len(comment)>0:
            outcomment = ": {}".format(comment)
        sys.stdout.write("\nOutput: {} {}\n".format(self.task, outcomment))
        sys.stdout.flush()