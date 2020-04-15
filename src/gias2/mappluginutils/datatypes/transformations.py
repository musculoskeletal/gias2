"""
FILE: transformations.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: Geometric transformations classes for mapclient plugins.

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import numpy as np


class Transform(object):
    transformType = None

    def __init__(self, T, P=None):
        self.T = T
        self.P = P

    def getT(self):
        return np.array(self.T)

    def setT(self, T):
        self.T = np.array(T)

    def getP(self):
        return np.array(self.P)

    def setP(self, P):
        self.P = np.array(P)


class AffineTransform(Transform):
    transformType = 'affine'


class RigidTransform(Transform):
    transformType = 'rigid'


class RigidTransformAboutPoint(Transform):
    transformType = 'rigid_about_point'


class RigidScaleTransform(Transform):
    transformType = 'rigidscale'


class RigidScaleTransformAboutPoint(Transform):
    transformType = 'rigidscale_about_point'


class RigidPCModesTransform(Transform):
    transformType = 'rigidpcmodes'
