"""
FILE: osim.py
LAST MODIFIED: 09-02-2016
DESCRIPTION: Module of wrappers and helper functions and classes for opensim
models

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import opensim
from opensim.opensim import ProbeSet
import numpy as np
import pdb

class Body(object):

    def __init__(self, b):
        self._osimBody = b

    @property
    def name(self):
        return self._osimBody.getName()

    @name.setter
    def name(self, name):
        self._osimBody.setName(name)

    @property
    def mass(self):
        return self._osimBody.getMass()

    @mass.setter
    def mass(self, m):
        self._osimBody.setMass(m)

    @property
    def massCenter(self):
        v = opensim.Vec3()
        self._osimBody.getMassCenter(v)
        return [v.get(i) for i in xrange(3)]

    @massCenter.setter
    def massCenter(self, x):
        v = opensim.Vec3(x[0], x[1], x[2])
        self._osimBody.setMassCenter(v)

    @property
    def inertia(self):
        m = opensim.Mat33()
        ma = np.zeros((3,3))
        self._osimBody.getInertia(m)
        for i in xrange(3):
            for j in xrange(3):
                ma[i,j] = m.get(i,j)
        return ma

    @inertia.setter
    def inertia(self, I):
        inertia = opensim.Inertia(I[0], I[1], I[2])
        self._osimBody.setInertia(inertia)

    def setDisplayGeometryFileName(self, filenames):
        geoset = self._osimBody.getDisplayer().getGeometrySet()
        nGeoms = geoset.getSize()

        # # remove existing geoms
        # for gi in xrange(nGeoms):
        #     geoset.remove(0)

        # # add new geoms
        # for fi, fn in enumerate(filenames):
        #     dgnew = opensim.DisplayGeometry()
        #     dgnew.setGeometryFile(fn)
        #     geoset.insert(fi, dgnew)

        # remove existing geoms
        if len(filenames)!=nGeoms:
            raise ValueError(
                'Expected {} filenames, got {}'.format(
                    nGeoms, len(filenames)
                    )
                )

        # add new geoms
        for fi, fn in enumerate(filenames):
            disp_geo = geoset.get(fi)
            disp_geo.setGeometryFile(fn)

        # if oldfilename is None:
        #     visibles.setGeometryFileName(0, filename)
        # else:
        #     for i in xrange(visibles.getNumGeometryFiles()):
        #         if oldfilename==visibles.getGeometryFileName(i):
        #             visibles.setGeometryFileName(i, filename)

class PathPoint(object):

    def __init__(self, p):
        self._osimPathPoint = p

    @property
    def name(self):
        return self._osimPathPoint.getName()

    @name.setter
    def name(self, name):
        self._osimPathPoint.setName(name)

    @property
    def location(self):
        return [self._osimPathPoint.getLocationCoord(i) for i in xrange(3)]

    @location.setter
    def location(self, x):
        self._osimPathPoint.setLocationCoord(0, x[0])
        self._osimPathPoint.setLocationCoord(1, x[1])
        self._osimPathPoint.setLocationCoord(2, x[2])
    
class Muscle(object):

    def __init__(self, m):
        self._osimMuscle = m

    @property
    def name(self):
        return self._osimMuscle.getName()

    @name.setter
    def name(self, name):
        self._osimMuscle.setName(name)

    def getPathPoint(self, i):
        gp = self._osimMuscle.getGeometryPath()
        pathPoints = gp.getPathPointSet()
        pp = pathPoints.get(i)
        return PathPoint(pp)

    def getAllPathPoints(self):
        pps = []
        gp = self._osimMuscle.getGeometryPath()
        pathPoints = gp.getPathPointSet()
        for i in xrange(pathPoints.getSize()):
            pp = pathPoints.get(i)
            pps.append(PathPoint(pp))

        return pps

class Joint(object):

    def __init__(self, j):
        self._osimJoint = j

    @property
    def locationInParent(self):
        v = opensim.Vec3()
        self._osimJoint.getLocationInParent(v)
        return [v.get(i) for i in xrange(3)]
    
    @locationInParent.setter
    def locationInParent(self, x):
        v = opensim.Vec3(x[0], x[1], x[2])
        self._osimJoint.setLocationInParent(v)

    @property
    def location(self):
        v = opensim.Vec3()
        self._osimJoint.getLocation(v)
        return [v.get(i) for i in xrange(3)]
    
    @location.setter
    def location(self, x):
        v = opensim.Vec3(x[0], x[1], x[2])
        self._osimJoint.setLocation(v)

    @property
    def orientationInParent(self):
        v = opensim.Vec3()
        self._osimJoint.getOrientationInParent(v)
        return [v.get(i) for i in xrange(3)]
    
    @orientationInParent.setter
    def orientationInParent(self, x):
        v = opensim.Vec3(x[0], x[1], x[2])
        self._osimJoint.setOrientationInParent(v)

    @property
    def orientation(self):
        v = opensim.Vec3()
        self._osimJoint.getOrientation(v)
        return [v.get(i) for i in xrange(3)]
    
    @orientation.setter
    def orientation(self, x):
        v = opensim.Vec3(x[0], x[1], x[2])
        self._osimJoint.setOrientation(v)

    @property
    def parentName(self):
        return self._osimJoint.getParentName()

    @parentName.setter
    def parentName(self, name):
        self._osimJoint.setParentName(name)
    
class Model(object):

    def __init__(self, filename=None):
        if filename is not None:
            self.load(filename)

    def load(self, filename):
        self._model = opensim.Model(filename)

    def save(self, filename):
        self._model.printToXML(filename)

    def getBody(self, bodyname):
        bodies = self._model.getBodySet()
        for bi in xrange(bodies.getSize()):
            b = bodies.get(bi)
            if bodyname==b.getName():
                return Body(b)

        raise ValueError('No bodies named {}'.format(bodyname))

    def getJoint(self, jointname):
        joints = self._model.getJointSet()
        for ji in xrange(joints.getSize()):
            j = joints.get(ji)
            if jointname==j.getName():
                return Joint(j)

        raise ValueError('No joints named {}'.format(jointname))

    def getMuscle(self, musclename):
        muscles = model.getMuscles()
        for mi in xrange(muscles.getSize()):
            m = muscles.get(ji)
            if musclename==m.getName():
                return Muscle(m)

        raise ValueError('No muscle named {}'.format(musclename))