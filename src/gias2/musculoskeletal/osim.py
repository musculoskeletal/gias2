"""
FILE: osim.py
LAST MODIFIED: 01-07-2016
DESCRIPTION: Module of wrappers and helper functions and classes for opensim
models

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""
import inspect

import numpy as np
import opensim

# import pdb

try:
    opensim_version = getattr(opensim, '__version__')
except AttributeError:
    opensim_version = None


class Body(object):

    def __init__(self, b):
        self._osimBody = b
        self._massScaleFactor = None
        self._inertialScaleFactor = None

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
        return np.array([v.get(i) for i in range(3)])

    @massCenter.setter
    def massCenter(self, x):
        v = opensim.Vec3(x[0], x[1], x[2])
        self._osimBody.setMassCenter(v)

    @property
    def inertia(self):
        m = opensim.Mat33()
        ma = np.zeros((3, 3))
        self._osimBody.getInertia(m)
        for i in range(3):
            for j in range(3):
                ma[i, j] = m.get(i, j)
        return ma

    @inertia.setter
    def inertia(self, I):
        _I = np.array(I)
        if len(_I.shape) == 1:
            inertia = opensim.Inertia(_I[0], _I[1], _I[2])
        else:
            inertia = opensim.Inertia(
                _I[0, 0], _I[1, 1], _I[2, 2],
                _I[0, 1], _I[0, 2], _I[1, 2],
            )
        self._osimBody.setInertia(inertia)

    @property
    def scaleFactors(self):
        v = opensim.Vec3()
        self._osimBody.getScaleFactors(v)
        return np.array([v.get(i) for i in range(3)])

    @scaleFactors.setter
    def scaleFactors(self, s):
        v = opensim.Vec3(s[0], s[1], s[2])
        self._osimBody.scale(v)

    def scale(self, scaleFactors, scaleMass=False):
        v = opensim.Vec3(scaleFactors[0], scaleFactors[1], scaleFactors[2])
        self._osimBody.scale(v, scaleMass)

    def scaleInertialProperties(self, scaleFactors, scaleMass=True):
        v = opensim.Vec3(scaleFactors[0], scaleFactors[1], scaleFactors[2])
        self._osimBody.scaleInertialProperties(v, scaleMass)

    def scaleMass(self, scaleFactor):
        self._osimBody.scaleMass(scaleFactor)

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
        if len(filenames) != nGeoms:
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
        self._isConditionalPathPoint = False
        self._isMovingPathPoint = False
        if p.getConcreteClassName() == 'MovingPathPoint':
            self._osimPathPoint = opensim.MovingPathPoint_safeDownCast(p)
            self._isMovingPathPoint = True
        elif p.getConcreteClassName() == 'ConditionalPathPoint':
            self._osimPathPoint = opensim.ConditionalPathPoint_safeDownCast(p)
            self._isConditionalPathPoint = True
        else:
            self._osimPathPoint = p

    @property
    def name(self):
        return self._osimPathPoint.getName()

    @name.setter
    def name(self, name):
        self._osimPathPoint.setName(name)

    @property
    def location(self):
        return np.array([self._osimPathPoint.getLocationCoord(i) for i in range(3)])

    @location.setter
    def location(self, x):
        self._osimPathPoint.setLocationCoord(0, x[0])
        self._osimPathPoint.setLocationCoord(1, x[1])
        self._osimPathPoint.setLocationCoord(2, x[2])

    @property
    def body(self):
        return Body(self._osimPathPoint.getBody())

    def scale(self, sf):
        raise (NotImplementedError, 'Consider using Muscle.scale.')
        # state = opensim.State()
        # scaleset = opensim.ScaleSet() # ???
        # scaleset.setScale([integer]) #???
        # mus._osimMuscle.scale(state, scaleset)

    @property
    def isMovingPathPoint(self):
        return self._isMovingPathPoint

    @property
    def isConditionalPathPoint(self):
        return self._isConditionalPathPoint

    def _getSimmSpline(self, axis):
        """
        Return the SimmSpline of a given axis (x, y, or z) if self.isMovingPathPoint
        """
        if axis == 'x':
            func = self._osimPathPoint.getXFunction()
        elif axis == 'y':
            func = self._osimPathPoint.getYFunction()
        elif axis == 'z':
            func = self._osimPathPoint.getZFunction()

        ss = opensim.SimmSpline_safeDownCast(func)
        if ss is None:
            raise TypeError('MovingPathPoint function not a simmspline, {} instead'.format(func.getConcreteClassName()))

        return ss

    def _getSimmSplineParams(self, axis):
        ss = self._getSimmSpline(axis)
        ss_x = np.array([ss.getX(i) for i in range(ss.getSize())])
        ss_y = np.array([ss.getY(i) for i in range(ss.getSize())])
        return np.array([ss_x, ss_y])

    def getSimmSplineParams(self):
        """
        Returns the SimmSpline parameters for the x, y, and z coordinates
        of this path point if it is a MovingPathPoint.

        inputs
        ======
        None

        returns
        =======
        x_params : 2 x n ndarray
            Array of SimmSpline parameters of the x coordinates. First
            row contains the x knot values, second row contains the 
            y knot values
        y_params : 2 x n ndarray
            Array of SimmSpline parameters of the y coordinates. First
            row contains the x knot values, second row contains the 
            y knot values
        z_params : 2 x n ndarray
            Array of SimmSpline parameters of the z coordinates. First
            row contains the x knot values, second row contains the 
            y knot values
        """
        if not self.isMovingPathPoint:
            raise TypeError('Not a MovingPathPoint')

        x_params = self._getSimmSplineParams('x')
        y_params = self._getSimmSplineParams('y')
        z_params = self._getSimmSplineParams('z')
        return x_params, y_params, z_params

    def _updateSimmSplineParams(self, axis, params):

        ss = self._getSimmSpline(axis)
        ssLength = ss.getSize()
        x, y = params
        if (len(x) != ssLength) or (len(y) != ssLength):
            raise (
                ValueError(
                    'Input parameters must be of length {}'.format(ssLength)
                )
            )
        for i in range(ssLength):
            ss.setX(i, x[i])
            ss.setY(i, y[i])

    def updateSimmSplineParams(self, x_params=None, y_params=None, z_params=None):
        """
        Update the SimmSpline parameters of the x, y, z coordinates of
        this path point if it is a MovingPathPoint.

        inputs
        ======
        x_params : 2 x n ndarray
            New x and y knot values for the x coordinate spline. Length must
            be the same as the existing spline.
        y_params : 2 x n ndarray
            New x and y knot values for the y coordinate spline. Length must
            be the same as the existing spline.
        z_params : 2 x n ndarray
            New x and y knot values for the z coordinate spline. Length must
            be the same as the existing spline.

        returns
        =======
        None
        """
        if not self.isMovingPathPoint:
            raise TypeError('Not a MovingPathPoint')

        if x_params is not None:
            self._updateSimmSplineParams('x', x_params)
        if y_params is not None:
            self._updateSimmSplineParams('y', y_params)
        if z_params is not None:
            self._updateSimmSplineParams('z', z_params)

    # def removeMultiplierFunction(self):
    #     """
    #     If pathpoint has a multiplierfunction for its X, Y, or Z, function,
    #     replace the multiplierfunction with the function it is multiplying.
    #     """
    #     if not self.isMovingPathPoint:
    #         raise TypeError('Not a MovingPathPoint')

    #     newfunc = self._osimPathPoint.getXFunction()
    #     if newfunc.getConcreteClassName()=='MultiplierFunction':
    #         oldfunc = opensim.MultiplierFunction_safeDownCast(newfunc).getFunction()
    #         owner.setFunction(oldfunc.clone())


class Muscle(object):

    def __init__(self, m):
        self._osimMuscle = m
        self.path_points = {}
        self._init_path_points()

    def _init_path_points(self):
        pps = self.getAllPathPoints()
        for pp in pps:
            self.path_points[pp.name] = pp

    @property
    def name(self):
        return self._osimMuscle.getName()

    @name.setter
    def name(self, name):
        self._osimMuscle.setName(name)

    @property
    def tendonSlackLength(self):
        return self._osimMuscle.getTendonSlackLength()

    @tendonSlackLength.setter
    def tendonSlackLength(self, tsl):
        self._osimMuscle.setTendonSlackLength(tsl)

    @property
    def optimalFiberLength(self):
        return self._osimMuscle.getOptimalFiberLength()

    @optimalFiberLength.setter
    def optimalFiberLength(self, tsl):
        self._osimMuscle.setOptimalFiberLength(tsl)

    def getPathPoint(self, i):
        gp = self._osimMuscle.getGeometryPath()
        pathPoints = gp.getPathPointSet()
        pp = pathPoints.get(i)
        return PathPoint(pp)

    def getAllPathPoints(self):
        pps = []
        gp = self._osimMuscle.getGeometryPath()
        pathPoints = gp.getPathPointSet()
        for i in range(pathPoints.getSize()):
            pp = pathPoints.get(i)
            pps.append(PathPoint(pp))

        return pps

    def preScale(self, state, *scales):
        """
        Scale a pathActuator for a given state by one or
        more Scale instances that define the scale factors
        on the inserted segments
        """
        scaleset = opensim.ScaleSet()
        for sc in scales:
            scaleset.cloneAndAppend(sc._osimScale)

        self._osimMuscle.preScale(state, scaleset)

    def scale(self, state, *scales):
        """
        Scale a pathActuator for a given state by one or
        more Scale instances that define the scale factors
        on the inserted segments
        """
        scaleset = opensim.ScaleSet()
        for sc in scales:
            scaleset.cloneAndAppend(sc._osimScale)

        self._osimMuscle.scale(state, scaleset)

    def postScale(self, state, *scales):
        """
        Scale a pathActuator for a given state by one or
        more Scale instances that define the scale factors
        on the inserted segments
        """
        scaleset = opensim.ScaleSet()
        for sc in scales:
            scaleset.cloneAndAppend(sc._osimScale)

        self._osimMuscle.postScale(state, scaleset)


class CoordinateSet(object):

    def __init__(self, cs):
        self._cs = cs
        self._defaultValue = None

    @property
    def defaultValue(self):
        return self._cs.getDefaultValue()

    @defaultValue.setter
    def defaultValue(self, x):
        self._cs.setDefaultValue(x)


class wrapObject(object):

    def __init__(self, WrObj):
        self._wrapObject = WrObj

    @property
    def name(self):
        return self._wrapObject.getName()

    @name.setter
    def name(self, name):
        self._wrapObject.setName(name)

    def getDimensions(self):
        return self._wrapObject.getDimensionsString()

    def scale(self, scaleFactors):
        v = opensim.Vec3(scaleFactors[0], scaleFactors[1], scaleFactors[2])
        self._wrapObject.scale(v)


class Joint(object):

    def __init__(self, j):
        if j.getConcreteClassName() == 'CustomJoint':
            self._osimJoint = opensim.CustomJoint_safeDownCast(j)
            self._isCustomJoint = True
        else:
            self._osimJoint = j
            self._isCustomJoint = False

        self._initCoordSets()
        if self.isCustomJoint:
            self._initSpatialTransform()
        else:
            self.spatialTransform = None

    def _initCoordSets(self):
        self.coordSets = {}
        cs = self._osimJoint.getCoordinateSet()
        for csi in range(cs.getSize()):
            _cs = cs.get(csi)
            self.coordSets[_cs.getName()] = CoordinateSet(_cs)

    def _initSpatialTransform(self):
        """
        Expose TransformAxes
        """
        self.spatialTransform = self._osimJoint.getSpatialTransform()

    @property
    def name(self):
        return self._osimJoint.getName()

    @name.setter
    def name(self, name):
        self._osimJoint.setName(name)

    def getSimmSplineParams(self, taxisname):
        """
        Returns the SimmSpline parameters for a given TransformAxis.

        inputs
        ======
        taxisname : str
            Name of the TransformAxis

        returns
        =======
        params : 2 x n ndarray
            Array of SimmSpline parameters.
        """

        _method_name = 'get_{}'.format(taxisname)
        _bound_methods = dict(
            inspect.getmembers(
                self.spatialTransform,
                lambda m: inspect.ismethod(m) and m.__func__ in m.im_class.__dict__.values()
            )
        )
        if _method_name not in _bound_methods:
            raise (ValueError('Unknown axis {}'.format(_method_name)))

        tfunc = _bound_methods[_method_name]().get_function()
        ss = opensim.SimmSpline_safeDownCast(tfunc)
        ss_x = np.array([ss.getX(i) for i in range(ss.getSize())])
        ss_y = np.array([ss.getY(i) for i in range(ss.getSize())])
        # ss_z = np.array([ss.getZ(i) for i in range(ss.getSize())])
        return np.array([ss_x, ss_y])

    def updateSimmSplineParams(self, taxisname, x, y):
        """
        Update the SimmSpline parameters for a given TransformAxis

        inputs
        ======
        taxisname : str
            Name of the TransformAxis
        x : 1d ndarray
            New SimmSpline x parameters. Length must be the same as the
            existing x.
        y : 1d ndarray
            New SimmSpline y parameters. Length must be the same as the
            existing y.

        returns
        =======
        None
        """
        _method_name = 'get_{}'.format(taxisname)
        _bound_methods = dict(
            inspect.getmembers(
                self.spatialTransform,
                lambda m: inspect.ismethod(m) and m.__func__ in m.im_class.__dict__.values()
            )
        )
        if _method_name not in _bound_methods:
            raise (ValueError('Unknown axis {}'.format(_method_name)))

        tfunc = _bound_methods[_method_name]().get_function()
        ss = opensim.SimmSpline_safeDownCast(tfunc)
        ssLength = ss.getSize()

        if (len(x) != ssLength) or (len(y) != ssLength):
            raise (
                ValueError(
                    'Input parameters must be of length {}'.format(ssLength)
                )
            )
        for i in range(ssLength):
            ss.setX(i, x[i])
            ss.setY(i, y[i])

    @property
    def isCustomJoint(self):
        return self._isCustomJoint

    @property
    def locationInParent(self):
        v = opensim.Vec3()
        self._osimJoint.getLocationInParent(v)
        return np.array([v.get(i) for i in range(3)])

    @locationInParent.setter
    def locationInParent(self, x):
        v = opensim.Vec3(x[0], x[1], x[2])
        self._osimJoint.setLocationInParent(v)

    @property
    def location(self):
        v = opensim.Vec3()
        self._osimJoint.getLocation(v)
        return np.array([v.get(i) for i in range(3)])

    @location.setter
    def location(self, x):
        v = opensim.Vec3(x[0], x[1], x[2])
        self._osimJoint.setLocation(v)

    @property
    def orientationInParent(self):
        v = opensim.Vec3()
        self._osimJoint.getOrientationInParent(v)
        return np.array([v.get(i) for i in range(3)])

    @orientationInParent.setter
    def orientationInParent(self, x):
        v = opensim.Vec3(x[0], x[1], x[2])
        self._osimJoint.setOrientationInParent(v)

    @property
    def orientation(self):
        v = opensim.Vec3()
        self._osimJoint.getOrientation(v)
        return np.array([v.get(i) for i in range(3)])

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

    def scale(self, *scales):
        """
        Scales joint parameters given one or more Scale instances
        which should define scale factors for joined segments.
        """

        # create ScaleSet
        scaleset = opensim.ScaleSet()
        for sc in scales:
            scaleset.cloneAndAppend(sc._osimScale)

        self._osimJoint.scale(scaleset)


class Scale(object):

    def __init__(self, sfactors=None, name=None, segname=None):

        if len(sfactors) != 3:
            raise (ValueError, 'sfactors must be of length 3')

        self._osimScale = opensim.Scale()
        if sfactors is not None:
            v = opensim.Vec3(
                sfactors[0],
                sfactors[1],
                sfactors[2],
            )
            self._osimScale.setScaleFactors(v)
        if segname is not None:
            self._osimScale.setSegmentName(segname)
        if name is not None:
            self._osimScale.setName(name)
        self._osimScale.setApply(True)

    @property
    def name(self):
        return self._osimScale.getName()

    @name.setter
    def name(self, name):
        self._osimScale.setName(name)

    @property
    def segmentName(self):
        return self._osimScale.getSegmentName()

    @segmentName.setter
    def segmentName(self, name):
        self._osimScale.setSegmentName(name)

    @property
    def scaleFactors(self):
        v = self._osimScale.getScaleFactors()
        return np.array([v.get(i) for i in range(3)])

    @scaleFactors.setter
    def scaleFactors(self, sfactors):
        v = opensim.Vec3(
            sfactors[0],
            sfactors[1],
            sfactors[2],
        )
        self._osimScale.setScaleFactors(v)

    def apply(self, isapply):
        self._osimScale.setApply(isapply)


class Marker(object):
    """
    Pythonic wrap of opensim's Marker class
    """

    def __init__(self, m=None, name=None, bodyname=None, offset=None):
        if m is None:
            self._osimMarker = opensim.Marker()
            self.bodyName = bodyname
            self.offset = offset
            self.name = name
        else:
            self._osimMarker = m

    @property
    def name(self):
        return self._osimMarker.getName()

    @name.setter
    def name(self, name):
        self._osimMarker.setName(str(name))

    @property
    def bodyName(self):
        return self._osimMarker.getBodyName()

    @bodyName.setter
    def bodyName(self, bodyName):
        self._osimMarker.setBodyName(bodyName)

    @property
    def offset(self):
        if opensim_version == 4.0:
            v = self._osimMarker.get_location(v)
        else:
            v = opensim.Vec3()
            self._osimMarker.getOffset(v)
        return np.array([v.get(i) for i in range(3)])

    @offset.setter
    def offset(self, x):
        v = opensim.Vec3(x[0], x[1], x[2])
        if opensim_version == 4.0:
            self._osimMarker.set_location(v)
        else:
            self._osimMarker.setOffset(v)

    # Same as location
    @property
    def offset(self):
        v = opensim.Vec3()
        self._osimMarker.getOffset(v)
        return np.array([v.get(i) for i in range(3)])

    # Same as location
    @offset.setter
    def offset(self, x):
        v = opensim.Vec3(x[0], x[1], x[2])
        self._osimMarker.setOffset(v)


class Model(object):

    def __init__(self, filename=None, model=None):
        self._model = None
        self.joints = {}
        self.bodies = {}
        self.muscles = {}
        self.wrapObjects = {}

        if filename is not None:
            self.load(filename)

        if model is not None:
            self._model = model
            self._init_model()

    def load(self, filename):
        self._model = opensim.Model(filename)
        self._init_model()

    def save(self, filename):
        self._model.printToXML(filename)

    def _init_model(self):
        self._init_joints()
        self._init_bodies()
        self._init_muscles()
        self._init_wrapObjects()

    def _init_joints(self):
        """
        Make a dict of all joints in model
        """
        joints = self._model.getJointSet()
        for ji in range(joints.getSize()):
            j = joints.get(ji)
            self.joints[j.getName()] = Joint(j)

    def _init_bodies(self):
        """
        Make a dict of all bodies in model
        """
        bodies = self._model.getBodySet()
        for bi in range(bodies.getSize()):
            b = bodies.get(bi)
            self.bodies[b.getName()] = Body(b)

    def _init_muscles(self):
        """
        Make a dict of all muscles in body
        """
        muscles = self._model.getMuscles()
        for mi in range(muscles.getSize()):
            m = muscles.get(mi)
            self.muscles[m.getName()] = Muscle(m)

    def _init_wrapObjects(self):
        """
        Make a dict of all wrapping objects in model
        """
        bodies = self._model.getBodySet()
        for bi in range(bodies.getSize()):
            b = bodies.get(bi)
            wObjects = b.getWrapObjectSet()
            if (wObjects.getSize() != 0):
                for wi in range(wObjects.getSize()):
                    w = wObjects.get(wi)
                    self.wrapObjects[w.getName()] = wrapObject(w)

    def scale(self, state, *scales):
        """
        Scale the entire model for a given state and one or more
        Scale instances that define the scale factors for different
        segments
        """

        scaleset = opensim.ScaleSet()
        for sc in scales:
            scaleset.cloneAndAppend(sc._osimScale)

        self._model.scale(state, scaleset)

    def view_init_state(self):
        self._model.setUseVisualizer(True)
        state = self._model.initSystem()
        v = self._model.updVisualizer()
        v.show(state)
        return v
