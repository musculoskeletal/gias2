"""
FILE: osim.py
LAST MODIFIED: 10-09-2021
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


class Body(object):
    """
    Pythonic wrap of OpenSim's Body class (4.2).
    """

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
        v = self._osimBody.getMassCenter()
        return np.array([v.get(i) for i in range(3)])

    @massCenter.setter
    def massCenter(self, x):
        v = opensim.Vec3(x[0], x[1], x[2])
        self._osimBody.setMassCenter(v)

    @property
    def inertia(self):
        ma = np.zeros((3, 3))
        m = self._osimBody.getInertia()

        # SimTK Inertia objects now seem to only output the diagonal of the
        # Inertia matrix rather than the whole matrix. Are the following values
        # sufficient?
        moments = m.getMoments()
        for i in range(3):
            ma[i, i] = moments[i]

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
        v = self._osimBody.getScaleFactors()
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
        """
        This method will currently only work if the parameters are given in the
        same order that they are listed in the model. It will break if the list
        of filenames is larger than the current attached_geometry list.
        """
        for index, file_name in enumerate(filenames):
            geometry = self._osimBody.upd_attached_geometry(index)
            opensim.Mesh.safeDownCast(geometry).set_mesh_file(file_name)


class PathPoint(object):
    """
    This class wraps PathPoint and MovingPathPoint objects.

    MovingPathPoint no longer inherits from PathPoint, we should re-name this
    to AbstractPathPoint to be consistent with the OpenSim 4.2 API.
    """

    def __init__(self, p):
        self._isConditionalPathPoint = False
        self._isMovingPathPoint = False

        if p.getConcreteClassName() == "PathPoint":
            self._osimPathPoint = opensim.PathPoint_safeDownCast(p)
        elif p.getConcreteClassName() == 'MovingPathPoint':
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
        if isinstance(self._osimPathPoint, opensim.simulation.PathPoint):
            v = self._osimPathPoint.get_location()
            return np.array([v[0], v[1], v[2]])
        else:
            return np.array([0.0, 0.0, 0.0])

    @location.setter
    def location(self, x):
        # MovingPathPoints no longer has this attribute, so should be skipped.
        if isinstance(self._osimPathPoint, opensim.simulation.PathPoint):
            v = opensim.Vec3(x[0], x[1], x[2])
            self._osimPathPoint.set_location(v)

    @property
    def body(self):
        return Body(self._osimPathPoint.getParentFrame())

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

    def get_concrete_class_name(self):
        return self._osimPathPoint.getConcreteClassName()

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
        Scale a pathActuator for a given state by one or more Scale instances
        that define the scale factors on the inserted segments.
        """
        scaleset = opensim.ScaleSet()
        for sc in scales:
            scaleset.cloneAndAppend(sc._osimScale)

        self._osimMuscle.preScale(state, scaleset)

    def scale(self, state, *scales):
        """
        Scale a pathActuator for a given state by one or more Scale instances
        that define the scale factors on the inserted segments.
        """
        scaleset = opensim.ScaleSet()
        for sc in scales:
            scaleset.cloneAndAppend(sc._osimScale)

        self._osimMuscle.scale(state, scaleset)

    def postScale(self, state, *scales):
        """
        Scale a pathActuator for a given state by one or more Scale instances
        that define the scale factors on the inserted segments.
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

    def get_translation(self):
        return self._wrapObject.get_translation()

    @name.setter
    def translation(self, translation):
        v = opensim.Vec3(translation)
        self._wrapObject.set_translation(v)

    def getDimensions(self):
        return self._wrapObject.getDimensionsString()

    def scale(self, scaleFactors):
        v = opensim.Vec3(scaleFactors[0], scaleFactors[1], scaleFactors[2])
        self._wrapObject.scale(v)


class Joint(object):
    """
    Pythonic wrap of OpenSim's Joint class (4.2).
    """

    def __init__(self, j):
        if j.getConcreteClassName() == 'CustomJoint':
            self._osimJoint = opensim.CustomJoint_safeDownCast(j)
            self._isCustomJoint = True

        else:
            self._osimJoint = j
            self._isCustomJoint = False

        self._initCoordSets()

        if self._isCustomJoint:
            self._initSpatialTransform()
        else:
            self.spatialTransform = None

    def _initCoordSets(self):
        self.coordSets = {}

        for i in range(self._osimJoint.numCoordinates()):
            _cs = self._osimJoint.get_coordinates(i)
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
                lambda m: inspect.ismethod(m) and hasattr(m, '__self__')
            )
        )
        if _method_name not in _bound_methods:
            raise (ValueError('Unknown axis {}'.format(_method_name)))

        tfunc = _bound_methods[_method_name]().get_function()
        ss = opensim.SimmSpline_safeDownCast(tfunc)
        ss_x = np.array([ss.getX(i) for i in range(ss.getSize())])
        ss_y = np.array([ss.getY(i) for i in range(ss.getSize())])
        # Why aren't we doing the z-axis?
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
                lambda m: inspect.ismethod(m) and hasattr(m, '__self__')
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
        v = self._osimJoint.get_frames(0).get_translation()
        return np.array([v.get(i) for i in range(3)])

    @locationInParent.setter
    def locationInParent(self, x):
        v = opensim.Vec3(x[0], x[1], x[2])
        self._osimJoint.upd_frames(0).set_translation(v)

    @property
    def location(self):
        v = self._osimJoint.get_frames(1).get_translation()
        return np.array([v.get(i) for i in range(3)])

    @location.setter
    def location(self, x):
        v = opensim.Vec3(x[0], x[1], x[2])
        self._osimJoint.upd_frames(1).set_translation(v)

    @property
    def orientationInParent(self):
        v = self._osimJoint.getOrientationInParent()
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

    def scale(self, state, *scales):
        """
        Scales joint parameters given one or more Scale instances
        which should define scale factors for joined segments.
        """

        # create ScaleSet
        scaleset = opensim.ScaleSet()
        for sc in scales:
            scaleset.cloneAndAppend(sc._osimScale)

        self._osimJoint.scale(state, scaleset)


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
    Pythonic wrap of OpenSim's Marker class (4.2).
    """

    def __init__(self, marker=None, name=None, frame_name=None, location=None):
        """
        This constructor can be used in multiple ways. Either the user can
        supply an OpenSim::Marker object, in which case the constructor simply
        wraps the Marker object; or the user can use the constructor to create
        a new Marker object by specifying the name, frame_name and location
        explicitly. name and frame_name should both be strings, location should
        be 3-dimensional tuple of integers (similarly, setting the marker's
        location should be done with a tuple).
        """
        if marker is None:
            self._osim_marker = opensim.Marker()
            self._osim_marker.setName(name)
            self._osim_marker.setParentFrameName(frame_name)
            self._osim_marker.set_location(opensim.Vec3(location))
        else:
            self._osim_marker = marker
            self._check_frame_name()

    def get_osim_marker(self):
        return self._osim_marker

    def _check_frame_name(self):
        if self._osim_marker.getParentFrameName()[:9] == "/bodyset/":
            self._osim_marker.setParentFrameName(
                self._osim_marker.getParentFrameName()[9:]
            )

    @property
    def name(self):
        return self._osim_marker.getName()

    @name.setter
    def name(self, name):
        self._osim_marker.setName(name)

    @property
    def frame_name(self):
        return self._osim_marker.getParentFrameName()

    @frame_name.setter
    def frame_name(self, frame_name):
        self._osim_marker.setParentFrameName(frame_name)

    @property
    def location(self):
        vector = self._osim_marker.get_location()
        return np.array([vector[0], vector[1], vector[2]])

    @location.setter
    def location(self, location):
        self._osim_marker.set_location(opensim.Vec3(location))


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
        print("Successfully Saved.")

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

    def init_system(self):
        return self._model.initSystem()

    def set_marker_set(self, marker_set):
        self._model.set_MarkerSet(marker_set)

    def update_marker_set(self, marker_set):
        self._model.updateMarkerSet(marker_set)

    def scale(self, state, *scales, preserve_mass_distribution, mass):
        """
        Scale the entire model for a given state and one or more
        Scale instances that define the scale factors for different
        segments
        """

        scaleset = opensim.ScaleSet()
        for sc in scales:
            scaleset.cloneAndAppend(sc._osimScale)

        self._model.scale(state, scaleset, preserve_mass_distribution, mass)

    def get_working_state(self):
        return self._model.getWorkingState()

    def get_muscles(self):
        return self._model.getMuscles()

    def view_init_state(self):
        self._model.setUseVisualizer(True)
        state = self._model.initSystem()
        v = self._model.updVisualizer()
        v.show(state)
        return v
