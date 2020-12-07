"""
FILE: geoprimitives.py
LAST MODIFIED: 24-12-2015
DESCRIPTION: geometric primitives like lines and planes

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""
import logging
from typing import Tuple, Union, List, Optional, Iterable

import numpy as np
from scipy.linalg import eig, inv
from scipy.optimize import leastsq, fmin
from scipy.spatial.distance import euclidean

from gias2.common import transform3D
from gias2.common.math import norm, mag, angle

""" 
functions and classes for fitting meshes to segmented data
"""

log = logging.getLogger(__name__)
PRECISION = 1e-16


# ======================================================================#
class NonInterceptError(Exception):
    pass


class InPlaneError(Exception):
    pass


class Line3D(object):
    """
    A line in 3D described by x = t*a + b
    """

    def __init__(self, a: Iterable, b: Iterable):
        """
        A line in 3D defined by its origin `b` and direction `a`
        :param a: length 3 array of line direction vector
        :param b: length 3 array of origin coordinates
        """
        self.a = None
        self.b = None
        self._a = None
        self.setAB(a, b)
        self.t0 = 0.0
        self.t1 = 1.0
        self._l = self

    def setAB(self, a: Iterable, b: Iterable) -> None:
        self.a = norm(np.array(a, dtype=float))
        self._a = self.a[:, np.newaxis]
        self.b = np.array(b, dtype=float)

    def eval(self, t: Union[float, np.ndarray]) -> np.ndarray:
        return (t * self._a).T.squeeze() + self.b

    def findClosest(self, p: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """ calc closest point on line to p
        """
        p = np.array(p, dtype=float)
        closest_t = np.dot((p - self.b), self.a)
        closest_p = self.eval(closest_t)
        return closest_p, closest_t

    def calcDistanceFromPoint(self, p: np.ndarray, retT: bool = False) -> Union[Tuple[float, np.ndarray], float]:
        """ calc closest distance to point p
        """
        p_line, t_line = self.findClosest(p)
        if retT:
            return euclidean(p, p_line), t_line
        else:
            return euclidean(p, p_line)

    def project_coordinate(self, x: np.ndarray, dim: int) -> float:
        """
        Calculate the parameter coordinate along the line that gives a point
        with a coordinate of x in dimension dim
        """
        t = (x - self.b[dim]) / self.a[dim]
        return t

    def checkParallel(self, line: 'Line3D') -> bool:
        """
        check if self is parallel to line l
        """
        u = self.a
        v = line.a
        a = np.dot(u, u)
        b = np.dot(u, v)
        c = np.dot(v, v)
        denom = (a * c - b * b)

        # if lines are parallel
        if denom < PRECISION:
            return True
        else:
            return False

    def checkCoincidence(self, line: 'Line3D') -> bool:
        """check if self is coincident to another infinite 3D line l.
        Checks by seeing if there are 2 common points
        """
        p1 = self.eval(0.0)
        p2 = self.eval(10.0)
        d1 = line.calcDistanceFromPoint(p1)
        d2 = line.calcDistanceFromPoint(p2)
        if (d1 < PRECISION) and (d2 < PRECISION):
            log.debug('coincident')
            return True
        else:
            log.debug('not coincident')
            return False

    def calcIntercept(self, line: 'Line3D') -> Tuple[np.ndarray, float, float]:
        """
        tries to calculate the intercept with line l

        If there is an intercept, returns
            - 3D coordinates of the point of intersection
            - parameter on self at the point of intersection
            - parameter on the input line at the point of intersection

        If there is no intercept, raises NonInterceptError
        """

        # c = self.a[0]/self.a[1]
        # Ix = l.b[0] - self.b[0]

        # # calculate parameters at intercept
        # ti = (Ix*(c-1))/(l.a[0]-c*l.a[1])
        # si = (l.a[0]*ti + Ix)/self.a[0]

        ex = line.b[0] - self.b[0]
        ey = self.b[1] - line.b[1]
        f = self.a[1] / self.a[0]

        ti = (f * ex + ey) / (line.a[1] - line.a[0] * f)
        si = (line.a[0] * ti + ex) / self.a[0]

        # check if parameter actually give an intercept
        p1 = self.eval(si)
        p2 = line.eval(ti)

        if np.sqrt(((p1 - p2) ** 2.0).sum()) < PRECISION:
            return (p1 + p2) / 2.0, si, ti
        else:
            raise NonInterceptError

    def calcClosestDistanceToLine(self, line: 'Line3D') -> Tuple[float, float, float]:
        """ Calculates the closest distance to another infinite 3D line.

        returns:
        d - closest distance
        sc - parameter of self at closest approach
        tc - parameter of input line at closest approach

        reference:
        http://geomalgorithms.com/a07-_distance.html
        """

        # first check if intersecting:
        try:
            pi, sc, tc = self.calcIntercept(line)
            d = 0.0
            log.debug('Intersecting lines...')
        except NonInterceptError:
            log.debug('Non intersecting lines...')
            u = self.a
            v = line.a
            w0 = self.b - line.b
            a = np.dot(u, u)
            b = np.dot(u, v)
            c = np.dot(v, v)
            d = np.dot(u, w0)
            e = np.dot(v, w0)
            denom = (a * c - b * b)
            if denom < 0.0:
                raise RuntimeError('negative denominator in closest approach calculation: {}'.format(denom))

            # if lines are parallel
            if denom < PRECISION:
                sc = self.t0
                d, tc = line.calcDistanceFromPoint(self.eval(sc), retT=True)
            else:
                sc = (b * e - c * d) / denom
                tc = (a * e - b * d) / denom

                wc = w0 + u * sc - v * tc
                d = np.sqrt((wc * wc).sum())

        return d, sc, tc

    def intersectSphere(self, o: np.ndarray, r: float) -> Tuple[List[np.ndarray], List[float]]:
        """
        Calculate intersection with a sphere with centre o and
        radius r. Returns a list of intercept points, and a list of 
        intercept distances. 

        Inputs:
        o: 3d coords of the sphere centre
        r: sphere radius

        Returns:
        xCoords: list of intercept coordinates, length of either 0,
            1, or 2 depending on number of intercepts.
        xDists: list of intercept distances, length of either 0,
            1, or 2 depending on number of intercepts.
        """
        TOL = 1e-12
        r = float(r)
        o = np.array(o, dtype=float)
        ob = self.b - o
        root = np.dot(self.a, ob) ** 2.0 - mag(ob) ** 2.0 + r ** 2.0
        if root < 0.0:
            return [], []
        else:
            const = -np.dot(self.a, ob)
            if abs(root) < TOL:
                return [self.eval(const), ], [const, ]
            else:
                d1 = const + np.sqrt(root)
                d2 = const - np.sqrt(root)
                return [self.eval(d1), self.eval(d2)], [d1, d2]

    def transformAffine(self, t: np.ndarray) -> None:
        self.b = np.dot(
            t,
            np.hstack([
                self.b,
                1.0
            ])
        )[:3]
        self.a = np.dot(
            t[:3, :3],
            self.a
        )


class LineOutOfBoundsError(Exception):
    pass


class LineSegment3D(Line3D):

    def __init__(self, a: Iterable, b: Iterable, t0: float, t1: float):
        """
        A line with finite length defined by `t0` and `t1` which are the distance
        of the two end points from `b` in direction `a`
        :param a: direction of line
        :param b: origin of line
        :param t0: end point 1 distance from b
        :param t1: end point 2 distance from b
        """
        self.t0 = t0
        self.t1 = t1
        self._l = Line3D(a, b)
        self.setAB(a, b)
        # super(LineSegment3D, self).__init__(a,b)
        self.p0 = self.eval(t0)
        self.p1 = self.eval(t1)

    def setAB(self, a: Iterable, b: Iterable) -> None:
        self.a = norm(np.array(a, dtype=float))
        self._a = self.a[:, np.newaxis]
        self.b = np.array(b, dtype=float)
        self.p0 = self.eval(self.t0)
        self.p1 = self.eval(self.t1)
        self._l.setAB(a, b)

    def _checkBound(self, t: float) -> bool:
        return (self.t0 <= t) & (t <= self.t1)

    def checkCoincidence(self, line_segment: 'LineSegment3D') -> bool:
        """check if self is coincident to another 3D line segment l.
        Checks by seeing if there are 2 common points
        """
        return self._l.checkCoincidence(line_segment._l)

    def eval(self, t: float, checkBound: bool = True) -> np.ndarray:
        """
        Evaluate the position on the line at distance `t` from its origin
        :param t:
        :param checkBound: if true, raises LineOutOfBoundsError if t is out beyond the line end points
        :return: length 3 array of position on line
        """

        x = (t * self._a).T.squeeze() + self.b
        if checkBound:
            if self._checkBound(t):
                return x
            else:
                raise LineOutOfBoundsError
        else:
            return x

    def findClosest(self, p: np.ndarray) -> Tuple[np.ndarray, float]:
        p = np.array(p, dtype=float)
        tc = float(np.dot((p - self.b), self.a))

        if isinstance(tc, float):
            if self._checkBound(tc):
                return self.eval(tc), tc
            elif tc < self.t0:
                return self.p0, self.t0
            else:
                return self.p1, self.t1
        else:
            # array
            bd = self._checkBound(tc)
            output_eval = self.eval(tc, checkBound=False)
            output_eval[tc < self.t0] = self.p0
            output_eval[tc > self.t1] = self.p1
            tc = float(np.clip(tc, min=self.t0, max=self.t1))
            return output_eval, tc

    def calcClosestDistanceToLine(self, line: Union[Line3D, 'LineSegment3D']) -> Tuple[float, float, float]:
        if isinstance(line, Line3D):
            return self._distanceToLine(line)
        elif isinstance(line, LineSegment3D):
            return self._distanceToLineSegment(line)

    def _distanceToLine(self, line: Line3D) -> Tuple[float, float, float]:
        """ closest distance to an infinite line
        """
        d, sc, tc = self._l.calcClosestDistanceToLine(line._l)
        if self._checkBound(sc):
            pass
        elif sc < self.t0:
            d, tc = line.calcDistanceFromPoint(self.p0, retT=True)
            sc = self.t0
        else:
            d, tc = line.calcDistanceFromPoint(self.p1, retT=True)
            sc = self.t1

        return d, sc, tc

    def _distanceToLineSegment(self, line: 'LineSegment3D') -> Tuple[float, float, float]:
        """ closest distance to another line segment
        """

        d, sc, tc = self._l.calcClosestDistanceToLine(line._l)

        if self._checkBound(sc):
            if line._checkBound(tc):
                # if both sc and tc are in bound, then all good
                pass
            else:
                # closest point is on self segment
                if tc < line.t0:
                    d, sc = self.calcDistanceFromPoint(line.p0, retT=True)
                    tc = line.t0
                else:
                    d, sc = self.calcDistanceFromPoint(line.p1, retT=True)
                    tc = line.t1
        else:
            if line._checkBound(tc):
                # closest point is on l
                if sc < self.t0:
                    d, sc = line.calcDistanceFromPoint(self.p0, retT=True)
                    sc = self.t0
                else:
                    d, sc = line.calcDistanceFromPoint(self.p1, retT=True)
                    sc = self.t1
            else:
                # closest point is beyond the limits of both segments
                # find the segment with end furtherest from closest point
                ds0 = abs(sc - self.t0)
                ds1 = abs(sc - self.t1)
                dt0 = abs(sc - line.t0)
                dt1 = abs(sc - line.t1)

                if ds0 < ds1:
                    if dt0 < dt1:
                        if ds0 < dt0:
                            # self.p0 is closest, then l.p0
                            tc = line.t0
                            d, sc = self.calcDistanceFromPoint(line.p0, retT=True)
                        else:
                            # l.p0 is closest, then self.p0
                            sc = self.t0
                            d, tc = line.calcDistanceFromPoint(self.p0, retT=True)
                    else:
                        if ds0 < dt1:
                            # self.p0 is closest, then l.p1
                            tc = line.t1
                            d, sc = self.calcDistanceFromPoint(line.p1, retT=True)
                        else:
                            # l.p1 is closest, then self.p0
                            sc = self.t0
                            d, tc = line.calcDistanceFromPoint(self.p0, retT=True)
                else:
                    if dt0 < dt1:
                        if ds1 < dt0:
                            # self.p1 is closest, then l.p0
                            tc = line.t0
                            d, sc = self.calcDistanceFromPoint(line.p0, retT=True)
                        else:
                            # l.p0 is closest, then self.p1
                            sc = self.t1
                            d, tc = line.calcDistanceFromPoint(self.p1, retT=True)
                    else:
                        if ds1 < dt1:
                            # self.p1 is closest, then l.p1
                            tc = line.t1
                            d, sc = self.calcDistanceFromPoint(line.p1, retT=True)
                        else:
                            # l.p1 is closest, then self.p1
                            sc = self.t1
                            d, tc = line.calcDistanceFromPoint(self.p1, retT=True)

        return d, sc, tc


class LineElement3D:
    """ X = (1-x)A + (x)B
    """

    def __init__(self, p1: np.ndarray, p2: np.ndarray):
        self.p1 = np.array(p1)
        self.p2 = np.array(p2)

    def eval(self, x):
        return (1 - x) * self.p1 + x * self.p2


class Plane(object):

    def __init__(
            self,
            origin: np.ndarray,
            normal: np.ndarray,
            x: Optional[np.ndarray] = None,
            y: Optional[np.ndarray] = None):
        """
        A plane in 3D
        :param origin: length 3 array of the origin point of the plane
        :param normal: length 3 array of the normal vector of the plane
        :param x: length 3 array of the X-axis in the plane
        :param y: length 3 array of the Y-axis of the plane
        """

        self.O = np.array(origin, dtype=float)
        self.N = norm(normal)
        self.X = None
        self.Y = None
        if x is not None:
            self.X = np.array(x, dtype=float)
        if y is not None:
            self.Y = np.array(y, dtype=float)

    @property
    def origin(self) -> np.ndarray:
        return self.O

    @property
    def normal(self) -> np.ndarray:
        return self.N

    @property
    def x_axis(self) -> Optional[np.ndarray]:
        return self.X

    @property
    def y_axis(self) -> Optional[np.ndarray]:
        return self.Y

    def calcDistanceToPlane(self, pts: np.ndarray) -> float:
        """
            Signed perpendicular distance from a point `pts` to this plane.
            Distance is positive if `pts` is on the side of the point the normal points to
        """
        d = ((pts - self.O) * self.N).sum(-1)
        return d

    def near_points(self, pts: np.ndarray, dmax: float) -> np.ndarray:
        dist = self.calcDistanceToPlane(pts)
        mask = abs(dist) <= dmax
        return np.array(pts[mask])

    def project2Plane3D(self, points: np.ndarray) -> np.ndarray:
        """
        returns the closest points to P on the plane, in 3D coordinates
        """
        d = self.calcDistanceToPlane(points)
        p = points - d * self.N
        return p

    def project2Plane2D(self, points: np.ndarray) -> np.ndarray:
        """
        returns the closest points to P on the plane, in 2D in-plane
        coordinates
        """
        if (self.X is None) or (self.Y is None):
            raise ValueError('plane X and Y vectors not set')

        u = ((points - self.O) * self.X).sum(-1)
        v = ((points - self.O) * self.Y).sum(-1)
        return np.array([u, v]).T

    def plane2Dto3D(self, points: np.ndarray) -> np.ndarray:
        """
        return 3D coordinates of 2D in-plane coordinates
        """
        if (self.X is None) or (self.Y is None):
            raise ValueError('plane X and Y vectors not set')

        p = points[:, 0, np.newaxis] * self.X + points[:, 1, np.newaxis] * self.Y + self.O
        return p

    def angleToPlane(self, plane: 'Plane') -> float:
        """
        calculates the angle between self normal and the normal
        or another plane p
        """
        return angle(self.N, plane.N)

    def angleToVector(self, vector: np.ndarray) -> float:
        """ calculate the angle between this plane and a vector v
        """
        # project vector onto plane
        v_proj = self.project2Plane3D(vector)

        # calc angle between v and v_proj
        return angle(v_proj, vector)

    def intersect_line(self, line: Line3D, ret_t: bool = False) -> Union[Tuple[np.ndarray, float], np.ndarray]:
        """
        Find the point of intersection between a line l and this plane
        """

        nom = np.dot((self.O - line.b), self.N)
        denom = np.dot(line.a, self.N)

        if abs(nom) < PRECISION:
            # line is in plane
            raise InPlaneError('line is in plane, infinite intersections')
        elif abs(denom) < PRECISION:
            # line is parallel to plane
            raise NonInterceptError('line is parallel to plane but out of plane, no intersections')
        else:
            t = nom / denom
            p = line.eval(t)
            if ret_t:
                return p, t
            else:
                return p

    def transformAffine(self, t: np.ndarray) -> None:
        self.O = np.dot(
            t,
            np.hstack([
                self.O,
                1.0
            ])
        )[:3]
        self.N = norm(
            np.dot(
                t[:3, :3],
                self.N
            )
        )
        if self.X is not None:
            self.X = norm(np.dot(t[:3, :3], self.X))
        if self.Y is not None:
            self.Y = norm(np.dot(t[:3, :3], self.Y))

    def drawPlane(self, mscene, l=100, acolor=(1, 0, 0), ascale=10.0, scolor=(0, 1, 0), sopacity=0.5):
        """ Draw the plane in a mayavi scene as a square and a normal vector arrow.
        Inputs:
        ========
        mscene: mayavi scene
        l: float, square side length
        acolor: 3-tuple, color of the normal arrow in (r,g,b)
        ascale: float, scale factor of the normal arrow
        scolor: 3-tuple, color of the plane in (r,g,b)
        sopacity: float, opacity of square, between 0 and 1.

        Returns:
        =========
        arrow: mayavi output for drawing the plane normal at the plane origin
        square: mayavi output for drawing the plane square
        """
        # plane data to draw
        l2 = 0.5 * l
        planeCorners2D = np.array([[-l2, -l2],
                                   [l2, -l2],
                                   [-l2, l2],
                                   [l2, l2]])
        p3D = self.plane2Dto3D(planeCorners2D)
        square = mscene.mlab.mesh([[p3D[0, 0], p3D[1, 0]],
                                   [p3D[2, 0], p3D[3, 0]]],
                                  [[p3D[0, 1], p3D[1, 1]],
                                   [p3D[2, 1], p3D[3, 1]]],
                                  [[p3D[0, 2], p3D[1, 2]],
                                   [p3D[2, 2], p3D[3, 2]]],
                                  color=scolor,
                                  opacity=sopacity)
        if self.X is None:
            arrow = mscene.mlab.quiver3d(
                [self.O[0]],
                [self.O[1]],
                [self.O[2]],
                [self.N[0]],
                [self.N[1]],
                [self.N[2]],
                mode='arrow',
                scale_factor=ascale,
                color=acolor)
        else:
            arrow = mscene.mlab.quiver3d(
                [self.O[0], self.O[0], self.O[0]],
                [self.O[1], self.O[1], self.O[1]],
                [self.O[2], self.O[2], self.O[2]],
                [self.X[0], self.Y[0], self.N[0], ],
                [self.X[1], self.Y[1], self.N[1], ],
                [self.X[2], self.Y[2], self.N[2], ],
                mode='arrow',
                scale_factor=ascale,
                scalars=[1, 2, 3],
                scale_mode='scalar',
                color=acolor)

        return arrow, square


# ===============================================================================#
def fitAxis3D(data: np.ndarray, axis: Line3D) -> Tuple[Line3D, np.ndarray, float]:
    xtol = 1e-5
    ftol = 1e-5
    maxfev = 6 * 1000
    data_com = data.mean(0)

    def obj(x):
        axis.setAB(x[0:3], x[3:6])
        axis_points = axis.findClosest(data)[0]
        # ~ axisPoints = np.array([axis.findClosest(d)[0] for d in data])
        com_dist = ((axis.b - data_com) ** 2.0).sum()
        d2 = ((data - axis_points) ** 2.0).sum(1)
        return np.hstack([d2, com_dist])

    x_opt = leastsq(obj, np.hstack([axis.a, axis.b]), xtol=xtol, ftol=ftol, maxfev=maxfev)[0]
    fitted_rmse = np.sqrt(obj(x_opt).mean())
    axis.setAB(x_opt[0:3], x_opt[3:6])

    return axis, x_opt, fitted_rmse


def fitPlaneLS(points: np.ndarray) -> Plane:
    # calc CoM
    com = points.mean(0)
    xc = points - com

    # eigen system
    A = np.zeros([3, 3])
    A[0, 0] = (xc[:, 0] ** 2.0).sum()
    A[1, 1] = (xc[:, 1] ** 2.0).sum()
    A[2, 2] = (xc[:, 2] ** 2.0).sum()
    A[0, 1] = A[1, 0] = (xc[:, 0] * xc[:, 1]).sum()
    A[0, 2] = A[2, 0] = (xc[:, 0] * xc[:, 2]).sum()
    A[1, 2] = A[2, 1] = (xc[:, 1] * xc[:, 2]).sum()

    w, v = eig(A)  # the right matrix v is returned by default
    v = v[:, w.argsort()]

    # project points onto plane
    plane = Plane(com, v[:, 0], v[:, 1], v[:, 2])

    return plane


def fitSphere(points: np.ndarray) -> Tuple[float, np.ndarray]:
    def obj(x):
        d = points - x[0:3]
        e = np.sqrt((d * d).sum(1)) - x[3]
        return e * e

    init_centre = np.mean(points, 0)
    init_r = np.mean(mag(points - init_centre))
    x0 = np.hstack((init_centre, init_r))
    try:
        x_opt = leastsq(obj, x0)[0]
    except TypeError:
        log.warning('probably not enough points to fit')
        return 0, np.array([0, 0, 0, 0])
    else:
        rms_opt = np.sqrt((obj(x_opt)).mean())
        return rms_opt, x_opt


def fitSphereAnalytic(points: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    ADAPTED FROM MATLAB SCRIPT:
    
    this fits a sphere to a collection of data using a closed form for the
    solution (opposed to using an array the size of the data set). 
    Minimizes Sum((x-xc)^2+(y-yc)^2+(z-zc)^2-r^2)^2
    x,y,z are the data, xc,yc,zc are the sphere's center, and r is the radius

    Assumes that points are not in a singular configuration, real numbers, ...
    if you have coplanar data, use a circle fit with svd for determining the
    plane, recommended Circle Fit (Pratt method), by Nikolai Chernov
    http://www.mathworks.com/matlabcentral/fileexchange/22643

    Input:
    X: n x 3 matrix of cartesian data
    Outputs:
    Center: Center of sphere 
    Radius: Radius of sphere
    Author:
    Alan Jennings, University of Dayton
    """
    a_matrix = np.zeros((3, 3), dtype=float)
    a_matrix[0, 0] = (points[:, 0] * (points[:, 0] - points[:, 0].mean())).mean()
    a_matrix[0, 1] = 2 * (points[:, 0] * (points[:, 1] - points[:, 1].mean())).mean()
    a_matrix[0, 2] = 2 * (points[:, 0] * (points[:, 2] - points[:, 2].mean())).mean()
    a_matrix[1, 0] = 0
    a_matrix[1, 1] = (points[:, 1] * (points[:, 1] - points[:, 1].mean())).mean()
    a_matrix[1, 2] = 2 * (points[:, 1] * (points[:, 2] - points[:, 2].mean())).mean()
    a_matrix[2, 0] = 0
    a_matrix[2, 1] = 0
    a_matrix[2, 2] = (points[:, 2] * (points[:, 2] - points[:, 2].mean())).mean()

    a_matrix = a_matrix + a_matrix.T

    b_matrix = np.zeros(3, dtype=float)
    b_matrix[0] = ((points[:, 0] ** 2.0 + points[:, 1] ** 2 + points[:, 2] ** 2) * (points[:, 0] - points[:, 0].mean())).mean()
    b_matrix[1] = ((points[:, 0] ** 2.0 + points[:, 1] ** 2 + points[:, 2] ** 2) * (points[:, 1] - points[:, 1].mean())).mean()
    b_matrix[2] = ((points[:, 0] ** 2.0 + points[:, 1] ** 2 + points[:, 2] ** 2) * (points[:, 2] - points[:, 2].mean())).mean()

    # Center=(A\B).';
    centre = np.dot(inv(a_matrix), b_matrix)

    radius = np.sqrt((np.vstack([points[:, 0] - centre[0],
                                 points[:, 1] - centre[1],
                                 points[:, 2] - centre[2]]) ** 2).sum(0).mean())

    return centre, radius


def fitBox(data: np.ndarray, centre: np.ndarray, axes: List) -> Tuple[np.ndarray, float, np.ndarray, List[Line3D]]:
    max_it = 10000
    # initialise axes
    x_line = Line3D(axes[0], centre)
    y_line = Line3D(axes[1], centre)
    z_line = Line3D(axes[2], centre)

    def obj(x):
        # update box axes
        new_b = x[:3]
        old_as = np.array([v.a.copy() for v in [x_line, y_line, z_line]])
        new_as = transform3D.transformRigid3D(old_as, np.array([0, 0, 0, x[3], x[4], x[5]]))

        x_line.setAB(new_as[0], new_b)
        y_line.setAB(new_as[1], new_b)
        z_line.setAB(new_as[2], new_b)

        # project data points
        p_x = x_line.findClosest(data)[1]
        p_y = y_line.findClosest(data)[1]
        p_z = z_line.findClosest(data)[1]

        # calc volume
        return (p_x.max() - p_x.min()) * (p_y.max() - p_y.min()) * (p_z.max() - p_z.min())

    x0 = np.hstack([centre, [0, 0, 0]])
    x_opt = fmin(obj, x0, maxiter=max_it)

    final_centre = x_opt[:3]
    final_volume = obj(x_opt)

    # calculate fitted box dimensions
    p_x = x_line.findClosest(data)[1]
    p_y = y_line.findClosest(data)[1]
    p_z = z_line.findClosest(data)[1]
    final_dim = np.array([p_x.max() - p_x.min(), p_y.max() - p_y.min(), p_z.max() - p_z.min()])

    return final_centre, final_volume, final_dim, [x_line, y_line, z_line]


def circumcentre3Points(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> Tuple[np.ndarray, float, np.ndarray]:
    """
    calculate the circum centre of 3 points and the circle radius,
    also calculates the normal of the circle plane.
    """
    a = np.array(a, dtype=float)
    b = np.array(b, dtype=float)
    c = np.array(c, dtype=float)

    # plane normal
    normal = norm(np.cross((a - b), (a - c)))

    # midpoints of ab and ac
    d = 0.5 * (a + b)
    e = 0.5 * (a + c)

    # calculate perpendicular bisectors
    do = norm(np.cross(normal, d - a))
    eo = norm(np.cross(normal, e - a))

    ldo = Line3D(do, d)
    leo = Line3D(eo, e)

    # find intercept of bisectors
    dist, sldo, sleo = ldo.calcClosestDistanceToLine(leo)
    origin = 0.5 * (ldo.eval(sldo) + leo.eval(sleo))

    # calc radius
    radius = float(np.mean([mag(a - origin), mag(b - origin), mag(c - origin)]))

    return origin, radius, normal
