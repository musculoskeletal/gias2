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

""" 
functions and classes for fitting meshes to segmented data
"""
import os
import sys
import time
import pickle
import itertools
import cProfile

import scipy
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist, euclidean
from scipy.optimize import leastsq, fmin
from scipy.linalg import eig, inv

from gias2.common import transform3D
from gias2.common.math import *

PRECISION = 1e-16

#======================================================================#
class NonInterceptError(Exception):
    pass

class Line3D( object ):
    """ X = t*a + b
    """
    def __init__(self, a, b ):
        self.setAB( a, b )
        self.t0 = 0.0
        self.t1 = 1.0
        self._l = self
        
    def setAB( self, a, b ):
        self.a = norm( scipy.array(a, dtype=float) )
        self._a = self.a[:,scipy.newaxis]
        self.b = scipy.array(b, dtype=float)
        
    def eval( self, t ):
        return (t*self._a).T.squeeze() + self.b
        # try:
        #   return self.a*t[:,scipy.newaxis] + self.b
        # except IndexError:
        #   return self.a*t + self.b

    def findClosest( self, p ):
        """ calc closest point on line to p
        """
        p = scipy.array(p, dtype=float)
        closestT = scipy.dot( (p-self.b), self.a )
        closestP = self.eval( closestT )
        return closestP, closestT
        
    def calcDistanceFromPoint( self, p, retT=False ):
        """ calc closest distance to point p
        """
        pLine, tLine = self.findClosest( p )
        #~ return scipy.sqrt( ( ( p - pLine )**2.0 ).sum() )
        if retT:
            return euclidean( p, pLine ), tLine
        else:
            return euclidean( p, pLine )

    def checkParallel(self, l):
        """
        check if self is parallel to l
        """
        u = self.a
        v = l.a
        a = scipy.dot(u, u)
        b = scipy.dot(u, v)
        c = scipy.dot(v, v)
        denom = (a*c - b*b)

        # if lines are parallel
        if denom < PRECISION:
            return True
        else:
            return False

    def checkCoincidence(self, l):
        """check if self is coincident to another infinite 3D line l.
        Checks by seeing if there are 2 common points
        """
        p1 = self.eval(0.0)
        p2 = self.eval(10.0)
        d1 = l.calcDistanceFromPoint(p1)
        d2 = l.calcDistanceFromPoint(p2)
        # print d1, d2
        if (d1<PRECISION) and (d2<PRECISION):
            # print 'coincident'
            return True
        else:
            # print 'not coincident'
            return False

    def calcIntercept(self, l):
        """
        tries to calculate the intercept with line l
        """

        # c = self.a[0]/self.a[1]
        # Ix = l.b[0] - self.b[0]

        # # calculate parameters at intercept
        # ti = (Ix*(c-1))/(l.a[0]-c*l.a[1])
        # si = (l.a[0]*ti + Ix)/self.a[0]

        ex = l.b[0] - self.b[0]
        ey = self.b[1] - l.b[1]
        f = self.a[1] / self.a[0]

        ti = (f*ex + ey)/(l.a[1] - l.a[0]*f)
        si = (l.a[0]*ti + ex)/self.a[0]

        # check if parameter actually give an intercept
        p1 = self.eval(si)
        p2 = l.eval(ti)

        if np.sqrt(((p1-p2)**2.0).sum())<PRECISION:
            return (p1+p2)/2.0, si, ti
        else:
            # print si, p1
            # print ti, p2
            raise NonInterceptError

    def calcClosestDistanceToLine(self, l):
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
            pi, sc, tc = self.calcIntercept(l)
            d = 0.0
            # print 'Intersecting lines...'
        except NonInterceptError:
            # print 'Non intersecting lines...'
            u = self.a
            v = l.a
            w0 = self.b - l.b
            a = scipy.dot(u, u)
            b = scipy.dot(u, v)
            c = scipy.dot(v, v)
            d = scipy.dot(u, w0)
            e = scipy.dot(v, w0)
            denom = (a*c - b*b)
            if denom < 0.0:
                raise RuntimeError('negative denominator in closest approach calculation: {}'.format(denom))

            # if lines are parallel
            if denom < PRECISION:
                sc = self.t0
                d, tc = l.calcDistanceFromPoint(self.eval(sc), retT=True)
                if self.checkCoincidence(l):
                    print('WARNING: coincident lines')
                else:
                    print('WARNING: parallel lines')
            else:
                sc = (b*e - c*d)/denom
                tc = (a*e - b*d)/denom

                wc = w0 + u*sc - v*tc
                d = scipy.sqrt((wc*wc).sum())

        return d, sc, tc

    def intersectSphere(self, o, r):
        """Calculate intersection with a sphere with centre o and
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
        o = scipy.array(o, dtype=float)
        ob = self.b - o
        root = scipy.dot(self.a, ob)**2.0 - mag(ob)**2.0 + r**2.0
        if root<0.0:
            return [], []
        else:
            const = -np.dot(self.a, ob)
            if abs(root)<TOL:
                return [self.eval(const),], [const,]
            else:
                d1 = const+scipy.sqrt(root)
                d2 = const-scipy.sqrt(root)
                return [self.eval(d1),self.eval(d2)], [d1, d2]

    def transformAffine(self, t):
        self.b = scipy.dot(
                    t,
                    scipy.hstack([
                        self.b,
                        1.0
                        ])
                    )[:3]
        self.a = scipy.dot(
                    t[:3,:3],
                    self.a
                    )

class LineOutOfBoundsError(Exception):
    pass

class LineSegment3D(Line3D):

    def __init__(self, a, b, t0, t1 ):
        super(LineSegment3D, self).__init__(a,b)
        self.t0 = t0
        self.t1 = t1
        self.p0 = self.eval(t0)
        self.p1 = self.eval(t1)
        self._l = Line3D(a,b)

    def _checkBound(self, t):
        return (self.t0<=t) & (t<=self.t1)

    def checkCoincidence(self, l):
        """check if self is coincident to another 3D line segment l.
        Checks by seeing if there are 2 common points
        """
        return self._l.checkCoincidence(l._l)

    def eval(self, t, checkBound=True):

        x = (t*self._a).T.squeeze() + self.b
        if checkBound:
            if self._checkBound(t):
                return x
            else:
                raise LineOutOfBoundsError
        else:
            return x

    def findClosest(self, p):
        p = scipy.array(p, dtype=float)
        tc = scipy.dot( (p-self.b), self.a )

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
            output_eval[tc<self.t0] = self.p0
            output_eval[tc>self.t1] = self.p1
            tc = tc.clip(min=self.t0, max=self.t1)
            return output_eval, tc

    def calcClosestDistanceToLine(self, l):
        if l.__class__==Line3D:
            return self._distanceToLine(l)
        elif l.__class__==LineSegment3D:
            return self._distanceToLineSegment(l)
        
    def _distanceToLine(self, l):
        """ closest distance to an infinite line
        """
        d, sc, tc = self._l.calcClosestDistanceToLine(l._l)
        if self._checkBound(sc):
            pass
        elif sc < self.t0:
            d, tc = l.calcDistanceFromPoint(self.p0, retT=True)
            sc = self.t0
        else:
            d, tc = l.calcDistanceFromPoint(self.p1, retT=True)
            sc = self.t1

        return d, sc, tc

    def _distanceToLineSegment(self, l):
        """ closest distance to another line segment
        """

        d, sc, tc = self._l.calcClosestDistanceToLine(l._l)
    
        if self._checkBound(sc):
            if l._checkBound(tc):
                # if both sc and tc are in bound, then all good
                pass 
            else:
                # closest point is on self segment
                if tc<l.t0:
                    d, sc = self.calcDistanceFromPoint(l.p0, retT=True)
                    tc = l.t0
                else:
                    d, sc = self.calcDistanceFromPoint(l.p1, retT=True)
                    tc = l.t1
        else:
            if l._checkBound(tc):
                # closest point is on l
                if sc<self.t0:
                    d, sc = l.calcDistanceFromPoint(self.p0, retT=True)
                    sc = self.t0
                else:
                    d, sc = l.calcDistanceFromPoint(self.p1, retT=True)
                    sc = self.t1
            else:
                # closest point is beyond the limits of both segments
                # find the segment with end furtherest from closest point
                ds0 = abs(sc - self.t0)
                ds1 = abs(sc - self.t1)
                dt0 = abs(sc - l.t0)
                dt1 = abs(sc - l.t1)

                if ds0 < ds1:
                    if dt0 < dt1:
                        if ds0 < dt0:
                            # self.p0 is closest, then l.p0
                            tc = l.t0
                            d, sc = self.calcDistanceFromPoint(l.p0, retT=True)
                        else:
                            # l.p0 is closest, then self.p0
                            sc = self.t0
                            d, tc = l.calcDistanceFromPoint(self.p0, retT=True)
                    else:
                        if ds0 < dt1:
                            # self.p0 is closest, then l.p1
                            tc = l.t1
                            d, sc = self.calcDistanceFromPoint(l.p1, retT=True)
                        else:
                            # l.p1 is closest, then self.p0
                            sc = self.t0
                            d, tc = l.calcDistanceFromPoint(self.p0, retT=True)
                else:
                    if dt0 < dt1:
                        if ds1 < dt0:
                            # self.p1 is closest, then l.p0
                            tc = l.t0
                            d, sc = self.calcDistanceFromPoint(l.p0, retT=True)
                        else:
                            # l.p0 is closest, then self.p1
                            sc = self.t1
                            d, tc = l.calcDistanceFromPoint(self.p1, retT=True)
                    else:
                        if ds1 < dt1:
                            # self.p1 is closest, then l.p1
                            tc = l.t1
                            d, sc = self.calcDistanceFromPoint(l.p1, retT=True)
                        else:
                            # l.p1 is closest, then self.p1
                            sc = self.t1
                            d, tc = l.calcDistanceFromPoint(self.p1, retT=True)

        return d, sc, tc

class LineElement3D(Line3D):
    """ X = (1-x)A + (x)B
    """
    def __init__( self, p1, p2):
        self.p1 = scipy.array(p1)
        self.p2 = scipy.array(p2)
        
    def eval( self, x ):
        return (1-x)*self.p1 + x*self.p2

class Plane( object ):
    
    def __init__( self, origin, normal, x=None, y=None ):
        
        self.O = scipy.array(origin, dtype=float)
        self.N = norm(normal)
        self.X = scipy.array(x, dtype=float)
        self.Y = scipy.array(y, dtype=float)
    
    def calcDistanceToPlane( self, P ):
        d = ((P - self.O) * self.N ).sum(-1)
        return d
        
    def project2Plane3D( self, P ):
        """
        returns the closest points to P on the plane, in 3D coordinates
        """
        d = self.calcDistanceToPlane( P )
        p = P - d*self.N
        return p
        
    def project2Plane2D( self, P ):
        """
        returns the closest points to P on the plane, in 2D in-plane
        coordinates
        """
        u = ((P - self.O) * self.X).sum(-1)
        v = ((P - self.O) * self.Y).sum(-1)
        return scipy.array([u,v]).T
    
    def plane2Dto3D( self, P ):
        """
        return 3D coordinates of from 2D in-plane coordinates
        """
        p = P[:,0,scipy.newaxis]*self.X + P[:,1,scipy.newaxis]*self.Y + self.O
        return p

    def angleToPlane(self, p):
        """
        calculates the angle between self normal and the normal
        or another plane p
        """
        return angle(self.N, p.N)

    def angleToVector(self, v):
        """ calcualte the angle between this plane and a vector v
        """
        # project vector onto plane
        v_proj = self.project2Plane3D(v)

        # calc angle between v and v_proj
        return angle(v_proj, v)

    def transformAffine(self, t):
        self.O = scipy.dot(
                    t,
                    scipy.hstack([
                        self.O,
                        1.0
                        ])
                    )[:3]
        self.N = scipy.dot(
                    t[:3,:3],
                    self.N
                    )
        self.X = scipy.dot(
                    t[:3,:3],
                    self.X
                    )
        self.Y = scipy.dot(
                    t[:3,:3],
                    self.Y
                    )

    def drawPlane(self, mscene, l=100, acolor=(1,0,0), ascale=10.0, scolor=(0,1,0), sopacity=0.5):
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
        l2 = 0.5*l
        planeCorners2D = np.array([[-l2,-l2],
                                   [ l2,-l2],
                                   [-l2, l2],
                                   [ l2, l2]])
        p3D = self.plane2Dto3D(planeCorners2D)
        arrow = mscene.mlab.quiver3d([self.O[0]],
                                    [self.O[1]],
                                    [self.O[2]],
                                    [self.N[0]],
                                    [self.N[1]],
                                    [self.N[2]],
                                    mode='arrow',
                                    scale_factor=ascale,
                                    color=acolor)

        square = mscene.mlab.mesh([[p3D[0,0], p3D[1,0]],
                                  [p3D[2,0], p3D[3,0]]],
                                 [[p3D[0,1], p3D[1,1]],
                                  [p3D[2,1], p3D[3,1]]],
                                 [[p3D[0,2], p3D[1,2]],
                                  [p3D[2,2], p3D[3,2]]],
                                  color=scolor,
                                  opacity=sopacity)

        return arrow, square
            
#===============================================================================#
def fitAxis3D( data, axis ):
    
    xtol = 1e-5
    ftol = 1e-5
    maxfev = 6*1000
    dataCoM = data.mean(0)
    
    def obj( x ):
         axis.setAB( x[0:3], x[3:6] )
         axisPoints = axis.findClosest(data)[0]
         #~ axisPoints = scipy.array([axis.findClosest(d)[0] for d in data])
         CoMDist = ((axis.b - dataCoM)**2.0).sum()
         d2 = ((data - axisPoints)**2.0).sum(1)
         return scipy.hstack( [d2, CoMDist] )
         
    xOpt = leastsq( obj, scipy.hstack([axis.a,axis.b]), xtol=xtol, ftol=ftol, maxfev=maxfev )[0]
    fittedRMSE = scipy.sqrt( obj(xOpt).mean() )
    axis.setAB( xOpt[0:3], xOpt[3:6] )
    
    return axis, xOpt, fittedRMSE

def fitPlaneLS( X ):
    
    # calc CoM
    CoM = X.mean(0)
    XC = X - CoM
    
    # eigen system
    A = scipy.zeros([3,3])
    A[0,0] = (XC[:,0]**2.0).sum()
    A[1,1] = (XC[:,1]**2.0).sum()
    A[2,2] = (XC[:,2]**2.0).sum()
    A[0,1] = A[1,0] = (XC[:,0]*XC[:,1]).sum()   
    A[0,2] = A[2,0] = (XC[:,0]*XC[:,2]).sum()   
    A[1,2] = A[2,1] = (XC[:,1]*XC[:,2]).sum()
    
    W, V = eig( A )
    V = V[:,W.argsort()]
    
    # project points onto plane
    P = Plane( CoM, V[:,0], V[:,1], V[:,2] ) 
    
    return P

def fitSphere( X ):

    def obj( x ):
        d = X - x[0:3]
        e = scipy.sqrt((d*d).sum(1)) - x[3]
        return e*e
    
    initCentre = scipy.mean(X, 0)
    initR = scipy.mean(mag(X - initCentre))
    x0 = scipy.hstack((initCentre, initR))
    try:
        xOpt = leastsq( obj, x0 )[0]
    except TypeError:
        print('WARNING: fitSphere: probably not enough points to fit')
        return 0, scipy.array([0,0,0,0])
    else:
        rmsOpt = scipy.sqrt((obj( xOpt )).mean())
        return rmsOpt, xOpt

def fitSphereAnalytic(X):
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
    A = scipy.zeros((3,3), dtype=float)
    A[0,0] =   (X[:,0]*(X[:,0]-X[:,0].mean())).mean()
    A[0,1] = 2*(X[:,0]*(X[:,1]-X[:,1].mean())).mean()
    A[0,2] = 2*(X[:,0]*(X[:,2]-X[:,2].mean())).mean()
    A[1,0] = 0
    A[1,1] =   (X[:,1]*(X[:,1]-X[:,1].mean())).mean()
    A[1,2] = 2*(X[:,1]*(X[:,2]-X[:,2].mean())).mean()
    A[2,0] = 0
    A[2,1] = 0
    A[2,2] =   (X[:,2]*(X[:,2]-X[:,2].mean())).mean()
    
    A = A + A.T

    B = scipy.zeros(3, dtype=float)
    B[0] = ((X[:,0]**2.0+X[:,1]**2+X[:,2]**2)*(X[:,0]-X[:,0].mean())).mean()
    B[1] = ((X[:,0]**2.0+X[:,1]**2+X[:,2]**2)*(X[:,1]-X[:,1].mean())).mean()
    B[2] = ((X[:,0]**2.0+X[:,1]**2+X[:,2]**2)*(X[:,2]-X[:,2].mean())).mean()

    # Center=(A\B).';
    Centre = scipy.dot(inv(A), B)
    
    Radius=scipy.sqrt((scipy.vstack([X[:,0]-Centre[0],
                                     X[:,1]-Centre[1],
                                     X[:,2]-Centre[2]])**2).sum(0).mean())

    return Centre, Radius


def fitBox( data, centre, axes ):
    
    maxIt = 10000
    # initialise axes
    X = Line3D( axes[0], centre )
    Y = Line3D( axes[1], centre )
    Z = Line3D( axes[2], centre )
    
    def obj( x ):
        
        # update box axes
        newB = x[:3]
        oldAs = scipy.array( [ v.a.copy() for v in [X,Y,Z] ] )
        newAs = transform3D.transformRigid3D( oldAs, [0,0,0,x[3],x[4],x[5]] )
        
        X.setAB( newAs[0], newB )
        Y.setAB( newAs[1], newB )
        Z.setAB( newAs[2], newB )
        
        # project data points
        pX = X.findClosest( data )[1]
        pY = Y.findClosest( data )[1]
        pZ = Z.findClosest( data )[1]
        
        # calc volume
        V = (pX.max() - pX.min())*(pY.max() - pY.min())*(pZ.max() - pZ.min())
        #~ print x
        print(V)
        return V
        
    x0 = scipy.hstack( [centre, [0,0,0]] )
    xOpt = fmin( obj, x0, maxiter=maxIt )

    finalCentre = xOpt[:3]
    finalAxes = scipy.array( [ v.a.copy() for v in [X,Y,Z] ] )
    finalVolume = obj( xOpt )
    
    # calculate fitted box dimensions
    pX = X.findClosest( data )[1]
    pY = Y.findClosest( data )[1]
    pZ = Z.findClosest( data )[1]
    finalDim = scipy.array( [pX.max() - pX.min(), pY.max() - pY.min(), pZ.max() - pZ.min()] )
    
    return finalCentre, finalVolume, finalDim,  [X,Y,Z]

def circumcentre3Points(a, b, c):
    """
    calculate the circum centre of 3 points and the circle radius,
    also calculates the normal of the circle plane.
    """
    a = scipy.array(a, dtype=float)
    b = scipy.array(b, dtype=float)
    c = scipy.array(c, dtype=float)
    
    # plane normal
    N = norm(scipy.cross((a-b), (a-c)))

    # midpoints of ab and ac
    d = 0.5*(a + b)
    e = 0.5*(a + c)

    # calculate perpendicular bisectors
    do = norm(np.cross(N, d-a))
    eo = norm(np.cross(N, e-a))

    ldo = Line3D(do, d)
    leo = Line3D(eo, e)

    # find intercept of bisectors
    dist, sldo, sleo = ldo.calcClosestDistanceToLine(leo)
    O = 0.5*(ldo.eval(sldo) + leo.eval(sleo))

    # calc radius
    R = scipy.mean([mag(a-O), mag(b-O), mag(c-O)])

    return O, R, N