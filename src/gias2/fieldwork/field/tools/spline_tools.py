"""
FILE: spline_tools.py
LAST MODIFIED: 24-12-2015
DESCRIPTION: spline fitting functions for generating meshes

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""
import logging
import shelve

import numpy
from scipy.interpolate import splprep, splev, splint
from scipy.linalg import svd
from scipy.optimize import fmin, leastsq
from scipy.spatial import cKDTree
from scipy.spatial.distance import euclidean

from gias2.common.geoprimitives import LineElement3D

# ======================================================================#
from gias2.common.math import norm

log = logging.getLogger(__name__)


class spline_parametric(object):
    """ Class for a spline curve fitted to a number of points in 3D for
    defining landmarks. Interpolates x,y,z given parametric variable u
    """

    def __init__(self, point_coords=None, tcku=None, smoothing=3.0, order=3, nknots=-1, per=1):
        """ point_coords is a list of N lists n long for n points in N-D 
        space.
        """
        if point_coords != None:
            self.tck, self.u = splprep(point_coords, s=smoothing, k=order, nest=nknots, per=per)
            self.order = order
        elif tcku != None:
            self.tck = tcku[0]
            self.u = tcku[1]
            self.order = self.tck[2]

        self.point_coords = point_coords
        self.smoothing = smoothing
        self.nknots = nknots
        self.type = 'parametric'

    # ==================================================================#
    def load(self, filename):
        try:
            s = shelve.open(filename, 'r')
        except:
            raise IOError(filename + ' not found')
        else:
            self.tck = s['tck']
            self.u = s['u']
            self.order = self.tck[2]
            s.close()

    def save(self, filename):
        """ save (tck, u) of spline S
        """
        s = shelve.open(filename + '.spl')
        s['tck'] = self.tck
        s['u'] = self.u
        s.close()

        return filename + '.spl'

    # ==================================================================#
    def eval(self, u, d=0):
        """ evaluate the spline g(u) at u. u is a 1D array of values.
        d is the derivative order
        """
        x = splev(u, self.tck, der=d)
        return numpy.array(x)

    # ==================================================================#
    def findClosest(self, p):
        """ returns coords and parameters u of point on the line closest
        to p
        """

        def obj(u):
            pLine = self.eval(u[0])
            d = (numpy.subtract(p, pLine) ** 2.0).sum()
            return d

        u0 = [0.5]
        uMin = fmin(obj, u0, disp=0)
        pClosest = self.eval(uMin)
        return pClosest, uMin

    # ==================================================================#
    def integrate(self, a, b):
        return splint(a, b, self.tck)

    # ==================================================================#
    # def draw( self, d, dp=0, figure=None, name='spline' ):
    #   x,y,z = self.eval(numpy.linspace(0.0,1.0,d))
    #   if dp:
    #       xp,yp,zp = self.eval(numpy.linspace(0.0,1.0,dp))

    #   if figure==None:
    #       figure = mlab.figure()

    #   mlab.plot3d( x,y,z, figure=figure, tube_radius=0.2, name=name )
    #   if dp:
    #       mlab.points3d( xp,yp,zp, figure=figure, scale_factor=1.0, name=name, color=(1.0,0.0,0.0) )

    #   return figure

    # ==================================================================#


def fitSplineCoeffDPEPObj(x, data, cShape, spline, nEP):
    spline.tck[1] = x.reshape(cShape)
    EP = spline.eval(numpy.linspace(0.0, 1.0, nEP)).T
    epTree = cKDTree(EP)
    d = epTree.query(list(data))[0]
    # ~ print 'rms', numpy.sqrt( (d*d).mean() )
    return d


def fitSplineCoeffEPDPObj(x, cShape, spline, nEP, dataTree):
    spline.tck[1] = x.reshape(cShape)
    EP = spline.eval(numpy.linspace(0.0, 1.0, nEP)).T
    d = dataTree.query(list(EP))[0]
    # ~ print 'rms', numpy.sqrt( (d*d).mean() )
    return d


def fitSplineCoeffPerDPEPObj(x, data, cShape, spline, nEP, k):
    temp = x.reshape(cShape)
    spline.tck[1] = numpy.hstack((temp, temp[:, :k]))
    EP = spline.eval(numpy.linspace(0.0, 1.0, nEP)).T
    epTree = cKDTree(EP)
    d = epTree.query(list(data))[0]
    # ~ print 'rms', numpy.sqrt( (d*d).mean() )
    return d


def fitSplineCoeffPerEPDPObj(x, cShape, spline, nEP, k, dataTree):
    temp = x.reshape(cShape)
    spline.tck[1] = numpy.hstack((temp, temp[:, :k]))
    EP = spline.eval(numpy.linspace(0.0, 1.0, nEP)).T
    d = dataTree.query(list(EP))[0]
    # ~ print 'rms', numpy.sqrt( (d*d).mean() )
    return d


def fitSplineCoeffPer(spline, data, xtol=1e-6):
    """ fit bspline coefficients to periodic data. The first k coeffs
    are the same as the last k, k being the spline order
    """

    k = spline.tck[2]
    # dont fit the last k coeffs
    x0 = numpy.array(spline.tck[1])[:, :-k].ravel()
    cShape = numpy.array(spline.tck[1])[:, :-k].shape

    def obj(x):
        temp = x.reshape(cShape)
        spline.tck[1] = numpy.hstack((temp, temp[:, :k]))
        P = numpy.array([spline.findClosest(d)[0] for d in data])
        err = ((data - P) ** 2.0).sum(1)
        # ~ print 'rms', numpy.sqrt( err.mean() )
        return err

    cOpt = leastsq(obj, x0, xtol=xtol)[0]
    temp = cOpt.reshape(cShape)
    spline.tck[1] = list(numpy.hstack((temp, temp[:, :k])))
    return spline


def fitSplineCoeffPerDPEP(spline, data, xtol=1e-6):
    """ fit bspline coefficients to periodic data. The first k coeffs
    are the same as the last k, k being the spline order
    """
    nEP = 100

    # dont fit the last k coeffs
    k = spline.tck[2]
    x0 = numpy.array(spline.tck[1])[:, :-k].ravel()
    cShape = numpy.array(spline.tck[1])[:, :-k].shape

    cOpt = leastsq(fitSplineCoeffPerDPEPObj, x0, args=(data, cShape, spline, nEP, k), xtol=xtol)[0]
    temp = cOpt.reshape(cShape)
    spline.tck[1] = list(numpy.hstack((temp, temp[:, :k])))
    return spline


def fitSplineCoeffPerEPDP(spline, data, xtol=1e-6):
    """ fit bspline coefficients to periodic data. The first k coeffs
    are the same as the last k, k being the spline order
    """
    nEP = 100
    dataTree = cKDTree(data)

    # dont fit the last k coeffs
    k = spline.tck[2]
    x0 = numpy.array(spline.tck[1])[:, :-k].ravel()
    cShape = numpy.array(spline.tck[1])[:, :-k].shape

    cOpt = leastsq(fitSplineCoeffPerEPDPObj, x0, args=(cShape, spline, nEP, k, dataTree), xtol=xtol)[0]
    temp = cOpt.reshape(cShape)
    spline.tck[1] = list(numpy.hstack((temp, temp[:, :k])))
    return spline


# ======================================================================#
def sortPointsCircular(points, startI=None):
    """ given a cloud of points in arbitrary order, orders the points by
    theta about the cloud CoM, from the startI'th point if given, or the
    1st point
    """

    # calc CoM
    CoM = points.mean(0)
    # find data plane axes from principal axes
    pAxes = calcPrincipalMomentsOfInertia(points, 1.0)[1]
    planeV = pAxes[:, [1, 2]].T
    if planeV[0][abs(planeV[0]).argmax()] < 0.0:
        planeV[0] *= -1.0

    planeV[1] = numpy.cross(planeV[0], pAxes[:, 0])

    # project points into plane
    planeX = projectDataOntoPlane(points, CoM, planeV[0], planeV[1]).T
    # find theta of points on the plane
    t = numpy.arctan2(planeX[1], planeX[0])
    t = numpy.where(t < 0.0, 2.0 * numpy.pi + t, t)
    # rank points by theta
    tArg = numpy.argsort(t)
    return points[tArg]


# ======================================================================#
def calcPrincipalMomentsOfInertia(d, mass):
    com = d.mean(0)
    P = d - com

    I11 = ((P[:, 1] * P[:, 1] + P[:, 2] * P[:, 2]) * mass).sum()
    I22 = ((P[:, 0] * P[:, 0] + P[:, 2] * P[:, 2]) * mass).sum()
    I33 = ((P[:, 1] * P[:, 1] + P[:, 0] * P[:, 0]) * mass).sum()
    I12 = -(P[:, 0] * P[:, 1] * mass).sum()
    I13 = -(P[:, 0] * P[:, 2] * mass).sum()
    I23 = -(P[:, 1] * P[:, 2] * mass).sum()

    I = numpy.array([[I11, I12, I13], [I12, I22, I23], [I13, I23, I33]])
    u, s, vh = svd(I)
    return (s.real, u.real)


# ======================================================================#
def projectDataOntoPlane(d, O, v0, v1):
    return numpy.array([(numpy.dot(di - O, v0), numpy.dot(di - O, v1)) for di in d])


# ======================================================================#

class pointCrawlerImage(object):
    """ traces shortest path in a binary image between 2 point
    """
    blocksize = 15
    blockShape = numpy.array([blocksize, blocksize, blocksize])
    tol = 4.0
    maxIt = 50
    smooth = 5.0

    def __init__(self, image):
        self.image = image

    def trace(self, p1, p2, debug=0):

        # initialise
        self.p2 = p2
        initialPaths = self._initialPathSearch(p1)
        line = None
        initialPath = 0
        goodPaths = []

        # get the path
        if debug:
            log.debug('tracing between', p1, ' and', self.p2)

        # starting-paths loop
        while (initialPath < len(initialPaths)):
            path = [numpy.array(p1), ]
            p = initialPaths[initialPath]
            d = euclidean(p, self.p2)
            it = 0

            # trace loop for each starting path
            while (d > self.tol) and (it < self.maxIt):
                path.append(numpy.array(p))
                pOld = p
                p = self._step(pOld)
                if p != None:
                    d = euclidean(p, self.p2)
                    it += 1
                    if debug:
                        log.debug(it, p, d)
                else:
                    # failed
                    break

            if d < self.tol:
                path.append(numpy.array(self.p2))
                goodPaths.append(numpy.array(path))

            initialPath += 1

        # find shortest path
        pathLengths = [self._calcPathLength(p) for p in goodPaths]
        if len(pathLengths) == 0:
            log.debug('ERROR: pointCrawler: trace unsuccessful, aborted')
            return None, None
        else:
            shortestPathI = numpy.argmin([self._calcPathLength(p) for p in goodPaths])
            bestPath = goodPaths[shortestPathI]
            if len(bestPath) > 3:
                line = spline_parametric(numpy.transpose(bestPath), smoothing=self.smooth, order=3)
                return line, bestPath
            elif len(bestPath) == 3:
                line = spline_parametric(numpy.transpose(bestPath), smoothing=self.smooth, order=2)
                return line, bestPath
            elif len(bestPath) == 2:
                line = LineElement3D(bestPath[0], path[1])
                return line, bestPath
            else:
                log.debug('ERROR: pointCrawler: length 1 path')
                return None, None

    def _initialPathSearch(self, centre, n=10):
        """ from the starting location, find the n best next points
        """
        # line from centre to target point
        l = LineElement3D(centre, self.p2)
        V, O = self._getSubVolume(centre)
        coords = numpy.array(V.nonzero()).transpose() + O
        if len(coords) == 0:
            log.debug('warning: no non-zero pixels')
            return None
        else:

            # get point inside a sphere of r = blocksize/2
            points = []
            distances = [euclidean(centre, p) for p in coords]
            for i in range(len(distances)):
                if distances[i] < self.blocksize / 2.0:
                    points.append(coords[i])

            # index sort distances to get n best points
            bestPoints = []
            bestPointIndices = numpy.argsort([l.findClosest(p)[1] for p in points])[::-1][:n]
            for i in bestPointIndices:
                bestPoints.append(points[i])

            return bestPoints

    def _step(self, centre):
        """ given position of a centre point, find the point in the 
        search block thats closest to self.p2
        """
        # line from centre to target point
        l = LineElement3D(centre, self.p2)
        V, O = self._getSubVolume(centre)
        coords = numpy.array(V.nonzero()).transpose() + O
        if len(coords) == 0:
            log.debug('warning: no non-zero pixels')
            return None
        else:

            # get point inside a sphere of r = blocksize/2
            points = []
            distances = [euclidean(centre, p) for p in coords]
            for i in range(len(distances)):
                if distances[i] < self.blocksize / 2.0:
                    points.append(coords[i])

            # find index of the point with the largest project on line
            i = numpy.argmax([l.findClosest(p)[1] for p in points])
            return points[i]

    def _getSubVolume(self, centre):
        corner = numpy.subtract(centre, self.blocksize / 2.0)
        corner = numpy.where(corner < 0.0, 0.0, corner)
        corner = corner.round().astype(int)
        return self.image[corner[0]: corner[0] + self.blocksize, \
               corner[1]: corner[1] + self.blocksize, \
               corner[2]: corner[2] + self.blocksize, ], corner

    def _calcPathLength(self, path):

        s = 0.0
        for i in range(len(path) - 1):
            s += euclidean(path[i], path[i + 1])

        return s


# ======================================================================#
class ridgeFitter(object):
    sliceShape = numpy.array((30, 30))

    def __init__(self, image, S, ):
        """ image: scan object with getSliceNormal method. S: spline
        object
        """
        self.image = image
        # ~ self.S = S

        self.splines = [S, ]

    def fit(self, it=3, knotsN=3, order=3, smoothing=3.0):
        """ samples original spline at knotsN positions evenly between 
        0.0 and 1,0 exclusive. Fit each of these positions to the apex
        of the closest ridge. iterative
        """

        for i in range(it):
            self.S = self.splines[-1]
            self.p0 = self.S.eval([0.0])
            self.p1 = self.S.eval([1.0])

            # get a list of coords and direction vectors for each new knot
            initial = self._sampleSpline(knotsN)

            newKnots = [self.p0, ]
            # fit each new knot to ridge apex
            for k in initial:
                newKnots.append(self._findApex(k[0], k[1]))

            newKnots.append(self.p1)

            # create new split interpolating new knots
            self.splines.append(spline_parametric(numpy.transpose(newKnots), smoothing=smoothing, order=order))
            # ~ self.sliceShape = (self.sliceShape*0.8).astype(int)

        return self.splines[-1]

    def _sampleSpline(self, n):
        """ returns the positions and direction vector of n points
        between 0.0 and 1.0 along the initial spline
        """
        T = numpy.linspace(0.0, 1.0, n + 2)[1:-1]
        ret = []
        for t in T:
            value = numpy.array(self.S.eval(t, d=0))
            dir = norm(self.S.eval(t, d=1))
            ret.append((value, dir))

        return ret

    def _findApex(self, initCentre, initNormal):

        slice = self.image.createSubSliceNormal(initCentre, initNormal, sliceShape=self.sliceShape)
        slice.fitQuadratic()
        sliceApex = slice.getQuadraticMidpoint()
        apex = slice.slice2StackCS(sliceApex)

        # ~ # get non zero pixels to fit a quadratic spline
        # ~ data = numpy.array( slice.I.nonzero() )
        # ~
        # ~ # sort by x (ind=1)
        # ~ sortInd = numpy.argsort( data[0] )
        # ~ newData = numpy.zeros( data.shape )
        # ~ for i in range( len(sortInd) ):
        # ~ newData[:,i] = data[:,sortInd[i]]
        # ~
        # ~ s = spline_parametric( newData, order=2, smoothing=1.0 )
        # ~
        # ~ v = s.eval( numpy.linspace(0.0,1.0,100), d=0 )
        # ~ d = s.eval( numpy.linspace(0.0,1.0,100), d=1 )

        # ~ print 'v: ',v
        # ~ print 'd: ',d
        # ~ plot.figure()
        # ~ plot.hold(True)
        # ~ plot.plot(v[0], v[1])
        # ~ plot.plot(d[0], d[1])

        # debug
        slice.viewSlice()
        return apex


# ======================================================================#

class pointCrawlerData(object):
    """ traces shortest path in a data cloud between 2 point
    """
    nNeighbours = 50
    rMaxNeighbours = 10.0
    leafSize = 50
    tol = 3.0
    maxIt = 100
    smooth = 5.0

    def __init__(self, data):
        self.data = data
        self.dataTree = cKDTree(self.data, self.leafSize)

    def trace(self, p1, p2, debug=0):

        # initialise
        self.p2 = p2
        initialPaths = self._initialPathSearch(p1)
        line = None
        initialPath = 0
        goodPaths = []
        goodPathClouds = []

        # get the path
        if debug:
            log.debug('tracing between', p1, ' and', self.p2)

        # starting-paths loop
        while (initialPath < len(initialPaths)):
            path = [numpy.array(p1), ]
            pathCloud = []
            p = initialPaths[initialPath]
            d = euclidean(p, self.p2)
            it = 0

            # trace loop for each starting path
            while (d > self.tol) and (it < self.maxIt):
                path.append(numpy.array(p))
                pOld = p
                p, neighb = self._step(pOld)
                pathCloud.append(neighb)
                if p != None:
                    d = euclidean(p, self.p2)
                    it += 1
                    if debug:
                        log.debug(it, p, d)
                else:
                    # failed
                    break

            if d < self.tol:
                path.append(numpy.array(self.p2))
                goodPaths.append(numpy.array(path))
                goodPathClouds.append(numpy.vstack(pathCloud))

            initialPath += 1

        # find shortest path
        pathLengths = [self._calcPathLength(p) for p in goodPaths]
        if len(pathLengths) == 0:
            raise RuntimeError('ERROR: pointCrawler: trace unsuccessful, aborted')
        else:
            shortestPathI = numpy.argmin([self._calcPathLength(p) for p in goodPaths])
            bestPath = numpy.array(goodPaths[shortestPathI])
            # ~ bestPath = numpy.array(goodPaths[0])

            # get unique points from the bestPatchCloud
            bestPathCloudUnique = []
            for p in goodPathClouds[shortestPathI]:
                # ~ pdb.set_trace()
                if tuple(p) not in bestPathCloudUnique:
                    bestPathCloudUnique.append(tuple(p))

            # if debug:
            #   from mayavi import mlab
            #   mlab.plot3d(bestPath[:,0], bestPath[:,1], bestPath[:,2])

            return numpy.array(bestPath), bestPathCloudUnique

    def _initialPathSearch(self, centre):
        """ from the starting location, find the n best next points
        """
        n = 10
        # get n closest neighbours
        nD, nI = self.dataTree.query(centre, n)
        nI = nI[numpy.where(nD < self.rMaxNeighbours)]
        try:
            nCoords = self.data[nI]
        except IndexError:
            log.debug('warning: no neighbours at starting point')
            return None
        else:
            # line from centre to target point
            l = LineElement3D(centre, self.p2)
            # sort neighbour points by largest projecting on l
            nRankedI = numpy.argsort([l.findClosest(p)[1] for p in nCoords])[::-1]
            return nCoords[nRankedI]
            # ~ return nCoords[[nRankedI[0]]]

    def _step(self, centre):
        """ given position of a centre point, find the point in the 
        search block that is closest to self.p2
        """
        # get n closest neighbours
        nD, nI = self.dataTree.query(centre, self.nNeighbours)
        nI = nI[numpy.where(nD < self.rMaxNeighbours)]
        try:
            nCoords = self.data[nI]
        except IndexError:
            log.debug('warning: no neighbours at starting point')
            return None
        else:
            # line from centre to target point
            l = LineElement3D(centre, self.p2)
            # find index of the point with the largest project on line
            i = numpy.argmax([l.findClosest(p)[1] for p in nCoords])
            return nCoords[i], nCoords

    def _calcPathLength(self, path):

        s = 0.0
        for i in range(len(path) - 1):
            s += euclidean(path[i], path[i + 1])

        return s


# ======================================================================#
def fitCurveEPDP(curve, data, pInit, debug=0):
    """ fit a 1D lagrange curve to data, fix ends. p0 should of the
    shape (dims,nNodes,paramsPerNode)
    """
    xtol = 1e-8
    ftol = 1e-8
    eps = 1e-6
    evalD = 60

    dim = pInit.shape[0]
    x0 = pInit[:, 1:-1].ravel()
    p = pInit.copy()
    nParams = pInit.shape[1]
    if debug:
        log.debug('dim:', dim)
        log.debug('p shape:', p.shape)
        log.debug('x0 shape:', x0.shape)
        log.debug('nParams:', nParams)

    def obj(x):
        # ~ pdb.set_trace()
        # set evaluate with new parameters
        p[:, 1:-1, :] = x.reshape((dim, nParams - 2, -1))
        ep = []
        for pI in p:
            ep.append(curve.evaluate_field_in_mesh(evalD, parameters=pI))
        ep = numpy.array(ep).T

        # for each ep, find dist to closest dp
        err = numpy.array([((data - epI) ** 2.0).sum(1).min() for epI in ep])
        if debug:
            log.debug('curve fitting rms:', numpy.sqrt(err.mean()))
        return err

    def objWithSmoothing(x):
        wD1 = 1.0e-5
        wD2 = 1.0e-5
        # ~ pdb.set_trace()
        # set evaluate with new parameters
        p[:, 1:-1, :] = x.reshape((dim, nParams - 2, -1))
        ep = []
        D1 = []
        D2 = []
        for pI in p:
            v, (d1, d2) = curve.evaluate_field_in_mesh(evalD, parameters=pI, derivs=-1)
            ep.append(v)
            D1.append(d1)
            D2.append(d2)

        ep = numpy.array(ep).T
        D1 = numpy.array(D1).T.squeeze()
        D2 = numpy.array(D2).T.squeeze()
        penalty = (D1 * D1).sum(1) * wD1 + (D2 * D2).sum(1) * wD2

        # for each ep, find dist to closest dp
        err = numpy.array([((data - epI) ** 2.0).sum(1).min() for epI in ep])
        err += penalty


        if debug:
            log.debug('curve fitting rms:', numpy.sqrt(err.mean()))
        return err

    # ~ xOpt = leastsq( obj, x0, ftol=ftol, xtol=xtol, epsfcn=eps )[0]
    # ~ xOpt = leastsq( obj, x0, ftol=ftol, xtol=xtol )[0]
    xOpt = leastsq(objWithSmoothing, x0, ftol=ftol, xtol=xtol)[0]
    p[:, 1:-1, :] = xOpt.reshape((dim, nParams - 2, -1))

    return p


# ======================================================================#
def fitCurveDPEP(curve, data, pInit, debug=0):
    """ fit a 1D lagrange curve to data, fix ends. p0 should of the
    shape (dims,nNodes,paramsPerNode)
    """
    xtol = 1e-3
    ftol = 1e-3
    eps = 1e-6
    evalD = 50

    dim = pInit.shape[0]
    x0 = pInit[:, 1:-1].ravel()
    p = pInit.copy()
    nParams = pInit.shape[1]
    if debug:
        log.debug('dim:', dim)
        log.debug('p shape:', p.shape)
        log.debug('x0 shape:', x0.shape)
        log.debug('nParams:', nParams)

    def obj(x):
        # ~ pdb.set_trace()
        # set evaluate with new parameters
        p[:, 1:-1, :] = x.reshape((dim, nParams - 2, -1))
        ep = []
        for pI in p:
            ep.append(curve.evaluate_field_in_mesh(evalD, parameters=pI))
        ep = numpy.array(ep).T.squeeze()
        epTree = cKDTree(ep)
        # for each ep, find dist to closest dp
        err, i = epTree.query(list(data))
        if debug:
            log.debug('curve fitting rms:', numpy.sqrt(err.mean()))
        return err

    def objWithSmoothing(x):
        wD1 = 1.0e-3
        wD2 = 1.0e-3
        # ~ pdb.set_trace()
        # set evaluate with new parameters
        p[:, 1:-1, :] = x.reshape((dim, nParams - 2, -1))
        ep = []
        D1 = []
        D2 = []
        for pI in p:
            v, (d1, d2) = curve.evaluate_field_in_mesh(evalD, parameters=pI, derivs=-1)
            ep.append(v)
            D1.append(d1)
            D2.append(d2)

        ep = numpy.array(ep).T.squeeze()
        D1 = numpy.array(D1).T.squeeze()
        D2 = numpy.array(D2).T.squeeze()
        penalty = (D1 * D1).sum(1) * wD1 + (D2 * D2).sum(1) * wD2

        # for each ep, find dist to closest dp
        # ~ pdb.set_trace()
        epTree = cKDTree(ep)
        err, i = epTree.query(list(data))
        err += penalty[i]

        if debug:
            log.debug('curve fitting rms:', numpy.sqrt(err.mean()))
        return err

    # ~ xOpt = leastsq( obj, x0, ftol=ftol, xtol=xtol, epsfcn=eps )[0]
    # ~ xOpt = leastsq( obj, x0, ftol=ftol, xtol=xtol )[0]
    xOpt = leastsq(objWithSmoothing, x0, ftol=ftol, xtol=xtol)[0]
    p[:, 1:-1, :] = xOpt.reshape((dim, nParams - 2, -1))

    return p
