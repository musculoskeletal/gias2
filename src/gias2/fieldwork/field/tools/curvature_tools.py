"""
FILE: curvature_tools.py
LAST MODIFIED: 24-12-2015
DESCRIPTION: functions and classes for evaluating curvature on 
fieldwork meshes.
    
===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""
import logging

from numpy import array, hstack, where, argwhere, sort, histogram, digitize, zeros, exp, unique
from scipy.spatial import cKDTree

log = logging.getLogger(__name__)


def filterCurv(C, cMin, cMax):
    """ filters data curvature so that values outside cMin and cMax
    are capped to the limits
    """
    return C.clip(cMin, cMax)


def normalise(x):
    # should scale -ves between -1 and 0, and +ves between 0 and 1
    xP = x[where(x > 0.0)]
    xN = x[where(x < 0.0)]
    # scale positives
    try:
        x = where(x > 0.0, (x - xP.min()) / (xP.max() - xP.min()), x)
    except:
        pass
    # scale negatives
    try:
        x = where(x < 0.0, (x - xN.max()) / (xN.max() - xN.min()), x)
    except:
        pass
    return x


def drawBins3d(D, C, bins=None, nBins=None):
    from mayavi import mlab
    if bins == None:
        bins = histogram(C, nBins)[1]
    else:
        nBins = len(bins) - 1

    CBini = digitize(C, bins)
    f = mlab.figure(bgcolor=(0.0, 0.0, 0.0), size=(800, 800))
    mlab.points3d(D[:, 0], D[:, 1], D[:, 2], CBini, mode='point', figure=f, vmin=1.0, vmax=nBins, scale_factor=0.2, )


def drawPolyMesh(v, f, C, bins=None, nBins=None):
    from mayavi import mlab
    if bins == None:
        bins = histogram(C, nBins)[1]
    else:
        nBins = len(bins) - 1

    CBini = digitize(C, bins)
    fig = mlab.figure(bgcolor=(0.0, 0.0, 0.0), size=(800, 800))
    mlab.triangular_mesh(v[:, 0], v[:, 1], v[:, 2], f, figure=fig, scalars=CBini, vmin=1.0, vmax=nBins)


def matchBins(x, yPD):
    """ find the bins of x (xBins) to match the prob density func of y
    (yPD )
    
    Instead of matching prob density, match n/lenX for each bin.
    """
    xSorted = sort(x)
    lenX = float(xSorted.shape[0])
    tol = 1e-3
    xBinsNew = []
    xI = 0
    n = 0
    lBin = xSorted[0]
    for ypd in yPD[:-1]:
        xBinsNew.append(xSorted[xI])
        n += 1
        xI += 1
        # ~ pd = n / (lenX * (xSorted[xI] - lBin) )
        pd = n / lenX
        while abs(pd - ypd) > tol:
            # ~ while pd < ypd:
            n += 1
            xI += 1
            # ~ pd = n / (lenX * (xSorted[xI] - lBin) )
            pd = n / lenX

        log.debug(pd, lBin, ypd, xI)
        n = 0
        lBin = xSorted[xI]

    xBinsNew.append(xSorted[xI])
    xBinsNew.append(xSorted[-1])
    return xBinsNew


def matchBins2(x, yCDF):
    """ find the bins of x (xBins) to match the normalised cumulative 
    distribution function of y (yCDF)
    """
    xSorted = sort(x)
    xCDF = (yCDF * len(x)).astype(int) - 1
    log.debug(xCDF)
    xBinsNew = xSorted[xCDF]

    return hstack((xSorted[0], xBinsNew))


def CDF(x, bins):
    """ calculates the cumulative distribution function of x given a
    sequence of bin edges
    """
    xLen = len(x)
    xSorted = sort(x)
    xCDF = []
    for i in range(1, len(bins) - 1):
        xCDF.append(where(xSorted < bins[i])[0].shape[0])
    xCDF.append(xLen)
    return array(xCDF), array(xCDF) / float(xLen)


def assignBins(aScalar, bScalar, nBins):
    # make bins for aScalar
    aBins = histogram(aScalar, nBins)[1]
    bBins = histogram(bScalar, nBins)[1]
    # digitise aScalar into aBins
    aBinIndices = digitize(aScalar, aBins)
    # digitise bScalar into bBins
    bBinIndices = digitize(bScalar, bBins)
    # handle out-of-bound bScalars
    # ~ bBinIndices = where( bBinIndices==0, 1, bBinIndices )
    # ~ bBinIndices = where( bBinIndices==nBins+1, nBins, bBinIndices )
    # get indices of a that belong to the bin of each bScalar
    bBinsA = [argwhere(aBinIndices == bBin).squeeze() for bBin in bBinIndices]
    for i, b in enumerate(bBinsA):
        if b.shape == ():
            bBinsA[i] = array([b])

    return bBinsA


def assignBins2(aScalar, aBins, bScalar, bBins):
    """ assign each point in b the corresponding bin index of a
    """

    # digitise aScalar into aBins
    aBinInd = digitize(aScalar, aBins)
    # digitise bScalar into bBins
    bBinInd = digitize(bScalar, bBins)
    # get indices of a that belong to the bin of each bScalar
    bBinsA = [argwhere(aBinInd == bBin).squeeze() for bBin in bBinInd]
    for i, b in enumerate(bBinsA):
        if b.shape == ():
            bBinsA[i] = array([b])

    return bBinsA


def smoothCurvField1(points, curvature):
    """ ckdtree
    """

    neighSize = 5
    dNeighMax = 10.0

    leafSize = 20
    sigma = 2.0  # for distance weight
    newCurv = zeros(curvature.shape)

    tree = cKDTree(points, leafSize)
    # for each point find the closest neighSize points
    DNeigh, PiNeigh = tree.query(list(points), neighSize)

    # for each points neighbours
    for i, dNeigh in enumerate(DNeigh):
        # get their curvature
        cNeigh = curvature[PiNeigh[i]]
        # calculate weight wrt to distance from point
        W = exp(- dNeigh / sigma)
        # normalise weights
        W = W / W.sum()
        # calculate smoothed curvature at point - weighted average
        newCurv[i] = (cNeigh * W).sum()
        log.debug(dNeigh)

    return newCurv


def smoothCurvField2(points, curvature):
    """ filter out duplicate points
    """
    neighSize = 20
    dNeighMax = 10.0
    sigma = 10.0  # for distance weight

    leafSize = 20
    newCurv = zeros(curvature.shape)

    tree = cKDTree(points, leafSize)
    # for each point find the closest neighSize points
    DNeigh, PiNeigh = tree.query(list(points), neighSize)

    # for each points neighbours
    for i, dNeigh in enumerate(DNeigh):
        dNeighU = unique(dNeigh)
        cNeigh = curvature[[PiNeigh[i][u] for u in [where(dNeigh == x)[0][0] for x in dNeighU]]]

        # get their curvature
        # ~ cNeigh = curvature[ PiNeigh[i] ]
        # calculate weight wrt to distance from point
        W = exp(- dNeighU / sigma)
        # normalise weights
        W = W / W.sum()
        # calculate smoothed curvature at point - weighted average
        newCurv[i] = (cNeigh * W).sum()

    return newCurv
