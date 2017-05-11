"""
FILE: rbfreg.py
LAST MODIFIED: 18-04-2017 
DESCRIPTION:
Module for geometric registration using rbfs

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import sys
import numpy as np
from scipy.spatial import cKDTree
import RBF

def generateBBoxPointsGrid(points, spacing=None, padding=None):
    """
    generate a grid of points internal to the bounding box of a set of points
    with spacing specified by tuple spacing.
    """
    points = np.array(points)
    # sample surface to get bounding box
    if padding is None:
        bboxMin = points.min(0)
        bboxMax = points.max(0)
    else:
        bboxMin = points.min(0) - padding
        bboxMax = points.max(0) + padding
    
    if spacing is None:
        spacing = (bboxMax - bboxMin)/3.0

    N = (bboxMax - bboxMin) / spacing
    for ni, n in enumerate(N):
        if n<2:
            N[ni] = 2
    
    # generate grid of points in bounding box
    PAll = np.array([[x,y,z] for z in np.linspace(bboxMin[2], bboxMax[2], N[2])\
                             for y in np.linspace(bboxMin[1], bboxMax[1], N[1])\
                             for x in np.linspace(bboxMin[0], bboxMax[0], N[0])
                             ])
        
    return PAll

def rbfreg(knots, source, target, basistype, basisargs, disttype):
    # create a RBF Coordinates field
    rcf = RBF.RBFComponentsField(3)

    # set radial basis function
    rcf.makeBasis(basistype, basisargs)

    # create a knot at every 10 data point
    rcf.setCentres(knots)

    # find correspondence
    if disttype=='st':
        # find closest target point to each source point
        # print('using ts')
        targetTree = cKDTree(target)
        closestInds = targetTree.query(source)[1]
        X = source
        Y = target[closestInds]
    elif disttype=='ts':
        # find closest source point to each target point
        # print('using ts')
        sourceTree = cKDTree(source)
        closestInds = sourceTree.query(target)[1]
        X = source[closestInds]
        Y = target
    else:
        raise ValueError("disttype must be 'st' or 'ts'")

    # fit knot weights
    rcf.fitData(X, Y)

    # evaluate registered source points
    sourceReg = rcf.evalMany(source).T

    # calculated RMS distance
    d = np.sqrt(((X - Y)**2.0).sum(1))
    dRms = np.sqrt((d*d).mean())
    print('RMS distance: {}'.format(dRms))

    return sourceReg, dRms, rcf, d


def _checkTermination(it, cost1, cost0, nknots, xtol, max_it, max_knots):

    if it>max_it:
        print('terminating because max iterations reached')
        return True

    if nknots>max_knots:
        print('terminating because max knots reached')
        return True

    if (abs(cost1-cost0)/cost0)<xtol:
        print(abs(cost1-cost0)/cost0)
        print('terminating because xtol reached')
        return True

    return False

def rbfRegIterative(source, target, distmode='ts', knots=None,
    basisType='gaussianNonUniformWidth', basisArgs=None, xtol=1e-3,
    minKnotDist=5.0, maxIt=50, maxKnots=500, maxKnotsPerIt=20):

    """Iterative RBF registration with greedy knots adding per iteration.

    inputs
    ------
    source: list of source point coordinates.
    target: list of target point coordinates.
    distmode: how source to target distance is calculated.
        'st': source to target - between each source point and closest target
              point.
        'ts': target to source - between each target point and closest source
              point.
    knots: list of knot coordinates.
    basisType: Radial basis type.
    basisArgs: dictionary of arguments fro the basis type.
    xtol: relative change in error for termination.
    maxIt: max number of iterations.
    minKnotDist: minimum distance between knots.
    maxKnotsPerIt: max number of knots to add per iteration.

    returns
    -------
    sourceCurrent: final morphed source point coordinates.
    rms: final RMS distance between morphed source and target points.
    rcf: final RBF deformation field.
    history: fitting results from each iteration. Dict containing
        'rms': rms error at each iteration,
        'ssdist': sum of squared distance at each iteration,
        'nknots': number knots at each iteration.
    """

    if basisArgs is None:
        basisArgs = {'s':10.0, 'scaling':50.0}

    if knots is None:
        knots = generateBBoxPointsGrid(source)

    terminate = False
    it = 0
    sourceCurrent = source
    ssdistCurrent = np.inf
    history = {
        'rms': [],
        'ssdist': [],
        'nknots': [],
        }  
    while not terminate:

        # perform fit
        sourceNew, rms, rcf, dist = rbfreg(
                                        knots,
                                        sourceCurrent,
                                        target,
                                        basisType,
                                        basisArgs,
                                        'ts', # 'ts' doesnt currently work
                                        )
        ssdistNew = (dist*dist).sum()

        # check if should terminate
        terminate = _checkTermination(
                        it, ssdistNew, ssdistCurrent, knots.shape[0], 
                        xtol, maxIt, maxKnots,
                        )

        # add knot
        if not terminate:
            print('\niteration {}'.format(it))
            # find source locations with highest errors
            sourceTree = cKDTree(sourceNew)
            tsDist, tsInds = sourceTree.query(target, k=1)
            sourceMaxDistInds = tsInds[np.argsort(tsDist)[::-1]]

            # go through source points from highest error and find
            # first one that is more than min_knot_dist from an
            # existing knot
            nKnotsAdded = 0
            for maxInd in sourceMaxDistInds:
                knotsTree  = cKDTree(knots)
                closestKnotDist = knotsTree.query(sourceNew[maxInd])[0]
                if closestKnotDist>minKnotDist:
                    knots = np.vstack([knots, sourceNew[maxInd]])
                    nKnotsAdded += 1
                
                if nKnotsAdded==maxKnotsPerIt:
                    break

            if nKnotsAdded==0:
                terminate = True
                print('terminating because no new knots can be added')

        sourceCurrent = sourceNew
        ssdistCurrent = ssdistNew
        history['rms'].append(rms)
        history['ssdist'].append(ssdistCurrent)
        history['nknots'].append(len(knots))
        it += 1

    return sourceCurrent, rms, rcf, history 
        
if __name__=='__main__':

    from gias2.registration import alignment_fitting as af
    from gias2.mesh import vtktools

    source_points_file = 'data/2007_5028_l.wrl'
    source_surf = vtktools.loadpoly(source_points_file)
    source_points = source_surf.v

    target_points_file = 'data/2008_0006_l.wrl'
    target_surf = vtktools.loadpoly(target_points_file)
    target_points = target_surf.v

    #=============================================================#
    # rigidly register source points to target points
    reg1_T, source_points_reg1, reg1_errors = af.fitDataRigidDPEP(
                                                source_points,
                                                target_points,
                                                xtol=1e-6,
                                                sample=1000,
                                                t0=np.deg2rad((0,0,0,0,0,0)),
                                                outputErrors=1
                                                )

    # add isotropic scaling to rigid registration
    reg2_T, source_points_reg2, reg2_errors = af.fitDataRigidScaleDPEP(
                                                source_points,
                                                target_points,
                                                xtol=1e-6,
                                                sample=1000,
                                                t0=np.hstack([reg1_T, 1.0]),
                                                outputErrors=1
                                                )

    #=============================================================#

    source_points_reg3, regRms, regRcf, regHist = rbfRegIterative(
        source_points_reg2, target_points, distmode='ts',
        basisType='gaussianNonUniformWidth',
        )

    knots = regRcf.C
    source_surf.v = source_points_reg3
    #=============================================================#
    # view
    try:
        from gias2.visualisation import fieldvi
        has_mayavi = True
    except ImportError:
        has_mayavi = False

    if has_mayavi:
        v = fieldvi.Fieldvi()
        v.addData('target points', target_points, renderArgs={'mode':'point', 'color':(1,0,0)})
        v.addTri('target', target_surf, renderArgs={'color':(1,0,0)})
        v.addTri('source morphed', source_surf, renderArgs={'color':(0,1,0)})
        v.addData('source points', source_points, renderArgs={'mode':'point'})
        v.addData('source points reg 1', source_points_reg1, renderArgs={'mode':'point'})
        v.addData('source points reg 2', source_points_reg2, renderArgs={'mode':'point'})
        v.addData('source points reg 3', source_points_reg3, renderArgs={'mode':'point', 'color':(0.5,0.5,1.0)})
        v.addData('knots', knots, renderArgs={'mode':'sphere', 'color':(0,1.0,0), 'scale_factor':5.0})
        v.scene.background=(0,0,0)
        v.configure_traits()
        


