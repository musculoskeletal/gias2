"""
FILE: reg_test.py
LAST MODIFIED: 17-04-2017 
DESCRIPTION:
Radial Basis Function non-rigid registration of point clouds

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import numpy as np
import sys
import itertools
import copy
import time

from scipy.spatial import cKDTree

from gias2.common import transform3D
from gias2.registration import alignment_fitting as af
from gias2.mesh import vtktools

import RBF

try:
    from gias2.visualisation import fieldvi
    has_mayavi = True
except ImportError:
    has_mayavi = False

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
        
        target_tree = cKDTree(target)
        closest_inds = target_tree.query(source)[1]
        X = source
        Y = target[closest_inds]
    elif disttype=='ts':
        # find closest source point to each target point
        source_tree = cKDTree(source)
        closest_inds = source_tree.query(target)[1]
        X = source_points[closest_inds]
        Y = target
    else:
        raise ValueError("disttype must be 'st' or 'ts'")

    # fit knot weights
    rcf.fitData(X, Y)

    # evaluate registered source points
    source_reg = rcf.evalMany(source).T

    # calculated RMS distance
    d = np.sqrt(((X - Y)**2.0).sum(1))
    d_rms = np.sqrt((d*d).mean())
    print('RMS distance: {}'.format(d_rms))

    return source_reg, d_rms, rcf, d

def generateBBoxPointsGrid(points, spacing, padding=None):
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


#=============================================================================#
# Affine transformed target
# source points for registration
# source_points_file = 'data/BN00105_E15006_S01_tpm 7_001_deci.xyz'
# source_points = np.loadtxt(source_points_file)[::2]*1000.0 # scale to mm

# # target point cloud (for this example, affine transform of source points for fitting)
# target_points = transform3D.transformAffine(
#                     source_points,
#                     np.array([
#                         [1.2, 0.01, 0.02, 0.01],
#                         [0.01, 0.9, 0.01, -0.01],
#                         [0.02, 0.01, 1.1, 0.005],
#                         [0.0, 0.0, 0.0,  1],
#                         ])
#                     )

#=============================================================================#
# Completely different target
# source points for registration
source_points_file = 'data/2007_5028_l.wrl'
# source_points_file = 'data/2008_0051_pelvis.wrl'
source_surf = vtktools.loadpoly(source_points_file)
source_points = source_surf.v

# target points for registration
target_points_file = 'data/2008_0006_l.wrl'
# target_points_file = 'data/2008_0058_pelvis.wrl'
target_surf = vtktools.loadpoly(target_points_file)
target_points = target_surf.v

t0 = time.time()

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

t1 = time.time()

#=============================================================#
# do 1st fit
# knots1 = generateBBoxPointsGrid(source_points_reg2, [50,50,50], [10,10,10])
# basis_type = 'gaussianNonUniformWidth'
# basis_args =  {'s':10.0, 'scaling':50.0}
# source_points_reg3, rms1, rcf1, d1 = rbfreg(
#                                     # source_points_reg2[::100],
#                                     knots1, 
#                                     source_points_reg2,
#                                     target_points,
#                                     basis_type,
#                                     basis_args,
#                                     'st'
#                                     )

# # do 2nd fit
# knots2 = generateBBoxPointsGrid(source_points_reg3, [25,25,25], [10,10,10])
# basis_type = 'gaussianNonUniformWidth'
# basis_args =  {'s':10.0, 'scaling':25.0}
# source_points_reg4, rms2, rcf2, d2 = rbfreg(
#                                     # source_points_reg3[::10],
#                                     knots2,
#                                     source_points_reg3,
#                                     target_points,
#                                     basis_type,
#                                     basis_args,
#                                     'st'
#                                     )
#=============================================================#
# Iterative fit, adding knots in each iteration

def check_termination(it, cost1, cost0, nknots, xtol, max_it, max_knots):

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

basis_type = 'gaussianNonUniformWidth'
basis_args =  {'s':10.0, 'scaling':50.0}
xtol = 1e-3  # relative change in error for termination
min_knot_dist = 5.0  # minimum distance between knots
max_it = 20  # max number of iterations
max_knots = 1000  # max number of knots
max_knots_to_add_per_it = 20  # max number of knots to add per iteration

knots = generateBBoxPointsGrid(source_points_reg2, [50,50,50], [5,5,5])

terminate = False
it = 0
source_points_current = source_points_reg2
ssdist_current = np.inf
while not terminate:

    # perform fit
    source_points_new, rms_new, rcf, dist = rbfreg(
                                    knots,
                                    source_points_current,
                                    target_points,
                                    basis_type,
                                    basis_args,
                                    'st', # 'ts' doesnt currently work
                                    )
    ssdist_new = (dist*dist).sum()

    # check if should terminate
    terminate = check_termination(
                    it, ssdist_new, ssdist_current, knots.shape[0], 
                    xtol, max_it, max_knots,
                    )

    # add knot
    if not terminate:
        print('\niteration {}'.format(it))
        # find source locations with highest errors
        source_tree = cKDTree(source_points_new)
        ts_dist, ts_inds = source_tree.query(target_points, k=1)
        source_max_dist_inds = ts_inds[np.argsort(ts_dist)[::-1]]

        # go through source points from highest error and find
        # first one that is more than min_knot_dist from an
        # existing knot
        n_knots_added = 0
        for max_ind in source_max_dist_inds:
            knots_tree  = cKDTree(knots)
            closest_knot_dist = knots_tree.query(source_points_new[max_ind])[0]
            if closest_knot_dist>min_knot_dist:
                knots = np.vstack([knots, source_points_new[max_ind]])
                n_knots_added += 1
            
            if n_knots_added==max_knots_to_add_per_it:
                break

        if n_knots_added==0:
            terminate = True
            print('terminating because no new knots can be added')

    source_points_current = source_points_new
    ssdist_current = ssdist_new
    it += 1

source_points_reg3 = source_points_current
source_points_reg4 = source_points_new
knots1 = knots
knots2 = knots
source_surf.v = source_points_reg4

t2 = time.time()

#=============================================================#
# view
if has_mayavi:
    v = fieldvi.Fieldvi()
    v.addData('target points', target_points, renderArgs={'mode':'point', 'color':(1,0,0)})
    v.addTri('target', target_surf, renderArgs={'color':(1,0,0)})
    v.addTri('source morphed', source_surf, renderArgs={'color':(0,1,0)})
    v.addData('source points', source_points, renderArgs={'mode':'point'})
    v.addData('source points reg 1', source_points_reg1, renderArgs={'mode':'point'})
    v.addData('source points reg 2', source_points_reg2, renderArgs={'mode':'point'})
    v.addData('source points reg 3', source_points_reg3, renderArgs={'mode':'point', 'color':(0.5,0.5,1.0)})
    v.addData('source points reg 4', source_points_reg4, renderArgs={'mode':'point', 'color':(0,1.0,0)})
    v.addData('knots 1', knots1, renderArgs={'mode':'sphere', 'color':(0,1.0,0), 'scale_factor':5.0})
    v.addData('knots 2', knots2, renderArgs={'mode':'sphere', 'color':(0,1.0,0), 'scale_factor':5.0})
    v.scene.background=(0,0,0)
    v.configure_traits()
    v._drawTriSurface('target')
    v._drawTriSurface('source morphed')

#=============================================================#
print('Fitting done')
print('Final RMS error: {:6.2f}'.format(rms_new))
print('Total time: {:6.2f} s'.format(t2-t0))
print('Rigid registration time: {:6.2f} s'.format(t1-t0))
print('RBF registration time: {:6.2f} s'.format(t2-t1))
