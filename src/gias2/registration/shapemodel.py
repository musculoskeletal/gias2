"""
FILE: shapemodel.py
LAST MODIFIED: 31-7-2017 
DESCRIPTION: Shape model-based non-rigid registration.

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
from scipy.optimize import leastsq
from gias2.common import transform3D

#=============================================================================#
# utility functions
#=============================================================================#
def _sampleData(data, N):
    """
    Pick N evenly spaced points from data
    """

    if N<1:
        raise ValueError('N must be > 1')
    elif N>len(data):
        return data
    else:
        i = np.linspace(0,len(data)-1,N).astype(int)
        return data[i,:]

def mahalanobis( x ):
    return np.sqrt(np.multiply(x,x).sum()) 

#=============================================================================#
# reconstruction functions
#=============================================================================#
def r2c13(x_recon):
    """
    Reconstruction in to nx3 array of coordinates
    """
    return x_recon.reshape((-1,3))

def r2c31(x_recon):
    """
    Reconstruction in to 3xn array of coordinates
    """
    return x_recon.reshape((3,-1))

#=============================================================================#
# Main registration function
#=============================================================================#
def fitSSMTo3DPoints(data, ssm, fit_comps, fit_mode, fit_inds=None, mw=0.0,
    init_t=None, fit_scale=False, ftol=1e-6, sample=None,
    ldmk_targs=None, ldmk_evaluator=None, ldmk_weights=None,
    recon2coords=None, verbose=False
    ):
    
    """
    Fit a shape model to a set of non-correspondent points by optimising
    translation, rotation, and PCA scores.

    Rigid and optionally scaling registration should be performed before
    running this.
    
    arguments
    ---------
    data: nx3 array of target point coordinates.
    ssm: a gias2.learning.PCA.PrincipalComponents object
    fit_comps: a list of PC modes to fit, e.g. [0,1,2] to fit the 1st 3 modes.
    fit_mode: {'st'|'ts'|'corr'} source to target, target to source, or
        corresponding fitting. Use sptp if target datacloud covers more of the
        object than the shape model. Use tpsp if the target data cloud covers
        less of the object than the shape model. Use corr if the target
        datacloud is correspondent with the points in the shape model.
    
    keyword arguments
    -----------------
    fit_inds: [list of ints] restrict fitting to certain points in the 
        shape model
    mw: [float] weighting on the mahalanobis fitting penalty. Reasonable value
        is 0.1 to 1.0
    init_t: [list] initial [tx,ty,tz,rx,ry,rz,[s]]
    fit_scale: [bool] fit for scaling, default is False
    ftol: [float] relative error desired in sum of squared error
    sample: [int] number of points to sample from the target data
    ldmk_targ: [mx3 array] additional target landmark points to use during
        fitting
    ldmk_evaluator: functional that evaluates corresponding landmark
        coordinates from reconstructed data points. Should take as input a nx3
        array of reconstructed coordinates from the shape model and output a
        mx3 array of landmark coordinates.
    ldmk_weights: [mx1 float array] the fitting weight for each landmark.
    recon2coords: A function for reconstructing point coordinates from shape
        model data. e.g. r2c13 and r2c31 in this module.
    verbose: [bool] extra info during fit
    
    Returns
    -------
    xOpt: array of optimised parameters. First 3 elements are translation, 
        then 3 for rotation, and the rest are PC scores in terms of standard
        deviations.
    recon_data_opt: fitted shape model points
    (err_opt, dist_opt_rms, mdist_opt): final error, final rms distance, final
        mahalanobis distance
    """

    if recon2coords is None:
        # Function to convert ssm data into point coordinates. Default is for
        # nx3 point clouds.
        recon2coords = r2c13

    print('fitting ssm to points')
    if init_t is None:
        init_t = np.array([0.0,0.0,0.0, 0.0,0.0,0.0])
        if fit_scale:
            init_t = np.hstack([init_t, 1.0])
    else:
        init_t = np.array(init_t)

    if fit_scale:
        assert(len(init_t)==7)
    else:
        assert(len(init_t)==6)
    
    if fit_inds is None:
        print('fit_inds shape: None')
    else:
        print('fit_inds shape:', fit_inds.shape)

    if sample is not None:
        data = _sampleData(data, sample)

    #-------------------------------------------------------------------------#
    # define reconstruction functions
    #-------------------------------------------------------------------------#
    def _recon_no_scale(X):
        recon = ssm.reconstruct(
                    ssm.getWeightsBySD(fit_comps, X[6:]), fit_comps
                    )
        # reconstruct rigid transform
        recon_pts = transform3D.transformRigid3DAboutCoM(
                        recon2coords(recon), X[:6]
                        )
        mahalanobis_dist = mahalanobis(X[6:])
        return recon_pts, mahalanobis_dist

    def _recon_scale(X):
        recon = ssm.reconstruct(
                    ssm.getWeightsBySD(fit_comps, X[7:]), fit_comps
                    )
        # reconstruct rigid transform
        recon_pts = transform3D.transformRigidScale3DAboutCoM(
                        recon2coords(recon), X[:7]
                        )
        mahalanobis_dist = mahalanobis(X[7:])
        return recon_pts, mahalanobis_dist

    if fit_scale:
        _recon = _recon_scale
    else:
        _recon = _recon_no_scale

    #-------------------------------------------------------------------------#
    # define distance error functions
    #-------------------------------------------------------------------------#
    targ_tree = cKDTree(data)
    def _dist_sptp(recon_pts, m):
        return targ_tree.query(recon_pts)[0] + mw*m

    def _dist_tpsp(recon_pts, m):
        recon_tree = cKDTree(recon_pts)
        return recon_tree.query(data)[0] + mw*m

    def _dist_corr(recon_pts, m):
        return np.sqrt(((data - recon_pts)**2.0).sum(1))

    fit_modes_map = {
        'st': _dist_sptp,
        'ts': _dist_tpsp,
        'corr': _dist_corr,
    }

    try:
        _dist = fit_modes_map[fit_mode]
    except KeyError:
        raise ValueError('invalid fit mode {}'.format(fit_mode))

    #-------------------------------------------------------------------------#
    # define objective functions
    #-------------------------------------------------------------------------#
    def _obj_no_ldmks(X):
        # reconstruct data points
        recon_data, mdist = _recon(X)

        # calc error
        err = _dist(recon_data, mdist)

        if verbose:
            sys.stdout.write( '\robj rms:'+str(np.sqrt(err.mean())) )
            sys.stdout.flush()

        return err

    def _obj_ldmks(X):
        # reconstruct data points
        recon_data, mdist = _recon(X)
        recon_ldmks = ldmk_evaluator(recon_data.T.ravel())

        # calc error
        err_data = _dist(recon_data, mdist)
        
        err_ldmks = ((ldmk_targs - recon_ldmks)**2.0).sum(1)*ldmk_weights
        err = np.hstack([err_data, err_ldmks])

        if verbose:
            sys.stdout.write(
                '\rPC fit rmse: %6.3f (data: %6.3f) (landmarks: %6.3f)'%\
                (np.sqrt(err.mean()), np.sqrt(err_data.mean()), np.sqrt(err_ldmks.mean()))
                )
            sys.stdout.flush()

        return err
    
    if ldmk_targs is None:
        _obj = _obj_no_ldmks
    else:
        _obj = _obj_ldmks

    #-------------------------------------------------------------------------#
    # fit
    #-------------------------------------------------------------------------#
    x0 = np.hstack([init_t, np.zeros(len(fit_comps), dtype=float)])
    
    if verbose:
        recon_data_init, mdist_init = _recon(x0)
        err_init = _obj(x0)
        dist_init_rms = np.sqrt((_dist(recon_data_init, 0.0)**2.0).mean())
        print('\ninitial rms distance: {}'.format(dist_init_rms))

    x_opt = leastsq(_obj, x0, ftol=ftol)[0]
    print('')
        
    recon_data_opt, mdist_opt = _recon(x_opt)
    err_opt = _obj(x_opt)
    dist_opt_rms = np.sqrt((_dist(recon_data_opt, 0.0)**2.0).mean())

    if verbose:
        print('\nfinal rms distance: {}'.format(dist_opt_rms))

    return x_opt, recon_data_opt, (err_opt, dist_opt_rms, mdist_opt)