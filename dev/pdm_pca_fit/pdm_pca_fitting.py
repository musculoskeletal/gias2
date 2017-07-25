"""
Dev script testing PCA fitting using a point distribution model, e.g. created
by rbf registration.
"""
from os import path
import sys
gias2_dir = path.split(__file__)[0]
gias2_dir = path.join(gias2_dir, '..', '..', 'src')
sys.path.append(gias2_dir)

import numpy as np
from scipy.spatial import cKDTree
from scipy.optimize import leastsq
from gias2.common import transform3D
from gias2.mesh import vtktools
from gias2.registration import alignment_fitting as af
from gias2.learning import PCA
from gias2.learning import PCA_fitting

def r2c13(x_recon):
	return x_recon.reshape((-1,3))

def mahalanobis( x ):
    return np.sqrt(np.multiply(x,x).sum()) 

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

def fitSSMTo3DPoints(data, ssm, fit_comps, fit_mode, fit_inds=None, mw=0.0,
    init_t=None, fit_scale=False, xtol=1e-6, sample=None,
    ldmk_targs=None, ldmk_evaluator=None, ldmk_weights=None,
    recon2coords=None, verbose=False
    ):
    
    """
    Fit a shape model to a set of non-correspondent points by optimising
    translation, rotation, and PCA scores.

    Rigid and optionally scaling registration should be performed before
    running this.
    
    inputs
    ------
    data: nx3 array of target point coordinates. N Points must correspond to
        points in the shape model
    ssm: a gias2.learning.PCA.PrincipalComponents object
    fit_comps: a list of PC modes to fit, e.g. [0,1,2] to fit the 1st 3 modes.
	fit_mode: {'sptp'|'tpsp'} source to target or target to source fitting
    fit_inds: [optional] restrict fitting to certain points in the 
        shape model
    mw: [float, optional] mahalanobis weight
    init_t: [list, optional] initial [tx,ty,tz,rx,ry,rz,[s]]
    fit_scale: [bool, False] fit for scaling
    
    Returns
    -------
    xOpt: array of optimised parameters. First 3 elements are translation, 
        then 3 for rotation, and the rest are PC scores in terms of standard
        deviations.
    recon_data_opt_t: fitted shape model points aligned to the mean shape
    data_t: target data points aligned to the mean shape
    recon_data_opt:  fitted shape model points aligned to the target data
    """

    if recon2coords is None:
        # Function to convert ssm data into point coordinates. Default is for
        # fieldwork models
        def recon2coords(xr):
            return xr.reshape((3,-1)).T

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

    # fit has to be done on a subset of recon data is projection is on a subset of variables
    # if fit_inds!=None:
    #     fit_inds = np.array(fit_inds, dtype=int)
    #     mean_data = recon2coords(ssm.getMean())[fit_inds,:]
    # else:
    #     mean_data = recon2coords(ssm.getMean())
    
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

    if fit_mode=='sptp':
        _dist = _dist_sptp
    elif fit_mode=='tpsp':
        _dist = _dist_tpsp
    else:
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
    x_opt = leastsq(_obj, x0, xtol=xtol)[0]
    print('')
        
    recon_data_opt, mdist_opt = _recon(x_opt)
    err_opt = _obj(x_opt)
    dist_opt_rms = np.sqrt((_dist(recon_data_opt, 0.0)**2.0).mean())

    if verbose:
        print('\nfinal rms distance: {}'.format(dist_opt_rms))

    return x_opt, recon_data_opt, (err_opt, dist_opt_rms, mdist_opt)

#=============================================================================#
targ_data_path = 'data/2008_0010_l.wrl'
mean_mesh_path = 'data/femur_mean.ply'
fitted_mesh_corr_path = 'data/pca_fitted_corr.ply'
fitted_mesh_nocorr_path = 'data/pca_fitted_nocorr.ply'
ssm_path = 'data/femur.pc.npz'

ssm = PCA.loadPrincipalComponents(ssm_path)
mean_mesh = vtktools.loadpoly(mean_mesh_path)
targ_data = vtktools.loadpoly(targ_data_path).v
fitted_mesh = vtktools.loadpoly(mean_mesh_path)

#=============================================================================#
# Correspondent fitting
#=============================================================================#
print('## CORR FITTING ##')
# xopt, recon_t, data_t, fitted_data = PCA_fitting.fitSSMTo3DPoints(
# 	targ_data, ssm, [0,1,2], mWeight=1.0, recon2coords=recon2coords, verbose=True
#     )

# fitted_mesh.v = fitted_data
# vtktools.savepoly(fitted_mesh, fitted_mesh_corr_path)

#=============================================================================#
# non correspondent
#=============================================================================#
print('## NON-CORR FITTING ##')

t0, src_data_aligned = af.fitDataRigidDPEP(
    mean_mesh.v, targ_data, xtol=1e-3, sample=5000
    )
xopt, fitted_pts, fit_errs = fitSSMTo3DPoints(
    targ_data, ssm, [0,1,2,3,4], 'sptp', init_t=t0, mw=1.0, xtol=1e-6, sample=10000,
    recon2coords=r2c13, verbose=True
    )
fitted_mesh.v = fitted_pts
vtktools.savepoly(fitted_mesh, fitted_mesh_nocorr_path)
