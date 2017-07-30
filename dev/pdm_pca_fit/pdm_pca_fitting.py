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
from gias2.mesh import vtktools
from gias2.registration import alignment_fitting as af
from gias2.registration import shapemodel as smreg
from gias2.learning import PCA

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
# print('## CORR FITTING ##')
# xopt, recon_t, data_t, fitted_data = PCA_fitting.fitSSMTo3DPoints(
#     targ_data, ssm, [0,1,2,3,4], mWeight=0.1, recon2coords=smreg.r2c13, verbose=True
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
xopt, fitted_pts, fit_errs = smreg.fitSSMTo3DPoints(
    targ_data, ssm, [0,1,2,3,4], 'tpsp', init_t=t0, mw=0.1, ftol=1e-6, sample=10000,
    recon2coords=smreg.r2c13, verbose=True
    )
fitted_mesh.v = fitted_pts
vtktools.savepoly(fitted_mesh, fitted_mesh_nocorr_path)
