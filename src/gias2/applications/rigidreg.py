"""
FILE: rigidreg.py
LAST MODIFIED: 24/05/17
DESCRIPTION:
Script for rigid-body registration of one model to another using ICP.

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

from os import path
import sys
import argparse
import numpy as np
from scipy.spatial import cKDTree
import copy
import logging

from gias2.registration import alignment_fitting as af
from gias2.registration import RBF
from gias2.mesh import vtktools

reg_methods = {
              'corr_r': AF.fitRigid,
              'corr_rs': AF.fitRigidScale,
              'corr_a': AF.fitAffine,
              'icp_r_st': AF.fitDataRigidEPDP,
              'icp_r_ts': AF.fitDataRigidDPEP,
              'icp_rs_st': AF.fitDataRigidScaleEPDP,
              'icp_rs_ts': AF.fitDataRigidScaleDPEP,
             }

def _makeX0(reg_method, source_pts, target_pts, t0=None, r0=None, s0=None):

        # auto initialise translation
        if t0 is None:
            t0 = target_pts.mean(0) - source_pts.mean(0)

        if r0 is None:
            r0 = np.array([0.0,0.0,0.0])
        else:
            r0 = np.deg2rad(r0)

        if s0 is None:
            s0 = 1.0

        if reg_method=='corr_a':
            return None
        elif 'rs' in reg_method:
            return np.hstack([t0, r0, s0])
        elif '_r' in reg_method:
            return np.hstack([t0, r0])
        else:
            return None

def register(reg_method, source, target, init_trans, init_rot, init_s,
    xtol=1e-3, samples=10000, pts_only=False, out=None, view=False):
    
    if pts_only:
        source_pts = source
        target_pts = target
    else:
        source_pts = source.v
        target_pts = target.v

    x0 = _makeX0(reg_method, source_pts, target_pts, init_trans, init_rot, init_s)
    print 'T0:', x0
    
    reg = reg_methods[reg_method]    
    if x0 is None:
        T, source_pts_reg, (rmse0, RMSE) = reg(
            source_pts, target_pts, xtol=xtol, sample=samples, outputErrors=True
            )
    else:
        T, source_pts_reg, (rmse0, RMSE) = reg(
            source_pts, target_pts, t0=x0, xtol=xtol, sample=samples,
            outputErrors=True
            )

    print 'Registered...'
    print 'RMSE:', RMSE
    print 'T:', T
    
    #=============================================================#
    # create regstered mesh
    if not pts_only:
        reg = copy.deepcopy(source)
        reg.v = source_pts_reg

    if out:
        if not pts_only:
            writer = vtktools.Writer(v=reg.v, f=reg.f)
            writer.write(args.out)
        else:
            n = np.arange(1,len(source_pts_reg))
            _out = np.hstack([n[:,np.newaxis], source_pts_reg])
            np.savetxt(
                args.out, _out, delimiter=',',
                fmt=['%6d', '%10.6f', '%10.6f', '%10.6f'],
                header='# rigid-body registered points'
                )

    #=============================================================#
    # view
    if view:
        try:
            from gias2.visualisation import fieldvi
            has_mayavi = True
        except ImportError:
            has_mayavi = False

        if has_mayavi:
            v = fieldvi.Fieldvi()
            if pts_only:
                v.addData('target points', target_pts, renderArgs={'mode':'point', 'color':(1,0,0)})
                v.addData('source points', source_pts, renderArgs={'mode':'point', 'color':(0,1,0)})
                v.addData('registered points', source_pts_reg, renderArgs={'mode':'point', 'color':(0.3,0.3,1)})
            else:
                v.addTri('target', target, renderArgs={'color':(1,0,0)})
                v.addTri('source', source, renderArgs={'color':(0,1,0)})
                v.addTri('registered', reg, renderArgs={'color':(0.3,0.3,1)})

            v.scene.background=(0,0,0)
            v.configure_traits()
        else:
            print('Visualisation error: cannot import mayavi')

    return source_pts_aligned, RMSE

def main(args):
    if args.points_only:
        source = np.loadtxt(args.source, skiprows=1, use_cols=(1,2,3))
    else:
        source = vtktools.loadpoly(args.source)
    
    if args.points_only:
        target = np.loadtxt(args.target, skiprows=1, use_cols=(1,2,3))
    else:
        target = vtktools.loadpoly(args.target)
    
    reg, rms = register(
        args.reg_method, source, target, args.init_trans, args.init_rot, args.init_scale,
        xtol=1e-3, samples=10000, pts_only=args.points_only, out=args.out, view=args.view,
        )

    logging.info('{}, rms: {}'.format(path.split(args.target)[1], rms))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Rigid-body registration of one model to another.')
    parser.add_argument(
        'reg-method',
        choices=[
            'corr_r', 'corr_rs', 'corr_a',
            'icp_r_st', 'icp_r_ts', 'icp_rs_st', 'icp_rs_ts',
            ],
        help='''Registration method. The choices are
corr_r: rigid on correspondent models
corr_rs: rigid plus scaling on correspondent models
corr_a: affine on correspondent models
icp_r_st: rigid using ICP, source to target distance minimisation 
icp_r_ts: rigid using ICP, target to source distance minimisation
icp_rs_st: rigid plus scaling using ICP, source to target distance minimisation 
icp_rs_ts: rigid plus scaling using ICP, target to source distance minimisation 
        '''
        )
    parser.add_argument(
        '-s', '--source',
        help='file path of the source model.'
        )
    parser.add_argument(
        '-t', '--target',
        help='file path of the target model.'
        )
    parser.add_argument(
        '-o', '--out',
        help='file path of the output registered model.'
        )
    parser.add_argument(
        '-p', '--points-only',
        help='Model are point clouds only. Expected file format is 1 header line, then n,x,y,z on each line after'
        )
    parser.add_argument(
        '--t0', nargs=3, type=float, default=None,
        help='initial translation, e.g. 100 0 0'
        )
    parser.add_argument(
        '--r0', nargs=3, type=float, default=None,
        help='initial rotation in degrees, e.g. 180 90 0'
        )
    parser.add_argument(
        '--s0', nargs=1, type=float, default=None,
        help='initial scaling, e.g. 1.1'
        )
    parser.add_argument(
        '-b', '--batch',
        help='file path of a list of model paths to fit. 1st model on list will be the source.'
        )
    parser.add_argument(
        '-d', '--outdir',
        help='directory path of the output registered models when using batch mode.'
        )
    parser.add_argument(
        '-v', '--view',
        action='store_true',
        help='Visualise measurements and model in 3D'
        )
    parser.add_argument(
        '-l', '--log',
        help='log file'
        )
    args = parser.parse_args()

    # start logging
    if args.log:
        log_fmt = '%(levelname)s - %(asctime)s: %(message)s'
        log_level = logging.INFO

        logging.basicConfig(
            filename=args.log,
            level=log_level,
            format=log_fmt,
            )
        logging.info(
            'Starting rigid-body registration',
            )

    if args.batch is None:
        main(args)
    else:
        model_paths = np.loadtxt(args.batch, dtype=str)
        args.source = model_paths[0]
        out_dir = args.outdir
        for i, mp in enumerate(model_paths[1:]):
            args.target = mp
            _p, _ext = path.splitext(path.split(mp)[1])
            args.out = path.join(out_dir, _p+'_rigidreg'+_ext)
            main(args)



                


