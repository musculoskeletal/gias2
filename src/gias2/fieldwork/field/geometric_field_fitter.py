"""
FILE: geometric_field_fitter.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: functions and classes for fitting geometric_fields.

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""
import logging
import pdb

from scipy import sparse
from scipy.optimize import leastsq
from scipy.spatial import cKDTree as KDTree

from gias2.common import math
from gias2.fieldwork.field import ensemble_field_function as EFF
from gias2.fieldwork.field import geometric_field
from gias2.fieldwork.field.tools import curvature_tools as CT
from gias2.fieldwork.field.topology import element_types

from numpy import array, newaxis, ones, sqrt, mean, dot, cos, sin, hstack, where, inf, digitize, linspace, zeros, cross, \
    dstack

log = logging.getLogger(__name__)


class geometryFit(object):
    """ object for fitting the geometry of a mesh
    """
    leaf_size = 20
    ftol = 1.0e-5
    xtol = 1.0e-5
    epsfcn = 1.0e-5
    cMin = -0.2
    cMax = 0.2
    projMethod = 'per_call'
    smoothing = False
    nBins = 5

    def __init__(self, G, data, eval_d, fitMode='geometry', projectionDirection='EPDP', projectionFreq='percall',
                 smoothing=False, sD=20, sNW=1000.0, sCW=1.0, dataCurvature=None):
        """ sN: smooth by minimising normal difference across element boundaries
        sNW: weight for sN
        sC: smooth by minimising curvature
        sCW: weight for sC
        sCMode: mean or gaussian curvature
        """

        self.T = array([[0.0, 0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0, 0.0],
                        [1.0, 1.0, 1.0, 1.0]])
        self.initH = None
        self.initK = None
        self.initK1 = None
        self.initK2 = None
        self.initMeshCurvature = None
        self.dataProjection = None

        self.G = G
        self.data = data
        self.dataCurvature = dataCurvature
        self.dataTree = KDTree(self.data, self.leaf_size)
        self.eval_d = eval_d
        self.fitMode = fitMode
        self.projectionDirection = projectionDirection
        self.projectionFreq = projectionFreq

        if smoothing:
            self.smoothing = True
            self.normSmoother = normalSmoother(self.G, sNW, sCW, sD)

        self.set_p0(self.G.field_parameters.copy())
        self._setObj()

    def setData(self, d):
        self.data = d
        self.dataTree = KDTree(self.data, self.leaf_size)

    def _setObj(self):

        if self.fitMode == 'geometry':
            self._obj = self._objGeometry
            if self.projectionDirection == 'EPDP':
                if self.projectionFreq == 'percall':
                    self.findClosestErr = self._findClosestErrEPDPPerCall
                elif self.projectionFreq == 'fixed':
                    self.findClosestErr = self._findClosestErrEPDPFixed
                    self.projectEPDP()
            elif self.projectionDirection == 'DPEP':
                if self.projectionFreq == 'percall':
                    self.findClosestErr = self._findClosestErrDPEPPerCall
                elif self.projectionFreq == 'fixed':
                    self.findClosestErr = self._findClosestErrDPEPFixed
                    self.projectDPEP()
            elif self.projectionDirection == 'both':
                if self.projectionFreq == 'percall':
                    self.findClosestErr = self._findClosestErrPerCallBoth
                elif self.projectionFreq == 'fixed':
                    self.findClosestErr = self._findClosestErrFixedBoth
                    self.projectDPEP()
        elif self.fitMode == 'curvature':
            self.dataCurvatureTree = KDTree(self.dataCurvature[:, newaxis], self.leaf_size)
            self.dataCurvatureTree = KDTree(self.dataCurvature[:, newaxis], self.leaf_size)
            self.normalise_curvature = 0
            self.cScale = 1.0

            # set _obj manually, so commented out below
            # ~ if self.projectionDirection=='EPDP':
            # ~ if self.projectionFreq=='percall':
            # ~ self._obj = self._objCurvatureDistancePerCall1EPDP
            # ~ elif self.projectionFreq=='fixed':
            # ~ self._obj= self._objCurvatureDistanceFixed1
            # ~ elif self.projectionDirection=='DPEP':
            # ~ if self.projectionFreq=='percall':
            # ~ self._obj = self._objCurvatureDistancePerCall1DPEP
            # ~ elif self.projectionFreq=='fixed':
            # ~ self._obj = self._objCurvatureDistanceFixed1
            # ~ elif self.projectionDirection=='both':
            # ~ if self.projectionFreq=='percall':
            # ~ self._obj = self._objCurvatureDistancePerCall1Both
            # ~ elif self.projectionFreq=='fixed':
            # ~ self._obj = self._objCurvatureDistanceFixed1

    def calc_initial_curvature(self):
        K, H, k1, k2 = self.G.evaluate_curvature_in_mesh(self.eval_d)
        self.initH = H
        self.initK = K
        self.initK1 = k1
        self.initK2 = k2
        self.initMeshCurvature = self.initH

    def preprocessCurvature(self, mode, cMin=None, cMax=None):
        if mode == 'filter':
            self.dataCurvature = CT.filterCurv(self.dataCurvature, cMin, cMax)
            self.initMeshCurvature = CT.filterCurv(self.initMeshCurvature, cMin, cMax)
        elif mode == 'normalise':
            self.dataCurvature = CT.normalise(self.dataCurvature)
            self.initMeshCurvature = CT.normalise(self.initMeshCurvature)
        elif mode == 'filternormalise':
            self.dataCurvature = CT.normalise(CT.filterCurv(self.dataCurvature, cMin, cMax))
            self.initMeshCurvature = CT.normalise(CT.filterCurv(self.initMeshCurvature, cMin, cMax))
        return

    def set_p0(self, p0):
        self.p0 = p0.squeeze().ravel()
        p = self.p0.reshape(3, -1)
        self.px0 = p[0]
        self.py0 = p[1]
        self.pz0 = p[2]
        self.nodes0 = p.T
        self.p0Affine = array([self.px0, self.py0, self.pz0, ones(self.px0.shape[0])])
        # ~ self._setObj()

    # def plotCurvatureHistograms( self, norm=False ):
    #   f1 = plot.figure()
    #   bData = plot.hist( self.dataCurvature, self.nBins, figure=f1, label='data', normed=norm )
    #   f2 = plot.figure()
    #   bMesh = plot.hist( self.initMeshCurvature, self.nBins, figure=f2, label='mesh', normed=norm )
    #   plot.show()
    #   return (bData[0], bData[1], bMesh[0], bMesh[1] )

    def reshapeParams(self, params):
        p = params.reshape((3, -1))
        p = array([p[0][:, newaxis], p[1][:, newaxis], p[2][:, newaxis]])
        return p

    def meshFit(self, it=0, errorOutput=False, verbose=True):
        self.it = 0
        if verbose:
            log.debug('mesh fit...')
        # ~ pdb.set_trace()
        maxFEval = it * self.p0.shape[0]

        initRMS = sqrt(self.objMesh(self.p0).mean())
        if verbose:
            log.debug('initial rms:', initRMS)

        output = leastsq(self.objMesh, self.p0, ftol=self.ftol, xtol=self.xtol, epsfcn=self.epsfcn, maxfev=maxFEval)

        finalRMS = sqrt(self.objMesh(output[0]).mean())
        if verbose:
            log.debug('final rms:', finalRMS)

        self.G.set_field_parameters(output[0].reshape((self.G.dimensions, -1, 1)))
        ep_coord = self.G.evaluate_geometric_field(self.eval_d)

        finalRMSDist = sqrt(mean(self.findClosestErr(ep_coord.T)))
        if verbose:
            log.debug('final rms distance:', finalRMSDist)

        if errorOutput:
            return self.reshapeParams(output[0]), output, (initRMS, finalRMS, finalRMSDist)
        else:
            return self.reshapeParams(output[0]), output

    def affineFit(self, it=0):
        log.debug(' affine fit...')
        self.it = 0
        T0 = array([[1.0, 0.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0, 0.0]]).ravel()

        maxFEval = it * 12
        output = leastsq(self.objAffine, T0, ftol=self.ftol, xtol=self.xtol, epsfcn=self.epsfcn, maxfev=maxFEval)
        TOpt = output[0].reshape((3, 4))
        self.T[:3, :] = TOpt
        # ~ print T
        pOpt = self.reshapeParams(dot(self.T, self.p0Affine)[:3, :].ravel())
        return (pOpt, output)

    def rigidScaleFit(self):
        self.it = 0
        log.debug('rigid and scaling fit...')
        # tx,ty,tz,rx,ry,rz,sx,sy,sz
        T0 = array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
        output = leastsq(self.objRigidScale, T0, ftol=self.ftol, xtol=self.xtol, epsfcn=self.epsfcn)

        TOpt = output[0]
        T = array([[1.0, 0.0, 0.0, TOpt[0]],
                   [0.0, 1.0, 0.0, TOpt[1]],
                   [0.0, 0.0, 1.0, TOpt[2]],
                   [1.0, 1.0, 1.0, 1.0]])

        S = array([[TOpt[-3], 0.0, 0.0],
                   [0.0, TOpt[-2], 0.0],
                   [0.0, 0.0, TOpt[-1]]])

        Rx = array([[1.0, 0.0, 0.0],
                    [0.0, cos(TOpt[3]), -sin(TOpt[3])],
                    [0.0, sin(TOpt[3]), cos(TOpt[3])]])

        Ry = array([[cos(TOpt[4]), 0.0, sin(TOpt[4])],
                    [0.0, 1.0, 0.0],
                    [-sin(TOpt[4]), 0.0, cos(TOpt[4])]])

        Rz = array([[cos(TOpt[5]), -sin(TOpt[5]), 0.0],
                    [sin(TOpt[5]), cos(TOpt[5]), 0.0],
                    [0.0, 0.0, 1.0]])

        T[:3, :3] = dot(dot(dot(Rx, Ry), Rz), S)
        pNew = self.reshapeParams(dot(T, self.p0Affine)[:3, :].ravel())
        return (pNew, output)

    def rigidFit(self, it=0):
        self.it = 0
        log.debug('rigid fit...')
        # tx,ty,tz,rx,ry,rz
        T0 = array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        maxFEval = it * 6
        output = leastsq(self.objRigid, T0, ftol=self.ftol, xtol=self.xtol, epsfcn=self.epsfcn, maxfev=maxFEval)

        TOpt = output[0]
        T = array([[1.0, 0.0, 0.0, TOpt[0]],
                   [0.0, 1.0, 0.0, TOpt[1]],
                   [0.0, 0.0, 1.0, TOpt[2]],
                   [1.0, 1.0, 1.0, 1.0]])

        Rx = array([[1.0, 0.0, 0.0],
                    [0.0, cos(TOpt[3]), -sin(TOpt[3])],
                    [0.0, sin(TOpt[3]), cos(TOpt[3])]])

        Ry = array([[cos(TOpt[4]), 0.0, sin(TOpt[4])],
                    [0.0, 1.0, 0.0],
                    [-sin(TOpt[4]), 0.0, cos(TOpt[4])]])

        Rz = array([[cos(TOpt[5]), -sin(TOpt[5]), 0.0],
                    [sin(TOpt[5]), cos(TOpt[5]), 0.0],
                    [0.0, 0.0, 1.0]])

        T[:3, :3] = dot(dot(Rx, Ry), Rz)
        pNew = self.reshapeParams(dot(T, self.p0Affine)[:3, :].ravel())
        return (pNew, output)

    def hostMeshFit(self, host, slaveObj, slaveXi=None, maxIt=0, smoothD=10, smoothW=1.0e-5):
        """ host mesh fit self.G using host (geometric_field) as the 
        host mesh and slaveObj as the objective function to minimise
        """
        log.debug('host mesh fit...')
        # calc slave node xi in host
        if slaveXi is None:
            log.debug('calculating slave xi...')
            slaveXi = array([host.findXi(0, node) for node in self.G.field_parameters[:, :, 0].T])
            # ~ savetxt( 'host_mesh_fitting/slaveXi.txt', slaveXi )

        # calc host basis values at slaveXis
        hostElem = host.ensemble_field_function.mesh.elements[0]
        evaluator = host.ensemble_field_function.evaluators[hostElem.type]
        basisFunction = host.ensemble_field_function.basis[hostElem.type]
        slaveBasis = basisFunction.eval(slaveXi.T)

        # initialise smoothing for host mesh
        hostParam0 = host.field_parameters.copy()
        smoother = elementDotPenalty3D(host.ensemble_field_function, smoothD, hostParam0, smoothW)
        smoothErr = array([])

        # hostmesh obj function
        def hostMeshObj(hostParams):
            hostParams = hostParams.reshape(3, -1, 1)
            host.set_field_parameters(hostParams)
            slaveParams = array([evaluator(slaveBasis, p) for p in hostParams]).ravel()

            slaveErr = slaveObj(slaveParams)
            # ~ smoothErr = smoother.obj( hostParams[:,:,newaxis] )
            Err = hstack((slaveErr, smoothErr))
            # ~ print 'it:', self.it, ' slave:', sqrt( slaveErr.mean() ), ' smoothing:', sqrt( smoothErr.mean() ), 'combined:', sqrt( Err.mean() )
            self.it += 1
            return slaveErr

        self.it = 0
        maxf = maxIt * (host.get_number_of_points() * 3)

        log.debug('initial rms:', sqrt(hostMeshObj(hostParam0).mean()))
        # do fit
        hostParamsOpt = \
            leastsq(hostMeshObj, hostParam0.ravel(), xtol=self.xtol, ftol=self.ftol, maxfev=maxf, epsfcn=self.epsfcn)[
                0].reshape(3, -1)[:, :, newaxis]
        host.set_field_parameters(hostParamsOpt)
        slaveParamsOpt = host.evaluate_geometric_field_at_element_points(0, slaveXi)[:, :, newaxis]
        self.G.set_field_parameters(slaveParamsOpt)

        finalHostRMS = sqrt(hostMeshObj(hostParamsOpt.ravel().copy()).mean())
        finalSlaveRMS = sqrt(slaveObj(slaveParamsOpt.ravel().copy()).mean())
        log.debug('final host rms:', finalHostRMS, ' final slave rms:', finalSlaveRMS)

        return hostParamsOpt, slaveParamsOpt, slaveXi, finalSlaveRMS

    def objRigidScale(self, t):

        T = array([[0.0, 0.0, 0.0, t[0]],
                   [0.0, 0.0, 0.0, t[1]],
                   [0.0, 0.0, 0.0, t[2]],
                   [1.0, 1.0, 1.0, 1.0]])

        S = array([[t[-3], 0.0, 0.0],
                   [0.0, t[-2], 0.0],
                   [0.0, 0.0, t[-1]]])

        Rx = array([[1.0, 0.0, 0.0],
                    [0.0, cos(t[3]), -sin(t[3])],
                    [0.0, sin(t[3]), cos(t[3])]])

        Ry = array([[cos(t[4]), 0.0, sin(t[4])],
                    [0.0, 1.0, 0.0],
                    [-sin(t[4]), 0.0, cos(t[4])]])

        Rz = array([[cos(t[5]), -sin(t[5]), 0.0],
                    [sin(t[5]), cos(t[5]), 0.0],
                    [0.0, 0.0, 1.0]])

        T[:3, :3] = dot(dot(dot(Rx, Ry), Rz), S)
        pNew = dot(T, self.p0Affine)[:3, :].ravel()

        d = self._obj(pNew)
        log.debug('it', self.it, '\trms:', sqrt((d).mean()))
        self.it += 1

        return d

    def objRigid(self, t):

        T = array([[1.0, 0.0, 0.0, t[0]],
                   [0.0, 1.0, 0.0, t[1]],
                   [0.0, 0.0, 1.0, t[2]],
                   [1.0, 1.0, 1.0, 1.0]])

        Rx = array([[1.0, 0.0, 0.0],
                    [0.0, cos(t[3]), -sin(t[3])],
                    [0.0, sin(t[3]), cos(t[3])]])

        Ry = array([[cos(t[4]), 0.0, sin(t[4])],
                    [0.0, 1.0, 0.0],
                    [-sin(t[4]), 0.0, cos(t[4])]])

        Rz = array([[cos(t[5]), -sin(t[5]), 0.0],
                    [sin(t[5]), cos(t[5]), 0.0],
                    [0.0, 0.0, 1.0]])

        T[:3, :3] = dot(dot(Rx, Ry), Rz)
        pNew = dot(T, self.p0Affine)[:3, :].ravel()
        d = self._obj(pNew)
        log.debug('it', self.it, '\trms:', sqrt((d).mean()))
        self.it += 1
        return d

    def objAffine(self, T):
        T = T.reshape((3, 4))
        self.T[:3, :] = T
        # ~ print T
        pNew = dot(self.T, self.p0Affine)[:3, :].ravel()
        d = self._obj(pNew)
        log.debug('it', self.it, '\trms:', sqrt((d).mean()))
        self.it += 1

        return d

    def objMesh(self, params):
        # ~ print 'params shape', params.shape
        # ~ dGeom = array([])
        dGeom = self._obj(params)

        if self.smoothing:
            dNorm, dCurv = self.normSmoother.smoothObj(self.reshapeParams(params))
            D = hstack((dGeom, dNorm, dCurv))

            log.debug('it', self.it,
                  ' geom:', sqrt((dGeom).mean()),
                  ' norm:', sqrt((dNorm).mean()),
                  ' curv:', sqrt((dCurv).mean()),
                  ' combined:', sqrt((D).mean()))
            self.it += 1
            return D
        else:
            log.debug('it', self.it, '\trms:', sqrt((dGeom).mean()))
            self.it += 1
            return dGeom

    def _objGeometry(self, params):
        # get new ep positions
        self.G.set_field_parameters(self.reshapeParams(params))
        ep_coord = self.G.evaluate_geometric_field(self.eval_d)
        # calculate error
        d = self.findClosestErr(ep_coord.T)
        return d

    # ~ def _objGeometrySobelov( self, params ):
    # ~ # get new ep positions
    # ~ p = self.reshapeParams(params)
    # ~ self.G.set_field_parameters( p )
    # ~ ep_coord = self.G.evaluate_geometric_field( self.eval_d, derivs=-1 )
    # ~ # calculate error
    # ~ epTree = KDTree( ep, self.leaf_size )
    # ~ projection = epTree.query( list(self.data), 1 )[1]
    # ~ d = self.data - ep[projection]

    def _objGeometrySobelov(self, params):
        # set new params
        p = self.reshapeParams(params)
        self.G.set_field_parameters(p)
        ep_coord = self.G.evaluate_geometric_field(self.eval_d)
        # calculate error
        dG = self.findClosestErr(ep_coord.T)
        # ~ pdb.set_trace()
        dS = self.sobelov2D.obj(p)

        # ~ print 'it', self.it, '\trms:', sqrt( (dG).mean() )
        self.it += 1
        try:
            return dG + dS
        except ValueError:
            return hstack((dG, dS))

    def _objGeometryAreaPenalty(self, params):
        # get new ep positions
        p = self.reshapeParams(params)
        self.G.set_field_parameters(p)
        ep_coord = self.G.evaluate_geometric_field(self.eval_d)
        # calculate error
        d = self.findClosestErr(ep_coord.T)
        dPenalty = self.areaPenalty(p)
        return hstack((d, dPenalty))

    def _objGeometryDotPenalty(self, params):
        # get new ep positions
        p = self.reshapeParams(params)
        self.G.set_field_parameters(p)
        ep_coord = self.G.evaluate_geometric_field(self.eval_d).T
        # calculate error
        d = self.findClosestErr(ep_coord)
        dPenalty = self.dotPenalty(p)
        log.debug('geom error:', sqrt(d.mean()), '\tpenalty:', sqrt(dPenalty.mean()))
        return hstack((d, dPenalty))

    def _objGeometryDotPenalty2(self, params):
        # get new ep positions
        p = self.reshapeParams(params)
        self.G.set_field_parameters(p)
        ep_coord = self.G.evaluate_geometric_field(self.eval_d).T
        # calculate error
        epTree = KDTree(ep_coord, self.leaf_size)
        d, d_i = epTree.query(list(self.data), 1)
        d = d * d

        dPenalty = self.dotPenalty(p)[d_i]
        log.debug('geom error:', sqrt(d.mean()), '\tpenalty:', sqrt(dPenalty.mean()))
        return d + dPenalty

    def _objGeometryDotAreaPenalty(self, params):
        # get new ep positions
        p = self.reshapeParams(params)
        self.G.set_field_parameters(p)
        ep_coord = self.G.evaluate_geometric_field(self.eval_d)
        # calculate error
        d = self.findClosestErr(ep_coord.T)
        dDotPenalty = self.dotPenalty(p)
        dAreaPenalty = self.areaPenalty(p)
        return hstack((d, dDotPenalty, dAreaPenalty))

    def _objCurvatureEPDP(self, params):
        """ finds the closest data point to each element point, and 
        calcs the sum of squared diff of datapoint curvature and
        element point initial curvature
        """
        # get new ep positions
        self.G.set_field_parameters(self.reshapeParams(params))
        ep_coord = self.G.evaluate_geometric_field(self.eval_d).T
        # find indices of closest data points
        # ~ leastArgs = self.findClosestPoint( ep_coord )[1]
        leastArgs = self.dataTree.query(list(ep_coord), 1, p=2)[1]
        # get corresponding data curvatures
        dataC = self.dataCurvature[leastArgs]
        # calculated sum of squared diff
        d = (self.initMeshCurvature - dataC) ** 2.0

        # add penalty for devation for p0
        # ~ dparams = self.p0Penalty( params ) * 0.01
        # ~ d = hstack( (d, dparams) )

        return d

    def _objCurvatureDPEP(self, params, plotProjection=0, fig=None):
        """ finds the closest element point to each data point, and 
        calcs the sum of squared diff of datapoint curvature and
        element point initial curvature
        """
        # get new ep positions
        self.G.set_field_parameters(self.reshapeParams(params))
        ep_coord = self.G.evaluate_geometric_field(self.eval_d).T
        # find indices of closest data points
        # ~ leastArgs = self.findClosestPoint( ep_coord )[1]
        epTree = KDTree(ep_coord)
        leastArgs = epTree.query(list(self.data), 1, p=2)[1]
        # get corresponding initial ep curvature
        epC = self.initMeshCurvature[leastArgs]
        # calculated sum of squared diff
        d = (self.dataCurvature - epC) ** 2.0

        # if plotProjection:
        #   from enthought.mayavi import mlab
        #   proj = ( ep_coord[leastArgs] - self.data )
        #   mlab.quiver3d( self.data[:,0], self.data[:,1], self.data[:,2], proj[:,0], proj[:,1], proj[:,2], line_width=0.0, figure=fig )

        # add penalty for devation for p0
        # ~ dparams = self.p0Penalty( params ) * 0.01
        # ~ d = hstack( (d, dparams) )

        return d

    def _objCurvatureDistancePerCall1Both(self, params):
        dEPDP = self._objCurvatureDistancePerCall1EPDP(params)
        dDPEP = self._objCurvatureDistancePerCall1DPEP(params)
        return hstack((dEPDP, dDPEP))

    def _objCurvatureDistancePerCall1EPDP(self, params):
        """ update projection every call
        """
        neighSize = 20
        neighRadius = 10.0

        # calculate new ep positions
        self.G.set_field_parameters(self.reshapeParams(params))
        ep_coord = self.G.evaluate_geometric_field(self.eval_d).T

        # get data neighbourhood for each ep
        neighbourhoods = self.dataTree.query(list(ep_coord), neighSize, distance_upper_bound=neighRadius)[1]
        # remove missing neighbours
        neighbourhoods = [where(n == self.dataTree.n, n[0], n) for n in neighbourhoods]

        # get curvatures for the points in each neighbourhood
        try:
            nCurv = array([self.dataCurvature[n] for n in neighbourhoods])
        except IndexError:
            # some eps have no dp in neighbourhood radius
            nCurv = []
            for n_i, n in enumerate(neighbourhoods):
                if self.dataTree.n in n:
                    log.debug('WARNING: no data points found in neighbourhood for element point', n_i)
                    nCurv.append(ones(neighSize) * inf)
                else:
                    nCurv.append(self.dataCurvature[n])

            nCurv = array(nCurv)

        # ~ pdb.set_trace()
        # the index in each neighbourhood for the point with most similar curvature
        temp = nCurv - self.initH[:, newaxis]
        nClosest_i = (temp * temp).argmin(1)
        # get data index of the closest dp in each neighbourhood
        d_i = []
        for i, n in enumerate(neighbourhoods):
            d_i.append(n[nClosest_i[i]])

        d = ep_coord - self.data[d_i]
        d = (d * d).sum(1)
        # weight by curvature
        w = self.initMeshCurvature * self.dataCurvature[d_i]
        d = d * w

        # ~ dPenalty = self.areaPenalty( self.reshapeParams(params)  )
        # from enthought.mayavi import mlab
        # proj = ( self.data[d_i] - ep_coord ) * w[:,newaxis]
        # mlab.quiver3d( ep_coord[:,0], ep_coord[:,1], ep_coord[:,2], proj[:,0], proj[:,1], proj[:,2], line_width=0.0 )
        # pdb.set_trace()

        return d
        # ~ return hstack((d, dPenalty))

    def _objCurvatureDistancePerCall1DPEP(self, params, plotProjection=0, fig=None):
        """ update projection every call
        """
        neighSize = 20
        neighRadius = 20.0

        # calculate new ep positions
        self.G.set_field_parameters(self.reshapeParams(params))
        ep_coord = self.G.evaluate_geometric_field(self.eval_d).T

        epTree = KDTree(ep_coord, self.leaf_size)
        # get ep neighbourhood for each dp
        neighbourhoods = epTree.query(list(self.data), neighSize, distance_upper_bound=neighRadius)[1]
        # remove missing neighbours
        neighbourhoods = [where(n == epTree.n, n[0], n) for n in neighbourhoods]

        # get curvatures for the points in each neighbourhood
        try:
            nCurv = array([self.initMeshCurvature[n] for n in neighbourhoods])
        except IndexError:
            # some eps have no dp in neighbourhood radius
            nCurv = []
            for n_i, n in enumerate(neighbourhoods):
                if self.epTree.n in n:
                    log.debug('WARNING: no element points found in neighbourhood for data point', n_i)
                    nCurv.append(ones(neighSize) * inf)
                else:
                    nCurv.append(self.initMeshCurvature[n])

            nCurv = array(nCurv)

        # ~ pdb.set_trace()
        # the index in each neighbourhood for the point with most similar curvature
        temp = nCurv - self.dataCurvature[:, newaxis]
        nClosest_i = (temp * temp).argmin(1)
        # get ep index of the closest dp in each neighbourhood
        ep_i = []
        for i, n in enumerate(neighbourhoods):
            ep_i.append(n[nClosest_i[i]])

        # calculate distance
        d = self.data - ep_coord[ep_i]
        d = (d * d).sum(1) * self.dataCurvature * self.initMeshCurvature[ep_i]

        # weight by curvature

        # ~ dPenalty = self.areaPenalty( self.reshapeParams(params)  )
        # if plotProjection:
        #   from enthought.mayavi import mlab
        #   #~ proj = ( ep_coord[ep_i] - self.data ) * ( self.dataCurvature * self.initMeshCurvature[ep_i] )[:,newaxis]
        #   proj = ( ep_coord[ep_i] - self.data )
        #   mlab.quiver3d( self.data[:,0], self.data[:,1], self.data[:,2], proj[:,0], proj[:,1], proj[:,2], line_width=0.0, figure=fig )
        # ~ pdb.set_trace()

        return d
        # ~ return hstack((d, dPenalty))

    def _objCurvatureDistancePerCall2EPDP(self, params):

        neighSize = 100
        neighRadius = 0.02

        # calculate new ep positions
        self.G.set_field_parameters(self.reshapeParams(params))
        ep_coord = self.G.evaluate_geometric_field(self.eval_d).T

        try:
            similarCurvDPCoord = self.similarCurvDPCoord
        except AttributeError:
            # 1. get indices of dp with similar curvature
            neighbourhoods = \
                self.dataCurvatureTree.query(list(self.initH[:, newaxis]), neighSize, distance_upper_bound=neighRadius)[
                    1]
            # remove missing neighbours
            neighbourhoods = [where(n == self.dataCurvatureTree.n, n[0], n) for n in neighbourhoods]
            # get coordinates of dps in the curv neighb of each ep

            # get curvatures for the points in each neighbourhood
            try:
                self.similarCurvDPCoord = array([self.data[n] for n in neighbourhoods])
            except IndexError:
                # some eps have no dp in neighbourhood radius
                self.similarCurvDPCoord = []
                for n_i, n in enumerate(neighbourhoods):
                    if self.dataCurvatureTree.n in n:
                        log.debug('WARNING: no data points found in neighbourhood for element point', n_i)
                        self.similarCurvDPCoord.append(ones((neighSize, 3)) * inf)
                    else:
                        self.similarCurvDPCoord.append(self.data[n])

                self.similarCurvDPCoord = array(similarCurvDPCoord)

            similarCurvEPCoord = self.similarCurvDPCoord

        # for each ep, calc the distance to the closest point in the neighbourhood
        d = similarCurvDPCoord - ep_coord[:, newaxis]
        d = (d * d).sum(2).min(1)
        # ~ pdb.set_trace()

        # ~ # get index in neighbourhood
        # ~ d_i_neigh = ( ( self.similarCurvDPCoord - ep_coord[:,newaxis] )**2.0 ).sum(2).argmin(1)
        # ~ # get datapoints
        # ~ D = []
        # ~ for i, d_i_n in enumerate( d_i_neigh ):
        # ~ D.append( self.data[ neighbourhoods[i][d_i_n] ] )

        # ~ from enthought.mayavi import mlab
        # ~ for i, p in enumerate(D):
        # ~ if not mod(i,10):
        # ~ mlab.plot3d( [p[0], ep_coord[i,0]], [p[1], ep_coord[i,1]], [p[2], ep_coord[i,2]], tube_radius=None )
        # ~
        # ~ pdb.set_trace()

        return d

    def _objCurvatureDistancePerCall2DPEP(self, params, plotProjection=0, fig=None):

        neighSize = 100
        neighRadius = 0.05

        # calculate new ep positions
        self.G.set_field_parameters(self.reshapeParams(params))
        ep_coord = self.G.evaluate_geometric_field(self.eval_d).T
        # ~ epTree = KDTree( ep_coord )

        try:
            neighbourhoods = self.neighbourhoods
        except AttributeError:
            epCurvatureTree = KDTree(self.initMeshCurvature[:, newaxis], self.leaf_size)
            # 1. get indices of dp with similar curvature
            self.neighbourhoods = \
                epCurvatureTree.query(list(self.dataCurvature[:, newaxis]), neighSize,
                                      distance_upper_bound=neighRadius)[1]
            # remove missing neighbours
            self.neighbourhoods = [where(n == epCurvatureTree.n, n[0], n) for n in self.neighbourhoods]
            neighbourhoods = self.neighbourhoods

        # get coordinates for the eps in each dp curvature neighbourhood
        try:
            similarCurvEPCoord = array([ep_coord[n] for n in neighbourhoods])
        except IndexError:
            # some eps have no dp in neighbourhood radius
            similarCurvEPCoord = []
            for n_i, n in enumerate(neighbourhoods):
                if self.initMeshCurvature.shape[0] in n:
                    log.debug('WARNING: no element points found in neighbourhood for data point', n_i)
                    similarCurvEPCoord.append(ones((neighSize, 3)) * inf)
                else:
                    similarCurvEPCoord.append(ep_coord[n])

            similarCurvEPCoord = array(similarCurvEPCoord)

        # for each ep, calc the distance to the closest point in the neighbourhood
        dtemp = similarCurvEPCoord - self.data[:, newaxis]
        d = (dtemp * dtemp).sum(2).min(1)

        # if plotProjection:
        #   from enthought.mayavi import mlab
        #   e_i_neigh = ( ( similarCurvEPCoord - self.data[:,newaxis] )**2.0 ).sum(2).argmin(1)
        #   E = []
        #   for i, e_i_n in enumerate( e_i_neigh ):
        #       E.append( ep_coord[ neighbourhoods[i][e_i_n] ] )  

        #   proj = E - self.data
        #   #~ proj = ( ep_coord[ep_i] - self.data )
        #   mlab.quiver3d( self.data[:,0], self.data[:,1], self.data[:,2], proj[:,0], proj[:,1], proj[:,2], line_width=0.0, figure=fig )
        #   #~ pdb.set_trace()

        return d

    def assignBinsDPEP(self, dataBins, meshBins):

        # Get ep indices for each ep bin
        meshBinI = digitize(self.initMeshCurvature, meshBins)
        self.epIBins = [where(meshBinI == i)[0] for i in range(1, len(meshBins))]
        # get datacoords for each bin
        dataBinI = digitize(self.dataCurvature, dataBins)
        self.dataCoordBins = [self.data[where(dataBinI == i)[0]] for i in range(1, len(dataBins))]

    def _objCurvatureDistancePerCall2DPEP2(self, params, plotProjection=0, fig=None):
        # calculate new ep positions
        self.G.set_field_parameters(self.reshapeParams(params))
        ep_coord = self.G.evaluate_geometric_field(self.eval_d).T

        # build kdtrees for each bin of element points
        epTrees = [KDTree(ep_coord[i], self.leaf_size) for i in self.epIBins]

        D = []
        epBinInd = []
        # for each bin of datapoints, find their closest ep in epTrees
        for i, dBin in enumerate(self.dataCoordBins):
            d, ep_i = epTrees[i].query(list(dBin), 1)
            D.append(d)
            epBinInd.append(ep_i)

        D = hstack(D)

        # if plotProjection:
        #   from enthought.mayavi import mlab
        #   for i, d in enumerate( self.dataCoordBins ):
        #       proj = ep_coord[self.epIBins[i][epBinInd[i]]] - d
        #       mlab.quiver3d( d[:,0], d[:,1], d[:,2], proj[:,0], proj[:,1], proj[:,2], line_width=0.0, figure=fig )

        # ~ for di, dp in enumerate( self.dataCoordBins[3] ):
        # ~ if not mod(di,10):
        # ~ try:
        # ~ ep = ep_coord[self.epIBins[3][epBinInd[3][di]] ]
        # ~ mlab.plot3d( [dp[0],ep[0]], [dp[1],ep[1]], [dp[2],ep[2]], line_width=0.0, figure=fig)
        # ~ except:
        # ~ pdb.set_trace()
        # ~ try:
        # ~ # get coords of eps for each datapoint bin
        # ~ DPEPBinsCoords = [ ep_coord[EPs] for EPs in self.DPEPBins ]
        # ~ except AttributeError:
        # ~ self.DPEPBins = CT.assignBins2( self.initMeshCurvature, self.dataCurvature, self.nBins )
        # ~ DPEPBinsCoords = [ ep_coord[EPs] for EPs in self.DPEPBins ]
        # ~
        # ~ pdb.set_trace()
        # ~ d = array( [ min( ( (DPEPBinsCoords[i] - self.data[i])**2.0).sum(1) )  for i in xrange(self.data.shape[0]) ] )
        # ~ d = []
        # ~ for i in xrange(self.data.shape[0]):
        # ~ print i
        # ~ d.append( min( ( (DPEPBinsCoords[i] - self.data[i])**2.0).sum(1) ) )
        # ~
        # ~ if plotProjection:
        # ~ from enthought.mayavi import mlab
        # ~ d_i_bins = [ argmin( ( (DPEPBinsCoords[i] - self.data[i])**2.0).sum(1) )  for i in xrange(self.data.shape[0]) ]
        # ~ proj = array( [ DPEPBinsCoords[i][d_i_bins[i]] - self.data[i] for i in xrange(self.data.shape[0]) ] )
        # ~ proj = ( ep_coord[ep_i] - self.data )
        # ~ mlab.quiver3d( self.data[:,0], self.data[:,1], self.data[:,2], proj[:,0], proj[:,1], proj[:,2], line_width=0.0, figure=fig )

        return D * D

    def _objCurvatureDistanceFixed1(self, params):
        """ use initial projection
        """
        neighSize = 20
        neighRadius = 50.0

        # calculate new ep positions
        self.G.set_field_parameters(self.reshapeParams(params))
        ep_coord = self.G.evaluate_geometric_field(self.eval_d)

        try:
            d = ep_coord.T - self.data[self._d_i]
            d = (d * d).sum(1)
        except AttributeError:
            # get data neighbourhood for each ep
            neighbourhoods = self.dataTree.query(list(ep_coord.T), neighSize, distance_upper_bound=neighRadius)[1]
            # remove missing neighbours
            neighbourhoods = [where(n == self.dataTree.n, n[0], n) for n in neighbourhoods]

            # get curvatures for the points in each neighbourhood
            nCurv = array([self.dataCurvature[n] for n in neighbourhoods])
            # ~ pdb.set_trace()
            # the index in each neighbourhood for the point with most similar curvature
            temp = nCurv - self.initH[:, newaxis]
            nClosest_i = (temp * temp).argmin(1)
            # get data index of the closest dp in each neighbourhood
            self._d_i = []
            for i, n in enumerate(neighbourhoods):
                self._d_i.append(n[nClosest_i[i]])

            d = ep_coord.T - self.data[self._d_i]
            d = (d * d).sum(1)

        return d

    def _objCurvatureDistanceFixed2(self, params):

        neighSize = 100
        neighRadius = 0.01

        # calculate new ep positions
        self.G.set_field_parameters(self.reshapeParams(params))
        ep_coord = self.G.evaluate_geometric_field(self.eval_d).T

        try:
            # for each ep, calc the distance to the closest point in the neighbourhood
            d = self.similarCurvDPCoord - ep_coord[:, newaxis]
            d = (d * d).sum(2).min(1)
            # ~ pdb.set_trace()
        except AttributeError:
            # 1. get indices of dp with similar curvature
            neighbourhoods = \
                self.dataCurvatureTree.query(list(self.initH[:, newaxis]), neighSize, distance_upper_bound=neighRadius)[
                    1]
            # remove missing neighbours
            neighbourhoods = [where(n == self.dataCurvatureTree.n, n[0], n) for n in neighbourhoods]
            # get coordinates of dps in the curv neighb of each ep
            self.similarCurvDPCoord = array([self.data[n] for n in neighbourhoods])
            # for each ep, calc the distance to the closest point in the neighbourhood
            d = self.similarCurvDPCoord - ep_coord[:, newaxis]
            d = (d * d).sum(2).min(1)

        return d

    def _findClosestErrPerCallBoth(self, ep):
        dEPDP = self._findClosestErrEPDPPerCall(ep)
        dDPEP = self._findClosestErrDPEPPerCall(ep)
        return hstack((dEPDP, dDPEP))

    def _findClosestErrFixedBoth(self, ep):
        dEPDP = self._findClosestErrEPDPFixed(ep)
        dDPEP = self._findClosestErrDPEPFixed(ep)
        return hstack((dEPDP, dDPEP))

    def _findClosestErrEPDPFixed(self, ep):
        d = ep - self.EPDPProjection
        return (d * d).sum(1)

    def _findClosestErrDPEPFixed(self, ep):
        d = self.data - ep[self.DPEPProjectionI]
        return (d * d).sum(1)

    def _findClosestErrEPDPPerCall(self, ep):
        # ~ projection = self.dataTree.query( list(ep), 1, p=2 )[1]
        # ~ d = ep - self.data[projection]
        # ~ return (d*d).sum(1)

        d = self.dataTree.query(list(ep), 1, p=2)[0]
        return d * d

    def _findClosestErrDPEPPerCall(self, ep):
        epTree = KDTree(ep, self.leaf_size)
        # ~ projection = epTree.query( list(self.data), 1 )[1]
        # ~ d = self.data - ep[projection]
        d = epTree.query(list(self.data), 1)[0]

        # ~ from enthought.mayavi import mlab
        # ~ proj = ( ep[projection] - self.data )
        # ~ mlab.quiver3d( self.data[:,0], self.data[:,1], self.data[:,2], proj[:,0], proj[:,1], proj[:,2], line_width=0.0 )
        # ~ pdb.set_trace()

        # ~ return ( d*d ).sum(1)
        return d * d

    def projectDPEP(self):
        """ for each data point, find the index of the closest element
        point
        """
        ep = self.G.evaluate_geometric_field(self.eval_d).T
        # ~ epTree = KDTree( ep, self.leaf_size )
        # ~ self.DPEPProjectionI = epTree.query( list(self.data), 1 )[1]

        # ~ self.DPEPProjectionI = cdist( self.data, ep, 'euclidean' ).argmin(1)

        I = []
        for d in self.data:
            I.append(((ep - d) ** 2.0).sum(1).argmin())

        self.DPEPProjectionI = array(I)

        return

    def projectEPDP(self):
        """ for each element point, find the closest data point
        """
        ep = self.G.evaluate_geometric_field(self.eval_d).T
        leastArgs = self.dataTree.query(list(ep), 1, p=2)[1]
        self.EPDPProjection = self.data[leastArgs]
        return

    def p0Penalty(self, params):
        p = params.reshape(3, -1).T
        return ((self.nodes0 - p) ** 2.0).sum(1)


# ======================================================================#
def makeObjEPEP(G, data, evalD, dataWeights=None, evaluator=None, nClosestPoints=None, treeArgs=None, epIndex=None,
                epXi=None, matPoints=None):
    if evaluator is None:
        evaluator = geometric_field.makeGeometricFieldEvaluatorSparse(G, evalD, epIndex=epIndex, epXi=epXi,
                                                                      matPoints=matPoints)

    data = array(data)

    if dataWeights is None:
        def obj(p):
            ep = evaluator(p).T
            err = ((ep - data) ** 2.0).sum(1)
            # err = ne.evaluate( 'sum((ep - data)**2.0, axis=1)', local_dict={'data':data, 'ep':ep} )
            return err
    else:
        data = array(data.T)

        def obj(p):
            ep = evaluator(p)
            err = (dataWeights * ((ep - data) ** 2.0)).sum(0)
            return err

    return obj


def makeObjDPEP(G, data, evalD, dataWeights=None, evaluator=None, nClosestPoints=1, treeArgs={}, epIndex=None):
    if evaluator is None:
        evaluator = geometric_field.makeGeometricFieldEvaluatorSparse(G, evalD)

    data1 = list(data)

    if epIndex is None:
        if nClosestPoints > 1:
            if dataWeights is None:
                def obj(p):
                    ep = evaluator(p).T
                    err = mean(KDTree(ep).query(data1, k=nClosestPoints, **treeArgs)[0], 1)
                    return err * err
            else:
                def obj(p):
                    ep = evaluator(p).T
                    err = mean(KDTree(ep).query(data1, k=nClosestPoints, **treeArgs)[0], 1)
                    return err * err * dataWeights
        else:
            if dataWeights is None:
                def obj(p):
                    ep = evaluator(p).T
                    err = KDTree(ep).query(data1, k=1, **treeArgs)[0]
                    return err * err
            else:
                def obj(p):
                    ep = evaluator(p).T
                    err = KDTree(ep).query(data1, k=1, **treeArgs)[0]
                    return err * err * dataWeights
    else:
        if nClosestPoints > 1:
            if dataWeights is None:
                def obj(p):
                    ep = evaluator(p).T[epIndex]
                    err = mean(KDTree(ep).query(data1, k=nClosestPoints, **treeArgs)[0], 1)
                    return err * err
            else:
                def obj(p):
                    ep = evaluator(p).T[epIndex]
                    err = mean(KDTree(ep).query(data1, k=nClosestPoints, **treeArgs)[0], 1)
                    return err * err * dataWeights
        else:
            if dataWeights is None:
                def obj(p):
                    ep = evaluator(p).T[epIndex]
                    err = KDTree(ep).query(data1, k=1, **treeArgs)[0]
                    return err * err
            else:
                def obj(p):
                    ep = evaluator(p).T[epIndex]
                    err = KDTree(ep).query(data1, k=1, **treeArgs)[0]
                    return err * err * dataWeights

    return obj


def makeObjEPDP(G, data, evalD, dataWeights=None, evaluator=None, nClosestPoints=1, treeArgs={}, epIndex=None):
    if evaluator is None:
        evaluator = geometric_field.makeGeometricFieldEvaluatorSparse(G, evalD)

    dataTree = KDTree(data)

    if epIndex is None:
        if nClosestPoints > 1:
            if dataWeights is None:
                def obj(p):
                    ep = evaluator(p).T
                    err = mean(dataTree.query(list(ep), k=nClosestPoints, **treeArgs)[0], 1)
                    return err * err
            else:
                def obj(p):
                    ep = evaluator(p).T
                    d, i = dataTree.query(list(ep), k=nClosestPoints, **treeArgs)
                    d = mean(d, 1)
                    w = dataWeights[i]
                    return d * d * w
        else:
            if dataWeights is None:
                def obj(p):
                    ep = evaluator(p).T
                    err = dataTree.query(list(ep), k=1, **treeArgs)[0]
                    return err * err
            else:
                def obj(p):
                    ep = evaluator(p).T
                    d, i = dataTree.query(list(ep), k=1, **treeArgs)
                    w = dataWeights[i]
                    return d * d * w
    else:
        if nClosestPoints > 1:
            if dataWeights is None:
                def obj(p):
                    ep = evaluator(p).T[epIndex]
                    err = mean(dataTree.query(list(ep), k=nClosestPoints, **treeArgs)[0], 1)
                    return err * err
            else:
                def obj(p):
                    ep = evaluator(p).T[epIndex]
                    d, i = dataTree.query(list(ep), k=nClosestPoints, **treeArgs)
                    d = mean(d, 1)
                    w = dataWeights[i]
                    return d * d * w
        else:
            if dataWeights is None:
                def obj(p):
                    ep = evaluator(p).T[epIndex]
                    err = dataTree.query(list(ep), k=1, **treeArgs)[0]
                    return err * err
            else:
                def obj(p):
                    ep = evaluator(p).T[epIndex]
                    d, i = dataTree.query(list(ep), k=1, **treeArgs)
                    w = dataWeights[i]
                    return d * d * w

    return obj


def makeObj2Way(G, data, evalD, dataWeights=None, evaluator=None, nClosestPoints=1, treeArgs={}):
    objEPDP = makeObjEPDP(G, data, evalD, dataWeights, evaluator, nClosestPoints, treeArgs)
    objDPEP = makeObjDPEP(G, data, evalD, dataWeights, evaluator, nClosestPoints, treeArgs)

    def obj(x):
        err = hstack([objEPDP(x), objDPEP(x)])
        return err

    return obj


# ======================================================================#
class normalSmoother2(object):
    """ penalises surface normal mismatch at 2D element boundaries
    """

    def __init__(self, F):

        self.F = F
        self.en2el, self.el2en = self.F.get_mapping()  # ens2elem, elem2ens maps
        # self._procMap()
        self._findCommonEdges()

    def _procMap(self):
        """ find element points of shared nodes
        """
        self.edgePoints = []
        self.vertexPoints = []  # not implemented

        # for each ensemble point
        for p in list(self.en2el.keys()):
            # edge nodes
            if len(self.en2el[p]) == 2:
                E = list(self.en2el[p].items())  # element points mapped to by ensemble point p
                self.edgePoints.append(((E[0][0], list(E[0][1].keys())[0]), (
                    E[1][0], list(E[1][1].keys())[0])))  # ((element1, elementnode1), (element2, elementnode2))
            # vertex nodes
            elif len(self.en2el[p]) > 2:
                pass

        return

    def _findCommonEdges(self):
        """
        Create a list of tuples containing pairs of element edge instances that are overlapped
        """

        def _findCommonElem(c):
            """Given a list contain lists of elemnum and elem node number connected to 
            edge nodes, find the element that is connect to all edge points

            e.g. c = [[(1,0),(2,1)], [(1,1)], [(1,2),(3,2)]] will return 1 since nodes
            from elem 1 is connect in all cases.
            """

            # create sets of elements connects to each point
            elemSets = []
            for ci in c:
                elemSets.append(set([t[0] for t in ci]))

            if len(elemSets) == 0:
                return None, None

            # find the intersection of the sets
            sCommon = elemSets[0]
            for s in elemSets[1:]:
                sCommon = sCommon.intersection(s)

            if len(sCommon) == 1:
                # the only elem left is the one that is connected
                connElemNum = sCommon.pop()

                # get the connected element's edge points
                connPoints = []
                for ci in c:
                    connPoints += [t[1] for t in ci if t[0] == connElemNum]

                return connElemNum, connPoints

            if len(sCommon) == 0:
                return None, None

            if len(sCommon) > 1:
                raise RuntimeError('Edge shared between more than 2 elements.')

        # list of shared edges
        # each tuple contains (elem1, edge1, elem2, edge2, [1|-1]) where the
        # last element is 1 if direction is align or -1 if direction is
        # opposite 
        self.commonEdges = []
        doneEdges = set()

        # loop through each element
        for elemNum, elem in self.F.mesh.elements.items():
            # loop through each set of edge points
            for edgeInd, edgePoints in enumerate(elem.edge_points):
                if (elemNum, edgeInd) not in doneEdges:
                    edge = elem.edges[edgeInd]
                    # get elem number and elem node number connected to each edge point
                    connections = [self.F.mesh.connectivity[(elemNum, ep)] for ep in edgePoints]
                    # find the element and its edge points that are connected to this edge
                    connElemNum, connPoints = _findCommonElem(connections)

                    if connElemNum is not None:
                        # find the edge that overlaps with this edge
                        connElem = self.F.mesh.elements[connElemNum]
                        connEdge, direction = connElem.get_edge_by_edge_points(connPoints)

                        self.commonEdges.append(
                            (elemNum, edge, connElemNum, connEdge, direction)
                        )
                        doneEdges.add((elemNum, edgeInd))
                        doneEdges.add((connElemNum, connEdge.edge_number))

    def _procEdge(self, D):
        """ for each pair in edgePoints, find the element and edge shared
        and generate eval points along the shared edges
        """
        self.edgeEvalBasis = []
        self.edgeEvalPoints = []
        self.nPairs = 0  # number of pairs of edge points
        doneEdges = set()  # list of element edges that have been done
        # for each pair
        for enum1, edge1, enum2, edge2, direction in self.commonEdges:
            # get elements
            e1 = self.F.mesh.elements[enum1]
            e2 = self.F.mesh.elements[enum2]

            # check if element edges are reversed
            if direction < 0:
                eval1 = edge1.get_elem_coord(linspace(0.0, 1.0, D))
                eval2 = edge2.get_elem_coord(linspace(1.0, 0.0, D))
            else:
                eval1 = edge1.get_elem_coord(linspace(0.0, 1.0, D))
                eval2 = edge2.get_elem_coord(linspace(0.0, 1.0, D))

            # get basis values for these
            basis1 = [self.F.basis[e1.type].eval_derivatives(eval1.T, d).T for d in ((1, 0), (0, 1))]
            basis2 = [self.F.basis[e2.type].eval_derivatives(eval2.T, d).T for d in ((1, 0), (0, 1))]

            self.edgeEvalPoints.append((enum1, eval1, enum2, eval2))  # ( element1, ep1, element2, ep2 )
            self.edgeEvalBasis.append((enum1, basis1, enum2, basis2))  # ( element1, basis1, element2, basis2 )
            self.nPairs += eval1.shape[0]

    def _procEdgeOld(self, D):
        """ for each pair in edgePoints, find the element and edge shared
        and generate eval points along the shared edges
        """
        self.edgeEvalBasis = []
        self.edgeEvalPoints = []
        self.nPairs = 0  # number of pairs of edge points
        doneEdges = set()  # list of element edges that have been done
        # for each pair
        for ep1, ep2 in self.edgePoints:
            # get elements
            e1 = self.F.mesh.elements[ep1[0]]
            e2 = self.F.mesh.elements[ep2[0]]
            # get edges for element points
            edgeList1 = e1.get_point_edge(ep1[1])
            edgeList2 = e2.get_point_edge(ep2[1])

            # ignore corners
            if len(edgeList1) > 1 or len(edgeList2) > 2:
                pass
            else:
                try:
                    # get edge number, point number (on edge), edge instance
                    (ei1, pi1, edge1) = edgeList1[0]
                    (ei2, pi2, edge2) = edgeList2[0]
                except IndexError:
                    pdb.set_trace()
                # ignore repeat element edges
                if ((ep1[0], ei1) in doneEdges) and ((ep2[0], ei2)) in doneEdges:
                    pass
                else:
                    # check if element edges are reversed
                    ### assumes no hanging nodes and same number of nodes per edge!!! ###

                    # get next node along edge
                    if pi1 != pi2:
                        eval1 = edge1.get_elem_coord(linspace(0.0, 1.0, D))
                        eval2 = edge2.get_elem_coord(linspace(1.0, 0.0, D))
                    else:
                        eval1 = edge1.get_elem_coord(linspace(0.0, 1.0, D))
                        eval2 = edge2.get_elem_coord(linspace(0.0, 1.0, D))

                    # get basis values for these
                    basis1 = [self.F.basis[e1.type].eval_derivatives(eval1.T, d).T for d in ((1, 0), (0, 1))]
                    basis2 = [self.F.basis[e2.type].eval_derivatives(eval2.T, d).T for d in ((1, 0), (0, 1))]

                    self.edgeEvalPoints.append((ep1[0], eval1, ep2[0], eval2))  # ( element1, ep1, element2, ep2 )
                    self.edgeEvalBasis.append(
                        (ep1[0], basis1, ep2[0], basis2))  # ( element1, basis1, element2, basis2 )
                    doneEdges.add((ep1[0], ei1))
                    doneEdges.add((ep2[0], ei2))
                    self.nPairs += eval1.shape[0]

        return

    def makeObj(self, D):
        """ make a lagrange multiplier element edge smoothing objective 
        function with each edge discretised at D
        """

        # calculate edge point basis values
        self._procEdge(D)

        # make derivatives evaluation matrices
        A1dxi1 = zeros((self.nPairs, self.F.get_number_of_ensemble_points()), dtype=float)
        A1dxi2 = zeros((self.nPairs, self.F.get_number_of_ensemble_points()), dtype=float)
        A2dxi1 = zeros((self.nPairs, self.F.get_number_of_ensemble_points()), dtype=float)
        A2dxi2 = zeros((self.nPairs, self.F.get_number_of_ensemble_points()), dtype=float)
        row1 = 0
        row2 = 0

        for e1, b1, e2, b2 in self.edgeEvalBasis:
            emap1 = self.el2en[e1]
            emap2 = self.el2en[e2]

            # make dxi1 and dxi2 matrices for points on one side
            for n in range(b1[0].shape[1]):
                A1dxi1[row1:row1 + b1[0].shape[0], emap1[n][0][0]] = b1[0][:, n]
                A1dxi2[row1:row1 + b1[1].shape[0], emap1[n][0][0]] = b1[1][:, n]

            row1 += b1[0].shape[0]

            # make dxi1 and dxi2 matrices for points on the other side
            for n in range(b2[0].shape[1]):
                A2dxi1[row2:row2 + b2[0].shape[0], emap2[n][0][0]] = b2[0][:, n]
                A2dxi2[row2:row2 + b2[1].shape[0], emap2[n][0][0]] = b2[1][:, n]

            row2 += b2[0].shape[0]

        # sparsify matrices
        sA1dxi1 = sparse.csc_matrix(A1dxi1)
        sA1dxi2 = sparse.csc_matrix(A1dxi2)
        sA2dxi1 = sparse.csc_matrix(A2dxi1)
        sA2dxi2 = sparse.csc_matrix(A2dxi2)

        def obj(x):
            P = x.reshape((3, -1)).T

            # evaluate normal on one side
            # evaluate  dxi1
            d1dxi1 = sA1dxi1 * P
            # evaluate  dxi2
            d1dxi2 = sA1dxi2 * P
            # cross product and normalise
            n1 = math.norms(cross(d1dxi1, d1dxi2))
            # n10, n11, n12 = n1.T

            # evaluation normal of the other side
            # evaluate  dxi1
            d2dxi1 = sA2dxi1 * P
            # evaluate  dxi2
            d2dxi2 = sA2dxi2 * P
            # cross product and normalise
            n2 = math.norms(cross(d2dxi1, d2dxi2))
            # n20, n21, n22 = n2.T

            # dot product normals
            err = 1.0 - (n1[:, 0] * n2[:, 0] + n1[:, 1] * n2[:, 1] + n1[:, 2] * n2[:, 2])
            # err = arccos(n1[:,0]*n2[:,0] + n1[:,1]*n2[:,1] + n1[:,2]*n2[:,2])
            # err = ne.evaluate( '1.0 - (n1[:,0]*n2[:,0] + n1[:,1]*n2[:,1] + n1[:,2]*n2[:,2])' )
            # err = ne.evaluate( '1.0 - (n10*n20 + n11*n21 + n12*n22)' )

            return err

        return obj


# ======================================================================#
class tangentSmoother2(object):
    """ penalises surface tangent (normal to element boundary) mismatch 
    at 2D element boundaries
    """

    def __init__(self, F):

        self.F = F
        self.en2el, self.el2en = self.F.get_mapping()  # ens2elem, elem2ens maps
        self._procMap()

    def _procMap(self):
        """ find element points of shared nodes
        """
        self.edgePoints = []
        self.vertexPoints = []  # not implemented

        # for each ensemble point
        for p in list(self.en2el.keys()):
            # edge nodes
            if len(self.en2el[p]) == 2:
                E = list(self.en2el[p].items())  # element points mapped to by ensemble point p
                self.edgePoints.append(((E[0][0], list(E[0][1].keys())[0]), (
                    E[1][0], list(E[1][1].keys())[0])))  # ((element1, elementnode1), (element2, elementnode2))
            # vertex nodes
            elif len(self.en2el[p]) > 2:
                pass

        return

    def _whichDeriv(self, elem, edge):
        """
        determine the derivative normal to the edge to calculate for a
        given element and edge.
        """
        quadMap = {0: ((0, 1), 1), 1: ((1, 0), -1), 2: ((0, 1), -1), 3: ((0, 1), 1)}
        triMap = {0: ((1, 0, 0), 1), 1: ((0, 1, 0), 1),
                  2: ((0, 0, 1), 1)}  # implement d wrt to area coordinates in basis

        if 'quad' in elem.type:
            return quadMap[edge]
        elif 'tri' in elem.type:
            return triMap[edge]

    def _procEdge(self, D):
        """ for each pair in edgePoints, find the element and edge shared
        and generate eval points along the shared edges, also record with
        element direction derivative needs to be calculated
        """
        self.edgeEvalBasis = []
        self.edgeEvalPoints = []
        self.edgeEvalXiDirection = []
        self.nPairs = 0  # number of pairs of edge points
        doneEdges = []  # list of element edges that have been done
        # for each pair
        for ep1, ep2 in self.edgePoints:
            # get elements
            e1 = self.F.mesh.elements[ep1[0]]
            e2 = self.F.mesh.elements[ep2[0]]
            # get edges for element points
            edgeList1 = e1.get_point_edge(ep1[1])
            edgeList2 = e2.get_point_edge(ep2[1])

            # ignore corners
            if len(edgeList1) > 1 or len(edgeList2) > 2:
                pass
            else:
                try:
                    (ei1, pi1, edge1) = edgeList1[0]
                    (ei2, pi2, edge2) = edgeList2[0]
                except IndexError:
                    pdb.set_trace()
                # ignore repeat element edges
                if ((ep1[0], ei1) in doneEdges) and ((ep2[0], ei2)) in doneEdges:
                    pass
                else:
                    # check if element edges are reversed
                    ### assumes no hanging nodes and same number of nodes per edge!!! ###
                    if pi1 != pi2:
                        eval1 = edge1.get_elem_coord(linspace(0.0, 1.0, D))
                        eval2 = edge2.get_elem_coord(linspace(1.0, 0.0, D))
                    else:
                        eval1 = edge1.get_elem_coord(linspace(0.0, 1.0, D))
                        eval2 = edge2.get_elem_coord(linspace(0.0, 1.0, D))

                    # get basis values for these
                    basis1 = [self.F.basis[e1.type].eval_derivatives(eval1.T, d).T for d in ((1, 0), (0, 1))]
                    basis2 = [self.F.basis[e2.type].eval_derivatives(eval2.T, d).T for d in ((1, 0), (0, 1))]

                    self.edgeEvalPoints.append((ep1[0], eval1, ep2[0], eval2))  # ( element1, ep1, element2, ep2 )
                    self.edgeEvalBasis.append(
                        (ep1[0], basis1, ep2[0], basis2))  # ( element1, basis1, element2, basis2 )
                    doneEdges.append((ep1[0], ei1))
                    doneEdges.append((ep2[0], ei2))
                    self.nPairs += eval1.shape[0]

        return

    def makeObj(self, D):
        """ make a lagrange multiplier element edge smoothing objective 
        function with each edge discretised at D
        """

        # calculate edge point basis values
        self._procEdge(D)

        # make derivatives evaluation matrices
        A1dxi1 = zeros((self.nPairs, self.F.get_number_of_ensemble_points()), dtype=float)
        A1dxi2 = zeros((self.nPairs, self.F.get_number_of_ensemble_points()), dtype=float)
        A2dxi1 = zeros((self.nPairs, self.F.get_number_of_ensemble_points()), dtype=float)
        A2dxi2 = zeros((self.nPairs, self.F.get_number_of_ensemble_points()), dtype=float)
        row1 = 0
        row2 = 0

        for e1, b1, e2, b2 in self.edgeEvalBasis:
            emap1 = self.el2en[e1]
            emap2 = self.el2en[e2]

            # make dxi1 and dxi2 matrices for points on one side
            for n in range(b1[0].shape[1]):
                A1dxi1[row1:row1 + b1[0].shape[0], emap1[n][0][0]] = b1[0][:, n]
                A1dxi2[row1:row1 + b1[1].shape[0], emap1[n][0][0]] = b1[1][:, n]

            row1 += b1[0].shape[0]

            # make dxi1 and dxi2 matrices for points on the other side
            for n in range(b2[0].shape[1]):
                A2dxi1[row2:row2 + b2[0].shape[0], emap2[n][0][0]] = b2[0][:, n]
                A2dxi2[row2:row2 + b2[1].shape[0], emap2[n][0][0]] = b2[1][:, n]

            row2 += b2[0].shape[0]

        # sparsify matrices
        sA1dxi1 = sparse.csc_matrix(A1dxi1)
        sA1dxi2 = sparse.csc_matrix(A1dxi2)
        sA2dxi1 = sparse.csc_matrix(A2dxi1)
        sA2dxi2 = sparse.csc_matrix(A2dxi2)

        def obj(x):
            P = x.reshape((3, -1)).T

            # evaluate normal on one side
            # evaluate  dxi1
            d1dxi1 = sA1dxi1 * P
            # evaluate  dxi2
            d1dxi2 = sA1dxi2 * P
            # cross product and normalise
            n1 = math.norms(cross(d1dxi1, d1dxi2))
            n10, n11, n12 = n1.T

            # evaluation normal of the other side
            # evaluate  dxi1
            d2dxi1 = sA2dxi1 * P
            # evaluate  dxi2
            d2dxi2 = sA2dxi2 * P
            # cross product and normalise
            n2 = math.norms(cross(d2dxi1, d2dxi2))
            n20, n21, n22 = n2.T

            # dot product normals
            # err = 1.0 - (n1[:,0]*n2[:,0] + n1[:,1]*n2[:,1] + n1[:,2]*n2[:,2])
            err = 1.0 - (n1 * n2).sum(1)
            # ~ err = ne.evaluate( '1.0 - (n1[:,0]*n2[:,0] + n1[:,1]*n2[:,1] + n1[:,2]*n2[:,2])' )
            # err = ne.evaluate( '1.0 - (n10*n20 + n11*n21 + n12*n22)' )

            return err

        return obj


# ======================================================================#
class tangentSmoother(object):
    """ for 1D curves
    """

    def __init__(self, F):

        self.F = F
        self.en2el, self.el2en = self.F.get_mapping()  # ens2elem, elem2ens maps
        self._procMap()

    def _procMap(self):
        """ find element points of shared nodes
        """
        self.edgePoints = []
        self.vertexPoints = []  # not implemented

        # for each ensemble point
        for p in list(self.en2el.keys()):
            # edge nodes
            if len(self.en2el[p]) == 2:
                E = list(self.en2el[p].items())  # element points mapped to by ensemble point p
                self.edgePoints.append(((E[0][0], list(E[0][1].keys())[0]), (
                    E[1][0], list(E[1][1].keys())[0])))  # ((element1, elementnode1), (element2, elementnode2))
            # vertex nodes
            elif len(self.en2el[p]) > 2:
                pass

        return

    def _procEdge(self):
        """ for each pair in edgePoints, find the element and calculate
        the basis for evaluation at the points
        """
        self.edgeEvalBasis = []
        self.edgeEvalPoints = []
        self.nPairs = 0  # number of pairs of edge points
        doneEdges = set()  # list of element edges that have been done
        # for each pair
        for ep1, ep2 in self.edgePoints:
            # skip edges already processed
            if (ep1, ep2) not in doneEdges:
                # get elements
                e1 = self.F.mesh.elements[ep1[0]]
                e2 = self.F.mesh.elements[ep2[0]]
                # get element coordinates of points
                try:
                    eval1 = e1.node_coordinates[ep1[1]]
                    eval2 = e2.node_coordinates[ep2[1]]
                except IndexError:
                    pdb.set_trace()

                # check direction of element coordinates
                if (eval1 - eval2) == 0.0:
                    coeff = -1.0
                else:
                    coeff = 1.0

                # calculate basis value at element coordinates
                basis1 = self.F.basis[e1.type].eval_derivatives(eval1, (1,))
                basis2 = coeff * self.F.basis[e2.type].eval_derivatives(eval2, (1,))

                self.edgeEvalPoints.append((ep1[0], eval1, ep2[0], eval2))  # ( element1, ep1, element2, ep2 )
                self.edgeEvalBasis.append((ep1[0], basis1, ep2[0], basis2))  # ( element1, basis1, element2, basis2 )

                doneEdges.add((ep1, ep2))
                self.nPairs += 1

        return

    def makeObj(self):
        """ make a lagrange multiplier element edge smoothing objective 
        function with each edge discretised at D
        """

        # calculate edge point basis values
        self._procEdge()

        if self.nPairs == 0:
            def obj(x):
                return 0.0

            return obj

        # make derivatives evaluation matrices
        A1dxi1 = zeros((self.nPairs, self.F.get_number_of_ensemble_points()), dtype=float)
        A2dxi1 = zeros((self.nPairs, self.F.get_number_of_ensemble_points()), dtype=float)
        row1 = 0
        row2 = 0

        for e1, b1, e2, b2 in self.edgeEvalBasis:
            emap1 = self.el2en[e1]
            emap2 = self.el2en[e2]

            # make dxi1 matrix for points on one side
            for n in range(b1.shape[0]):
                A1dxi1[row1, emap1[n][0][0]] = b1[n]
            row1 += 1

            # make dxi1 and dxi2 matrices for points on the other side  
            for n in range(b2.shape[0]):
                A2dxi1[row2, emap2[n][0][0]] = b2[n]
            row2 += 1

        # sparsify matrices
        sA1dxi1 = sparse.csc_matrix(A1dxi1)
        sA2dxi1 = sparse.csc_matrix(A2dxi1)

        def obj(x):
            P = x.reshape((3, -1)).T

            # evaluate tangent on one side
            # evaluate  dxi1
            d1dxi1 = math.norms(sA1dxi1 * P)

            # evaluation normal of the other side
            # evaluate  dxi1
            d2dxi1 = math.norms(sA2dxi1 * P)

            # dot product normals
            err = 1.0 - (d1dxi1 * d2dxi1).sum(1)

            return err

        return obj


class normalSmoother(object):
    """ 
    DEPRECATED
    """
    simplex_cubic_lagrange_2d_element_node = {0: (0.0, 0.0), 1: (1 / 3.0, 0.0), 2: (2 / 3.0, 0.0), 3: (1.0, 0.0),
                                              4: (0.0, 1 / 3.0), 5: (1 / 3.0, 1 / 3.0), 6: (2 / 3.0, 1 / 3.0),
                                              7: (0.0, 2 / 3.0), 8: (1 / 3.0, 2 / 3.0),
                                              9: (0.0, 1.0)}
    simplex_eval_points = array([(0.0, 0.0), (1 / 3.0, 0.0), (2 / 3.0, 0.0), (1.0, 0.0),
                                 (0.0, 1 / 3.0), (1 / 3.0, 1 / 3.0), (2 / 3.0, 1 / 3.0),
                                 (0.0, 2 / 3.0), (1 / 3.0, 2 / 3.0),
                                 (0.0, 1.0)])

    quad_eval_points = array([(0.0, 0.0), (1 / 3.0, 0.0), (2 / 3.0, 0.0), (1.0, 0.0),
                              (0.0, 1 / 3.0), (1 / 3.0, 1 / 3.0), (2 / 3.0, 1 / 3.0), (1.0, 1 / 3.0),
                              (0.0, 2 / 3.0), (1 / 3.0, 2 / 3.0), (2 / 3.0, 2 / 3.0), (1.0, 2 / 3.0),
                              (0.0, 1.0), (1 / 3.0, 1.0), (2 / 3.0, 1.0), (1.0, 1.0)])

    def __init__(self, G, nW, cW, d):
        """ old method: calculate normals at all edge element nodes,
        then use self.edgePoints to calculate the dot products of share
        nodes. 
        
        new method: first generate self.edgeEvalBasis using self.edgePoints
        (do once). Then per iteration, use self.edgeEvalBasis to calculate
        normals and dot product of normal at shared points along shared
        edges.
        """
        self.G = G
        self.F = self.G.ensemble_field_function
        self.nW = nW
        self.cW = cW
        self.d = d
        # ~ self.basis_values = {'tri10':  [ F.basis['tri10'].eval_derivatives( self.simplex_eval_points.T, d ).T for d in ((1,0),(0,1)) ],\
        # ~ 'quad44': [ F.basis['quad44'].eval_derivatives( self.quad_eval_points.T, d ).T for d in ((1,0),(0,1)) ]}
        self.map = self.F.get_mapping()[0]
        self._procMap()
        self._procEdge()

    def _procMap(self):
        """ find element points of shared nodes
        """
        self.edgePoints = []
        self.vertexPoints = []  # not implemented

        for p in list(self.map.keys()):
            # edge nodes
            if len(self.map[p]) == 2:
                E = list(self.map[p].items())
                self.edgePoints.append(((E[0][0], list(E[0][1].keys())[0]), (E[1][0], list(E[1][1].keys())[0])))
            elif len(self.map[p]) > 2:
                pass

        return

    def _procEdge(self):
        """ for each pair in edgePoints, find the element and edge shared
        and generate eval points along the shared edges
        """
        self.edgeEvalBasis = []
        self.edgeEvalPoints = []
        doneEdges = []  # list of element edges that have been done
        # for each pair
        for ep1, ep2 in self.edgePoints:
            # get elements
            e1 = self.F.mesh.elements[ep1[0]]
            e2 = self.F.mesh.elements[ep2[0]]
            # get edges for element points
            edgeList1 = e1.get_point_edge(ep1[1])
            edgeList2 = e2.get_point_edge(ep2[1])

            # ignore corners
            if len(edgeList1) > 1 or len(edgeList2) > 2:
                pass
            else:
                try:
                    (ei1, pi1, edge1) = edgeList1[0]
                    (ei2, pi2, edge2) = edgeList2[0]
                except IndexError:
                    pdb.set_trace()
                # ignore repeat element edges
                if ((ep1[0], ei1) in doneEdges) and ((ep2[0], ei2)) in doneEdges:
                    pass
                else:
                    # check if element edges are reversed
                    ### assumes no hanging nodes and same number of nodes per edge!!! ###
                    if pi1 != pi2:
                        eval1 = edge1.get_elem_coord(linspace(0.0, 1.0, self.d))
                        eval2 = edge2.get_elem_coord(linspace(1.0, 0.0, self.d))
                    else:
                        eval1 = edge1.get_elem_coord(linspace(0.0, 1.0, self.d))
                        eval2 = edge2.get_elem_coord(linspace(0.0, 1.0, self.d))

                    # get basis values for these
                    basis1 = [self.F.basis[e1.type].eval_derivatives(eval1.T, d).T for d in ((1, 0), (0, 1))]
                    basis2 = [self.F.basis[e2.type].eval_derivatives(eval2.T, d).T for d in ((1, 0), (0, 1))]

                    self.edgeEvalPoints.append((ep1[0], eval1, ep2[0], eval2))
                    self.edgeEvalBasis.append((ep1[0], basis1, ep2[0], basis2))
                    doneEdges.append((ep1[0], ei1))
                    doneEdges.append((ep2[0], ei2))

        return

    def smoothObj(self, params):
        """ add penalty for difference in normal direction at shared nodes
        between elements
        """
        # ~ n = self._calcNormal( params )
        # ~ d = self._compareNormals( n ) * self.W
        # ~ print 'smoothing rms:', sqrt( (d**2.0).mean() )
        dNorm = None
        dCurvature = None
        DD1 = []
        DD2 = []
        # evaluate derivatives of coordinate fields
        for P in params:
            self.F.set_parameters(P)
            D1 = None
            D2 = None
            for e1, ep1, e2, ep2 in self.edgeEvalPoints:
                d1 = self.F.evaluate_field_at_element_point(e1, ep1, parameters=P, derivs=-1)[1]
                d2 = self.F.evaluate_field_at_element_point(e2, ep2, parameters=P, derivs=-1)[1]
                try:
                    D1 = hstack((D1, d1))
                    D2 = hstack((D2, d2))
                except ValueError:
                    D1 = d1.copy()
                    D2 = d2.copy()

            DD1.append(D1)
            DD2.append(D2)

        DD1 = array(DD1)
        # ~ print 'dd1 ', DD1.shape
        DD2 = array(DD2)
        # ~ print 'dd2 ', DD2.shape
        # ~ print 'input ', dstack( (D1,D2) ).shape

        # penalise diff in curvature on edges
        # ~ H1 = self.G._calculate_curvature( DD1 )[1]
        # ~ H2 = self.G._calculate_curvature( DD2 )[1]
        # ~ dCurvature = ( H1 - H2 )**2.0 * self.cW
        # penalise high curavture on edges
        if self.cW:
            dCurvature = self.G._calculate_curvature(dstack((DD1, DD2)))[1] ** 2.0 * self.cW
        else:
            dCurvature = array([0.])
        # penalise diff in normal on edge
        if self.nW:
            dNorm = self._calcNormal3(DD1[:, 0, :], DD1[:, 1, :], DD2[:, 0, :], DD2[:, 1, :]) * self.nW
        else:
            dNorm = array([0.])

        return dNorm, dCurvature
        # ~ return dNorm, array([0.0,])

    def _calcCurvatureDiff(self, params):
        dH = []
        self.G.set_field_parameters(params)
        for e1, ep1, e2, ep2 in self.edgeEvalPoints:
            H1 = self.G.evaluate_curvature_at_element_points(e1, ep1)[1]
            H2 = self.G.evaluate_curvature_at_element_points(e2, ep2)[1]
            dH.append(H1 - H2)

        d = array(dH).ravel()
        return d * d

    def _calcNormal3(self, D1dx1, D1dx2, D2dx1, D2dx2):
        # calculate normals
        N1 = math.norms(cross(D1dx1.T, D1dx2.T))
        N2 = math.norms(cross(D2dx1.T, D2dx2.T))

        D = (N1[:, 0] * N2[:, 0] + N1[:, 1] * N2[:, 1] + N1[:, 2] * N2[:, 2]) - 1.0

        # calc dot of normals
        # ~ N = array((N1,N2)).swapaxes(0,1)
        # ~ D = array( [ 1.0-dot(n1, n2) for n1,n2 in N ] )**2.0
        return D * D

    def _calcNormal2(self, params):
        # evaluate normals at precalculated points on shared edges

        # calculate derivatives
        D1dx1 = []  # store derivatives for all p
        D1dx2 = []
        D2dx1 = []
        D2dx2 = []
        # cycle params last to take advantage of params caching
        for p in params:
            self.F.set_parameters(p)
            d1 = []  # store derivatives for current p
            d2 = []  # store derivatives for current p
            for e1, b1, e2, b2 in self.edgeEvalBasis:
                eParams1 = self.F._get_element_parameters(e1)
                eParams2 = self.F._get_element_parameters(e2)
                d1.append(dot(b1, eParams1))
                d2.append(dot(b2, eParams2))

            d1 = array(d1)
            d2 = array(d2)
            D1dx1.append(d1[:, 0, :].ravel())
            D1dx2.append(d1[:, 1, :].ravel())
            D2dx1.append(d2[:, 0, :].ravel())
            D2dx2.append(d2[:, 1, :].ravel())

        D1dx1 = array(D1dx1)
        D1dx2 = array(D1dx2)
        D2dx1 = array(D2dx1)
        D2dx2 = array(D2dx2)

        # calculate normals
        N1 = math.norms(cross(D1dx1.T, D1dx2.T))
        N2 = math.norms(cross(D2dx1.T, D2dx2.T))
        N = array((N1, N2)).swapaxes(0, 1)

        # calc dot of normals
        D = array([1.0 - dot(n1, n2) for n1, n2 in N])

        return D * D

    def _calcNormal(self, params):
        # evaluate normals at element boundary nodes
        normals = []
        for e in list(self.F.mesh.elements.keys()):
            basis = self.basis_values[self.F.mesh.elements[e].type]
            D = []
            for p in params:
                self.F.set_parameters(p)
                e_params = self.F._get_element_parameters(e)
                D.append(dot(basis, e_params))

            D = array(D)
            normals.append(math.norms(cross(D[:, 0].T, D[:, 1].T)))

        # ~ pdb.set_trace()
        # ~ n = array(normals)

        return normals

    def _compareNormals(self, N):

        D = zeros(len(self.edgePoints))
        i = 0
        for (p1, p2) in self.edgePoints:
            D[i] = 1.0 - dot(N[p1[0]][p1[1]], N[p2[0]][p2[1]])
            i += 1

        return D ** 2.0


# def norm(x):
#     return x/(sqrt((x**2.0).sum(1))[:,newaxis])
# ~ return ne.evaluate('x/sqrt(sum(x*x, axis=1))')[:,newaxis]

# penalty functions for meshing fitting obj

# element area ratio
class elementAreaPenalty1(object):

    def __init__(self, F, evalD, p0, w):

        self.F = F
        self.evalD = evalD
        self.calcAR0(p0)
        self.w = w

    def calcAR0(self, p0):
        # evaluate 1st derivatives on eval_grid over mesh
        D1P1 = self.F.evaluate_derivatives_in_mesh(self.evalD, parameters=p0[0], derivs=(1, 0), unpack=False)
        D2P1 = self.F.evaluate_derivatives_in_mesh(self.evalD, parameters=p0[0], derivs=(0, 1), unpack=False)
        D1P2 = self.F.evaluate_derivatives_in_mesh(self.evalD, parameters=p0[1], derivs=(1, 0), unpack=False)
        D2P2 = self.F.evaluate_derivatives_in_mesh(self.evalD, parameters=p0[1], derivs=(0, 1), unpack=False)
        D1P3 = self.F.evaluate_derivatives_in_mesh(self.evalD, parameters=p0[2], derivs=(1, 0), unpack=False)
        D2P3 = self.F.evaluate_derivatives_in_mesh(self.evalD, parameters=p0[2], derivs=(0, 1), unpack=False)

        aElem = zeros(D1P1.__len__())
        for e in range(aElem.shape[0]):
            aElem[e] = cross(array([D1P1[e], D1P2[e], D1P3[e]]).T, array([D2P1[e], D2P2[e], D2P3[e]]).T).sum()

        self.aR0 = aElem / aElem.max()

        log.debug('max element area:', aElem.max())
        return

    def obj(self, params):
        """ for each element, calculate its area, and the ratio of its area
        to the largest element area/mean element area. Returns the square
        difference of the current ratio to the initial ratio as error to 
        minimise
        """

        D1P1 = self.F.evaluate_derivatives_in_mesh(self.evalD, parameters=params[0], derivs=(1, 0), unpack=False)
        D2P1 = self.F.evaluate_derivatives_in_mesh(self.evalD, parameters=params[0], derivs=(0, 1), unpack=False)
        D1P2 = self.F.evaluate_derivatives_in_mesh(self.evalD, parameters=params[1], derivs=(1, 0), unpack=False)
        D2P2 = self.F.evaluate_derivatives_in_mesh(self.evalD, parameters=params[1], derivs=(0, 1), unpack=False)
        D1P3 = self.F.evaluate_derivatives_in_mesh(self.evalD, parameters=params[2], derivs=(1, 0), unpack=False)
        D2P3 = self.F.evaluate_derivatives_in_mesh(self.evalD, parameters=params[2], derivs=(0, 1), unpack=False)

        aElem = zeros(D1P1.__len__())
        for e in range(aElem.shape[0]):
            aElem[e] = cross(array([D1P1[e], D1P2[e], D1P3[e]]).T, array([D2P1[e], D2P2[e], D2P3[e]]).T).sum()

        # get largest area and divide the rest by it
        d = self.aR0 - aElem / aElem.max()
        d = self.w * d * d
        log.debug('area rms', sqrt(d.sum()))

        return d


# element vertex angles
class elementDotPenalty2D(object):
    """ penalise deviation from the initial mesh angle field (dot(dxi1,dxi2)
    """

    def __init__(self, F, evalD, p0, w):
        self.F = F
        self.evalD = evalD
        self.calcDot0(p0)
        self.w = w

    def calcDot0(self, p0):
        # evaluate 1st derivatives on eval_grid over mesh
        D1P1 = self.F.evaluate_derivatives_in_mesh(self.evalD, parameters=p0[0], derivs=(1, 0), unpack=True)
        D2P1 = self.F.evaluate_derivatives_in_mesh(self.evalD, parameters=p0[0], derivs=(0, 1), unpack=True)
        D1P2 = self.F.evaluate_derivatives_in_mesh(self.evalD, parameters=p0[1], derivs=(1, 0), unpack=True)
        D2P2 = self.F.evaluate_derivatives_in_mesh(self.evalD, parameters=p0[1], derivs=(0, 1), unpack=True)
        D1P3 = self.F.evaluate_derivatives_in_mesh(self.evalD, parameters=p0[2], derivs=(1, 0), unpack=True)
        D2P3 = self.F.evaluate_derivatives_in_mesh(self.evalD, parameters=p0[2], derivs=(0, 1), unpack=True)

        self.dot0 = D1P1 * D2P1 + D1P2 * D2P2 + D1P3 * D2P3
        return

    def obj(self, params):
        # evaluate 1st derivatives on eval_grid over mesh
        D1P1 = self.F.evaluate_derivatives_in_mesh(self.evalD, parameters=params[0], derivs=(1, 0), unpack=True)
        D2P1 = self.F.evaluate_derivatives_in_mesh(self.evalD, parameters=params[0], derivs=(0, 1), unpack=True)
        D1P2 = self.F.evaluate_derivatives_in_mesh(self.evalD, parameters=params[1], derivs=(1, 0), unpack=True)
        D2P2 = self.F.evaluate_derivatives_in_mesh(self.evalD, parameters=params[1], derivs=(0, 1), unpack=True)
        D1P3 = self.F.evaluate_derivatives_in_mesh(self.evalD, parameters=params[2], derivs=(1, 0), unpack=True)
        D2P3 = self.F.evaluate_derivatives_in_mesh(self.evalD, parameters=params[2], derivs=(0, 1), unpack=True)

        dot0 = D1P1 * D2P1 + D1P2 * D2P2 + D1P3 * D2P3
        d = self.dot0 - dot0
        d = self.w * d * d
        return d


def makeElementDotPenalty2D(G, evalD, p0, w):
    gDEval = geometric_field.makeGeometricFieldDerivativesEvaluatorSparse(G, evalD)

    D = gDEval(p0.ravel())
    D1P10, D2P10, D11P10, D22P10, D12P10 = D[0]
    D1P20, D2P20, D11P20, D22P20, D12P20 = D[1]
    D1P30, D2P30, D11P30, D22P30, D12P30 = D[2]
    dot0 = D1P10 * D2P10 + D1P20 * D2P20 + D1P30 * D2P30

    def obj(p):
        D = gDEval(p)
        D1P1, D2P1, D11P1, D22P1, D12P1 = D[0]
        D1P2, D2P2, D11P2, D22P2, D12P2 = D[1]
        D1P3, D2P3, D11P3, D22P3, D12P3 = D[2]
        dotp = D1P1 * D2P1 + D1P2 * D2P2 + D1P3 * D2P3
        d = dot0 - dotp
        return d * d * w

    return obj


class elementDotPenalty3D(object):

    def __init__(self, F, evalD, p0, w):
        self.F = F
        self.evalD = evalD
        self.w = w
        self.calcDot0(p0)

    def calcDot0(self, p):
        # evaluate derivatives on eval_grid over mesh
        self.F.set_parameters(p[0])
        D1P1 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(1, 0, 0), unpack=True)
        D2P1 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(0, 1, 0), unpack=True)
        D3P1 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(0, 0, 1), unpack=True)
        self.F.set_parameters(p[1])
        D1P2 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(1, 0, 0), unpack=True)
        D2P2 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(0, 1, 0), unpack=True)
        D3P2 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(0, 0, 1), unpack=True)
        self.F.set_parameters(p[2])
        D1P3 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(1, 0, 0), unpack=True)
        D2P3 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(0, 1, 0), unpack=True)
        D3P3 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(0, 0, 1), unpack=True)

        self.dot012 = D1P1 * D2P1 + D1P2 * D2P2 + D1P3 * D2P3
        self.dot013 = D1P1 * D3P1 + D1P2 * D3P2 + D1P3 * D3P3
        self.dot023 = D2P1 * D3P1 + D2P2 * D3P2 + D2P3 * D3P3

    def obj(self, p):
        # evaluate derivatives on eval_grid over mesh
        self.F.set_parameters(p[0])
        D1P1 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(1, 0, 0), unpack=True)
        D2P1 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(0, 1, 0), unpack=True)
        D3P1 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(0, 0, 1), unpack=True)
        self.F.set_parameters(p[1])
        D1P2 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(1, 0, 0), unpack=True)
        D2P2 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(0, 1, 0), unpack=True)
        D3P2 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(0, 0, 1), unpack=True)
        self.F.set_parameters(p[2])
        D1P3 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(1, 0, 0), unpack=True)
        D2P3 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(0, 1, 0), unpack=True)
        D3P3 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(0, 0, 1), unpack=True)

        dot12 = D1P1 * D2P1 + D1P2 * D2P2 + D1P3 * D2P3
        dot13 = D1P1 * D3P1 + D1P2 * D3P2 + D1P3 * D3P3
        dot23 = D2P1 * D3P1 + D2P2 * D3P2 + D2P3 * D3P3

        d1 = self.w * (self.dot012 - dot12)
        d2 = self.w * (self.dot013 - dot13)
        d3 = self.w * (self.dot023 - dot23)

        return d1 * d1 + d2 * d2 + d3 * d3


class sobelovPenalty3D(object):

    def __init__(self, F, evalD, p0, w):
        self.F = F
        self.evalD = evalD
        self.w = w

    def obj(self, p):
        # evaluate derivatives on eval_grid over mesh
        self.F.set_parameters(p[0])
        D1P1 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(1, 0, 0), unpack=True)
        D11P1 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(2, 0, 0), unpack=True)
        D2P1 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(0, 1, 0), unpack=True)
        D22P1 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(0, 2, 0), unpack=True)
        D3P1 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(0, 0, 1), unpack=True)
        D33P1 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(0, 0, 2), unpack=True)
        self.F.set_parameters(p[1])
        D1P2 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(1, 0, 0), unpack=True)
        D11P2 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(2, 0, 0), unpack=True)
        D2P2 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(0, 1, 0), unpack=True)
        D22P2 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(0, 2, 0), unpack=True)
        D3P2 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(0, 0, 1), unpack=True)
        D33P2 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(0, 0, 2), unpack=True)
        self.F.set_parameters(p[2])
        D1P3 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(1, 0, 0), unpack=True)
        D11P3 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(2, 0, 0), unpack=True)
        D2P3 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(0, 1, 0), unpack=True)
        D22P3 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(0, 2, 0), unpack=True)
        D3P3 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(0, 0, 1), unpack=True)
        D33P3 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=(0, 0, 2), unpack=True)


class sobelovPenalty2D(object):

    def __init__(self, F, evalD, w):
        self.F = F
        self.evalD = evalD
        self.w = w

    def obj(self, p):
        # evaluate derivatives on eval_grid over mesh
        self.F.set_parameters(p[0])
        D1P1, D2P1, D11P1, D22P1, D12P1 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=-1, unpack=True)
        self.F.set_parameters(p[1])
        D1P2, D2P2, D11P2, D22P2, D12P2 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=-1, unpack=True)
        self.F.set_parameters(p[2])
        D1P3, D2P3, D11P3, D22P3, D12P3 = self.F.evaluate_derivatives_in_mesh(self.evalD, derivs=-1, unpack=True)

        S = self.w[0] * (D1P1 * D1P1 + D1P2 * D1P2 + D1P3 * D1P3) + \
            self.w[1] * (D11P1 * D11P1 + D11P2 * D11P2 + D11P3 * D11P3) + \
            self.w[2] * (D2P1 * D2P1 + D2P2 * D2P2 + D2P3 * D2P3) + \
            self.w[3] * (D22P1 * D22P1 + D22P2 * D22P2 + D22P3 * D22P3) + \
            self.w[4] * (D12P1 * D12P1 + D12P2 * D12P2 + D12P3 * D12P3)

        return S


def makeSobelovPenalty3D(G, evalD, w):
    gDEval = geometric_field.makeGeometricFieldDerivativesEvaluatorSparse(G, evalD)

    def obj(p):
        D = gDEval(p)
        # ~ D1P1, D2P1, D11P1, D22P1, D12P1 = D[0]
        # ~ D1P2, D2P2, D11P2, D22P2, D12P2 = D[1]
        # ~ D1P3, D2P3, D11P3, D22P3, D12P3 = D[2]

        # ~ pdb.set_trace()
        S = dot(w, (D ** 2.0).sum(0))

        # ~ S = w[0]*(D1P1*D1P1 + D1P2*D1P2 + D1P3*D1P3) +\
        # ~ w[1]*(D11P1*D11P1 + D11P2*D11P2 + D11P3*D11P3) +\
        # ~ w[2]*(D2P1*D2P1 + D2P2*D2P2 + D2P3*D2P3) +\
        # ~ w[3]*(D22P1*D22P1 + D22P2*D22P2 + D22P3*D22P3) +\
        # ~ w[4]*(D12P1*D12P1 + D12P2*D12P2 + D12P3*D12P3)

        return S

    return obj


def makeSobelovPenalty2D(G, evalD, w):
    gDEval = geometric_field.makeGeometricFieldDerivativesEvaluatorSparse(G, evalD)

    def obj(p):

        D = gDEval(p)
        # ~ D1P1, D2P1, D11P1, D22P1, D12P1 = D[0]
        # ~ D1P2, D2P2, D11P2, D22P2, D12P2 = D[1]
        # ~ D1P3, D2P3, D11P3, D22P3, D12P3 = D[2]

        # ~ S = w[0]*(D1P1*D1P1 + D1P2*D1P2 + D1P3*D1P3) +\
        # ~ w[1]*(D11P1*D11P1 + D11P2*D11P2 + D11P3*D11P3) +\
        # ~ w[2]*(D2P1*D2P1 + D2P2*D2P2 + D2P3*D2P3) +\
        # ~ w[3]*(D22P1*D22P1 + D22P2*D22P2 + D22P3*D22P3) +\
        # ~ w[4]*(D12P1*D12P1 + D12P2*D12P2 + D12P3*D12P3)

        try:
            S = (w * ((D * D).sum(0).T)).sum(1)
            # S = (w*D0*D0 + w*D1*D1 + w*D2*D2).sum(1)
            # S = ne.evaluate( 'sum(w*D0*D0 + w*D1*D1 + w*D2*D2, axis=1)', local_dict={'w':w, 'D0':D[0].T, 'D1':D[1].T, 'D2':D[2].T} )
        except ValueError:
            pdb.set_trace()

        return S

    return obj


def makeSobelovPenalty1D(G, evalD, w):
    gDEval = geometric_field.makeGeometricFieldDerivativesEvaluatorSparse(G, evalD)

    def obj(p):
        D = gDEval(p)
        D1P1, D11P1 = D[0]
        D1P2, D11P2 = D[1]
        D1P3, D11P3 = D[2]

        S = w[0] * (D1P1 * D1P1 + D1P2 * D1P2 + D1P3 * D1P3) + \
            w[1] * (D11P1 * D11P1 + D11P2 * D11P2 + D11P3 * D11P3)

        return S

    return obj


# ======================================================================#
# Host mesh fitting functions                                          #
# ======================================================================#
def makeQuadraticCube(xRange, yRange, zRange):
    F = EFF.ensemble_field_function('host', 3, debug=0)
    F.set_basis({'quad333': 'quad_L2_L2_L2'})
    F.set_new_mesh('host')
    F.create_elements('quad333', 1)
    F.map_parameters()
    G = geometric_field.geometric_field('host', 3, ensemble_field_function=F)

    x = list(linspace(xRange[0], xRange[1], 3)) * 9
    x = array(x)

    y = list(linspace(yRange[0], yRange[1], 3).repeat(3)) * 3
    y = array(y)

    z = linspace(zRange[0], zRange[1], 3).repeat(9)

    G.set_field_parameters([x[:, newaxis], y[:, newaxis], z[:, newaxis]])
    return G


def makeQuadraticCubeMulti(xRange, yRange, zRange, discretisation):
    F = EFF.ensemble_field_function('host', 3, debug=0)
    F.set_basis({'quad333': 'quad_L2_L2_L2'})
    F.set_new_mesh('host')
    G = geometric_field.geometric_field('host', 3, ensemble_field_function=F)

    # generate global nodes
    x = linspace(xRange[0], xRange[1], 3 + (discretisation[0] - 1) * 2)
    y = linspace(yRange[0], yRange[1], 3 + (discretisation[1] - 1) * 2)
    z = linspace(zRange[0], zRange[1], 3 + (discretisation[2] - 1) * 2)

    for ex in range(discretisation[0]):
        for ey in range(discretisation[1]):
            for ez in range(discretisation[2]):
                elem = element_types.create_element('quad333')
                if ex == 0:
                    elemCoordsX = list(x[0:3]) * 9
                else:
                    elemCoordsX = list(x[ex * 2:ex * 2 + 3]) * 9

                if ey == 0:
                    elemCoordsY = list(y[0:3].repeat(3)) * 3
                else:
                    elemCoordsY = list(y[ey * 2:ey * 2 + 3].repeat(3)) * 3

                if ez == 0:
                    elemCoordsZ = z[0:3].repeat(9)
                else:
                    elemCoordsZ = z[ez * 2:ez * 2 + 3].repeat(9)

                elemParams = array([elemCoordsX,
                                    elemCoordsY,
                                    elemCoordsZ])[:, :, newaxis]

                G.add_element_with_parameters(elem, elemParams)

    return G


def makeCubicCube(xRange, yRange, zRange):
    F = EFF.ensemble_field_function('host', 3, debug=0)
    F.set_basis({'quad444': 'quad_L3_L3_L3'})
    F.set_new_mesh('host')
    F.create_elements('quad444', 1)
    F.map_parameters()
    G = geometric_field.geometric_field('host', 3, ensemble_field_function=F)

    x = list(linspace(xRange[0], xRange[1], 4)) * 16
    x = array(x)

    y = list(linspace(yRange[0], yRange[1], 4).repeat(4)) * 4
    y = array(y)

    z = linspace(zRange[0], zRange[1], 4).repeat(16)

    G.set_field_parameters([x[:, newaxis], y[:, newaxis], z[:, newaxis]])
    return G


def makeCubicCubeMulti(xRange, yRange, zRange, discretisation):
    F = EFF.ensemble_field_function('host', 3, debug=0)
    F.set_basis({'quad444': 'quad_L3_L3_L3'})
    F.set_new_mesh('host')
    G = geometric_field.geometric_field('host', 3, ensemble_field_function=F)

    # generate global nodes
    x = linspace(xRange[0], xRange[1], 4 + (discretisation[0] - 1) * 3)
    y = linspace(yRange[0], yRange[1], 4 + (discretisation[1] - 1) * 3)
    z = linspace(zRange[0], zRange[1], 4 + (discretisation[2] - 1) * 3)

    for ex in range(discretisation[0]):
        for ey in range(discretisation[1]):
            for ez in range(discretisation[2]):
                elem = element_types.create_element('quad444')
                if ex == 0:
                    elemCoordsX = list(x[0:4]) * 16
                else:
                    elemCoordsX = list(x[ex * 3:ex * 3 + 4]) * 16

                if ey == 0:
                    elemCoordsY = list(y[0:4].repeat(4)) * 4
                else:
                    elemCoordsY = list(y[ey * 3:ey * 3 + 4].repeat(4)) * 4

                if ez == 0:
                    elemCoordsZ = list(z[0:4].repeat(16))
                else:
                    elemCoordsZ = list(z[ez * 3:ez * 3 + 4].repeat(16))

                elemParams = array([elemCoordsX,
                                    elemCoordsY,
                                    elemCoordsZ])[:, :, newaxis]

                G.add_element_with_parameters(elem, elemParams)

    return G


def makeQuarticCube(xRange, yRange, zRange):
    F = EFF.ensemble_field_function('host', 3, debug=0)
    F.set_basis({'quad555': 'quad_L4_L4_L4'})
    F.set_new_mesh('host')
    F.create_elements('quad555', 1)
    F.map_parameters()
    G = geometric_field.geometric_field('host', 3, ensemble_field_function=F)

    x = list(linspace(xRange[0], xRange[1], 5)) * 25
    x = array(x)

    y = list(linspace(yRange[0], yRange[1], 5).repeat(5)) * 5
    y = array(y)

    z = linspace(zRange[0], zRange[1], 5).repeat(25)

    G.set_field_parameters([x[:, newaxis], y[:, newaxis], z[:, newaxis]])
    return G


def makeQuarticCubeMulti(xRange, yRange, zRange, discretisation):
    F = EFF.ensemble_field_function('host', 3, debug=0)
    F.set_basis({'quad555': 'quad_L4_L4_L4'})
    F.set_new_mesh('host')
    G = geometric_field.geometric_field('host', 3, ensemble_field_function=F)

    # generate global nodes
    x = linspace(xRange[0], xRange[1], 5 + (discretisation[0] - 1) * 4)
    y = linspace(yRange[0], yRange[1], 5 + (discretisation[1] - 1) * 4)
    z = linspace(zRange[0], zRange[1], 5 + (discretisation[2] - 1) * 4)

    for ex in range(discretisation[0]):
        for ey in range(discretisation[1]):
            for ez in range(discretisation[2]):
                elem = element_types.create_element('quad555')
                if ex == 0:
                    elemCoordsX = list(x[0:5]) * 25
                else:
                    elemCoordsX = list(x[ex * 4:ex * 4 + 5]) * 25

                if ey == 0:
                    elemCoordsY = list(y[0:5].repeat(5)) * 5
                else:
                    elemCoordsY = list(y[ey * 4:ey * 4 + 5].repeat(5)) * 5

                if ez == 0:
                    elemCoordsZ = z[0:5].repeat(25)
                else:
                    elemCoordsZ = z[ez * 4:ez * 4 + 5].repeat(25)

                elemParams = array([elemCoordsX,
                                    elemCoordsY,
                                    elemCoordsZ])[:, :, newaxis]

                G.add_element_with_parameters(elem, elemParams)

    return G


def makeHostMesh(P, pad, elemType):
    hosts = {'quad333': makeQuadraticCube,
             'quad444': makeCubicCube,
             'quad555': makeQuarticCube}

    xmin = P.squeeze()[0].min() - pad
    xmax = P.squeeze()[0].max() + pad
    ymin = P.squeeze()[1].min() - pad
    ymax = P.squeeze()[1].max() + pad
    zmin = P.squeeze()[2].min() - pad
    zmax = P.squeeze()[2].max() + pad
    return hosts[elemType]([xmin, xmax], [ymin, ymax], [zmin, zmax])


def makeHostMeshMulti(P, pad, elemType, discretisation):
    hosts = {'quad333': makeQuadraticCubeMulti,
             'quad444': makeCubicCubeMulti,
             'quad555': makeQuarticCubeMulti}

    xmin = P.squeeze()[0].min() - pad
    xmax = P.squeeze()[0].max() + pad
    ymin = P.squeeze()[1].min() - pad
    ymax = P.squeeze()[1].max() + pad
    zmin = P.squeeze()[2].min() - pad
    zmax = P.squeeze()[2].max() + pad
    return hosts[elemType]([xmin, xmax], [ymin, ymax], [zmin, zmax], discretisation)


def viewHostMesh(m, d, figure):
    scale = 0.5
    radius = 1.0
    edges = m.ensemble_field_function.mesh.elements[0].edge_points
    v = m.evaluate_geometric_field(d)
    if figure is None:
        figure = geometric_field.mlab.figure()

    geometric_field.mlab.points3d(v[0], v[1], v[2], scale_factor=scale, mode='point', figure=figure)
    p = m.field_parameters[:, :, 0].T
    geometric_field.mlab.plot3d(p[edges[0], 0], p[edges[0], 1], p[edges[0], 2], tube_radius=radius, tube_sides=3,
                                figure=figure)
    geometric_field.mlab.plot3d(p[edges[1], 0], p[edges[1], 1], p[edges[1], 2], tube_radius=radius, tube_sides=3,
                                figure=figure)
    geometric_field.mlab.plot3d(p[edges[2], 0], p[edges[2], 1], p[edges[2], 2], tube_radius=radius, tube_sides=3,
                                figure=figure)
    geometric_field.mlab.plot3d(p[edges[3], 0], p[edges[3], 1], p[edges[3], 2], tube_radius=radius, tube_sides=3,
                                figure=figure)
    geometric_field.mlab.plot3d(p[edges[4], 0], p[edges[4], 1], p[edges[4], 2], tube_radius=radius, tube_sides=3,
                                figure=figure)
    geometric_field.mlab.plot3d(p[edges[5], 0], p[edges[5], 1], p[edges[5], 2], tube_radius=radius, tube_sides=3,
                                figure=figure)
    geometric_field.mlab.plot3d(p[edges[6], 0], p[edges[6], 1], p[edges[6], 2], tube_radius=radius, tube_sides=3,
                                figure=figure)
    geometric_field.mlab.plot3d(p[edges[7], 0], p[edges[7], 1], p[edges[7], 2], tube_radius=radius, tube_sides=3,
                                figure=figure)
    geometric_field.mlab.plot3d(p[edges[8], 0], p[edges[8], 1], p[edges[8], 2], tube_radius=radius, tube_sides=3,
                                figure=figure)
    geometric_field.mlab.plot3d(p[edges[9], 0], p[edges[9], 1], p[edges[9], 2], tube_radius=radius, tube_sides=3,
                                figure=figure)
    geometric_field.mlab.plot3d(p[edges[10], 0], p[edges[10], 1], p[edges[10], 2], tube_radius=radius, tube_sides=3,
                                figure=figure)
    geometric_field.mlab.plot3d(p[edges[11], 0], p[edges[11], 1], p[edges[11], 2], tube_radius=radius, tube_sides=3,
                                figure=figure)

    return figure
