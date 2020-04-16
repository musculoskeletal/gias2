"""
FILE: mesh_fitter.py
LAST MODIFIED: 24-12-2015
DESCRIPTION: Fitting class that combines multiple fitting methods

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""
import logging
import traceback

import scipy

from gias2.fieldwork.field import geometric_field
from gias2.fieldwork.field import geometric_field_fitter as GFF
from gias2.fieldwork.field.tools import fitting_tools
from gias2.learning import PCA_fitting

log = logging.getLogger(__name__)


# fitting class, does rigid, HMF, and normal fits, not specific to pelvis mesh
class log(object):

    def __init__(self, filename):
        self.filename = filename

    def logFitterSettings(self, settingNames, F):
        with open(self.filename, 'w') as f:
            for S in settingNames:
                f.write(str(S) + ': ' + str(getattr(F, S)) + '\n')
            f.write('\n')

    def logFit(self, jobName, errors):
        with open(self.filename, 'a') as f:
            f.write(
                '%(f)9s' % {'f': jobName} + ' || ' + ' | '.join(['%(err)8.6f' % {'err': err} for err in errors]) + '\n')

    def logError(self, jobName, error):
        with open(self.filename, 'a') as f:
            f.write('%(f)9s' % {'f': jobName} + ' || ' + '%(err)s' % {'err': error} + '\n')

        with open(self.filename + '.error.' + str(jobName), 'w') as errlog:
            traceback.print_exc(file=errlog)


class MeshFitter(object):
    """Class for fitting a template mesh to a data cloud via a variety of
    methods. Fitting methods are:
    - rigid-body transformations
    - host-mesh fitting
    - shape-model principal components fitting
    - nodal fitting
    """

    fitSettings = ('templateGFFilename', 'simplemeshDir', 'simplemeshSuffix',
                   'HMFMaxIt', 'HMFHostElemType', 'pcFitMaxIt', 'pcFitNModes',
                   'pcFitmW', 'epD', 'sobW', 'nD', 'nW', 'meshFitMaxIt',
                   'meshFitNClosestPoints', 'meshFitTreeArgs',
                   )

    templateGFFilename = None
    templateEnsFilename = None
    templateMeshFilename = None
    templatePath = None

    gfSaveFileStr = None

    alignObjMode = 'EPDP'
    alignEPD = [10, 10]

    HMFObjMode = 'EPDP'
    HMFMaxIt = 10
    HMFMaxItPerIt = 2
    HMFHostElemType = 'quad444'
    HMFHostElemDisc = [2, 2, 2]
    HMFHostSobD = [4, 4, 4]
    HMFHostSobW = 1e-5
    HMFSlaveEPD = [10, 10]
    HMFSlaveSobD = [5, 5]
    HMFSlaveSobW = scipy.array([1, 1, 1, 1, 2]) * 1e-8
    HMFSlaveNormD = 5
    HMFSlaveNormW = 100.0
    HMFXtol = 1e-6
    HMFTreeArgs = {'distance_upper_bound': 50.0}
    HMFFixedSlaveNodes = None
    HMFCallback = None

    PCFObjMode = 'EPDP'
    PCFEPD = [10, 10]
    PCFInitTrans = None
    PCFInitRot = None
    PCFitMaxIt = 10
    PCFitmW = 2.0
    PCFitNModes = [1, 2, 3, 4]
    PCFitXtol = 1e-6

    meshFitObjMode = 'EPDP'
    meshFitEPD = 5.0
    meshFitSobD = [9, 9]
    meshFitSobW = scipy.array((1e-5, 1e-5, 1e-5, 1e-5, 2e-5)) * 0.1
    meshFitND = 9
    meshFitNW = 20.0
    meshFitXtol = 1e-5
    meshFitFixedNodes = None
    meshFitMaxIt = 5
    meshFitMaxItperIt = 2
    meshFitNClosestPoints = 1
    meshFitTreeArgs = {'distance_upper_bound': 50.0}

    fitter = None

    def __init__(self, jobName):
        self.jobName = jobName
        self.templateGF = None
        self.data = None
        self.dataWeights = None
        self.meshFitError = -1.0
        self.fitterMeshFitError = -1.0
        self.PCFitError = -1.0
        self.HMFError = -1.0

    def setData(self, data, dataWeights=None):
        self.data = data
        self.dataWeights = dataWeights

    def loadTemplateMesh(self):
        self.templateGF = geometric_field.load_geometric_field(
            self.templateGFFilename,
            self.templateEnsFilename,
            self.templateMeshFilename,
            path=self.templatePath
        )

    def setTemplateMesh(self, g):
        self.templateGF = g

    def align(self, initTranslation=None, initRotation=None, initScale=None, sampleData=200):
        """align template mesh to data using translation, rotation and scale
        """
        log.debug('aligning...')

        gfEP = self.templateGF.evaluate_geometric_field(self.alignEPD).T

        if initTranslation is None:
            gfEPCoM = self.templateGF.calc_CoM_2D(self.alignEPD)
            dataCoM = self.data.mean(0)
            initTranslation = dataCoM - gfEPCoM
        if initRotation is None:
            initRotation = [0.0, 0.0, 0.0]
        if initScale is None:
            initScale = 1.0

        rigidScaleT0 = scipy.hstack([initTranslation, initRotation, initScale])
        if self.alignObjMode == 'DPEP':
            rigidScaleTOpt, rigidScaleXOpt = fitting_tools.fitDataRigidScaleDPEP(
                gfEP, self.data, xtol=1e-6,
                maxfev=0, t0=rigidScaleT0,
                sample=sampleData
            )
        elif self.alignObjMode == 'EPDP':
            rigidScaleTOpt, rigidScaleXOpt = fitting_tools.fitDataRigidScaleEPDP(
                gfEP, self.data, xtol=1e-6,
                maxfev=0, t0=rigidScaleT0,
                sample=sampleData
            )

        rigidScalePOpt = fitting_tools.transform3D.transformRigidScale3DAboutCoM(
            self.templateGF.get_field_parameters().squeeze().T,
            rigidScaleTOpt
        )
        self.templateGF.set_field_parameters(rigidScalePOpt.T[:, :, scipy.newaxis])
        self.rigidScaleTOpt = rigidScaleTOpt

    def alignRigid(self, initTranslation=None, initRotation=None, sampleData=200):
        """align template mesh to data using translation, rotation and scale
        """
        log.debug('aligning...')

        gfEP = self.templateGF.evaluate_geometric_field(self.alignEPD).T

        if initTranslation is None:
            gfEPCoM = self.templateGF.calc_CoM_2D(self.alignEPD)
            dataCoM = self.data.mean(0)
            initTranslation = dataCoM - gfEPCoM
        if initRotation is None:
            initRotation = [0.0, 0.0, 0.0]

        rigidT0 = scipy.hstack([initTranslation, initRotation])
        if self.alignObjMode == 'DPEP':
            rigidTOpt, rigidXOpt = fitting_tools.fitDataRigidDPEP(
                gfEP, self.data, xtol=1e-6,
                maxfev=0, t0=rigidT0,
                sample=sampleData
            )
        elif self.alignObjMode == 'EPDP':
            rigidTOpt, rigidXOpt = fitting_tools.fitDataRigidEPDP(
                gfEP, self.data, xtol=1e-6,
                maxfev=0, t0=rigidT0,
                sample=sampleData
            )

        rigidPOpt = fitting_tools.transform3D.transformRigid3DAboutCoM(
            self.templateGF.get_field_parameters().squeeze().T,
            rigidTOpt
        )
        self.templateGF.set_field_parameters(rigidPOpt.T[:, :, scipy.newaxis])
        self.rigidTOpt = rigidTOpt

    def HMF(self):
        """Host mesh fit template mesh to data using a single element host mesh.
        """
        log.debug('host-mesh fitting...')

        slaveGF = self.templateGF
        slaveP0 = slaveGF.get_field_parameters()
        hostGF = GFF.makeHostMesh(slaveP0, 5.0, self.HMFHostElemType)

        # make slave obj
        # squared distance between slaveGF boundary nodes and boundary curve nodes

        if self.HMFObjMode == 'DPEP':
            slaveGObj = GFF.makeObjDPEP(slaveGF, self.data, self.HMFSlaveEPD)
        elif self.HMFObjMode == 'EPDP':
            slaveGObj = GFF.makeObjEPDP(slaveGF, self.data, self.HMFSlaveEPD)
        elif self.HMFObjMode == '2way':
            slaveGObj = GFF.makeObj2Way(slaveGF, self.data, self.HMFSlaveEPD)

        slaveSobObj = GFF.makeSobelovPenalty2D(slaveGF, self.HMFSlaveSobD, self.HMFSlaveSobW)
        slaveNormalSmoother = GFF.normalSmoother2(slaveGF.ensemble_field_function.flatten()[0])
        slaveNormObj = slaveNormalSmoother.makeObj(self.HMFSlaveNormD)

        def slaveObj(x):
            errSurface = slaveGObj(x)
            errSob = slaveSobObj(x)
            errNorm = slaveNormObj(x) * self.HMFSlaveNormW
            return scipy.hstack([errSurface, errSob, errNorm])

        hostParamsOpt, slaveParamsOpt, \
        slaveXi, self.HMFError = fitting_tools.hostMeshFit(
            hostGF, slaveGF, slaveObj,
            maxIt=self.HMFMaxIt,
            sobD=self.HMFHostSobD,
            sobW=self.HMFHostSobW,
            fixedSlaveNodes=self.HMFFixedSlaveNodes
        )
        self.templateGF.set_field_parameters(slaveParamsOpt)

        self.HMFHostGF = hostGF
        self.HMFGFParams = slaveParamsOpt.copy()
        log.debug('HMF fit rms: %(rms)6.4f' % {'rms': self.HMFError})

        return self.HMFError

    def HMFMulti(self):
        """Host mesh fit template mesh to data using a multi-element host mesh.
        """

        log.debug('host-mesh fitting with multiple host elements...')

        slaveGF = self.templateGF
        slaveP0 = slaveGF.get_field_parameters()
        hostGF = GFF.makeHostMeshMulti(slaveP0, 5.0, self.HMFHostElemType, self.HMFHostElemDisc)

        # make slave obj
        # squared distance between slaveGF boundary nodes and boundary curve nodes

        if self.HMFObjMode == 'DPEP':
            slaveGObj = GFF.makeObjDPEP(slaveGF, self.data, self.HMFSlaveEPD)
        elif self.HMFObjMode == 'EPDP':
            slaveGObj = GFF.makeObjEPDP(slaveGF, self.data, self.HMFSlaveEPD)
        elif self.HMFObjMode == '2way':
            slaveGObj = GFF.makeObj2Way(slaveGF, self.data, self.HMFSlaveEPD)

        slaveSobObj = GFF.makeSobelovPenalty2D(slaveGF, self.HMFSlaveSobD, self.HMFSlaveSobW)
        slaveNormalSmoother = GFF.normalSmoother2(slaveGF.ensemble_field_function.flatten()[0])
        slaveNormObj = slaveNormalSmoother.makeObj(self.HMFSlaveNormD)

        def slaveObj(x):
            errSurface = slaveGObj(x)
            errSob = slaveSobObj(x)
            errNorm = slaveNormObj(x) * self.HMFSlaveNormW
            return scipy.hstack([errSurface, errSob, errNorm])

        hostParamsOpt, slaveParamsOpt, \
        slaveXi, HMFError = fitting_tools.hostMeshFitMulti(
            hostGF, slaveGF, slaveObj,
            maxIt=self.HMFMaxIt,
            sobD=self.HMFHostSobD,
            sobW=self.HMFHostSobW,
            fixedSlaveNodes=self.HMFFixedSlaveNodes
        )
        self.templateGF.set_field_parameters(slaveParamsOpt)

        self.HMFHostGF = hostGF
        self.HMFGFParams = slaveParamsOpt.copy()
        self.HMFSlaveXi = slaveXi
        self.HMFError = HMFError
        log.debug('HMF fit rms: %(rms)6.4f' % {'rms': self.HMFError})

        return self.HMFError

    def HMFMultiPerItSearch(self):
        """Host mesh fit template mesh to data using a multi-element host mesh
        and closest-point updates per n iterations where n is defined by
        self.HMFMaxItPerIt.
        """
        log.debug('host-mesh fitting with multiple host elements and closes point search per n iterations...')

        slaveGF = self.templateGF
        slaveP0 = slaveGF.get_field_parameters()
        hostGF = GFF.makeHostMeshMulti(slaveP0, 5.0, self.HMFHostElemType, self.HMFHostElemDisc)

        fitOutput = fitting_tools.hostMeshFitMultiPerItSearch(
            self.data, hostGF, slaveGF, self.HMFObjMode,
            self.HMFSlaveEPD,
            self.HMFSlaveSobD, self.HMFSlaveSobW,
            self.HMFSlaveNormD, self.HMFSlaveNormW,
            hostSobD=self.HMFHostSobD, hostSobW=self.HMFHostSobW,
            dataWeights=self.dataWeights, slaveXi=None,
            xtol=self.HMFXtol, maxIt=self.HMFMaxIt, maxItPerIt=self.HMFMaxItPerIt,
            fixedSlaveNodes=self.HMFFixedSlaveNodes,
            treeArgs=self.HMFTreeArgs, fitOutputCallback=self.HMFCallback,
            verbose=True)

        hostParamsOpt, slaveParamsOpt, slaveXi, HMFError = fitOutput
        self.templateGF.set_field_parameters(slaveParamsOpt)

        self.HMFHostGF = hostGF
        self.HMFGFParams = slaveParamsOpt.copy()
        self.HMFSlaveXi = slaveXi
        self.HMFError = HMFError
        log.debug('HMF fit rms: %(rms)6.4f' % {'rms': self.HMFError})

        return self.HMFError

    def pcFit(self, pc, mode0Offset=0.0):
        """Fit template mesh to data by deformation along principal components
        """

        self.templateGF.set_field_parameters(pc.getMean().reshape((3, -1, 1)))
        if self.PCFInitTrans is None:
            gfCoM = self.templateGF.calc_CoM_2D(self.PCFEPD)
            dataCoM = self.data.mean(0)
            self.PCFInitTrans = dataCoM - gfCoM
        if self.PCFInitRot is None:
            self.PCFInitRot = [0.0, 0.0, 0.0]

        rigidMode0T0 = scipy.hstack([self.PCFInitTrans, self.PCFInitRot])
        # ~ rigidMode0T0 = scipy.hstack( [[0.0,0.0,0.0], [0.0,0.0,0.0]] )

        if self.PCFObjMode == 'DPEP':
            gObj = GFF.makeObjDPEP(self.templateGF, self.data, self.PCFEPD)
        elif self.PCFObjMode == 'EPDP':
            gObj = GFF.makeObjEPDP(self.templateGF, self.data, self.PCFEPD)
        elif self.PCFObjMode == '2way':
            gObj = GFF.makeObj2Way(self.templateGF, self.data, self.PCFEPD)

        # ~ gObj = GFF.makeObjEPDP( self.templateGF, self.data, self.HMFSlaveEPD, dataWeights=None ) # same as region fit
        # ~ gObj = GFF.makeObjDPEP( self.templateGF, self.simplemesh.v, self.epD, dataWeights=None )

        pcFitter = PCA_fitting.PCFit(pc=pc)
        pcFitter.xtol = self.PCFitXtol

        log.debug('\npc fitting...')
        # rigid fit template gf
        rigidOpt, rigidPOpt = pcFitter.rigidFit(gObj, x0=rigidMode0T0)
        self.templateGF.set_field_parameters(rigidPOpt.reshape((3, -1, 1)))

        # rigid+mode0 fit template gf
        rigidMode0X0 = scipy.hstack([rigidOpt, mode0Offset])
        rigidMode0Opt, rigidMode0POpt = pcFitter.rigidMode0Fit(
            gObj, x0=rigidMode0X0,
            mWeight=self.PCFitmW
        )
        self.templateGF.set_field_parameters(rigidMode0POpt.reshape((3, -1, 1)))

        # rigid+modeN fit template gf
        rigidModeNOpt, rigidModeNPOpt = pcFitter.rigidModeNFit(
            gObj, modes=self.PCFitNModes,
            x0=None, mWeight=self.PCFitmW,
            maxfev=0
        )
        self.templateGF.set_field_parameters(rigidModeNPOpt.reshape((3, -1, 1)))
        self.PCFitError = scipy.sqrt(gObj(rigidModeNPOpt).mean())  # obj returns squared error

        self.PCGFParams = rigidModeNPOpt.reshape((3, -1, 1)).copy()
        log.debug('pc fit rms: %(pcrms)6.4f' % {'pcrms': self.PCFitError})
        log.debug('mode weights: ' + ' '.join(['%(0)5.3f' % {'0': i} for i in rigidModeNOpt[6:]]))

        return self.PCFitError

    def meshFit(self):
        """Fit template mesh to data by optimisation of nodal parameters.
        """
        log.debug('fitting...')
        self.templateGF, gfFitPOpt, \
        self.meshFitError = fitting_tools.fitSurfacePerItSearch(
            self.meshFitObjMode, self.templateGF,
            self.data, self.meshFitEPD,
            self.meshFitSobD, self.meshFitSobW,
            self.meshFitND, self.meshFitNW,
            fixedNodes=self.meshFitFixedNodes,
            xtol=self.meshFitXtol,
            itMax=self.meshFitMaxIt,
            itMaxPerIt=self.meshFitMaxItperIt,
            nClosestPoints=self.meshFitNClosestPoints,
            treeArgs=self.meshFitTreeArgs
        )

        self.meshFitGFParams = self.templateGF.get_field_parameters()
        log.debug('mesh fit rms: %(rms)6.4f' % {'rms': self.meshFitError})
        return self.meshFitError

    def fitterMeshFit(self, drms=0.0, output=True):
        self.fitter, self.fitterMeshFitError = fitting_tools.fitterFit(self.templateGF,
                                                                       self.epD,
                                                                       self.data,
                                                                       fit=self.fitter,
                                                                       maxIt=self.meshFitMaxIt,
                                                                       drms=drms,
                                                                       output=output,
                                                                       doFit=True
                                                                       )

        # ~ self.meshFitError = 0.0
        self.fitterMeshFitGFParams = self.templateGF.get_field_parameters()
        log.debug('mesh fit rms: %(rms)6.4f' % {'rms': self.fitterMeshFitError})
        return self.fitterMeshFitError

    def logFitErrors(self, logger):

        # ~ if self.pcFitError==-1.0:
        # ~ err1 = self.HMFError
        # ~ else:
        # ~ err1 = self.pcFitError
        errors = [self.PCFitError, self.HMFError, self.meshFitError]
        logger.logFit(self.jobName, errors)

    def saveGF(self, filename=None):
        if filename is None:
            filename = self.gfSaveFileStr % {'jobName': self.jobName}

        log.debug('...saving GF to... ' + filename)
        self.templateGF.save_geometric_field(filename)
