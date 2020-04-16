"""
FILE: fw_pelvis_measurements.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: functions and classes for taking anthropometric measurements from
femur geometric meshes

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import copy
import logging
import pickle

import numpy

from gias2.common import geoprimitives as FT
from gias2.musculoskeletal import fw_model_landmarks as fml
from gias2.musculoskeletal import fw_pelvis_model_data as pmd
from gias2.musculoskeletal import model_alignment
from gias2.musculoskeletal import pelvis_hjc_estimation as HJC

log = logging.getLogger(__name__)


class Measurement(object):

    def __init__(self, name, value=None):
        self.name = name
        self.value = value


class PelvisMeasurements(object):
    regionName2ElementMap = {'LH': 1, 'RH': 0, 'sac': 2}
    leftAcetabulumElems = ((regionName2ElementMap['LH'], pmd.hemiPelvisAcetabulumElements),)
    rightAcetabulumElems = ((regionName2ElementMap['RH'], pmd.hemiPelvisAcetabulumElements),)
    landmarks = ('LASIS', 'RASIS', 'LPSIS', 'RPSIS', 'Sacral', 'LPS', 'RPS', 'LIS', 'RIS',
                 'LIT', 'RIT', 'LHJC', 'RHJC',
                 )
    acs = 'isb'  # or app
    epD = [10, 10]

    def __init__(self, GF=None, acs='isb'):

        self.acs = acs
        self.measurements = {
            'left_acetabulum_diameter': None,
            'right_acetabulum_diameter': None,
            'left_HJC_Bell': None,
            'left_HJC_Tylkowski': None,
            'left_HJC_Andriacchi': None,
            'left_HJC_mesh': None,
            'right_HJC_Bell': None,
            'right_HJC_Tylkowski': None,
            'right_HJC_Andriacchi': None,
            'right_HJC_mesh': None,
            'landmarks_unaligned': None,
            'landmarks_ACS': None,
            'max_height': None,
            'max_width': None,
            'ilial_depth': None,
            'inter_ASIS_distance': None,
            'inter_PSIS_distance': None,
            'inter_PS_distance': None,
        }

        if GF != None:
            self.GF = copy.deepcopy(GF)
        else:
            self.GF = None

        self.GFACS = None  # anatomic CS
        self.EP = None
        self.EPACS = None
        self.EPMap = None

        if self.GF != None:
            self._alignToAnatomicCS()
            self._getLandmarks()

    def saveMeasurements(self, filename):
        with open(filename, 'w') as f:
            pickle.dump(self.measurements, f, protocol=2)

    def loadMeasurements(self, filename):
        with open(filename, 'r') as f:
            self.measurements = pickle.load(f)

    def calcMeasurements(self):
        self.calcAcetabulumDiameters()
        # self.calcACSAcetabulum()
        self.calcInterHJCDistance()
        self.calcInterASISDistance()
        self.calcInterPSISDistance()
        self.calcInterPSDistance()
        self.calcIlialDepth()
        self.calcMaxHeight()
        self.calcMaxWidth()

    def printMeasurements(self):
        m = list(self.measurements.keys())
        m.sort()
        for mi in m:
            if self.measurements[mi] is not None:
                log.debug(mi, ':', self.measurements[mi].value)

    def _alignToAnatomicCS(self):

        if self.acs == 'isb':
            alignedParams = model_alignment.alignWholePelvisMeshParametersAnatomic([self.GF, ])[0]
        elif self.acs == 'app':
            alignedParams = model_alignment.alignWholePelvisMeshParametersAnatomicAPP([self.GF, ])[0]
        else:
            raise ValueError('Unknown anatomic coordinate system {}'.format(self.acs))

        self.GFACS = copy.deepcopy(self.GF)
        self.GFACS.set_field_parameters(alignedParams)

    def _getLandmarks(self):

        # unaligned landmarks
        landmarksUnaligned = {}
        for ln in self.landmarks:
            e = fml.makeLandmarkEvaluator('pelvis-' + ln, self.GF)
            landmarksUnaligned[ln] = e(self.GF.field_parameters)

        self.measurements['landmarks_unaligned'] = Measurement(
            'landmarks_unaligned',
            landmarksUnaligned
        )

        # aligned landmarks
        landmarksACS = {}
        for ln in self.landmarks:
            e = fml.makeLandmarkEvaluator('pelvis-' + ln, self.GFACS)
            landmarksACS[ln] = e(self.GFACS.field_parameters)

        self.measurements['landmarks_ACS'] = Measurement(
            'landmarks_ACS',
            landmarksACS
        )

        # nodesUnaligned = self.GF.get_all_point_positions()
        # landmarksUnaligned = {}
        # for landmarkName, node in pmd.landmarksNodes.items():
        #   landmarksUnaligned[landmarkName] = nodesUnaligned[node]

        # landmarksUnaligned['PS'] = 0.5*(landmarksUnaligned['LPS']+landmarksUnaligned['RPS'])

        # self.measurements['landmarks_unaligned'] = Measurement('landmarks_unaligned',\
        #                                                      landmarksUnaligned
        #                                                      )

        # nodesACS = self.GFACS.get_all_point_positions()
        # landmarksACS = {}
        # for landmarkName, node in pmd.landmarksNodes.items():
        #   landmarksACS[landmarkName] = nodesACS[node]

        # landmarksACS['PS'] = 0.5*(landmarksACS['LPS']+landmarksACS['RPS'])

        # self.measurements['landmarks_ACS'] = Measurement('landmarks_ACS',\
        #                                                landmarksACS )

    def _evaluateGF(self):
        self.EPACS = self.GFACS.evaluate_geometric_field(self.epD).T
        self.EP = self.GF.evaluate_geometric_field(self.epD).T
        self.EPMap = self.GF.getElementPointI(
            self.epD,
            'all'
        )
        # self.EPMap = self.GF.getElementPointINested(
        #               self.epD,
        #               numpy.sort(self.regionName2ElementMap.values())
        #               )

    def calcAcetabulumDiameters(self):

        lhjcEval = fml.makeLandmarkEvaluator('pelvis-LHJC', self.GF, radius=True)
        lhjcCenter, lhjcRadius = lhjcEval(self.GF.field_parameters)
        LAm = Measurement('left_acetabulum_diameter', lhjcRadius * 2.0)
        LAm.centre = numpy.array([lhjcCenter])
        self.measurements['left_acetabulum_diameter'] = LAm

        rhjcEval = fml.makeLandmarkEvaluator('pelvis-RHJC', self.GF, radius=True)
        rhjcCenter, rhjcRadius = rhjcEval(self.GF.field_parameters)
        RAm = Measurement('right_acetabulum_diameter', rhjcRadius * 2.0)
        RAm.centre = numpy.array([rhjcCenter])
        self.measurements['right_acetabulum_diameter'] = RAm

        #######
        # OLD #
        #######
        # LHGF = self.GF.makeGFFromElements('LH', [self.regionName2ElementMap['LH'],],\
        #                                 pmd.pelvisCubicBasisTypes )
        # LHGF.flatten_ensemble_field_function()
        # leftAcetabEPs = LHGF.evaluate_geometric_field_in_elements(self.epD,\
        #               self.leftAcetabulumElems[0][1]).T

        # # LArms, [LAcx,LAcy,LAcz,LAr] = FT.fitSphere( leftAcetabEPs )
        # (LAcx, LAcy, LAcz), LAr = FT.fitSphereAnalytic( leftAcetabEPs )

        # LAm = Measurement( 'left_acetabulum_diameter', LAr*2.0 )
        # LAm.centre = numpy.array([LAcx, LAcy, LAcz])
        # self.measurements['left_acetabulum_diameter'] = LAm

        # RHGF = self.GF.makeGFFromElements('RH', [self.regionName2ElementMap['RH'],],\
        #                                 pmd.pelvisCubicBasisTypes )
        # RHGF.flatten_ensemble_field_function()
        # rightAcetabEPs = RHGF.evaluate_geometric_field_in_elements(self.epD,\
        #                self.rightAcetabulumElems[0][1]).T

        # # RArms, [RAcx,RAcy,RAcz,RAr] = FT.fitSphere( rightAcetabEPs )
        # (RAcx, RAcy, RAcz), RAr = FT.fitSphereAnalytic( rightAcetabEPs )

        # RAm = Measurement( 'right_acetabulum_diameter', RAr*2.0 )
        # RAm.centre = numpy.array([RAcx, RAcy, RAcz])
        # self.measurements['right_acetabulum_diameter'] = RAm

        # # add HJC to landmarks
        # self.measurements['landmarks_unaligned'].value['LHJC'] = LAm.centre.copy()
        # self.measurements['landmarks_unaligned'].value['RHJC'] = RAm.centre.copy()

        return LAm.value, RAm.value

    def calcACSAcetabulum(self):
        """ Calculate acetabular anteversion and abduction angles using the ACS
        aligned model.

        Needs to update to use flat pelvis meshes.
        """

        if self.acs == 'isb':
            _calcAnteversionACS = _calcAnteversionISB
            _calcAbductionACS = _calcAbductionISB
        elif self.acs == 'app':
            _calcAnteversionACS = _calcAnteversionAPP
            _calcAbductionACS = _calcAbductionAPP
        else:
            raise ValueError('Unknown ACS {}'.format(self.acs))

        # left
        # lhjcEval = fml.makeLandmarkEvaluator('pelvis-LHJC', self.GFACS)
        # lhjcCenter = lhjcEval(self.GFACS.field_parameters)
        # LHJC = numpy.array([lhjcCenter])

        self.LHGF = self.GFACS.makeGFFromElements('LH',
                                                  [self.regionName2ElementMap['LH'], ],
                                                  pmd.pelvisCubicBasisTypes)
        self.LHGF.flatten_ensemble_field_function()
        leftAcetabEPs = self.LHGF.evaluate_geometric_field_in_elements(self.epD,
                                                                       self.leftAcetabulumElems[0][1]).T
        (LAcx, LAcy, LAcz), LAr = FT.fitSphereAnalytic(leftAcetabEPs)
        LHJC = numpy.array([LAcx, LAcy, LAcz])
        self.measurements['landmarks_ACS'].value['LHJC'] = LHJC

        self.LHJCPlane = self._calcAcetabulumPlaneACS(
            LHJC,
            self.LHGF.get_all_point_positions()[pmd.hemiPelvisAcetabularCupRimNodes, :],
            leftAcetabEPs
        )
        self.measurements['landmarks_ACS'].value['left_acetabulum_plane'] = self.LHJCPlane.N
        self.measurements['left_anteversion'] = Measurement(
            'left_anteversion',
            _calcAnteversionACS(
                self.LHJCPlane.N,
                'left'
            )
        )
        self.measurements['left_abduction'] = Measurement(
            'left_abduction',
            _calcAbductionACS(
                self.LHJCPlane.N,
                'left'
            )
        )

        # right
        self.RHGF = self.GFACS.makeGFFromElements('RH', [self.regionName2ElementMap['RH'], ],
                                                  pmd.pelvisCubicBasisTypes)
        self.RHGF.flatten_ensemble_field_function()
        rightAcetabEPs = self.RHGF.evaluate_geometric_field_in_elements(self.epD,
                                                                        self.rightAcetabulumElems[0][1]).T
        (RAcx, RAcy, RAcz), RAr = FT.fitSphereAnalytic(rightAcetabEPs)
        RHJC = numpy.array([RAcx, RAcy, RAcz])
        self.measurements['landmarks_ACS'].value['RHJC'] = LHJC

        self.RHJCPlane = self._calcAcetabulumPlaneACS(
            RHJC,
            self.RHGF.get_all_point_positions()[pmd.hemiPelvisAcetabularCupRimNodes, :],
            rightAcetabEPs
        )
        self.measurements['landmarks_ACS'].value['right_acetabulum_plane'] = self.RHJCPlane.N
        self.measurements['right_anteversion'] = Measurement(
            'right_anteversion',
            _calcAnteversionACS(
                self.RHJCPlane.N,
                'right'
            )
        )
        self.measurements['right_abduction'] = Measurement(
            'right_abduction',
            _calcAbductionACS(
                self.RHJCPlane.N,
                'right'
            )
        )

        return LHJC, RHJC

    def _calcACSHJCMesh(self):
        """
        Calculate HJC from sphere fit to mesh
        """
        lhjcEval = fml.makeLandmarkEvaluator('pelvis-LHJC', self.GFACS, radius=True)
        lhjcCenter, lhjcRadius = lhjcEval(self.GFACS.field_parameters)

        rhjcEval = fml.makeLandmarkEvaluator('pelvis-RHJC', self.GFACS, radius=True)
        rhjcCenter, rhjcRadius = rhjcEval(self.GFACS.field_parameters)

        return lhjcCenter, rhjcCenter

    def calcHJCPredictions(self, popClass):

        L = self.measurements['landmarks_ACS']
        LHJC_T, RHJC_T, ASIS2ASIS = HJC.HJCTylkowski(L.value['LASIS'],
                                                     L.value['RASIS'],
                                                     popClass)
        LHJC_A, RHJC_A, LO, RO = HJC.HJCAndriacchi(L.value['LASIS'],
                                                   L.value['RASIS'],
                                                   0.5 * (L.value['LPS'] + L.value['RPS']),
                                                   popClass)
        LHJC_B, RHJC_B, ASIS2ASIS, LO, RO = HJC.HJCBell(L.value['LASIS'],
                                                        L.value['RASIS'],
                                                        0.5 * (L.value['LPS'] + L.value['RPS']),
                                                        popClass)
        LHJC_D, RHJC_D, ASIS2ASIS, H, D = HJC.HJCSeidel(L.value['LASIS'],
                                                        L.value['RASIS'],
                                                        L.value['LPSIS'],
                                                        L.value['RPSIS'],
                                                        0.5 * (L.value['LPS'] + L.value['RPS']),
                                                        popClass)
        LHJC_H, RHJC_H, PW, PD = HJC.HJCHarrington(L.value['LASIS'],
                                                   L.value['RASIS'],
                                                   L.value['LPSIS'],
                                                   L.value['RPSIS'],
                                                   popClass)
        LHJC_mesh, RHJC_mesh = self._calcACSHJCMesh()

        self.measurements['left_HJC_Bell'] = Measurement('left_HJC_Bell', LHJC_B)
        self.measurements['left_HJC_Tylkowski'] = Measurement('left_HJC_Tylkowski', LHJC_T)
        self.measurements['left_HJC_Andriacchi'] = Measurement('left_HJC_Andriacchi', LHJC_A)
        self.measurements['left_HJC_Seidel'] = Measurement('left_HJC_Seidel', LHJC_D)
        self.measurements['left_HJC_Harrington'] = Measurement('left_HJC_Harrington', LHJC_H)
        self.measurements['left_HJC_mesh'] = Measurement('left_HJC_mesh', LHJC_mesh)

        self.measurements['right_HJC_Bell'] = Measurement('right_HJC_Bell', RHJC_B)
        self.measurements['right_HJC_Tylkowski'] = Measurement('right_HJC_Tylkowski', RHJC_T)
        self.measurements['right_HJC_Andriacchi'] = Measurement('right_HJC_Andriacchi', RHJC_A)
        self.measurements['right_HJC_Seidel'] = Measurement('right_HJC_Seidel', RHJC_D)
        self.measurements['right_HJC_Harrington'] = Measurement('right_HJC_Harrington', RHJC_H)
        self.measurements['right_HJC_mesh'] = Measurement('right_HJC_mesh', RHJC_mesh)

        self.measurements['_ASIS2ASIS'] = Measurement('_ASIS2ASIS', ASIS2ASIS)
        self.measurements['_LO'] = Measurement('_LO', LO)
        self.measurements['_RO'] = Measurement('_RO', RO)
        self.measurements['_H'] = Measurement('_H', H)
        self.measurements['_D'] = Measurement('_D', D)

    def calcInterHJCDistance(self):
        """ Distance between LHJC and RHJC
        """
        x1 = self.measurements['landmarks_unaligned'].value['LHJC']
        x2 = self.measurements['landmarks_unaligned'].value['RHJC']
        d = numpy.sqrt(((x1 - x2) ** 2.0).sum())
        self.measurements['inter_HJC_distance'] = Measurement('inter_HJC_distance', d)

    def calcInterASISDistance(self):
        """ Distance between LASIS and RASIS
        """
        x1 = self.measurements['landmarks_unaligned'].value['LASIS']
        x2 = self.measurements['landmarks_unaligned'].value['RASIS']
        d = numpy.sqrt(((x1 - x2) ** 2.0).sum())
        self.measurements['inter_ASIS_distance'] = Measurement('inter_ASIS_distance', d)

    def calcInterPSISDistance(self):
        """ Distance between LPSIS and RPSIS
        """
        x1 = self.measurements['landmarks_unaligned'].value['LPSIS']
        x2 = self.measurements['landmarks_unaligned'].value['RPSIS']
        d = numpy.sqrt(((x1 - x2) ** 2.0).sum())
        self.measurements['inter_PSIS_distance'] = Measurement('inter_PSIS_distance', d)

    def calcInterPSDistance(self):
        """ Distance between LPSIS and RPSIS
        """
        x1 = self.measurements['landmarks_unaligned'].value['LPS']
        x2 = self.measurements['landmarks_unaligned'].value['RPS']
        d = numpy.sqrt(((x1 - x2) ** 2.0).sum())
        self.measurements['inter_PS_distance'] = Measurement('inter_PS_distance', d)

    def calcIlialDepth(self):
        """ Distance between ASIS and PSIS
        """
        x1l = self.measurements['landmarks_unaligned'].value['LASIS']
        x2l = self.measurements['landmarks_unaligned'].value['LPSIS']
        dl = numpy.sqrt(((x1l - x2l) ** 2.0).sum())
        self.measurements['left_ilial_depth'] = Measurement('left_ilial_depth', dl)

        x1r = self.measurements['landmarks_unaligned'].value['RASIS']
        x2r = self.measurements['landmarks_unaligned'].value['RPSIS']
        dr = numpy.sqrt(((x1r - x2r) ** 2.0).sum())
        self.measurements['right_ilial_depth'] = Measurement('right_ilial_depth', dr)

        d = 0.5 * (dl + dr)
        self.measurements['ilial_depth'] = Measurement('ilial_depth', d)

    def calcMaxHeight(self):
        """ Distance between highest and lowest points in anatomic
        coordinate system
        """
        if self.EPACS is None:
            self._evaluateGF()

        if self.acs == 'app':
            d = self.EPACS[:, 2].max() - self.EPACS[:, 2].min()
        elif self.acs == 'isb':
            d = self.EPACS[:, 1].max() - self.EPACS[:, 1].min()

        self.measurements['max_height'] = Measurement('max_height', d)

    def calcMaxWidth(self):
        """ Distance between leftest and rightest points in anatomic
        coordinate system
        """
        if self.EPACS is None:
            self._evaluateGF()

        if self.acs == 'app':
            d = self.EPACS[:, 0].max() - self.EPACS[:, 0].min()
        elif self.acs == 'isb':
            d = self.EPACS[:, 2].max() - self.EPACS[:, 2].min()

        self.measurements['max_width'] = Measurement('max_width', d)

    def _calcAcetabulumPlaneACS(self, hjc, cupRimPoints, cupPoints):
        # calculate orientation based on plane of acetabular rim
        cupRimPlane = FT.fitPlaneLS(cupRimPoints)
        hjcProj = cupRimPlane.project2Plane3D(hjc)
        cupRimPlane.O = hjcProj

        # plane normal should point towards bone
        offsetPlus = cupRimPlane.O + 10 * cupRimPlane.N
        offsetMinus = cupRimPlane.O - 10 * cupRimPlane.N
        offsetPlusMaxD = ((cupPoints - offsetPlus) ** 2.0).sum(1).max()
        offsetMinusMaxD = ((cupPoints - offsetMinus) ** 2.0).sum(1).max()
        if offsetPlusMaxD > offsetMinusMaxD:
            cupRimPlane.N *= -1.0

        return cupRimPlane


def _calcAnteversionISB(N, side):
    # calc anteversion angle given cup rim plane normal
    # anteversion is measured in the transverse (axial) plane, angle btw left-right
    # axis and transverse plane projection (dims 0,1) of cup plane normal
    a = -numpy.arctan(N[0] / N[2]) * 180 / numpy.pi
    if side == 'left':
        a *= -1
    return a


def _calcAbductionISB(N, side):
    # calc abduction angle given cup rim plane normal.
    # abduction is measured in the coronal plane, angle btw sup-inf axis and
    # coronal plane projection (dims 0,2) of cup plane normal
    a = numpy.arctan(N[2] / N[1]) * 180 / numpy.pi
    if side == 'left':
        a *= -1
    return a


def _calcAnteversionAPP(N, side):
    # calc anteversion angle given cup rim plane normal
    # anteversion is measured in the transverse (axial) plane, angle btw left-right
    # axis and transverse plane projection (dims 0,1) of cup plane normal
    a = numpy.arctan(N[1] / N[0]) * 180 / numpy.pi
    if side == 'left':
        a *= -1
    return a


def _calcAbductionAPP(N, side):
    # calc abduction angle given cup rim plane normal.
    # abduction is measured in the coronal plane, angle btw sup-inf axis and
    # coronal plane projection (dims 0,2) of cup plane normal
    a = -numpy.arctan(N[0] / N[2]) * 180 / numpy.pi
    if side == 'left':
        a *= -1
    return a
