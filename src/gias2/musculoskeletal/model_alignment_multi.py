"""
FILE: model_alignment_multi.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: functions for alignining multiple fieldwork models, e.g. joints

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import copy

import numpy

from gias2.common import transform3D
from gias2.musculoskeletal import fw_femur_measurements
from gias2.musculoskeletal import fw_femur_model_data as fmd
from gias2.musculoskeletal import fw_model_landmarks as fml
from gias2.musculoskeletal import fw_pelvis_measurements
from gias2.musculoskeletal import model_alignment


def normaliseVector(v):
    return v / numpy.linalg.norm(v)


def alignTibiaFibulaMeshParametersAnatomicSingle(tib, fib):
    """ given a list of femur geometric fields, align them geometrically.
    returns the aligned field parameters
    """
    tibNodes = tib.get_all_point_positions()
    fibNodes = fib.get_all_point_positions()
    LM = fibNodes[fml._fibulaLMNode]
    MM = tibNodes[fml._tibiaMMNode]
    LC = tibNodes[fml._tibiaLCNode]
    MC = tibNodes[fml._tibiaMCNode]

    alignedTibParams, T = model_alignment.alignAnatomicTibiaFibulaGroodSuntay(tibNodes, MM, LM, MC, LC, True)
    alignedFibParams, T = model_alignment.alignAnatomicTibiaFibulaGroodSuntay(fibNodes, MM, LM, MC, LC, True)
    alignedTibParams = alignedTibParams.T[:, :, numpy.newaxis]
    alignedFibParams = alignedFibParams.T[:, :, numpy.newaxis]
    return alignedTibParams, alignedFibParams, T


def alignTibiaFibulaMeshParametersAnatomic(tibs, fibs):
    """ given a list of femur geometric fields, align them geometrically.
    returns the aligned field parameters
    """
    alignedTibParamsAll = []
    alignedFibParamsAll = []
    for tib, fib in zip(tibs, fibs):
        alignedTibParams, alignedFibParams, T = alignTibiaFibulaMeshParametersAnatomicSingle(tib, fib)
        alignedTibParamsAll.append(alignedTibParams)
        alignedFibParamsAll.append(alignedFibParams)

    return alignedTibParamsAll, alignedFibParamsAll


def calcAngle(v1, v2):
    return numpy.arccos(numpy.dot(v1, v2) / (numpy.linalg.norm(v1) * numpy.linalg.norm(v2)))


def alignPelvisRightFemurAnatomic(pelvisG, rightFemurG):
    """
    aligns meshes of the pelvis and femur to the pelvis
    anatomic coordinate system, and aligns the femur z 
    to the pelvis y axis (inf-sup) by rotation about the hip joint 
    centre. The HJC is the midpoint between the femoral head 
    centre and the acetabulum centre as estimated by 
    sphere fitting.

    The input mesh is expected to have 4 submeshes - 
    RH, LH, sac, right femur
    """

    # first transform meshes to pelvis anatomic CS
    pelvisG = copy.deepcopy(pelvisG)
    rightFemurG = copy.deepcopy(rightFemurG)

    pelvisAlignedParams, pelvisAnatAlignT = model_alignment.alignWholePelvisMeshParametersAnatomicSingle(pelvisG)
    pelvisG.set_field_parameters(pelvisAlignedParams)

    rightFemurG.transformAffine(pelvisAnatAlignT)

    # calculate femoral head and acetabulum centres 
    pelvisM = fw_pelvis_measurements.PelvisMeasurements(pelvisG)
    pelvisM.calcAcetabulumDiameters()
    rAcetab = pelvisM.measurements['right_acetabulum_diameter']

    femurM = fw_femur_measurements.FemurMeasurements(rightFemurG)
    femurM.calcHeadDiameter()
    rFH = femurM.measurements['head_diameter']

    # average to get joint rotation centre
    HJC = (rFH.centre + rAcetab.centre) / 2.0

    # estimate joint spacing
    jointSpacing = rAcetab.value - rFH.value

    """
    # align femur z's X-Z projection to the Z axis - rotate about Y
    rFz = _calcFemurZ(rightFemurG)
    thetaY = calcAngle( numpy.array([rFz[2], rFz[0]]), numpy.array([1.0,0.0]))
    rightFemurG.transformRotateAboutP( numpy.array([0.0,thetaY,0.0]), HJC )

    # calculate angle between femur z's global Y-Z projection and the Y axis - rotate about X
    rFz = _calcFemurZ(rightFemurG)
    thetaX = calcAngle( numpy.array([rFz[2], rFz[1]]), numpy.array([0.0,1.0]))
    rightFemurG.transformRotateAboutP( numpy.array([thetaX,0.0,0.0]), HJC )
    """

    # align femur anatomic coord system to global
    d = (10, 10)
    head = rightFemurG.calc_CoM_2D(d, elem=fmd.assemblyElementsNumbers['head'])
    lc = rightFemurG.calc_CoM_2D(d, elem=fmd.assemblyElementsNumbers['lateralcondyle'])
    mc = rightFemurG.calc_CoM_2D(d, elem=fmd.assemblyElementsNumbers['medialcondyle'])
    oF = (mc + lc) / 2.0
    zF = normaliseVector(head - oF)
    yF = normaliseVector(numpy.cross(zF, (lc - oF)))
    xF = normaliseVector(numpy.cross(yF, zF))

    oF = rFH.centre

    u = numpy.array([oF, oF - yF, oF - zF, oF + xF])
    ut = numpy.array([oF, oF + [1.0, 0, 0], oF + [0, 1.0, 0], oF + [0, 0, 1.0]])
    femurAlignT = transform3D.directAffine(u, ut)
    rightFemurG.transformAffine(femurAlignT)

    # check new joint spacing
    # pelvisM2 = fw_pelvis_measurements.PelvisMeasurements(pelvisG)
    # pelvisM2.calcAcetabulumDiameters()
    # rAcetab2 = pelvisM2.measurements['right_acetabulum_diameter']

    femurM2 = fw_femur_measurements.FemurMeasurements(rightFemurG)
    femurM2.calcHeadDiameter()
    rFH2 = femurM2.measurements['head_diameter']
    jointSpacing2 = rAcetab.value - rFH2.value
    HJC2 = (rFH2.centre + rAcetab.centre) / 2.0

    # translate femur to correct for changes in joint spacing
    tCorrect = (rAcetab.centre - rFH2.centre) * (jointSpacing2 - jointSpacing)
    rightFemurG.transformTranslate(tCorrect)

    return pelvisG, rightFemurG


def _calcFemurZ(femurG):
    d = [10, 10]
    head = femurG.calc_CoM_2D(d, elem=fmd.assemblyElementsNumbers['head'])
    lc = femurG.calc_CoM_2D(d, elem=fmd.assemblyElementsNumbers['lateralcondyle'])
    mc = femurG.calc_CoM_2D(d, elem=fmd.assemblyElementsNumbers['medialcondyle'])
    rFOrigin = (mc + lc) / 2.0
    rFz = normaliseVector(head - rFOrigin)
    if rFz[2] < 0.0:
        rFz *= -1.0

    return rFz


def alignAnatomicTibiaFibulaGroodSuntay(X, MM, LM, MC, LC, returnT=False):
    IC = (MC + LC) / 2.0
    IM = (MM + LM) / 2.0

    z = normaliseVector(IC - IM)
    y = normaliseVector(numpy.cross(z, MC - LC))
    x = normaliseVector(numpy.cross(y, z))

    u = numpy.array([IC, IC + x, IC + y, IC + z])
    ut = numpy.array([[0, 0, 0],
                      [1, 0, 0],
                      [0, 1, 0],
                      [0, 0, 1]])

    t = transform3D.directAffine(u, ut)
    if returnT:
        return transform3D.transformAffine(X, t), t
    else:
        return transform3D.transformAffine(X, t)
