"""
FILE: fw_model_landmark.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: Functions for creating evaluator functions for anatomic landmarks
on fieldwork models

You can get a list of all landmarks by

fw_model_landmarks.landmarkNames()

To evaluate a landmark (e.g. femoral head centre), we first generate a function
to evaluate that landmark:

femur_HC_evaluator = fw_model_landmarks.makeLandmarkEvaluator('femur-HC', gf)

where gf is the geometric_field of a femur. Then to evaluate the the landmark,
we call the generated function with the gf parameters:

femur_HC_coords = femur_HC_evaluator(gf.field_parameters)

It is implemented this way so that landmark coordinates can be evaluated quickly
during fitting optimisations by simpling providing the gf parameters. The
femur_HC_evaluator function can now be called in your script with new gf
parameters to get the new landmark coordinates if the parameters change (e.g.
if the gf has been fitted or transformed).

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import numpy as np

from gias2.common import geoprimitives
from gias2.fieldwork.field import geometric_field

# ===========================================================================#
# femur landmark variables
_femurHeadElem = 0
_femurHeadElems = [0, 1, 2, 3, 4, 5, 6, 7]
_femurShaftElems = [23, 24, 25, 40, 41, 42]
_femurNeckElems = [13, 14, 15, 16]
_femurNeckLongElems = [1, 2, 4, 6, 13, 14, 15, 16]
_femurMedialCondyleElems = [48, 49, 50, 51, 52, 53]
_femurLateralCondyleElems = [43, 44, 45, 46, 47]
_femurGreaterTrochanterElems = [8, 9, 10, 11, 12]
_femurProximalElems = [17, 18, 19, 20, 21, 22] + _femurHeadElems + \
                      _femurGreaterTrochanterElems
_femurDistalElems = [28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39] + \
                    _femurLateralCondyleElems + \
                    _femurMedialCondyleElems
_femurSubtrochanterNodes = [259, 260, 261, 262, 263, 291, 292, 293, 294, 307, 308, 309]
_femurMidShaftNodes = [323, 334, 335, 336, 337, 346, 347, 348, 319, 320, 321, 322, 323]
_femurMedialEpicondyleElems = [53, ]
_femurLateralEpicondyleElems = [45, ]
_femurCondyleAlignmentNodes = [546, 633]
_femurMECNode = 630  # was 633
_femurRightMECNode = 628
_femurLECNode = 550  # was 546
_femurRightLECNode = 550
_femurGTNode = 172
_femurRightGTNode = 174
_femurLTNode = 276
_femurRightLTNode = 276


# _femurLandmarkEvaluators = {
#                   'Head Centre': makeEvaluatorFemurHeadCentre,
#                   'MEC': makeEvaluatorFemurMedialEpicondyle,
#                   'LEC': makeEvaluatorFemurLateralEpicondyle,
#                   'GT': makeEvaluatorFemurGreaterTrochanter,
#                   }

# def getFemurLandmarkNames():
#   return _femurLandmarkEvaluators.keys()

# def makeFemurLandmarkEvaluator(gf, landmarkName):
#   try:
#       return _femurLandmarkEvaluators[landmarkName]
#   except KeyError:
#       raise ValueError, 'Unknown femur landmark name '+landmarkName

def makeEvaluatorFemurHeadCentre(gf, flattened=False, side='left'):
    if flattened:
        headNodes = []
        for e in _femurHeadElems:
            headNodes += list(gf.ensemble_field_function.mapper._element_to_ensemble_map[e].keys())
            headNodes = list(set(headNodes))
    else:
        headNodes = list(gf.ensemble_field_function.mapper._element_to_ensemble_map[_femurHeadElem].keys())

    def evalFemurHeadCentre(meshParams):
        return geoprimitives.fitSphereAnalytic(meshParams[:, headNodes].squeeze().T)[0]

    return evalFemurHeadCentre


def makeEvaluatorFemurMedialEpicondyle(gf, side='left'):
    if side == 'left':
        def evalFemurMedialEpicondyle(meshParams):
            return meshParams[:, _femurMECNode].squeeze()
    elif side == 'right':
        def evalFemurMedialEpicondyle(meshParams):
            return meshParams[:, _femurRightMECNode].squeeze()

    return evalFemurMedialEpicondyle


def makeEvaluatorFemurLateralEpicondyle(gf, side='left'):
    if side == 'left':
        def evalFemurLateralEpicondyle(meshParams):
            return meshParams[:, _femurLECNode].squeeze()
    elif side == 'right':
        def evalFemurLateralEpicondyle(meshParams):
            return meshParams[:, _femurRightLECNode].squeeze()
    return evalFemurLateralEpicondyle


def makeEvaluatorFemurGreaterTrochanter(gf, side='left'):
    if side == 'left':
        def evalFemurGreaterTrochanter(meshParams):
            return meshParams[:, _femurGTNode].squeeze()
    elif side == 'right':
        def evalFemurGreaterTrochanter(meshParams):
            return meshParams[:, _femurRightGTNode].squeeze()

    return evalFemurGreaterTrochanter


def makeEvaluatorFemurLesserTrochanter(gf, side='left'):
    if side == 'left':
        def evalFemurLesserTrochanter(meshParams):
            return meshParams[:, _femurLTNode].squeeze()
    elif side == 'right':
        def evalFemurLesserTrochanter(meshParams):
            return meshParams[:, _femurRightLTNode].squeeze()

    return evalFemurLesserTrochanter


def makeEvaluatorFemurKneeCentre(gf, side='left'):
    evalMC = makeEvaluatorFemurMedialEpicondyle(gf, side=side)
    evalLC = makeEvaluatorFemurLateralEpicondyle(gf, side=side)

    def evalFemurKneeCentre(meshParams):
        mc = evalMC(meshParams)
        lc = evalLC(meshParams)
        return (mc + lc) * 0.5

    return evalFemurKneeCentre


# ===========================================================================#
# hemi pelvis landmark variables
_hemiPelvisAcetabulumElements = [36, 35, 38, 39, 40, 41, 42]
_hemiPelvisASISNode = 464  # left
_hemiPelvisRightASISNode = 466  # right
_hemiPelvisPSISNode = 384  # left
_hemiPelvisRightPSISNode = 384  # right
_hemiPelvisPSNode = 90  # left
_hemiPelvisRightPSNode = 92  # right
_hemiPelvisPTNode = 102  # left
_hemiPelvisRightPTNode = 103  # right
_hemiPelvisISNode = 233  # left
_hemiPelvisRightISNode = 233  # right
_hemiPelvisITNode = 26  # left
_hemiPelvisRightITNode = 25  # right
_hemiPelvisANNode = 287  # left
_hemiPelvisRightANNode = 289  # right


def makeEvaluatorHemiPelvisAcetabularCentre(gf):
    m = gf.ensemble_field_function.mapper._element_to_ensemble_map
    acNodes = []
    for e in _hemiPelvisAcetabulumElements:
        acNodes.extend(list(m[e].keys()))

    def evalAcetabularCentre(meshParams):
        return geoprimitives.fitSphereAnalytic(meshParams[:, acNodes].squeeze().T)[0]

    return evalAcetabularCentre


def makeEvaluatorHemiPelvisASIS(gf, side='left'):
    if side == 'left':
        def evalHemiPelvisASIS(meshParams):
            return meshParams[:, _hemiPelvisASISNode].squeeze()
    elif side == 'right':
        def evalHemiPelvisASIS(meshParams):
            return meshParams[:, _hemiPelvisRightASISNode].squeeze()
    return evalHemiPelvisASIS


def makeEvaluatorHemiPelvisPSIS(gf, side='left'):
    if side == 'left':
        def evalHemiPelvisPSIS(meshParams):
            return meshParams[:, _hemiPelvisPSISNode].squeeze()
    elif side == 'right':
        def evalHemiPelvisPSIS(meshParams):
            return meshParams[:, _hemiPelvisRightPSISNode].squeeze()
    return evalHemiPelvisPSIS


def makeEvaluatorHemiPelvisPubisSymphysis(gf, side='left'):
    if side == 'left':
        def evalHemiPelvisPubisSymphysis(meshParams):
            return meshParams[:, _hemiPelvisPSNode].squeeze()
    elif side == 'right':
        def evalHemiPelvisPubisSymphysis(meshParams):
            return meshParams[:, _hemiPelvisRightPSNode].squeeze()
    return evalHemiPelvisPubisSymphysis


def makeEvaluatorHemiPelvisPubisTubercle(gf, side='left'):
    if side == 'left':
        def evalHemiPelvisPubisTubercle(meshParams):
            return meshParams[:, _hemiPelvisPTNode].squeeze()
    elif side == 'right':
        def evalHemiPelvisPubisTubercle(meshParams):
            return meshParams[:, _hemiPelvisRightPTNode].squeeze()
    return evalHemiPelvisPubisTubercle


def makeEvaluatorHemiPelvisIlialSpine(gf, side='left'):
    if side == 'left':
        def evalHemiPelvisIlialSpine(meshParams):
            return meshParams[:, _hemiPelvisISNode].squeeze()
    elif side == 'right':
        def evalHemiPelvisIlialSpine(meshParams):
            return meshParams[:, _hemiPelvisRightISNode].squeeze()
    return evalHemiPelvisIlialSpine


def makeEvaluatorHemiPelvisIschialTuberosity(gf, side='left'):
    if side == 'left':
        def evalHemiPelvisIschialTuberosity(meshParams):
            return meshParams[:, _hemiPelvisITNode].squeeze()
    elif side == 'right':
        def evalHemiPelvisIschialTuberosity(meshParams):
            return meshParams[:, _hemiPelvisRightITNode].squeeze()
    return evalHemiPelvisIschialTuberosity


def makeEvaluatorHemiPelvisAcetabularNotch(gf, side='left'):
    if side == 'left':
        def evalHemiPelvisAcetabularNotch(meshParams):
            return meshParams[:, _hemiPelvisANNode].squeeze()
    elif side == 'right':
        def evalHemiPelvisAcetabularNotch(meshParams):
            return meshParams[:, _hemiPelvisRightANNode].squeeze()
    return evalHemiPelvisAcetabularNotch


# ===========================================================================#
# whole pelvis landmarks
_pelvisLASISNode = 1004
_pelvisRASISNode = 466
_pelvisLPSISNode = 924
_pelvisRPSISNode = 384
_pelvisLPTNode = 642
_pelvisRPTNode = 103
_pelvisLPSNode = 630
_pelvisRPSNode = 92
_pelvisLISNode = 773
_pelvisRISNode = 233
_pelvisLITNode = 566
_pelvisRITNode = 25
_pelvisSacPlatNode = 1300  # centre-posterior-most point on the vertebral plateau on the sacrum
_pelvisLHJCElems = [109, 111, 112, 113, 114, 115]  # in flattened mesh
_pelvisRHJCElems = [36, 38, 39, 40, 41, 42]  # in flattened mesh
_pelvisRHElems = list(range(0, 73))
_pelvisLHElems = list(range(73, 146))
_pelvisSacElems = list(range(146, 260))


def makeEvaluatorPelvisLASIS(gf, **kwargs):
    def evalPelvisLASIS(meshParams):
        return meshParams[:, _pelvisLASISNode].squeeze()

    return evalPelvisLASIS


def makeEvaluatorPelvisRASIS(gf, **kwargs):
    def evalPelvisRASIS(meshParams):
        return meshParams[:, _pelvisRASISNode].squeeze()

    return evalPelvisRASIS


def makeEvaluatorPelvisLPSIS(gf, **kwargs):
    def evalPelvisLPSIS(meshParams):
        return meshParams[:, _pelvisLPSISNode].squeeze()

    return evalPelvisLPSIS


def makeEvaluatorPelvisRPSIS(gf, **kwargs):
    def evalPelvisRPSIS(meshParams):
        return meshParams[:, _pelvisRPSISNode].squeeze()

    return evalPelvisRPSIS


def makeEvaluatorPelvisLPT(gf, **kwargs):
    def evalPelvisLPT(meshParams):
        return meshParams[:, _pelvisLPTNode].squeeze()

    return evalPelvisLPT


def makeEvaluatorPelvisRPT(gf, **kwargs):
    def evalPelvisRPT(meshParams):
        return meshParams[:, _pelvisRPTNode].squeeze()

    return evalPelvisRPT


def makeEvaluatorPelvisLPS(gf, **kwargs):
    def evalPelvisLPS(meshParams):
        return meshParams[:, _pelvisLPSNode].squeeze()

    return evalPelvisLPS


def makeEvaluatorPelvisRPS(gf, **kwargs):
    def evalPelvisRPS(meshParams):
        return meshParams[:, _pelvisRPSNode].squeeze()

    return evalPelvisRPS


def makeEvaluatorPelvisLIS(gf, **kwargs):
    def evalPelvisLIS(meshParams):
        return meshParams[:, _pelvisLISNode].squeeze()

    return evalPelvisLIS


def makeEvaluatorPelvisRIS(gf, **kwargs):
    def evalPelvisRIS(meshParams):
        return meshParams[:, _pelvisRISNode].squeeze()

    return evalPelvisRIS


def makeEvaluatorPelvisLIT(gf, **kwargs):
    def evalPelvisLIT(meshParams):
        return meshParams[:, _pelvisLITNode].squeeze()

    return evalPelvisLIT


def makeEvaluatorPelvisRIT(gf, **kwargs):
    def evalPelvisRIT(meshParams):
        return meshParams[:, _pelvisRITNode].squeeze()

    return evalPelvisRIT


def makeEvaluatorPelvisSacPlat(gf, **kwargs):
    def evalPelvisSacPlat(meshParams):
        return meshParams[:, _pelvisSacPlatNode].squeeze()

    return evalPelvisSacPlat


def makeEvaluatorPelvisSacral(gf, **kwargs):
    """Mid-point of PSISes
    """

    def evalPelvisSacral(meshParams):
        s = 0.5 * (meshParams[:, _pelvisLPSISNode].squeeze() +
                   meshParams[:, _pelvisRPSISNode].squeeze()
                   )
        return s

    return evalPelvisSacral


def makeEvaluatorPelvisLHJC(gf, disc=5.0, radius=False, side=None):
    # make evaluator for left acetabulum elements
    acetabElemEval = geometric_field.makeGeometricFieldElementsEvaluatorSparse(
        gf, _pelvisLHJCElems, disc
    )
    if radius:
        def evalPelvisLHJC(meshParams):
            acetabPoints = acetabElemEval(meshParams).T
            return geoprimitives.fitSphereAnalytic(acetabPoints)
    else:
        def evalPelvisLHJC(meshParams):
            acetabPoints = acetabElemEval(meshParams).T
            return geoprimitives.fitSphereAnalytic(acetabPoints)[0]
    return evalPelvisLHJC


def makeEvaluatorPelvisRHJC(gf, disc=5.0, radius=False, side=None):
    # make evaluator for left acetabulum elements
    acetabElemEval = geometric_field.makeGeometricFieldElementsEvaluatorSparse(
        gf, _pelvisRHJCElems, disc
    )
    if radius:
        def evalPelvisRHJC(meshParams):
            acetabPoints = acetabElemEval(meshParams).T
            return geoprimitives.fitSphereAnalytic(acetabPoints)
    else:
        def evalPelvisRHJC(meshParams):
            acetabPoints = acetabElemEval(meshParams).T
            return geoprimitives.fitSphereAnalytic(acetabPoints)[0]
    return evalPelvisRHJC


# ===========================================================================#
# tibia fibula combined landmarks

# aligned with opensim tibia model
_tibiaFibulaLCNode = 257  # 256
_tibiaFibulaRightLCNode = 258
_tibiaFibulaMCNode = 236  # s235
_tibiaFibulaRightMCNode = 235  # 237
_tibiaFibulaLMNode = 528
_tibiaFibulaRightLMNode = 527  # 535
_tibiaFibulaMMNode = 150
_tibiaFibulaRightMMNode = 151

_tibiaFibulaTTNode = 203
_tibiaFibulaRightTTNode = 203
_tibiaFibulaKneeCentreOffset = 50.0  # 389.55 from KC to ankle average of 45 subjects from Mousa K.
_tibiaElements = list(range(0, 46))
_fibulaElements = list(range(46, 88))


def makeEvaluatorTibiaFibulaLC(gf, side='left'):
    if side == 'left':
        def evalTibiaFibulaLC(meshParams):
            return meshParams[:, _tibiaFibulaLCNode].squeeze()
    elif side == 'right':
        def evalTibiaFibulaLC(meshParams):
            return meshParams[:, _tibiaFibulaRightLCNode].squeeze()
    return evalTibiaFibulaLC


def makeEvaluatorTibiaFibulaMC(gf, side='left'):
    if side == 'left':
        def evalTibiaFibulaMC(meshParams):
            return meshParams[:, _tibiaFibulaMCNode].squeeze()
    elif side == 'right':
        def evalTibiaFibulaMC(meshParams):
            return meshParams[:, _tibiaFibulaRightMCNode].squeeze()
    return evalTibiaFibulaMC


def makeEvaluatorTibiaFibulaMM(gf, side='left'):
    if side == 'left':
        def evalTibiaFibulaMM(meshParams):
            return meshParams[:, _tibiaFibulaMMNode].squeeze()
    elif side == 'right':
        def evalTibiaFibulaMM(meshParams):
            return meshParams[:, _tibiaFibulaRightMMNode].squeeze()
    return evalTibiaFibulaMM


def makeEvaluatorTibiaFibulaLM(gf, side='left'):
    if side == 'left':
        def evalTibiaFibulaLM(meshParams):
            return meshParams[:, _tibiaFibulaLMNode].squeeze()
    elif side == 'right':
        def evalTibiaFibulaLM(meshParams):
            return meshParams[:, _tibiaFibulaRightLMNode].squeeze()
    return evalTibiaFibulaLM


def makeEvaluatorTibiaFibulaTT(gf, side='left'):
    if side == 'left':
        def evalTibiaFibulaTT(meshParams):
            return meshParams[:, _tibiaFibulaTTNode].squeeze()
    elif side == 'right':
        def evalTibiaFibulaTT(meshParams):
            return meshParams[:, _tibiaFibulaRightTTNode].squeeze()
    return evalTibiaFibulaTT


def makeEvaluatorTibiaFibulaKneeCentre(gf, side='left'):
    evalLC = makeEvaluatorTibiaFibulaLC(gf, side=side)
    evalMC = makeEvaluatorTibiaFibulaMC(gf, side=side)
    evalLM = makeEvaluatorTibiaFibulaLM(gf, side=side)
    evalMM = makeEvaluatorTibiaFibulaMM(gf, side=side)

    def evalTibiaFibulaKneeCentre(meshParams):
        lc = evalLC(meshParams)
        mc = evalMC(meshParams)
        lm = evalLM(meshParams)
        mm = evalMM(meshParams)

        # calc tibfib ACS
        ic = (mc + lc) / 2.0
        im = (mm + lm) / 2.0  # origin
        # superiorly, IM to IC
        y = geoprimitives.norm(ic - im)
        # anteriorly, normal to plane of IM, LC and MC
        x = geoprimitives.norm(np.cross(lc - im, mc - im))
        # right
        z = geoprimitives.norm(np.cross(x, y))

        # estimate knee centre
        kc = ic + y * _tibiaFibulaKneeCentreOffset
        return kc

    return evalTibiaFibulaKneeCentre


# ===========================================================================#
# patella landmarks
_patellaInfNode = 29
_patellaRightInfNode = 29
_patellaSupNode = 59
_patellaRightSupNode = 58
_patellaLatNode = 72
_patellaRightLatNode = 72


def makeEvaluatorPatellaInf(gf, side='left'):
    if side == 'left':
        def evalPatellaInf(meshParams):
            return meshParams[:, _patellaInfNode].squeeze()
    elif side == 'right':
        def evalPatellaInf(meshParams):
            return meshParams[:, _patellaRightInfNode].squeeze()
    return evalPatellaInf


def makeEvaluatorPatellaSup(gf, side='left'):
    if side == 'left':
        def evalPatellaSup(meshParams):
            return meshParams[:, _patellaSupNode].squeeze()
    elif side == 'right':
        def evalPatellaSup(meshParams):
            return meshParams[:, _patellaRightSupNode].squeeze()
    return evalPatellaSup


def makeEvaluatorPatellaLat(gf, side='left'):
    if side == 'left':
        def evalPatellaLat(meshParams):
            return meshParams[:, _patellaLatNode].squeeze()
    elif side == 'right':
        def evalPatellaLat(meshParams):
            return meshParams[:, _patellaRightLatNode].squeeze()
    return evalPatellaLat


# ===========================================================================#
_landmarkEvaluators = {
    'femur-HC': makeEvaluatorFemurHeadCentre,
    'femur-MEC': makeEvaluatorFemurMedialEpicondyle,
    'femur-LEC': makeEvaluatorFemurLateralEpicondyle,
    'femur-GT': makeEvaluatorFemurGreaterTrochanter,
    'femur-LT': makeEvaluatorFemurLesserTrochanter,
    'femur-kneecentre': makeEvaluatorFemurKneeCentre,
    'hpelvis-ASIS': makeEvaluatorHemiPelvisASIS,
    'hpelvis-PSIS': makeEvaluatorHemiPelvisPSIS,
    'hpelvis-PS': makeEvaluatorHemiPelvisPubisSymphysis,
    'hpelvis-PT': makeEvaluatorHemiPelvisPubisTubercle,
    'hpelvis-IS': makeEvaluatorHemiPelvisIlialSpine,
    'hpelvis-IT': makeEvaluatorHemiPelvisIschialTuberosity,
    'hpelvis-AN': makeEvaluatorHemiPelvisAcetabularNotch,
    'pelvis-LASIS': makeEvaluatorPelvisLASIS,
    'pelvis-RASIS': makeEvaluatorPelvisRASIS,
    'pelvis-LPSIS': makeEvaluatorPelvisLPSIS,
    'pelvis-RPSIS': makeEvaluatorPelvisRPSIS,
    'pelvis-Sacral': makeEvaluatorPelvisSacral,
    'pelvis-LPS': makeEvaluatorPelvisLPS,
    'pelvis-RPS': makeEvaluatorPelvisRPS,
    'pelvis-LIS': makeEvaluatorPelvisLIS,
    'pelvis-RIS': makeEvaluatorPelvisRIS,
    'pelvis-LIT': makeEvaluatorPelvisLIT,
    'pelvis-RIT': makeEvaluatorPelvisRIT,
    'pelvis-LHJC': makeEvaluatorPelvisLHJC,
    'pelvis-RHJC': makeEvaluatorPelvisRHJC,
    'pelvis-SacPlat': makeEvaluatorPelvisSacPlat,
    'tibiafibula-LC': makeEvaluatorTibiaFibulaLC,
    'tibiafibula-MC': makeEvaluatorTibiaFibulaMC,
    'tibiafibula-LM': makeEvaluatorTibiaFibulaLM,
    'tibiafibula-MM': makeEvaluatorTibiaFibulaMM,
    'tibiafibula-TT': makeEvaluatorTibiaFibulaTT,
    'tibiafibula-kneecentre': makeEvaluatorTibiaFibulaKneeCentre,
    'patella-inf': makeEvaluatorPatellaInf,
    'patella-sup': makeEvaluatorPatellaSup,
    'patella-lat': makeEvaluatorPatellaLat,
}

validLandmarks = sorted(_landmarkEvaluators.keys())


def landmarkNames():
    """Return a list of implemented landmarks
    """
    return list(_landmarkEvaluators.keys())


def makeLandmarkEvaluator(name, gf, **args):
    """
    Generate a function to evaluate the named landmark on the given
    geometric_field. Call landmarkNames to get a list of possible landmark
    names.

    inputs
    ------
    name : str of landmark name.
    gf : geometric_field instance on which to evaluate the landmark

    returns
    -------
    func : function that evaluates the named landmark given the field
        parameters of gf. e.g. ldmk_coords = func(gf.field_parameters)
    """
    if args is None:
        args = {}
    try:
        return _landmarkEvaluators[name](gf, **args)
    except KeyError:
        raise ValueError('Unknown landmark name ' + name)
