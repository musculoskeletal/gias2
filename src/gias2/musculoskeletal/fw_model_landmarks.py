"""
FILE: fw_model_landmark.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: Functions for creating evaluator functions for anatomic landmarks
on fieldwork models

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

#===========================================================================#
# femur landmark variables
_femurHeadElem = 0
_femurHeadElems = [0,1,2,3,4,5,6,7]
_femurShaftElems = [23,24,25,40,41,42]
_femurNeckElems = [13,14,15,16]
_femurNeckLongElems = [1,2,4,6,13,14,15,16]
_femurMedialCondyleElems = [48,49,50,51,52,53]
_femurLateralCondyleElems = [43,44,45,46,47]
_femurGreaterTrochanterElems = [8,9,10,11,12]
_femurProximalElems = [17,18,19,20,21,22]+_femurHeadElems+\
					  _femurGreaterTrochanterElems
_femurDistalElems = [28,29,30,31,32,33,34,35,36,37,38,39]+\
					_femurLateralCondyleElems+\
					_femurMedialCondyleElems
_femurSubtrochanterNodes = [259,260,261,262,263,291,292,293,294,307,308,309]
_femurMidShaftNodes = [323,334,335,336,337,346,347,348,319, 320,321,322,323]
_femurMedialEpicondyleElems = [53,]
_femurLateralEpicondyleElems = [45,]
_femurCondyleAlignmentNodes = [546, 633] 
_femurMECNode = 633
_femurLECNode = 546
_femurGTNode = 172
# _femurLandmarkEvaluators = {
# 					'Head Centre': makeEvaluatorFemurHeadCentre,
# 					'MEC': makeEvaluatorFemurMedialEpicondyle,
# 					'LEC': makeEvaluatorFemurLateralEpicondyle,
# 					'GT': makeEvaluatorFemurGreaterTrochanter,
# 					}

# def getFemurLandmarkNames():
# 	return _femurLandmarkEvaluators.keys()

# def makeFemurLandmarkEvaluator(gf, landmarkName):
# 	try:
# 		return _femurLandmarkEvaluators[landmarkName]
# 	except KeyError:
# 		raise ValueError, 'Unknown femur landmark name '+landmarkName

def makeEvaluatorFemurHeadCentre(gf, flattened=False):
	if flattened:
		headNodes = []
		for e in _femurHeadElemFlat:
			headNodes += gf.ensemble_field_function.mapper._element_to_ensemble_map[e].keys()
			headNodes = list(set(headNodes))
	else:
		headNodes = gf.ensemble_field_function.mapper._element_to_ensemble_map[_femurHeadElem].keys()
	def evalFemurHeadCentre(meshParams):
		return geoprimitives.fitSphereAnalytic(meshParams[:,headNodes].squeeze().T)[0]
	return evalFemurHeadCentre

def makeEvaluatorFemurMedialEpicondyle(gf):
	def evalFemurMedialEpicondyle(meshParams):
		return meshParams[:,_femurMECNode].squeeze()
	return evalFemurMedialEpicondyle

def makeEvaluatorFemurLateralEpicondyle(gf):
	def evalFemurLateralEpicondyle(meshParams):
		return meshParams[:,_femurLECNode].squeeze()
	return evalFemurLateralEpicondyle

def makeEvaluatorFemurGreaterTrochanter(gf):
	def evalFemurGreaterTrochanter(meshParams):
		return meshParams[:,_femurGTNode].squeeze()
	return evalFemurGreaterTrochanter

def makeEvaluatorFemurKneeCentre(gf):
	evalMC = makeEvaluatorFemurMedialEpicondyle(gf)
	evalLC = makeEvaluatorFemurLateralEpicondyle(gf)
	def evalFemurKneeCentre(meshParams):
		mc = evalMC(meshParams)
		lc = evalLC(meshParams)
		return (mc+lc)*0.5
	return evalFemurKneeCentre

#===========================================================================#
# hemi pelvis landmark variables
_hemiPelvisAcetabulumElements = [36,35,38,39,40,41,42]
_hemiPelvisASISNode = 464 # left
_hemiPelvisRightASISNode = 466 # right
_hemiPelvisPSISNode = 384 # left
_hemiPelvisRightPSISNode = 384 # right
_hemiPelvisPSNode = 90 # left
_hemiPelvisRightPSNode = 92 # right
_hemiPelvisPTNode = 102 # left
_hemiPelvisRightPTNode = 103 # right
_hemiPelvisISNode = 233 # left
_hemiPelvisRightISNode = 233 # right
_hemiPelvisITNode = 26 # left
_hemiPelvisRightITNode = 25 # right
_hemiPelvisANNode = 287 # left
_hemiPelvisRightANNode = 289 # right

def makeEvaluatorHemiPelvisAcetabularCentre(gf):
	m = gf.ensemble_field_function.mapper._element_to_ensemble_map
	acNodes = []
	for e in _hemiPelvisAcetabulumElements:
		acNodes.extend(m[e].keys())
	def evalAcetabularCentre(meshParams):
		return geoprimitives.fitSphereAnalytic(meshParams[:,acNodes].squeeze().T)[0]
	return evalAcetabularCentre

def makeEvaluatorHemiPelvisASIS(gf):
	def evalHemiPelvisASIS(meshParams):
		return meshParams[:,_hemiPelvisASISNode].squeeze()
	return evalHemiPelvisASIS

def makeEvaluatorHemiPelvisPSIS(gf):
	def evalHemiPelvisPSIS(meshParams):
		return meshParams[:,_hemiPelvisPSISNode].squeeze()
	return evalHemiPelvisPSIS

def makeEvaluatorHemiPelvisPubisSymphysis(gf):
	def evalHemiPelvisPubisSymphysis(meshParams):
		return meshParams[:,_hemiPelvisPSNode].squeeze()
	return evalHemiPelvisPubisSymphysis

def makeEvaluatorHemiPelvisPubisTubercle(gf):
	def evalHemiPelvisPubisTubercle(meshParams):
		return meshParams[:,_hemiPelvisPTNode].squeeze()
	return evalHemiPelvisPubisTubercle

def makeEvaluatorHemiPelvisIlialSpine(gf):
	def evalHemiPelvisIlialSpine(meshParams):
		return meshParams[:,_hemiPelvisISNode].squeeze()
	return evalHemiPelvisIlialSpine

def makeEvaluatorHemiPelvisIschialTuberosity(gf):
	def evalHemiPelvisIschialTuberosity(meshParams):
		return meshParams[:,_hemiPelvisITNode].squeeze()
	return evalHemiPelvisIschialTuberosity

def makeEvaluatorHemiPelvisAcetabularNotch(gf):
	def evalHemiPelvisAcetabularNotch(meshParams):
		return meshParams[:,_hemiPelvisANNode].squeeze()
	return evalHemiPelvisAcetabularNotch


#===========================================================================#
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
_pelvisLHJCElems = [109,111,112,113,114,115] # in flattened mesh
_pelvisRHJCElems = [36,38,39,40,41,42] # in flattened mesh
_pelvisRHElems = range(0, 73)
_pelvisLHElems = range(73,146)
_pelvisSacElems = range(146, 260)

def makeEvaluatorPelvisLASIS(gf):
	def evalPelvisLASIS(meshParams):
		return meshParams[:,_pelvisLASISNode].squeeze()
	return evalPelvisLASIS

def makeEvaluatorPelvisRASIS(gf):
	def evalPelvisRASIS(meshParams):
		return meshParams[:,_pelvisRASISNode].squeeze()
	return evalPelvisRASIS

def makeEvaluatorPelvisLPSIS(gf):
	def evalPelvisLPSIS(meshParams):
		return meshParams[:,_pelvisLPSISNode].squeeze()
	return evalPelvisLPSIS

def makeEvaluatorPelvisRPSIS(gf):
	def evalPelvisRPSIS(meshParams):
		return meshParams[:,_pelvisRPSISNode].squeeze()
	return evalPelvisRPSIS

def makeEvaluatorPelvisLPT(gf):
	def evalPelvisLPT(meshParams):
		return meshParams[:,_pelvisLPTNode].squeeze()
	return evalPelvisLPT

def makeEvaluatorPelvisRPT(gf):
	def evalPelvisRPT(meshParams):
		return meshParams[:,_pelvisRPTNode].squeeze()
	return evalPelvisRPT

def makeEvaluatorPelvisLPS(gf):
	def evalPelvisLPS(meshParams):
		return meshParams[:,_pelvisLPSNode].squeeze()
	return evalPelvisLPS

def makeEvaluatorPelvisRPS(gf):
	def evalPelvisRPS(meshParams):
		return meshParams[:,_pelvisRPSNode].squeeze()
	return evalPelvisRPS

def makeEvaluatorPelvisLIS(gf):
	def evalPelvisLIS(meshParams):
		return meshParams[:,_pelvisLISNode].squeeze()
	return evalPelvisLIS

def makeEvaluatorPelvisRIS(gf):
	def evalPelvisRIS(meshParams):
		return meshParams[:,_pelvisRISNode].squeeze()
	return evalPelvisRIS

def makeEvaluatorPelvisLIT(gf):
	def evalPelvisLIT(meshParams):
		return meshParams[:,_pelvisLITNode].squeeze()
	return evalPelvisLIT

def makeEvaluatorPelvisRIT(gf):
	def evalPelvisRIT(meshParams):
		return meshParams[:,_pelvisRITNode].squeeze()
	return evalPelvisRIT

def makeEvaluatorPelvisSacral(gf):
	"""Mid-point of PSISes
	"""
	def evalPelvisSacral(meshParams):
		s = 0.5*(meshParams[:,_pelvisLPSISNode].squeeze()+
				 meshParams[:,_pelvisRPSISNode].squeeze()
				 )
		return s
	return evalPelvisSacral

def makeEvaluatorPelvisLHJC(gf, disc=5.0, radius=False):
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

def makeEvaluatorPelvisRHJC(gf, disc=5.0, radius=False):
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

#===========================================================================#
# tibia fibula combined landmarks
_tibiaFibulaLCNode = 256
_tibiaFibulaMCNode = 235
_tibiaFibulaLMNode = 528
_tibiaFibulaMMNode = 150
_tibiaFibulaTTNode  = 203
_tibiaFibulaKneeCentreOffset = 50.0 # 389.55 from KC to ankle average of 45 subjects from Mousa K.
_tibiaElements = range(0, 46)
_fibulaElements = range(46, 88)

def makeEvaluatorTibiaFibulaLC(gf):
	def evalTibiaFibulaLC(meshParams):
		return meshParams[:,_tibiaFibulaLCNode].squeeze()
	return evalTibiaFibulaLC

def makeEvaluatorTibiaFibulaMC(gf):
	def evalTibiaFibulaMC(meshParams):
		return meshParams[:,_tibiaFibulaMCNode].squeeze()
	return evalTibiaFibulaMC

def makeEvaluatorTibiaFibulaMM(gf):
	def evalTibiaFibulaMM(meshParams):
		return meshParams[:,_tibiaFibulaMMNode].squeeze()
	return evalTibiaFibulaMM

def makeEvaluatorTibiaFibulaLM(gf):
	def evalTibiaFibulaLM(meshParams):
		return meshParams[:,_tibiaFibulaLMNode].squeeze()
	return evalTibiaFibulaLM

def makeEvaluatorTibiaFibulaTT(gf):
	def evalTibiaFibulaTT(meshParams):
		return meshParams[:,_tibiaFibulaTTNode].squeeze()
	return evalTibiaFibulaTT

def makeEvaluatorTibiaFibulaKneeCentre(gf):
	evalLC = makeEvaluatorTibiaFibulaLC(gf)
	evalMC = makeEvaluatorTibiaFibulaMC(gf)
	evalLM = makeEvaluatorTibiaFibulaLM(gf)
	evalMM = makeEvaluatorTibiaFibulaMM(gf)
	def evalTibiaFibulaKneeCentre(meshParams):
		lc = evalLC(meshParams)
		mc = evalMC(meshParams)
		lm = evalLM(meshParams)
		mm = evalMM(meshParams)

		# calc tibfib ACS
		ic = (mc + lc)/2.0
		im = (mm + lm)/2.0 # origin
		# superiorly, IM to IC
		y = geoprimitives.norm(ic-im)
		# anteriorly, normal to plane of IM, LC and MC
		x = geoprimitives.norm(np.cross(lc-im, mc-im))
		# right
		z = geoprimitives.norm(np.cross(x,y))

		# estimate knee centre
		kc = ic + y*_tibiaFibulaKneeCentreOffset
		return kc

	return evalTibiaFibulaKneeCentre

#===========================================================================#
# patella landmarks
_patellaInfNode = 29
_patellaSupNode = 59
_patellaLatNode = 72

def makeEvaluatorPatellaInf(gf):
	def evalPatellaInf(meshParams):
		return meshParams[:,_patellaInfNode].squeeze()
	return evalPatellaInf

def makeEvaluatorPatellaSup(gf):
	def evalPatellaSup(meshParams):
		return meshParams[:,_patellaSupNode].squeeze()
	return evalPatellaSup

def makeEvaluatorPatellaLat(gf):
	def evalPatellaLat(meshParams):
		return meshParams[:,_patellaLatNode].squeeze()
	return evalPatellaLat

#===========================================================================#
_landmarkEvaluators = {
				  'femur-HC': makeEvaluatorFemurHeadCentre,
				  'femur-MEC': makeEvaluatorFemurMedialEpicondyle,
				  'femur-LEC': makeEvaluatorFemurLateralEpicondyle,
				  'femur-GT': makeEvaluatorFemurGreaterTrochanter,
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

def landmarkNames():
	return _landmarkEvaluators.keys()

def makeLandmarkEvaluator(name, gf, **args):
	if args is None:
		args = {}
	try:
		return _landmarkEvaluators[name](gf, **args)
	except KeyError:
		raise ValueError, 'Unknown landmark name '+name

