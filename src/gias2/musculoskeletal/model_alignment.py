"""
FILE: model_alignment.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: functions for alignining fieldwork models of individual bones

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import copy
import scipy
from gias2.common import geoprimitives
from gias2.registration import alignment_fitting
from gias2.registration import alignment_analytic
from gias2.common import transform3D
from gias2.learning import PCA_fitting
from gias2.musculoskeletal import fw_femur_model_data as fmd
from gias2.musculoskeletal import fw_femur_measurements
from gias2.musculoskeletal import fw_model_landmarks

def normaliseVector(v):
	return v/scipy.linalg.norm(v)

def _makeLandmarkObj(targ, evaluator):
    def obj(P):
        # print targ
        return ((targ - evaluator(P))**2.0).sum()

    return obj

#=======================================================#
# general alignment                                     #
#=======================================================#
def alignMeshParametersRigid( gFields, targetGF=None, retTransforms=False ):
	# procrustes alignment of gFields to the 1st gField

	# evaluate points from each g
	# d = 5
	# X = [ g.evaluate_geometric_field( d ).T for g in gFields ]
	XNodes = [ g.get_all_point_positions() for g in gFields ]
	
	if targetGF==None:
		# use the first one
		targetGF = gFields[0]
		
	# targ = targetGF.evaluate_geometric_field( d ).T
	targNodes = targetGF.get_all_point_positions()
	targetCoM = targetGF.calc_CoM()

	# rigid fit each to X[targI] 
	alignedParams = []
	Ts = []
	for i in xrange(len(gFields)):

		CoMTrans = targetCoM - gFields[i].calc_CoM()
		x0 = scipy.hstack([CoMTrans,0,0,0])
		# fit nodes
		tOpt = alignment_fitting.fitRigid( XNodes[i], targNodes, xtol=1e-6, verbose=1 )[0]
		# fit surface data
		# tOpt = alignment_fitting.fitRigid( X[i], targ, xtol=1e-5, verbose=1 )[0]

		# apply transform to gfield parameters
		gFieldNodes = gFields[i].get_field_parameters().squeeze().T
		alignedParams.append( transform3D.transformRigid3DAboutCoM( gFieldNodes, tOpt ).T[:,:,scipy.newaxis] )
		Ts.append(tOpt)
	
	if retTransforms:
		return alignedParams, Ts
	else:
		return alignedParams

def alignMeshParametersProcrustes( gFields, targetGF=None, retTransforms=False ):
	# procrustes alignment of gFields to the 1st gField

	# evaluate points from each g
	# d = 5
	# X = [ g.evaluate_geometric_field( d ).T for g in gFields ]
	XNodes = [ g.get_all_point_positions() for g in gFields ]
	
	if targetGF==None:
		# use the first one
		targetGF = gFields[0]
		
	# targ = targetGF.evaluate_geometric_field( d ).T
	targNodes = targetGF.get_all_point_positions()
	
	# rigid fit each to X[targI] 
	sizes = []
	alignedParams = []
	Ts = []
	for i in xrange(len(gFields)):
		# fit nodes
		tOpt = alignment_fitting.fitRigidScale( XNodes[i], targNodes, xtol=1e-3, verbose=1 )[0]
		# fit surface data
		# tOpt = alignment_fitting.fitRigidScale( X[i], targ, xtol=1e-5 )[0]
		
		# apply transform to gfield parameters
		gFieldNodes = gFields[i].get_field_parameters().squeeze().T
		alignedParams.append( transform3D.transformRigidScale3DAboutCoM( gFieldNodes, tOpt ).T[:,:,scipy.newaxis] )
		sizes.append( tOpt[-1] )
		Ts.append(tOpt)

	if retTransforms:
		return alignedParams, scipy.array(sizes), Ts
	else:	
		return alignedParams, scipy.array(sizes)

#=======================================================#
# femur alignment                                       #
#=======================================================#
def createFemurACS(head, mc, lc):
	o = ( mc + lc )/2.0
	z = normaliseVector( head - o )
	y = normaliseVector( scipy.cross( z,  ( o - lc ) ) )
	x = normaliseVector( scipy.cross( y, z ) )
	u = scipy.array([ o, o+x, o+y, o+z])
	return u

def createFemurACSISB(head, mc, lc, side='left'):
	"""Axes: x-anterior, y-superior, z-right
	origin: midpoint of epicondyles
	"""
	#  origin - midpoint of epicondyles
	o = (mc + lc)/2.0
	# y - origin to head
	y = normaliseVector(head - o)
	# z - right in plane of head, mc, lc
	n1 = normaliseVector(scipy.cross(mc-head, lc-head))
	z = normaliseVector(scipy.cross(n1, y))
	if side=='right':
		z *= -1.0
	# x - anteriorly 
	x = normaliseVector(scipy.cross(y, z))
	return o, x, y, z

def alignAnatomicFemur( X, head, mc, lc, returnT=False ):
	""" aligns points X, with head CoM, mc CoM, lc CoM, to the origin
	and global axes (femur only). Grood and Suntay 1983 system.
	"""

	o = ( mc + lc )/2.0
	z = normaliseVector( head - o )
	y = normaliseVector( scipy.cross( z,  ( o - lc ) ) )
	x = normaliseVector( scipy.cross( y, z ) )
	
	# o = X.mean(0)
	
	u = scipy.array([ o, o+x, o+y, o+z])
				
	ut = scipy.array([[0, 0, 0],\
				[1, 0, 0],\
				[0, 1, 0],\
				[0, 0, 1]])
				
	t = transform3D.directAffine( u, ut )
	if t.shape==(3,4):
		t = scipy.vstack([t, [0,0,0,1]])
	
	if returnT:
		return transform3D.transformAffine( X, t ), t
	else:
		return transform3D.transformAffine( X, t )

femurLandmarkNodes = {'MEC': 633,
					  'LEC': 546,
					  'FGT': 172,
					  }

def alignFemurLandmarksRigidScale(gf, landmarks, t0=None, r0=None, s0=None):
	"""
	landmarks: a list of tuples [(landmark name, landmark coords),...]
	valid landmark names: FHC, MEC, LEC

	"""
	targetLandmarks = []
	sourceLandmarks = []
	for ldName, ldTarg in landmarks:
		if ldName is 'FHC':
			ldName = 'HC'
		evaluator = fw_model_landmarks.makeLandmarkEvaluator('femur-'+ldName, gf)
		sourceLandmarks.append(evaluator(gf.get_field_parameters()))
		targetLandmarks.append(ldTarg)

	targetLandmarks = scipy.array(targetLandmarks)
	sourceLandmarks = scipy.array(sourceLandmarks)
	T0 = scipy.zeros(7)
	if s0 is None:
		T0[6] = 1.0
	else:
		T0[6] = s0

	if t0 is None:
		T0[:3] = targetLandmarks.mean(0) - sourceLandmarks.mean(0)
	else:
		T0[:3] = t0[:]

	if r0 is not None:
		T0[3:6] = r0[:]


	TOpt, fittedLandmarks, (rms0, rmsOpt) = alignment_fitting.fitRigidScale(
												sourceLandmarks,
								   				targetLandmarks,
								   				t0=T0,
								   				xtol=1e-9,
								   				outputErrors=1)

	gf.transformRigidScaleRotateAboutP(TOpt, sourceLandmarks.mean(0))

	return gf, (rms0, rmsOpt), TOpt 

def alignFemurLandmarksPC(gf, pc, landmarks, GFParamsCallback=None, mw0=1.0, mwn=1.0):
	"""
	landmarks: a list of tuples [(landmark name, landmark coords),...]
	valid landmark names: FHC, MEC, LEC

	"""
	headElem = 0

	sourceGF = copy.deepcopy(gf)
	headNodes = sourceGF.ensemble_field_function.mapper._element_to_ensemble_map[headElem].keys()
	hasFHC = False
	targetFHC = None
	targetLandmarks = []
	landmarkNodes = []
	for ln, lc in landmarks:
		if ln=='FHC':
			hasFHC = True
			targetFHC = lc
		elif ln in femurLandmarkNodes:
			targetLandmarks.append(lc)
			landmarkNodes.append(femurLandmarkNodes[ln])
		else:
			print 'WARNING: landmark %s unsupported'

	if hasFHC:
		targetLandmarks.append(targetFHC)

	def obj(p):
		P = p.reshape((3,-1))
		x = scipy.array([P[:,l] for l in landmarkNodes])
		if hasFHC:
			headC = geoprimitives.fitSphereAnalytic(P[:,headNodes].T)[0]
			x = scipy.vstack([x, headC])
		# e = scipy.sqrt((((x - targetLandmarks)**2.0).sum(1)).mean())
		sse = (((x - targetLandmarks)**2.0).sum(1)).sum()
		return sse

	pcFitter = PCA_fitting.PCFit(pc=pc)
	pcFitter.useFMin = True
	pcFitter.ftol = 1e-3

	P0 = pc.getMean().reshape((3,-1))
	n0 = P0[:,landmarkNodes[0]]
	x0 = scipy.hstack([targetLandmarks[0] - n0, 0, 0, 0])

	# targetCoM = (targetHC + ((targetMEC+targetLEC)/2.0))/2.0
	# x0 = scipy.hstack([targetCoM - sourceGF.calc_CoM(), 0, 0, 0])
	
	rigidT, rigidP = pcFitter.rigidFit(obj, x0=x0)
	sourceGF.set_field_parameters(rigidP.reshape((3,-1,1)))
	rigidSSE = obj(rigidP)
	if GFParamsCallback is not None:
		GFParamsCallback(rigidP)

	rigidMode0T, rigidMode0P = pcFitter.rigidMode0Fit(obj, mWeight=mw0)
	sourceGF.set_field_parameters(rigidMode0P.reshape((3,-1,1)))
	rigidMode0SSE = obj(rigidMode0P)
	if GFParamsCallback is not None:
		GFParamsCallback(rigidMode0P)
	
	rigidModeNT, rigidModeNP = pcFitter.rigidModeNFit(obj, modes=[1,2], mWeight=mwn)
	sourceGF.set_field_parameters(rigidModeNP.reshape((3,-1,1)))
	rigidModeNSSE = obj(rigidModeNP)
	if GFParamsCallback is not None:
		GFParamsCallback(rigidModeNP)

	return sourceGF, (rigidSSE, rigidMode0SSE, rigidModeNSSE), rigidModeNT 

def alignAnatomicFemurOrthoload( X, head, p1, p2, lcdorsal, mcdorsal, returnT=False ):
	""" aligns points X, with head CoM, mc CoM, lc CoM, to the origin
	and global axes (femur only). Grood and Suntay 1983 system.
	"""

	o = head
	z = normaliseVector( p1 - p2 )
	y = normaliseVector( scipy.cross( z,  ( lcdorsal - mcdorsal ) ) )
	x = normaliseVector( scipy.cross( y, z ) )
	
	# o = X.mean(0)
	
	u = scipy.array([ o, o+x, o+y, o+z])
				
	ut = scipy.array([[0, 0, 0],\
				[1, 0, 0],\
				[0, 1, 0],\
				[0, 0, 1]])
				
	t = transform3D.directAffine( u, ut )
	if t.shape==(3,4):
		t = scipy.vstack([t, [0,0,0,1]])
	
	if returnT:
		return transform3D.transformAffine( X, t ), t
	else:
		return transform3D.transformAffine( X, t )

def alignFemurMeshParametersOrtholoadSingle( femurModel ):
	""" given a femur geometric field, align it geometrically.
	returns the aligned field parameters
	""" 

	# first align to standard ACS
	femurParamsACS, femurACST = alignFemurMeshParametersAnatomicSingle(femurModel)
	femurModel.set_field_parameters(femurParamsACS)
	FM = fw_femur_measurements.FemurMeasurements(femurModel)
	FM.calcMeasurements()

	o = FM.measurements['head_diameter'].centre
	p1 = FM.shaftAxis.a
	p2 = scipy.array([0,0,0])

	# condyle dorsal vector
	lcondX = femurModel.evaluate_geometric_field_in_elements([10,10],
				[fmd.assemblyElementsNumbers['lateralcondyle']]).T
	mcondX = femurModel.evaluate_geometric_field_in_elements([10,10],
				[fmd.assemblyElementsNumbers['medialcondyle']]).T
	mcDorsal = mcondX[mcondX[:,1].argmin()]
	lcDorsal = lcondX[lcondX[:,1].argmin()]

	alignedParams, T = alignAnatomicFemurOrthoload( femurModel.get_field_parameters().squeeze().T,
													o, p1, p2, mcDorsal, lcDorsal, returnT=True )

	alignedParams = alignedParams.T[:,:,scipy.newaxis]
	return alignedParams, T

def alignFemurMeshParametersAnatomicSingle( g ):
	""" given a femur geometric field, align it geometrically.
	returns the aligned field parameters
	""" 
	d = (10,10)
	head = g.calc_CoM_2D( d, elem=fmd.assemblyElementsNumbers['head'] ) 
	lc = g.calc_CoM_2D( d, elem=fmd.assemblyElementsNumbers['lateralcondyle'] ) 
	mc = g.calc_CoM_2D( d, elem=fmd.assemblyElementsNumbers['medialcondyle'] ) 

	alignedParams, T = alignAnatomicFemur( g.get_field_parameters().squeeze().T, head, mc, lc, returnT=True )
	alignedParams = alignedParams.T[:,:,scipy.newaxis]
	return alignedParams, T

def alignFemurMeshParametersAnatomic( Gs ):
	""" given a list of femur geometric fields, align them geometrically.
	returns the aligned field parameters
	""" 
	alignedParams = []
	d = (10,10)
	for g in Gs:
		head = g.calc_CoM_2D( d, elem=fmd.assemblyElementsNumbers['head'] ) 
		lc = g.calc_CoM_2D( d, elem=fmd.assemblyElementsNumbers['lateralcondyle'] ) 
		mc = g.calc_CoM_2D( d, elem=fmd.assemblyElementsNumbers['medialcondyle'] ) 
		alignedParams.append( alignment_analytic.alignAnatomic( g.get_field_parameters().squeeze().T, head, mc, lc ).T[:,:,scipy.newaxis] )
		
	return alignedParams

#=======================================================#
# pelvis alignment                                      #
#=======================================================#
def createPelvisACSISB(lasis, rasis, lpsis, rpsis):
	"""Calculate the ISB pelvis anatomic coordinate system
	axes: x-anterior, y-superior, z-right
	"""
	oa = (lasis + rasis)/2.0
	op = (lpsis + lpsis)/2.0
	# right
	z = normaliseVector(rasis - lasis)
	# anterior, in plane of op, rasis, lasis
	n1 = normaliseVector(scipy.cross(rasis-op, lasis-op))
	x = normaliseVector(scipy.cross(n1,z))
	# superior
	y = normaliseVector(scipy.cross(z,x))
	return oa, x, y, z

def createPelvisACSAPP(lasis, rasis, lpt, rpt):
	"""Calculate the anterior pelvic plane anatomic
	coordinate system: x-right, y-anterior, z-superior
	"""
	# lasis = scipy.array(lasis)
	# rasis = scipy.array(rasis)
	# lpt = scipy.array(lpt)
	# rpt = scipy.array(rpt)
	
	o = 0.5*(lasis+rasis) 
	pt = 0.5*(lpt + rpt)
	x = normaliseVector(rasis-lasis)
	y = normaliseVector(scipy.cross(x, pt-o))
	# y = normaliseVector(scipy.cross(x, lpt-lasis))
	z = normaliseVector(scipy.cross(x, y))
	return o, x, y, z

def alignAnatomicPelvis( X, lasis, rasis, lpsis, rpsis, returnT=False ):
	
	# oa = ( lasis + rasis )/2.0
	# op = ( lpsis + lpsis )/2.0
	# z = normaliseVector( rasis - lasis )
	# y = normaliseVector( scipy.cross( z, op - rasis ) )
	# x = normaliseVector( scipy.cross( y, z ) )

	o, x, y, z = createPelvisACSISB(lasis, rasis, lpsis, rpsis)
	
	u = scipy.array([ o, o+x, o+y, o+z])		
	ut = scipy.array([[0, 0, 0],\
						[1, 0, 0],\
						[0, 1, 0],\
						[0, 0, 1]])
				
	t = transform3D.directAffine( u, ut )
	if t.shape==(3,4):
		t = scipy.vstack([t, [0,0,0,1]])

	if returnT:
		return transform3D.transformAffine( X, t ), t
	else:
		return transform3D.transformAffine( X, t )

def alignAnatomicPelvisAPP(X, lasis, rasis, lpt, rpt, returnT=False):
	"""
	Align to the Anterior Pelvic Plane (APP) coordinate system commonly 
	used in hip surgery.

	The APP is defined by the LASIS, RPSIS, and the midpoint between left and
	right pubic tubercles. The x axis is parallel to LASIS-RASIS, the y axis is
	normal to the APP, and the z axis is normal to the x and y axes. The
	origin is the midpoint between LASIS and RASIS.
	"""
	# lasis = scipy.array(lasis)
	# rasis = scipy.array(rasis)
	# lpt = scipy.array(lpt)
	# rpt = scipy.array(rpt)
	
	# o = 0.5*(lasis+rasis) 
	# pt = 0.5*(lpt + rpt)
	# x = normaliseVector(rasis-lasis)
	# y = normaliseVector(scipy.cross(x, pt-o))
	# z = normaliseVector(scipy.cross(x, y))

	o, x, y, z = createPelvisACSAPP(lasis, rasis, lpt, rpt)

	u = scipy.array([ o, o+x, o+y, o+z])		
	ut = scipy.array([[0, 0, 0],
					  [1, 0, 0],
					  [0, 1, 0],
					  [0, 0, 1]])
				
	t = transform3D.directAffine(u, ut)
	if t.shape==(3,4):
		t = scipy.vstack([t, [0,0,0,1]])

	if returnT:
		return transform3D.transformAffine( X, t ), t
	else:
		return transform3D.transformAffine( X, t )

def alignAnatomicLH( X, lasis, lpsis, FHC ):
	
	y = normaliseVector( lpsis - FHC )
	x = normaliseVector( scipy.cross( y, lasis - FHC ) )
	z = normaliseVector( scipy.cross( x, y ) )
	
	u = scipy.array([ FHC, FHC+x, FHC+y, FHC+z])		
	ut = scipy.array([[0, 0, 0],\
						[1, 0, 0],\
						[0, 1, 0],\
						[0, 0, 1]])
				
	t = transform3D.directAffine( u, ut )
	if t.shape==(3,4):
		t = scipy.vstack([t, [0,0,0,1]])
	return transform3D.transformAffine( X, t )	

# pelvisLandmarkNodes = {'lasis': 1005,
# 					  'rasis':  465,
# 					  'lpsis':  924,
# 					  'rpsis':  384,
# 					  }

pelvisLandmarkNodes = {'lasis': fw_model_landmarks._pelvisLASISNode,
					   'rasis': fw_model_landmarks._pelvisRASISNode,
					   'lpsis': fw_model_landmarks._pelvisLPSISNode,
					   'rpsis': fw_model_landmarks._pelvisRPSISNode,
					   'lpt': fw_model_landmarks._pelvisLPTNode,
					   'rpt': fw_model_landmarks._pelvisRPTNode,
					  }

LHLandmarkNodes = {'lasis': 466,
					'lpsis': 384,
					}

LHAcetabulumElements = [38,39,40,41,42]

def alignLHMeshParametersAnatomic( Gs ):
	"""
	three landmarks are the CoMs of the three pelvis bones
	"""
	
	alignedParams = []
	d = (10,10)
	for g in Gs:
		nodeCoords = g.get_all_point_positions()
		lasis = nodeCoords[LHLandmarkNodes['lasis']]
		lpsis = nodeCoords[LHLandmarkNodes['lpsis']]
		acetEP = g.evaluate_geometric_field_in_elements( d, LHAcetabulumElements ).T
		FHC, FHRadius = alignment_fitting.fitSphere( acetEP )
		alignedParams.append( alignAnatomicLH( g.get_field_parameters().squeeze().T, lasis, lpsis, FHC ).T[:,:,scipy.newaxis] )
		
	return alignedParams

def alignWholePelvisMeshParametersAnatomicSingle( g ):
	nodeCoords = g.get_all_point_positions()
	lasis = nodeCoords[pelvisLandmarkNodes['lasis']]
	rasis = nodeCoords[pelvisLandmarkNodes['rasis']]
	lpsis = nodeCoords[pelvisLandmarkNodes['lpsis']]
	rpsis = nodeCoords[pelvisLandmarkNodes['rpsis']]
	
	alignedParams, t = alignAnatomicPelvis( g.get_field_parameters().squeeze().T, lasis, rasis, lpsis, rpsis, returnT=True )
	alignedParams = alignedParams.T[:,:,scipy.newaxis]

	return alignedParams, t

def alignWholePelvisMeshParametersAnatomic( Gs ):
	"""
	three landmarks are the CoMs of the three pelvis bones
	"""
	
	alignedParams = []
	for g in Gs:
		alignedParams.append( alignWholePelvisMeshParametersAnatomicSingle( g )[0] )
		
	return alignedParams

def alignWholePelvisMeshParametersAnatomicAPPSingle( g ):
	nodeCoords = g.get_all_point_positions()
	lasis = nodeCoords[pelvisLandmarkNodes['lasis']]
	rasis = nodeCoords[pelvisLandmarkNodes['rasis']]
	lpt = nodeCoords[pelvisLandmarkNodes['lpt']]
	rpt = nodeCoords[pelvisLandmarkNodes['rpt']]
	
	alignedParams, t = alignAnatomicPelvisAPP(g.get_field_parameters().squeeze().T,
											  lasis, rasis, lpt, rpt,
											  returnT=True
											  )
	alignedParams = alignedParams.T[:,:,scipy.newaxis]

	return alignedParams, t

def alignWholePelvisMeshParametersAnatomicAPP( Gs ):
	"""
	three landmarks are the CoMs of the three pelvis bones
	"""
	
	alignedParams = []
	for g in Gs:
		alignedParams.append(alignWholePelvisMeshParametersAnatomicAPPSingle(g)[0])
		
	return alignedParams

def alignPelvisLandmarksPC(gf, pc, landmarks, weights=1.0, GFParamsCallback=None, mw0=1.0, mwn=1.0):
	"""
	landmarks: a list of tuples [(landmark name, landmark coords),...]
	valid landmark names: LASIS, RASIS, LPSIS, RPSIS, Sacral
	"""

	##############
	sourceGF = copy.deepcopy(gf)
	targetLandmarks = []
	ldObjs = []
	for ldName, ldTarg in landmarks:
		targetLandmarks.append(ldTarg)
		evaluator = fw_model_landmarks.makeLandmarkEvaluator('pelvis-'+ldName, sourceGF)
		ldObjs.append(_makeLandmarkObj(ldTarg, evaluator))

	def obj(P):
		P3 = P.reshape((3,-1))
		se = scipy.array([f(P3) for f in ldObjs])
		# print se
		sse = (se*weights).sum()
		return sse

	pcFitter = PCA_fitting.PCFit(pc=pc)
	pcFitter.useFMin = True
	pcFitter.ftol = 1e-3

	P0 = pc.getMean().reshape((3,-1))
	n0 = ldObjs[0](P0)
	x0 = scipy.hstack([targetLandmarks[0] - n0, 0, 0, 0])

	# targetCoM = (targetHC + ((targetMEC+targetLEC)/2.0))/2.0
	# x0 = scipy.hstack([targetCoM - sourceGF.calc_CoM(), 0, 0, 0])
	
	rigidT, rigidP = pcFitter.rigidFit(obj, x0=x0)
	sourceGF.set_field_parameters(rigidP.reshape((3,-1,1)))
	rigidSSE = obj(rigidP)
	if GFParamsCallback is not None:
		GFParamsCallback(rigidP)

	rigidMode0T, rigidMode0P = pcFitter.rigidMode0Fit(obj, mWeight=mw0)
	sourceGF.set_field_parameters(rigidMode0P.reshape((3,-1,1)))
	rigidMode0SSE = obj(rigidMode0P)
	if GFParamsCallback is not None:
		GFParamsCallback(rigidMode0P)
	
	rigidModeNT, rigidModeNP = pcFitter.rigidModeNFit(obj, modes=[1,2], mWeight=mwn)
	sourceGF.set_field_parameters(rigidModeNP.reshape((3,-1,1)))
	rigidModeNSSE = obj(rigidModeNP)
	if GFParamsCallback is not None:
		GFParamsCallback(rigidModeNP)

	return sourceGF, (rigidSSE, rigidMode0SSE, rigidModeNSSE), rigidModeNT 

#========================#
# Tibia fibula alignment #
#========================#
def createTibiaFibulaACSGroodSuntay(MM, LM, MC, LC, side='left'):
	"""Axes: medial, anterior, proximal
	"""
	IC = (MC + LC)/2.0
	IM = (MM + LM)/2.0
	
	z = normaliseVector(IC - IM)
	if side=='right':
		z *= -1.0
	y = normaliseVector(scipy.cross(z, MC-LC))
	x = normaliseVector(scipy.cross(y, z))

	return IC, x, y, z

def createTibiaFibulaACSISB(MM, LM, MC, LC, side='left'):
	"""Axes: x-anterior, y-superior, z-right. Calcaneus CS
	"""
	IC = (MC + LC)/2.0
	IM = (MM + LM)/2.0 # origin
	
	# superiorly, IM to IC
	y = normaliseVector(IC-IM)

	# anteriorly, normal to plane of IM, LC and MC
	x = normaliseVector(scipy.cross(LC-IM, MC-IM))

	# right
	z = normaliseVector(scipy.cross(x,y))
	if side=='right':
		z *= -1.0
		x *= -1.0

	return IM, x, y, z

def alignAnatomicTibiaFibulaGroodSuntay(X, MM, LM, MC, LC, returnT=False):

	# IC = (MC + LC)/2.0
	# IM = (MM + LM)/2.0
	
	# z = normaliseVector(IC - IM)
	# y = normaliseVector(scipy.cross(z, MC-LC))
	# x = normaliseVector(scipy.cross(y, z))

	IC, x, y, z = createTibiaFibulaACSGroodSuntay(MM, LM, MC, LC)
	
	u = scipy.array([IC, IC+x, IC+y, IC+z])		
	ut = scipy.array([[0, 0, 0],\
					  [1, 0, 0],\
					  [0, 1, 0],\
					  [0, 0, 1]])
				
	t = transform3D.directAffine(u, ut)
	if t.shape==(3,4):
		t = scipy.vstack([t, [0,0,0,1]])
	if returnT:
		return transform3D.transformAffine(X, t), t
	else:
		return transform3D.transformAffine(X, t)

def alignTibiaFibulaLandmarksPC(gf, pc, landmarks, weights=1.0, GFParamsCallback=None, mw0=1.0, mwn=1.0):
	"""
	landmarks: a list of tuples [(landmark name, landmark coords),...]
	valid landmark names: LM, MM, TT, kneecentre
	"""

	##############
	sourceGF = copy.deepcopy(gf)
	targetLandmarks = []
	ldObjs = []
	for ldName, ldTarg in landmarks:
		targetLandmarks.append(ldTarg)
		evaluator = fw_model_landmarks.makeLandmarkEvaluator('tibiafibula-'+ldName, sourceGF)
		ldObjs.append(_makeLandmarkObj(ldTarg, evaluator))

	def obj(P):
		P3 = P.reshape((3,-1))
		se = scipy.array([f(P3) for f in ldObjs])
		# print se
		sse = (se*weights).sum()
		return sse

	pcFitter = PCA_fitting.PCFit(pc=pc)
	pcFitter.useFMin = True
	pcFitter.ftol = 1e-3

	P0 = pc.getMean().reshape((3,-1))
	n0 = ldObjs[0](P0)
	x0 = scipy.hstack([targetLandmarks[0] - n0, 0, 0, 0])

	# targetCoM = (targetHC + ((targetMEC+targetLEC)/2.0))/2.0
	# x0 = scipy.hstack([targetCoM - sourceGF.calc_CoM(), 0, 0, 0])
	
	rigidT, rigidP = pcFitter.rigidFit(obj, x0=x0)
	sourceGF.set_field_parameters(rigidP.reshape((3,-1,1)))
	rigidSSE = obj(rigidP)
	if GFParamsCallback is not None:
		GFParamsCallback(rigidP)

	rigidMode0T, rigidMode0P = pcFitter.rigidMode0Fit(obj, mWeight=mw0)
	sourceGF.set_field_parameters(rigidMode0P.reshape((3,-1,1)))
	rigidMode0SSE = obj(rigidMode0P)
	if GFParamsCallback is not None:
		GFParamsCallback(rigidMode0P)
	
	rigidModeNT, rigidModeNP = pcFitter.rigidModeNFit(obj, modes=[1,2], mWeight=mwn)
	sourceGF.set_field_parameters(rigidModeNP.reshape((3,-1,1)))
	rigidModeNSSE = obj(rigidModeNP)
	if GFParamsCallback is not None:
		GFParamsCallback(rigidModeNP)

	return sourceGF, (rigidSSE, rigidMode0SSE, rigidModeNSSE), rigidModeNT 

#===================#
# Patella alignment #
#===================#
def createPatellaACSTest(sup, inf, lat, side='left'):
	"""Axes: x-anterior, y-superior, z-right
	"""
	o = (sup + inf)/2.0
	y = normaliseVector(sup - inf)
	x = normaliseVector(scipy.cross(lat-inf, sup-inf))
	z = normaliseVector(scipy.cross(x,y))
	if side=='right':
		z *= -1.0
	return o, x, y, z