"""
FILE: PCA_fitting.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: functions and classes for mesh/data fitting using principal
components.

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import scipy
from scipy.optimize import leastsq, fmin
from gias2.common import transform3D
from gias2.registration import alignment_fitting
import sys

def rigidObj( x, p, func, funcArgs=() ):
	""" x = [tx,ty,tz,rx,ry,rz]
	"""
	# print x
	# apply transform to parameters
	#~ pT = transform3D.transformRigid3D( p, x )
	pT = transform3D.transformRigid3DAboutCoM( p, x )
	# calc error
	#~ pdb.set_trace()
	return func( pT.T.ravel(), *funcArgs )
	
def rigidScaleObj( x, p, func, funcArgs=() ):
	""" x = [tx,ty,tz,rx,ry,rz,s]
	"""
	# apply transform to parameters
	pT = transform3D.transformRigidScale3DAboutCoM( p, x )
	# calc error
	return func( pT.T.ravel(), *funcArgs )

def rigidMode0Obj( x, func, pca, mWeight, funcArgs=() ):
	# reconstruct shape using 1st mode weight
	p = pca.reconstruct( pca.getWeightsBySD( [0,], [x[6],] ), [0,] )
	# apply rigid transform
	#~ p = transform3D.transformRigid3D( p.reshape((3,-1)).T, x[:6] )
	p = transform3D.transformRigid3DAboutCoM( p.reshape((3,-1)).T, x[:6] )
	# calc error
	return func( p.T.ravel(), *funcArgs ) + mahalanobis( [x[6],] )*mWeight

def rigidModeNObj( x, func, pca, modes, mWeight, funcArgs ):
	# reconstruct shape using modes
	p = pca.reconstruct( pca.getWeightsBySD( modes, x[6:] ), modes )
	# reconstruct rigid transform
	#~ p = transform3D.transformRigid3D( p.reshape((3,-1)).T, x[:6] )
	p = transform3D.transformRigid3DAboutCoM( p.reshape((3,-1)).T, x[:6] )
	# calc error
	funcErr = func( p.T.ravel(), *funcArgs )
	MErr = mahalanobis( x[6:] )*mWeight
	Err = funcErr + MErr
	
	#~ sys.stdout.write( '\rfunc rms:'+str(scipy.sqrt(funcErr.mean())) )
	# sys.stdout.write( '\rfunc rms:'+str(scipy.sqrt(Err.mean())) )
	# sys.stdout.flush()
	
	return funcErr + MErr
	#~ return func( p.T.ravel(), *funcArgs ) + mahalanobis( x[6:] )*mWeight

def rigidModeNRotateAboutCoMObj( x, func, pca, modes, mWeight, funcArgs ):

	# print 'rigidModeNRotateAboutCoMObj x:', x 

	# reconstruct shape using modes
	p = pca.reconstruct( pca.getWeightsBySD( modes, x[6:] ), modes )
	# reconstruct rigid transform
	p = transform3D.transformRigid3DAboutCoM( p.reshape((3,-1)).T, x[:6] )
	# calc error
	funcErr = func( p.T.ravel(), *funcArgs )
	MErr = mahalanobis( x[6:] )*mWeight
	Err = funcErr + MErr
	
	#~ sys.stdout.write( '\rfunc rms:'+str(scipy.sqrt(funcErr.mean())) )
	sys.stdout.write( '\rfunc rms:'+str(scipy.sqrt(Err.mean())) )
	sys.stdout.flush()
	
	return funcErr + MErr
	#~ return func( p.T.ravel(), *funcArgs ) + mahalanobis( x[6:] )*mWeight

def rigidScaleModeNObj( x, func, pca, modes, mWeight, funcArgs ):
	# reconstruct shape using modes
	p = pca.reconstruct( pca.getWeightsBySD( modes, x[7:] ), modes )
	# reconstruct rigid transform
	p = transform3D.transformRigidScale3DAboutCoM( p.reshape((3,-1)).T, x[:7] )
	# calc error
	funcErr = func( p.T.ravel(), *funcArgs )
	MErr = mahalanobis( x[7:] )*mWeight
	Err = funcErr + MErr
	
	#~ sys.stdout.write( '\rfunc rms:'+str(scipy.sqrt(funcErr.mean())) )
	sys.stdout.write( '\rfunc rms:'+str(scipy.sqrt(Err.mean())) )
	sys.stdout.flush()
	
	return funcErr + MErr

def modeNObj( x, func, pca, modes, rigidT, funcArgs ):
	# reconstruct shape using modes
	p = pca.reconstruct( pca.getWeightsBySD( modes, x ), modes )
	# reconstruct rigid transform
	#~ p = transform3D.transformRigid3D( p.reshape((3,-1)).T, rigidT )
	p = transform3D.transformRigid3DAboutCoM( p.reshape((3,-1)).T, rigidT )
	# calc error
	return func( p.T.ravel(), *funcArgs ) + mahalanobis( x[6:] )
	
def mahalanobis( x ):
	return scipy.sqrt( scipy.multiply(x,x).sum() ) 
	
class PCFit( object ):
	
	xtol = 1e-6
	ftol = 1e-6
	maxiter = None
	maxfev = None
	useFMin = False

	def __init__(self, pc=None):
		self.pc = pc

		self.rigidOpt = None
		self.rigidScaleOut = None
		self.rigidMode0Opt = None
		self.rigidModeNOpt = None
		self.rigidScaleModeNOpt = None
	
	def setPC( self, pc ):
		self.pc = pc
	
	def rigidFit( self, func, x0=None, funcArgs=(), p0=None ):
		if x0==None:
			x0=scipy.array([0.0,0.0,0.0,0.0,0.0,0.0])
		if p0==None:
			p0 = self.pc.getMean().reshape( (3,-1) ).T
				
		x0 = scipy.array(x0)
				
		if self.useFMin:
			self.rigidOpt = fmin(rigidObj, x0,
								 args=(p0, func, funcArgs),
								 xtol=self.xtol,
								 ftol=self.ftol,
								 maxiter=self.maxiter,
								 maxfun=self.maxfev,
								 )
		else:
			self.rigidOpt = leastsq(rigidObj, x0,
									args=(p0, func, funcArgs),
									xtol=self.xtol,
									ftol=self.ftol,
									epsfcn=1e-5
									)[0]

		#~ pOpt = transform3D.transformRigid3D( p0, self.rigidOpt ).T
		pOpt = transform3D.transformRigid3DAboutCoM( p0, self.rigidOpt ).T
		return self.rigidOpt, pOpt.ravel()
	
	def rigidScaleFit( self, func, x0=None, funcArgs=(), p0=None ):
		if x0==None:
			x0=[0.0,0.0,0.0,0.0,0.0,0.0,1.0]
		if p0==None:
			p0 = self.pc.getMean().reshape( (3,-1) ).T	
			
		x0 = scipy.array(x0)
		p0 = self.pc.getMean().reshape( (3,-1) ).T

		if self.useFMin:
			self.rigidScaleOpt = fmin(rigidScaleObj, x0,
								 args=(p0, func, funcArgs),
								 xtol=self.xtol,
								 ftol=self.ftol,
								 maxiter=self.maxiter,
								 maxfun=self.maxfev,
								 )
		else:
			self.rigidScaleOpt = leastsq(rigidScaleObj, x0,
										 args=(p0, func, funcArgs),
										 xtol=self.xtol,
										 ftol=self.ftol
										 )[0]
		
		pOpt = transform3D.transformRigidScale3DAboutCoM( p0, self.rigidScaleOpt ).T
		return self.rigidScaleOpt, pOpt.ravel()
	
	def rigidMode0Fit( self, func, x0=None, mWeight=0.0, funcArgs=() ):
		if x0==None:
			x0 = scipy.hstack( (self.rigidOpt, 0.0) )
		elif len(x0)!=7:
			raise ValueError, 'x0 must be 7 long'

		if self.useFMin:
			self.rigidMode0Opt = fmin(rigidMode0Obj, x0,
								 	  args=(func, self.pc, mWeight, funcArgs),
								 	  xtol=self.xtol,
								 	  ftol=self.ftol,
								 	  maxiter=self.maxiter,
								 	  maxfun=self.maxfev,
								 	  )
		else:	
			self.rigidMode0Opt = leastsq(rigidMode0Obj, x0,
										 args=(func, self.pc, mWeight, funcArgs),
										 xtol=self.xtol,
										 ftol=self.ftol
										 )[0]
		
		pOpt = self.pc.reconstruct( self.pc.getWeightsBySD( [0,], [self.rigidMode0Opt[6],] ), [0,] )
		#~ pOpt = transform3D.transformRigid3D( pOpt.reshape((3,-1)).T, self.rigidMode0Opt[:6] )
		pOpt = transform3D.transformRigid3DAboutCoM( pOpt.reshape((3,-1)).T, self.rigidMode0Opt[:6] )
		return self.rigidMode0Opt, pOpt.T.ravel()
		
	def rigidModeNFit( self, func, modes=None, x0=None, mWeight=0.0, maxfev=None, funcArgs=() ):
		# fit modes beyond mode 0
		if x0==None:
			x0 = scipy.hstack( ( self.rigidMode0Opt, scipy.zeros( len(modes) ) ) )
		if modes==None:
			modes = [1,2,3]
			
		modes = scipy.hstack( (0,modes) )
		# print '\nx0:', x0
		# print 'mWeight:', mWeight

		if self.useFMin:
			if maxfev is None:
				maxfev = self.maxfev
			self.rigidModeNOpt = fmin(rigidModeNObj, x0,
								 	  args=(func, self.pc, modes, mWeight, funcArgs),
								 	  xtol=self.xtol,
								 	  ftol=self.ftol,
								 	  maxfun=maxfev,
								 	  maxiter=self.maxiter,
								 	  )
		else:
			if maxfev is None:
				self.rigidModeNOpt = leastsq(rigidModeNObj, x0,
										 args=(func, self.pc, modes, mWeight, funcArgs),
										 xtol=self.xtol,
										 ftol=self.ftol,
										 )[0]
			else:
				self.rigidModeNOpt = leastsq(rigidModeNObj, x0,
											 args=(func, self.pc, modes, mWeight, funcArgs),
											 xtol=self.xtol,
											 ftol=self.ftol,
											 maxfev=maxfev
											 )[0]
		
		pOpt = self.pc.reconstruct( self.pc.getWeightsBySD( modes, self.rigidModeNOpt[6:] ), modes )
		pOpt = transform3D.transformRigid3DAboutCoM( pOpt.reshape((3,-1)).T, self.rigidModeNOpt[:6] )
		# print '\nxOpt:', self.rigidModeNOpt
		
		return self.rigidModeNOpt, pOpt.T.ravel()
		
	def rigidModeNRotateAboutCoMFit( self, func, modes=None, x0=None, mWeight=0.0, maxfev=None, funcArgs=() ):
		# fit modes beyond mode 0
		if x0==None:
			x0 = scipy.hstack( ( self.rigidMode0Opt, scipy.zeros( len(modes) ) ) )
		if modes==None:
			modes = [1,2,3]
			
		modes = scipy.hstack( (0,modes) )
		# print '\nrigidModeNRotateAboutCoMFit x0:', x0
		# print 'mWeight:', mWeight
		if self.useFMin:
			if maxfev is None:
				maxfev = self.maxfev
			self.rigidModeNOpt = fmin(rigidModeNRotateAboutCoMObj, x0,
								 	  args=(func, self.pc, modes, mWeight, funcArgs),
								 	  xtol=self.xtol,
								 	  ftol=self.ftol,
								 	  maxiter=self.maxiter,
								 	  maxfun=maxfev,
								 	  )
		else:
			if maxfev is None:
				self.rigidModeNOpt = leastsq(rigidModeNRotateAboutCoMObj, x0,
											 args=(func, self.pc, modes, mWeight, funcArgs),
											 xtol=self.xtol,
											 ftol=self.ftol,
											 )[0]
			else:
				self.rigidModeNOpt = leastsq(rigidModeNRotateAboutCoMObj, x0,
											 args=(func, self.pc, modes, mWeight, funcArgs),
											 xtol=self.xtol,
											 ftol=self.ftol,
											 maxfev=maxfev
											 )[0]
		
		pOpt = self.pc.reconstruct( self.pc.getWeightsBySD( modes, self.rigidModeNOpt[6:] ), modes )
		pOpt = transform3D.transformRigid3DAboutCoM( pOpt.reshape((3,-1)).T, self.rigidModeNOpt[:6] )
		# print '\nxOpt:', self.rigidModeNOpt
		
		return self.rigidModeNOpt, pOpt.T.ravel()

	def rigidScaleModeNFit( self, func, modes=None, x0=None, mWeight=0.0, maxfev=None, funcArgs=() ):
		# fit modes
		if x0==None:
			x0 = scipy.hstack( ( self.rigidScaleOpt, scipy.zeros( len(modes) ) ) )
		if modes==None:
			modes = [0,1,2]
			
		if self.useFMin:
			if maxfev is None:
				maxfev = self.maxfev
			self.rigidScaleModeNOpt = fmin(rigidScaleModeNObj, x0,
								 	  	   args=(func, self.pc, modes, mWeight, funcArgs),
								 	  	   xtol=self.xtol,
								 	  	   ftol=self.ftol,
								 	  	   maxfun=maxfev,
								 	  	   maxiter=self.maxiter,
								 	  	   )
		else:
			if maxfev is None:
				self.rigidScaleModeNOpt = leastsq(rigidScaleModeNObj, x0,
												  args=(func, self.pc, modes, mWeight, funcArgs),
												  xtol=self.xtol,
												  ftol=self.ftol,
												  )[0]
			else:
				self.rigidScaleModeNOpt = leastsq(rigidScaleModeNObj, x0,
												  args=(func, self.pc, modes, mWeight, funcArgs),
												  xtol=self.xtol,
												  ftol=self.ftol,
												  maxfev=maxfev
												  )[0]
		
		pOpt = self.pc.reconstruct( self.pc.getWeightsBySD( modes, self.rigidScaleModeNOpt[7:] ), modes )
		pOpt = transform3D.transformRigidScale3DAboutCoM( pOpt.reshape((3,-1)).T, self.rigidScaleModeNOpt[:7] )
		return self.rigidScaleModeNOpt, pOpt.T.ravel()

	def modeNFit( self, func, modes=None, x0=None, funcArgs=() ):
		""" fit mode weights only
		"""
		if x0==None:
			x0 = scipy.hstack( (self.rigidMode0Opt[-1], scipy.zeros( len(modes) ) ) )
		if modes==None:
			modes = [1,2,3]
			
		modes = scipy.hstack( (0,modes) )
		rigidT = self.rigidMode0Opt[:6]
		
		if self.useFMin:
			if maxfev is None:
				maxfev = self.maxfev
			self.modeNOpt = fmin(modeNObj, x0,
					 	  	     args=(func, self.pc, modes, mWeight, funcArgs),
					 	  	     xtol=self.xtol,
					 	  	     ftol=self.ftol,
					 	  	     maxfun=maxfev,
					 	  	     maxiter=self.maxiter,
					 	  	     )
		else:
			self.modeNOpt = leastsq(modeNObj, x0,
									args=(func, self.pc, modes, rigidT, funcArgs),
									xtol=self.xtol,
									ftol=self.ftol
									)[0]
		
		pOpt = self.pc.reconstruct( self.pc.getWeightsBySD( modes, self.modeNOpt ), modes )
		#~ pOpt = transform3D.transformRigid3D( pOpt.reshape((3,-1)).T, rigidT[:6] )
		pOpt = transform3D.transformRigid3DAboutCoM( pOpt.reshape((3,-1)).T, rigidT[:6] )
		return self.modeNOpt, pOpt.T.ravel()

def project3DPointsToSSM(data, SSM, projectModes, projectVariables=None, 
						 landmarkIs=None, initRotation=None, doScale=False, 
						 verbose=False, retT=False):
	
	# rigid align data to SSM mean data
	if initRotation==None:
		initRotation = scipy.array([0.0,0.0,0.0])

	# landmarkIs = scipy.array(landmarkIs, dtype=int)

	# fit has to be done on a subset of recon data is projection is on a subset of variables
	if projectVariables!=None:
		# meanData = SSM.getMean()[projectVariables].reshape((3,-1)).T
		meanData = SSM.getMean().reshape((3,-1)).T[landmarkIs,:]
	else:
		meanData = SSM.getMean().reshape((3,-1)).T
	
	print 'data shape', data.shape
	print 'meandata shape', meanData.shape
	if landmarkIs is not None:
		print 'landmarkIs shape', landmarkIs.shape

	if doScale:
		align1X0 = scipy.hstack([ meanData.mean(0)-data.mean(0), initRotation, 1.0 ])
		T1, dataT = alignment_fitting.fitRigidScale( data, meanData, align1X0, verbose=verbose )
	else:
		align1X0 = scipy.hstack([ meanData.mean(0)-data.mean(0), initRotation ])
		# T1, dataT = alignment_fitting.fitRigid( data, meanData, align1X0, verbose=verbose )
		T1, dataT = alignment_fitting.fitRigid( data, meanData, align1X0, xtol=1e-6, epsfcn=1e-9, verbose=verbose )

	# project data onto principal components
	pcWeights = SSM.project( (dataT.T-meanData.T).ravel(), modes=projectModes, variables=projectVariables )
	# reconstruction data using projected weights
	reconData = SSM.reconstruct(pcWeights, projectModes).reshape((3,-1)).T

	# fit has to be done on a subset of recon data is projection is on a subset of variables
	if projectVariables!=None:
		# reconDataFit = reconData.ravel()[projectVariables].reshape((3,-1)).T
		reconDataFit = reconData.T.ravel().reshape((3,-1)).T[landmarkIs,:]
	else:
		reconDataFit = reconData
	
	# inverse transform back to image
	if doScale:
		align2X0 = scipy.hstack([ reconDataFit.mean(0)-data.mean(0), -initRotation, 1.0 ])
		T2, reconDataFitT = alignment_fitting.fitRigidScale( reconDataFit, data, align2X0, verbose=verbose )
	else:
		align2X0 = scipy.hstack([ reconDataFit.mean(0)-data.mean(0), -initRotation ])
		# T2, reconDataT = alignment_fitting.fitRigid( reconData, data, align2X0, verbose=verbose )
		T2, reconDataFitT = alignment_fitting.fitRigid( reconDataFit, data, align2X0, xtol=1e-6, epsfcn=1e-9, verbose=verbose )

	# transform full recon data to data space if data is a subset of variables
	if projectVariables!=None:
		reconDataT = transform3D.transformRigid3DAboutCoM(reconData, T2)
	else:
		reconDataT = reconDataFitT

	if retT:
		return pcWeights, reconDataT, dataT, reconData, T2
	else:
		return pcWeights, reconDataT, dataT, reconData

def fitSSMTo3DPoints(data, SSM, fitModes, fitPointIndices=None, mWeight=0.0,
					 initRotation=None, doScale=False, 
					 landmarkTargets=None, landmarkEvaluator=None, landmarkWeights=None,
					 verbose=False):
	
	print 'fitting SSM to points'
	# rigid align data to SSM mean data
	if initRotation==None:
		initRotation = scipy.array([0.0,0.0,0.0])
	else:
		initRotation = scipy.array(initRotation)

	# fit has to be done on a subset of recon data is projection is on a subset of variables
	if fitPointIndices!=None:
		fitPointIndices = scipy.array(fitPointIndices, dtype=int)
		meanData = SSM.getMean().reshape((3,-1)).T[fitPointIndices,:]
	else:
		meanData = SSM.getMean().reshape((3,-1)).T
	
	# print 'data shape:', data.shape
	# print 'meandata shape:', meanData.shape
	if fitPointIndices==None:
		print 'fitPointIndices shape: None'
	else:
		print 'fitPointIndices shape:', fitPointIndices.shape

	if doScale:
		align1X0 = scipy.hstack([ meanData.mean(0)-data.mean(0), initRotation, 1.0 ])
		T1, dataT = alignment_fitting.fitRigidScale( data, meanData, align1X0, verbose=verbose )
		if landmarkTargets!=None:
			landmarkTargetsT = transform3D.transformRigidScale3DAboutP(landmarkTargets, T1,
																	   data.mean(0)
																	   )
	else:
		align1X0 = scipy.hstack([ meanData.mean(0)-data.mean(0), initRotation ])
		# T1, dataT = alignment_fitting.fitRigid( data, meanData, align1X0, verbose=verbose )
		T1, dataT = alignment_fitting.fitRigid( data, meanData, align1X0, xtol=1e-6, epsfcn=1e-9, verbose=verbose )
		if landmarkTargets!=None:
			landmarkTargetsT = transform3D.transformRigid3DAboutP(landmarkTargets, T1,
																  data.mean(0)
																  )

	def _obj(X):
		recon = SSM.reconstruct( SSM.getWeightsBySD( fitModes, X[6:] ), fitModes )
		# reconstruct rigid transform
		reconData = transform3D.transformRigid3DAboutCoM( recon.reshape((3,-1)).T, X[:6] )

		# print reconData[466]

		# calc error
		if fitPointIndices!=None:
			E = ((reconData[fitPointIndices,:] - dataT)**2.0).sum(1) + mahalanobis( X[6:] )*mWeight
		else:
			E = ((reconData - dataT)**2.0).sum(1) + mahalanobis( X[6:] )*mWeight

		sys.stdout.write( '\robj rms:'+str(scipy.sqrt(E.mean())) )
		sys.stdout.flush()

		return E

	def _objLandmarks(X):
		recon = SSM.reconstruct( SSM.getWeightsBySD( fitModes, X[6:] ), fitModes )
		# reconstruct rigid transform
		reconData = transform3D.transformRigid3DAboutCoM( recon.reshape((3,-1)).T, X[:6] )
		# reconLandmarks = transform3D.transformRigid3DAboutP( landmarkEvaluator(recon), X[:6], recon.reshape((3,-1)).T.mean(0) )
		reconLandmarks = landmarkEvaluator(reconData.T.ravel())

		# print reconLandmarks[0]
		# print reconData[466]
		# print landmarkTargets[0]

		# calc error
		if fitPointIndices!=None:
			EData = ((reconData[fitPointIndices,:] - dataT)**2.0).sum(1) + mahalanobis( X[6:] )*mWeight
		else:
			EData = ((reconData - dataT)**2.0).sum(1) + mahalanobis( X[6:] )*mWeight

		ELandmarks = ((landmarkTargetsT - reconLandmarks)**2.0).sum(1)*landmarkWeights
		E = scipy.hstack([EData, ELandmarks])

		sys.stdout.write('\rPC fit rmse: %6.3f (data: %6.3f) (landmarks: %6.3f)'%\
						 (scipy.sqrt(E.mean()), scipy.sqrt(EData.mean()),scipy.sqrt(ELandmarks.mean()))
						)
		sys.stdout.flush()

		return E
	
	
	# PC Fit
	x0 = scipy.hstack([[0,0,0], initRotation, scipy.zeros(len(fitModes), dtype=float)])

	if landmarkTargets==None:
		print 'using non-landmark obj func'
		xOpt = leastsq( _obj, x0, xtol=1e-6 )[0]
	else:
		print 'using landmark obj func'
		if verbose:
			print 'landmarks'
			print landmarkTargets
		xOpt = leastsq( _objLandmarks, x0, xtol=1e-6 )[0]
	print ' '
		
	reconOpt = SSM.reconstruct( SSM.getWeightsBySD( fitModes, xOpt[6:] ), fitModes )
	reconDataOpt = transform3D.transformRigid3DAboutCoM( reconOpt.reshape((3,-1)).T, xOpt[:6] )

	# fit has to be done on a subset of recon data is projection is on a subset of variables
	if fitPointIndices!=None:
		# reconDataFit = reconData.ravel()[projectVariables].reshape((3,-1)).T
		reconDataFit = scipy.array(reconDataOpt.T.ravel().reshape((3,-1)).T[fitPointIndices,:])
	else:
		reconDataFit = scipy.array(reconDataOpt)
	
	# inverse transform back to image
	if doScale:
		align2X0 = scipy.hstack([ reconDataFit.mean(0)-data.mean(0), -initRotation, 1.0 ])
		T2, reconDataFitT = alignment_fitting.fitRigidScale( reconDataFit, data, align2X0, verbose=verbose )
	else:
		align2X0 = scipy.hstack([ reconDataFit.mean(0)-data.mean(0), -initRotation ])
		# T2, reconDataT = alignment_fitting.fitRigid( reconData, data, align2X0, verbose=verbose )
		T2, reconDataFitT = alignment_fitting.fitRigid( reconDataFit, data, align2X0, xtol=1e-6, epsfcn=1e-9, verbose=verbose )

	# transform full recon data to data space if data is a subset of variables
	if fitPointIndices!=None:
		reconDataOptT = transform3D.transformRigid3DAboutCoM(reconDataOpt, T2)
	else:
		reconDataOptT = reconDataFitT

	xOpt[:6] += T2

	return xOpt, reconDataOptT, dataT, reconDataOpt