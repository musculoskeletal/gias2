"""
FILE: transform3D.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: 3D geometric transformations

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import scipy
from scipy.linalg import inv
from numpy.linalg import svd, det
from gias2.common import math

def calcAffineMatrix( scale=None, trans=None, shear=None, rot=None ):

	T = scipy.array([[1.0, 0.0, 0.0, 0.0],\
					 [0.0, 1.0, 0.0, 0.0],\
					 [0.0, 0.0, 1.0, 0.0],\
					 [0.0, 0.0, 0.0, 1.0]])
	
	if scale!=None:
		T[0,0] = scale[0]
		T[1,1] = scale[1]
		T[2,2] = scale[2]
	           
	if trans!=None:
		T[:3,3] = trans
		
	if shear!=None:
		T[0,1] = shear[0]	# xy
		T[0,2] = shear[1]	# xz
		T[1,0] = shear[2]	# yx
		T[1,2] = shear[3]	# yz
		T[2,0] = shear[4]	# zx
		T[2,1] = shear[5]	# zy
	
	if rot!=None:        
		Rx = scipy.array( [[1.0, 0.0, 0.0],\
						   [0.0, scipy.cos(rot[0]), -scipy.sin(rot[0])],\
						   [0.0, scipy.sin(rot[0]),  scipy.cos(rot[0])]] )
					 
		Ry = scipy.array( [[ scipy.cos(rot[1]), 0.0, scipy.sin(rot[1])],\
						   [0.0, 1.0, 0.0],\
						   [-scipy.sin(rot[1]), 0.0, scipy.cos(rot[1])]] )
					 
		Rz = scipy.array( [[scipy.cos(rot[2]), -scipy.sin(rot[2]), 0.0],\
						   [scipy.sin(rot[2]),  scipy.cos(rot[2]), 0.0],\
						   [0.0, 0.0, 1.0]] )
		
		#~ print rot
		#~ print Rx 
		#~ print Ry
		#~ print Rz
		#~ print scipy.dot( scipy.dot( Rx,Ry ),Rz )
		T[:3,:3] = scipy.dot( scipy.dot( scipy.dot( Rx,Ry ),Rz ), T[:3,:3] )
		
	return T

def calcRigidAffineMatrix(t, com=None):
	"""
	calculates an affine transformation matrix that applies a
	rigid transform equivalent to that produced by the rigid
	transform functions here using their output translation
	and rotation parameters.
	"""

	if com is None:
		com = [0.0, 0.0, 0.0]

	# rotations are about the CoM of the source data cloud
	# first translate CoM to origin
	T1 = scipy.array([[1.0, 0.0, 0.0, -com[0]],
				   [0.0, 1.0, 0.0, -com[1]],
				   [0.0, 0.0, 1.0, -com[2]],
				   [0.0, 0.0, 0.0, 1.0]
				   ])
	# then rotate and translate
	Rx = scipy.array([[1.0, 0.0, 0.0, 1.0],
				   [0.0, scipy.cos(t[3]), -scipy.sin(t[3]), 1.0],
				   [0.0, scipy.sin(t[3]),  scipy.cos(t[3]), 1.0],
				   [0.0, 0.0, 0.0, 1.0]
				   ])

	Ry = scipy.array([[scipy.cos(t[4]), 0.0, scipy.sin(t[4]), 1.0],
	              [0.0, 1.0, 0.0, 1.0],
	              [-scipy.sin(t[4]), 0.0, scipy.cos(t[4]), 1.0],
	              [0.0, 0.0, 0.0, 1.0]
	              ])

	Rz = scipy.array([[scipy.cos(t[5]), -scipy.sin(t[5]), 0.0, 1.0],
	              [scipy.sin(t[5]), scipy.cos(t[5]), 0.0, 1.0],
	              [0.0, 0.0, 1.0, 1.0],
	              [0.0, 0.0, 0.0, 1.0]
	              ])
	T2 = scipy.dot(scipy.dot(Rx,Ry),Rz)
	T2[0:3,3] = t[:3]	# translation elements

	# then translate back to com
	T3 = scipy.array([[1.0, 0.0, 0.0, com[0]],
	               [0.0, 1.0, 0.0, com[1]],
	               [0.0, 0.0, 1.0, com[2]],
	               [0.0, 0.0, 0.0, 1.0]
	               ])

	# T3 x T2 x T1
	T = scipy.dot(scipy.dot(T3, T2), T1)
	return T

def transformRigid3D( x, t ):
	""" applies a rigid transform to list of points x.
	T = (tx,ty,tz,rx,ry,rz)
	"""
	X = scipy.vstack( (x.T, scipy.ones(x.shape[0]) ) )
	T = scipy.array([[1.0, 0.0, 0.0, t[0]],\
	           [0.0, 1.0, 0.0, t[1]],\
	           [0.0, 0.0, 1.0, t[2]],\
	           [0.0, 0.0, 0.0, 1.0]])
	            
	Rx = scipy.array( [[1.0, 0.0, 0.0],\
				 [0.0, scipy.cos(t[3]), -scipy.sin(t[3])],\
				 [0.0, scipy.sin(t[3]),  scipy.cos(t[3])]] )
				 
	Ry = scipy.array( [[scipy.cos(t[4]), 0.0, scipy.sin(t[4])],\
				 [0.0, 1.0, 0.0],\
				 [-scipy.sin(t[4]), 0.0, scipy.cos(t[4])]] )
				 
	Rz = scipy.array( [[scipy.cos(t[5]), -scipy.sin(t[5]), 0.0],\
				 [scipy.sin(t[5]), scipy.cos(t[5]), 0.0],\
				 [0.0, 0.0, 1.0]] )
	
	T[:3,:3] = scipy.dot( scipy.dot( Rx,Ry ),Rz )
	return scipy.dot( T, X )[:3,:].T

def transformRigid3DAboutCoM( x, t ):
	""" applies a rigid transform to list of points x.
	T = (tx,ty,tz,rx,ry,rz), rotation is about center of mass. Rotates
	then translates.
	"""	
	com = x.mean(0)
	xO = x - com
	xOT = transformRigid3D( xO, t )
	return xOT + com

def transformRigid3DAboutP( x, t, P ):
	""" applies a rigid transform to list of points x.
	T = (tx,ty,tz,rx,ry,rz), rotation is about point P. Rotates
	then translates.
	"""	
	xO = x - P
	xOT = transformRigid3D( xO, t )
	return xOT + P

def transformRigidScale3DAboutCoM( x, t ):
	""" applies a rigid + scale transform to list of points x.
	T = (tx,ty,tz,rx,ry,rz,sx,sy,sz), rotation and scaling is about center of mass. 
	Scales, rotates, then translates.
	"""	
	com = x.mean(0)
	xO = x - com
	xOS = transformScale3D( xO, t[6:])
	xOT = transformRigid3D( xOS, t[:6] )
	return xOT + com

def transformRigidScale3DAboutP( x, t, P ):
	""" applies a rigid + scale transform to list of points x.
	T = (tx,ty,tz,rx,ry,rz,sx,sy,sz), rotation and scaling is about point P. 
	Scales, rotates, then translates.
	"""	
	xO = x - P
	xOS = transformScale3D( xO, t[6:])
	xOT = transformRigid3D( xOS, t[:6] )
	return xOT + P

def transformScale3D( x, S ):
	""" applies scaling to a list of points x. S = (sx,sy,sz)
	"""
	return scipy.multiply( x, S )
	
def transformRigidScale3D( x, t ):
	return transformScale3D( transformRigid3D(x, t[:6]), t[6] )	

def transformAffine( x, t ):
	""" applies affine transform t (shape = (3,4) or (4,4)) to list of points x
	"""
	return scipy.dot( t[:3,:], scipy.vstack( (x.T, scipy.ones(x.shape[0])) ) )[:3,:].T 
	
def transformRotateAboutP( x, r, P ):
	"""
	rotate points x about P by r=[rx,ry,rz]
	"""
	# P to origin
	xO = x - P
	# rotate
	xOR = transformRigid3D( xO, scipy.hstack([[0,0,0],r]) )
	# move to P
	xR = xOR + P
	return xR	

def directAffine(u, ut):
	""" calculate the affine transformation using least squares direct
	method (Kumar)
	"""

	# transformation matrix
	A = scipy.zeros( (4,4) )
	b0 = scipy.zeros( 4 )
	b1 = scipy.zeros( 4 )
	b2 = scipy.zeros( 4 )
	
	for d in range(len(u)):
		A[0,0]=A[0,0]+u[d,0]*u[d,0];
		A[0,1]=A[0,1]+u[d,0]*u[d,1];
		A[0,2]=A[0,2]+u[d,0]*u[d,2];
		A[0,3]=A[0,3]+u[d,0]*1;    
		A[1,0]=A[1,0]+u[d,1]*u[d,0];
		A[1,1]=A[1,1]+u[d,1]*u[d,1];
		A[1,2]=A[1,2]+u[d,1]*u[d,2];
		A[1,3]=A[1,3]+u[d,1]*1;
		A[2,0]=A[2,0]+u[d,2]*u[d,0];
		A[2,1]=A[2,1]+u[d,2]*u[d,1];
		A[2,2]=A[2,2]+u[d,2]*u[d,2];
		A[2,3]=A[2,3]+u[d,2]*1;
		A[3,0]=A[3,0]+1*u[d,0];
		A[3,1]=A[3,1]+1*u[d,1];
		A[3,2]=A[3,2]+1*u[d,2];
		A[3,3]=A[3,3]+1*1;
			
		b0[0]=b0[0]+u[d,0]*ut[d,0];
		b0[1]=b0[1]+u[d,1]*ut[d,0];
		b0[2]=b0[2]+u[d,2]*ut[d,0];
		b0[3]=b0[3]+     1*ut[d,0];
		
		b1[0]=b1[0]+u[d,0]*ut[d,1];
		b1[1]=b1[1]+u[d,1]*ut[d,1];
		b1[2]=b1[2]+u[d,2]*ut[d,1];
		b1[3]=b1[3]+     1*ut[d,1];
		
		b2[0]=b2[0]+u[d,0]*ut[d,2];
		b2[1]=b2[1]+u[d,1]*ut[d,2];
		b2[2]=b2[2]+u[d,2]*ut[d,2];
		b2[3]=b2[3]+     1*ut[d,2];        

	B = inv(A)
	
	t0= scipy.dot( B, b0 )
	t1 = scipy.dot( B, b1 )
	t2 = scipy.dot( B, b2 )
	
	T = scipy.vstack( [ t0, t1, t2 ] )

	return T

def transformRotateAboutAxisOld(x, theta, p0, p1, p2=None, p3=None):

	# create orthogonal vectors with p0 and p1
	p0 = scipy.array(p0)
	v1 = math.norm(scipy.subtract(p1,p0))
	v2 = math.norm(scipy.cross(v1, v1*[2.0,0,0]))
	v3 = math.norm(scipy.cross(v1, v2))

	# if p2 is None:
	# 	p2 = p0 + math.norm(scipy.cross(p1-p0, scipy.multiply(p1,[2.0,0,0])-p0))
	# 	p3 = p0 + math.norm(scipy.cross(p1-p0, p2-p0))

	# transform v to global x - T0
	cs_global = scipy.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]], dtype=float)
	cs_local = scipy.array([p0, p0+v1, p0+v2, p0+v3], dtype=float)
	T0 = directAffine(cs_local, cs_global)
	x0 = transformAffine(x, T0)

	# rotate about x with theta
	x1 = transformRigid3DAboutP(x0, [0,0,0,theta,0,0], (0,0,0))

	# apply inverse of T0
	x2 = transformAffine(x1, inv(scipy.vstack([T0,[0,0,0,1]])))

	return x2

def transformRotateAboutAxis(x, theta, p0, p1, retmat=False):
	"""
	Rotate point x by angle theta about an axis defined by 2 points p0, p1.
	http://paulbourke.net/geometry/rotate/
	"""

	# translate to axis is at origin
	p = (x - p0).T
	# rotation unit vector
	n = math.norm(p1 - p0)

	# matrix common factors
	c = scipy.cos(theta)
	t = 1.0-scipy.cos(theta)
	s = scipy.sin(theta)

	# matrix
	d11 = t*n[0]*n[0] + c
	d12 = t*n[0]*n[1] - s*n[2]
	d13 = t*n[0]*n[2] + s*n[1]
	d21 = t*n[0]*n[1] + s*n[2]
	d22 = t*n[1]*n[1] + c
	d23 = t*n[1]*n[2] - s*n[0]
	d31 = t*n[0]*n[2] - s*n[1]
	d32 = t*n[1]*n[2] + s*n[0]
	d33 = t*n[2]*n[2] + c

	# matrix M*|p|
	q = scipy.array([d11*p[0] + d12*p[1] + d13*p[2],
				 	 d21*p[0] + d22*p[1] + d23*p[2],
				 	 d31*p[0] + d32*p[1] + d33*p[2]])

	if retmat:
		M1 = scipy.array([
				[1,0,0,-p0[0]],
				[0,1,0,-p0[1]],
				[0,0,1,-p0[2]],
				[0,0,0,1],
			])
		M2 = scipy.array([
				[d11,d12,d13,p0[0]],
				[d21,d22,d23,p0[1]],
				[d31,d32,d33,p0[2]],
				[0,0,0,1],
			])
		M = M2.dot(M1)
		return q.T + p0, M
	else:
		return q.T + p0


def transformRotateAboutCartCS(x, r, o, v1, v2, v3):
	"""
	rotate about an arbitrary cartesian coordinate system.

	Inputs:
	x [array]: nx3 array of points to transform
	r [list]: a list of rotation angles about the CS
	o [array]: coordinates of the CS origin
	v1, v2, v3 [arrays]: 3 orthogonal unit vectors defining the axes of the CS
	""" 

	# transform v to global x - T0
	cs_global = scipy.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]], dtype=float)
	cs_local = scipy.array([o, o+v1, o+v2, o+v3], dtype=float)
	T0 = directAffine(cs_local, cs_global)
	x0 = transformAffine(x, T0)

	# rotate about x with theta
	x1 = transformRigid3DAboutP(x0, [0,0,0,r[0],r[1],r[2]], (0,0,0))

	# apply inverse of T0
	x2 = transformAffine(x1, inv(scipy.vstack([T0,[0,0,0,1]])))

	return x2

def calcAffineMatrixSVD(A, B):
	"""Calculate rigid transformation between two list of corresponding
	points using SVD
	"""
	assert len(A) == len(B)

	N = A.shape[0]; # total points

	centroid_A = scipy.mean(A, axis=0)
	centroid_B = scipy.mean(B, axis=0)

	# centre the points
	AA = A - scipy.tile(centroid_A, (N, 1))
	BB = B - scipy.tile(centroid_B, (N, 1))

	# dot is matrix multiplication for array
	H = scipy.dot(scipy.transpose(AA), BB)

	U, S, Vt = svd(H)

	R = Vt.T * U.T

	# special reflection case
	if det(R) < 0:
	   # print "Reflection detected"
	   Vt[2,:] *= -1
	   R = Vt.T * U.T

	t = scipy.dot(-R, centroid_A.T) + centroid_B.T

	T = scipy.eye(4)
	T[:3,:3] = R
	T[:3,3] = t.squeeze()

	return T