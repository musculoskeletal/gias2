"""
FILE: alignment_fitting.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: Iterative alignment of points.

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import scipy
from scipy.spatial import cKDTree
from scipy.optimize import leastsq, fmin
from scipy.linalg import lstsq
from gias2.common import transform3D

def _sampleData(data, N):
    """
    Pick N evenly spaced points from data
    """

    if N<1:
        raise ValueError('N must be > 1')
    elif N>len(data):
        return data
    else:
        i = scipy.linspace(0,len(data)-1,N).astype(int)
        return data[i,:]


#======================================================================#
# correspondent fitting data fitting                                   #
#======================================================================#
def fitAffine( data, target, xtol=1e-5, maxfev=0, sample=None, verbose=0, outputErrors=0 ):
    if len(data)!=len(target):
        raise ValueError('data and target points must have same number of points')

    rms0 = scipy.sqrt(((data - target)**2.0 ).sum(1).mean())

    if sample is not None:
        D = _sampleData(data,sample)
        T = _sampleData(target,sample)
    else:
        D = data
        T = target
    
    t = transform3D.directAffine(D, T)
    dataFitted = transform3D.transformAffine( data, t )
    rmsOpt = scipy.sqrt(((dataFitted - target)**2.0 ).sum(1).mean())
    if verbose:
        print('initial RMS:', rms0)
        print('final RMS:', rmsOpt)

    if outputErrors:
        return t, dataFitted, (rms0, rmsOpt)
    else:
        return t, dataFitted    

def fitTranslation( data, target, xtol=1e-5, maxfev=0, sample=None, verbose=0, outputErrors=0 ):
    """ fits for tx,ty for transforms points in data to points
    in target. Points in data and target are assumed to correspond by
    order
    """ 

    if sample is not None:
        D = _sampleData(data,sample)
        T = _sampleData(target,sample)
    else:
        D = data
        T = target
    
    def obj( x ):
        DT = D + x
        d = ( ( DT - T )**2.0 ).sum(1)
        return d
    
    t0 = target.mean(0) - data.mean(0)
    # t0 = scipy.array([ 0.0, 0.0, 0.0 ])
    
    rms0 = scipy.sqrt( obj( t0 ).mean() )
    if verbose:
        print('initial RMS:', rms0)
        
    xOpt = leastsq( obj, t0, xtol=xtol, maxfev=maxfev )[0]
    
    rmsOpt = scipy.sqrt( obj(xOpt).mean() )
    if verbose:
        print('final RMS:', rmsOpt)
    
    dataFitted = data + xOpt
    if outputErrors:
        return xOpt, dataFitted, (rms0, rmsOpt)
    else:
        return xOpt, dataFitted 
    
def fitRigid(data, target, t0=None, xtol=1e-3, rotcentre=None, maxfev=None,
    maxfun=None, sample=None, verbose=0, epsfcn=0, outputErrors=0 ):
    """ fits for tx,ty,tz,rx,ry,rz to transform points in data to points
    in target. Points in data and target are assumed to correspond by
    order
    """ 

    if sample is not None:
        D = _sampleData(data,sample)
        T = _sampleData(target,sample)
    else:
        D = data
        T = target

    if t0 is None:
        t0 = scipy.array([0.0,0.0,0.0,0.0,0.0,0.0])
    else:
        t0 = scipy.array(t0)
        
    if rotcentre is None:
        rotcentre = D.mean(0)

    if data.shape[0]>=t0.shape[0]:
        def obj( x ):
            # DT = transform3D.transformRigid3DAboutCoM( D, x )
            DT = transform3D.transformRigid3DAboutP( D, x, rotcentre )
            d = ( ( DT - T )**2.0 ).sum(1)
            return d
    else:
        def obj( x ):
            # DT = transform3D.transformRigid3DAboutCoM( D, x )
            DT = transform3D.transformRigid3DAboutP( D, x, rotcentre )
            d = ( ( DT - T )**2.0 ).sum(1)
            return d.sum()
        
    t0 = scipy.array(t0)
    rms0 = scipy.sqrt( obj( t0 ).mean() )
    if verbose:
        print('initial RMS:', rms0)
        
    if data.shape[0]>=t0.shape[0]:
        if maxfev is None:
            maxfev = 0
        xOpt = leastsq( obj, t0, xtol=xtol, maxfev=maxfev, epsfcn=epsfcn )[0]
    else:
        xOpt = fmin( obj, t0, xtol=xtol, maxiter=maxfev, maxfun=maxfun, disp=verbose )
    
    rmsOpt = scipy.sqrt( obj(xOpt).mean() )
    if verbose:
        print('final RMS:', rmsOpt)
    
    # dataFitted = transform3D.transformRigid3DAboutCoM( data, xOpt )
    dataFitted = transform3D.transformRigid3DAboutP( data, xOpt, rotcentre )
    if outputErrors:
        return xOpt, dataFitted, (rms0, rmsOpt)
    else:
        return xOpt, dataFitted 
    
def fitRigidFMin( data, target, t0=None, xtol=1e-3, maxfev=0, sample=None, verbose=0, outputErrors=0 ):
    """ fits for tx,ty,tz,rx,ry,rz to transform points in data to points
    in target. Points in data and target are assumed to correspond by
    order
    """ 

    if sample is not None:
        D = _sampleData(data,sample)
        T = _sampleData(target,sample)
    else:
        D = data
        T = target

    if t0 is None:
        t0 = [0.0,0.0,0.0,0.0,0.0,0.0]
        
    def obj( x ):
        DT = transform3D.transformRigid3DAboutCoM( D, x )
        d = ( ( DT - T )**2.0 ).sum(1)
        rmsd = scipy.sqrt(d.mean())
        return rmsd
        
    t0 = scipy.array(t0)
    rms0 = scipy.sqrt( obj( t0 ).mean() )
    if verbose:
        print('initial RMS:', rms0)
        
    xOpt = fmin( obj, t0, xtol=xtol, maxiter=maxfev )
    
    rmsOpt = scipy.sqrt( obj(xOpt).mean() )
    if verbose:
        print('final RMS:', rmsOpt)
    
    dataFitted = transform3D.transformRigid3DAboutCoM( data, xOpt )
    if outputErrors:
        return xOpt, dataFitted, (rms0, rmsOpt)
    else:
        return xOpt, dataFitted 

def fitRigidScale( data, target, t0=None, xtol=1e-3, maxfev=None, sample=None, verbose=0, outputErrors=0 ):
    """ fits for tx,ty,tz,rx,ry,rz,s to transform points in data to points
    in target. Points in data and target are assumed to correspond by
    order
    """

    if sample is not None:
        D = _sampleData(data,sample)
        T = _sampleData(target,sample)
    else:
        D = data
        T = target

    if t0 is None:
        t0 = scipy.array([0.0,0.0,0.0,0.0,0.0,0.0,1.0])
    else:
        t0 = scipy.array(t0)
    
    if data.shape[0]>=t0.shape[0]:
        def obj( x ):
            DT = transform3D.transformRigidScale3DAboutCoM( D, x )
            d = ( ( DT - T )**2.0 ).sum(1)
            return d
    else:
        def obj( x ):
            DT = transform3D.transformRigidScale3DAboutCoM( D, x )
            d = ( ( DT - T )**2.0 ).sum(1)
            return d.sum()
        
    t0 = scipy.array(t0)
    rms0 = scipy.sqrt( obj( t0 ).mean() )
    if verbose: 
        print('initial RMS:', rms0)
    
    if data.shape[0]>=t0.shape[0]:
        if maxfev is None:
            maxfev = 0
        xOpt = leastsq( obj, t0, xtol=xtol, maxfev=maxfev )[0]
    else:
        xOpt = fmin( obj, t0, xtol=xtol, maxiter=maxfev )
    
    rmsOpt = scipy.sqrt( obj(xOpt).mean() )
    if verbose:
        print('final RMS:', rmsOpt)
    
    dataFitted = transform3D.transformRigidScale3DAboutCoM( data, xOpt )
    if outputErrors:
        return xOpt, dataFitted, (rms0, rmsOpt)
    else:
        return xOpt, dataFitted 

#======================================================================#
# Non correspondent fitting data fitting                               #
#======================================================================#
def fitDataRigidEPDP( data, target, xtol=1e-5, maxfev=0, t0=None, sample=None, outputErrors=0 ):
    """ fit list of points data to list of points target by minimising
    least squares distance between each point in data and closest neighbour
    in target
    """

    if sample is not None:
        D = _sampleData(data,sample)
        T = _sampleData(target,sample)
    else:
        D = data
        T = target

    if t0 is None:
        t0 = scipy.array([0.0,0.0,0.0,0.0,0.0,0.0])
        
    TTree = cKDTree( T )
    D = scipy.array(D)

    def obj( t ):
        DT = transform3D.transformRigid3DAboutCoM( D, t )
        d = TTree.query( DT )[0]
        # print d.mean()
        return d*d

    initialRMSE = scipy.sqrt(obj(t0).mean())
    tOpt = leastsq( obj, t0, xtol=xtol, maxfev=maxfev )[0]
    dataFitted = transform3D.transformRigid3DAboutCoM( data, tOpt )
    finalRMSE = scipy.sqrt(obj( tOpt ).mean())
    # print 'fitDataRigidEPDP finalRMSE:', finalRMSE
    
    if outputErrors:
        return tOpt, dataFitted, (initialRMSE, finalRMSE)
    else:
        return tOpt, dataFitted

def fitDataTranslateEPDP( data, target, xtol=1e-5, maxfev=0, t0=None, sample=None, outputErrors=0 ):
    """ fit list of points data to list of points target by minimising
    least squares distance between each point in data and closest neighbour
    in target
    """

    if sample is not None:
        D = _sampleData(data,sample)
        T = _sampleData(target,sample)
    else:
        D = data
        T = target

    if t0 is None:
        t0 = scipy.array([0.0,0.0,0.0])
        
    TTree = cKDTree( T )
    D = scipy.array(D)
    
    def obj( t ):
        DT = transform3D.transformRigid3DAboutCoM( D, scipy.hstack( (t,[0.0,0.0,0.0])) )
        d = TTree.query( list(DT) )[0]
        #~ print d.mean()
        return d*d
    
    initialRMSE = scipy.sqrt(obj(t0).mean())
    tOpt = leastsq( obj, t0, xtol=xtol, maxfev=maxfev )[0]
    dataFitted = transform3D.transformRigid3DAboutCoM( data, scipy.hstack((tOpt,[0.0,0.0,0.0])) )
    finalRMSE = scipy.sqrt(obj( tOpt ).mean())
    
    if outputErrors:
        return tOpt, dataFitted, (initialRMSE, finalRMSE)
    else:
        return tOpt, dataFitted
    
def fitDataRigidDPEP( data, target, xtol=1e-5, maxfev=0, t0=None, sample=None, outputErrors=0 ):
    """ fit list of points data to list of points target by minimising
    least squares distance between each point in target and closest neighbour
    in data
    """

    if sample is not None:
        D = _sampleData(data,sample)
        T = _sampleData(target,sample)
    else:
        D = data
        T = target

    if t0 is None:
        t0 = scipy.array([0.0,0.0,0.0,0.0,0.0,0.0])
        
    D = scipy.array(D)

    def obj( t ):
        DT = transform3D.transformRigid3DAboutCoM( D, t )
        DTTree = cKDTree( DT )
        d = DTTree.query( list(T) )[0]
        #~ print d.mean()
        return d*d
    
    initialRMSE = scipy.sqrt(obj(t0).mean())
    tOpt = leastsq( obj, t0, xtol=xtol, maxfev=maxfev )[0]
    dataFitted = transform3D.transformRigid3DAboutCoM( data, tOpt )
    finalRMSE = scipy.sqrt(obj( tOpt ).mean())

    if outputErrors:
        return tOpt, dataFitted, (initialRMSE, finalRMSE)
    else:
        return tOpt, dataFitted
    
def fitDataRigidScaleEPDP( data, target, xtol=1e-5, maxfev=0, t0=None, sample=None, outputErrors=0, scaleThreshold=None ):
    """ fit list of points data to list of points target by minimising
    least squares distance between each point in data and closest neighbour
    in target
    """

    if sample is not None:
        D = _sampleData(data,sample)
        T = _sampleData(target,sample)
    else:
        D = data
        T = target

    if t0 is None:
        t0 = scipy.array([0.0,0.0,0.0,0.0,0.0,0.0,1.0])
            
    TTree = cKDTree( T )
    D = scipy.array(D)
    
    if scaleThreshold is not None:
        # print 'scale penalty on'
        def obj( t ):
            DT = transform3D.transformRigidScale3DAboutCoM( D, t )
            d = TTree.query( list(DT) )[0]
            s = max(t[-1], 1.0/t[-1])
            if s>scaleThreshold:
                sw = 1000.0*s
            else:
                sw = 1.0
            return d*d + sw
    else:
        def obj( t ):
            DT = transform3D.transformRigidScale3DAboutCoM( D, t )
            d = TTree.query( list(DT) )[0]
            return d*d
    
    initialRMSE = scipy.sqrt(obj(t0).mean())
    tOpt = leastsq( obj, t0, xtol=xtol, maxfev=maxfev )[0]
    dataFitted = transform3D.transformRigidScale3DAboutCoM( data, tOpt )
    finalRMSE = scipy.sqrt(obj( tOpt ).mean())

    if outputErrors:
        return tOpt, dataFitted, (initialRMSE, finalRMSE)
    else:
        return tOpt, dataFitted
    
def fitDataRigidScaleDPEP( data, target, xtol=1e-5, maxfev=0, t0=None, sample=None, outputErrors=0, scaleThreshold=None ):
    """ fit list of points data to list of points target by minimising
    least squares distance between each point in target and closest neighbour
    in data
    """

    if sample is not None:
        D = _sampleData(data,sample)
        T = _sampleData(target,sample)
    else:
        D = data
        T = target

    if t0 is None:
        t0 = scipy.array([0.0,0.0,0.0,0.0,0.0,0.0,1.0])
        
    D = scipy.array(D)

    if scaleThreshold is not None:
        def obj( t ):
            DT = transform3D.transformRigidScale3DAboutCoM( D, t )
            DTTree = cKDTree( DT )
            d = DTTree.query(T)[0]
            if s>scaleThreshold:
                sw = 1000.0*s
            else:
                sw = 1.0
            return d*d + sw
    else:
        def obj( t ):
            DT = transform3D.transformRigidScale3DAboutCoM( D, t )
            DTTree = cKDTree( DT )
            d = DTTree.query( T )[0]
            return d*d
    
    initialRMSE = scipy.sqrt(obj(t0).mean())    
    tOpt = leastsq( obj, t0, xtol=xtol, maxfev=maxfev )[0]
    dataFitted = transform3D.transformRigidScale3DAboutCoM( data, tOpt )
    finalRMSE = scipy.sqrt(obj( tOpt ).mean())
    
    if outputErrors:
        return tOpt, dataFitted, (initialRMSE, finalRMSE)
    else:
        return tOpt, dataFitted
#===========================================================================#


def fitSphere( X ):
    """
    least squares fits the sphere centre and radius to a cloud of points X
    """ 
    B = (X**2.0).sum(1)
    A = scipy.hstack( [2.0*X, scipy.ones(X.shape[0])[:,scipy.newaxis]] )
    x, res, rank, s = lstsq( A, B )
    
    a, b, c, m = x
    r = scipy.sqrt( m + a*a + b*b + c*c )
    return (a,b,c), r
