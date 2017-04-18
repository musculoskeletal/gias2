"""
functions and classes for radial basis functions
"""

import sys
import scipy as sp
from scipy.spatial.distance import pdist,cdist,squareform
from scipy.linalg import lu_factor, lu_solve, svd, lstsq, solve, inv, qr
from scipy.spatial import cKDTree, KDTree
from scipy.interpolate import Rbf as scipyRBF
# import alglibRBF
from scipy.optimize import leastsq
import cPickle
import pdb

#=============================================================================#
# Basis functions
#=============================================================================#
def polyCubic3D( x, y, z ):
    
    X0 = sp.ones(x.shape[0])
    X1 = x
    X2 = y
    X3 = z
    X4 = x*x
    X5 = y*y
    X6 = z*z
    X7 = x*y
    X8 = x*z
    X9 = y*z
    X10 = X4*x  #x^3
    X11 = X5*y  #y^3
    X12 = X6*z  #z^3
    X13 = X4*y  #x^2y
    X14 = X5*x  #xy^2
    X15 = X4*z  #x^2z
    X16 = X6*x  #xz^2
    X17 = X5*z  #y^2z
    X18 = X6*y  #yz^2
    X19 = x*y*z
    
    return sp.array([X0,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X19])

def polyLinear3D( x, y, z ):

    X0 = sp.ones(x.shape[0])
    X1 = x
    X2 = y
    X3 = z
    
    return sp.array([X0,X1,X2,X3])

def polyConst3D( x, y, z ):

    X0 = sp.ones(x.shape[0])
    
    return sp.array([X0,])

def makeBasisBiharmonic():
    def b( r ):
        return r
        
    return b
    
def makeBasisGaussian( s ): 
    """
    Gaussian basis function, normalised to integrate to 1
    """
    def b( r ):
        #~ return 1.0/sp.sqrt(2.0*sp.pi*s) * sp.exp( -0.5*r*r/(s*s) )
        return 1.0/sp.sqrt(2.0*sp.pi*s) * sp.exp( -0.5*r*r/(s*s) ) +1.0
        
    return b
    
def makeBasisGaussianNonUniformWidth( s, scaling ): 
    """
    Gaussian basis function, normalised to integrate to 1
    """
    # print 'making gaussian non-uniform width RBF'
    S = s*scaling
    def b( r ):
        return sp.exp( -r*r/(S*S) )
        
    return b

def makeBasisNormalisedGaussianNonUniformWidth( s, scaling ):   
    """
    Gaussian basis function, normalised to integrate to 1
    """
    # print 'making normalised gaussian non-uniform width RBF'
    
    S = s*scaling
    def b( r ):
        return 1.0/sp.sqrt(2.0*sp.pi*S) * sp.exp( -0.5*r*r/(S*S) ) + 1.0
        
    return b

#~ def makeBasisGaussianBilinear( s ):
    #~ 
    #~ def b( r, x, y ):
        #~ g = 1.0/sp.sqrt(2.0*sp.pi*s) * sp.exp( -0.5*r*r/(s*s) ) +1.0
        #~ return g, x, y
#~ 
#~ def makeBasisGaussianBasis( s ):
    #~ 
    #~ def b( r, x, y, z ):
        #~ g = 1.0/sp.sqrt(2.0*sp.pi*s) * sp.exp( -0.5*r*r/(s*s) ) +1.0
        #~ return g, x, y, z
        
def makeBasisMultiQuadric( c ):
    def b( r ):
        return sp.sqrt( c*c + r*r )
        
    return b

def makeBasisTPS():
    def b( r ):
        x = r*r*sp.log(r)
        return sp.where( sp.isfinite(x), x, 0.0 )

    return b

def makeBasisWendlandsC31NonUniformWidth(s, scaling):
    """
    see Fornefett, Rohr, Steihl (1999) elastic registration of medical
    images using RBFs with compact support
    
    3-dimensional, c2 continuous
    """
    
    def b(r):
        rNorm = r/(s*scaling)
        return sp.where( rNorm<1.0, (4.0*rNorm+1)*(1-rNorm)**4.0, 0.0 )
        #~ return (4.0*rNorm+1)*(1-rNorm)**4.0
        
    return b

def makeBasisWendlandsC32NonUniformWidth(s, scaling):
    """
    see Fornefett, Rohr, Steihl (1999) elastic registration of medical
    images using RBFs with compact support
    
    3-dimensional, c4 continuous
    """
    
    def b(r):
        rNorm = r/(s*scaling)
        return sp.where( rNorm<1.0, (35.0*rNorm**2.0 + 18.0*rNorm + 3.0)*(1.0-rNorm)**6.0, 0.0 )
        #~ return (4.0*rNorm+1)*(1-rNorm)**4.0
        
    return b
    
RBFBases = {'gaussian': makeBasisGaussian,
            'gaussianNonUniformWidth': makeBasisGaussianNonUniformWidth,
            'normalisedgaussianNonUniformWidth': makeBasisNormalisedGaussianNonUniformWidth,
            'TPS': makeBasisTPS,
            'multiquadric': makeBasisMultiQuadric,
            'biharmonic': makeBasisBiharmonic,
            'WC31NonUniformWidth': makeBasisWendlandsC31NonUniformWidth,
            'WC32NonUniformWidth': makeBasisWendlandsC32NonUniformWidth,
            }

#=============================================================================#
# Util functions
#=============================================================================#
def estimateNonUniformWidth( X, k=2 ):
    
    XTree = cKDTree( X )
    XNNd = XTree.query( list(X), k=k+1 )[0][:,1:] # 1st NN is always point itself
    
    s = sp.sqrt( (XNNd**2.0).sum(1) ) / k
    return s

def xDist( x, X ):
    """ calculates the distance of x to each point in X
    """
    v = ( X - x )
    return sp.sqrt( (v*v).sum(1) )
    
def xDist1D( x, X ):
    return X - x

def fitData( C, basis, dataX, dataU ):
    
    print('fitting {} knots to {} data points'.format(len(C), len(dataX)))
    # split distance matrix calculation in groups to save memory
    groupSize = 5000
    #~ pdb.set_trace()
    if len(dataX)>groupSize:
        A = sp.zeros((len(dataX), len(C)))
        for i in xrange( 0, len(dataX), groupSize ):
            # print(str(i)+' - '+str(i+groupSize))
            A[i:i+groupSize,:] = basis(cdist(dataX[i:i+groupSize,:], C))
    else:
        A = basis( cdist( dataX, C ) )
    
    #~ pdb.set_trace()
    #~ A = sp.array( [ basis( xDist( d, C ) ) for d in dataX ] )
    W, residual, rank, singVal = lstsq( A, dataU.copy(), overwrite_a=True, overwrite_b=True )
    W = W.T
    
    return W, (residual, rank, singVal)

def fitDataPoly3D( C, basis, dataX, dataU, poly ):
    
    print('fitting {} knots to {} data points'.format(len(C), len(dataX)))
    # split distance matrix calculation in groups to save memory
    groupSize = 5000
    #~ pdb.set_trace()
    if len(dataX)>groupSize:
        A = sp.zeros((len(dataX), len(C)))
        for i in xrange( 0, len(dataX), groupSize ):
            # print(str(i)+' - '+str(i+groupSize))
            A[i:i+groupSize,:] = basis( cdist( dataX[i:i+groupSize,:], C ) )
    else:
        A = basis( cdist( dataX, C ) )
    
    # append polynomial
    # polynomial terms in 2nd dim, points in 1st dim
    P = poly( dataX[:,0],dataX[:,1],dataX[:,2] ).T
    
    AAUpper = sp.hstack( [A,P] )
    
    #~ pdb.set_trace()
    #~ AALower = sp.hstack( [P.T, sp.zeros([P.shape[1], AAUpper.shape[1]-P.shape[0]])] )
    #~ AA = sp.vstack([AAUpper, AALower])
    
    AA = AAUpper
    
    W, residual, rank, singVal = lstsq( AA, dataU.copy(), overwrite_a=True, overwrite_b=True )
    W = W.T
    
    rbfW = W[:C.shape[0]]
    polyCoeff = W[C.shape[0]:]
    
    return rbfW, polyCoeff, (residual, rank, singVal)

def fitDataQR( C, basis, dataX, dataU ):
    
    print('fitting {} knots to {} data points'.format(len(C), len(dataX)))
    # split distance matrix calculation in groups to save memory
    groupSize = 5000
    #~ pdb.set_trace()
    if len(dataX)>groupSize:
        A = sp.zeros((len(dataX), len(C)))
        for i in xrange( 0, len(dataX), groupSize ):
            # print str(i)+' - '+str(i+groupSize)
            A[i:i+groupSize,:] = basis( cdist( dataX[i:i+groupSize,:], C ) )
    else:
        A = basis( cdist( dataX, C ) )
    
    Q,R = qr(A)
    P = sp.dot(Q.T, dataU.copy())
    pdb.set_trace()
    W = sp.dot(inv(R[:R.shape[1],:]), P)    # not working, need a mapping matrix to match matrix shapes
    residual = dataU - sp.dot(A,W)
    rank = -1
    singVal = -1
    W = W.T
    
    return W, (residual, rank, singVal)


def fitComponentFieldKnotCoordWidth( RBFF, dataX, dataU, fullOutput=False, xtol=1e-3, maxfev=0 ):
    """
    nonlinear least squares fit of knot coordinates and width
    """
    
    nKnots = len( RBFF.C )
    
    def obj( X ):
        
        knotCoords = X[:(nKnots*RBFF.nComponents)].reshape((nKnots, RBFF.nComponents))
        knotWidths = X[(nKnots*RBFF.nComponents):]
        RBFF.setCentres( knotCoords, width=widths )
        W, (res, rank, singVal) = RBFF.fitData( dataX, dataU, fullOutput=True )
        sys.stdout.write('rbf fit rmse: %6.4f\r'%(sp.sqrt((res**2.0).mean())))
        sys.stdout.flush()
        return res
    
    x0 = sp.hstack( [RBFF.C.ravel(), RBFF.CWidth] )
    
    xOpt, cov_x, fitInfo, mesg, ier = leastsq( obj, x0, xtol=xtol, maxfev=maxfev, full_output=1 )
    knotCoords = xOpt[:(nKnots*RBFF.nComponents)].reshape((nKnots, RBFF.nComponents))
    knotWidths = xOpt[(nKnots*RBFF.nComponents):]
    RBFF.setCentres( knotCoords, width=knotWidths )
    
    fitRMSE = sp.sqrt((fitInfo['fvec']**2.0).mean())
    
    if not fullOutput:
        return RBFF
    else:
        return RBFF, fitRMSE, fitInfo

def fitComponentFieldKnotWidth( RBFF, dataX, dataU, fullOutput=False, xtol=1e-3, maxfev=0 ):
    """
    nonlinear least squares fit of knot width
    """
    
    nKnots = len( RBFF.C )
    knotCoords = sp.array( RBFF.C )
    
    def obj( widths ):
        RBFF.setCentres( knotCoords, width=widths )
        W, (res, rank, singVal) = RBFF.fitData( dataX, dataU, fullOutput=True )
        err = sp.sqrt(((RBFF.evalMany( dataX ).T - dataU)**2.0).sum(1))
        sys.stdout.write('rbf fit rmse: %6.4f\r'%(sp.sqrt((err**2.0).mean())))
        sys.stdout.flush()
        return err
    
    x0 = sp.array(RBFF.CWidth)
    try:
        xOpt, cov_x, fitInfo, mesg, ier = leastsq( obj, x0, xtol=xtol, maxfev=maxfev, full_output=1 )
    except:
        pdb.set_trace()
        
    RBFF.setCentres( knotCoords, width=xOpt )
    fitRMSE = sp.sqrt((fitInfo['fvec']**2.0).mean())
    
    if not fullOutput:
        return RBFF
    else:
        return RBFF, fitRMSE, fitInfo

polynomials = {0:polyConst3D,
               1:polyLinear3D,
               3:polyCubic3D,
               }

#=============================================================================#
# Main RBF Classes
#=============================================================================#
class RBFField( object ):
    
    usePoly = -1
    CWidthNN = 3
    
    def __init__( self ):
        self.W = None
        self.C = None
        self.U = None
        self.polyOrder = None
        self.polyCoeffs = None
        self.poly = None
        self.basis = None
        self.basisType = None
        self.basisArgs = {}
        self.CWidth = None
    
    def save( self, filename ):
        
        d = {'CWidthNN': self.CWidthNN,
             'C': self.C,
             'U': self.U,
             'basisType': self.basisType,
             'basisArgs': self.basisArgs,
             'W': self.W,
             'polyCoeffs': self.polyCoeffs,
             'polyOrder': self.polyCoeffs,
             }
        
        with open(filename, 'w') as f:
            cPickle.dump(d, f, protocol=2)
    
    def load( self, filename ):
        with open(filename, 'r') as f:
            d = cPickle.load(f)
        
        self.CWidthNN = d['CWidthNN']
        self.W = d['W']
        self.C = d['C']
        self.U = d['U']
        self.basisType = d['basisType']
        self.basisArgs = d['basisArgs']
        print self.basisType
        if self.basisArgs!=None:
            if isinstance(self.basisArgs, list) and self.basisType=='gaussian':
                self.basisArgs = {'s':self.basisArgs[0]}
                
            self.setBasis( RBFBases[self.basisType](**self.basisArgs) )
        else:
            self.setBasis( RBFBases[self.basisType]() )
        
        self.setCentres( self.C )
        
        try:
            self.polyOrder = d['polyOrder']
            self.polyCoeffs = d['polyCoeffs']
        except KeyError:
            self.polyCoeffs = None
            self.polyOrder = None
        else:
            if self.polyOrder != None:
                self.setPolynomial( self.polyOrder ) 
    
    def setPolynomial( self, polyOrder ):
        self.poly = polynomials[polyOrder]
        self.usePoly = 1
        
    def setCentres( self, C ):
        self.C = C
        M = squareform( pdist( self.C ) )
        self.CWidth = estimateNonUniformWidth( C, k=self.CWidthNN )
        
        if 'NonUniformWidth' in self.basisType:
            self.basisArgs['s'] = self.CWidth
            self.makeBasis( self.basisType, self.basisArgs )
            
        self.M = self.basis( M )
        
    def setValues( self, U ):
        self.U = U
        
    def setBasis (self, phi ):
        self.basis = phi

    def makeBasis( self, basisType, basisArgs=None ):
        self.basisType = basisType
        self.basisArgs = basisArgs
        
        if self.basisArgs!=None:
            self.setBasis( RBFBases[self.basisType](**self.basisArgs) )
        else:
            self.setBasis( RBFBases[self.basisType]() )
            
    def calcWeights( self ):
        print 'interpolating '+str(len(self.C))+' knots'
        print 'calculating distances...'
        R = pdist( self.C )
        RS = squareform( R )
        A = self.basis( RS )
        
        print 'solving system...'
        self.W = lu_solve( lu_factor(A), self.U.copy(), overwrite_b=True )
        
    def eval( self, x ):
        """
        evaluate the value of the field at point x
        """
        r = xDist( x, self.C )
        B = self.basis(r)
        y = ( self.W * B ).sum()
        return y
    
    def evalMany( self, x ):
        """
        evaluate the value of the field at points x
        """
        print 'evaluating at '+str(len(x))+' points'
        # calculate distance from each x to each rbf centre
        # row i of D contains the distances of each rbf centre to point
        # x[i] 
        
        # split distance matrix calculation in groups to save memory
        groupSize = 10000
        if len(x)>groupSize:
            y = sp.zeros(len(x))
            for i in xrange( 0, len(x), groupSize ):
                #~ print str(i)+' - '+str(i+groupSize)
                r = cdist( x[i:i+groupSize,:], self.C, 'euclidean' )
                B = self.basis( r )
                y[i:i+groupSize] = (self.W * B).sum(1)
            
        else:
            r = cdist( x, self.C, 'euclidean' )
            B = self.basis( r )
            y = (self.W * B).sum(1)
        
        
        return y
        
    def evalManyPoly3D( self, x ):
        """
        evaluate the value of the field at points x
        """
        print 'evaluating at '+str(len(x))+' points'
        # calculate distance from each x to each rbf centre
        # row i of D contains the distances of each rbf centre to point
        # x[i] 
        
        # split distance matrix calculation in groups to save memory
        groupSize = 10000
        if len(x)>groupSize:
            y = sp.zeros(len(x))
            for i in xrange( 0, len(x), groupSize ):
                #~ print str(i)+' - '+str(i+groupSize)
                xGroup = x[i:i+groupSize,:]
                r = cdist( xGroup, self.C, 'euclidean' )
                B = self.basis( r )
                try:
                    y[i:i+groupSize] = (self.W * B).sum(1) + sp.dot( self.poly(xGroup[:,0],xGroup[:,1],xGroup[:,2]).T, self.polyCoeffs )
                except ValueError:
                    pdb.set_trace()
        else:
            r = cdist( x, self.C, 'euclidean' )
            B = self.basis( r )
            y = (self.W * B).sum(1) + sp.dot( self.poly(x[:,0],x[:,1],x[:,2]).T, self.polyCoeffs )
        
        
        return y
        
    def fitData( self, dataX, dataU, fullOutput=False ):
        """
        calculate weights for self.C to fit to field sampled at dataX
        with field values dataU
        """
        #~ print 'fitting data...'
        self.W, extraInfo = fitData( self.C, self.basis, dataX, dataU )
        #~ self.W, extraInfo = fitDataQR( self.C, self.basis, dataX, dataU )
        if fullOutput:
            return self.W, extraInfo
        else:
            return self.W
            
    def fitDataPoly3D( self, dataX, dataU, fullOutput=False ):
        """
        calculate weights for self.C to fit to field sampled at dataX
        with field values dataU
        """
        #~ print 'fitting data...'
        self.W, self.polyCoeffs, extraInfo = fitDataPoly3D( self.C, self.basis, dataX, dataU, self.poly )
        #~ self.W, extraInfo = fitDataQR( self.C, self.basis, dataX, dataU )
        if fullOutput:
            return self.W, self.polyCoeffs, extraInfo
        else:
            return self.W, self.polyCoeffs
            
    
        
class RBFComponentsField( object ):
    """
    Multivarite RBF field

    The values of the field has nComponent number of components.
    """
    
    CWidthNN = 5
    
    def __init__( self, nComponents ):
        self.W = None
        self.C = None
        self.U = None
        self.basis = None
        self.M = None
        self.nComponents = nComponents
        self.basisType = None
        self.basisArgs = {}
        self.CWidth = None
        self.polyOrder = None
        self.polyCoeffs = None
        self.poly = None
    
    def save( self, filename ):
        
        d = {'CWidthNN': self.CWidthNN,
             'C': self.C,
             'U': self.U,
             'basisType': self.basisType,
             'basisArgs': self.basisArgs,
             'W': self.W,
             'nComponents': self.nComponents,
             'polyCoeffs': self.polyCoeffs,
             'polyOrder': self.polyCoeffs,
             }
        
        with open(filename, 'w') as f:
            cPickle.dump(d, f, protocol=2)
    
    def load( self, filename ):
        with open(filename, 'r') as f:
            d = cPickle.load(f)
        
        self.CWidthNN = d['CWidthNN']
        self.nComponents = d['nComponents']
        self.W = d['W']
        self.C = d['C']
        self.U = d['U']
        self.basisType = d['basisType']
        self.basisArgs = d['basisArgs']
        print self.basisType
        if self.basisArgs!=None:
            if isinstance(self.basisArgs, list) and self.basisType=='gaussian':
                self.basisArgs = {'s':self.basisArgs[0]}
                
            self.setBasis( RBFBases[self.basisType](**self.basisArgs) )
        else:
            self.setBasis( RBFBases[self.basisType]() )
        
        self.setCentres( self.C )
        
        try:
            self.polyOrder = d['polyOrder']
            self.polyCoeffs = d['polyCoeffs']
        except KeyError:
            self.polyCoeffs = None
            self.polyOrder = None
        else:
            if self.polyOrder != None:
                self.setPolynomial( self.polyOrder ) 
        
    def setCentres( self, C, width=None ):
        self.C = C
        M = squareform( pdist( self.C ) )
        if width!=None:
            self.CWidth = width
        else:
            self.CWidth = estimateNonUniformWidth( C, k=self.CWidthNN )
        
        if 'NonUniformWidth' in self.basisType:
            self.basisArgs['s'] = self.CWidth
            self.makeBasis( self.basisType, self.basisArgs )
            
        self.M = self.basis( M )
    
    def setPolynomial( self, polyOrder ):
        self.poly = polynomials[polyOrder]
        self.usePoly = 1
        
    def setCentreValues( self, U ):
        if U.shape[0]!=self.nComponents:
            raise ValueError, 'incorrect number of components'
        else:
            self.U = U
        
    def setBasis (self, phi ):
        self.basis = phi
        
    def makeBasis( self, basisType, basisArgs ):
        self.basisType = basisType
        self.basisArgs = basisArgs
        
        if self.basisArgs!=None:
            self.setBasis( RBFBases[self.basisType](**self.basisArgs) )
        else:
            self.setBasis( RBF.RBFBases[self.basisType]() )

    def calcWeights( self ):
        #~ print 'solving system...'
        #~ LUA = lu_factor(A)
        #~ self.W = sp.array([lu_solve( LUA, u ) for u in self.U])

        self.W = solve( self.M, self.U.T ).T
        
    def eval( self, x ):
        """
        evaluate the value of the field at point x
        """
        r = xDist( x, self.C )
        B = self.basis( r )
        y = ( self.W * B ).sum(1)
        return y
    
    def evalMany( self, x ):
        """
        evaluate the value of the field at points x
        """
        
        #~ # calculate distance from each x to each rbf centre
        #~ # row i of D contains the distances of each rbf centre to point
        #~ # x[i] 
        #~ D = cdist( x, self.C, 'euclidean' )
        #~ B = self.basis( D )
        #~ y = sp.dot(self.W, B.T)
        #~ return y
        if self.polyOrder!=None:
            return self.evalManyPoly3D(x)
            
        print 'evaluating at '+str(len(x))+' points'
        # calculate distance from each x to each rbf centre
        # row i of D contains the distances of each rbf centre to point
        # x[i] 
        
        # split distance matrix calculation in groups to save memory
        groupSize = 10000
        if len(x)>groupSize:
            y = sp.zeros( (self.nComponents, len(x)), dtype=float )
            for i in xrange( 0, len(x), groupSize ):
                #~ print str(i)+' - '+str(i+groupSize)
                r = cdist( x[i:i+groupSize,:], self.C, 'euclidean' )
                B = self.basis( r )
                y[:,i:i+groupSize] = sp.dot(self.W, B.T)
            
        else:
            r = cdist( x, self.C, 'euclidean' )
            B = self.basis( r )
            y = sp.dot(self.W, B.T)
        
        return y
        
    def evalManyPoly3D( self, x ):
        """
        evaluate the value of the field at points x
        """
        print 'evaluating at '+str(len(x))+' points poly'
        # calculate distance from each x to each rbf centre
        # row i of D contains the distances of each rbf centre to point
        # x[i] 
        
        # split distance matrix calculation in groups to save memory
        groupSize = 10000
        if len(x)>groupSize:
            y = sp.zeros( (self.nComponents, len(x)), dtype=float )
            for i in xrange( 0, len(x), groupSize ):
                #~ print str(i)+' - '+str(i+groupSize)
                xGroup = x[i:i+groupSize,:]
                r = cdist( xGroup, self.C, 'euclidean' )
                B = self.basis( r )
                try:
                    y[:,i:i+groupSize] = sp.dot(self.W, B.T) + sp.dot( self.poly(xGroup[:,0],xGroup[:,1],xGroup[:,2]).T, self.polyCoeffs )
                except ValueError:
                    pdb.set_trace()
        else:
            r = cdist( x, self.C, 'euclidean' )
            B = self.basis( r )
            #~ pdb.set_trace()
            y = sp.dot(self.W, B.T) + sp.dot( self.poly(x[:,0],x[:,1],x[:,2]).T, self.polyCoeffs )
        
        return y
        
    def fitData( self, dataX, dataU, fullOutput=True ):
        """
        calculate weights for self.C to fit to field sampled at dataX
        with field values dataU
        """
        if self.polyOrder!=None:
            return self.fitDataPoly3D(dataX, dataU, fullOutput)
        
        if dataU.shape[1]!=self.nComponents:
            #~ pdb.set_trace()
            raise ValueError, 'incorrect number of components in data'
            
        else:
            # print 'fitting data...'
            self.W, extraInfo = fitData( self.C, self.basis, dataX, dataU )
            if fullOutput:
                return self.W, extraInfo
            else:
                return self.W

    def fitDataPoly3D( self, dataX, dataU, fullOutput=False ):
        """
        calculate weights for self.C to fit to field sampled at dataX
        with field values dataU
        """
        if dataU.shape[1]!=self.nComponents:
            #~ pdb.set_trace()
            raise ValueError, 'incorrect number of components in data'
            
        else:
            # print 'fitting data poly...'
            self.W, self.polyCoeffs, extraInfo = fitDataPoly3D( self.C, self.basis, dataX, dataU, self.poly )
            #~ self.W, extraInfo = fitDataQR( self.C, self.basis, dataX, dataU )
            if fullOutput:
                return self.W, self.polyCoeffs, extraInfo
            else:
                return self.W, self.polyCoeffs  

