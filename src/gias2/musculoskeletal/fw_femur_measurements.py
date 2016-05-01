"""
FILE: fw_femur_measurements.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION:
functions and classes for taking anthropometric measurements from 
the fieldwork femur model

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import copy
import pickle
import scipy
from scipy.spatial import distance
from scipy.optimize import fmin, leastsq

from gias2.common import geoprimitives as FT
from gias2.registration import alignment_analytic as alignment
from gias2.musculoskeletal import fw_model_landmarks as fml

def _norm( v ):
    return scipy.divide( v, scipy.sqrt((scipy.array(v)**2.0).sum()) )

def _mag( v ):
    return scipy.sqrt((scipy.array(v)**2.0).sum())

class _Plane3D( object ):
    """
    ax + by + cz + d = 0
    """
    
    def __init__(self, abcd=None, NP=None):
        if abcd!=None:
            self.a, self.b, self.c, self.d = abcd
            self.N = scipy.array([self.a, self.b, self.c])
            self.P = None
        if NP!=None:
            self.N, self.P = NP
            self.a, self.b, self.c = self.N
            self.d = -self.a*self.P[0] - self.b*self.P[1] - self.c*self.P[2]
            
    def distanceToPoint( self, p ):
        
        D = ( self.a*p[0] + self.b*p[1] + self.c*p[2] + self.d )\
            / (self.a**2.0 + self.b**2.0 + self.c**2.0)
            
        pointOnPlane = p - self.N*D
        
        return abs(D), pointOnPlane

    def distanceToPoints( self, p ):
        
        D = ( self.a*p[:,0] + self.b*p[:,1] + self.c*p[:,2] + self.d )\
            / (self.a**2.0 + self.b**2.0 + self.c**2.0)
        
        # print p.shape
        # print (self.N*D).shape
        pointOnPlane = p - self.N*D[:,scipy.newaxis]
        
        return abs(D), pointOnPlane
        
    def projectClosePoints( self, P, maxDist ):
        
        # distances = scipy.zeros(P.shape[0], dtype=float)
        # projections = scipy.zeros([P.shape[0],2], dtype=float)

        distances, projections = self.distanceToPoints(P)

        # for i, p in enumerate(P):
        #   distances[i], projections[i] = self.distanceToPoint( p )
        
        closeIndices = scipy.where(distances<maxDist)[0]
        closeProjections = projections[closeIndices]
        
        return closeProjections, closeIndices

class measurement( object ):
    
    def __init__( self, name, value=None ):
        self.name = name
        self.value = value

class FemurMeasurements( object ):
    
    measurements = { 
                     'shaft_axis': None,
                     'neck_axis': None,
                     'head_diameter': None,
                     'neck_shaft_angle': None,
                     'neck_width': None,
                     'femoral_axis_length': None,
                     'length_long': None,
                     'length_short': None,
                     'subtrochanteric_diameter': None,
                     'subtrochanteric_width': None,
                     'midshaft_diameter': None,
                     'midshaft_width': None,
                     'epicondylar_width': None,
                     }
                     
    headElems = fml._femurHeadElems
    shaftElems = fml._femurShaftElems
    neckElems = fml._femurNeckElems
    neckLongElems = fml._femurNeckLongElems
    subtrochanterNodes = fml._femurSubtrochanterNodes
    midshaftNodes = fml._femurMidShaftNodes
    epiCondylarMCondElems = fml._femurMedialEpicondyleElems
    epiCondylarLCondElems = fml._femurLateralEpicondyleElems
    condyleAlignmentNodes = fml._femurCondyleAlignmentNodes
    gTrocElems = fml._femurGreaterTrochanterElems
    mCondElems = fml._femurMedialCondyleElems
    lCondElems = fml._femurLateralCondyleElems
    
    epD = [20,20]
                     
    def __init__( self, GF=None, side='left' ):
        if GF!=None:
            self.GF = copy.deepcopy( GF )
            if not self.GF.ensemble_field_function.is_flat():
                self.GF.flatten_ensemble_field_function()
        else:
            self.GF = None
        
        self.side = side
        self.EP = None
        self.EPMap = None
        self.neckAxis = None
        self.shaftAxis = None
        self.xAxis = None
        self.yAxis = None
        self.epicondylarAxis = None
        
        self.measurements = {}
        self.shaftAligned = False
        
        if self.GF != None:
            self._evaluateGF()
            self.calcShaftAxis()
            self.calcNeckAxis()
            self.calcXYAxes( alignGF=False )
        
        #~ self.alignToShaftAxis()
        #~ 
        #~ self._evaluateAlignedGF()
    
    def saveMeasurements( self, filename ):
        with open( filename, 'w' ) as f:
            pickle.dump( self.measurements, f, protocol=2 )
            
    def loadMeasurements( self, filename ):
        with open( filename, 'r' ) as f:
            self.measurements = pickle.load( f )
            
        if self.measurements['shaft_axis'] != None:
            self.shaftAxis = FT.Line3D(self.measurements['shaft_axis'].value[0],
                                       self.measurements['shaft_axis'].value[1])
        
        if self.measurements['neck_axis'] != None:
            self.neckAxis = FT.Line3D(self.measurements['neck_axis'].value[0],
                                      self.measurements['neck_axis'].value[1])
        
        #~ if self.measurements['epicondylar_axis'] != None:
            #~ self.neckAxis = FT.Line3D(self.measurements['epicondylar_axis'].value[0],
                                      #~ self.measurements['epicondylar_axis'].value[1])
    
    def calcMeasurements( self ):
        self.calcHeadDiameter()
        self.calcNeckWidth()
        self.calcNeckShaftAngle()
        self.calcFemoralAxisLength()
        self.calcFemoralLengthLong()
        self.calcFemoralLengthShort()
        self.calcSubtrochantericDiameter()
        self.calcSubtrochantericWidth()
        self.calcMidshaftDiameter()
        self.calcMidshaftWidth()
        self.calcEpicondylarWidth()
        self.calcAnteversion()

    def getMeasurementsDict(self, measurementNames):
        d = {}
        for n in measurementNames:
            d[n] = self.measurements[n].value

        return d
    
    def printMeasurements( self ):
        m = list(self.measurements.keys())
        m.sort()
        for mi in m:
            print('{} : {}'.format(mi, self.measurements[mi].value))
        
    def _evaluateGF( self ):
        self.EP = self.GF.evaluate_geometric_field( self.epD ).T
        self.EPMap = self.GF.getElementPointI(self.epD,
                        sorted(self.GF.ensemble_field_function.mesh.elements.keys()))
    
    def _evaluateAlignedGF( self ):
        self.EPAligned = self.GF.evaluate_geometric_field( self.epD ).T
    
    def calcShaftAxis( self ):
        
        shaftEP = scipy.vstack([self.EP[self.EPMap[e],:] for e in self.shaftElems])
        
        # for elems in self.shaftElems:
        #   shaftEPi.append( scipy.hstack([self.EPMap[e] for e in elems]) )
        
        # shaftEP = self.EP[ scipy.hstack( shaftEPi ) ]
        
        # initialise shaft axis
        self.shaftAxis = FT.Line3D( [0.0,0.0,1.0], shaftEP.mean(0) )
        
        # fit shaft axis
        self.shaftAxis, fittedAxesParams, fittedRMSE = FT.fitAxis3D( shaftEP, self.shaftAxis )
        if self.shaftAxis.a[2] < 0:
            self.shaftAxis.a *= -1.0
            
        m = measurement( 'shaft_axis', ( self.shaftAxis.a, self.shaftAxis.b ) )
        self.measurements['shaft_axis'] = m
        
    def calcNeckAxis( self ):
        """
        fits axis through neck elements and head elements
        """

        neckEP = scipy.vstack([self.EP[self.EPMap[e],:] for e in self.neckElems])
        headEP = scipy.vstack([self.EP[self.EPMap[e],:] for e in self.headElems])
        EP = scipy.vstack([neckEP, headEP])
        
        # neckEPi = []
        # for submesh, elems in self.neckElems:
        #   neckEPi.append( scipy.hstack([self.EPMap[submesh][e] for e in elems]) )
        
        # neckEP = self.EP[ scipy.hstack( neckEPi ) ]
        # headEP = self.EP[ scipy.hstack(self.EPMap[self.regionName2ElementMap['head']]) ]
        # EP = scipy.vstack( [neckEP, headEP] )
        
        # initialise neck axis
        initialNeckAxis = FT.Line3D( headEP.mean(0)-neckEP.mean(0), neckEP.mean(0) )
        
        # rough fit neck axis
        initialNeckAxis, fittedAxesParams, fittedRMSE = FT.fitAxis3D( EP, initialNeckAxis )
        # m = measurement( 'neck_axis', ( self.neckAxis.a, self.neckAxis.b ) )
        # self.measurements['neck_axis'] = m

        # fine fit by finding plane of min xsection area
        neckPlaneNormal, neckPlaneOrigin, neckPlanePoints, searchP0, searchP1 = self._findNeckPlaneMinXSection(initialNeckAxis)

        self.neckAxis = FT.Line3D(neckPlaneNormal, neckPlaneOrigin)
        self.neckAxis._planePoints = neckPlanePoints
        self.neckAxis._searchP0 = searchP0
        self.neckAxis._searchP1 = searchP1
        m = measurement( 'neck_axis', ( neckPlaneNormal, neckPlaneOrigin ) )
        self.measurements['neck_axis'] = m

    def _findNeckPlaneMinXSection(self, initialNeckAxis):   
        # use neck EPs
        neckEP = scipy.vstack([self.EP[self.EPMap[e],:] for e in self.neckLongElems])

        # neckEPi = []
        # for submesh, elems in self.neckLongElems:
        #   neckEPi.append( scipy.hstack([self.EPMap[submesh][e] for e in elems]) )         
        # neckEP = self.EP[ scipy.hstack( neckEPi ) ]
        
        axisPoints = []
        axisPointsT = []
        for d in neckEP:
            p,t = initialNeckAxis.findClosest(d)
            axisPoints.append(p)
            axisPointsT.append(t)
            
        #~ axisPoints = scipy.array([self.neckAxis.findClosest(d)[0] for d in neckEP])
        initPI = scipy.sqrt(((neckEP - axisPoints)**2.0).sum(1)).argmin()
        initP = axisPoints[ initPI ]
        initT = axisPointsT[ initPI ]
        searchP0 = initialNeckAxis.eval( min(axisPointsT)+0.0 )
        searchP1 = initialNeckAxis.eval( max(axisPointsT)-5.0 )
        
        neckXSectioner = NeckMinXSectionFit( neckEP )
        # initial search constrained to neck axis
        planePoints1, finalP1 = neckXSectioner.findNeckMinAlongLine( searchP0, searchP1, 0.05 )
        # fine search
        # planePoints2, finalP2, finalN = neckXSectioner.findNeckMin( finalP1, initialNeckAxis.a )

        
        neckEP = scipy.vstack([self.EP[self.EPMap[e],:] for e in self.neckElems])

        # neckEPi = []
        # for submesh, elems in self.neckElems:
        #   neckEPi.append( scipy.hstack([self.EPMap[submesh][e] for e in elems]) )
        # neckEP = self.EP[ scipy.hstack( neckEPi ) ]

        planePoints2, finalP2, finalN = neckXSectioner.findNeckMin( neckEP.mean(0), initialNeckAxis.a )
        # planePoints2, finalP2, finalN = neckXSectioner.findNeckMin( neckEP.mean(0), self.neckAxis.a )
        # planePoints2 = []
        
        if len(planePoints2)==0:
            neckPlanePoints = planePoints1
            neckPlaneNormal = scipy.array(initialNeckAxis.a)
            neckPlaneOrigin = finalP1
        else:
            neckPlanePoints = planePoints2
            neckPlaneNormal = finalN
            neckPlaneOrigin = finalP2

        return neckPlaneNormal, neckPlaneOrigin, neckPlanePoints, searchP0, searchP1
            
    def calcXYAxes( self, alignGF=False ):
        """ calculate X and Y axes normal to the shaft axis
        """

        nodeCoords = self.GF.get_field_parameters()[:,:,0].T
        # lateral epicondyle - medial epicondyle
        v1 = nodeCoords[self.condyleAlignmentNodes[0]] - nodeCoords[self.condyleAlignmentNodes[1]]
        
        z = self.shaftAxis.a
        y = _norm( scipy.cross( _norm(v1), _norm(z) ) )
        x = _norm( scipy.cross( y, _norm(z) ) )
        CoM = self.GF.calc_CoM_2D([5,5])
        
        self.xAxis = FT.Line3D( x, CoM )
        self.yAxis = FT.Line3D( y, CoM )
        
        if alignGF:
            pAxes = scipy.array( [x,y,z] )
            
            targetCoM = scipy.zeros(3)
            targetPAxes = scipy.array([[1,0,0],[0,1,0],[0,0,1]], dtype=float)
            
            T = alignment.calcAffine( (CoM, pAxes), ( targetCoM, targetPAxes ) )
            self.GF.transformAffine( scipy.vstack((T,scipy.ones(4))) )
            self.shaftAligned = True

    def calcHeadDiameter( self ):
        headEP = scipy.vstack([self.EP[self.EPMap[e],:] for e in self.headElems])
        rms, [cx,cy,cz,r] = FT.fitSphere( headEP )
        m = measurement( 'head_diameter', r*2.0 )
        m.centre = scipy.array([cx, cy, cz])
        self.measurements['head_diameter'] = m
        return m.value
    
    def calcNeckShaftAngle( self ):
        """
        angle between shaftAxis and neckAxis
        """
        nom = scipy.dot( self.neckAxis.a, self.shaftAxis.a )
        denom = _mag( self.neckAxis.a ) * _mag( self.shaftAxis.a )
        a = scipy.rad2deg(scipy.arccos(nom/denom))
        if a<90.0:
        	a = 180.0 - a
        m = measurement('neck_shaft_angle', a)
        self.measurements['neck_shaft_angle'] = m
        return m.value
    
    def calcNeckWidth( self ):
        
        # use neck EPs
        neckEP = scipy.vstack([self.EP[self.EPMap[e],:] for e in self.neckLongElems])
        
        # neckEPi = []
        # for submesh, elems in self.neckLongElems:
        #   neckEPi.append( scipy.hstack([self.EPMap[submesh][e] for e in elems]) )         
        # neckEP = self.EP[ scipy.hstack( neckEPi ) ]
        
        axisPoints = []
        axisPointsT = []
        for d in neckEP:
            p,t = self.neckAxis.findClosest(d)
            axisPoints.append(p)
            axisPointsT.append(t)
            
        #~ axisPoints = scipy.array([self.neckAxis.findClosest(d)[0] for d in neckEP])
        initPI = scipy.sqrt(((neckEP - axisPoints)**2.0).sum(1)).argmin()
        initP = axisPoints[ initPI ]
        initT = axisPointsT[ initPI ]
        searchP0 = self.neckAxis.eval( min(axisPointsT)+5.0 )
        searchP1 = self.neckAxis.eval( max(axisPointsT)-5.0 )
        
        neckXSectioner = NeckMinXSectionFit( neckEP )
        # initial search constrained to neck axis
        planePoints1, finalP1 = neckXSectioner.findNeckMinAlongLine( searchP0, searchP1, 0.05 )
        # fine search - doesnt work well
        planePoints2, finalP2, finalN = neckXSectioner.findNeckMin( finalP1, self.neckAxis.a )
        # planePoints2, finalP2, finalN = neckXSectioner.findNeckMin( neckEP.mean(0), self.neckAxis.a )
        # planePoints2 = []
        
        if len(planePoints2)==0:
            planePoints = planePoints1
            finalP = finalP1
        else:
            planePoints = planePoints2
            finalP = finalP2
        
        planePoints = self.neckAxis._planePoints

        # diameter = distance between highest and lowest plane points
        neckSup = planePoints[ planePoints[:,2].argmax() ]
        neckInf = planePoints[ planePoints[:,2].argmin() ]
        neckW = scipy.sqrt( ((neckSup - neckInf)**2.0).sum() )
        
        m = measurement( 'neck_width', neckW )
        m.centre = self.neckAxis.b
        m.interceptSup = neckSup
        m.interceptInf = neckInf
        m.searchMin = self.neckAxis._searchP0
        m.searchMax = self.neckAxis._searchP1
        self.measurements[ 'neck_width' ] = m
        return m.value
    
    def calcNeckDiameter( self ):
        """
        as a first approximation, finds the closest EP in neck elements
        to the neck axis
        """
        neckEP = scipy.vstack([self.EP[self.EPMap[e],:] for e in self.neckElems])
        # neckEPi = []
        # for submesh, elems in self.neckElems:
        #   neckEPi.append( scipy.hstack([self.EPMap[submesh][e] for e in elems]) )
        # neckEP = self.EP[ scipy.hstack( neckEPi ) ]
        
        axisPoints = scipy.array([self.neckAxis.findClosest(d)[0] for d in neckEP])
        dMin = scipy.sqrt(((neckEP - axisPoints)**2.0).sum(1)).min()
        
        m = measurement( 'neck_diameter', dMin*2.0 )
        self.measurements[ 'neck_diameter' ] = m
        return m.value
        
    def calcFemoralAxisLength( self ):
        
        # calc closest point to axis in head
        headEPi = scipy.hstack([self.EPMap[e] for e in self.headElems])
        headEP = scipy.vstack([self.EP[self.EPMap[e],:] for e in self.headElems])
        # headEPi = scipy.hstack(self.EPMap[ self.regionName2ElementMap['head']])
        # headEP = self.EP[ headEPi ]
        
        headD = scipy.array([ self.neckAxis.calcDistanceFromPoint( p ) for p in headEP ])
        # headIntercept = headEP[ headD.argmin() ]
        headIntercept = self.neckAxis.findClosest( headEP[ headD.argmin() ] )[0]
        
        # calc closest point to axis in gtroc
        gTrocEPi = scipy.hstack([self.EPMap[e] for e in self.gTrocElems])
        gTrocEP = scipy.vstack([self.EP[self.EPMap[e],:] for e in self.gTrocElems])
        # gTrocEPi = scipy.hstack(self.EPMap[ self.regionName2ElementMap['greatertrochanter']])
        # gTrocEP = self.EP[ gTrocEPi ]
        
        gTrocD = scipy.array([ self.neckAxis.calcDistanceFromPoint( p ) for p in gTrocEP ])
        # gTrocIntercept = gTrocEP[ gTrocD.argmin() ]
        gTrocIntercept = self.neckAxis.findClosest( gTrocEP[ gTrocD.argmin() ] )[0]
        
        m = measurement( 'femoral_axis_length', scipy.sqrt( ((headIntercept - gTrocIntercept)**2.0).sum() ) )
        m.headIntercept = ( headEPi[headD.argmin()], scipy.array( headIntercept ) )
        m.gTrocIntercept = ( gTrocEPi[gTrocD.argmin()], scipy.array( gTrocIntercept ) )
        self.measurements['femoral_axis_length'] = m
        return m.value
        
    def calcFemoralLengthLong( self ):
        
        headEPi = scipy.hstack([self.EPMap[e] for e in self.headElems])
        headEP = scipy.vstack([self.EP[self.EPMap[e],:] for e in self.headElems])
        mCondEPi = scipy.hstack([self.EPMap[e] for e in self.mCondElems])
        mCondEP = scipy.vstack([self.EP[self.EPMap[e],:] for e in self.mCondElems])

        # headEPi = scipy.hstack(self.EPMap[ self.regionName2ElementMap['head']])
        # headEP = self.EP[ headEPi ]
        # mCondEPi = scipy.hstack(self.EPMap[ self.regionName2ElementMap['medialcondyle']])
        # mCondEP = self.EP[ mCondEPi ]
        
        # project eps on the shaft axis
        headP = self.shaftAxis.findClosest( headEP )[1]
        mCondP = self.shaftAxis.findClosest( mCondEP )[1]
        
        if headP.mean()>mCondP.mean():
            length = headP.max() - mCondP.min()
            headPointI = headEPi[headP.argmax()]
            headPoint = headEP[headP.argmax()]
            mCondPointI = mCondEPi[mCondP.argmin()]
            mCondPoint = mCondEP[mCondP.argmin()]
        else:
            length = mCondP.max() - headP.min()
            headPointI = headEPi[headP.argmin()]
            headPoint = headEP[headP.argmin()]
            mCondPointI = mCondEPi[mCondP.argmax()]
            mCondPoint = mCondEP[mCondP.argmax()]

        # length = max(abs(mCondP.max() - headP.min()), abs(headP.max() - mCondP.min())) 
        
        # headPointI = headEPi[ headP.argmin() ]
        # headPoint = headEP[ headP.argmin() ]
        # mCondPointI = mCondEPi[ mCondP.argmin() ]
        # mCondPoint = mCondEP[ mCondP.argmin() ]
        
        m = measurement( 'length_long', length )
        m.headPoint = ( headPointI, headPoint )
        m.mCondPoint = ( mCondPointI, mCondPoint )
        self.measurements['length_long'] = m
        return m.value
    
    def calcFemoralLengthShort( self ):
        
        gTrocEPi = scipy.hstack([self.EPMap[e] for e in self.gTrocElems])
        gTrocEP = scipy.vstack([self.EP[self.EPMap[e],:] for e in self.gTrocElems])
        lCondEPi = scipy.hstack([self.EPMap[e] for e in self.lCondElems])
        lCondEP = scipy.vstack([self.EP[self.EPMap[e],:] for e in self.lCondElems])

        # gTrocEPi = scipy.hstack(self.EPMap[ self.regionName2ElementMap['greatertrochanter']])
        # gTrocEP = self.EP[ gTrocEPi ]
        # lCondEPi = scipy.hstack(self.EPMap[ self.regionName2ElementMap['lateralcondyle']])
        # lCondEP = self.EP[ lCondEPi ]
        
        # project eps on the shaft axis
        gTrocP = self.shaftAxis.findClosest( gTrocEP )[1]
        lCondP = self.shaftAxis.findClosest( lCondEP )[1]

        if gTrocP.mean()>lCondP.mean():
            length = gTrocP.max() - lCondP.min()
            gTrocPointI = gTrocEPi[gTrocP.argmax()]
            gTrocPoint = gTrocEP[gTrocP.argmax()]
            lCondPointI = lCondEPi[lCondP.argmin()]
            lCondPoint = lCondEP[lCondP.argmin()]
        else:
            length = lCondP.max() - gTrocP.min()
            gTrocPointI = gTrocEPi[gTrocP.argmin()]
            gTrocPoint = gTrocEP[gTrocP.argmin()]
            lCondPointI = lCondEPi[lCondP.argmax()]
            lCondPoint = lCondEP[lCondP.argmax()]
        
        # length = lCondP.max() - gTrocP.min() 
        
        # gTrocPointI = gTrocEPi[ gTrocP.argmin() ]
        # gTrocPoint = gTrocEP[ gTrocP.argmin() ]
        # lCondPointI = lCondEPi[ lCondP.argmin() ]
        # lCondPoint = lCondEP[ lCondP.argmin() ]
        
        m = measurement( 'length_short', length )
        m.gTrocPoint = ( gTrocPointI, gTrocPoint )
        m.lCondPoint = ( lCondPointI, lCondPoint )
        self.measurements['length_short'] = m
        return m.value
        
    def calcSubtrochantericDiameter( self ):
        # get line mesh of element boundary below the ltroc
        S = self.GF.makeLineElementsFromPoints( self.subtrochanterNodes, 5, {'line5l':'line_L4'} )
        
        # evaluate line EPs
        ep = S.evaluate_geometric_field( [self.epD[0]] ).T
        
        # fit sphere to EPs
        rms, [cx,cy,cz,r] = FT.fitSphere( ep )
        
        m = measurement( 'subtrochanteric_diameter', r*2.0 )
        m.centre = scipy.array( [cx, cy, cz] )
        self.measurements['subtrochanteric_diameter'] = m
        return m.value
        
    def calcSubtrochantericWidth( self ):
        # get line mesh of element boundary below the ltroc
        S = self.GF.makeLineElementsFromPoints( self.subtrochanterNodes, 5, {'line5l':'line_L4'} )
        
        # evaluate line EPs
        ep = S.evaluate_geometric_field( [self.epD[0]] ).T

        # project on x axis
        projX = self.xAxis.findClosest( ep )[1]
        width = projX.max() - projX.min()
        
        m = measurement( 'subtrochanteric_width', width )
        m.p1 = scipy.array( ep[ projX.argmax() ] )
        m.p2 = scipy.array( ep[ projX.argmin() ] )
        self.measurements['subtrochanteric_width'] = m
        return m.value
        
    def calcMidshaftDiameter( self ):
        # get line mesh of element boundary below the ltroc
        S = self.GF.makeLineElementsFromPoints( self.midshaftNodes, 5, {'line5l':'line_L4'} )
        
        # evaluate line EPs
        ep = S.evaluate_geometric_field( [self.epD[0]] ).T
        
        # fit sphere to EPs
        rms, [cx,cy,cz,r] = FT.fitSphere( ep )
        
        m = measurement( 'midshaft_diameter', r*2.0 )
        m.centre = scipy.array( [cx, cy, cz] )
        self.measurements['midshaft_diameter'] = m
        return m.value
    
    def calcMidshaftWidth( self ):
        """
        a caliper-like measurement
        """

        # get line mesh of element boundary below the ltroc
        S = self.GF.makeLineElementsFromPoints( self.midshaftNodes, 5, {'line5l':'line_L4'} )
        
        # evaluate line EPs
        ep = S.evaluate_geometric_field( [self.epD[0]] ).T
        
        # project on x axis
        projX = self.xAxis.findClosest( ep )[1]
        width = projX.max() - projX.min()
        
        m = measurement( 'midshaft_width', width )
        m.p1 = scipy.array( ep[ projX.argmax() ] )
        m.p2 = scipy.array( ep[ projX.argmin() ] )
        self.measurements['midshaft_width'] = m
        return m.value
    
    def calcEpicondylarWidth( self ):
        #~ return self.calcEpicondylarWidthCaliper()
        #~ return self.calcEpicondylarWidthFitBox()
        return self.calcEpicondylarWidthByNode()
                
    def calcEpicondylarWidthCaliper( self ):
        """Finds the most distant points between the medial and lateral epicondyles in the
        shaft-XY CS.
        """
        if not self.shaftAligned:
            self.alignToShaftAxis()

        lCondEPi = scipy.hstack([self.EPMap[e] for e in self.lCondElems])
        lCondEP = scipy.vstack([self.EP[self.EPMap[e],:] for e in self.lCondElems])
        mCondEPi = scipy.hstack([self.EPMap[e] for e in self.mCondElems])
        mCondEP = scipy.vstack([self.EP[self.EPMap[e],:] for e in self.mCondElems])
        
        # mCondEPi = []
        # for submesh, elems in self.epiCondylarMCondElems:
        #   mCondEPi.append( scipy.hstack([self.EPMap[submesh][e] for e in elems]) )
        
        # mCondEPi = scipy.hstack( mCondEPi )
        # mCondEP = self.EP[ mCondEPi ]
        medialPoint = mCondEP[ mCondEP[:,0].argmax() ]
        
        # lCondEPi = []
        # for submesh, elems in self.epiCondylarLCondElems:
        #   lCondEPi.append( scipy.hstack([self.EPMap[submesh][e] for e in elems]) )
        
        # lCondEPi = scipy.hstack( lCondEPi )
        # lCondEP = self.EP[ lCondEPi ] 
        lateralPoint = lCondEP[ lCondEP[:,0].argmin() ]
        
        m = measurement( 'epicondylar_width', scipy.sqrt(((medialPoint - lateralPoint)**2.0).sum()) )
        m.p1 = ( mCondEPi[ mCondEP[:,0].argmax() ], scipy.array( medialPoint ) )
        m.p2 = ( lCondEPi[ lCondEP[:,0].argmin() ], scipy.array( lateralPoint ) )
        self.measurements['epicondylar_width'] = m
        
        self.epicondylarAxis = FT.Line3D( medialPoint-lateralPoint, (medialPoint+lateralPoint)/2.0 )
        mAxis = measurement( 'epicondylar_axis', ( self.epicondylarAxis.a, self.epicondylarAxis.b ) )
        self.measurements['epicondylar_axis'] = mAxis
        
        return m.value
        
    def calcEpicondylarWidthDistance( self ):
        """
        finds most distant pair of points in predefined condyle elements in the shaft-XY CS.
        """
        if not self.shaftAligned:
            self.alignToShaftAxis()
        
        lCondEPi = scipy.hstack([self.EPMap[e] for e in self.lCondElems])
        lCondEP = scipy.vstack([self.EP[self.EPMap[e],:] for e in self.lCondElems])
        mCondEPi = scipy.hstack([self.EPMap[e] for e in self.mCondElems])
        mCondEP = scipy.vstack([self.EP[self.EPMap[e],:] for e in self.mCondElems])

        # mCondEPi = []
        # for submesh, elems in self.epiCondylarMCondElems:
        #   mCondEPi.append( scipy.hstack([self.EPMap[submesh][e] for e in elems]) )
        
        # mCondEPi = scipy.hstack( mCondEPi )
        # mCondEP = self.EP[ mCondEPi ]
        
        # lCondEPi = []
        # for submesh, elems in self.epiCondylarLCondElems:
        #   lCondEPi.append( scipy.hstack([self.EPMap[submesh][e] for e in elems]) )
        
        # lCondEPi = scipy.hstack( lCondEPi )
        # lCondEP = self.EP[ lCondEPi ] 
        
        d = distance.cdist( mCondEP, lCondEP, 'euclidean' )
        dInd = d.argmax()
        
        di = scipy.floor( dInd / d.shape[0] )
        dj = dInd - di*d.shape[0]
        
        medialPoint = mCondEP[ di ]
        medialPointI = mCondEPi[ di ]
        lateralPoint = lCondEP[ dj ]
        lateralPointI = lCondEPi[ dj ]
        
        m = measurement( 'epicondylar_width', scipy.sqrt(((medialPoint - lateralPoint)**2.0).sum()) )
        m.p1 = ( medialPointI, scipy.array( medialPoint ) )
        m.p2 = ( lateralPointI, scipy.array( lateralPoint ) )
        self.measurements['epicondylar_width'] = m
        
        self.epicondylarAxis = FT.Line3D( medialPoint-lateralPoint, (medialPoint+lateralPoint)/2.0 )
        mAxis = measurement( 'epicondylar_axis', ( self.epicondylarAxis.a, self.epicondylarAxis.b ) )
        self.measurements['epicondylar_axis'] = mAxis
        
        return m.value

    def calcEpicondylarWidthFitBox( self ):
        """
        fit box to condyle EPs, find max width
        """
        
        def getMinMax( data ):
            dataMin = ( data.argmin(), data.min() )
            dataMax = ( data.argmax(), data.max() )
            return dataMin, dataMax
        
        if not self.shaftAligned:
            self.alignToShaftAxis()

        lCondEPi = scipy.hstack([self.EPMap[e] for e in self.lCondElems])
        lCondEP = scipy.vstack([self.EP[self.EPMap[e],:] for e in self.lCondElems])
        mCondEPi = scipy.hstack([self.EPMap[e] for e in self.mCondElems])
        mCondEP = scipy.vstack([self.EP[self.EPMap[e],:] for e in self.mCondElems])
        
        # mCondEPi = scipy.hstack(self.EPMap[ self.regionName2ElementMap['medialcondyle']])
        # mCondEP = self.EP[ mCondEPi ]
        
        # lCondEPi = scipy.hstack(self.EPMap[ self.regionName2ElementMap['lateralcondyle']])
        # lCondEP = self.EP[ lCondEPi ]
        
        # condyleEPi = scipy.hstack( [lCondEPi, mCondEPi] )
        condyleEP = scipy.vstack( [lCondEP, mCondEP] )
        
        boxCoM0 = condyleEP.mean(0)
        boxCoMOpt, boxVolume, boxDim, boxAxes = FT.fitBox( condyleEP, boxCoM0, [self.xAxis.a, self.yAxis.a, self.shaftAxis.a] )

        # find the widest axis
        p = [ axis.findClosest( condyleEP )[1] for axis in boxAxes ]
        
        pMin, pMax = getMinMax( p[ boxDim.argmax() ] )
            
        epMin = condyleEP[ pMin[0] ]
        epMax = condyleEP[ pMax[0] ]

        width = scipy.sqrt(((epMin - epMax )**2.0).sum())

        m = measurement( 'epicondylar_width', width )
        m.p1 = ( condyleEPi[pMin[0]], condyleEP[ pMin[0] ] )
        m.p2 = ( condyleEPi[pMax[0]], condyleEP[ pMax[0] ] )
        self.measurements['epicondylar_width'] = m
        
        self.epicondylarAxis = boxAxes[ boxDim.argmax() ]
        
        return m.value

    def calcEpicondylarWidthByNode( self ):
        nodeCoords = self.GF.get_field_parameters()[:,:,0].T
        p1 = nodeCoords[self.condyleAlignmentNodes[0]] 
        p2 = nodeCoords[self.condyleAlignmentNodes[1]]
        
        width = scipy.sqrt((( p1 - p2 )**2.0).sum())
        
        epTree = FT.cKDTree( self.EP )
        p1i, p2i = epTree.query( [p1,p2] )[1]
        
        m = measurement( 'epicondylar_width', width )
        m.p1 = ( p1i, p1 )
        m.p2 = ( p2i, p2 )
        self.measurements['epicondylar_width'] = m
        
        self.epicondylarAxis = FT.Line3D( p1-p2, (p1+p2)/2.0 )
        
        return m.value

    def calcAnteversion(self):

        # create plane normal to shaft axis
        px = FT.norm(scipy.cross(self.shaftAxis.a, [1,0,0]))
        py = FT.norm(scipy.cross(self.shaftAxis.a, px))
        p = FT.Plane([0,0,0], self.shaftAxis.a, x=px, y=py)
        # p = FT.Plane([0,0,0], [0,0,1], x=[1,0,0], y=[0,1,0])

        # get vector connecting posterior extremes of condyles
        lcEP = scipy.vstack([self.EP[self.EPMap[e],:] for e in self.lCondElems])
        mcEP = scipy.vstack([self.EP[self.EPMap[e],:] for e in self.mCondElems])
        # mcEP = self.EP[scipy.hstack(self.EPMap[self.regionName2ElementMap['medialcondyle']])]
        # lcEP = self.EP[scipy.hstack(self.EPMap[self.regionName2ElementMap['lateralcondyle']])]
        mcCP, mcCT = self.yAxis.findClosest(mcEP)
        mcP = mcEP[mcCT.argmin(),:]
        lcCP, lcCT = self.yAxis.findClosest(lcEP)
        lcP = lcEP[lcCT.argmin(),:]
        cv = FT.Line3D(mcP-lcP, (mcP+lcP)/2.0)

        # project this vector on plane
        v1 = p.project2Plane2D(cv.a)
        # v1 = scipy.array([1,0,0], dtype=float)

        # project neckAxis onto plane
        v2 = p.project2Plane2D(self.neckAxis.a)
        
        # calculate angle between projections
        a = FT.angle(v1,v2)*180.0/scipy.pi
        if a>90.0:
        	a = 180.0 - a

        m = measurement( 'anteversion_angle', a )
        m.posteriorCondyleAxis = cv
        self.measurements['anteversion_angle'] = m

        return m.value

#=========================================================================#
class NeckMinXSectionFit( object ):
    """ Class for finding the plane cutting the narrowest part of the
    femoral neck.
    
    Initial guess should be a point halfway between CoM of the proximal 
    shaft and the head centre.
    
    Returns the slice, slice CoM and normal
    """
    
    thresh = 50.0
    planeProjectionMaxDist = 1.0
    useCallback = False
    
    def __init__( self, neckDataPoints ):
        
        self.dataPoints = neckDataPoints
        self.k = 0
        
    def findNeckMin( self, initialP, initialN ):
        """ find CoM and normal of xsection slice at with minimum area
        in neck
        """
        
        #~ self.a = scipy.array( initialN )
        #~ self.b = scipy.array( initialP )
        #~ t = 0
        self.k = -1
        xtol=1e-6
        ftol=1e-6
        
        #~ X0 = scipy.hstack( (t, initialN) )
        print('initP:', initialP)
        print('initN:', initialN)
        
        X0 = scipy.hstack( (initialP, initialN) )
        print('\nStarting neck search')
        if self.useCallback:
            callBack = self._callback
            callback( X0 )
        else:
            callBack = None
        
        X = fmin( self._objFunc, X0, callback=callBack, xtol=xtol, ftol=ftol  )
        # X = leastsq( self._objFuncLeastsq, X0, xtol=xtol, ftol=ftol  )[0]
        
        finalP = X[:3]
        finalN = X[3:]
        
        finalPlane = _Plane3D( NP=[finalN,finalP] )
        planePoints, planePointsI = finalPlane.projectClosePoints( self.dataPoints, self.planeProjectionMaxDist )
        
        return (planePoints, finalP, finalN)
    
    def findNeckMinFixedN( self, initialP, N ):
        """ find CoM of xsection slice at with minimum area
        in neck, given a fixed normal direction
        """
        
        #~ self.a = scipy.array( initialN )
        #~ self.b = scipy.array( initialP )
        #~ t = 0
        self.k = -1
        xtol=0.001
        ftol=0.001
        
        #~ X0 = scipy.hstack( (t, initialN) )
        self.planeN = N
        X0 = initialP
        print('\nStarting neck search')
        if self.useCallback:
            callBack = self._callbackFixedN
            callback( X0 )
        else:
            callBack = None
        
        finalP = fmin( self._objFuncFixedN, X0, callback=callBack, xtol=xtol, ftol=ftol  )
        
        finalPlane = _Plane3D( NP=[self.planeN,finalP] )
        planePoints, planePointsI = finalPlane.projectClosePoints( self.dataPoints, self.planeProjectionMaxDist )
        
        return planePoints, finalP
    
    def findNeckMinAlongLineOld( self, p0, p1, res ):
        """ find CoM and normal of xsection slice at with minimum area
        in neck
        """
        
        T = scipy.arange(0.0, 1.0+res, res)
        line = FT.lineElement3D( p0, p1 )
        self.planeN = _norm(p1 - p0)
        
        areas = []
        for t in T:
            P = line.eval(t)
            areas.append( self._objFuncFixedN( P ) )
            
        bestT = T[ scipy.argmin(areas) ]
        bestP = line.eval(bestT)
        
        finalPlane = _Plane3D( NP=[self.planeN,bestP] )
        planePoints, planePointsI = finalPlane.projectClosePoints( self.dataPoints, self.planeProjectionMaxDist )
        
        return planePoints, bestP
    
    def findNeckMinAlongLine( self, p0, p1, res ):
        """ find CoM and normal of xsection slice at with minimum area
        in neck
        """
        
        line = FT.LineElement3D( p0, p1 )
        self.planeN = _norm(p1 - p0)

        def obj( x ):
            P = line.eval(x)
            return self._objFuncFixedN( P )
        
        xtol=0.001
        ftol=0.001
        bestT = fmin( obj, 0.5, xtol=xtol, ftol=ftol  )
        bestP = line.eval(bestT)
        
        finalPlane = _Plane3D( NP=[self.planeN,bestP] )
        planePoints, planePointsI = finalPlane.projectClosePoints( self.dataPoints, self.planeProjectionMaxDist )
        
        return planePoints, bestP
    
    def _lineFunc( self, t ):
        return  t*self.a + self.b
        
    def _objFunc( self, X ):
        """ objective function to minimise: the xsectional area of the
        slice.
        """
        #~ P = self._lineFunc( X[0] )
        #~ N = X[1:]
        P = X[:3]
        N = X[3:]
        plane = _Plane3D( NP=[N,P] )
        planePoints, planePointsI = plane.projectClosePoints( self.dataPoints, self.planeProjectionMaxDist )
        
        
        method = 4
        if method == 1:
            # Approximate using sum of pixel radii from CoM
            
            sr = slice.calculateSlicePolar()[0]
            area = sr.sum()
        
        elif method == 2:
            # approximate area as area of ellipse characterised by principal
            # axes
            [b,a]=slice.calculatePrincipalAxes ()[1]
            area = scipy.pi * (0.5*scipy.sqrt(a)) * (0.5*scipy.sqrt(b))
        
        elif method == 3:
            # Approximate using sum of pixel radii from centre of slice
            
            # get all nonzero point indices
            xi, yi = scipy.nonzero( slice.I )
            # center
            xi = xi - slice.I.shape[0]/2.0
            yi = yi - slice.I.shape[1]/2.0
            nTotal = xi.shape[0]
            
            # convert all point indices into radial coordinates
            R = scipy.sqrt( xi**2.0 + yi**2.0 )
            area = R.sum()
            
        elif method == 4:
            # approximate using sum of datapoint distance from plane P
            R = scipy.sqrt(((planePoints - plane.P)**2.0).sum(1))
            # area = R.sum()
            area = R.mean() + scipy.sqrt(((plane.P - planePoints.mean(0))**2.0).sum())
            # print 'plane points:', len(planePoints), 'area:', area    
            
        
        return area

    def _objFuncLeastsq( self, X ):
        """ objective function to minimise: the xsectional area of the
        slice.
        """
        #~ P = self._lineFunc( X[0] )
        #~ N = X[1:]
        P = X[:3]
        N = X[3:]
        plane = _Plane3D( NP=[N,P] )
        planePoints, planePointsI = plane.projectClosePoints( self.dataPoints, self.planeProjectionMaxDist )
        
        
        method = 4
        if method == 1:
            # Approximate using sum of pixel radii from CoM
            
            sr = slice.calculateSlicePolar()[0]
            area = sr.sum()
        
        elif method == 2:
            # approximate area as area of ellipse characterised by principal
            # axes
            [b,a]=slice.calculatePrincipalAxes ()[1]
            area = scipy.pi * (0.5*scipy.sqrt(a)) * (0.5*scipy.sqrt(b))
        
        elif method == 3:
            # Approximate using sum of pixel radii from centre of slice
            
            # get all nonzero point indices
            xi, yi = scipy.nonzero( slice.I )
            # center
            xi = xi - slice.I.shape[0]/2.0
            yi = yi - slice.I.shape[1]/2.0
            nTotal = xi.shape[0]
            
            # convert all point indices into radial coordinates
            R = scipy.sqrt( xi**2.0 + yi**2.0 )
            area = R.sum()
            
        elif method == 4:
            # approximate using sum of datapoint distance from plane P
            R = scipy.sqrt(((planePoints - plane.P)**2.0).sum(1))
            # area = R.sum()
            # area = R.mean() + scipy.sqrt(((plane.P - planePoints.mean(0))**2.0).sum())
            # print 'plane points:', len(planePoints), 'area:', area    
            area = scipy.hstack([R, scipy.sqrt(((plane.P - planePoints.mean(0))**2.0).sum())])
        
        return area
    
    def _objFuncFixedN( self, X ):
        """ objective function to minimise: the xsectional area of the
        slice.
        """
        #~ P = self._lineFunc( X[0] )
        #~ N = X[1:]
        P = X[:3]
        plane = _Plane3D( NP=[self.planeN,P] )
        planePoints, planePointsI = plane.projectClosePoints( self.dataPoints, self.planeProjectionMaxDist )
        
        method = 4
        if method == 1:
            # Approximate using sum of pixel radii from CoM
            
            sr = slice.calculateSlicePolar()[0]
            area = sr.sum()
        
        elif method == 2:
            # approximate area as area of ellipse characterised by principal
            # axes
            [b,a]=slice.calculatePrincipalAxes ()[1]
            area = scipy.pi * (0.5*scipy.sqrt(a)) * (0.5*scipy.sqrt(b))
        
        elif method == 3:
            # Approximate using sum of pixel radii from centre of slice
            
            # get all nonzero point indices
            xi, yi = scipy.nonzero( slice.I )
            # center
            xi = xi - slice.I.shape[0]/2.0
            yi = yi - slice.I.shape[1]/2.0
            nTotal = xi.shape[0]
            
            # convert all point indices into radial coordinates
            R = scipy.sqrt( xi**2.0 + yi**2.0 )
            area = R.sum()
            
        elif method == 4:
            # approximate using sum of datapoint distance from plane P
            R = scipy.sqrt(((planePoints - plane.P)**2.0).sum(1))
            area = R.sum()
            print(area)
            #~ pdb.set_trace()
            
        return area
        
    def _callback( self, xk ):
        self.k += 1
        print('it:%(it)i\t P: %(p0)8.5f %(p1)8.5f %(p2)8.5f\t N: %(n0)8.5f %(n1)8.5f %(n2)8.5f'\
              %{'it':self.k, 'p0':xk[0], 'p1':xk[1], 'p2':xk[2], 'n0':xk[3], 'n1':xk[4], 'n2':xk[5]})

    def _callbackFixedN( self, xk ):
        self.k += 1
        print('it:%(it)i\t P: %(p0)8.5f %(p1)8.5f %(p2)8.5f\t'\
              %{'it':self.k, 'p0':xk[0], 'p1':xk[1], 'p2':xk[2]})

    #~ def calculateSlicePolar( self ):
        #~ """ calculates the polar coordinates of pixels in slice centered
        #~ at the CoM
        #~ """
    #~ 
        #~ # get all nonzero point indices
        #~ xi, yi = scipy.nonzero( self.I )
        #~ # center
        #~ xi = xi - self.CoM[0]
        #~ yi = yi - self.CoM[1]
        #~ nTotal = xi.shape[0]
        #~ 
        #~ # convert all point indices into radial coordinates
        #~ R = scipy.sqrt( xi**2.0 + yi**2.0 )      # r
        #~ theta = scipy.zeros( nTotal )
        #~ for i in range( nTotal ):
            #~ theta[i] = calcTheta( xi[i], yi[i] ) # theta
        #~ 
        #~ self.theta = theta
        #~ self.R = R
            #~ 
        #~ return ( self.R, self.theta )    
