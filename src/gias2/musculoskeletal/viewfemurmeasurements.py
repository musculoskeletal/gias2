"""
FILE: viewfemurmeasurements.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: function to render femur measurements in a Fieldvi instance.

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import numpy

from gias2.fieldwork.field import geometric_field
from gias2.visualisation import fieldvi


def viewMeasurements(M, G, onCloseCallback=None, F=None):
    if F is None:
        F = fieldvi.fieldvi()
        if onCloseCallback == None:
            F.addOnCloseCallback(onCloseCallback)
        F.GFD = [8, 8]
        F.displayGFNodes = False
        F.scene.background = (0, 0, 0)
        F.addGeometricField(G.name, G, geometric_field.makeGeometricFieldEvaluatorSparse(G, F.GFD))
        F._drawGeometricField(G.name)
    # print 'dingdingding'
    # EP = G.evaluate_geometric_field(M.epD).T

    # draw axes
    _drawAxes(F, M)

    # plot head sphere
    _drawHead(F, M)

    # draw femoral axis length intercepts
    _drawFemoralAxisLength(F, M)

    # draw neck width and tube
    _drawNeckWidth(F, M, drawTube=False)

    # draw subtrochanteric width
    _drawSubTrochantericWidth(F, M)

    # draw midshaft width and tube
    _drawMidshaftWidth(F, M, drawTube=False)

    # draw epicondyle intercepts
    _drawEpicondyleWidth(F, M)

    # add text of measurements
    if 0:
        textMeasurements = (
        'head_diameter', 'neck_width', 'neck_shaft_angle', 'femoral_axis_length', 'subtrochanteric_width')
        tx = 0.02
        ty = 0.02
        tspacing = 0.05
        charWidth = 0.01
        for m in textMeasurements:
            value = M.measurements[m].value
            mString = '{m}: {v:5.2f}'.format(m=m, v=value)
            F.scene.mlab.text(tx, ty, mString, width=len(mString) * charWidth)
            ty += tspacing

    F.configure_traits()
    return F


def _addText3D(F, name, value, unit, mOrigin, offset):
    charWidth = 0.01
    lineWidth = 0.2
    textOrigin = numpy.array(mOrigin) + numpy.array(offset)
    textLine = numpy.array([mOrigin, textOrigin]).T
    mStr = '{}: {:5.2f} {}'.format(name, value, unit)
    F.scene.mlab.text(textOrigin[0], textOrigin[1], mStr, z=textOrigin[2], width=len(mStr) * charWidth,
                      name='text_' + name)
    F.scene.mlab.plot3d(textLine[0], textLine[1], textLine[2], tube_radius=lineWidth, name='textline_' + name)


def _drawAxes(F, M):
    saPoints = M.shaftAxis.eval(numpy.array([-300, 300])).T
    F.scene.mlab.plot3d(saPoints[0], saPoints[1], saPoints[2], name='axis_shaft', tube_radius=1.0)
    naPoints = M.neckAxis.eval(numpy.array([-100, 100])).T
    F.scene.mlab.plot3d(naPoints[0], naPoints[1], naPoints[2], name='axis_neck', tube_radius=1.0)
    ecPoints = M.epicondylarAxis.eval(numpy.array([-100, 100])).T
    F.scene.mlab.plot3d(ecPoints[0], ecPoints[1], ecPoints[2], name='axis_epicondylar', tube_radius=1.0)
    pcPoints = M.measurements['anteversion_angle'].posteriorCondyleAxis.eval(numpy.array([-100, 100])).T
    F.scene.mlab.plot3d(pcPoints[0], pcPoints[1], pcPoints[2], name='axis_posteriorcondyle', tube_radius=1.0)


def _drawHead(F, M):
    headM = M.measurements['head_diameter']
    C = headM.centre
    F.scene.mlab.points3d(C[0], C[1], C[2], mode='sphere', scale_factor=headM.value, resolution=16,
                          name='glyph_headSphere', color=(0.0, 1.0, 0.0), opacity=0.3)
    _addText3D(F, 'head diameter', headM.value, 'mm', C, [-50.0, 0, -50])


def _drawNeckWidth(F, M, drawTube=False):
    # width
    NW = M.measurements['neck_width']
    NWC = NW.centre
    NWSup = NW.interceptSup
    NWInf = NW.interceptInf
    # NWMin = NW.searchMin
    # NWMax = NW.searchMax
    # NWPoints = numpy.array([NWSup, NWC, NWInf, NWMin, NWMax]).T
    NWPoints = numpy.array([NWSup, NWInf]).T
    F.scene.mlab.points3d(NWPoints[0], NWPoints[1], NWPoints[2], name='glyph_neckWidthPoints', mode='sphere',
                          scale_factor=5, resolution=16, color=(1.0, 0.0, 0.0))
    NWLinePoints = numpy.array([NWSup, NWInf]).T
    F.scene.mlab.plot3d(NWLinePoints[0], NWLinePoints[1], NWLinePoints[2], name='glyph_neckWidthLine', tube_radius=1.0)
    _addText3D(F, 'neck width', NW.value, 'mm', NWC, [0.0, 0.0, -100])

    # tube
    if drawTube:
        neckRadiusM = M.measurements['neck_width']
        # neckEnds = M.neckAxis.eval(numpy.array([-50,10])).T
        # NW = M.measurements['neck_width']
        # neckEnds = numpy.array([NW.searchMin, NW.searchMax]).T
        neckEnds = M.neckAxis.eval(numpy.array([-30, 20])).T
        F.scene.mlab.plot3d(neckEnds[0], neckEnds[1], neckEnds[2], name='glyph_neckWidthTube',
                            tube_radius=neckRadiusM.value / 2.0, tube_sides=16, color=(0.0, 0.0, 1.0), opacity=0.3)


def _drawFemoralAxisLength(F, M):
    FAL = M.measurements['femoral_axis_length']
    H = FAL.headIntercept[1]
    G = FAL.gTrocIntercept[1]
    F.scene.mlab.points3d([H[0], G[0]], [H[1], G[1]], [H[2], G[2]], name='glyph_FALPoints', mode='sphere',
                          scale_factor=5, resolution=16, color=(1.0, 0.0, 0.0))
    _addText3D(F, 'femoral axis length', FAL.value, 'mm', G, [200.0, 0.0, -150.0])


def _drawNeckShaftAngle(F, M):
    pass
    # NSA = M.measurements['neck_shaft_angle']
    # O = 
    # _addText3D(F, 'neck shaft angle', NSA.value, 'degrees', O, [-100.0,0.0,0.0])


def _drawSubTrochantericWidth(F, M):
    sTW = M.measurements['subtrochanteric_width']
    points = numpy.array([sTW.p1, sTW.p2]).T
    centre = (sTW.p1 + sTW.p2) * 0.5
    F.scene.mlab.points3d(points[0], points[1], points[2], name='glyph_sTWPoints', mode='sphere', scale_factor=5,
                          resolution=16, color=(1.0, 0.0, 0.0))
    F.scene.mlab.plot3d(points[0], points[1], points[2], name='glyph_sTWLine', tube_radius=1.0)
    _addText3D(F, 'subtrochanteric width', sTW.value, 'mm', centre, [-100.0, 0.0, 0.0])


def _drawMidshaftWidth(F, M, drawTube=False):
    mSW = M.measurements['midshaft_width']
    points = numpy.array([mSW.p1, mSW.p2]).T
    centre = (mSW.p1 + mSW.p2) * 0.5
    F.scene.mlab.points3d(points[0], points[1], points[2], name='glyph_midshaftWidthPoints', mode='sphere',
                          scale_factor=5, resolution=16, color=(1.0, 0.0, 0.0))
    F.scene.mlab.plot3d(points[0], points[1], points[2], name='glyph_midshaftWidthLine', tube_radius=1.0)
    _addText3D(F, 'midshaft width', mSW.value, 'mm', centre, [-100.0, 0.0, 0.0])

    # draw midshaft tube
    if drawTube:
        midshaftEnds = M.shaftAxis.eval(numpy.array([-20, 20])).T
        F.scene.mlab.plot3d(midshaftEnds[0], midshaftEnds[1], midshaftEnds[2], name='glyph_midshaftWidthTube',
                            tube_radius=mSW.value / 2.0, tube_sides=16, color=(0.0, 0.0, 1.0), opacity=0.3)


def _drawEpicondyleWidth(F, M):
    ECW = M.measurements['epicondylar_width']
    # l = EP[ ECW.p1[0] ]
    l = ECW.p1[1]
    # m = EP[ ECW.p2[0] ]
    m = ECW.p2[1]
    c = (l + m) * 0.5
    F.scene.mlab.points3d([l[0], m[0]], [l[1], m[1]], [l[2], m[2]], name='glyph_epicondylarWidthPoints', mode='sphere',
                          scale_factor=5, resolution=16, color=(1.0, 0.0, 0.0))
    # F.scene.mlab.plot3d( [l[0], m[0]], [l[1], m[1]], [l[2], m[2]], name='ECW-axis', tube_radius=1.0, tube_sides=16 )
    _addText3D(F, 'epicondylar width', ECW.value, 'mm', c, [-50.0, 0.0, 100.0])
