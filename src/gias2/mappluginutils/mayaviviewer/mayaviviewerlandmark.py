"""
FILE: mayaviviewerlandmark.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: Container class for landmark points in a mayavi scene.

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import numpy as np

from gias2.mappluginutils.mayaviviewer import MayaviViewerSceneObject, MayaviViewerObject


class MayaviViewerLandmarkSceneObject(MayaviViewerSceneObject):
    typeName = 'landmark'

    def __init__(self, name, sceneObject=None):
        self.name = name
        if sceneObject == None:
            self.sceneObject = {}
        else:
            self.sceneObject = sceneObject

    def addSceneObject(self, name, obj):
        self.sceneObject[name] = obj

    def setVisibility(self, visible):
        for obj in list(self.sceneObject.values()):
            obj.visible = visible

    def remove(self):
        for name in list(self.sceneObject.keys()):
            self.sceneObject[name].remove()
            del self.sceneObject[name]


# object for landmarks
class MayaviViewerLandmark(MayaviViewerObject):
    typeName = 'landmark'
    textLineRadius = 0.5
    charWidth = 0.01
    textColour = (1, 1, 1)
    _renderArgs = dict(mode='sphere',
                       scale_factor=2.0,
                       resolution=16,
                       color=(0.0, 1.0, 0.0),
                       opacity=1.0
                       )
    add3DText = False

    def __init__(self, name, landmarkCoords, drawWidthTubes=False, text2d=False, renderArgs=None):
        self.name = name
        self.coords = landmarkCoords
        self._drawWidthTubes = drawWidthTubes
        self._text2d = text2d
        self.sceneObject = None

        if renderArgs is not None:
            self._renderArgs = renderArgs

        self._renderArgs['name'] = self.name

    def setRenderArgs(self, args):
        self._renderArgs = args
        if 'name' not in list(self._renderArgs.keys()):
            self._renderArgs['name'] = self.name

    def setVisibility(self, visible):
        self.sceneObject.setVisibility(visible)

    def remove(self):

        if self.sceneObject:
            self.sceneObject.remove()
            self.sceneObject = None

        self._M = None

    def draw(self, scene):

        # print 'DRAWING landmark '+self.name
        self.sceneObject = MayaviViewerLandmarkSceneObject(self.name)

        # draw axes
        self._drawPoint(scene)
        # add text of landmark
        if self._text2d:
            self._drawText2D(scene)

    def _drawText2D(self, scene):
        tx = 0.02
        ty = 0.02
        tspacing = 0.05
        lString = self.name + ': ' + ', '.join(['{:5.2f}'.format(ci) for ci in self.coords])
        sObj = scene.mlab.text(tx, ty, lString,
                               width=len(lString) * self.charWidth,
                               name='text2d_' + self.name,
                               color=self.textColour)
        self.sceneObject.addSceneObject('text2d_' + self.name, sObj)

    def _addText3D(self, scene, origin, offset):

        # print 'DRAWING 3D TEXT'

        textOrigin = np.array(origin) + np.array(offset)
        textLine = np.array([origin, textOrigin]).T
        lString = self.name + ': ' + ', '.join(['{:5.2f}'.format(ci) for ci in self.coords])
        text = scene.mlab.text(textOrigin[0], textOrigin[1], lString,
                               z=textOrigin[2],
                               width=len(lString) * self.charWidth,
                               name='text3d_' + self.name,
                               color=self.textColour)
        self.sceneObject.addSceneObject('text3d_' + self.name, text)
        line = scene.mlab.plot3d(textLine[0], textLine[1], textLine[2],
                                 tube_radius=self.textLineRadius,
                                 name='text3dline_' + self.name)
        self.sceneObject.addSceneObject('text3dline_' + self.name, line)

    def _drawPoint(self, scene):
        C = self.coords
        name = 'landmark point ' + self.name
        point = scene.mlab.points3d(C[0], C[1], C[2],
                                    **self._renderArgs)
        self.sceneObject.addSceneObject(name, point)
        if self.add3DText:
            self._addText3D(scene, C, [20.0, 0.0, 0.0])

    def updateGeometry(self, coords, scene):
        self.coords = coords
        if self.sceneObject == None:
            self.draw(scene)
        else:
            self.sceneObject.remove()
            self.sceneObject = None
            self.draw(scene)
