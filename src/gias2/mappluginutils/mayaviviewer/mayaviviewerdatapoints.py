"""
FILE: mayaviviewerdatapoints.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: Container class for datapoints in a mayavi scene.

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

from gias2.mappluginutils.mayaviviewer import MayaviViewerSceneObject, MayaviViewerObject, colours


class MayaviViewerDataPointsSceneObject(MayaviViewerSceneObject):
    typeName = 'datapoints'

    def __init__(self, name, points):
        self.name = name
        self.points = points

    def setVisibility(self, visible):
        self.points.visible = visible

    def remove(self):
        self.points.remove()
        self.points = None


class MayaviViewerDataPoints(MayaviViewerObject):
    typeName = 'datapoints'

    def __init__(self, name, coords, scalars=None, renderArgs=None):
        self.name = name
        self.coords = coords

        self.scalarName = 'None'
        if scalars == None:
            self.scalars = {}
        else:
            self.scalars = scalars

        if renderArgs == None:
            self.renderArgs = {}
        else:
            self.renderArgs = renderArgs

        self.sceneObject = None
        self.defaultColour = colours['bone']
        if 'color' not in list(self.renderArgs.keys()):
            self.renderArgs['color'] = self.defaultColour

    def setScalarSelection(self, scalarName):
        self.scalarName = scalarName

    def setVisibility(self, visible):
        self.sceneObject.setVisibility(visible)

    def remove(self):
        self.sceneObject.remove()
        self.sceneObject = None

    def draw(self, scene):
        scene.disable_render = True
        d = self.coords
        s = self.scalars.get(self.scalarName)
        renderArgs = self.renderArgs
        if s != None:
            self.sceneObject = MayaviViewerDataPointsSceneObject(self.name, \
                                                                 scene.mlab.points3d(d[:, 0], d[:, 1], d[:, 2], s,
                                                                                     **renderArgs))
        else:
            self.sceneObject = MayaviViewerDataPointsSceneObject(self.name, \
                                                                 scene.mlab.points3d(d[:, 0], d[:, 1], d[:, 2],
                                                                                     **renderArgs))

        scene.disable_render = False
        return self.sceneObject

    def updateGeometry(self, coords, scene):

        if self.sceneObject == None:
            self.coords = coords
            self.draw(scene)
        else:
            self.sceneObject.points.mlab_source.set(x=coords[:, 0], y=coords[:, 1], z=coords[:, 2])
            self.coords = coords

    def updateScalar(self, scalarName, scene):
        self.setScalarSelection(scalarName)
        if self.sceneObject == None:
            self.draw(scene)
        else:
            d = self.coords
            s = self.scalars.get(self.scalarName)
            renderArgs = self.renderArgs
            if s != None:
                self.sceneObject.points.actor.mapper.scalar_visibility = True
                self.sceneObject.points.mlab_source.reset(x=d[:, 0], y=d[:, 1], z=d[:, 2], s=s)
            else:
                if 'color' not in renderArgs:
                    color = self.defaultColour
                else:
                    color = renderArgs['color']

                self.sceneObject.points.actor.mapper.scalar_visibility = False
                self.sceneObject.points.actor.property.specular_color = color
                self.sceneObject.points.actor.property.diffuse_color = color
                self.sceneObject.points.actor.property.ambient_color = color
                self.sceneObject.points.actor.property.color = color
                self.sceneObject.points.mlab_source.reset(x=d[:, 0], y=d[:, 1], z=d[:, 2])
