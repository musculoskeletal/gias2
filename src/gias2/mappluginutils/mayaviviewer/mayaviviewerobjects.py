"""
FILE: mayaviviewerlandmark.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: Generic container class for a mayavi scene.

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import os

os.environ['ETS_TOOLKIT'] = 'qt4'

colours = {'bone': (0.84705882, 0.8, 0.49803922)}


class MayaviViewerObject(object):

    def __init__(self):
        pass

    def draw(self, scene):
        pass

    def setScalarSelection(self, scalarName):
        self.scalarName = scalarName

    def setVisibility(self, visible):
        pass

    def updateGeometry(self, params):
        pass

    def updateScalar(self, scalarName):
        pass

    def remove(self):
        pass


class MayaviViewerSceneObject(object):

    def __init__(self):
        pass


class MayaviViewerObjectsContainer(object):
    """
    stores objects to be rendered in the viewer
    """

    def __init__(self):
        self._objects = {}

    def addObject(self, name, obj):
        # if name in self._objects.keys():
        #     raise ValueError, 'name must be unique'

        if not isinstance(obj, MayaviViewerObject):
            raise TypeError('obj must a MayaviViewerObject')

        self._objects[name] = obj

    def getObjectAll(self, name):
        return self._objects[name]

    def getObjectType(self, name):
        return self._objects[name].typeName

    def getObject(self, name):
        return self._objects[name]

    def getObjectNamesOfType(self, typeName):
        ret = []
        for name, o in list(self._objects.items()):
            if o.typeName == typeName:
                ret.append(name)

        return ret

    def getObjectNames(self):
        return list(self._objects.keys())

    def getNumberOfObjects(self):
        return len(list(self._objects.keys()))

    def removeObject(self, name):
        if name in self._objects:
            o = self._objects[name]
            o.remove()
            del self._objects[name]
        else:
            raise ValueError('Object {} not in container'.format(name))
