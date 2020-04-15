"""
FILE: mayaviscenewidget.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: Class for a mayavi 3D widget.

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import os

os.environ['ETS_TOOLKIT'] = 'qt4'

# from pyface.qt import QtGui, QtCore
from PySide import QtGui

from traits.api import HasTraits, Instance
from traitsui.api import View, Item
from mayavi.core.ui.api import MayaviScene, MlabSceneModel, \
    SceneEditor


################################################################################
# The actual visualization
class Visualisation(HasTraits):
    scene = Instance(MlabSceneModel, ())

    # @on_trait_change('scene.activated')
    # def update_plot(self):
    #     # This function is called when the view is opened. We don't
    #     # populate the scene when the view is not yet open, as some
    #     # VTK features require a GLContext.

    #     # We can do normal mlab calls on the embedded scene.
    #     self.scene.mlab.test_points3d()

    # the layout of the dialog screated
    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=250, width=300, show_label=False),
                resizable=True  # We need this to resize with the parent widget
                )


################################################################################
# The QWidget containing the visualization, this is pure PyQt4 code.
class MayaviSceneWidget(QtGui.QWidget):

    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        layout = QtGui.QVBoxLayout(self)
        # layout.setMargin(0)
        # layout.setSpacing(0)
        self.visualisation = Visualisation()
        # print '############'
        # print self.visualisation.__dict__
        # print '############'

        # If you want to debug, beware that you need to remove the Qt
        # input hook.
        # QtCore.pyqtRemoveInputHook()
        # import pdb ; pdb.set_trace()
        # QtCore.pyqtRestoreInputHook()

        # The edit_traits call will generate the widget to embed.
        self.ui = self.visualisation.edit_traits(parent=self,
                                                 kind='subpanel').control
        # self.ui = self.visualisation.edit_traits(parent=self,
        #                                          kind='modal').control
        layout.addWidget(self.ui)
        self.ui.setParent(self)
        self.setLayout(layout)
