"""
FILE: point_picker_test.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: testing a mayavi widget for picking points in 3D

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

################################################################################
# Create some data
import logging

import numpy

log = logging.getLogger(__name__)

elements = 3
height = 10
radius = [10, 11, 12]
theta = numpy.linspace(0.0, 2.0 * numpy.pi, elements * 2 + 1)
zcoords = numpy.linspace(0.0, height, 3)
xparam = []
yparam = []
zparam = []

for row in range(3):
    xcoords = radius[row] * numpy.cos(theta)
    ycoords = radius[row] * numpy.sin(theta)
    for i in range(len(xcoords) - 1):
        xparam.append([xcoords[i]])
        yparam.append([ycoords[i]])
        zparam.append([float(zcoords[row])])

################################################################################
# Plot the data
from mayavi import mlab

# A first plot in 3D
fig = mlab.figure(1)
mlab.clf()
point_plot = mlab.points3d(xparam, yparam, zparam, mode='sphere', scale_factor=1.0, color=(1, 0, 0))
labels = list(range(len(xparam)))
for i in range(len(labels)):
    l = mlab.text(xparam[i][0], yparam[i][0], str(labels[i]), z=zparam[i][0], line_width=0.01,
                  width=0.005 * len(str(labels[i])) ** 1.2)

################################################################################
# Some logic to select 'mesh' and the data index when picking.
# observer for the figure's point picker object
from tvtk.api import tvtk


# ~ def picker_callback(picker_obj, evt):
# ~ picker_obj = tvtk.to_tvtk(picker_obj)
# ~ picked = picker_obj.actors
# ~
# ~ print 'point number', picker_obj.point_id
# ~ print 'point coordinates', picker_obj.pick_position
# ~
# ~ fig.scene.picker.pointpicker.add_observer('EndPickEvent', picker_callback)

################################################################################
# Some logic to pick on click but no move
class MvtPicker(object):
    mouse_mvt = False

    def __init__(self, picker, points_plot):
        self.picker = picker
        self.points_plot = points_plot
        self.picked_points = []

    def on_button_press(self, obj, evt):
        # ~ print 'down'
        self.mouse_mvt = False

    def on_mouse_move(self, obj, evt):
        # ~ print 'move'
        self.mouse_mvt = True

    def on_button_release(self, obj, evt):
        log.debug('release!')
        if not self.mouse_mvt:
            x, y = obj.GetEventPosition()
            self.picker.pick((x, y, 0), fig.scene.renderer)

            # calculate distance to all data points from picked point
            tmp = (self.points_plot.mlab_source.points - self.picker.pick_position) ** 2.0
            # find nearest data point
            nearest_point = tmp.sum(axis=1).argmin()
            log.debug('nearest point:', nearest_point)
            self.picked_points.append(nearest_point)

        self.mouse_mvt = False

    def on_p(self, obj, evt):
        obj = tvtk.to_tvtk(obj)
        picked = obj.actors

        log.debug('point number', obj.point_id)
        log.debug('point coordinates', obj.pick_position)
        # calculate distance to all data points from picked point
        tmp = (self.points_plot.mlab_source.points - obj.pick_position) ** 2.0
        # find nearest data point
        nearest_point = tmp.sum(axis=1).argmin()
        log.debug('nearest point:', nearest_point)
        self.picked_points.append(nearest_point)

        return None


mvt_picker = MvtPicker(fig.scene.picker.pointpicker, point_plot)

fig.scene.interactor.add_observer('LeftButtonPressEvent',
                                  mvt_picker.on_button_press)
fig.scene.interactor.add_observer('MouseMoveEvent',
                                  mvt_picker.on_mouse_move)
fig.scene.interactor.add_observer('LeftButtonReleaseEvent',
                                  mvt_picker.on_button_release)
fig.scene.picker.pointpicker.add_observer('EndPickEvent', mvt_picker.on_p)

mlab.show()
