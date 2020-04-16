"""
FILE: geometric_field_interactor.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: classes and functions to interact with geometric_fields in mayavi2

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""
import logging

import numpy as np
from mayavi import mlab
from tvtk.api import tvtk

from gias2.fieldwork.field.tools import spline_tools

log = logging.getLogger(__name__)


class point_picker(object):
    """observers methods for a mayavi scene interactor. 
    Applied logic to make it record the picked position only
    if the mouse has not moved between lmb press and release.
    picks up to a user defined number of points before notifying a 
    parent class.
    """

    mouse_mvt = False

    def __init__(self, picker, points_plot, figure, parent):
        self.picker = picker
        self.points_plot = points_plot
        self.parent = parent
        self.figure = figure
        self.picked_points = []

    def set_number_of_points_to_pick(self, n):
        """ clears the list of picked points and defines how many points
        to pick before calling the parent's _picking_done method """

        self.picked_points = []
        self.n = n

    def on_button_press(self, obj, evt):
        self.mouse_mvt = False

    def on_mouse_move(self, obj, evt):
        self.mouse_mvt = True

    def on_button_release(self, obj, evt):
        """ if left mouse button has not moved, gets picked position
        from scene pointpicker and calculates the closest data point.
        records this data point until the number of recorded data points
        reaches self.n, where it calls callback """

        if not self.mouse_mvt:
            x, y = obj.GetEventPosition()
            self.picker.pick((x, y, 0), self.figure.scene.renderer)

            # calculate distance to all data points from picked point
            tmp = (self.points_plot.mlab_source.points - self.picker.pick_position) ** 2.0
            # find nearest data point
            nearest_point = tmp.sum(axis=1).argmin()
            log.debug('point:', nearest_point)
            self.picked_points.append(nearest_point)
            log.debug('points picked:', self.picked_points)

            if len(self.picked_points) == self.n:
                log.debug('adding element...')
                self.parent._picking_done()
                self.mouse_mvt = False
                return

            self.mouse_mvt = False


class point_picker_2(object):
    """observers methods for a mayavi scene interactor. 
    Applied logic to make it record the picked position only
    if the mouse has not moved between lmb press and release.
    """

    mouse_mvt = False

    def __init__(self, figure, points_plot, callback, name=None):
        self.points_plot = points_plot
        self.callback = callback
        self.figure = figure
        self.figure.scene.picker.pointpicker.add_observer('EndPickEvent', self.on_p)

    def set_number_of_points_to_pick(self, n):
        """ clears the list of picked points and defines how many points
        to pick before calling the parent's _picking_done method """

        self.picked_points = []
        self.n = n

    def on_p(self, obj, evt):
        """ when p is pressed, gets picked position
        from scene pointpicker and calculates the closest data point,
        callback with the point number
        """

        ## new
        obj = tvtk.to_tvtk(obj)
        self.callback(obj.picked_position)

        return None


class point_picker_3(object):
    """observers methods for a mayavi scene interactor. 
    Applied logic to make it record the picked position only
    if the mouse has not moved between lmb press and release.
    """

    mouse_mvt = False
    pointScaleFactor = 0.2
    pointColour = (1.0, 0.0, 0.0)

    def __init__(self, figure, points_plot, callback, name=None):
        self.points_plot = points_plot
        self.callback = callback
        self.figure = figure
        self.figure.scene.picker.pointpicker.add_observer('EndPickEvent', self.on_p)
        self.do_callback = 0
        self.callbackArgs = None

    def set_callback_mode(self, args):
        self.callbackArgs = args

    def set_number_of_points_to_pick(self, n):
        """ clears the list of picked points and defines how many points
        to pick before calling the parent's _picking_done method """

        self.picked_points = np.zeros((n, 3))
        self.picked_plot = None
        self.n = n
        self.pickedn = 0
        return

    def on_p(self, obj, evt):
        """ when p is pressed, gets picked position
        from scene pointpicker and calculates the closest data point,
        adds to picked_points list. calls back with list when enough
        points are picked
        """

        if self.do_callback:

            p = obj.GetPickPosition()
            log.debug(p)
            self.picked_points[self.pickedn, :] = p
            self.pickedn += 1

            self.plot_picked_points()

            # ~ self.picked_point_plots.append( mlab.points3d( [p[0]], [p[1]], [p[2]] , mode='sphere', color=(1,0,0), figure=self.figure) )
            if self.pickedn == self.n:
                self.callback(tuple(self.picked_points), self.callbackArgs)
                self.picked_points = None
                self.picked_plot.remove()
                self.picked_plot = None
                self.pickedn = 0

        return None

    def plot_picked_points(self):
        view = self.figure.scene.mlab.view()
        roll = self.figure.scene.mlab.roll()
        self.figure.scene.disable_render = True
        try:
            self.picked_plot.remove()
        except AttributeError:
            pass

        p = self.picked_points[:self.pickedn, :]
        self.picked_plot = self.figure.scene.mlab.points3d(p[:, 0], p[:, 1], p[:, 2], mode='sphere',
                                                           scale_factor=self.pointScaleFactor, color=self.pointColour)

        self.figure.scene.mlab.view(view[0], view[1], view[2], view[3])
        self.figure.scene.mlab.roll(roll)
        self.figure.scene.disable_render = False

    def undo_pick(self):
        self.picked_points[self.pickedn, :] = 0.0
        self.pickedn -= 1
        self.plot_picked_points()
        return


class point_picker_3D(object):
    """observers methods for a mayavi scene interactor. 
    Applied logic to make it record the picked position only
    if the mouse has not moved between lmb press and release.
    """

    mouse_mvt = False

    def __init__(self, data, N, callback, mode='points', name=None):
        """ coords: point coordinates ( x, y, z ), N: total number of points to pick,
        callback: callback function
        """
        self.name = name
        self.N = N
        self.pickedPoints = []
        self.callback = callback
        self.figure = mlab.figure(bgcolor=(0, 0, 0), size=(600, 600))
        if mode == 'points':
            self.plot = mlab.points3d(data[0], data[1], data[2], mode='point')
        elif mode == 'contour':
            self.plot = mlab.contour3d(input)

        self.figure.scene.picker.pointpicker.add_observer('EndPickEvent', self.on_p)
        mlab.show()

    def on_p(self, obj, evt):
        """ when p is pressed, gets picked position and stores. If
        required number of points have been picked, calls back, and
        resets pickedPoints
        """

        obj = tvtk.to_tvtk(obj)
        self.pickedPoints.append(np.array(obj.pick_position))
        mlab.points3d([obj.pick_position[0]], [obj.pick_position[1]], [obj.pick_position[2]], mode='sphere',
                      color=(1, 1, 0), scale_factor=0.2)

        if len(self.pickedPoints) == self.N:
            if self.name != None:
                self.callback(list(self.pickedPoints), self.name)
            else:
                self.callback(list(self.pickedPoints))

            self.pickedPoints = []

        return

    def undoPick(self):
        del self.pickedPoints[-1]
        del self.figure.children[-1]


class element_adder(object):
    """ Class for interactively adding an element to a geometric field 
    by clicking on the geometric points the new element's element points 
    will connect to. Initialise with the geometric_field object.
    """

    def __init__(self, geometric_field):
        self.gf = geometric_field
        self.element = None

        # initialise plot
        self.figure = mlab.figure(1)
        self.points_plot = self.gf._plot_points()
        self.field_plot = self.gf._plot_field(10, glyph='sphere')

        # initialise point_picker
        self.picker = point_picker(self.figure.scene.picker.pointpicker,
                                   self.points_plot, self.figure, self)
        self.figure.scene.interactor.add_observer('LeftButtonPressEvent',
                                                  self.picker.on_button_press)
        self.figure.scene.interactor.add_observer('MouseMoveEvent',
                                                  self.picker.on_mouse_move)
        self.figure.scene.interactor.add_observer('LeftButtonReleaseEvent',
                                                  self.picker.on_button_release)

    def add_element(self, element):
        """ sets the point_picker for the user to select the required
        number of points """

        self.element = element
        # get number of points needed
        element_points = self.element.get_number_of_ensemble_points()
        log.debug('element requires', element_points, 'points.')

        self.picker.set_number_of_points_to_pick(element_points)
        mlab.show()

    def _picking_done(self):
        """ once picking is done, add the element with the picked point
        numbers. Re-draws the field plot with the new element """

        if self.gf.add_element(self.element, self.picker.picked_points):
            log.debug('element_added')
            self.field_plot.remove()
            self.field_plot = self.gf._plot_field(10, glyph='sphere')

        return 1


class element_adder_data(object):
    """ Class for interactively adding an element to a geometric field 
    by clicking on the geometric points the new element's element points 
    will connect to. Initialise with the geometric_field object.
    """
    evalD = 15

    def __init__(self, geometric_field, data=None, data_mesh=None, data_scalar=None):
        self.gf = geometric_field
        self.element = None
        self.plots = {}

        # initialise plot
        self.figure = mlab.figure(1, bgcolor=(0, 0, 0), size=(800, 800))

        if self.gf.field_parameters != None:
            self.field_plot = self.gf._draw_surface_curvature(self.evalD, 'mean', self.gf.name, figure=self.figure)
            # ~ K,H,k1,k2 = self.gf.evaluate_curvature_in_mesh( self.evalD )
            # ~ self.field_plot = self.gf._draw_surface( 10, figure=self.figure, scalar=H, name=self.gf.name)
            self.plots['field'] = self.figure.children[-1]

        self.point_plot = self.gf._plot_points(scale=0.5, figure=self.figure)
        if self.point_plot:
            self.plots['points'] = self.figure.children[-1]

        if data_mesh:
            self.data_plot = mlab.triangular_mesh(data_mesh[0][:, 0], data_mesh[0][:, 1], data_mesh[0][:, 2],
                                                  data_mesh[1], scalars=data_scalar, figure=self.figure)
            self.plots['data'] = self.figure.children[-1]
        else:
            if data != None:
                self.data = data
                self.crawler = spline_tools.pointCrawlerData(self.data.T)
                if data_scalar == None:
                    self.data_scalar = data_scalar
                    self.data_plot = mlab.points3d(data[0], data[1], data[2], mode='point', figure=self.figure,
                                                   scalars=data_scalar)
                else:
                    self.data_plot = mlab.points3d(data[0], data[1], data[2], data_scalar, mode='point',
                                                   figure=self.figure)

                self.plots['data'] = self.figure.children[-1]
            else:
                self.plots['data'] = self.figure.children[-1]
                self.data_plot = self.point_plot

        if self.gf.splines:
            for s in list(self.gf.splines.values()):
                self.gf._draw_spline(s, 200)

        # initialise point_picker
        self.picker = point_picker_3(self.figure, self.data_plot, self._picking_done)

        # ~ self.point_picker = point_picker_2( self.figure, self.point_plot, self._point_picker_callback )
        # ~ self.figure.scene.picker.pointpicker.add_observer('EndPickEvent', self.picker.on_p)

        # ~ self.figure.scene.interactor.add_observer('LeftButtonPressEvent',
        # ~ self.picker.on_button_press)
        # ~ self.figure.scene.interactor.add_observer('MouseMoveEvent',
        # ~ self.picker.on_mouse_move)
        # ~ self.figure.scene.interactor.add_observer('LeftButtonReleaseEvent',
        # ~ self.picker.on_button_release)
        # ~ mlab.show()

    def add_element(self, element):
        """ sets the point_picker for the user to select the required
        number of points """

        self.element = element
        # get number of points needed
        element_points = self.element.get_number_of_ensemble_points()
        log.debug('element requires', element_points, 'points.')
        self.picker.set_number_of_points_to_pick(element_points)
        self.picker.set_callback_mode('add_element')
        self.picker.do_callback = 1

    def add_points_by_curve(self, curve):
        """ prompts user to pick 2 points (p0,p1), to define the start 
        and end of a curve. uses the pointcrawler to trace a path 
        between p0 and p1. Fits curve to the path, and adds the fitted
        node positions as geometric points. Curve should be an EFF with
        one line element.
        """
        self.curve = curve
        log.debug('pick curve start and end positions')
        self.picker.set_number_of_points_to_pick(2)
        self.picker.set_callback_mode('add_curve')
        self.picker.do_callback = 1

    def undo_pick(self):
        self.picker.undo_pick()

    def _add_element_callback(self, picked_points):
        p = np.array(picked_points).transpose()[:, :, np.newaxis]
        if self.gf.add_element_with_parameters(self.element, p):
            log.debug('element added')
        return 1

    def _modify_point_callback(self, picked_points):
        p = np.array(picked_points[0])[:, np.newaxis]
        self.gf.modify_geometric_point(self.mod_point, p)
        return 1

    def _add_curve_callback(self, picked_points):

        # check for close points to existing geometric points
        proxTol = 0.5
        picked_points = np.array(picked_points)
        gp = np.array(self.gf.get_all_point_positions())
        if len(gp) > 0:
            for i, p in enumerate(picked_points):
                dist = ((gp - p) ** 2.0).sum(1)
                if dist.min() <= proxTol:
                    picked_points[i] = gp[np.argmin(dist)]
        try:
            pathPoints, pathCloud = self.crawler.trace(picked_points[0], picked_points[1], debug=0)
        except RuntimeError:
            log.debug('could not find path between points', picked_points)
            return
        else:
            # initialise curve p0
            nNodes = self.curve.get_number_of_ensemble_points()
            # assume lagrange
            if len(pathPoints) > nNodes:
                # put initial nodes along path
                p0 = np.zeros((3, nNodes))
                for i in range(nNodes):
                    p0[:, i] = pathPoints[round((i / float(nNodes - 1)) * (len(pathPoints) - 1))]

                p0 = p0[:, :, np.newaxis]
            else:
                # evenly distribute nodes between 2 ends, ignoring path         
                p0x = np.linspace(picked_points[0][0], picked_points[1][0], nNodes)
                p0y = np.linspace(picked_points[0][1], picked_points[1][1], nNodes)
                p0z = np.linspace(picked_points[0][2], picked_points[1][2], nNodes)
                p0 = np.array([p0x, p0y, p0z])[:, :, np.newaxis]

            # fit curve
            # ~ fittedNodes = spline_tools.fitCurveEPDP( self.curve, pathPoints, p0, debug=1 )
            fittedNodes = spline_tools.fitCurveDPEP(self.curve, pathCloud, p0, debug=1)
            # plot curve
            ep = []
            for nI in fittedNodes:
                ep.append(self.curve.evaluate_field_in_mesh(50, parameters=nI).ravel())
            mlab.plot3d(ep[0], ep[1], ep[2])
            # ~ print 'fitted curve nodes:', fittedPoints
            # add points to geometric field
            for i in range(nNodes):
                self.gf.add_geometric_point(fittedNodes[:, i, :])

            return

    def _picking_done(self, picked_points, mode):
        """ once picking is done, add the element with the picked point
        numbers. Re-draws the field plot with the new element """

        callbacks = {'add_element': self._add_element_callback,
                     'mod_point': self._modify_point_callback,
                     'add_curve': self._add_curve_callback,
                     }

        r = callbacks.get(mode)(picked_points)
        self._refresh()
        return r

    def modify_point(self, point):
        self.mod_point = point
        self.picker.set_number_of_points_to_pick(1)
        self.picker.set_callback_mode('mod_point')
        self.picker.do_callback = 1

    def moveNode(self, node, p):
        self.gf.modify_geometric_point(node, np.array(p)[:, np.newaxis])
        self._refresh()

    def _refresh(self):
        self.figure.scene.disable_render = True
        try:
            self.field_plot.remove()
        except AttributeError:
            pass
        else:
            self.plots['field'].remove()

        try:
            self.point_plot.remove()
        except AttributeError:
            pass
        else:
            self.plots['points'].remove()

        view = mlab.view()
        roll = mlab.roll()

        # ~ K,H,k1,k2 = self.gf.evaluate_curvature_in_mesh( self.evalD )
        # ~ self.field_plot = self.gf._draw_surface( self.evalD, figure=self.figure, scalar=H)
        self.field_plot = self.gf._draw_surface_curvature(self.evalD, 'mean', self.gf.name, figure=self.figure)
        if self.field_plot:
            self.plots['field'] = self.figure.children[-1]

        self.point_plot = self.gf._plot_points(scale=0.5, figure=self.figure)
        if self.point_plot:
            self.plots['points'] = self.figure.children[-1]

        self.picker.do_callback = 0
        mlab.view(view[0], view[1], view[2], view[3])
        mlab.roll(roll)
        self.figure.scene.disable_render = False

    # ~ def _point_picker_callback( self, pos ):
    # ~
    # ~ print obj.scalars
