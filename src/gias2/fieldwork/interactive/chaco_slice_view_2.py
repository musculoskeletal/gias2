"""
FILE: chaco_slice_view2.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: A chaco widget for viewing image slices.

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

from enthought.traits.api import Bool, Any, HasTraits, Instance, Int, Array, Range, on_trait_change, Property, Float
from enthought.traits.ui.api import View, Group, Item
from enthought.enable.component_editor import ComponentEditor
from enthought.chaco.api import Plot, ArrayPlotData, bone, BaseTool
from enthought.chaco.tools.api import DataPrinter, ZoomTool, PanTool, ImageInspectorTool, ImageInspectorOverlay

from scipy import arange, zeros

class ClickPositionTool(BaseTool):
    """ A tool to get mouse click position
    """

    # This callback will be called with the index into self.component's
    # index and value:
    #     callback(tool, x_index, y_index)
    # where *tool* is a reference to this tool instance.  The callback
    # can then use tool.token.
    callback = Any()

    # This callback (if it exists) will be called with the integer number
    # of mousewheel clicks
    wheel_cb = Any()

    # This token can be used by the callback to decide how to process
    # the event.
    token = Any()

    # Whether or not to update the slice info; we enter select mode when
    # the left mouse button is pressed and exit it when the mouse button
    # is released
    # FIXME: This is not used right now.
    select_mode = Bool(False)
            
    def normal_left_down(self, event):
        #~ self._update_slices(event)
        self._update_position(event)

    #~ def normal_right_down(self, event):
        #~ self._update_slices(event)
#~ 
    #~ def normal_mouse_move(self, event):
        #~ if event.left_down or event.right_down:
            #~ self._update_slices(event)

    def _update_position(self, event):
        plot = self.component
        xy = plot.map_data((event.x, event.y))
        if xy != None:
            self.callback(xy)


class slice_viewer(HasTraits):

    #~ start_img = zeros((100,100))
    plotdata = ArrayPlotData( imagedata = None )
    plot = Plot( plotdata )
    #~ imgplot = plot.img_plot( 'imagedata', xbounds=arange(100), ybounds=arange(100), colormap=bone )[0]

    image = Array
    slicehigh = Property( depends_on=['image'] )
    slicelow = Property
    slice_range = Range( low='slicelow', high ='slicehigh', mode='slider' )
    slice_view = View( Item( 'plot', editor=ComponentEditor(), show_label=False),\
                       Item( name='slice_range'),\
                       width=800, height=600, resizable=True, title='Slice Viewer' )
                       
    
    def __init__(self, image):
        super( slice_viewer, self ).__init__()
        
        self.collecting_points = False      # toggle for whether click positions are recorded
        self.image = image
        self.plotdata.set_data( 'imagedata', self.image[ self.slice_range, :, :] )
        imgplot = self.plot.img_plot( 'imagedata', xbounds=arange(self.image.shape[2]),\
                                      ybounds=arange(self.image.shape[1]), colormap=bone )[0]
        
        self.add_tools( self.plot )
        # image inspector
        #~ self.imgtool = ImageInspectorTool(imgplot)
        #~ imgplot.tools.append(self.imgtool)
        #~ overlay = ImageInspectorOverlay(component=imgplot, image_inspector=self.imgtool,
                                    #~ bgcolor="white", border_visible=True)
        #~ imgplot.overlays.append(overlay)
    
    def add_tools( self, plot ):
        """ add ClickPositionTool, Zoom and Pan tools to the input plot
        """
        
        zoom = ZoomTool( plot, tool_mode='box' )
        plot.overlays.append(zoom)
        pan = PanTool(plot, drag_button="right", constrain_key="shift")
        plot.tools.append(pan)
        datapicker = ClickPositionTool(component=self.plot, callback=self._index_callback)
        self.plot.tools.append( datapicker )
        return
        
    def collect_n_points( self, npoints, callback ):
        """ prepares object for recording the positions of clicks.
        Records npoints points, and calls the callback when done
        """
        self.collecting_points = True
        self.points = []
        self.npoints = npoints
        self.collected_points_callback = callback
        #~ self.click_callback = callback
        return
    
    def _index_callback( self, xy ):
        """ callback functions for the ClickPositionTool. If collecting_points,
        positions of clicks are appended to self.points and when number
        of collected points == self.npoints, calls the collected_points_callback
        function
        """
        print('position:', [self.slice_range, xy[1], xy[0]])
        #~ self.click_callback( [self.slice_range, xy[1], xy[0]] )
        #~ return
        
        if self.collecting_points:
            self.points.append( [self.slice_range, xy[1], xy[0]] )
            if len(self.points) == self.npoints:
                points = self.points
                self.points = None
                self.collecting_points = False
                self.collected_points_callback( points )
                return
    
    def _get_slicehigh( self ):
        return self.image.shape[0] - 1
        
    def _get_slicelow( self ):
        return 0

    @on_trait_change('slice_range') 
    def calculate_slice_image(self):
        #~ print self.slice_range
        self.plotdata.set_data( 'imagedata', self.image[self.slice_range, :, :] )

