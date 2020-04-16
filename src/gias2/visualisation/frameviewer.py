"""
Module for viewing framesets in a matplotlib widget
"""

import json
import logging

import matplotlib.pyplot as plt
from matplotlib.image import AxesImage
from matplotlib.widgets import Slider, Button

log = logging.getLogger(__name__)


class FrameViewer(object):
    # parameters
    pick_marker = 'x'
    pick_out_path = 'picked_points.json'
    cmap = 'Greys'

    def __init__(self, frames=None, nodes=None, title=None, cmap=None):
        """
        Visualise target frames with nodes and sensor paths using a slider widget
        """
        self.frames = frames
        self.nodes = nodes
        self.title = title
        self.cmap = cmap

        self.n_frames = None
        if self.frames is not None:
            self.n_frames = len(self.frames)
        elif self.nodes is not None:
            self.n_frames = len(self.nodes)
        else:
            raise ValueError('Need at least 1 piece of input data')

        self.init_fig()

    def init_fig(self):

        # create figure
        self.fig, self.ax = plt.subplots(figsize=(10, 8))
        plt.subplots_adjust(left=0.1, bottom=0.25)
        if self.title is not None:
            self.ax.set_title(self.title)

        # create slider
        self.axslider = plt.axes([0.125, 0.125, 0.7, 0.03])
        self.slider = Slider(
            self.axslider, 'Frame', 1, self.n_frames, valinit=1
        )

        # draw frame
        if self.frames is not None:
            self.img = self.ax.imshow(
                self.frames[0], picker=True, cmap=self.cmap,
                vmin=self.frames.min(), vmax=self.frames.max()
            )
            self.ax.axis(
                (0, self.frames[0].shape[1], self.frames[0].shape[0], 0)
            )

        # draw nodes if we have them
        if self.nodes is not None:
            self.nodes_plot = self.ax.plot(
                self.nodes[0, :, 0], self.nodes[0, :, 1], '.g'
            )[0]

        # picker
        self.picked_points = {}
        for i in range(1, self.n_frames + 1):
            self.picked_points[i] = []

        self.picked_scatter = None

        # add undo button        left bot  width height
        self.ax_undo = plt.axes([0.85, 0.8, 0.08, 0.05])
        self.undo_button = Button(self.ax_undo, 'Undo')
        self.undo_button.on_clicked(self._on_undo)

        # add save button
        self.ax_save = plt.axes([0.85, 0.75, 0.08, 0.05])
        self.save_button = Button(self.ax_save, 'Save')
        self.save_button.on_clicked(self._on_save)

        # add close button
        self.ax_close = plt.axes([0.85, 0.7, 0.08, 0.05])
        self.close_button = Button(self.ax_close, 'Close')
        self.close_button.on_clicked(self._on_close)

        # start up
        self.slider.on_changed(self._update_slider)
        self.fig.canvas.mpl_connect('pick_event', self._on_img_pick)
        self.fig.show()

    def _frame_number(self):
        return int(round(self.slider.val))

    # what happens when a point is picked
    def _on_img_pick(self, event):
        artist = event.artist
        mevent = event.mouseevent
        if isinstance(artist, AxesImage):
            img_x = mevent.xdata
            img_y = mevent.ydata
            img_frame_num = self._frame_number()
            log.debug(
                'image {} picked at {}, {}'.format(
                    img_frame_num, img_x, img_y
                )
            )
            self.picked_points[img_frame_num].append((img_x, img_y))
            self._update_picked()

    # draw picked points
    def _update_picked(self):
        # destroy any previously drawn picked points
        if self.picked_scatter:
            try:
                self.picked_scatter.remove()
            except ValueError:
                pass

        # draw any picked points
        f = self._frame_number()
        if len(self.picked_points[f]):
            p_x = [p[0] for p in self.picked_points[f]]
            p_y = [p[1] for p in self.picked_points[f]]
            self.picked_scatter = self.ax.scatter(p_x, p_y, marker=self.pick_marker)

        self.fig.canvas.draw_idle()

    # what happens when slider is moved
    def _update_slider(self, val):
        f = self._frame_number()
        if self.frames is not None:
            self.img.set_data(self.frames[f - 1])

        # draw transformed nodes and sensor paths
        if self.nodes is not None:
            self.nodes_plot.set_xdata(self.nodes[f - 1, :, 0])
            self.nodes_plot.set_ydata(self.nodes[f - 1, :, 1])

        self._update_picked()
        # self.fig.canvas.draw_idle()
        self.fig.canvas.draw()

    def _on_undo(self, event):
        f = self._frame_number()
        if len(self.picked_points[f]):
            last_p = self.picked_points[f].pop()
            log.debug('removed point {},{}'.format(last_p[0], last_p[1]))
            self._update_picked()

    def _on_save(self, event):
        with open(self.pick_out_path, 'w') as f:
            json.dump(self.picked_points, f, indent=4, sort_keys=True)
        log.debug('picked points saved to {}'.format(self.pick_out_path))

    def _on_close(self, event):
        log.debug('Use the window close button')
        # plt.close(self.fig)
