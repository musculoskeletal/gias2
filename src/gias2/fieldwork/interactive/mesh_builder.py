"""
FILE: mesh_builder.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: functions and classes for interactively creating fieldwork meshes

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""
import logging

import numpy
from mayavi.core.ui.mayavi_scene import MayaviScene
from mayavi.tools.mlab_scene_model import MlabSceneModel
from traits.api import HasTraits, Button, Str, List, Instance
from traitsui.api import View, Item, HGroup, VGroup, EnumEditor, \
    TextEditor
from tvtk.pyface.scene_editor import SceneEditor

from gias2.fieldwork.field import ensemble_field_function as EFF
from gias2.fieldwork.field import geometric_field
from gias2.fieldwork.field.tools import spline_tools, fitting_tools
from gias2.fieldwork.field.topology import element_types
from gias2.fieldwork.interactive.geometric_field_interactor import point_picker_3

log = logging.getLogger(__name__)


class Viewer(HasTraits):
    initButton = Button()
    refreshButton = Button()
    drawCurvesButton = Button()
    drawDataMesh = False

    addCurveButton = Button()
    curveToAdd = 'line5l'
    curveTypes = List(['line2l', 'line3', 'line4l', 'line5l'])

    addElementButton = Button()
    elementToAdd = list(element_types.element_types.keys())[0]
    elementTypes = List(list(element_types.element_types.keys()))

    undoPickButton = Button()
    repositionNodeButton = Button()
    fitButton = Button()

    saveGFFilename = Str('')
    saveEnsFilename = Str('')
    saveMeshFilename = Str('')
    saveVersionNumber = Str('')
    saveGFButton = Button()
    loadGFButton = Button()

    view = View(HGroup(
        VGroup(
            # general operations
            Item('initButton', show_label=False, label='initialise', springy=True),
            Item('refreshButton', show_label=False, label='refresh', springy=True),
            Item('drawCurvesButton', show_label=False, label='draw curves', springy=True),

            # element adder
            VGroup(
                Item('elementToAdd', show_label=True, label='element types', editor=EnumEditor(name='elementTypes')),
                Item('addElementButton', show_label=False, label='add element', springy=True)
            ),

            # curve adder
            VGroup(
                Item('curveToAdd', show_label=True, label='curve types', editor=EnumEditor(name='curveTypes')),
                Item('addCurveButton', show_label=False, label='add curve', springy=True)
            ),

            Item('undoPickButton', show_label=False, label='undo pick', springy=True),
            Item('repositionNodeButton', show_label=False, label='reposition node', springy=True),
            Item('fitButton', show_label=False, label='fit mesh', springy=True),

            VGroup(
                Item('saveGFFilename', show_label=True, label='GF', editor=TextEditor()),
                Item('saveEnsFilename', show_label=True, label='ens', editor=TextEditor()),
                Item('saveMeshFilename', show_label=True, label='mesh', editor=TextEditor()),
                Item('saveVersionNumber', show_label=True, label='version', editor=TextEditor()),
                Item('saveGFButton', show_label=False, label='save GF', springy=True),
                Item('loadGFButton', show_label=False, label='load GF', springy=True),
            ),

            label='Mesh Builder'
        ),

        Item('scene', editor=SceneEditor(scene_class=MayaviScene), height=600, width=800, show_label=False)
    ),
        resizable=True,
    )

    scene = Instance(MlabSceneModel, ())
    nodeScaleFactor = 0.2
    dataScaleFactor = 0.1
    curveRadius = 2.0
    dataGlyph = '2dcross'

    def __init__(self):
        HasTraits.__init__(self)
        self.GF = None
        self.curves = []
        self.data = None
        self.dataScalar = None
        self.gfD = [20, 20]
        self.curveD = [20]
        self.gfScalar = 'mean'
        self.sceneObjectGF = None
        self.sceneObjectGFPoints = None
        self.sceneObjectCurves = []
        self.sceneObjectCurvesPoints = []
        self.picker = None
        self.sceneObjectData = None

    def setMeshBuilder(self, MB):
        self.MB = MB
        self.setGF(self.MB.GF)
        self.setSurfData(self.MB.data, scalar=self.MB.dataScalar)
        self.setCurves(self.MB.boundaryCurves)
        self.saveGFFilename = self.MB.GFFilename
        self.saveEnsFilename = self.MB.ensFilename
        self.saveMeshFilename = self.MB.meshFilename
        self.saveVersionNumber = self.MB.GFVersion

    def setGF(self, GF):
        self.GF = GF

    def setGFEvalD(self, d):
        self.gfD = d

    def setSurfData(self, data, scalar=None):
        self.data = data
        self.dataScalar = scalar

    def setCurves(self, curves):
        self.curves = curves

    def _initButton_fired(self):
        self._drawCurves()
        self._drawData()
        self.picker = point_picker_3(self, self.sceneObjectData, self._picking_done)
        self.picker.pointScaleFactor = self.nodeScaleFactor
        self.scene.background = (0.0, 0.0, 0.0)

    def _drawCurvesButton_fired(self):
        self._drawCurves()

    def _refreshButton_fired(self):

        self.scene.disable_render = True
        view = self.scene.mlab.view()
        roll = self.scene.mlab.roll()

        try:
            self.sceneObjectGF.remove()
            self.sceneObjectGFPoints.remove()
        except AttributeError:
            pass
        self._drawGeometricField()

        for i, c in enumerate(self.sceneObjectCurves):
            c.remove()
            self.sceneObjectCurvesPoints[i].remove()

        self.sceneObjectCurves = []
        self.sceneObjectCurvesPoints = []
        self._drawCurves()

        self.scene.mlab.view(view[0], view[1], view[2], view[3])
        self.scene.mlab.roll(roll)
        self.scene.disable_render = False

    def _saveGFButton_fired(self):
        self.MB.GFFilename = self.saveGFFilename
        self.MB.ensFilename = self.saveEnsFilename
        self.MB.meshFilename = self.saveMeshFilename
        self.MB.GFVersion = self.saveVersionNumber
        self.MB.saveGF()

    def _loadGFButton_fired(self):
        gf = self.saveGFFilename + '_' + self.saveVersionNumber + '.geof'
        ens = self.saveEnsFilename + '_' + self.saveVersionNumber + '.ens'
        mesh = self.saveMeshFilename + '_' + self.saveVersionNumber + '.mesh'
        self.MB.loadGF(gf, ens, mesh)
        self.setGF(self.MB.GF)
        self._refreshButton_fired()

        self.MB.GFFilename = self.saveGFFilename
        self.MB.ensFilename = self.saveEnsFilename
        self.MB.meshFilename = self.saveMeshFilename
        self.MB.GFVersion = self.saveVersionNumber

    def _undoPickButton_fired(self):
        self.picker.undo_pick()

    def _repositionNodeButton_fired(self):
        log.debug('pick node to move, then a point for its new position')
        self.picker.set_number_of_points_to_pick(2)
        self.picker.set_callback_mode('reposition')
        self.picker.do_callback = 1

    def _fitButton_fired(self):
        if self.MB.fitter != None:
            self.MB.fitMesh()
            self.setGF(self.MB.GF)
            self._refreshButton_fired()

    def _addElementButton_fired(self):
        self.element = element_types.create_element(self.elementToAdd)
        # get number of points needed
        element_points = self.element.get_number_of_ensemble_points()
        log.debug('element requires', element_points, 'points.')
        self.picker.set_number_of_points_to_pick(element_points)
        self.picker.set_callback_mode('addElement')
        self.picker.do_callback = 1

    def _addCurveButton_fired(self):

        # create curve GF
        cEns = EFF.ensemble_field_function('curve', 1, debug=0)
        cEns.set_basis({self.curveToAdd: self.MB.curveElemBasis[self.curveToAdd]})
        cEns.set_new_mesh('curve')
        cEns.create_elements(self.curveToAdd, 1)
        cEns.map_parameters()
        c = geometric_field.geometric_field('curve', 3, ensemble_field_function=cEns)

        self.curve = c
        log.debug('pick curve start and end positions')
        self.picker.set_number_of_points_to_pick(2)
        self.picker.set_callback_mode('addCurve')
        self.picker.do_callback = 1

    def _drawGeometricField(self):
        # evaluate mesh
        try:
            P = self.GF.evaluate_geometric_field(self.gfD)
        except:
            return
        else:
            # calc scalar
            S = None
            if self.gfScalar != 'none':
                K, H, k1, k2 = self.GF.evaluate_curvature_in_mesh(self.gfD)
                if self.gfScalar == 'mean':
                    S = H
                elif self.gfScalar == 'gaussian':
                    S = K

            T = self.GF.triangulator._triangulate(self.gfD)

            # draw
            # ~ pdb.set_trace()
            self.sceneObjectGF = self.scene.mlab.triangular_mesh(
                P[0], P[1], P[2], T,
                scalars=S,
                name=self.GF.name
            )
            points = self.GF.get_all_point_positions().T
            self.sceneObjectGFPoints = self.scene.mlab.points3d(
                points[0], points[1], points[2],
                numpy.arange(points.shape[1]),
                mode='sphere', scale_mode='none',
                scale_factor=self.nodeScaleFactor,
                name=self.GF.name
            )

    def _drawCurves(self):

        for c in self.curves:
            P = c.evaluate_geometric_field(self.curveD)

            # remove every density points to avoid overlap
            allPoints = P.T
            points = []
            for i, p in enumerate(allPoints):
                if numpy.mod(i, self.curveD[0]) != 0 or i == 0:
                    points.append(p)

            points = numpy.transpose(points)

            nodes = c.get_all_point_positions().T
            self.sceneObjectCurves.append(
                self.scene.mlab.plot3d(
                    points[0], points[1], points[2],
                    tube_radius=self.curveRadius
                )
            )
            self.sceneObjectCurvesPoints.append(
                self.scene.mlab.points3d(
                    nodes[0], nodes[1], nodes[2],
                    mode='sphere',
                    scale_factor=self.nodeScaleFactor,
                    name=c.name
                )
            )

            # display point labels
            nodeLabels = [str(i) for i in range(nodes.shape[1])]
            for i, l in enumerate(nodeLabels):
                self.scene.mlab.text3d(nodes[0, i], nodes[1, i], nodes[2, i], l, line_width=0.01)

    def _drawData(self):

        if self.dataScalar == None:
            self.sceneObjectData = self.scene.mlab.points3d(
                self.data[:, 0], self.data[:, 1], self.data[:, 2],
                mode=self.dataGlyph,
                scale_factor=self.dataScaleFactor
            )
        else:
            self.sceneObjectData = self.scene.mlab.points3d(
                self.data[:, 0], self.data[:, 1], self.data[:, 2],
                self.dataScalar,
                mode=self.dataGlyph,
                scale_factor=self.dataScaleFactor
            )

        if self.drawDataMesh and self.MB.simplemesh != None:

            if self.dataScalar == None:
                renderArgs = {'color': (0.84705882, 0.8, 0.49803922)}
                self.sceneObjectDataMesh = self.scene.mlab.triangular_mesh(
                    self.MB.simplemesh.v[:, 0],
                    self.MB.simplemesh.v[:, 1],
                    self.MB.simplemesh.v[:, 2],
                    self.MB.simplemesh.f,
                    **renderArgs
                )
            else:
                self.sceneObjectDataMesh = self.scene.mlab.triangular_mesh(
                    self.MB.simplemesh.v[:, 0],
                    self.MB.simplemesh.v[:, 1],
                    self.MB.simplemesh.v[:, 2],
                    self.MB.simplemesh.f,
                    scalars=self.dataScalar
                )

    def _picking_done(self, pickedPoints, mode):
        """ once picking is done, add the element with the picked point
        numbers. Re-draws the field plot with the new element """

        if mode == 'addElement':
            self.MB.addElementByPoints(self.element, pickedPoints)
            self.element = None
        elif mode == 'addCurve':
            self.MB.addCurve(self.curve, pickedPoints)
            self.curve = None
        elif mode == 'reposition':
            # get node number of 1st picked point
            p0 = pickedPoints[0]
            nodeCoords = self.MB.GF.get_all_point_positions()
            nodeN = numpy.sqrt(((nodeCoords - p0) ** 2.0).sum(1)).argmin()

            # change its coordinates to the 2nd picked point
            newParams = pickedPoints[1][:, numpy.newaxis]

            log.debug('moving node ' + str(nodeN) + ' to ' + str(newParams))
            self.MB.modifyPoint(nodeN, newParams)

        self.picker.do_callback = 0
        self._refreshButton_fired()
        return


class MeshBuilder(object):
    triMesh = None
    simplemesh = None
    boundaryCurves = []
    GF = None
    elementAdder = None
    fitter = None

    def __init__(self):
        pass

    def setSurfaceData(self, data, scalar=None):
        """ load a data cloud or trimesh on which the mesh will be built
        on
        """
        self.data = data
        self.dataScalar = scalar
        self.crawler = spline_tools.pointCrawlerData(self.data)

    def setSurfaceSimplemesh(self, mesh, scalar=None):
        self.simplemesh = mesh
        self.data = self.simplemesh.v
        self.dataScalar = scalar
        self.crawler = spline_tools.pointCrawlerData(self.data)

    def loadBoundaryCurve(self, gf, ens, mesh):
        self.boundaryCurves.append(geometric_field.load_geometric_field(gf, ens, mesh))

    def loadGF(self, gf, ens, mesh):
        self.GF = geometric_field.load_geometric_field(gf, ens, mesh)

    def initialiseGF(self, name, fieldDim, GFDim, elemBasis):
        self.GF = geometric_field.geometric_field(
            name, GFDim, field_dimensions=fieldDim, field_basis=elemBasis
        )

    def setCurveElemBasis(self, cEB):
        self.curveElemBasis = cEB

    def initVisualiser(self):
        pass

    def setGFFilenames(self, gf, ens, mesh):
        self.GFFilename = gf
        self.ensFilename = ens
        self.meshFilename = mesh

    def setGFVersion(self, v):
        self.GFVersion = v

    def saveGF(self):
        self.GF.save_geometric_field(
            self.GFFilename + '_' + self.GFVersion,
            self.ensFilename + '_' + self.GFVersion,
            self.meshFilename + '_' + self.GFVersion
        )

    def setFitter(self, fitter):
        self.fitter = fitter

    def fitMesh(self):
        self.fitter(self.GF)

    def addElementByPoints(self, e, points):
        p = numpy.array(points).transpose()[:, :, numpy.newaxis]
        if self.GF.add_element_with_parameters(e, p):
            log.debug('element added')
        return

    def modifyPoint(self, pointN, newParams):
        # ~ p = numpy.array( newParams[0] )[:,numpy.newaxis]
        p = newParams
        self.GF.modify_geometric_point(pointN, p)
        return

    def movePoint(self, pointN, move):
        p0 = self.GF.get_point_position(pointN)
        p1 = p0 + move
        self.GF.modify_geometric_point(pointN, p1[:, numpy.newaxis])
        return

    def addCurve(self, curve, picked_points):

        # check for close points to existing geometric points
        proxTol = 0.5
        picked_points = numpy.array(picked_points)
        gp = numpy.array(self.GF.get_all_point_positions())
        if len(gp) > 0:
            for i, p in enumerate(picked_points):
                dist = ((gp - p) ** 2.0).sum(1)
                if dist.min() <= proxTol:
                    picked_points[i] = gp[numpy.argmin(dist)]
        try:
            pathPoints, pathCloud = self.crawler.trace(picked_points[0], picked_points[1], debug=0)
        except RuntimeError:
            log.debug('could not find path between points', picked_points)
            return
        else:
            # initialise curve p0
            nNodes = curve.ensemble_field_function.get_number_of_ensemble_points()
            # assume lagrange
            if len(pathPoints) > nNodes:
                # put initial nodes along path
                p0 = numpy.zeros((3, nNodes))
                for i in range(nNodes):
                    p0[:, i] = pathPoints[round((i / float(nNodes - 1)) * (len(pathPoints) - 1))]

                p0 = p0[:, :, numpy.newaxis]
            else:
                # evenly distribute nodes between 2 ends, ignoring path         
                p0x = numpy.linspace(picked_points[0][0], picked_points[1][0], nNodes)
                p0y = numpy.linspace(picked_points[0][1], picked_points[1][1], nNodes)
                p0z = numpy.linspace(picked_points[0][2], picked_points[1][2], nNodes)
                p0 = numpy.array([p0x, p0y, p0z])[:, :, numpy.newaxis]

            curve.set_field_parameters(p0)

            # fit curve
            # ~ fittedNodes = splinet_tools.fitCurveEPDP( self.curve, pathPoints, p0, debug=1 )
            # ~ fittedNodes = splinet_tools.fitCurveDPEP( curve, pathCloud, p0, debug=1 )
            sobW = [1e-4, 1e-4]
            tangentW = 0.0
            fixedNodes = [0, nNodes - 1]
            curve, Opt, finalErr = fitting_tools.fitBoundaryCurve2WayFixNodes(
                curve, pathCloud, [50], sobW, tangentW,
                fixedNodes, itMax=20, fitVerbose=False
            )

            # ~ pdb.set_trace()
            # ~ curve = geometric_field.geometric_field( 'curve', 3, ensemble_field_function=curve )
            # ~ curveGF.set_field_parameters( fittedNodes )
            self.boundaryCurves.append(curve)
            # ~ print 'fitted curve nodes:', fittedPoints
            # add points to geometric field
            for p in curve.get_all_point_positions():
                self.GF.add_geometric_point(p[:, numpy.newaxis])

            return
