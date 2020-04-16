"""
FILE: fieldvi.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: Mayavi-based widget for viewing 3D data.
 
===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""
import logging

import scipy
from mayavi.core.ui.mayavi_scene import MayaviScene
from mayavi.tools.mlab_scene_model import MlabSceneModel
from traits.api import HasTraits, Range, Bool, Button, Str, Enum, List, Instance
from traitsui.api import View, Item, HGroup, HSplit, VGroup, Tabbed, EnumEditor, \
    TextEditor
from tvtk.api import tvtk
from tvtk.pyface.scene_editor import SceneEditor

try:
    from mayavi import mlab
except ImportError:
    raise ImportError('Mayavi not installed')

log = logging.getLogger(__name__)


class Fieldvi(HasTraits):
    renderAll = Button()
    renderImagePlane = Button()
    renderGeometricFields = Button()
    renderTri = Button()
    renderData = Button()

    imageList0 = 'None'
    imageList = List([])
    imagePlane = Enum(('x_axes', 'y_axes', 'z_axes'))
    imageVisible = Bool(False)

    dataList0 = 'None'
    dataList = List([])
    dataVisible = Bool(True)
    dataUpdate = Button()

    triList0 = 'None'
    triList = List([])
    triScalarList0 = 'None'
    triScalarList = List([])
    triVisible = Bool(True)
    triUpdate = Button()
    triCutPlane = Button()
    refreshTriScalars = Button()

    GFList0 = 'None'
    GFList = List([])
    GFScalarList0 = 'None'
    GFScalarList = List([])
    GFVisible = Bool(True)
    GFUpdate = Button()
    GFCutPlane = Button()
    refreshGFScalars = Button()

    PCList0 = 'None'
    PCList = List([])
    PCGeomList0 = 'None'
    PCGeomList = List([])
    modeIndex = Range(0, 13)
    mode1 = Range(-2., 2., 0.0)
    modeQuiver = Button()

    SFModeIndex = Range(0, 13)
    SFModeWeight = Range(-2., 2., 0.0)
    SFModeQuiver = Button()

    GFSFModeIndex = Range(0, 13)
    GFSFModeWeight = Range(-2., 2., 0.0)
    GFSFModeQuiver = Button()

    # LL Atlas
    LLPC0 = Range(-2., 2., 0.0)
    LLPC1 = Range(-2., 2., 0.0)
    LLPC2 = Range(-2., 2., 0.0)
    LLPC3 = Range(-2., 2., 0.0)
    LLPC4 = Range(-2., 2., 0.0)
    LLHJC1 = Range(-150.0, +150.0, 0.0)
    LLHJC2 = Range(-90.0, +90.0, 0.0)
    LLHJC3 = Range(-90.0, +90.0, 0.0)
    LLKnee1 = Range(-150.0, 0.0, 0.0)

    # 2 rigid-body
    RBAX0 = Range(-180., +180., 0.0)
    RBAX1 = Range(-180., +180., 0.0)
    RBAX2 = Range(-180., +180., 0.0)

    # screenshot
    saveImageFilename = Str('screenshot.jpeg')
    saveImageWidth = Str('2000')
    saveImageLength = Str('1500')
    saveImage = Button()

    scene = Instance(MlabSceneModel, ())

    mergeGFVertices = False
    displayGFNodes = True
    GFScalar = 'none'  # 'none', 'mean', or 'gaussian', or 'custom'
    triScalar = 'none'

    colours = {'bone': (0.84705882, 0.8, 0.49803922)}
    defaultColor = colours['bone']  # bone
    # ~ defaultColor = (0.5,0.5,0.5) # bone

    view = View(HSplit(
        Tabbed(
            VGroup(
                Item('renderAll', show_label=False, label='render all objects', springy=True),
                # ~ Item('renderImagePlane', show_label=False, label='render image plane', springy=True),

                Item('renderGeometricFields', show_label=False, label='render geometric fields', springy=True),
                Item('renderData', show_label=False, label='render data clouds', springy=True),

                VGroup(
                    Item('imageList0', show_label=True, label='images', editor=EnumEditor(name='imageList')),
                    HGroup(
                        Item('imageVisible', show_label=True, label='visible', springy=True),
                        Item('imagePlane', show_label=True, label='image plane', springy=True),
                    )
                ),

                VGroup(
                    Item('dataList0', show_label=True, label='Dataclouds', editor=EnumEditor(name='dataList')),
                    HGroup(
                        Item('dataVisible', show_label=True, label='visible', springy=True),
                        Item('dataUpdate', show_label=False, label='update', springy=True)
                    )
                ),

                VGroup(
                    Item('GFList0', show_label=True, label='GFs', editor=EnumEditor(name='GFList')),
                    HGroup(
                        Item('GFScalarList0', show_label=True, label='GFScalars',
                             editor=EnumEditor(name='GFScalarList')),
                        Item('refreshGFScalars', show_label=False, label='refresh', springy=True)
                    ),
                    HGroup(
                        Item('GFVisible', show_label=True, label='visible', springy=True),
                        Item('GFCutPlane', show_label=False, label='cut plane', springy=False),
                        Item('GFUpdate', show_label=False, label='update', springy=True)
                    )
                ),

                VGroup(
                    Item('triList0', show_label=True, label='tri surfaces', editor=EnumEditor(name='triList')),
                    HGroup(
                        Item('triScalarList0', show_label=True, label='triScalars',
                             editor=EnumEditor(name='triScalarList')),
                        Item('refreshTriScalars', show_label=False, label='refresh', springy=True)
                    ),
                    HGroup(
                        Item('triVisible', show_label=True, label='visible', springy=True),
                        Item('triCutPlane', show_label=False, label='cut plane', springy=False),
                        Item('triUpdate', show_label=False, label='update', springy=True)
                    )
                ),

                VGroup(
                    Item('saveImageFilename', show_label=True, label='filename', editor=TextEditor()),
                    HGroup(
                        Item('saveImageWidth', show_label=True, label='W', editor=TextEditor()),
                        Item('saveImageLength', show_label=True, label='L', editor=TextEditor()),
                    ),
                    Item('saveImage', show_label=False, label='saveImage', springy=True),
                ),

                label='Render'
            ),
            VGroup(
                Item('PCList0', show_label=True, label='PC Models', editor=EnumEditor(name='PCList')),
                Item('PCGeomList0', show_label=True, label='PC Geometry', editor=EnumEditor(name='PCGeomList')),
                Item('modeIndex', style='custom', springy=True),
                Item('mode1', show_label=False, label='Mode 1', springy=True),
                Item('modeQuiver', show_label=False, label='mode vectors', springy=False),
                label='statistical shape model'
            ),

            VGroup(
                Item('SFModeIndex', style='custom', springy=True),
                Item('SFModeWeight', show_label=False, label='Mode Weight', springy=True),
                Item('SFModeQuiver', show_label=False, label='Mode Vectors', springy=False),
                label='statistical surface field model'
            ),

            VGroup(
                Item('GFSFModeIndex', style='custom', springy=True),
                Item('GFSFModeWeight', show_label=False, label='Mode Weight', springy=True),
                Item('GFSFModeQuiver', show_label=False, label='Mode Vectors', springy=False),
                label='statistical geometry/surface field model'
            ),

            VGroup(
                Item('LLPC0', style='custom', springy=True, label='PC 1'),
                Item('LLPC1', style='custom', springy=True, label='PC 2'),
                Item('LLPC2', style='custom', springy=True, label='PC 3'),
                Item('LLPC3', style='custom', springy=True, label='PC 4'),
                Item('LLPC4', style='custom', springy=True, label='PC 5'),
                Item('LLHJC1', style='custom', springy=True, label='Hip Flexion'),
                Item('LLHJC2', style='custom', springy=True, label='Hip Rotation'),
                Item('LLHJC3', style='custom', springy=True, label='Hip Adduction'),
                Item('LLKnee1', style='custom', springy=True, label='Knee Flexion'),
                label='Lower limb atlas'
            ),

            VGroup(
                Item('RBAX0', style='custom', springy=True, label='Axis 0'),
                Item('RBAX1', style='custom', springy=True, label='Axis 1'),
                Item('RBAX2', style='custom', springy=True, label='Axis 2'),
                label='Rigid Bodies'
            ),
            springy=False),

        Item('scene', editor=SceneEditor(scene_class=MayaviScene), height=600, width=800, show_label=False),
    ),
        resizable=True,
    )

    def __init__(self):
        HasTraits.__init__(self)

        self.I = None

        self.images = {}
        self.imageCounter = 0
        self.sceneObjectImages = {}
        self.imageRenderArgs = {}

        self.data = {}
        self.dataScalar = {}
        self.dataCounter = 0
        self.sceneObjectData = {}
        self.dataRenderArgs = {}

        self.triCounter = 0
        self.triSurface = {}
        self.triScalarData = {}
        self.sceneObjectTri = {}
        self.triRenderArgs = {}

        self.geometricFields = {}
        self.GFCounter = 0
        self.GFD = [10, 10]
        self.GFDSpecific = {}
        self.GFEvaluators = {}
        self.GFScalarData = {}
        self.scalarData = {}
        self.scalarType = {}
        self.sceneObjectGF = {}
        self.sceneObjectGFPoints = {}
        self.sceneObjectGFCutPlanes = {}
        self.GFRenderArgs = {}
        self.bCurves = {}

        self.sceneObjectImageVolume = None
        self.modeVectorFlip = 1.0

        self.pc = None
        self.pcGFName = None
        self.SFPC = None
        self.SFPCFieldName = None
        self.GFSFPC = None
        self.GFSFPCFieldName = None

        self.PCs = {}

        self.view = None
        self.roll = None

        self.onCloseCallback = None

        # image plane widget attributes
        self._ipw_picked_obj = None
        self._ipw_picked_points = []

    def start(self):
        self.edit_traits()

    def addOnCloseCallback(self, callback):
        self.onCloseCallback = callback

    def closed(self, info, is_ok):
        if self.onCloseCallback is not None:
            self.onCloseCallback()

    def addTri(self, name, tri, renderArgs=None):

        if name not in list(self.triSurface.keys()):
            self.triCounter += 1
            self.triList.append(name)
            self.PCGeomList.append(name + '::tri')

        self.triSurface[name] = tri
        self.triScalarData[name] = {}
        self.addTriScalarData(name, 'none', None)
        self.triList0 = name

        if renderArgs is None:
            self.triRenderArgs[name] = {}
        else:
            self.triRenderArgs[name] = renderArgs

    def addTriScalarData(self, triName, scalarName, scalarData):
        self.triScalarData[triName][scalarName] = scalarData
        self.triScalarList0 = scalarName

    def addGeometricField(self, name, G, evaluator, GD=None, renderArgs=None):

        if name not in list(self.geometricFields.keys()):
            self.GFCounter += 1
            self.GFList.append(name)
            self.PCGeomList.append(name + '::gf')

        self.geometricFields[name] = G
        self.GFEvaluators[name] = evaluator
        self.GFScalarData[name] = []
        if GD is None:
            self.GFDSpecific[name] = self.GFD
        else:
            self.GFDSpecific[name] = GD
        self.addGeometricFieldScalarData(name, 'none', None)
        self.GFList0 = name

        if renderArgs is None:
            self.GFRenderArgs[name] = {}
        else:
            self.GFRenderArgs[name] = renderArgs

    def addGeometricFieldScalarData(self, GFNames, scalarName, scalarData, dataType='data'):
        """
        dataType is either 'data' or 'field'
        """

        if isinstance(GFNames, str):
            GFNames = [GFNames]

        for GFName in GFNames:
            self.GFScalarData[GFName].append(scalarName)

        self.GFScalarList0 = scalarName
        self.scalarType[scalarName] = dataType
        self.scalarData[scalarName] = scalarData
        self.GFScalar = 'custom'

    def addImageVolume(self, I, name=None, renderArgs=None):
        if name is None:
            name = str(self.imageCounter)

        if renderArgs is None:
            self.imageRenderArgs[name] = {}
        else:
            self.imageRenderArgs[name] = renderArgs
        if name not in list(self.images.keys()):
            self.imageCounter += 1
            self.imageList.append(name)

        self.imageList0 = name
        self.images[name] = I

    def addData(self, name, d, scalar=None, renderArgs=None):

        if name not in list(self.data.keys()):
            self.dataCounter += 1
            self.dataList.append(name)
            self.PCGeomList.append(name + '::data')

        self.data[name] = d
        if renderArgs is None:
            self.dataRenderArgs[name] = {}
        else:
            self.dataRenderArgs[name] = renderArgs
        if scalar is not None:
            self.dataScalar[name] = scalar
        else:
            self.dataScalar[name] = None

        self.dataList0 = name

    def removeData(self, name):
        try:
            del self.data[name]
            del self.dataRenderArgs[name]
            del self.dataScalar[name]
        except KeyError:
            log.debug('data', name, 'does not exist')
        else:
            self.dataList.remove(name)
            self.PCGeomList.remove[name + '::data']
            self.dataCounter -= 1
            try:
                self.sceneObjectData[name].remove()
                del self.sceneObjectData[name]
            except KeyError:
                pass

            self.dataList0 = self.dataList[-1]

    def addPC(self, pcName, pc):
        self.PCs[pcName] = pc
        self.PCList.append(pcName)
        self.PCList0 = pcName

    # def addPCToGF(self, pc, GFName):
    #   self.PCsGF[GFName] = pc

    # def addPCToData(self, pc, dataName):
    #   self.PCsData[dataName] = pc

    # def addPCToTriMesh(self, pc, triName):
    #   self.PCsTri[triName] = pc

    def addSFPCA(self, pc, SFName):
        self.SFPC = pc
        self.SFPCFieldName = SFName

    def addGFSFPCA(self, pc, GFName, SFName, GFRows, SFRows, comb=True):
        self.GFSFPC = pc
        self.GFSFGFName = GFName
        self.GFSFSFName = SFName
        self.GFSFGFRows = GFRows
        self.GFSFSFRows = SFRows
        self.GFSFPCComb = comb

    def _saveImage_fired(self):
        self.scene.mlab.savefig(str(self.saveImageFilename), size=(int(self.saveImageWidth), int(self.saveImageLength)))

    def _drawAxes(self):
        o = scipy.array([[0, 0, 0], ])
        o2 = scipy.array([0, 0, 0])
        self.addData('origin', o)
        a = scipy.array([[15, 0, 0], [0, 10, 0], [0, 0, 5]])
        self.addData('axes', a)

        self.scene.disable_render = True
        self.sceneObjectData['origin'] = self.scene.mlab.points3d(o[:, 0], o[:, 1], o[:, 2], mode='sphere',
                                                                  scale_factor=4.0)
        self.sceneObjectData['axes'] = self.scene.mlab.quiver3d(o2, o2, o2,
                                                                a.T[0], a.T[1], a.T[2],
                                                                mode='arrow',
                                                                line_width=2.0,
                                                                scale_factor=2.0)
        self.scene.disable_render = False

    def _renderAll_fired(self):
        """ redraws all objects
        """

        try:
            self.scene.mlab.clf()
        except:
            pass

        self.sceneObjectTri = {}
        self.sceneObjectGF = {}
        self.sceneObjectData = {}
        self.sceneObjectImageVolume = None

        for i in list(self.triSurface.keys()):
            self._drawTriSurface(i)

        for i in list(self.geometricFields.keys()):
            self._drawGeometricField(i)

        for i in list(self.data.keys()):
            self._drawData(i)

        for i in list(self.images.keys()):
            self._drawImage(i)

    def _drawImage(self, name):
        self.scene.disable_render = True

        I = self.images[name]
        renderArgs = self.imageRenderArgs[name]
        if 'colormap' not in renderArgs:
            renderArgs['colormap'] = 'black-white'

        ISrc = mlab.pipeline.scalar_field(I, colormap=renderArgs['colormap'])
        self.sceneObjectImages[name] = self.scene.mlab.pipeline.image_plane_widget(ISrc,
                                                                                   plane_orientation=self.imagePlane,
                                                                                   slice_index=0,
                                                                                   # figure=self.scene,
                                                                                   **renderArgs
                                                                                   )
        mlab.outline()
        self.scene.disable_render = False

        # attach pickers
        self.sceneObjectImages[name].ipw.add_observer(
            'EndWindowLevelEvent',
            self._ipw_pick_callback
        )

    def imagePlaneSliceIndex(self, sliceIndex=None, plane=0):
        if sliceIndex is None:
            return self.sceneObjectImages[self.imageList0].widgets[plane].slice_index
        else:
            self.sceneObjectImages[self.imageList0].widgets[plane].slice_index = sliceIndex

    def _imagePlane_changed(self):
        self.sceneObjectImages[self.imageList0].widgets[0].set(plane_orientation=self.imagePlane)

    def _imageVisible_changed(self):
        try:
            self.sceneObjectImages[self.imageList0].visible = self.imageVisible
        except KeyError:
            if self.imageVisible:
                self._drawImage(self.imageList0)

    def _dataVisible_changed(self):
        try:
            self.sceneObjectData[self.dataList0].visible = self.dataVisible
        except KeyError:
            if self.dataVisible:
                self._drawData(self.dataList0)

    def _dataUpdate_fired(self):
        self.updateData(self.dataList0)

    def updateData(self, name, coords=None):
        if coords is None:
            d = self.data[name]
        else:
            d = coords
        s = self.dataScalar.get(name)
        renderArgs = self.dataRenderArgs[name]

        try:
            if s is not None:
                self.sceneObjectData[name].actor.mapper.scalar_visibility = True
                self.sceneObjectData[name].mlab_source.reset(x=d[:, 0], y=d[:, 1], z=d[:, 2], s=s)
            else:
                if 'color' not in renderArgs:
                    color = self.defaultColor
                else:
                    color = renderArgs['color']

                self.sceneObjectData[name].actor.mapper.scalar_visibility = False
                self.sceneObjectData[name].actor.property.specular_color = color
                self.sceneObjectData[name].actor.property.diffuse_color = color
                self.sceneObjectData[name].actor.property.ambient_color = color
                self.sceneObjectData[name].actor.property.color = color
                self.sceneObjectData[name].mlab_source.reset(x=d[:, 0], y=d[:, 1], z=d[:, 2])
        except KeyError:
            self._drawData(name)

    def _renderData_fired(self):
        """ redraw all data clouds
        """
        for i in list(self.sceneObjectData.values()):
            i.remove()

        self.sceneObjectData = {}
        for i in list(self.data.keys()):
            self._drawData(i)

    def _drawData(self, name):
        self.scene.disable_render = True
        d = self.data[name]
        s = self.dataScalar.get(name)
        renderArgs = self.dataRenderArgs[name]
        if s is not None:
            # ~ self.sceneObjectData[name] = self.scene.mlab.points3d( d[:,0], d[:,1], d[:,2], s,
            # ~ mode='point', scale_factor=0.5,
            # ~ scale_mode='none',
            # ~ name=name )
            self.sceneObjectData[name] = self.scene.mlab.points3d(d[:, 0], d[:, 1], d[:, 2], s, **renderArgs)
        else:
            # ~ self.sceneObjectData[name] = self.scene.mlab.points3d( d[:,0], d[:,1], d[:,2],
            # ~ mode='point', scale_factor=0.5,
            # ~ color=(0.0,1.0,0.0),
            # ~ name=name )
            self.sceneObjectData[name] = self.scene.mlab.points3d(d[:, 0], d[:, 1], d[:, 2], **renderArgs)

        self.scene.disable_render = False

    def _refreshGFScalars_fired(self):
        log.debug('refresh scalars!')
        self.GFScalarList = self.GFScalarData[self.GFList0]

    def _refreshTriScalars_fired(self):
        log.debug('refresh scalars!')
        self.triScalarList = list(self.triScalarData[self.triList0].keys())

    def _triVisible_changed(self):
        try:
            self.sceneObjectTri[self.triList0].visible = self.triVisible
        except KeyError:
            if self.triVisible:
                self._drawTriSurface(self.triList0)

    def _GFVisible_changed(self):
        try:
            scnObj = self.sceneObjectGF[self.GFList0]
            if isinstance(scnObj, (list, tuple)):
                for _scnObj in scnObj:
                    _scnObj.visible = self.GFVisible
            else:
                self.sceneObjectGF[self.GFList0].visible = self.GFVisible
            if self.displayGFNodes:
                self.sceneObjectGFPoints[self.GFList0].visible = self.GFVisible
        except KeyError:
            if self.GFVisible:
                self._drawGeometricField(self.GFList0)

    def _triUpdate_fired(self):
        self.updateTriSurface(self.triList0)

    def _GFUpdate_fired(self):
        self.updateGeometricField(self.GFList0)

    def _renderTriSurfaces_fired(self):
        """ redraw all geometric fields
        """
        for i in list(self.sceneObjectTri.values()):
            i.remove()

        self.sceneObjectTri = {}
        for i in list(self.triSurface.keys()):
            self._drawTriSurface(i)

    def _renderGeometricFields_fired(self):
        """ redraw all geometric fields
        """
        for i in list(self.sceneObjectGF.values()):
            i.remove()

        self.sceneObjectGF = {}
        for i in list(self.geometricFields.keys()):
            self._drawGeometricField(i)

    def _triCutPlane_fired(self):
        """ create a vector cut plane for the current selected geometric
        field. If one already exists, remove it
        """
        if self.sceneObjectTriCutPlanes.get(self.triList0):
            self.sceneObjectTriCutPlanes[self.triList0].remove()
            del self.sceneObjectTriCutPlanes[self.triList0]
        else:
            self.sceneObjectTriCutPlanes[self.triList0] = mlab.pipeline.vector_cut_plane(
                self.sceneObjectTri[self.triList0],
                mask_points=1,
                mode='sphere',
                scale_factor=3)

    def _GFCutPlane_fired(self):
        """ create a vector cut plane for the current selected geometric
        field. If one already exists, remove it
        """
        if self.sceneObjectGFCutPlanes.get(self.GFList0):
            self.sceneObjectGFCutPlanes[self.GFList0].remove()
            del self.sceneObjectGFCutPlanes[self.GFList0]
        else:
            self.sceneObjectGFCutPlanes[self.GFList0] = mlab.pipeline.vector_cut_plane(self.sceneObjectGF[self.GFList0],
                                                                                       mask_points=1,
                                                                                       mode='sphere',
                                                                                       scale_factor=1.0,
                                                                                       )
            self.sceneObjectGFCutPlanes[self.GFList0].handle_size = 0.001
            try:
                self.sceneObjectGFCutPlanes[self.GFList0].origin = [0, 0,
                                                                    self.sceneObjectImages['0'].widgets[0].slice_index]
            except KeyError:
                pass

    def updateTriSurface(self, name, coords=None):
        if coords is None:
            V = self.triSurface[name].v
        else:
            V = coords
        renderArgs = self.triRenderArgs.get(name)
        scalar = self.triScalarData.get(name).get(self.triScalarList0)
        # log.debug('scalar:', scalar)
        try:
            if scalar is None:
                if 'color' not in renderArgs:
                    color = self.defaultColor
                else:
                    color = renderArgs['color']

                self.sceneObjectTri[name].actor.mapper.scalar_visibility = False
                self.sceneObjectTri[name].actor.property.specular_color = color
                self.sceneObjectTri[name].actor.property.diffuse_color = color
                self.sceneObjectTri[name].actor.property.ambient_color = color
                self.sceneObjectTri[name].actor.property.color = color
                self.sceneObjectTri[name].mlab_source.set(x=V[:, 0], y=V[:, 1], z=V[:, 2], **renderArgs)
            else:
                if 'color' in renderArgs:
                    del renderArgs['color']

                self.sceneObjectTri[name].actor.mapper.scalar_visibility = True
                self.sceneObjectTri[name].mlab_source.set(x=V[:, 0], y=V[:, 1], z=V[:, 2],
                                                          scalars=self.triScalarData[name][self.triScalarList0],
                                                          **renderArgs)

        except KeyError:
            self._drawTriSurface(name)

    def updateGeometricField(self, name, params=None):

        if self.sceneObjectGF.get(name) is None:
            self._drawGeometricField(name)
        else:
            if params is None:
                params = self.geometricFields[name].get_field_parameters().ravel()

            V = self.GFEvaluators[name](params)
            p = params.reshape((3, -1))

            if self.mergeGFVertices:
                V = V[:, self.uniqueVertexIndices]

            scalar = self._getGFScalarData(self.GFScalarList0, name)
            renderArgs = self.GFRenderArgs[name]

            if scalar is None:
                if 'color' not in renderArgs:
                    color = self.defaultColor
                else:
                    color = renderArgs['color']

                self.sceneObjectGF[name].actor.mapper.scalar_visibility = False
                self.sceneObjectGF[name].actor.property.specular_color = color
                self.sceneObjectGF[name].actor.property.diffuse_color = color
                self.sceneObjectGF[name].actor.property.ambient_color = color
                self.sceneObjectGF[name].actor.property.color = color
                self.sceneObjectGF[name].mlab_source.set(x=V[0], y=V[1], z=V[2])
            else:
                if self.mergeGFVertices:
                    if scalar.shape[0] != V.shape[1]:
                        scalar = scalar[self.uniqueVertexIndices]
                self.sceneObjectGF[name].mlab_source.set(x=V[0], y=V[1], z=V[2], scalars=scalar)
                self.sceneObjectGF[name].actor.mapper.scalar_visibility = True

            if self.displayGFNodes:
                self.sceneObjectGFPoints[name].mlab_source.set(x=p[0], y=p[1], z=p[2])

    def _drawTriSurface(self, name):

        self.scene.disable_render = True

        t = self.triSurface[name]
        P = t.v.T

        # triangulate vertices
        T = t.f

        # calc scalar
        S = self.triScalarData[name][self.triScalarList0]

        # render args
        renderArgs = self.triRenderArgs.get(name)
        if renderArgs is None:
            renderArgs = {}

        # draw
        if S is None:
            if 'color' not in renderArgs:
                renderArgs['color'] = self.defaultColor
            self.sceneObjectTri[name] = self.scene.mlab.triangular_mesh(P[0], P[1], P[2], T, name=name, **renderArgs)
        else:
            self.sceneObjectTri[name] = self.scene.mlab.triangular_mesh(P[0], P[1], P[2], T, scalars=S, name=name,
                                                                        **renderArgs)

        self.scene.disable_render = False

    def _evaluateScalarField(self, scalarFieldName, GFName):

        sf = self.scalarData[scalarFieldName]
        try:
            GFD = self.GFDSpecific[GFName]
        except KeyError:
            GFD = self.GFD

        scalarValues = sf.evaluateField(GFD)
        return scalarValues

    def _getGFScalarData(self, scalarName, GFName):

        if scalarName in list(self.scalarData.keys()):
            if self.scalarType[scalarName] == 'field':
                d = self._evaluateScalarField(scalarName, GFName)
            else:
                d = self.scalarData[scalarName]
            return d
        else:
            log.debug('ERROR: no scalar data called', scalarName)
            return None

    def drawElementBoundaries(self, name, GD, evaluatorMaker, nNodesElemMap, elemBasisMap, renderArgs):
        g = self.geometricFields[name]

        self.bCurves[name] = {}
        for elemN in list(g.ensemble_field_function.mesh.elements.keys()):
            self.bCurves[name][name + '_elem_' + str(elemN)] = g.makeElementBoundaryCurve(elemN, nNodesElemMap,
                                                                                          elemBasisMap)

        for b in self.bCurves[name]:
            evaluator = evaluatorMaker(self.bCurves[name][b], GD)
            self.addGeometricField(b, self.bCurves[name][b], evaluator, GD, renderArgs)
            self._drawGeometricField(b)

    def hideElementBoundaries(self, name):
        for b in self.bCurves[name]:
            SOb = self.sceneObjectGF.get(b)
            if SOb is not None:
                SOb.visible = False

    def showElementBoundaries(self, name):
        for b in self.bCurves[name]:
            SOb = self.sceneObjectGF.get(b)
            if SOb is not None:
                SOb.visible = True

    def _drawGeometricField(self, name):

        self.scene.disable_render = True

        g = self.geometricFields[name]
        evaluator = self.GFEvaluators[name]
        # evaluate mesh at triangle vertices
        # evaluate field at the desired element descritisation
        # P[0] = x, P[1] = y, P[2] = z
        P = evaluator(g.get_field_parameters().ravel())
        GFD = self.GFDSpecific[name]
        renderArgs = self.GFRenderArgs[name]
        if 'color' not in list(renderArgs.keys()):
            renderArgs['color'] = self.defaultColor

        # calc scalar
        S = None
        if self.GFScalar != 'none':
            if self.GFScalar == 'mean':
                K, H, k1, k2 = g.evaluate_curvature_in_mesh(GFD)
                S = H
            elif self.GFScalar == 'gaussian':
                K, H, k1, k2 = g.evaluate_curvature_in_mesh(GFD)
                S = K
            elif self.GFScalar == 'custom':
                S = self._getGFScalarData(self.GFScalarList0, name)

        # draw
        if g.ensemble_field_function.dimensions == 2:
            # triangulate vertices
            T = g.triangulator._triangulate(GFD)

            if self.mergeGFVertices:
                P, T, self.uniqueVertexIndices, vertMap = g.triangulator._mergePoints2(P.T)
                P = P.T

                if S is not None:
                    if S.shape[0] != P.shape[1]:
                        S = S[self.uniqueVertexIndices]

            if (S is None) or (S == 'none'):
                log.debug('S = None')
                self.sceneObjectGF[name] = self.scene.mlab.triangular_mesh(P[0], P[1], P[2], T, name=name, **renderArgs)
            else:
                # print S
                self.sceneObjectGF[name] = self.scene.mlab.triangular_mesh(P[0], P[1], P[2], T, scalars=S, name=name,
                                                                           **renderArgs)
        elif g.ensemble_field_function.dimensions == 1:
            self.sceneObjectGF[name] = g._draw_curve([GFD[0]], name=name, **renderArgs)

        self.sceneObjectGFPoints[name] = g._plot_points(glyph='sphere', label=None, scale=1.0, figure=None)
        if not self.displayGFNodes:
            self.sceneObjectGFPoints[name].visible = False

        self.scene.disable_render = False

    def drawGeometricFieldElementNumbers(self, name, textScale=5.0, textColor=(0, 0, 0)):
        g = self.geometricFields[name]
        ei, ex = g.get_element_numbers(coordinates=True)
        gElemLabels = [
            self.scene.mlab.text3d(ex[i][0], ex[i][1], ex[i][2], str(ei[i]), color=textColor, scale=textScale) for i in
            range(len(ei))]
        return gElemLabels

    def drawGeometricFieldNodeNumbers(self, name, textScale=5.0, textColor=(0, 0, 1)):
        P = self.geometricFields[name].get_all_point_positions()
        nodeLabels = [
            self.scene.mlab.text3d(
                P[i, 0], P[i, 1], P[i, 2], str(i), color=textColor, scale=textScale
            ) for i in range(len(P))
        ]
        return nodeLabels

    def _modeIndex_changed(self):
        q = self.sceneObjectGF.get('popQuiver')
        if q:
            lengthScale = 5.0
            # get mode and convert into vectors for each node
            m = self.pc.getMode(self.modeIndex).reshape((3, -1)).T

            # normalise magnitudes of vectors against the largest
            mag = scipy.sqrt((m * m).sum())
            magNorm = mag / mag.max()
            mNorm = (lengthScale * m / magNorm).T * self.modeVectorFlip

            nodeCoord = self.pc.reconstruct(self.pc.getWeightsBySD([self.modeIndex, ], [self.mode1, ]),
                                            [self.modeIndex, ]).reshape((3, -1))
            q.mlab_source.set(x=nodeCoord[0], y=nodeCoord[1], z=nodeCoord[2], u=mNorm[0], v=mNorm[1], w=mNorm[2])

    def _mode1_changed_old(self):
        # reconstruct field parameters
        params = self.pc.reconstruct(self.pc.getWeightsBySD([self.modeIndex, ], [self.mode1, ]), [self.modeIndex, ])
        # ~ V = self.GEval( params )
        # ~ # update vertices
        # ~ self.sceneObjectGF[self.pcGFName].mlab_source.set( x=V[0], y=V[1], z=V[2] )

        self.updateGeometricField(self.pcGFName, params=params)

        q = self.sceneObjectGF.get('popQuiver')
        if q:
            nodeCoord = params.reshape((3, -1))
            q.mlab_source.set(x=nodeCoord[0], y=nodeCoord[1], z=nodeCoord[2])

    def _mode1_changed(self):
        pc = self.PCs[self.PCList0]
        # reconstruct field parameters
        params = pc.reconstruct(pc.getWeightsBySD([self.modeIndex, ], [self.mode1, ]), [self.modeIndex, ])

        # get geometric object
        geomName, geomType = self.PCGeomList0.split('::')
        if geomType == 'gf':
            self.updateGeometricField(geomName, params=params)
            q = self.sceneObjectGF.get('popQuiver')
            if q:
                nodeCoord = params.reshape((3, -1))
                q.mlab_source.set(x=nodeCoord[0], y=nodeCoord[1], z=nodeCoord[2])
        elif geomType == 'data':
            self.updateData(geomName, coords=params.reshape([-1, 3]))
        elif geomType == 'tri':
            self.updateTriSurface(geomName, coords=params.reshape([-1, 3]))

    def _modeQuiver_fired(self):
        """ add quiver plot to show the direction of displacement at each
        node for the given mode
        """
        lengthScale = 2.0

        q = self.sceneObjectGF.get('popQuiver')
        if q:
            q.remove()
            del self.sceneObjectGF['popQuiver']
        else:
            # get mode and convert into vectors for each node
            m = self.pc.getMode(self.modeIndex).reshape((3, -1)).T

            # normalise magnitudes of vectors against the largest
            mag = scipy.sqrt((m * m).sum())
            magNorm = mag / mag.max()
            mNorm = (lengthScale * m / magNorm).T * self.modeVectorFlip

            # get current node coordinates
            nodeCoord = self.pc.reconstruct(self.pc.getWeightsBySD([self.modeIndex, ], [self.mode1, ]),
                                            [self.modeIndex, ]).reshape((3, -1))

            # make quiver plot
            # ~ self.sceneObjectGF['popQuiver'] = mlab.quiver3d( nodeCoord[0], nodeCoord[1], nodeCoord[2], mNorm[0], mNorm[1], mNorm[2] )
            self.sceneObjectGF['popQuiver'] = mlab.quiver3d(nodeCoord[0], nodeCoord[1], nodeCoord[2], mNorm[0],
                                                            mNorm[1], mNorm[2], scale_mode='none', scale_factor=5.0,
                                                            mode='arrow', color=(0, 0, 0))

    def _SFModeWeight_changed(self):
        # reconstruct field parameters
        params = self.SFPC.reconstruct(self.SFPC.getWeightsBySD([self.SFModeIndex, ], [self.SFModeWeight, ]),
                                       [self.SFModeIndex, ])
        # ~ V = self.GEval( params )
        # ~ # update vertices
        # ~ self.sceneObjectGF[self.pcGFName].mlab_source.set( x=V[0], y=V[1], z=V[2] )

        self.scalarData[self.SFPCFieldName].setFieldParameters(params[:, scipy.newaxis])
        s = self._evaluateScalarField(self.SFPCFieldName, self.GFList0)
        self.sceneObjectGF[self.GFList0].mlab_source.set(scalars=s)

    def _GFSFModeWeight_changed(self):
        # reconstruct combined pc weights
        weightsComb = self.GFSFPC.reconstruct(
            self.GFSFPC.getWeightsBySD([self.GFSFModeIndex, ], [self.GFSFModeWeight, ]), [self.GFSFModeIndex, ])
        # ~ V = self.GEval( params )
        # ~ # update vertices
        # ~ self.sceneObjectGF[self.pcGFName].mlab_source.set( x=V[0], y=V[1], z=V[2] )

        GFPCWeights = weightsComb[self.GFSFGFRows]
        SFPCWeights = weightsComb[self.GFSFSFRows]

        if self.GFSFPCComb:
            GFParams = GFPCWeights
        else:
            GFParams = self.pc.reconstruct(GFPCWeights, list(range(len(self.GFSFGFRows))))
        self.updateGeometricField(self.GFSFGFName, params=GFParams)

        if self.GFSFPCComb:
            SFParams = SFPCWeights
        else:
            SFParams = self.SFPC.reconstruct(SFPCWeights, list(range(len(self.GFSFSFRows))))
        self.scalarData[self.GFSFSFName].setFieldParameters(SFParams[:, scipy.newaxis])
        s = self._evaluateScalarField(self.GFSFSFName, self.GFList0)
        self.sceneObjectGF[self.GFList0].mlab_source.set(scalars=s)

    def saveImagesRevolve(self, nAngles, flipRoll=False, elevation=None,
                          distance=None, focalPoint=None, save=1):

        fileSuffix = str(self.saveImageFilename).split('.')[-1]
        filePrefix = '.'.join(str(self.saveImageFilename).split('.')[:-1])

        # if extra parameters not provided, use the current
        azimuth, e, d, f = self.scene.mlab.view()
        if elevation is None:
            elevation = e
        if distance is None:
            distance = d
        if focalPoint is None:
            focalPoint = f

        azimuths = azimuth + scipy.linspace(0.0, 360.0, nAngles)
        azimuths = scipy.where(azimuths > 360.0, azimuths - 360.0, azimuths)

        for i, a in enumerate(azimuths):
            self.scene.mlab.view(
                a, elevation, distance, focalPoint, reset_roll=True
            )
            if flipRoll:
                self.scene.mlab.roll(-self.scene.mlab.roll())
            if save:
                filename = filePrefix + '_%(i)05i' % {'i': i} + '.' + fileSuffix
                self.scene.mlab.savefig(
                    filename,
                    size=(int(self.saveImageWidth), int(self.saveImageLength))
                )

    def saveImagesModalDeformation(self, modeN, nImages, sdRange, modeType, save=1, startNumber=0, loop=False):

        fileSuffix = str(self.saveImageFilename).split('.')[-1]
        filePrefix = '.'.join(str(self.saveImageFilename).split('.')[:-1])

        if modeType == 'shape':
            self.modeIndex = modeN
        elif modeType == 'scalar':
            self.SFModeIndex = modeN
        elif modeType == 'shapescalar':
            self.GFSFModeIndex = modeN
        else:
            raise ValueError('unknown modeType')

        # self.modeIndex = modeN
        # ~ self._modeIndex_changed()

        if not loop:
            sdList = scipy.linspace(sdRange[0], sdRange[1], nImages)
        else:
            sdList = scipy.hstack([scipy.linspace(sdRange[0], sdRange[1], nImages),
                                   scipy.linspace(sdRange[1], sdRange[0], nImages)[1:]])

        for i, sd in enumerate(sdList):

            if modeType == 'shape':
                self.mode1 = sd
            elif modeType == 'scalar':
                self.SFModeWeight = sd
            elif modeType == 'shapescalar':
                self.GFSFModeWeight = sd
            else:
                raise ValueError('unknown modeType')
            # ~ self._mode1_changed()

            if save:
                filename = filePrefix + '_%(i)05i' % {'i': i + startNumber} + '.' + fileSuffix
                self.scene.mlab.savefig(filename, size=(int(self.saveImageWidth), int(self.saveImageLength)))

    def saveImagesGeneric(self, func, save=1):
        """
        func is a generator function that changes something in the scene
        and returns a number that will be used in the filename
        """
        fileSuffix = str(self.saveImageFilename).split('.')[-1]
        filePrefix = '.'.join(str(self.saveImageFilename).split('.')[:-1])

        for i in func:
            if save:
                filename = filePrefix + '_%(i)05i' % {'i': i} + '.' + fileSuffix
                self.scene.mlab.savefig(filename, size=(int(self.saveImageWidth), int(self.saveImageLength)))

    def scrollSlicesGenerator(self, slices, plane):

        for i, s in enumerate(slices):
            self.imagePlaneSliceIndex(s, plane)
            yield i

    def storeView(self):
        self.view = self.scene.mlab.view()
        self.roll = self.scene.mlab.roll()
        self.move = self.scene.mlab.move()

    def restoreView(self):
        # ~ self.scene.mlab.view( self.view[0],self.view[1],self.view[2],self.view[3] )
        self.scene.mlab.view(*self.view)
        self.scene.mlab.roll(self.roll)

    def setView(self, view, roll):
        self.view = view
        self.roll = roll
        self.restoreView()

    def setViewFocalPoint(self, f):
        if self.view is None:
            self.storeView()

        viewNew = list(self.view)
        viewNew[3] = f

        self.view = viewNew

    def saveSceneVRML(self, filename):
        renWin = self.scene.render_window
        vrmlEx = tvtk.VRMLExporter()
        vrmlEx.input = renWin
        vrmlEx.file_name = filename
        vrmlEx.write()

    # =========================================================================#
    # Lower limb atlas
    # =========================================================================#
    def setLowerLimbAtlas(self, LL, gfEvalMaker):
        self.LL = LL
        self.LLParams = list(LL._neutral_params)
        self.LLParams[0] = [0.0, ] * 5
        self.LLParams[1] = scipy.arange(5)

        for mn, m in list(self.LL.models.items()):
            mgfeval = gfEvalMaker(m.gf, self.GFD)
            self.addGeometricField('LL_' + mn, m.gf, mgfeval, self.GFD)

    def _LLPC0_changed(self):
        self.LLParams[0][0] = self.LLPC0
        self.LL.update_all_models(*self.LLParams)
        self._updateLLGF()

    def _LLPC1_changed(self):
        self.LLParams[0][1] = self.LLPC1
        self.LL.update_all_models(*self.LLParams)
        self._updateLLGF()

    def _LLPC2_changed(self):
        self.LLParams[0][2] = self.LLPC2
        self.LL.update_all_models(*self.LLParams)
        self._updateLLGF()

    def _LLPC3_changed(self):
        self.LLParams[0][3] = self.LLPC3
        self.LL.update_all_models(*self.LLParams)
        self._updateLLGF()

    def _LLPC4_changed(self):
        self.LLParams[0][4] = self.LLPC4
        self.LL.update_all_models(*self.LLParams)
        self._updateLLGF()

    def _LLHJC1_changed(self):
        self.LLParams[3][0] = scipy.deg2rad(self.LLHJC1)
        self.LL.update_femur(self.LLParams[3])
        self.LL.update_tibiafibula(self.LLParams[4])
        self.LL.update_patella()
        self._updateLLGF()

    def _LLHJC2_changed(self):
        self.LLParams[3][1] = scipy.deg2rad(self.LLHJC2)
        self.LL.update_femur(self.LLParams[3])
        self.LL.update_tibiafibula(self.LLParams[4])
        self.LL.update_patella()
        self._updateLLGF()

    def _LLHJC3_changed(self):
        self.LLParams[3][2] = scipy.deg2rad(self.LLHJC3)
        self.LL.update_femur(self.LLParams[3])
        self.LL.update_tibiafibula(self.LLParams[4])
        self.LL.update_patella()
        self._updateLLGF()

    def _LLKnee1_changed(self):
        self.LLParams[4][0] = scipy.deg2rad(self.LLKnee1)
        self.LL.update_tibiafibula(self.LLParams[4])
        self.LL.update_patella()
        self._updateLLGF()

    def _updateLLGF(self):
        for mn in self.LL.models:
            self.updateGeometricField('LL_' + mn)

    # =========================================================================#
    # 2 rigid-body
    # =========================================================================#
    def setRigidBodies(self, meshparent, meshchild, *joints):
        self.rb_parent = meshparent
        self.rb_child = meshchild
        self.joints = joints

        self.addTri('Parent body', self.rb_parent)
        self.addTri('Child body', self.rb_child)

    def _updateRigidBodies(self):
        self.updateTriSurface('Parent body')
        self.updateTriSurface('Child body')

    def _RBAX0_changed(self):
        if len(self.joints) > 0:
            self.joints[0].angle = scipy.deg2rad(self.RBAX0)
            self._updateRigidBodies()
        else:
            log.debug('Axis 0 not defined')

    def _RBAX1_changed(self):
        if len(self.joints) > 1:
            self.joints[1].angle = scipy.deg2rad(self.RBAX1)
            self._updateRigidBodies()
        else:
            log.debug('Axis 1 not defined')

    def _RBAX2_changed(self):
        if len(self.joints) > 2:
            self.joints[2].angle = scipy.deg2rad(self.RBAX2)
            self._updateRigidBodies()
        else:
            log.debug('Axis 2 not defined')

    # =========================================================================#
    # image picker
    # =========================================================================#
    def _ipw_pick_callback(self, obj, evt):
        img_coords = scipy.around(obj.GetCurrentCursorPosition()).astype(int)
        img_val = obj.GetCurrentImageValue()
        log.debug(('picked location: {}, value: {}'.format(img_coords, img_val)))
        self._ipw_picked_obj = obj
        self._ipw_picked_points.append(img_coords)

    def clear_ipw_picked_points(self):
        self._ipw_picked_points = []

    def get_ipw_picked_points(self):
        return list(self._ipw_picked_points)

    def addText3D(self, name, text, origin, offset, charWidth=0.01, lineWidth=0.2):
        """
        Adds text with a line to a point
        """
        textOrigin = scipy.array(origin) + scipy.array(offset)
        textLine = scipy.array([origin, textOrigin]).T
        self.scene.mlab.text(textOrigin[0], textOrigin[1], text, z=textOrigin[2], width=len(text) * charWidth,
                             name='text_' + name)
        self.scene.mlab.plot3d(textLine[0], textLine[1], textLine[2], tube_radius=lineWidth, name='textline_' + name)


def find_shaft_frame_gen(V, nFrames, direction='up', shaftMode='reveal'):
    if direction == 'up':
        slices = scipy.linspace(V.I.shape[2] - 1, 0, nFrames).astype(int)
    elif direction == 'down':
        slices = scipy.linspace(0, V.I.shape[2] - 1, nFrames).astype(int)

    a0, e, d, f = V.scene.mlab.view()
    azimuths = a0 + scipy.linspace(0.0, 360.0, nFrames)
    azimuths = scipy.where(azimuths > 360.0, azimuths - 360.0, azimuths)

    if shaftMode == 'reveal':
        shaft0 = V.data['shaft0'].copy()
        shaft1 = V.data['shaft1'].copy()

    for i in range(nFrames):
        # set slice index
        V.imagePlaneSliceIndex(slices[i], plane=0)

        # set shaft data
        if shaftMode == 'reveal':
            V.storeView()
            V.removeData('shaft0')
            V.removeData('shaft1')

            if direction == 'up':
                V.addData('shaft0', shaft0[scipy.where(shaft0[:, 2] > slices[i])[0]])
                V.addData('shaft1', shaft1[scipy.where(shaft1[:, 2] > slices[i])[0]])
            elif direction == 'down':
                V.addData('shaft0', shaft0[scipy.where(shaft0[:, 2] < slices[i])[0]])
                V.addData('shaft1', shaft1[scipy.where(shaft1[:, 2] < slices[i])[0]])

            V._drawData('shaft0')
            V._drawData('shaft1')
            V.restoreView()

        # set scene rotation
        V.scene.mlab.view(azimuths[i], e, d, f, reset_roll=True)

        yield i
