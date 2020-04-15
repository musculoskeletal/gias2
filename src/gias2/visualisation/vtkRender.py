"""
FILE: vtkRender.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: classes and functions for rendering numpy/scipy arrays using VTK

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import vtk
from numpy import uint8, array

raise DeprecationWarning('Use classes in gias2.mesh.vtktools')


class Colours:
    def __init__(self):
        self.colours = dict()

        red = vtk.vtkProperty()
        red.SetColor(1.0, 0.0, 0.0);
        self.colours['red'] = red

        green = vtk.vtkProperty()
        green.SetColor(0.0, 1.0, 0.0);
        self.colours['green'] = green

        blue = vtk.vtkProperty()
        blue.SetColor(0.0, 0.0, 1.0);
        self.colours['blue'] = blue

        magenta = vtk.vtkProperty()
        magenta.SetColor(1.0, 0.0, 1.0);
        self.colours['magenta'] = magenta

        yellow = vtk.vtkProperty()
        yellow.SetColor(1.0, 1.0, 0.0);
        self.colours['yellow'] = yellow

        cyan = vtk.vtkProperty()
        cyan.SetColor(0.0, 1.0, 1.0);
        self.colours['cyan'] = cyan

    def getColour(self, colourStr):
        return self.colours[colourStr]


class VtkImageVolumeRenderer:
    def __init__(self, image=None):

        self.renderWindowSize = 400
        self.bgColour = [0.0, 0.0, 0.0]
        self.colours = Colours()
        self.imageImporter = vtk.vtkImageImport()
        self.image = None

        if image != None:
            # get image into right format
            # ~ self.image = array( image, dtype = uint8 )
            self.image = image

            # import image array into vtk
            self._importImage()

        # ~ self.volumeList = []
        # ~ self.actorList = []
        self.CoMActors = []
        self.PDActors = []
        self.nodeActors = []

    def _importImage(self):
        # import image into vtk
        imageString = self.image.astype(uint8).tostring()
        self.imageImporter.CopyImportVoidPointer(imageString, len(imageString))
        self.imageImporter.SetDataScalarTypeToUnsignedChar()
        self.imageImporter.SetNumberOfScalarComponents(1)

        # set imported image size
        S = self.image.shape
        self.imageImporter.SetDataExtent(0, S[2] - 1, 0, S[1] - 1, 0, S[0] - 1)
        self.imageImporter.SetWholeExtent(0, S[2] - 1, 0, S[1] - 1, 0, S[0] - 1)

    def setImage(self, image):
        # ~ self.image = array( image, dtype = uint8 )
        self.image = image
        self._importImage()
        self.CoMActors = []
        self.PDActors = []
        self.nodeActors = []

    # Set the position of a sphere to mark the image CoM, and the
    # principal directions (optional)
    # CoM: array-like eg [x,y,z]
    # PD: list or tuple of 3 unit-vectors
    # lineScale: a list of 3 line scaling factors
    def setCoM(self, inputCoM, PD=None, lineScale=None):

        CoM = list(inputCoM)
        CoM.reverse()
        CoMSphere = vtk.vtkSphereSource()
        CoMSphere.SetCenter(CoM)
        CoMSphere.SetRadius(2.0)
        CoMSphere.SetPhiResolution(16)
        CoMSphere.SetThetaResolution(16)
        CoMSphereMapper = vtk.vtkPolyDataMapper()
        CoMSphereMapper.SetInput(CoMSphere.GetOutput())
        CoMSphereActor = vtk.vtkActor()
        CoMSphereActor.SetProperty(self.colours.getColour('magenta'))
        CoMSphereActor.SetMapper(CoMSphereMapper)

        self.CoMActors.append(CoMSphereActor)

        if PD != None:

            # line scaling from unit length
            if lineScale == None:
                s = array([max(self.image.shape) * 0.5] * 3)
            else:
                s = array(lineScale) * max(self.image.shape) * (1 / max(lineScale))

            # ~ s = list( s )
            # ~ s.reverse()
            for d in range(0, len(PD)):
                v = list(PD[:, d])
                v.reverse()
                v = array(v)
                startPoint = CoM - s[d] * v
                endPoint = CoM + s[d] * v

                self.addLine(startPoint, endPoint)

    def addNode(self, coord, colourStr='red'):

        coord = list(coord)
        coord.reverse()
        node = vtk.vtkSphereSource()
        node.SetCenter(coord)
        node.SetRadius(1.5)
        node.SetPhiResolution(16)
        node.SetThetaResolution(16)
        nodeMapper = vtk.vtkPolyDataMapper()
        nodeMapper.SetInput(node.GetOutput())
        nodeActor = vtk.vtkActor()
        nodeActor.SetProperty(self.colours.getColour(colourStr))
        nodeActor.SetMapper(nodeMapper)

        self.nodeActors.append(nodeActor)

    def addLine(self, p1, p2, colourStr='red'):
        """ add a line to the scene, between points p1 and p2
        """

        line = vtk.vtkLineSource()
        line.SetPoint1(p1[0], p1[1], p1[2])
        line.SetPoint2(p2[0], p2[1], p2[2])
        line.SetResolution(21)
        lineMapper = vtk.vtkPolyDataMapper()
        lineMapper.SetInput(line.GetOutput())
        lineActor = vtk.vtkActor()
        lineActor.SetProperty(self.colours.getColour(colourStr))
        lineActor.SetMapper(lineMapper)

        self.PDActors.append(lineActor)

    def addPlane(self, origin, normal):

        plane = vtk.vtkPlaneSource()
        plane.SetOrigin((0.0, 0.0, 0.0))
        plane.SetPoint1((50.0, 0.0, 0.0))
        plane.SetPoint2((0.0, 50.0, 0.0))

        plane.SetCenter(origin[::-1])
        plane.SetNormal(normal[::-1])
        planeMapper = vtk.vtkPolyDataMapper()
        planeMapper.SetInputConnection(plane.GetOutputPort())
        planeActor = vtk.vtkActor()
        planeActor.SetMapper(planeMapper)

        self.PDActors.append(planeActor)

    def clearNodes(self):
        self.nodeActors = []

    def renderVolume(self, cRange=[0, 255], oRange=[0, 255]):
        # volume rendering

        # Volume mapper 
        volumeMapper = vtk.vtkVolumeRayCastMapper()
        volumeMapper.SetInput(self.imageImporter.GetOutput())
        compositeFunc = vtk.vtkVolumeRayCastCompositeFunction()
        volumeMapper.SetVolumeRayCastFunction(compositeFunc)

        # Colour transfer functions
        colorFunc = vtk.vtkColorTransferFunction()
        colorFunc.AddRGBPoint(cRange[0], 0.0, 0.0, 0.0)
        colorFunc.AddRGBPoint(cRange[1], 1.0, 1.0, 1.0)

        # Opacity transfer functions
        opacityFunc = vtk.vtkPiecewiseFunction()
        opacityFunc.AddPoint(oRange[0], 0.0)
        # ~ opacity_transfer_func.AddPoint( 99, 0.0 )
        # ~ opacity_transfer_func.AddPoint( 250, 0.0 )
        opacityFunc.AddPoint(oRange[1], 0.1)

        # Volume properties
        volumeProperties = vtk.vtkVolumeProperty()
        volumeProperties.SetColor(colorFunc)
        volumeProperties.SetScalarOpacity(opacityFunc)

        # VTK volume
        volume = vtk.vtkVolume()
        volume.SetMapper(volumeMapper)
        volume.SetProperty(volumeProperties)

        # ~ self.volumeList.append( volume )

        self._render(volumeList=[volume])

    def renderContour(self, contourValueList):
        # render polydata contour surfaces at iso values defined in
        # list contourValueList

        contourExtractor = vtk.vtkContourFilter()
        contourExtractor.SetInput(self.imageImporter.GetOutput())
        # set contour values
        for i in range(0, len(contourValueList)):
            contourExtractor.SetValue(i, contourValueList[i])

        contourExtractor.Update()

        map = vtk.vtkPolyDataMapper()
        map.SetInput(contourExtractor.GetOutput())
        actor = vtk.vtkActor()
        actor.SetMapper(map)
        actor.GetProperty().SetColor(0.9, 0.9, 0.7)

        # ~ self.actorList.append( actor )

        self._render(actorList=[actor])

    def _render(self, actorList=None, volumeList=None):

        # axes
        axes = vtk.vtkAxesActor()
        axes.SetAxisLabels(0)
        axes.SetTotalLength(50, 100, 150)
        axes.SetConeRadius(0.1)

        # bounding box
        outline = vtk.vtkOutlineFilter()
        outline.SetInput(self.imageImporter.GetOutput())
        outlineMapper = vtk.vtkPolyDataMapper()
        outlineMapper.SetInput(outline.GetOutput())
        outlineActor = vtk.vtkActor()
        outlineActor.SetMapper(outlineMapper)
        outlineActor.GetProperty().SetColor(0.0, 0.0, 1.0)

        # renderer
        renderer = vtk.vtkRenderer()
        renderer.SetBackground(self.bgColour[0], self.bgColour[1], self.bgColour[2])

        renderer.AddActor(axes)
        renderer.AddActor(outlineActor)

        # add other actors
        if actorList:
            for actor in actorList:
                renderer.AddActor(actor)

        # add other volumes
        if volumeList:
            for volume in volumeList:
                renderer.AddVolume(volume)

        # add node spheres
        if len(self.nodeActors) > 0:
            for node in self.nodeActors:
                renderer.AddActor(node)

        # add CoM spheres
        if len(self.CoMActors) > 0:
            for CoM in self.CoMActors:
                renderer.AddActor(CoM)

        # add principal direction lines
        if len(self.PDActors) > 0:
            for PD in self.PDActors:
                renderer.AddActor(PD)

        # render window
        renderWindow = vtk.vtkRenderWindow()
        renderWindow.AddRenderer(renderer)
        renderWindow.SetSize(self.renderWindowSize, self.renderWindowSize)

        # render window interactor
        interactor = vtk.vtkRenderWindowInteractor()
        interactor.SetRenderWindow(renderWindow)

        renderWindow.AddObserver("AbortCheckEvent", self._exitCheck)

        # render and start interaction.
        interactor.Initialize()
        renderWindow.Render()
        interactor.Start()

    def clearActors(self):
        self.clearCoMs()
        self.clearPDs()

        return

    def clearCoMs(self):
        """ removes all CoM markers """
        self.CoMActors = []
        return

    def clearPDs(self):
        """ removes all principal direction lines """
        self.PDActors = []

    def _exitCheck(self, obj, event):
        if obj.GetEventPending() != 0:
            obj.SetAbortRender(1)
