"""
FILE: vtktool.py
LAST MODIFIED: 05-08-2018
DESCRIPTION: Classes and functions for working with vtkpolydata

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import copy
import logging
import pickle
import warnings
from os import path

import sys
import vtk
from numpy import zeros, array, uint8, int16, ones, newaxis, ascontiguousarray
from vtk.util import numpy_support

from gias2.mesh import plywriter
from gias2.mesh import simplemesh

log = logging.getLogger(__name__)


class Writer(object):
    """Class for writing polygons to file formats supported by VTK.
    """

    def __init__(self, **kwargs):
        """Keyword arguments:
        filename: output filename
        polydata: vtkPolydata instance
        v: array of vertices coordinates
        f: list of faces composed of lists of vertex indices
        vn: array of vertex normals
        rw: vtkRenderWindow instance
        colour: 3-tuple of colour (only works for ply)
        vcolour: 3-tuple of colour for each vertex
        fcolour: 3-tuple of colour for each face
        ascii: boolean, write in ascii (True) or binary (False)
        """
        self.filename = kwargs.get('filename')
        if self.filename is not None:
            self._parse_format()
        self._polydata = kwargs.get('polydata')
        self._vertices = kwargs.get('v')
        self._faces = kwargs.get('f')
        self._vertex_normals = kwargs.get('vn')
        self._render_window = kwargs.get('rw')
        self._colour = kwargs.get('colour')
        self._vertex_colour = kwargs.get('vcolour')
        self._face_colour = kwargs.get('fcolour')
        # self._field_data = kwargs.get('field')
        self._write_ascii = kwargs.get('ascii')

    def setFilename(self, f):
        self.filename = f
        self._parse_format()

    def _parse_format(self):
        self.file_prefix, self.file_ext = path.splitext(self.filename)
        self.file_ext = self.file_ext.lower()

    def _make_polydata(self):
        self._polydata = polygons2Polydata(
            self._vertices,
            self._faces,
            vcolours=self._vertex_colour,
            fcolours=self._face_colour,
            vnormals=self._vertex_normals,
        )

    def _make_render_window(self):
        if self._polydata is None:
            self._make_polydata()
        ply_mapper = vtk.vtkPolyDataMapper()
        if vtk.VTK_MAJOR_VERSION < 6:
            ply_mapper.SetInput(self._polydata)
        else:
            ply_mapper.SetInputDataObject(self._polydata)
        ply_actor = vtk.vtkActor()
        ply_actor.SetMapper(ply_mapper)

        ren1 = vtk.vtkRenderer()
        self._render_window = vtk.vtkRenderWindow()
        self._render_window.AddRenderer(ren1)
        ren1.AddActor(ply_actor)

    def write(self, filename=None, ascenc=True):
        if filename is not None:
            self.filename = filename

        filePrefix, fileExt = path.splitext(self.filename)
        fileExt = fileExt.lower()
        if fileExt == '.obj':
            self.writeOBJ()
        elif fileExt == '.wrl':
            self.writeVRML()
        elif fileExt == '.stl':
            self.writeSTL(ascenc=ascenc)
        elif fileExt == '.ply':
            self.writePLY(ascenc=ascenc)
        elif fileExt == '.vtp':
            self.writeVTP(ascenc=ascenc)
        else:
            raise ValueError('unknown file extension')

    def writeOBJ(self, filename=None):
        if filename is not None:
            self.filename = filename
        if self._render_window is None:
            self._make_render_window()

        w = vtk.vtkOBJExporter()
        w.SetRenderWindow(self._render_window)
        w.SetFilePrefix(path.splitext(self.filename)[0])
        w.Write()

    def writePLY(self, filename=None, ascenc=True):
        if filename is not None:
            self.filename = filename

        # if there are vnormals, have to use the gias2 writer since the vtk
        # writer does not write normals
        if self._vertex_normals is not None:
            if not ascenc:
                warnings.warn(
                    'Using GIAS2 PLYWriter, binary encoding not supported.',
                    UserWarning
                )
            w = plywriter.PLYWriter(
                v=self._vertices, f=self._faces, filename=self.filename,
                vn=self._vertex_normals, vcolours=self._vertex_colour
            )
            w.write()
            return

        if self._polydata is None:
            self._make_polydata()

        w = vtk.vtkPLYWriter()
        if vtk.VTK_MAJOR_VERSION < 6:
            w.SetInput(self._polydata)
        else:
            w.SetInputDataObject(self._polydata)
        w.SetFileName(self.filename)
        if ascenc:
            w.SetFileTypeToASCII()
        else:
            w.SetFileTypeToBinary()
        w.SetDataByteOrderToLittleEndian()
        # w.SetColorModeToUniformCellColor()
        # w.SetColor(255, 0, 0)

        if self._vertex_colour is not None:
            w.SetArrayName('colours')

        w.Write()

    def writeSTL(self, filename=None, ascenc=True):
        if filename is not None:
            self.filename = filename
        if self._polydata is None:
            self._make_polydata()

        w = vtk.vtkSTLWriter()
        if vtk.VTK_MAJOR_VERSION < 6:
            w.SetInput(self._polydata)
        else:
            w.SetInputDataObject(self._polydata)
        w.SetFileName(self.filename)
        if ascenc:
            w.SetFileTypeToASCII()
        else:
            w.SetFileTypeToBinary()
        w.Write()

    def writeVRML(self, filename=None):
        if filename is not None:
            self.filename = filename
        if self._render_window is None:
            self._make_render_window()

        w = vtk.vtkVRMLExporter()
        w.SetRenderWindow(self._render_window)
        w.SetFileName(self.filename)
        w.Write()

    def writeVTP(self, filename=None, ascenc=True, xmlenc=True):
        if filename is not None:
            self.filename = filename
        if self._polydata is None:
            self._make_polydata()

        if xmlenc:
            w = vtk.vtkXMLPolyDataWriter()
        else:
            w = vtk.vtkPolyDataWriter()

        # if (vtk.VTK_MAJOR_VERSION<6) and not xmlenc:
        if (vtk.VTK_MAJOR_VERSION < 6):
            w.SetInput(self._polydata)
        else:
            w.SetInputDataObject(self._polydata)
        w.SetFileName(self.filename)

        if not xmlenc:
            if ascenc:
                w.SetFileTypeToASCII()
            else:
                w.SetFileTypeToBinary()

        w.Write()


class Reader(object):
    """Class for reading polygon files of various formats
    """

    def __init__(self, **kwargs):
        self.filename = kwargs.get('filename')
        self.verbose = kwargs.get('verbose', True)
        self._points = None
        self._triangles = None
        self._vertexNormals = None
        self._nPoints = None
        self._nFaces = None
        self._nVertexNormals = None
        self._dimensions = None
        self.polydata = None

    def setFilename(self, filename):
        self.filename = filename

    def read(self, filename=None):
        if filename is not None:
            self.filename = filename

        filePrefix, fileExt = path.splitext(self.filename)
        fileExt = fileExt.lower()
        if fileExt == '.obj':
            self.readOBJ()
        elif fileExt == '.wrl':
            self.readVRML()
        elif fileExt == '.stl':
            self.readSTL()
        elif fileExt == '.ply':
            self.readPLY()
        elif fileExt == '.vtp':
            self.readVTP()
        else:
            raise ValueError('unknown file extension in {}'.format(self.filename))

    def readVRML(self, filename=None, actor=0):
        if filename is not None:
            self.filename = filename
        r = vtk.vtkVRMLImporter()
        r.SetFileName(self.filename)
        r.Update()
        actors = r.GetRenderer().GetActors()
        actors.InitTraversal()
        i = 0
        while i != actor:
            actors.GetNextActor()
            i += 1

        self.polydata = actors.GetNextActor().GetMapper().GetInput()

        if self.polydata.GetPoints() == None:
            raise IOError('file not loaded {}'.format(self.filename))
        else:
            self._load()

    def readOBJ(self, filename=None):
        if filename is not None:
            self.filename = filename

        r = vtk.vtkOBJReader()
        r.SetFileName(self.filename)
        r.Update()
        self.polydata = r.GetOutput()

        if self.polydata.GetPoints() == None:
            raise IOError('file not loaded {}'.format(self.filename))
        else:
            self._load()

    def readPLY(self, filename=None):
        if filename is not None:
            self.filename = filename

        r = vtk.vtkPLYReader()
        r.SetFileName(self.filename)
        r.Update()
        self.polydata = r.GetOutput()

        if self.polydata.GetPoints() == None:
            raise IOError('file not loaded {}'.format(self.filename))
        else:
            self._load()

    def readSTL(self, filename=None):
        if filename is not None:
            self.filename = filename

        r = vtk.vtkSTLReader()
        r.SetFileName(self.filename)
        r.Update()
        self.polydata = r.GetOutput()

        if self.polydata.GetPoints() == None:
            raise IOError('file not loaded {}'.format(self.filename))
        else:
            self._load()

    def readVTP(self, filename=None):
        if filename is not None:
            self.filename = filename

        if self._isXML(self.filename):
            r = vtk.vtkXMLPolyDataReader()
        else:
            r = vtk.vtkPolyDataReader()
        r.SetFileName(self.filename)
        r.Update()
        self.polydata = r.GetOutput()

        if self.polydata.GetPoints() == None:
            raise IOError('file not loaded {}'.format(self.filename))
        else:
            self._load()

    def _isXML(self, f):
        """Check if file is an xml file
        """
        with open(f, 'r') as fp:
            l = fp.readline()

        if l[0] == '<':
            return True
        else:
            return False

    def _load(self):
        self._loadPoints()
        self._loadTriangles()
        self._loadVertexNormals()

    def _loadPoints(self):
        P = self.polydata.GetPoints().GetData()
        self._dimensions = P.GetNumberOfComponents()
        self._nPoints = P.GetNumberOfTuples()

        if self.verbose:
            log.debug('loading %(np)i points in %(d)i dimensions' % {'np': self._nPoints, 'd': self._dimensions})

        self._points = numpy_support.vtk_to_numpy(P)

        # if self._dimensions==1:
        #     self._points = array([P.GetTuple1(i) for i in range(self._nPoints)])
        # elif self._dimensions==2:
        #     self._points = array([P.GetTuple2(i) for i in range(self._nPoints)])
        # elif self._dimensions==3:
        #     self._points = array([P.GetTuple3(i) for i in range(self._nPoints)])
        # elif self._dimensions==4:
        #     self._points = array([P.GetTuple4(i) for i in range(self._nPoints)])
        # elif self._dimensions==9:
        #     self._points = array([P.GetTuple9(i) for i in range(self._nPoints)])

    def _loadTriangles(self):
        polyData = self.polydata.GetPolys().GetData()
        # X = [int(polyData.GetTuple1(i)) for i in range(polyData.GetNumberOfTuples())]
        X = numpy_support.vtk_to_numpy(polyData)

        # assumes that faces are triangular
        X = array(X).reshape((-1, 4))
        self._nFaces = X.shape[0]
        self._triangles = X[:, 1:].copy()

        if self.verbose:
            log.debug('loaded %(f)i faces' % {'f': self._nFaces})

    def _loadVertexNormals(self):
        ptsNormals = self.polydata.GetPointData().GetNormals()
        if ptsNormals is not None:
            self._vertexNormals = numpy_support.vtk_to_numpy(ptsNormals)
            self._nVertexNormals = self._vertexNormals.shape[0]
        else:
            self._nVertexNormals = 0
            self._vertexNormals = None

        if self.verbose:
            log.debug('loaded %(f)i vertex normals' % {'f': self._nVertexNormals})

    def getSimplemesh(self):
        S = simplemesh.SimpleMesh(self._points, self._triangles)
        # S.calcFaceProperties()

        if self._vertexNormals is not None:
            S.vertexNormals = array(self._vertexNormals)
            S.hasVertexNormals = True

        return S


def savepoly(sm, filename, ascenc=True):
    w = Writer(v=sm.v, f=sm.f, vn=sm.vertexNormals)
    w.write(filename, ascenc=ascenc)


def loadpoly(filename, verbose=False):
    r = Reader(verbose=verbose)
    r.read(filename)
    return r.getSimplemesh()


class PolydataReader:
    def __init__(self, filename):
        self.fileName = filename

    def loadData(self):
        log.debug("opening", self.fileName)
        self.file = open(self.fileName, "r")
        self.file.seek(0, 2)
        self.fileEnd = self.file.tell()

        # determine number of points in the data
        self.file.seek(0)
        foundNumPoints = False
        while (foundNumPoints == False) and (self.file.tell() < self.fileEnd):
            line = self.file.readline()
            if line.find("POINTS") > -1:
                foundNumPoints = True
                line = line.split()
                self.numPoints = int(line[1])
                log.debug(str(self.numPoints) + " points in file")

        if foundNumPoints == False:
            log.debug("Unable to find number of points!")

            # determine number of faces in the data
        self.file.seek(0)
        foundNumFaces = False
        while (foundNumFaces == False) and (self.file.tell() < self.fileEnd):
            line = self.file.readline()
            if line.find("POLYGONS") > -1:
                foundNumFaces = True
                line = line.split()
                self.numFaces = int(line[1])
                log.debug(str(self.numFaces) + " faces in file")

        if foundNumFaces == False:
            log.debug("Unable to find number of Faces!")

    def getPoints(self):

        try:
            self.numPoints
        except AttributeError:
            log.debug("number of points unknown. Run .loadData() first")
        else:
            self.file.seek(0)
            foundPoints = False
            while (foundPoints == False) and (self.file.tell() < self.fileEnd):

                # find points section header
                line = self.file.readline()
                if line.find("POINTS") > -1:
                    foundPoints = True
                    log.debug("getting points...")
                    self.points = zeros([self.numPoints, 3], dtype=float)
                    pointCounter = 0

                    while pointCounter < self.numPoints:
                        pointLine = array(self.file.readline().split(), dtype=float)

                        for a in range(0, len(pointLine) / 3):
                            self.points[pointCounter, :] = pointLine[0 + 3 * a:3 + 3 * a]
                            pointCounter += 1

                    log.debug("Got " + str(self.points.shape[0]) + " points")
                    return self.points

            if foundPoints == False:
                log.debug("No points found!")

    def getPointNormals(self):

        try:
            self.numPoints
        except AttributeError:
            log.debug("number of points unknown. Run .loadData() first")
        else:
            self.file.seek(0)
            foundNormals = False
            while (foundNormals == False) and (self.file.tell() < self.fileEnd):

                # find points section header
                line = self.file.readline()
                if line.find("NORMALS") > -1:
                    foundNormals = True
                    log.debug("getting point normals...")
                    self.pointNormals = zeros([self.numPoints, 3], dtype=float)
                    normalCounter = 0

                    while normalCounter < self.numPoints:
                        normalLine = array(self.file.readline().split(), dtype=float)

                        for a in range(0, len(normalLine) / 3):
                            self.pointNormals[normalCounter, :] = normalLine[0 + 3 * a:3 + 3 * a]
                            normalCounter += 1

                    log.debug("Got " + str(self.pointNormals.shape[0]) + " point normals")
                    return self.pointNormals

            if foundNormals == False:
                log.debug("No point Normals found!")

    def getFaces(self):
        try:
            self.numFaces
        except AttributeError:
            log.debug("number of faces unknown. Run .loadData() first")
        else:
            self.file.seek(0)
            foundFaces = False
            while (foundFaces == False) and (self.file.tell() < self.fileEnd):

                # find points section header
                line = self.file.readline()
                if line.find("POLYGONS") > -1:
                    foundFaces = True
                    log.debug("getting faces...")
                    self.faces = zeros([self.numFaces, 3], dtype=int)
                    faceCounter = 0

                    while faceCounter < self.numFaces:
                        faceLine = array(self.file.readline().split(), dtype=int)

                        self.faces[faceCounter, :] = faceLine[1:4]
                        faceCounter += 1

                    log.debug("Got " + str(self.faces.shape[0]) + " faces")
                    return self.faces

            if foundFaces == False:
                log.debug("No faces found!")

    def getCurvature(self):

        try:
            self.numPoints
        except AttributeError:
            log.debug("number of points unknown. Run .loadData() first")
        else:
            self.file.seek(0)
            foundCurv = False
            while (foundCurv == False) and (self.file.tell() < self.fileEnd):

                # find curvature section header
                line = self.file.readline()
                if line.find("Curvature") > -1:
                    foundCurv = True
                    log.debug("Getting curvature values...")

                    self.curv = zeros([self.numPoints], dtype=float)
                    curvCounter = 0
                    self.file.readline()

                    while curvCounter < self.numPoints:
                        curvLine = array(self.file.readline().split(), dtype=float)
                        for i in curvLine:
                            self.curv[curvCounter] = i
                            curvCounter += 1

                    self.curvMean = self.curv.mean()
                    self.curvSD = self.curv.std()
                    log.debug("Got " + str(self.curv.shape[0]) + " curvature values")
                    return self.curvMean, self.curvSD

            if foundCurv == False:
                log.debug("No curvature values found!")

    def getEdgePoints(self, sd):

        # gets datapoints with local curvature greater than sd standard
        # deviations away from the mean

        try:
            self.points
        except AttributeError:
            self.getPoints()
            self.getCurvature()
            self.getEdgePoints(sd)
        else:
            self.edgePoints = []
            limit = sd * self.curvSD
            counter = 0

            for i in range(0, self.numPoints):
                if abs(self.curv[i] - self.curvMean) > limit:
                    self.edgePoints.append(self.points[i])
                    counter += 1

            log.debug(str(counter) + " edge points found at " + str(sd) + " SD")
            return self.edgePoints


def renderPolyData(data):
    mapper = vtk.vtkPolyDataMapper()
    if vtk.VTK_MAJOR_VERSION < 6:
        mapper.SetInput(data)
    else:
        mapper.SetInputDataObject(data)
    Actor = vtk.vtkActor()
    Actor.SetMapper(mapper)
    Actor.GetProperty().SetColor(0.5, 0.5, 0.5)

    # Create the RenderWindow
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    renWin.SetSize(400, 400)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    ren.AddActor(Actor)

    # set the properties of the renderers
    ren.SetBackground(1, 1, 1)
    ren.SetViewport(0.0, 0.0, 1.0, 1.0)
    ren.GetActiveCamera().SetPosition(1, -1, 0)
    ren.ResetCamera()

    # Render the image and start interaction.
    iren.Initialize()
    renWin.Render()
    iren.Start()


def array2vtkImage(arrayImage, dtype, flipDim=False, retImporter=False, extent=None, pipeline=False):
    # import array image into vtk
    imageImporter = vtk.vtkImageImport()

    # just to be sure
    # arr = array(arrayImage, dtype=dtype)
    arr = array(arrayImage, dtype=dtype, order='C')

    if dtype == int16:
        # log.debug('setting data scalar to int16')
        imageImporter.SetDataScalarTypeToShort()
    elif dtype == uint8:
        # log.debug('setting data scalar to uint8')
        imageImporter.SetDataScalarTypeToUnsignedChar()
    else:
        raise ValueError('Unsupported datatype {}'.format(dtype))

    imageImporter.SetNumberOfScalarComponents(1)
    # set imported image size
    s = arrayImage.shape
    if extent is None:
        if flipDim:
            extent = [0, s[2] - 1, 0, s[1] - 1, 0, s[0] - 1]
        else:
            extent = [0, s[0] - 1, 0, s[1] - 1, 0, s[2] - 1]

    imageImporter.SetDataExtent(extent)
    imageImporter.SetWholeExtent(extent)

    if vtk.VTK_MAJOR_VERSION >= 6:
        imageImporter.CopyImportVoidPointer(arr, arr.nbytes)
        # imageImporter.SetImportVoidPointer(arr, 1)
    else:
        imageString = arr.tostring()
        imageImporter.CopyImportVoidPointer(imageString, len(imageString))

    imageImporter.Update()

    # in VTK 6+, not returning the importer results in erroneous image array values
    # in the vtkImage
    if retImporter or pipeline:
        return imageImporter
    else:
        if vtk.VTK_MAJOR_VERSION >= 6:
            warnings.warn(
                'You should return the importer (retImporter=True) in VTK6 and above. Otherwise, values in the vtkimage will likely be garbage.',
                UserWarning
            )
        return imageImporter.GetOutput()


def vtkImage2Array(vtkImage, dtype, flipDim=False):
    exporter = vtk.vtkImageExport()
    if vtk.VTK_MAJOR_VERSION < 6:
        exporter.SetInput(vtkImage)
    else:
        exporter.SetInputDataObject(vtkImage)
    s = array(exporter.GetDataDimensions())
    if flipDim:
        s = s[::-1]
        I = zeros(s, dtype=dtype)
        exporter.Export(I)
        return I.transpose((2, 1, 0))
    else:
        I = zeros(s, dtype=dtype)
        exporter.Export(I)
        return I


def tri2Polydata(V, T, normals=True, featureangle=60.0):
    points = vtk.vtkPoints()
    triangles = vtk.vtkCellArray()

    for v in V:
        points.InsertNextPoint(v)

    for t in T:
        triangle = vtk.vtkTriangle()
        triangle.GetPointIds().SetId(0, t[0]);
        triangle.GetPointIds().SetId(1, t[1]);
        triangle.GetPointIds().SetId(2, t[2]);
        triangles.InsertNextCell(triangle)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(triangles)

    if normals:
        normalFilter = vtk.vtkPolyDataNormals()
        if vtk.VTK_MAJOR_VERSION < 6:
            normalFilter.SetInput(polydata)
        else:
            normalFilter.SetInputDataObject(polydata)

        normalFilter.SetComputePointNormals(True)
        normalFilter.SetComputeCellNormals(True)
        if featureangle is not None:
            normalFilter.SetFeatureAngle(featureangle)
        else:
            normalFilter.SetSplitting(False)

        normalFilter.Update()
        polydata = normalFilter.GetOutput()

    return polydata


class NoPolyDataError(Exception):
    pass


def polyData2Tri(p, pipeline=False):
    if pipeline:
        p = p.GetOutput()

    if p.GetNumberOfPoints() == 0:
        # raise ValueError('no points in polydata')
        raise NoPolyDataError('no points in polydata')

    # get vertices
    V = array([p.GetPoint(i) for i in range(p.GetNumberOfPoints())])

    # get triangles
    T = []
    for i in range(p.GetNumberOfCells()):
        ids = p.GetCell(i).GetPointIds()
        T.append((ids.GetId(0), ids.GetId(1), ids.GetId(2)))

    T = array(T, dtype=int)

    # curvature

    # normals
    polydataNormals = p.GetPointData().GetNormals()
    if polydataNormals != None:
        s = polydataNormals.GetDataSize()
        N = zeros(s, dtype=float)
        for i in range(s):
            N[i] = polydataNormals.GetValue(i)

        N = N.reshape((int(s / 3), 3))
    else:
        N = None

    return V, T, N


def polygons2Polydata(vertices, faces, vcolours=None, fcolours=None, vnormals=None):
    """
    Create a vtkPolyData instance from a set of vertices and
    faces.

    Inputs:
    vertices: (nx3) array of vertex coordinates
    faces: list of lists of vertex indices for each face
    vcolour : list of 3-tuple, vertex colours. Assigned to a vtkPointData
        array named "colours".
    fcolour : list of 3-tuple, face colours [Not implemented]
    normals: vertex normals

    Returns:
    P: vtkPolyData instance
    """
    # define points
    points = vtk.vtkPoints()

    if sys.version_info.major == 3:
        points.SetData(numpy_support.numpy_to_vtk(
            ascontiguousarray(array(vertices))
        ))
    else:
        for x, y, z in vertices:
            points.InsertNextPoint(x, y, z)

    # create polygons
    polygons = vtk.vtkCellArray()
    for f in faces:
        polygon = vtk.vtkPolygon()
        polygon.GetPointIds().SetNumberOfIds(len(f))
        for fi, gfi in enumerate(f):
            polygon.GetPointIds().SetId(fi, gfi)
        polygons.InsertNextCell(polygon)

    # create polydata
    P = vtk.vtkPolyData()
    P.SetPoints(points)
    P.SetPolys(polygons)

    # assign vertex colours
    if vcolours is not None:
        colors = vtk.vtkUnsignedCharArray()
        colors.SetNumberOfComponents(3)
        colors.SetName("colours")
        for c in vcolours:
            if vtk.VTK_MAJOR_VERSION < 6:
                colors.InsertNextTupleValue(c)
            else:
                colors.InsertNextTuple3(*c)

        P.GetPointData().SetScalars(colors)
        P.Modified()

    # assign normals to points
    if vnormals is not None:
        if vnormals.shape != vertices.shape:
            raise ValueError('vnormals must have same shape as vertices')

        vtk_vnormals = numpy_support.numpy_to_vtk(
            ascontiguousarray(array(vnormals))
        )
        P.GetPointData().SetNormals(vtk_vnormals)
        P.Modified()
        log.debug('normals set')

    return P


def polygons2Tri(vertices, faces, clean=False, normals=False):
    """
    Uses vtkTriangleFilter to convert a set of polygons
    to triangles. 

    Inputs:
    vertices: (nx3) array of vertex coordinates
    faces: list of lists of vertex indices for each face
    clean: run vtkCleanPolyData
    normals: run vtkPolyDataNormals

    Returns:
    V: (mx3) array of triangulated vertex coordinates
    T: (px3) array of vertex indices for each triangle
    N: (px3) face normals (optional)
    """
    polydata = polygons2Polydata(vertices, faces)

    # triangle filter
    tri_filter = vtk.vtkTriangleFilter()
    if vtk.VTK_MAJOR_VERSION < 6:
        tri_filter.SetInput(polydata)
    else:
        tri_filter.SetInputDataObject(polydata)
    tri_filter.Update()
    getPreviousOutput = tri_filter.GetOutput

    # clean mesh
    if clean:
        log.debug("cleaning...")
        cleaner = vtk.vtkCleanPolyData()
        if vtk.VTK_MAJOR_VERSION < 6:
            cleaner.SetInput(getPreviousOutput())
        else:
            cleaner.SetInputDataObject(getPreviousOutput())
        cleaner.SetConvertLinesToPoints(1)
        cleaner.SetConvertStripsToPolys(1)
        cleaner.SetConvertPolysToLines(1)
        cleaner.SetPointMerging(True)
        cleaner.SetTolerance(0.0)
        cleaner.Update()
        getPreviousOutput = cleaner.GetOutput

    # filter normals
    if normals:
        log.debug("filtering normals...")
        normal = vtk.vtkPolyDataNormals()
        if vtk.VTK_MAJOR_VERSION < 6:
            normal.SetInput(getPreviousOutput())
        else:
            normal.SetInputDataObject(getPreviousOutput())
        normal.SetAutoOrientNormals(1)
        normal.SetComputePointNormals(1)
        normal.SetConsistency(1)
        normal.Update()
        getPreviousOutput = normal.GetOutput

    # get triangulated vertices and faces
    return polyData2Tri(getPreviousOutput())


class polydataFromImageParams(object):
    def __init__(self):
        self.smoothImage = 1
        self.imgSmthSD = 2.0
        self.imgSmthRadius = 1.5
        self.isoValue = 200.0
        self.smoothIt = 100
        self.smoothFeatureEdge = 0
        self.deciRatio = 0.5  # higher the ratio, more decimation
        self.deciPerserveTopology = 0
        self.clean = True
        self.cleanPointMerging = 1
        self.cleanTolerance = 0.0
        self.filterNormal = 1
        self.calcCurvature = 1

    def save(self, filename):
        f = open(filename + '.polyparams', 'w')
        pickle.dump(self, f)
        f.close()


class DummyFilter(object):

    def __init__(self, output_data_object):
        self.data_object = output_data_object

    def getOutputDataObject(self):
        return self.data_object

    def getOutput(self):
        return self.data_object


def polydataFromImage(vtkImage, params, disp=0, pipeline=False):
    if pipeline:
        previousFilter = vtkImage
    else:
        previousFilter = DummyFilter(vtkImage)

    # testing - gaussian smoothing to binary image
    if params.smoothImage:
        log.debug('smoothing image...')
        imageSmoother = vtk.vtkImageGaussianSmooth()
        imageSmoother.SetStandardDeviation(params.imgSmthSD)
        imageSmoother.SetRadiusFactor(params.imgSmthRadius)

        if vtk.VTK_MAJOR_VERSION < 6:
            imageSmoother.SetInput(previousFilter.getOutput())
        else:
            if pipeline:
                imageSmoother.SetInputConnection(previousFilter.GetOutputPort())
            else:
                imageSmoother.SetInputDataObject(previousFilter.getOutputDataObject())

        previousFilter = imageSmoother

    # triangulate image to create mesh  
    log.debug("extracting contour...")
    # contourExtractor = vtk.vtkContourFilter()
    # contourExtractor.GenerateTrianglesOn()
    contourExtractor = vtk.vtkMarchingCubes()  # causes artefact faces in the corner of the volume
    # contourExtractor = vtk.vtkImageMarchingCubes()  # causes artefact faces in the corner of the volume
    contourExtractor.ComputeNormalsOn()
    contourExtractor.ComputeScalarsOn()
    contourExtractor.SetValue(0, params.isoValue)

    if vtk.VTK_MAJOR_VERSION < 6:
        contourExtractor.SetInput(previousFilter.getOutput())
    else:
        if pipeline:
            contourExtractor.SetInputConnection(previousFilter.GetOutputPort())
        else:
            contourExtractor.SetInputDataObject(previousFilter.getOutputDataObject())
            contourExtractor.Update()

    previousFilter = contourExtractor

    # triangle filter
    log.debug("filtering triangles...")
    triFilter = vtk.vtkTriangleFilter()

    if vtk.VTK_MAJOR_VERSION < 6:
        triFilter.SetInput(previousFilter.getOutput())
    else:
        if pipeline:
            triFilter.SetInputConnection(previousFilter.GetOutputPort())
        else:
            triFilter.SetInputDataObject(previousFilter.getOutputDataObject())
            triFilter.Update()

    previousFilter = triFilter

    # smooth polydata
    if params.smoothIt:
        log.debug("smoothing...")
        smoother = vtk.vtkSmoothPolyDataFilter()
        smoother.SetNumberOfIterations(params.smoothIt)
        smoother.SetFeatureEdgeSmoothing(params.smoothFeatureEdge)

        if vtk.VTK_MAJOR_VERSION < 6:
            smoother.SetInput(previousFilter.getOutput())
        else:
            if pipeline:
                smoother.SetInputConnection(previousFilter.GetOutputPort())
            else:
                smoother.SetInputDataObject(previousFilter.getOutputDataObject())
                smoother.Update()

        previousFilter = smoother

    # decimate polydata
    if params.deciRatio:
        # print "decimating..."
        # decimator = vtk.vtkDecimatePro()
        # decimator.SetInput( getPreviousOutput() )
        # decimator.SetTargetReduction( params.deciRatio )
        # decimator.SetPreserveTopology( params.deciPerserveTopology )
        # decimator.SplittingOn()
        # decimator.Update()
        # getPreviousOutput = decimator.GetOutput
        # if disp:
        #   RenderPolyData( decimator.GetOutput() )

        log.debug("decimating using quadric...")
        decimator = vtk.vtkQuadricDecimation()
        decimator.SetTargetReduction(params.deciRatio)

        if vtk.VTK_MAJOR_VERSION < 6:
            decimator.SetInput(previousFilter.getOutput())
        else:
            if pipeline:
                decimator.SetInputConnection(previousFilter.GetOutputPort())
            else:
                decimator.SetInputDataObject(previousFilter.getOutputDataObject())
                decimator.Update()

        previousFilter = decimator

    # clean mesh
    if params.clean:
        log.debug("cleaning...")
        cleaner = vtk.vtkCleanPolyData()
        cleaner.SetConvertLinesToPoints(1)
        cleaner.SetConvertStripsToPolys(1)
        cleaner.SetConvertPolysToLines(1)
        cleaner.SetPointMerging(params.cleanPointMerging)
        cleaner.SetTolerance(params.cleanTolerance)

        if vtk.VTK_MAJOR_VERSION < 6:
            cleaner.SetInput(previousFilter.getOutput())
        else:
            if pipeline:
                cleaner.SetInputConnection(previousFilter.GetOutputPort())
            else:
                cleaner.SetInputDataObject(previousFilter.getOutputDataObject())
                cleaner.Update()

        previousFilter = cleaner

    # filter normals
    if params.filterNormal:
        log.debug("filtering normals...")
        normal = vtk.vtkPolyDataNormals()
        normal.SetAutoOrientNormals(1)
        normal.SetComputePointNormals(1)
        normal.SetConsistency(1)

        if vtk.VTK_MAJOR_VERSION < 6:
            normal.SetInput(previousFilter.getOutput())
        else:
            if pipeline:
                normal.SetInputConnection(previousFilter.GetOutputPort())
            else:
                normal.SetInputDataObject(previousFilter.getOutputDataObject())
                normal.Update()

        previousFilter = normal

    if params.calcCurvature:
        log.debug("calculating curvature...")
        curvature = vtk.vtkCurvatures()
        curvature.SetCurvatureTypeToMean()

        if vtk.VTK_MAJOR_VERSION < 6:
            curvature.SetInput(previousFilter.getOutput())
        else:
            if pipeline:
                curvature.SetInputConnection(previousFilter.GetOutputPort())
            else:
                curvature.SetInputDataObject(previousFilter.getOutputDataObject())
                curvature.Update()

        previousFilter = curvature

    previousFilter.Update()

    if pipeline:
        return previousFilter
    else:
        return previousFilter.getOutputDataObject()


def triSurface2BinaryMask(v, t, imageShape, outputOrigin=None, outputSpacing=None, extent=None):
    """Create a binary image mask from a triangulated surface.

    Inputs
    ------
    v : an nx3 array of vertex coordinates.
    t : an mx3 array of the vertex indices of triangular faces
    imageShape : a 3-tuple of the output binary image array shape
    outputOrigin : 3D coordinates of the origin of the output image array
    outputSpacing : Voxel spacing of the output image array

    Returns
    -------
    maskImageArray : binary image array
    gfPoly : vtkPolyData instance of the triangulated surface
    """

    if outputOrigin is None:
        outputOrigin = [0.0, 0.0, 0.0]
    if outputSpacing is None:
        outputSpacing = [1.0, 1.0, 1.0]

    imgDtype = uint8

    # make into vtkPolydata
    gfPoly = tri2Polydata(v, t)

    # create mask vtkImage
    maskImageArray = ones(imageShape, dtype=imgDtype)
    maskVTKImageImporter = array2vtkImage(
        maskImageArray, imgDtype, flipDim=False, extent=extent,
        retImporter=True
    )
    maskVTKImage = maskVTKImageImporter.GetOutput()

    # create stencil from polydata
    stencilMaker = vtk.vtkPolyDataToImageStencil()
    if vtk.VTK_MAJOR_VERSION < 6:
        stencilMaker.SetInput(gfPoly)
    else:
        stencilMaker.SetInputDataObject(gfPoly)
    stencilMaker.SetOutputOrigin(outputOrigin)
    stencilMaker.SetOutputSpacing(outputSpacing)
    stencilMaker.SetOutputWholeExtent(maskVTKImage.GetExtent())
    stencilMaker.Update()  # needed in VTK 6

    stencil = vtk.vtkImageStencil()
    if vtk.VTK_MAJOR_VERSION < 6:
        stencil.SetInput(maskVTKImage)
        stencil.SetStencil(stencilMaker.GetOutput())
    else:
        stencil.SetInputDataObject(maskVTKImage)
        stencil.SetStencilData(stencilMaker.GetOutput())
    stencil.SetBackgroundValue(0)
    stencil.ReverseStencilOff()
    stencil.Update()

    maskImageArray2 = vtkImage2Array(stencil.GetOutput(), imgDtype, flipDim=True)
    return maskImageArray2, gfPoly
    # return maskImageArray2, gfPoly, maskVTKImage


def _makeImageSpaceGF(scan, GF, negSpacing=False, zShift=True):
    """
    Transform a fieldwork geometric field from physical coords to image voxel indices
    """
    newGF = copy.deepcopy(GF)
    p = GF.get_all_point_positions()
    pImg = scan.coord2Index(p, negSpacing=negSpacing, zShift=zShift, roundInt=False)
    newGF.set_field_parameters(pImg.T[:, :, newaxis])

    return newGF


def gf2BinaryMask(gf, scan, xiD=None, negSpacing=False, zShift=True,
                  outputOrigin=(0, 0, 0), outputSpacing=(1, 1, 1)):
    """Create a binary image mask from a GeometricField instance.

    Inputs
    ------
    gf : GeometricField instance
    scan : a Scan instance of the image volume the mask should be created for
    xiD : discretisation of gf, default is [10,10]
    zShift : shift model in the Z axis by image volume height.
    negSpacing : mirror the model in the Z axis
    outputOrigin : 3D coordinates of the origin of the output image array
    outputSpacing : Voxel spacing of the output image array

    Returns
    -------
    maskImageArray : binary image array
    gfPoly : vtkPolyData instance of the triangulated surface
    """
    if xiD is None:
        xiD = [10, 10]
    imgDtype = int16
    gfImage = _makeImageSpaceGF(scan, gf, negSpacing, zShift)
    vertices, triangles = gfImage.triangulate(xiD, merge=True)
    return triSurface2BinaryMask(vertices, triangles, scan.I.shape, outputOrigin, outputSpacing)


def simplemesh2BinaryMask(sm, scan, zShift=True, negSpacing=False,
                          outputOrigin=(0, 0, 0), outputSpacing=(1, 1, 1)):
    """Create a binary image mask from a SimpleMesh instance.

    Inputs
    ------
    sm : SimpleMesh instance
    scan : a Scan instance of the image volume the mask should be created for
    zShift : shift model in the Z axis by image volume height.
    negSpacing : mirror the model in the Z axis
    outputOrigin : 3D coordinates of the origin of the output image array
    outputSpacing : Voxel spacing of the output image array

    Returns
    -------
    maskImageArray : binary image array
    gfPoly : vtkPolyData instance of the triangulated surface
    """
    vImage = scan.coord2Index(sm.v, zShift=zShift, negSpacing=negSpacing, roundInt=False)
    return triSurface2BinaryMask(vImage, sm.f, scan.I.shape, outputOrigin, outputSpacing)


def image2Simplemesh(imageArray, index2Coord, isoValue, deciRatio=None, smoothIt=200, zShift=True):
    """Convert an image array into a SimpleMesh surface.

    imageArray must be uint8.

    inputs
    ------
    imageArray : binary image array
    index2Coord : a function that given an array of image indices, returns the scanner coordinates
    isoValue : the image value at which an isosurface will be created to generate the surface mesh
    deciRatio : amount of decimation to apply to the surface, for vtkQuadricDecimation
    smoothIt : number of smoothing iterations
    zShift : shift surface in the Z axis by image volume height.

    """

    IMGDTYPE = uint8
    if imageArray.dtype != IMGDTYPE:
        raise ValueError('imageArray must be {}.'.format(IMGDTYPE))

    if vtk.VTK_MAJOR_VERSION < 6:
        vtkImage = array2vtkImage(imageArray, IMGDTYPE, flipDim=True, retImporter=False)
    else:
        vtkImage = array2vtkImage(imageArray, IMGDTYPE, flipDim=True, retImporter=True, pipeline=True)
        # vtkImage = array2vtkImage(imageArray, IMGDTYPE, flipDim=True, retImporter=True).GetOutput()

    params = polydataFromImageParams()
    params.smoothImage = False
    params.imgSmthRadius = 1.0
    params.imgSmthSD = 1.0
    params.isoValue = isoValue
    params.smoothIt = smoothIt
    params.smoothFeatureEdge = 0
    params.deciRatio = deciRatio
    params.deciPerserveTopology = 1
    params.clean = 1
    params.cleanPointMerging = 1
    params.cleanTolerance = 1e-9
    params.filterNormal = 1
    params.calcCurvature = 0
    if vtk.VTK_MAJOR_VERSION < 6:
        polydata = polydataFromImage(vtkImage, params, pipeline=False)
        V, T, N = polyData2Tri(polydata, pipeline=False)
    else:
        polydata = polydataFromImage(vtkImage, params, pipeline=True)
        V, T, N = polyData2Tri(polydata, pipeline=True)

    V = V[:, ::-1] + [0.0, 0.0, 1.0]
    SMImg = simplemesh.SimpleMesh(v=V, f=T)
    SM = simplemesh.SimpleMesh(v=index2Coord(V, zShift=zShift), f=T)
    SM.data = {'vertexnormal': N}
    log.debug('image-to-mesh done')
    log.debug('vertices: {}'.format(SM.v.shape[0]))
    log.debug('faces: {}'.format(SM.f.shape[0]))
    return SM, SMImg, polydata


def smoothMeshVTK(mesh, it, smoothboundary=False, smoothfeatures=False, relaxfactor=1.0, usewsinc=True):
    """
    Apply smoothing to a SimpleMesh instance using VTK's SmoothPolyDataFilter
    or WindowedSincPolyDataFilter.

    inputs
    ======
    mesh : SimpleMesh instance
        Mesh to be smoothed
    it : int
        Smoothing iterations to apply
    smoothboundary : bool
        Whether to smooth boundary vertices
    smoothfeatures : bool
        Whether to smooth mesh feature edges differently.
    relaxfactor : float
        Relaxation factor for vtkSmoothPolyDataFilter
    usewsinc : bool
        Use vtkWindowedSincPolyDataFilter instead of 
        vtkSmoothPolyDataFilter if True

    returns
    =======
    mesh_smooth : SimpleMesh instance
        A smoothed copy of the input mesh

    """
    poly = tri2Polydata(mesh.v, mesh.f, featureangle=None)

    if usewsinc:
        smoother = vtk.vtkWindowedSincPolyDataFilter()
    else:
        smoother = vtk.vtkSmoothPolyDataFilter()

    if vtk.VTK_MAJOR_VERSION < 6:
        smoother.SetInput(poly)
    else:
        smoother.SetInputDataObject(poly)
    smoother.SetNumberOfIterations(it)
    if smoothfeatures:
        smoother.FeatureEdgeSmoothingOn()
    else:
        smoother.FeatureEdgeSmoothingOff()
    if smoothboundary:
        smoother.BoundarySmoothingOn()
    else:
        smoother.BoundarySmoothingOff()
    if not usewsinc:
        smoother.SetRelaxationFactor(relaxfactor)
    smoother.Update()
    poly_smooth = smoother.GetOutput()

    v, t, n = polyData2Tri(poly_smooth)
    mesh_smooth = simplemesh.SimpleMesh(v=v, f=t)
    if mesh.has1Ring:
        mesh_smooth.vertices1Ring = dict(mesh.vertices1Ring)
        mesh_smooth.faces1Ring = dict(mesh.faces1Ring)
        mesh_smooth.has1Ring = True
    if mesh.hasNeighbourhoods:
        mesh_smooth.neighbourFaces = list(mesh.neighbourFaces)
        mesh_smooth.neighbourVertices = list(mesh.neighbourVertices)
        mesh_smooth.hasNeighbourhoods = True
    return mesh_smooth


def optimiseMesh(sm, deciratio, clean=False):
    """
    Optimise a triangle mesh by removing degenerate elements
    and decimation.

    Inputs
    ======
    sm : SimpleMesh instance
        Mesh to be optimised
    deciratio : float
        Decimation target, fraction of original number of faces to remove.
        Closer to 1 - more faces removed.

    Returns
    =======
    newSM : SimpleMesh instance
        Optimised mesh
    """

    poly = tri2Polydata(
        sm.v, sm.f,
        normals=True,
        featureangle=None
    )

    if clean:
        log.debug("cleaning...")
        cleaner = vtk.vtkCleanPolyData()
        if vtk.VTK_MAJOR_VERSION < 6:
            cleaner.SetInput(poly)
        else:
            cleaner.SetInputDataObject(poly)
        cleaner.SetConvertLinesToPoints(1)
        cleaner.SetConvertStripsToPolys(1)
        cleaner.SetConvertPolysToLines(1)
        cleaner.SetPointMerging(True)
        cleaner.SetTolerance(1e-9)
        cleaner.Update()
        getPreviousOutput = cleaner.GetOutput

    # decimate polydata
    log.debug("decimating using quadric...")
    decimator = vtk.vtkQuadricDecimation()
    if vtk.VTK_MAJOR_VERSION < 6:
        if not clean:
            decimator.SetInput(poly)
        else:
            decimator.SetInput(getPreviousOutput())
    else:
        if not clean:
            decimator.SetInputDataObject(poly)
        else:
            decimator.SetInputDataObject(getPreviousOutput())
    decimator.SetTargetReduction(deciratio)
    # decimator.SetPreserveTopology(True)
    # decimator.SplittingOn()
    decimator.Update()
    getPreviousOutput = decimator.GetOutput

    # convert back to sm
    v, f, N = polyData2Tri(getPreviousOutput())
    newSM = simplemesh.SimpleMesh(v=v, f=f)
    return newSM


# ====================================================#
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


def renderVtkImageVolume(vtkImage, cRange=[0, 255], oRange=[0, 255]):
    bgColour = [0.0, 0.0, 0.0]
    renderWindowSize = 800

    def _exitCheck(obj, event):
        if obj.GetEventPending() != 0:
            obj.SetAbortRender(1)

    # Volume mapper 
    volumeMapper = vtk.vtkVolumeRayCastMapper()
    if vtk.VTK_MAJOR_VERSION < 6:
        volumeMapper.SetInput(vtkImage)
    else:
        volumeMapper.SetInputDataObject(vtkImage)
    compositeFunc = vtk.vtkVolumeRayCastCompositeFunction()
    volumeMapper.SetVolumeRayCastFunction(compositeFunc)

    # Colour transfer functions
    colorFunc = vtk.vtkColorTransferFunction()
    colorFunc.AddRGBPoint(cRange[0], 0.0, 0.0, 0.0)
    colorFunc.AddRGBPoint(cRange[1], 1.0, 1.0, 1.0)

    # Opacity transfer functions
    opacityFunc = vtk.vtkPiecewiseFunction()
    opacityFunc.AddPoint(oRange[0], 0.0)
    opacityFunc.AddPoint(oRange[1], 0.1)

    # Volume properties
    volumeProperties = vtk.vtkVolumeProperty()
    volumeProperties.SetColor(colorFunc)
    volumeProperties.SetScalarOpacity(opacityFunc)

    # VTK volume
    volume = vtk.vtkVolume()
    volume.SetMapper(volumeMapper)
    volume.SetProperty(volumeProperties)

    # render axes
    axes = vtk.vtkAxesActor()
    axes.SetAxisLabels(0)
    axes.SetTotalLength(50, 100, 150)
    axes.SetConeRadius(0.1)

    # render bounding box
    outline = vtk.vtkOutlineFilter()
    if vtk.VTK_MAJOR_VERSION < 6:
        outline.SetInput(vtkImage)
    else:
        outline.SetInputDataObject(vtkImage)
    outlineMapper = vtk.vtkPolyDataMapper()
    if vtk.VTK_MAJOR_VERSION < 6:
        outlineMapper.SetInput(outline.GetOutput())
    else:
        outlineMapper.SetInputDataObject(outline.GetOutput())
    outlineActor = vtk.vtkActor()
    outlineActor.SetMapper(outlineMapper)
    outlineActor.GetProperty().SetColor(0.0, 0.0, 1.0)

    # renderer
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(bgColour[0], bgColour[1], bgColour[2])
    renderer.AddActor(axes)
    renderer.AddActor(outlineActor)
    renderer.AddVolume(volume)

    # render window
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindow.SetSize(renderWindowSize, renderWindowSize)

    # render window interactor
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renderWindow)

    renderWindow.AddObserver("AbortCheckEvent", _exitCheck)

    # render and start interaction.
    interactor.Initialize()
    renderWindow.Render()
    interactor.Start()


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

        self.imageImporter = array2vtkImage(self.image, uint8, retImporter=True)

        # imageImporter = vtk.vtkImageImport()
        # imageImporter.SetDataScalarTypeToShort()
        # imageString = arrayImage.astype(dtype).tostring()
        # imageImporter.CopyImportVoidPointer( imageString, len( imageString ) )
        # imageImporter.SetNumberOfScalarComponents(1)
        # # set imported image size
        # s = arrayImage.shape
        # # imageImporter.SetWholeExtent(0, s[2]-1, 0, s[1]-1, 0, s[0]-1)
        # imageImporter.SetWholeExtent(0, s[0]-1, 0, s[1]-1, 0, s[2]-1)
        # imageImporter.SetDataExtentToWholeExtent()    

        if 0:
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
        if vtk.VTK_MAJOR_VERSION < 6:
            CoMSphereMapper.SetInput(CoMSphere.GetOutput())
        else:
            CoMSphereMapper.SetInputDataObject(CoMSphere.GetOutput())
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
        if vtk.VTK_MAJOR_VERSION < 6:
            nodeMapper.SetInput(node.GetOutput())
        else:
            nodeMapper.SetInputDataObject(node.GetOutput())

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
        if vtk.VTK_MAJOR_VERSION < 6:
            lineMapper.SetInput(line.GetOutput())
        else:
            lineMapper.SetInputDataObject(line.GetOutput())
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
        if vtk.VTK_MAJOR_VERSION < 6:
            volumeMapper.SetInput(self.imageImporter.GetOutput())
        else:
            volumeMapper.SetInputDataObject(self.imageImporter.GetOutput())
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

    def renderPoly(self, poly):
        # render polydata

        # self.polydata.append(poly)

        mapper = vtk.vtkPolyDataMapper()
        if vtk.VTK_MAJOR_VERSION < 6:
            mapper.SetInput(poly)
        else:
            mapper.SetInputDataObject(poly)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(0.9, 0.9, 0.7)

        # ~ self.actorList.append( actor )

        self._render(actorList=[actor])

    def renderContour(self, contourValueList):
        # render polydata contour surfaces at iso values defined in
        # list contourValueList

        contourExtractor = vtk.vtkContourFilter()
        if vtk.VTK_MAJOR_VERSION < 6:
            contourExtractor.SetInput(self.imageImporter.GetOutput())
        else:
            contourExtractor.SetInputDataObject(self.imageImporter.GetOutput())
        # set contour values
        for i in range(0, len(contourValueList)):
            contourExtractor.SetValue(i, contourValueList[i])

        contourExtractor.Update()

        mapper = vtk.vtkPolyDataMapper()
        if vtk.VTK_MAJOR_VERSION < 6:
            mapper.SetInput(contourExtractor.GetOutput())
        else:
            mapper.SetInputDataObject(contourExtractor.GetOutput())
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
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
        if vtk.VTK_MAJOR_VERSION < 6:
            outline.SetInput(self.imageImporter.GetOutput())
        else:
            outline.SetInputDataObject(self.imageImporter.GetOutput())
        outlineMapper = vtk.vtkPolyDataMapper()
        if vtk.VTK_MAJOR_VERSION < 6:
            outlineMapper.SetInput(outline.GetOutput())
        else:
            outlineMapper.SetInputDataObject(outline.GetOutput())
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
