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
from typing import Optional, Callable, Tuple, List, Any, Type

import sys
import vtk
from numpy import zeros, array, uint8, int16, ones, newaxis, ascontiguousarray, ndarray
from vtk.util import numpy_support

from gias2.image_analysis.image_tools import Scan
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

    def setFilename(self, f: str) -> None:
        self.filename = f
        self._parse_format()

    def _parse_format(self) -> None:
        self.file_prefix, self.file_ext = path.splitext(self.filename)
        self.file_ext = self.file_ext.lower()

    def _make_polydata(self) -> None:
        self._polydata = polygons2Polydata(
            self._vertices,
            self._faces,
            vcolours=self._vertex_colour,
            fcolours=self._face_colour,
            vnormals=self._vertex_normals,
        )

    def _make_render_window(self) -> None:
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

    def write(self, filename: Optional[str] = None, ascenc: bool = True) -> None:
        if filename is not None:
            self.filename = filename

        file_prefix, file_ext = path.splitext(self.filename)
        file_ext = file_ext.lower()
        if file_ext == '.obj':
            self.writeOBJ()
        elif file_ext == '.wrl':
            self.writeVRML()
        elif file_ext == '.stl':
            self.writeSTL(ascenc=ascenc)
        elif file_ext == '.ply':
            self.writePLY(ascenc=ascenc)
        elif file_ext == '.vtp':
            self.writeVTP(ascenc=ascenc)
        else:
            raise ValueError('unknown file extension')

    def writeOBJ(self, filename: Optional[str] = None) -> None:
        if filename is not None:
            self.filename = filename
        if self._render_window is None:
            self._make_render_window()

        w = vtk.vtkOBJExporter()
        w.SetRenderWindow(self._render_window)
        w.SetFilePrefix(path.splitext(self.filename)[0])
        w.Write()

    def writePLY(self, filename: Optional[str] = None, ascenc: bool = True) -> None:
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

    def writeSTL(self, filename: Optional[str] = None, ascenc: bool = True) -> None:
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

    def writeVRML(self, filename: Optional[str] = None) -> None:
        if filename is not None:
            self.filename = filename
        if self._render_window is None:
            self._make_render_window()

        w = vtk.vtkVRMLExporter()
        w.SetRenderWindow(self._render_window)
        w.SetFileName(self.filename)
        w.Write()

    def writeVTP(self, filename: Optional[str] = None, ascenc: bool = True, xmlenc: bool = True) -> None:
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
        self.filename: str = kwargs.get('filename')
        self.verbose: bool = kwargs.get('verbose', True)
        self._points = None
        self._triangles = None
        self._vertexNormals = None
        self._nPoints = None
        self._nFaces = None
        self._nVertexNormals = None
        self._dimensions = None
        self.polydata = None

    def setFilename(self, filename: str) -> None:
        self.filename = filename

    def read(self, filename: Optional[str] = None) -> None:
        if filename is not None:
            self.filename = filename

        file_prefix, file_ext = path.splitext(self.filename)
        file_ext = file_ext.lower()
        if file_ext == '.obj':
            self.readOBJ()
        elif file_ext == '.wrl':
            self.readVRML()
        elif file_ext == '.stl':
            self.readSTL()
        elif file_ext == '.ply':
            self.readPLY()
        elif file_ext == '.vtp':
            self.readVTP()
        else:
            raise ValueError('unknown file extension in {}'.format(self.filename))

    def readVRML(self, filename: Optional[str] = None, actor: int = 0) -> None:
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

    def readOBJ(self, filename: Optional[str] = None) -> None:
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

    def readPLY(self, filename: Optional[str] = None) -> None:
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

    def readSTL(self, filename: Optional[str] = None) -> None:
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

    def readVTP(self, filename: Optional[str] = None) -> None:
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

    def _isXML(self, f: str) -> bool:
        """Check if file is an xml file
        """
        with open(f, 'r') as fp:
            l = fp.readline()

        if l[0] == '<':
            return True
        else:
            return False

    def _load(self) -> None:
        self._loadPoints()
        self._loadTriangles()
        self._loadVertexNormals()

    def _loadPoints(self) -> None:
        points = self.polydata.GetPoints().GetData()
        self._dimensions = points.GetNumberOfComponents()
        self._nPoints = points.GetNumberOfTuples()

        if self.verbose:
            log.debug('loading %(np)i points in %(d)i dimensions' % {'np': self._nPoints, 'd': self._dimensions})

        self._points = numpy_support.vtk_to_numpy(points)

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

    def _loadTriangles(self) -> None:
        poly_data = self.polydata.GetPolys().GetData()
        face_indices = numpy_support.vtk_to_numpy(poly_data)

        # assumes that faces are triangular
        face_indices = array(face_indices).reshape((-1, 4))
        self._nFaces = face_indices.shape[0]
        self._triangles = face_indices[:, 1:].copy()

        if self.verbose:
            log.debug('loaded %s faces', self._nFaces)

    def _loadVertexNormals(self) -> None:
        pts_normals = self.polydata.GetPointData().GetNormals()
        if pts_normals is not None:
            self._vertexNormals = numpy_support.vtk_to_numpy(pts_normals)
            self._nVertexNormals = self._vertexNormals.shape[0]
        else:
            self._nVertexNormals = 0
            self._vertexNormals = None

        if self.verbose:
            log.debug('loaded %s vertex normals', self._nVertexNormals)

    def getSimplemesh(self) -> simplemesh.SimpleMesh:
        mesh = simplemesh.SimpleMesh(self._points, self._triangles)

        if self._vertexNormals is not None:
            mesh.vertexNormals = array(self._vertexNormals)
            mesh.hasVertexNormals = True

        return mesh


def savepoly(sm: simplemesh.SimpleMesh, filename: str, ascenc: bool = True) -> None:
    w = Writer(v=sm.v, f=sm.f, vn=sm.vertexNormals)
    w.write(filename, ascenc=ascenc)


def loadpoly(filename: str, verbose: bool = False) -> simplemesh.SimpleMesh:
    r = Reader(verbose=verbose)
    r.read(filename)
    return r.getSimplemesh()


def renderPolyData(data: vtk.vtkPolyData):
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


def array2vtkImage(
        array_image: ndarray,
        dtype: Type,
        flipDim: bool = False,
        retImporter: bool = False,
        extent: Optional[List[int]] = None,
        pipeline: bool = False) -> vtk.vtkImageData:
    # import array image into vtk
    image_importer = vtk.vtkImageImport()

    # just to be sure
    arr = array(array_image, dtype=dtype, order='C')

    if dtype == int16:
        # log.debug('setting data scalar to int16')
        image_importer.SetDataScalarTypeToShort()
    elif dtype == uint8:
        # log.debug('setting data scalar to uint8')
        image_importer.SetDataScalarTypeToUnsignedChar()
    else:
        raise ValueError('Unsupported datatype {}'.format(dtype))

    image_importer.SetNumberOfScalarComponents(1)
    # set imported image size
    s = array_image.shape
    if extent is None:
        if flipDim:
            extent = [0, s[2] - 1, 0, s[1] - 1, 0, s[0] - 1]
        else:
            extent = [0, s[0] - 1, 0, s[1] - 1, 0, s[2] - 1]

    image_importer.SetDataExtent(*extent)
    image_importer.SetWholeExtent(*extent)

    if vtk.VTK_MAJOR_VERSION >= 6:
        image_importer.CopyImportVoidPointer(arr, arr.nbytes)
    else:
        image_string = arr.tostring()
        image_importer.CopyImportVoidPointer(image_string, len(image_string))

    image_importer.Update()

    # in VTK 6+, not returning the importer results in erroneous image array values
    # in the vtkImage
    if retImporter or pipeline:
        return image_importer
    else:
        if vtk.VTK_MAJOR_VERSION >= 6:
            warnings.warn(
                'You should return the importer (retImporter=True) in VTK6 and above. Otherwise, values in the vtkimage will likely be garbage.',
                UserWarning
            )
        return image_importer.GetOutput()


def vtkImage2Array(vtk_image: vtk.vtkImageData, dtype: Type, flipDim: bool = False) -> ndarray:
    exporter = vtk.vtkImageExport()
    if vtk.VTK_MAJOR_VERSION < 6:
        exporter.SetInput(vtk_image)
    else:
        exporter.SetInputDataObject(vtk_image)
    s = array(exporter.GetDataDimensions())
    if flipDim:
        s = s[::-1]
        img_arr = zeros(s, dtype=dtype)
        exporter.Export(img_arr)
        return img_arr.transpose((2, 1, 0))
    else:
        img_arr = zeros(s, dtype=dtype)
        exporter.Export(img_arr)
        return img_arr


def tri2Polydata(
        vertices: ndarray,
        tris: ndarray,
        normals: bool = True,
        featureangle: Optional[float] = 60.0) -> vtk.vtkPolyData:
    points = vtk.vtkPoints()
    triangles = vtk.vtkCellArray()

    for v in vertices:
        points.InsertNextPoint(v)

    for t in tris:
        triangle = vtk.vtkTriangle()
        triangle.GetPointIds().SetId(0, t[0])
        triangle.GetPointIds().SetId(1, t[1])
        triangle.GetPointIds().SetId(2, t[2])
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


def polyData2Tri(p: vtk.vtkPolyData, pipeline: bool = False) -> Tuple[ndarray, ndarray, ndarray]:
    if pipeline:
        p = p.GetOutput()

    if p.GetNumberOfPoints() == 0:
        raise NoPolyDataError('no points in polydata')

    # get vertices
    verts = array([p.GetPoint(i) for i in range(p.GetNumberOfPoints())])

    # get triangles
    tris = []
    for i in range(p.GetNumberOfCells()):
        ids = p.GetCell(i).GetPointIds()
        tris.append((ids.GetId(0), ids.GetId(1), ids.GetId(2)))

    tris = array(tris, dtype=int)

    # curvature

    # normals
    polydata_normals = p.GetPointData().GetNormals()
    if polydata_normals is not None:
        s = polydata_normals.GetDataSize()
        normals = zeros(s, dtype=float)
        for i in range(s):
            normals[i] = polydata_normals.GetValue(i)

        normals = normals.reshape((int(s / 3), 3))
    else:
        normals = None

    return verts, tris, normals


def polygons2Polydata(
        vertices: ndarray,
        faces: List[List[int]],
        vcolours: Optional[ndarray] = None,
        fcolours: Optional[ndarray] = None,
        vnormals: Optional[ndarray] = None) -> vtk.vtkPolyData:
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
            points.InsertNextPoint((x, y, z))

    # create polygons
    polygons = vtk.vtkCellArray()
    for f in faces:
        polygon = vtk.vtkPolygon()
        polygon.GetPointIds().SetNumberOfIds(len(f))
        for fi, gfi in enumerate(f):
            polygon.GetPointIds().SetId(fi, gfi)
        polygons.InsertNextCell(polygon)

    # create polydata
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(polygons)

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

        polydata.GetPointData().SetScalars(colors)
        polydata.Modified()

    # assign normals to points
    if vnormals is not None:
        if vnormals.shape != vertices.shape:
            raise ValueError('vnormals must have same shape as vertices')

        vtk_vnormals = numpy_support.numpy_to_vtk(
            ascontiguousarray(array(vnormals))
        )
        polydata.GetPointData().SetNormals(vtk_vnormals)
        polydata.Modified()
        log.debug('normals set')

    return polydata


def polygons2Tri(
        vertices: ndarray,
        faces: List[List[int]],
        clean: bool = False,
        normals: bool = False) -> Tuple[ndarray, ndarray, ndarray]:
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


class PolydataFromImageParams(object):
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

    def save(self, filename: str) -> None:
        f = open(filename + '.polyparams', 'w')
        pickle.dump(self, f)
        f.close()


class DummyFilter(object):

    def __init__(self, output_data_object: Any):
        self.data_object = output_data_object

    def getOutputDataObject(self) -> Any:
        return self.data_object

    def getOutput(self) -> Any:
        return self.data_object


def polydataFromImage(
        vtk_image: vtk.vtkImageData,
        params: PolydataFromImageParams,
        disp: int = 0,
        pipeline: bool = False) -> vtk.vtkPolyData:
    if pipeline:
        previous_filter = vtk_image
    else:
        previous_filter = DummyFilter(vtk_image)

    # testing - gaussian smoothing to binary image
    if params.smoothImage:
        log.debug('smoothing image...')
        image_smoother = vtk.vtkImageGaussianSmooth()
        image_smoother.SetStandardDeviation(params.imgSmthSD)
        image_smoother.SetRadiusFactor(params.imgSmthRadius)

        if vtk.VTK_MAJOR_VERSION < 6:
            image_smoother.SetInput(previous_filter.getOutput())
        else:
            if pipeline:
                image_smoother.SetInputConnection(previous_filter.GetOutputPort())
            else:
                image_smoother.SetInputDataObject(previous_filter.getOutputDataObject())

        previous_filter = image_smoother

    # triangulate image to create mesh  
    log.debug("extracting contour...")
    # contour_extractor = vtk.vtkContourFilter()
    # contour_extractor.GenerateTrianglesOn()
    contour_extractor = vtk.vtkMarchingCubes()  # causes artefact faces in the corner of the volume
    # contour_extractor = vtk.vtkImageMarchingCubes()  # causes artefact faces in the corner of the volume
    contour_extractor.ComputeNormalsOn()
    contour_extractor.ComputeScalarsOn()
    contour_extractor.SetValue(0, params.isoValue)

    if vtk.VTK_MAJOR_VERSION < 6:
        contour_extractor.SetInput(previous_filter.getOutput())
    else:
        if pipeline:
            contour_extractor.SetInputConnection(previous_filter.GetOutputPort())
        else:
            contour_extractor.SetInputDataObject(previous_filter.getOutputDataObject())
            contour_extractor.Update()

    previous_filter = contour_extractor

    # triangle filter
    log.debug("filtering triangles...")
    tri_filter = vtk.vtkTriangleFilter()

    if vtk.VTK_MAJOR_VERSION < 6:
        tri_filter.SetInput(previous_filter.getOutput())
    else:
        if pipeline:
            tri_filter.SetInputConnection(previous_filter.GetOutputPort())
        else:
            tri_filter.SetInputDataObject(previous_filter.getOutputDataObject())
            tri_filter.Update()

    previous_filter = tri_filter

    # smooth polydata
    if params.smoothIt:
        log.debug("smoothing...")
        smoother = vtk.vtkSmoothPolyDataFilter()
        smoother.SetNumberOfIterations(params.smoothIt)
        smoother.SetFeatureEdgeSmoothing(params.smoothFeatureEdge)

        if vtk.VTK_MAJOR_VERSION < 6:
            smoother.SetInput(previous_filter.getOutput())
        else:
            if pipeline:
                smoother.SetInputConnection(previous_filter.GetOutputPort())
            else:
                smoother.SetInputDataObject(previous_filter.getOutputDataObject())
                smoother.Update()

        previous_filter = smoother

    # decimate polydata
    if params.deciRatio:
        log.debug("decimating using quadric...")
        decimator = vtk.vtkQuadricDecimation()
        decimator.SetTargetReduction(params.deciRatio)

        if vtk.VTK_MAJOR_VERSION < 6:
            decimator.SetInput(previous_filter.getOutput())
        else:
            if pipeline:
                decimator.SetInputConnection(previous_filter.GetOutputPort())
            else:
                decimator.SetInputDataObject(previous_filter.getOutputDataObject())
                decimator.Update()

        previous_filter = decimator

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
            cleaner.SetInput(previous_filter.getOutput())
        else:
            if pipeline:
                cleaner.SetInputConnection(previous_filter.GetOutputPort())
            else:
                cleaner.SetInputDataObject(previous_filter.getOutputDataObject())
                cleaner.Update()

        previous_filter = cleaner

    # filter normals
    if params.filterNormal:
        log.debug("filtering normals...")
        normal = vtk.vtkPolyDataNormals()
        normal.SetAutoOrientNormals(1)
        normal.SetComputePointNormals(1)
        normal.SetConsistency(1)

        if vtk.VTK_MAJOR_VERSION < 6:
            normal.SetInput(previous_filter.getOutput())
        else:
            if pipeline:
                normal.SetInputConnection(previous_filter.GetOutputPort())
            else:
                normal.SetInputDataObject(previous_filter.getOutputDataObject())
                normal.Update()

        previous_filter = normal

    if params.calcCurvature:
        log.debug("calculating curvature...")
        curvature = vtk.vtkCurvatures()
        curvature.SetCurvatureTypeToMean()

        if vtk.VTK_MAJOR_VERSION < 6:
            curvature.SetInput(previous_filter.getOutput())
        else:
            if pipeline:
                curvature.SetInputConnection(previous_filter.GetOutputPort())
            else:
                curvature.SetInputDataObject(previous_filter.getOutputDataObject())
                curvature.Update()

        previous_filter = curvature

    previous_filter.Update()

    if pipeline:
        return previous_filter
    else:
        return previous_filter.getOutputDataObject()


def triSurface2BinaryMask(
        v: ndarray,
        t: ndarray,
        image_shape: Tuple[int, int, int],
        outputOrigin: List[float] = None,
        outputSpacing: List[float] = None,
        extent: List[int] = None):
    """
    Create a binary image mask from a triangulated surface.

    Inputs
    ------
    v : an nx3 array of vertex coordinates.
    t : an mx3 array of the vertex indices of triangular faces
    imageShape : a 3-tuple of the output binary image array shape
    outputOrigin : 3D coordinates of the origin of the output image array
    outputSpacing : Voxel spacing of the output image array

    Returns
    -------
    mask_image_array : binary image array
    surf_poly : vtkPolyData instance of the triangulated surface
    """

    if outputOrigin is None:
        outputOrigin = [0.0, 0.0, 0.0]
    if outputSpacing is None:
        outputSpacing = [1.0, 1.0, 1.0]

    img_dtype = uint8

    # make into vtkPolydata
    surf_poly = tri2Polydata(v, t)

    # create mask vtkImage
    mask_image_array = ones(image_shape, dtype=img_dtype)
    mask_vtk_image_importer = array2vtkImage(
        mask_image_array, img_dtype, flipDim=False, extent=extent,
        retImporter=True
    )
    mask_vtk_image = mask_vtk_image_importer.GetOutput()

    # create stencil from polydata
    stencil_maker = vtk.vtkPolyDataToImageStencil()
    if vtk.VTK_MAJOR_VERSION < 6:
        stencil_maker.SetInput(surf_poly)
    else:
        stencil_maker.SetInputDataObject(surf_poly)
    stencil_maker.SetOutputOrigin(outputOrigin)
    stencil_maker.SetOutputSpacing(outputSpacing)
    stencil_maker.SetOutputWholeExtent(mask_vtk_image.GetExtent())
    stencil_maker.Update()  # needed in VTK 6

    stencil = vtk.vtkImageStencil()
    if vtk.VTK_MAJOR_VERSION < 6:
        stencil.SetInput(mask_vtk_image)
        stencil.SetStencil(stencil_maker.GetOutput())
    else:
        stencil.SetInputDataObject(mask_vtk_image)
        stencil.SetStencilData(stencil_maker.GetOutput())
    stencil.SetBackgroundValue(0)
    stencil.ReverseStencilOff()
    stencil.Update()

    mask_image_array2 = vtkImage2Array(stencil.GetOutput(), img_dtype, flipDim=True)
    return mask_image_array2, surf_poly


def _makeImageSpaceGF(
        scan: Scan,
        GF,
        negSpacing=False,
        zShift=True):
    """
    Transform a fieldwork geometric field from physical coords to image voxel indices
    """
    new_gf = copy.deepcopy(GF)
    p = GF.get_all_point_positions()
    p_img = scan.coord2Index(p, negSpacing=negSpacing, zShift=zShift, roundInt=False)
    new_gf.set_field_parameters(p_img.T[:, :, newaxis])

    return new_gf


def gf2BinaryMask(
        gf,
        scan: Scan,
        xiD: List[int] = None,
        negSpacing: bool = False,
        zShift: bool = True,
        outputOrigin: Tuple[float, float, float] = (0, 0, 0),
        outputSpacing: Tuple[float, float, float] = (1, 1, 1)) -> ndarray:
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
    gf_image = _makeImageSpaceGF(scan, gf, negSpacing, zShift)
    vertices, triangles = gf_image.triangulate(xiD, merge=True)
    return triSurface2BinaryMask(vertices, triangles, scan.I.shape, outputOrigin, outputSpacing)


def simplemesh2BinaryMask(
        sm: simplemesh.SimpleMesh,
        scan: Scan,
        zShift: bool = True,
        negSpacing: bool = False,
        outputOrigin: Tuple[float, float, float] = (0, 0, 0),
        outputSpacing: Tuple[float, float, float] = (1, 1, 1)) -> ndarray:
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
    v_image = scan.coord2Index(sm.v, zShift=zShift, negSpacing=negSpacing, roundInt=False)
    return triSurface2BinaryMask(v_image, sm.f, scan.I.shape, outputOrigin, outputSpacing)


def image2Simplemesh(
        image_array: ndarray,
        index_2_coord: Callable,
        iso_value: int,
        deciRatio: Optional[float] = None,
        smoothIt: int = 200,
        zShift: bool = True) -> Tuple[simplemesh.SimpleMesh, simplemesh.SimpleMesh, vtk.vtkPolyData]:
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
    if image_array.dtype != IMGDTYPE:
        raise ValueError('imageArray must be {}.'.format(IMGDTYPE))

    if vtk.VTK_MAJOR_VERSION < 6:
        vtk_image = array2vtkImage(image_array, IMGDTYPE, flipDim=True, retImporter=False)
    else:
        vtk_image = array2vtkImage(image_array, IMGDTYPE, flipDim=True, retImporter=True, pipeline=True)

    params = PolydataFromImageParams()
    params.smoothImage = False
    params.imgSmthRadius = 1.0
    params.imgSmthSD = 1.0
    params.isoValue = iso_value
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
        polydata = polydataFromImage(vtk_image, params, pipeline=False)
        verts, tris, normals = polyData2Tri(polydata, pipeline=False)
    else:
        polydata = polydataFromImage(vtk_image, params, pipeline=True)
        verts, tris, normals = polyData2Tri(polydata, pipeline=True)

    verts = verts[:, ::-1] + [0.0, 0.0, 1.0]
    sm_img = simplemesh.SimpleMesh(v=verts, f=tris)
    sm = simplemesh.SimpleMesh(v=index_2_coord(verts, zShift=zShift), f=tris)
    sm.data = {'vertexnormal': normals}
    log.debug('image-to-mesh done')
    log.debug('vertices: %s', sm.v.shape[0])
    log.debug('faces: %s', sm.f.shape[0])
    return sm, sm_img, polydata


def smoothMeshVTK(
        mesh: simplemesh.SimpleMesh,
        it: int,
        smoothboundary: bool = False,
        smoothfeatures: bool = False,
        relaxfactor: float = 1.0,
        usewsinc: bool = True) -> simplemesh.SimpleMesh:
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


def optimiseMesh(sm: simplemesh.SimpleMesh, deciratio: float, clean: bool = False) -> simplemesh.SimpleMesh:
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
        get_previous_output = cleaner.GetOutput

    # decimate polydata
    log.debug("decimating using quadric...")
    decimator = vtk.vtkQuadricDecimation()
    if vtk.VTK_MAJOR_VERSION < 6:
        if not clean:
            decimator.SetInput(poly)
        else:
            decimator.SetInput(get_previous_output())
    else:
        if not clean:
            decimator.SetInputDataObject(poly)
        else:
            decimator.SetInputDataObject(get_previous_output())
    decimator.SetTargetReduction(deciratio)
    # decimator.SetPreserveTopology(True)
    # decimator.SplittingOn()
    decimator.Update()
    get_previous_output = decimator.GetOutput

    # convert back to sm
    v, f, normals = polyData2Tri(get_previous_output())
    new_sm = simplemesh.SimpleMesh(v=v, f=f)
    return new_sm


# ====================================================#
class Colours:
    def __init__(self):
        self.colours = dict()

        red = vtk.vtkProperty()
        red.SetColor(1.0, 0.0, 0.0)
        self.colours['red'] = red

        green = vtk.vtkProperty()
        green.SetColor(0.0, 1.0, 0.0)
        self.colours['green'] = green

        blue = vtk.vtkProperty()
        blue.SetColor(0.0, 0.0, 1.0)
        self.colours['blue'] = blue

        magenta = vtk.vtkProperty()
        magenta.SetColor(1.0, 0.0, 1.0)
        self.colours['magenta'] = magenta

        yellow = vtk.vtkProperty()
        yellow.SetColor(1.0, 1.0, 0.0)
        self.colours['yellow'] = yellow

        cyan = vtk.vtkProperty()
        cyan.SetColor(0.0, 1.0, 1.0)
        self.colours['cyan'] = cyan

    def getColour(self, colourStr):
        return self.colours[colourStr]


def renderVtkImageVolume(vtkImage, cRange=(0, 255), oRange=(0, 255)):
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
