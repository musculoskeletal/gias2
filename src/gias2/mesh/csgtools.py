"""
FILE: csgtools.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: Constructive Solid Geometry module based on PyCSG

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import numpy as np

import pyximport

pyximport.install(
    setup_args={"include_dirs": np.get_include()},
    reload_support=True,
    language_level=3
)
from gias2.mesh import cython_csg as CSG

# try:
#     from _cython_csg import CSG
#     HAS_CYTHON_CSG = True
# except ImportError:
#     try:
#         from csg.core import CSG
#     except ImportError:
#         raise ImportError('Unable to import csg or cython-csg')
#     HAS_CYTHON_CSG = False

# from csg import geom
from gias2.mesh import vtktools, simplemesh

vtk = vtktools.vtk


# class CSG_Pos(object):
#     """ A very simple implementation of PyCSG's Pos class
#     """
#     def __init__(self, x, y, z):
#         self.x = x
#         self.y = y
#         self.z = z

# def make_csg_vertex( x, y, z):
#     pos = CSG_Pos(x, y, z)
#     return geom.Vertex(pos)

def _unit(v):
    """
    return the unit vector of vector v
    """
    return v / np.sqrt((v ** 2.0).sum(-1))


def poly_2_csgeom(vertices, faces, normals=None):
    """
    Create a CSG geometry from a list of vertices and faces.

    Inputs:
    vertices: an nx3 nested list of vertices coordinates
    faces: an mxp nested list of faces

    Returns:
    geom: a csg geometry instance
    """
    if normals is None:
        normals = list(np.zeros((len(vertices), 3)))
    return CSG.poly_2_csg(vertices, faces, normals)

    # # instantiate csg vertices for all vertices
    # csg_vertices = [geom.Vertex(list(v)) for v in vertices]

    # # instantiate csg polygons for all faces
    # csg_polygons = []
    # for f in faces:
    #     face_vertices = [csg_vertices[i] for i in f]
    #     p = geom.Polygon(face_vertices)
    #     csg_polygons.append(p)

    # # create csg geom
    # return CSG.fromPolygons(csg_polygons)


def get_csg_polys(csgeom):
    """
    return the vertex coordinates and polygon vertex indices
    of a csg geometry
    """

    return CSG.csg_2_polys(csgeom)

    # polygons = csgeom.toPolygons()

    # # get vertices for each polygon
    # vertices = []
    # vertex_numbers = {}
    # faces = []
    # new_vertex_number = 0
    # for polygon in polygons:
    #     face_vertex_numbers = []
    #     for v in polygon.vertices:
    #         pos = (v.pos.x, v.pos.y, v.pos.z)
    #         vertex_number = vertex_numbers.get(pos)
    #         if vertex_number is None:
    #             vertices.append(pos)
    #             vertex_numbers[pos] = new_vertex_number
    #             vertex_number = new_vertex_number
    #             new_vertex_number += 1
    #         face_vertex_numbers.append(vertex_number)
    #     faces.append(face_vertex_numbers)

    # return vertices, faces


def get_csg_triangles(csgeom, clean=False, normals=False):
    """
    Return the vertex coordinates, triangle vertex indices, and point normals
    (if defined) of a triangulated csg geometry.

    inputs
    ======
    csgeom : CSG Solid instance
        CSG solid to be meshed
    clean : bool (default=False)
        Clean the mesh
    normals : bool (default=False)
        Calculated normals

    Returns
    =======
    v : nx3 array
        a list of vertex coordinates
    f : mx3 array
        a list of 3-tuples face vertex indices
    n : mx3 array
        a list of face normals if normals=True, else None.
    """
    vertices, faces = get_csg_polys(csgeom)
    if len(vertices) == 0:
        raise ValueError('no polygons in geometry')
    return vtktools.polygons2Tri(vertices, faces, clean, normals)


def csg2simplemesh(csgeom, clean=True):
    v, f, n = get_csg_triangles(csgeom, clean=clean, normals=False)
    return simplemesh.SimpleMesh(v=v, f=f)


def simplemesh2csg(sm):
    if sm.vertexNormals is not None:
        normals = sm.vertexNormals.tolist()
    else:
        normals = None
    return poly_2_csgeom(sm.v.tolist(), sm.f.tolist(), normals)


def cube(center=[0, 0, 0], radius=[1, 1, 1]):
    return CSG.cube(center=list(center), radius=list(radius))


def cup(centre, normal, ri, ro, slices=12, stacks=12):
    return CSG.cup(list(centre), list(normal), ri, ro, slices, stacks)


def cylinder_var_radius(**kwargs):
    """ Returns a cylinder with linearly changing radius between the two ends.
        
        Kwargs:
            start (list): Start of cylinder, default [0, -1, 0].
            
            end (list): End of cylinder, default [0, 1, 0].
            
            startr (float): Radius of cylinder at the start, default 1.0.
            
            enr (float): Radius of cylinder at the end, default 1.0.
            
            slices (int): Number of radial slices, default 16.

            stacks (int): Number of axial slices, default=2.
    """
    return CSG.cylinder_var_radius(**kwargs)

    # s = kwargs.get('start', np.array([0.0, -1.0, 0.0]))
    # e = kwargs.get('end', np.array([0.0, 1.0, 0.0]))
    # if isinstance(s, list):
    #     s = np.array(s)
    # if isinstance(e, list):
    #     e = np.array(e)
    # sr = kwargs.get('startr', 1.0)
    # er = kwargs.get('endr', 1.0)
    # slices = kwargs.get('slices', 16)
    # stacks = kwargs.get('stacks', 2)
    # stack_l = 1.0/stacks # length of each stack segment
    # ray = e - s

    # axisZ = _unit(ray)
    # isY = np.abs(axisZ[1])>0.5
    # axisX = _unit(np.cross([float(isY), float(not isY), 0], axisZ))
    # axisY = _unit(np.cross(axisX, axisZ))
    # start = geom.Vertex(list(s), list(-axisZ))
    # end = geom.Vertex(list(e), list(axisZ))
    # polygons = []
    # _verts = {}

    # def make_vert(stacki, slicei, normalBlend):
    #     stackr = stacki*stack_l
    #     slicer = slicei/float(slices)
    #     angle = slicer*np.pi*2.0
    #     out = axisX*np.cos(angle) + axisY*np.sin(angle)
    #     r = sr + stackr*(er-sr)
    #     pos = s + ray*stackr + out*r
    #     normal = out*(1.0 - np.abs(normalBlend)) + (axisZ*normalBlend)
    #     return geom.Vertex(list(pos), list(normal))  

    # def point(stacki, slicei, normalBlend):
    #     # wrap around
    #     if slicei==slices:
    #         slicei = 0

    #     # check if vertex already exists. Duplicated vertices may
    #     # cause self-intersection errors
    #     vert = _verts.get((stacki, slicei), None)
    #     if vert is None:
    #         vert = make_vert(stacki, slicei, normalBlend)
    #         _verts[(stacki, slicei)] = vert
    #     return vert

    # for i in range(0, stacks):
    #     for j in range(0, slices):
    #         # start side triangle
    #         if i==0:
    #             polygons.append(
    #                 geom.Polygon([
    #                     start,
    #                     point(i, j,   -1.), 
    #                     point(i, j+1, -1.)
    #                     ])
    #                 )
    #         # round side quad
    #         polygons.append(
    #             geom.Polygon([
    #                 point(i,   j+1, 0.),
    #                 point(i,   j,   0.),
    #                 point(i+1, j,   0.),
    #                 point(i+1, j+1, 0.)
    #                 ])
    #             )

    #         # end side triangle
    #         if i==(stacks-1):
    #             polygons.append(
    #                 geom.Polygon([
    #                     end,
    #                     point(i+1, j+1, 1.), 
    #                     point(i+1, j,   1.)
    #                     ])
    #                 )

    # return CSG.fromPolygons(polygons)
