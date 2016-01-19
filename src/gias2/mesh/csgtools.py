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
from csg.core import CSG
from csg import geom
from gias2.mesh import vtktools
from gias2.common import math
vtk = vtktools.vtk

class CSG_Pos(object):
    """ A very simple implementation of PyCSG's Pos class
    """
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

def make_csg_vertex( x, y, z):
    pos = CSG_Pos(x, y, z)
    return geom.Vertex(pos)

def poly_2_csgeom(vertices, faces):
    """
    Create a CSG geometry from a list of vertices and faces.

    Inputs:
    vertices: an nx3 array of vertices coordinates
    faces: an mxp array of faces

    Returns:
    geom: a csg geometry instance
    """

    # instantiate csg vertices for all vertices
    csg_vertices = [make_csg_vertex(v[0],v[1],v[2]) for v in vertices]
    
    # instantiate csg polygons for all faces
    csg_polygons = []
    for f in faces:
        face_vertices = [csg_vertices[i] for i in f]
        p = geom.Polygon(face_vertices)
        csg_polygons.append(p)

    # create csg geom
    return CSG.fromPolygons(csg_polygons)

def get_csg_polys(csgeom):
    """
    return the vertex coordinates and polygon vertex indices
    of a csg geometry
    """

    polygons = csgeom.toPolygons()

    # get vertices for each polygon
    vertices = []
    vertex_numbers = {}
    faces = []
    new_vertex_number = 0
    for polygon in polygons:
        face_vertex_numbers = []
        for v in polygon.vertices:
            pos = (v.pos.x, v.pos.y, v.pos.z)
            vertex_number = vertex_numbers.get(pos)
            if vertex_number is None:
                vertices.append(pos)
                vertex_numbers[pos] = new_vertex_number
                vertex_number = new_vertex_number
                new_vertex_number += 1
            face_vertex_numbers.append(vertex_number)
        faces.append(face_vertex_numbers)

    return vertices, faces

def get_csg_triangles(csgeom, clean=False, normals=False):
    """
    Return the vertex coordinates, triangle vertex indices, and point normals
    (if defined) of a triangulated csg geometry.

    Returns a list of vertex coordinates, a list of 3-tuples, and a list of 
    face normalse if normals=True, else last return variable is None.
    """
    vertices, faces = get_csg_polys(csgeom)
    return vtktools.polygons2Tri(vertices, faces, clean, normals)

def cup(centre, normal, ri, ro):
    centre = np.array(centre)
    normal = math.norm(np.array(normal))

    # create outer sphere
    sphere_out = CSG.sphere(center=list(centre), radius=ro)

    # create inner sphere
    sphere_in = CSG.sphere(center=list(centre), radius=ri)

    # create shell
    shell = sphere_out.subtract(sphere_in)
    shell_poly = shell.toPolygons()

    # create cylinder to cut shell
    cylinder = CSG.cylinder(start=list(centre),
                            end=list(centre-normal*(ro*1.5)),
                            radius=ro*1.5
                            )
    # create cup
    cup = shell.subtract(cylinder)

    return cup
    

