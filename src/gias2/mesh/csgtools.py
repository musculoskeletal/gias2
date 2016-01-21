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

import math
import numpy as np
from csg.core import CSG
from csg import geom
from gias2.mesh import vtktools
from gias2.common import math as gmath
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
    csg_vertices = [geom.Vertex(v) for v in vertices]
    
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
    if len(vertices)==0:
        raise ValueError('no polygons in geometry')
    return vtktools.polygons2Tri(vertices, faces, clean, normals)

def cup(centre, normal, ri, ro):
    centre = np.array(centre)
    normal = gmath.norm(np.array(normal))

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
    

def cylinder_var_radius(**kwargs):
    """ Returns a cylinder with linearly changing radius between the two ends.
        
        Kwargs:
            start (list): Start of cylinder, default [0, -1, 0].
            
            end (list): End of cylinder, default [0, 1, 0].
            
            startr (float): Radius of cylinder at the start, default 1.0.
            
            enr (float): Radius of cylinder at the end, default 1.0.
            
            slices (int): Number of slices, default 16.
    """
    s = kwargs.get('start', geom.Vector(0.0, -1.0, 0.0))
    e = kwargs.get('end', geom.Vector(0.0, 1.0, 0.0))
    if isinstance(s, list):
        s = geom.Vector(*s)
    if isinstance(e, list):
        e = geom.Vector(*e)
    sr = kwargs.get('startr', 1.0)
    er = kwargs.get('endr', 1.0)
    slices = kwargs.get('slices', 16)
    ray = e.minus(s)
    
    axisZ = ray.unit()
    isY = (math.fabs(axisZ.y) > 0.5)
    axisX = geom.Vector(float(isY), float(not isY), 0).cross(axisZ).unit()
    axisY = axisX.cross(axisZ).unit()
    start = geom.Vertex(s, axisZ.negated())
    end = geom.Vertex(e, axisZ.unit())
    polygons = []
    
    def point(stack, slice, normalBlend):
        angle = slice * math.pi * 2.0
        out = axisX.times(math.cos(angle)).plus(
            axisY.times(math.sin(angle)))
        if stack==0:
            r = sr
        else:
            r = er
        pos = s.plus(ray.times(stack)).plus(out.times(r))
        normal = out.times(1.0 - math.fabs(normalBlend)).plus(
            axisZ.times(normalBlend))
        return geom.Vertex(pos, normal)
        
    for i in range(0, slices):
        t0 = i / float(slices)
        t1 = (i + 1) / float(slices)
        # start side triangle
        polygons.append(geom.Polygon([start, point(0., t0, -1.), 
                                      point(0., t1, -1.)]))
        # round side quad
        polygons.append(geom.Polygon([point(0., t1, 0.), point(0., t0, 0.),
                                      point(1., t0, 0.), point(1., t1, 0.)]))
        # end side triangle
        polygons.append(geom.Polygon([end, point(1., t1, 1.), 
                                      point(1., t0, 1.)]))
    
    return CSG.fromPolygons(polygons)