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
from typing import List, Optional, Tuple, Iterable

import numpy as np
import pyximport

from gias2.mesh import vtktools, simplemesh
from gias2.mesh.simplemesh import SimpleMesh

pyximport.install(
    setup_args={"include_dirs": np.get_include()},
    language_level=3
)
from gias2.mesh import cython_csg as CSG

vtk = vtktools.vtk


def _unit(v: np.ndarray) -> np.ndarray:
    """
    return the unit vector of vector v
    """
    return v / np.sqrt((v ** 2.0).sum(-1))


def poly_2_csgeom(
        vertices: List[List[float]],
        faces: List[List[int]],
        normals: Optional[List[List[float]]] = None) -> CSG.CSG:
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


def get_csg_polys(csgeom: CSG.CSG) -> Tuple[List[List[float]], List[List[int]]]:
    """
    return the vertex coordinates and polygon vertex indices
    of a csg geometry
    """

    return CSG.csg_2_polys(csgeom)


def get_csg_triangles(
        csgeom: CSG.CSG,
        clean: bool = False,
        normals: bool = False) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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


def csg2simplemesh(csgeom: CSG.CSG, clean: bool = True) -> SimpleMesh:
    v, f, n = get_csg_triangles(csgeom, clean=clean, normals=False)
    return simplemesh.SimpleMesh(v=v, f=f)


def simplemesh2csg(sm: SimpleMesh) -> CSG.CSG:
    if sm.vertexNormals is not None:
        normals = sm.vertexNormals.tolist()
    else:
        normals = None
    return poly_2_csgeom(sm.v.tolist(), sm.f.tolist(), normals)


def cube(
        center: Tuple[float, float, float] = (0, 0, 0),
        radius: Tuple[float, float, float] = (1, 1, 1)) -> CSG.CSG:
    return CSG.cube(center=list(center), radius=list(radius))


def cup(
        centre: Iterable[float],
        normal: Iterable[float],
        ri: float,
        ro: float,
        slices: int = 12,
        stacks: int = 12) -> CSG.CSG:
    return CSG.cup(list(centre), list(normal), ri, ro, slices, stacks)


def cylinder_var_radius(**kwargs) -> CSG.CSG:
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
