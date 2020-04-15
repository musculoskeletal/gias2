"""
Functions for using the GNU GTS library
"""
import gts
import numpy as np

from gias2.mesh import simplemesh


def _unit(v):
    """
    return the unit vector of vector v
    """
    return v / np.sqrt((v ** 2.0).sum(-1))


def simplemesh2gtssurf(smesh):
    """
    Create a GTS surface from a SimpleMesh.

    input
    =====
    smesh : SimpleMesh instance

    returns
    =======
    gts_surf : gts.Surface instance
    """
    ## create vertices
    gts_verts = [gts.Vertex(*vi) for vi in smesh.v]

    ## create surface
    gts_surf = gts.Surface()
    for fi in smesh.f:
        ## create edges, anticlockwise by indices in f
        e1 = gts.Edge(gts_verts[fi[0]], gts_verts[fi[1]])
        e2 = gts.Edge(gts_verts[fi[1]], gts_verts[fi[2]])
        e3 = gts.Edge(gts_verts[fi[2]], gts_verts[fi[0]])

        ## create face
        gts_surf.add(gts.Face(e1, e2, e3))

    return gts_surf


def gtssurf2simplemesh(gts_surf):
    """
    Create a SimpleMesh from a GTS surface.

    input
    =====
    gts_surf : gts.Surface instance

    returns
    =======
    smesh : SimpleMesh instance
    """
    x, y, z, f = gts.get_coords_and_face_indices(gts_surf, True)
    v = np.array([x, y, z]).T
    return simplemesh.SimpleMesh(v=v, f=f)


def cylinder(**kwargs):
    """ Returns a cylinder with linearly changing radius between the two ends.
        
        Kwargs:
            start (list): Start of cylinder, default [0, -1, 0].
            
            end (list): End of cylinder, default [0, 1, 0].
            
            startr (float): Radius of cylinder at the start, default 1.0.
            
            endr (float): Radius of cylinder at the end, default 1.0.
            
            slices (int): Number of radial slices, default 16.

            stacks (int): Number of axial slices, default=2.
    """
    s = np.array(kwargs.get('start', (0.0, -1.0, 0.0)))
    e = np.array(kwargs.get('end', (0.0, 1.0, 0.0)))
    sr = kwargs.get('startr', 1.0)
    er = kwargs.get('endr', 1.0)
    slices = kwargs.get('slices', 16)
    stacks = kwargs.get('stacks', 2)
    stack_l = 1.0 / stacks  # length of each stack segment
    ray = e - s

    axisZ = _unit(ray)
    isY = np.abs(axisZ[1]) > 0.5
    axisX = _unit(np.cross([float(isY), float(not isY), 0], axisZ))
    axisY = _unit(np.cross(axisX, axisZ))
    start = gts.Vertex(*s)
    end = gts.Vertex(*e)
    surf = gts.Surface()
    _verts = {}

    def make_vert(stacki, slicei):
        stackr = stacki * stack_l
        slicer = slicei / float(slices)
        angle = slicer * np.pi * 2.0
        out = axisX * np.cos(angle) + axisY * np.sin(angle)
        r = sr + stackr * (er - sr)
        pos = s + ray * stackr + out * r
        return gts.Vertex(*pos)

    def point(stacki, slicei):
        # wrap around
        if slicei == slices:
            slicei = 0

        # check if vertex already exists. Duplicated vertices will
        # cause self-intersection errors
        vert = _verts.get((stacki, slicei), None)
        if vert is None:
            vert = make_vert(stacki, slicei)
            _verts[(stacki, slicei)] = vert
        return vert

    def add_triangle(v1, v2, v3):
        e1 = gts.Edge(v1, v2)
        e2 = gts.Edge(v2, v3)
        e3 = gts.Edge(v3, v1)
        surf.add(gts.Face(e1, e2, e3))

    for i in range(0, stacks):
        # through angles
        for j in range(0, slices):
            # start side triangle
            if i == 0:
                add_triangle(start, point(i, j), point(i, j + 1))

            # round side tris
            add_triangle(point(i, j + 1), point(i, j), point(i + 1, j))
            add_triangle(point(i, j + 1), point(i + 1, j), point(i + 1, j + 1))

            # end side triangle
            if i == (stacks - 1):
                add_triangle(end, point(i + 1, j + 1), point(i + 1, j))

    return surf
