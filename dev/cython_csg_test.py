
import numpy as np
import pyximport
pyximport.install(
    setup_args={
        'include_dirs':[np.get_include(),] 
    }
)

# import cython_csg_gias2 as CSG
from gias2.mesh import cython_csg as CSG

from gias2.mesh import vtktools, simplemesh

# Vector tests
vec1 = CSG.Vector(1,0,0)
vec2 = CSG.Vector(0,1,0)
vec1_clone = vec1.clone()
vec1_neg = vec1.negated()
vec3 = vec1+vec2
vec4 = vec1-vec2
vec5 = vec1*2.0
vec6 = vec1/2.0
dot12 = vec1.dot(vec2)
vec7 = vec1.lerp(vec2, 0.5)
l1 = vec1.length()
vec8 = vec1.unit()
vec9 = vec1.cross(vec2)

# Vertex tests
vert1 = CSG.Vertex(vec1, vec2)
vert2 = CSG.Vertex(vec2, vec1)
vert3 = vert1.clone()
vert1.flip()
vert4 = vert1.interpolate(vert2, 0.5)

# plane tests
p1 = CSG.Plane(vec1, 10.0)
p2 = p1.clone()
p1.flip()

# polygon tests
pv1 = CSG.Vertex(CSG.Vector(0,0,0), CSG.Vector(0,0,1))
pv2 = CSG.Vertex(CSG.Vector(1,0,0), CSG.Vector(0,0,1))
pv3 = CSG.Vertex(CSG.Vector(0,1,0), CSG.Vector(0,0,1))
poly1 = CSG.Polygon([pv1,pv2,pv3], False)
poly2 = poly1.clone()
poly1.flip()
cutplane = CSG.Plane(CSG.Vector(-1,1,0).unit(), 0.0)
cpfront = []
cpback = []
front = []
back = []
cutplane.splitPolygon(poly1, cpfront, cpback, front, back)

# BSPNode tests
a = CSG.sphere(center=[0.5, 0.5, 0.5], radius=0.5, slices=8, stacks=4)
print('sphere csg created')
b = CSG.cylinder(start=[0.,0.,0.], end=[1.,0.,0.], radius=0.3, slices=16)
print('cylinder csg created')
union_ab_1 = a.union(b)

# shape tests
cube = CSG.cube([10,10,10], [1,2,3])
cube2 = CSG.cube([10,10,10], [5,1,1])
cube_sub = cube.subtract(cube2)
cone = CSG.cone(start=[10,0,0], end=[20,0,0], radius=5, slices=8)
cup = CSG.cup([0,0,0], [0,0,1], 9, 10, 8, 8)
trunc_cone = CSG.cylinder_var_radius(start=[0,0,0], end=[10,0,0], startr=2.0, endr=4.0, slice=8, stacks=8)


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
    vertices, faces = CSG.csg_2_polys(csgeom)
    if len(vertices)==0:
        raise ValueError('no polygons in geometry')
    return vtktools.polygons2Tri(vertices, faces, clean, normals)

def csg2simplemesh(csgeom, clean=True):
    v, f, n = get_csg_triangles(csgeom, clean=clean, normals=False)
    return simplemesh.SimpleMesh(v=v, f=f)

sphere_sm = csg2simplemesh(a)
cylinder_sm = csg2simplemesh(b)
cube_sm = csg2simplemesh(cube)
cube2_sm = csg2simplemesh(cube2)
cube_sub_sm = csg2simplemesh(cube_sub)
cone_sm = csg2simplemesh(cone)
cup_sm = csg2simplemesh(cup)
trunc_cone_sm = csg2simplemesh(trunc_cone)