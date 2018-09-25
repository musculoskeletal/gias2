
import numpy as np
import pyximport
pyximport.install(
    setup_args={
        'include_dirs':[np.get_include(),] 
    }
)

import cython_csg_gias2 as CSG

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

bspa = CSG.BSPNode(a.clone().polygons)
bspb = CSG.BSPNode(b.clone().polygons)
bspa.clipTo(bspb)
bspb.clipTo(bspa)
bspb.invert()
bspb.clipTo(bspa)
bspb.invert()
bspa.build(bspb.allPolygons())
union_ab = CSG.csgFromPolygons(bspa.allPolygons())

# union_ab_1 = a.union(b)