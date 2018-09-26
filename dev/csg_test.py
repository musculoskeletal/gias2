
from _cython_csg import CSG, BSPNode

# BSPNode tests
a = CSG.sphere(center=[0.5, 0.5, 0.5], radius=0.5, slices=8, stacks=4)
print('sphere csg created')
b = CSG.cylinder(start=[0.,0.,0.], end=[1.,0.,0.], radius=0.3, slices=16)
print('cylinder csg created')

# bspa = BSPNode(a.clone().polygons)
# bspb = BSPNode(b.clone().polygons)
# bspa.clipTo(bspb)
# bspb.clipTo(bspa)
# bspb.invert()
# bspb.clipTo(bspa)
# bspb.invert()
# bspa.build(bspb.allPolygons())
# union_ab = CSG.fromPolygons(bspa.allPolygons())

union_ab_1 = a.union(b)