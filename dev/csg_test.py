
import pstats, cProfile

import numpy as np
import pyximport
pyximport.install(
    setup_args={
        'include_dirs':[np.get_include(),] 
    }
)
import cython_csg_gias2 as gCSG
from _cython_csg import CSG, BSPNode

# #=============================================================================#
# # naive cython
# #=============================================================================#
# # BSPNode tests
# a = CSG.sphere(center=[0.5, 0.5, 0.5], radius=0.5, slices=8, stacks=4)
# print('sphere csg created')
# b = CSG.cylinder(start=[0.,0.,0.], end=[1.,0.,0.], radius=0.3, slices=16)
# print('cylinder csg created')
# union_ab_1 = a.union(b)

# #=============================================================================#
# # gias2 cython
# #=============================================================================#
# BSPNode tests
ga = gCSG.sphere(center=[0.5, 0.5, 0.5], radius=0.5, slices=8, stacks=4)
print('sphere csg created')
gb = gCSG.cylinder(start=[0.,0.,0.], end=[1.,0.,0.], radius=0.3, slices=16)
print('cylinder csg created')
g_union_ab_1 = ga.union(gb)

def gunion():
	ga = gCSG.sphere(center=[0.5, 0.5, 0.5], radius=0.5, slices=8, stacks=4)
	print('sphere csg created')
	gb = gCSG.cylinder(start=[0.,0.,0.], end=[1.,0.,0.], radius=0.3, slices=16)
	print('cylinder csg created')
	g_union_ab_1 = ga.union(gb)

cProfile.runctx("gunion()", globals(), locals(), "gunion_profile.prof")
s = pstats.Stats("gunion_profile.prof")
s.strip_dirs().sort_stats("time").print_stats()