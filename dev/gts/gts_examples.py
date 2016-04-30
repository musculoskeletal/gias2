"""
Examples and tests scratch pad for using the GNU Triangulated Surface Library
"""
import numpy as np
from gias2.mesh import vtktools, simplemesh, gtstools
from gias2.visualisation import fieldvi
import gts

#=============================================================================#
# example surface
surf_file = '2008_1741_tibia_fibula.wrl'
reader = vtktools.Reader()
reader.read(surf_file)
smesh = reader.getSimplemesh()

# converting between a simplemesh and a GTS surface
gts_surface = gtstools.simplemesh2gtssurf(smesh)
smesh_gts = gtstools.gtssurf2simplemesh(gts_surface)

# generating simple shapes
# gts_cube = gtstools.gts.cube()
gts_sphere = gtstools.gts.sphere(4)

# transform shapes
target = smesh.v.mean(0)
gts_sphere.scale(50,50,50)
gts_sphere.translate(*target)
smesh_sphere = gtstools.gtssurf2simplemesh(gts_sphere)

# cylinder/trunc cone
gts_cylinder = gtstools.cylinder(
    start=smesh.v.min(0),
    end=smesh.v.max(0),
    startr=10.0,
    endr=10.0,
    slices=32,
    stacks=50,
    )
smesh_cylinder = gtstools.gtssurf2simplemesh(gts_cylinder)
vtktools.Writer(v=smesh_cylinder.v, f=smesh_cylinder.f).write('cylinder.ply')

# boolean
gts_surface_diff = gts_surface.difference(gts_cylinder)
gts_surface_union = gts_surface.union(gts_cylinder)
smesh_diff = gtstools.gtssurf2simplemesh(gts_surface_diff)
smesh_union = gtstools.gtssurf2simplemesh(gts_surface_union)

# view
v = fieldvi.Fieldvi()
v.addTri('original', smesh)
v.addTri('gts', smesh_gts)
v.addTri('gts sphere', smesh_sphere)
v.addTri('gts cylinder', smesh_cylinder)
v.addTri('diff', smesh_diff)
v.addTri('union', smesh_union)
v.configure_traits()
