"""
Examples and tests scratch pad for using the PyCSG Library
"""
import numpy as np
from gias2.mesh import vtktools, simplemesh, csgtools
from gias2.visualisation import fieldvi
from gias2.image_analysis import image_tools
from csg.core import CSG
from csg import geom

def mesh2img(mesh, pad, voxel_spacing, vol_max=None, vol_min=None, retextent=False):
    scan = image_tools.Scan('dummy')
    if vol_max is None:
        vol_max = mesh.v.max(0) + pad
    if vol_min is None:
        vol_min = mesh.v.min(0) - pad
    voxel_spacing = np.array(voxel_spacing)
    voxel_origin = vol_min
    img_size = ((vol_max - vol_min)/voxel_spacing).round().astype(int)
    img_empty = np.zeros(img_size, dtype=int)
    scan.setImageArray(img_empty, voxel_spacing, voxel_origin)

    img = vtktools.simplemesh2BinaryMask(
        mesh, scan, zShift=False, negSpacing=False
        )[0]
    scan.setImageArray(img, voxel_spacing, voxel_origin)
    if retextent:
        return scan, vol_min, vol_max
    else:
        return scan

#=============================================================================#
# example surface
surf_file = '2008_1741_tibia_fibula.wrl'
# surf_file = 'OCAI-001 5mm x 70mm Screw Boolean Rev 1.STL'
reader = vtktools.Reader()
reader.read(surf_file)
smesh = reader.getSimplemesh()
smesh.calcFaceProperties()

# converting between a simplemesh and a GTS surface
csg_surface = csgtools.simplemesh2csg(smesh)
# csg_surface.scale(0.5, 0.5, 0.5)
smesh_csg = csgtools.csg2simplemesh(csg_surface)
# with open('surface.csg', 'w') as f:
#     csg_surface.write(f)

# generating simple shapes
# csg_cube = csgtools.csg.cube()


# transform shapes
target = smesh.calcCoM()
csg_sphere = CSG.sphere(center=target.tolist(), radius=30, slices=16, stacks=16)
smesh_sphere = csgtools.csg2simplemesh(csg_sphere)
# with open('sphere.csg', 'w') as f:
    # csg_sphere.write(f)

# cylinder/trunc cone
csg_cylinder = csgtools.cylinder_var_radius(
    start=smesh.v.min(0),
    end=smesh.v.max(0),
    startr=3.0,
    endr=6.0,
    slices=16,
    stacks=64,
    )
smesh_cylinder = csgtools.csg2simplemesh(csg_cylinder)
# cylinder_face_centres = np.array([f.circumcenter().coords() for f in csg_cylinder.faces()])
# cylinder_face_normals = np.array([f.normal() for f in csg_cylinder.faces()])
# vtktools.Writer(v=smesh_cylinder.v, f=smesh_cylinder.f).write('cylinder.ply')
# with open('cylinder.csg', 'w') as f:
#     csg_cylinder.write(f)

# boolean
# csg_surface_diff = csg_surface.difference(csg_sphere)
# csg_surface_union = csg_surface.union(csg_sphere)
# csg_surface_intersect = csg_surface.intersection(csg_sphere)
# smesh_diff = csgtools.csgsurf2simplemesh(csg_surface_diff)
# smesh_union = csgtools.csgsurf2simplemesh(csg_surface_union)
# smesh_intersect = csgtools.csgsurf2simplemesh(csg_surface_intersect)

csg_surface_diff = csg_sphere.subtract(csg_cylinder)
csg_surface_union = csg_sphere.union(csg_cylinder)
csg_surface_intersect = csg_sphere.intersect(csg_cylinder)
smesh_diff = csgtools.csg2simplemesh(csg_surface_diff)
smesh_union = csgtools.csg2simplemesh(csg_surface_union)
smesh_intersect = csgtools.csg2simplemesh(csg_surface_intersect)

# isosurface
# sphere_img, vol_min, vol_max = mesh2img(smesh_sphere, 10, [2.0,2.0,2.0], retextent=True)
# csg_sphere_img = csgtools.csg.isosurface(
#     sphere_img.I, 0.5, method='cube',
#     extents=[
#         vol_min[0], vol_max[0],
#         vol_min[1], vol_max[1],
#         vol_min[2], vol_max[2],
#         ]
#     )
# smesh_sphere_img = csgtools.csgsurf2simplemesh(csg_sphere_img)

# view
v = fieldvi.Fieldvi()
v.addTri('original', smesh)
v.addTri('csg', smesh_csg)
v.addTri('csg sphere', smesh_sphere)
v.addTri('csg cylinder', smesh_cylinder)
v.addTri('diff', smesh_diff)
v.addTri('union', smesh_union)
v.addTri('intersect', smesh_intersect)
# v.addTri('sphere_img', smesh_sphere_img)
v.configure_traits()

# v.scene.mlab.quiver3d(
#     cylinder_face_centres[:,0],
#     cylinder_face_centres[:,1],
#     cylinder_face_centres[:,2],
#     cylinder_face_normals[:,0],
#     cylinder_face_normals[:,1],
#     cylinder_face_normals[:,2],
#     )
