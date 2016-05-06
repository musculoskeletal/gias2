"""
Examples and tests scratch pad for using the GNU Triangulated Surface Library
"""
import numpy as np
from gias2.mesh import vtktools, simplemesh, gtstools
from gias2.visualisation import fieldvi
from gias2.image_analysis import image_tools
import gts

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
# surf_file = '2008_1741_tibia_fibula.wrl'
surf_file = 'OCAI-001 5mm x 70mm Screw Boolean Rev 1.STL'
reader = vtktools.Reader()
reader.read(surf_file)
smesh = reader.getSimplemesh()
smesh.calcFaceProperties()

# converting between a simplemesh and a GTS surface
gts_surface = gtstools.simplemesh2gtssurf(smesh)
# gts_surface.scale(0.5, 0.5, 0.5)
smesh_gts = gtstools.gtssurf2simplemesh(gts_surface)
with open('surface.gts', 'w') as f:
    gts_surface.write(f)

# generating simple shapes
# gts_cube = gtstools.gts.cube()
gts_sphere = gtstools.gts.sphere(4)

# transform shapes
target = smesh.calcCoM()
gts_sphere.scale(10,10,10)
gts_sphere.translate(*target)
smesh_sphere = gtstools.gtssurf2simplemesh(gts_sphere)
with open('sphere.gts', 'w') as f:
    gts_sphere.write(f)

# cylinder/trunc cone
gts_cylinder = gtstools.cylinder(
    start=smesh.v.min(0),
    end=smesh.v.max(0),
    startr=5.0,
    endr=5.0,
    slices=32,
    stacks=32,
    )
smesh_cylinder = gtstools.gtssurf2simplemesh(gts_cylinder)
cylinder_face_centres = np.array([f.circumcenter().coords() for f in gts_cylinder.faces()])
cylinder_face_normals = np.array([f.normal() for f in gts_cylinder.faces()])
vtktools.Writer(v=smesh_cylinder.v, f=smesh_cylinder.f).write('cylinder.ply')
with open('cylinder.gts', 'w') as f:
    gts_cylinder.write(f)

# boolean
# gts_surface_diff = gts_surface.difference(gts_sphere)
# gts_surface_union = gts_surface.union(gts_sphere)
# gts_surface_intersect = gts_surface.intersection(gts_sphere)
# smesh_diff = gtstools.gtssurf2simplemesh(gts_surface_diff)
# smesh_union = gtstools.gtssurf2simplemesh(gts_surface_union)
# smesh_intersect = gtstools.gtssurf2simplemesh(gts_surface_intersect)

gts_surface_diff = gts_sphere.difference(gts_cylinder)
gts_surface_union = gts_sphere.union(gts_cylinder)
gts_surface_intersect = gts_sphere.intersection(gts_cylinder)
smesh_diff = gtstools.gtssurf2simplemesh(gts_surface_diff)
smesh_union = gtstools.gtssurf2simplemesh(gts_surface_union)
smesh_intersect = gtstools.gtssurf2simplemesh(gts_surface_intersect)

# isosurface
# sphere_img, vol_min, vol_max = mesh2img(smesh_sphere, 10, [2.0,2.0,2.0], retextent=True)
# gts_sphere_img = gtstools.gts.isosurface(
#     sphere_img.I, 0.5, method='cube',
#     extents=[
#         vol_min[0], vol_max[0],
#         vol_min[1], vol_max[1],
#         vol_min[2], vol_max[2],
#         ]
#     )
# smesh_sphere_img = gtstools.gtssurf2simplemesh(gts_sphere_img)

# view
v = fieldvi.Fieldvi()
v.addTri('original', smesh)
v.addTri('gts', smesh_gts)
v.addTri('gts sphere', smesh_sphere)
v.addTri('gts cylinder', smesh_cylinder)
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
