"""
Examples and tests scratch pad for using the GNU Triangulated Surface Library
"""

from gias2.mesh import vtktools
import gts

# example surface
surf_file = '2008_1741_tibia_fibula.wrl'
reader = vtktools.Reader()
reader.read(surf_file)
smesh = reader.getSimplemesh()

# converting between a simplemesh and a GTS surface

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

# read out vertices and faces
coords = np.array([v.coords() for f in gts_surf.vertices()])