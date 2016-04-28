
from gias2.fieldwork.field.topology import mesh
# reload(mesh.mesh_io)
# reload(mesh)

m1 = mesh.load_mesh('data/femur_right_quartic_flat.mesh')
m1.save_mesh('data/femur_right_quartic_flat_json')
m2 = mesh.load_mesh('data/femur_right_quartic_flat_json.mesh')

assert m1.name==m2.name
assert m1.dimensions==m2.dimensions
assert m1.element_counter==m2.element_counter
assert m1.number_of_points==m2.number_of_points
assert m1.is_element==m2.is_element
assert m1.elements==m2.elements
assert m1.connectivity==m2.connectivity
assert m1.submesh_counter==m2.submesh_counter