from gias2.fieldwork.field import ensemble_field_function as eff

# reload(eff.ensemble_field_function_io)
reload(eff)

m1 = eff.load_ensemble(
    'data/femur_right_quartic_flat.ens',
    'data/femur_right_quartic_flat.mesh'
)
m1.save_ensemble(
    'data/femur_right_quartic_flat_json',
    'data/femur_right_quartic_flat_json'
)
m2 = eff.load_ensemble(
    'data/femur_right_quartic_flat_json.ens',
    'data/femur_right_quartic_flat.mesh'
)

assert m1.name == m2.name
assert m1.dimensions == m2.dimensions
# assert m1.basis==m2.basis
assert m1.subfields == m2.subfields
assert m1.subfield_counter == m2.subfield_counter
assert m1.mapper._custom_ensemble_order == m2.mapper._custom_ensemble_order

assert m1.mesh.name == m2.mesh.name
assert m1.mesh.dimensions == m2.mesh.dimensions
assert m1.mesh.element_counter == m2.mesh.element_counter
assert m1.mesh.number_of_points == m2.mesh.number_of_points
assert m1.mesh.is_element == m2.mesh.is_element
assert m1.mesh.elements == m2.mesh.elements
assert m1.mesh.connectivity == m2.mesh.connectivity
assert m1.mesh.submesh_counter == m2.mesh.submesh_counter
