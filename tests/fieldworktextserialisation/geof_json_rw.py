from gias2.fieldwork.field import geometric_field as gf
from gias2.visualisation.fieldvi import Fieldvi

reload(gf)

m1 = gf.load_geometric_field(
    'data/femur_right_mean_cortex_outer_2012-05-14_July2011TS_2_AB_rigid.geof',
    'data/femur_right_quartic_flat.ens',
    'data/femur_right_quartic_flat.mesh'
)
m1.save_geometric_field(
    'data/femur_right_mean_cortex_outer_2012-05-14_July2011TS_2_AB_rigid_json',
    'data/femur_right_quartic_flat_json',
    'data/femur_right_quartic_flat_json'
)
m2 = gf.load_geometric_field(
    'data/femur_right_mean_cortex_outer_2012-05-14_July2011TS_2_AB_rigid_json.geof',
    'data/femur_right_quartic_flat_json.ens',
    'data/femur_right_quartic_flat_json.mesh'
)

assert m1.name == m2.name
assert m1.dimensions == m2.dimensions
assert m1.ensemble_point_counter == m2.ensemble_point_counter
assert (m1.field_parameters - m2.field_parameters).sum() == 0

assert m1.ensemble_field_function.name == m2.ensemble_field_function.name
assert m1.ensemble_field_function.dimensions == m2.ensemble_field_function.dimensions
assert list(m1.ensemble_field_function.basis.keys()) == list(m2.ensemble_field_function.basis.keys())
assert m1.ensemble_field_function.subfields == m2.ensemble_field_function.subfields
assert m1.ensemble_field_function.subfield_counter == m2.ensemble_field_function.subfield_counter
assert m1.ensemble_field_function.mapper._custom_ensemble_order == m2.ensemble_field_function.mapper._custom_ensemble_order

assert m1.ensemble_field_function.mesh.name == m2.ensemble_field_function.mesh.name
assert m1.ensemble_field_function.mesh.dimensions == m2.ensemble_field_function.mesh.dimensions
assert m1.ensemble_field_function.mesh.element_counter == m2.ensemble_field_function.mesh.element_counter
assert m1.ensemble_field_function.mesh.number_of_points == m2.ensemble_field_function.mesh.number_of_points
assert m1.ensemble_field_function.mesh.is_element == m2.ensemble_field_function.mesh.is_element
assert m1.ensemble_field_function.mesh.elements == m2.ensemble_field_function.mesh.elements
assert m1.ensemble_field_function.mesh.connectivity == m2.ensemble_field_function.mesh.connectivity
assert m1.ensemble_field_function.mesh.submesh_counter == m2.ensemble_field_function.mesh.submesh_counter

v = Fieldvi()
v.displayGFNodes = False
v.GFD = [8, 8]
v.addGeometricField('m1', m1, gf.makeGeometricFieldEvaluatorSparse(m1, v.GFD))
v.addGeometricField('m2', m2, gf.makeGeometricFieldEvaluatorSparse(m2, v.GFD))
v.configure_traits()
