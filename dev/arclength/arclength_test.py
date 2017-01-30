"""
Script to test the evaluation of arc length of a polynomial parametric
curve.

"""

import numpy as np
from scipy import integrate
import timeit
from gias2.fieldwork.field import geometric_field
from gias2.fieldwork.field import template_fields as tf

# create a template mesh
# 4 _7_ 6
#  |\  |  
# 5| 3 |8
#  |__\| 
# 0  1  2

ens = tf.two_quad_patch()
ens.flatten()
ens.mapper.has_custom_map = False
ens.mapper._custom_ensemble_order = None
ens.mapper._custom_ensemble_order_inverse = None
ens.map_parameters()

x = np.array([0, 1, 2, 1, 0, 0, 2, 1, 2])
y = np.array([0, 0, 0, 1, 2, 1, 2, 2, 1])
z = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0])

BASIS_TYPE = {'tri6':'simplex_L2_L2'}
gf = geometric_field.geometric_field(
            '2patch', 3,
            field_dimensions=2,
            field_basis=BASIS_TYPE,
            )
gf.add_element_with_parameters(ens, np.array([x,y,z])[:,:,np.newaxis])

# create a 1-d curve
curve_nodes = [[0,1,2],[2,3,4],[4,5,0],[6,7,4],[2,8,6]]
curve_gf = gf.makeLineElementsFromPointSets(curve_nodes, {3: 'line3l'}, {'line3l':'line_L2'})
l_eval = geometric_field.makeArclengthEvalDisc(curve_gf, 100)
l = l_eval(curve_gf.field_parameters)