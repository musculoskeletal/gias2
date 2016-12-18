"""
Script to test methods for evaluating the arc lengths of edges on a piecewise
parametric lagrange mesh.

"""

import numpy as np
from scipy import integrate
import timeit
from gias2.fieldwork.field import geometric_field
from gias2.fieldwork.field import template_fields as tf

#=============================================================================#
# create a mesh
#=============================================================================#
ens, x, y, z = tf.two_quad_ring(2,10,10)
BASIS_TYPE = {'tri6':'simplex_L2_L2'}
gf = geometric_field.geometric_field(
            'ring', 3,
            field_dimensions=2,
            field_basis=BASIS_TYPE,
            )
gf.add_element_with_parameters(ens, np.array([x,y,z]))

# define curve mesh
curves_nodes = (
    (0,1,2), (2,3,0),
    (0,4,8), (2,5,8), (2,6,10), (0,7,10),
    (8,9,10), (10,11,8),
    )

# create a 1-d curve
curves_gf = gf.makeLineElementsFromPointSets(
                curves_nodes,
                {3: 'line3l'},
                {'line3l': 'line_L2'}
                )

#=============================================================================#
# METHOD 1: sum of lengths of discretised linear segments
#=============================================================================#
def make_M1(c, d):
    """
    Evaluate all element arclengths by discretisation
    """

    ceval = geometric_field.makeGeometricFieldEvaluatorSparse(c, [d,])
    p = np.array(c.field_parameters)
    n_elems = c.ensemble_field_function.mesh.get_number_of_true_elements()

    def M1():
        x = ceval(p).T
        
        # separate x into element points per element
        _x = x.reshape((n_elems,-1,3))
        return np.sum(
            np.sqrt(
                ((_x[:,1:,:] - _x[:,:-1,:])**2.0).sum(2)
                ), 1)

    return M1

m1 = make_M1(curves_gf, 1000)
m1_res = m1()
m1_time = timeit.timeit(m1, number=1000)
print('Method 1 - l: {}, time: {}'.format(m1_res, m1_time))

#=============================================================================#
# METHOD 5: SciPy quad with ctypes
#=============================================================================#
import ctypes

lib = ctypes.CDLL('arclengthlib.so')   # Use absolute path to testlib
func = lib.quad   # Assign specific function to name func (for simplicity)
func.restype = ctypes.c_double
func.argtypes = (ctypes.c_int, ctypes.c_double)

def make_M5(c):
    X = c.field_parameters[0]
    Y = c.field_parameters[1]
    Z = c.field_parameters[2]
    mapper = c.ensemble_field_function.mapper
    # get tuples of the parameters of each element
    elem_params = []
    for e in c.ensemble_field_function.mesh.elements:
        elem_params.append(
            tuple(
                np.hstack(
                    [mapper.get_element_parameters(e, p) for p in c.field_parameters]
                    )
                )
            )


    def M5():
        return [integrate.quad(func, 0.0, 1.0, args=args)[0] for args in elem_params]

    return M5

m5 = make_M5(curves_gf)
m5_res = m5()
m5_time = timeit.timeit(m5, number=1000)
print('Method 5 - l: {}, time: {}'.format(m5_res, m5_time))
