"""
Script to test methods for evaluating the arc length of a polynomial parametric
curve.

"""

import numpy as np
from scipy import integrate
import timeit
from gias2.fieldwork.field import geometric_field
from gias2.fieldwork.field import template_fields as tf

# create a template mesh
ens, x, y, z = tf.two_quad_ring(2,10,10)
BASIS_TYPE = {'tri6':'simplex_L2_L2'}
gf = geometric_field.geometric_field(
            'ring', 3,
            field_dimensions=2,
            field_basis=BASIS_TYPE,
            )
gf.add_element_with_parameters(ens, np.array([x,y,z]))

# create a 1-d curve
curve = gf.makeLineElementsFromPoints([0,1,2], 3, {'line3l':'line_L2'})

#=============================================================================#
# METHOD 1: sum of lengths of discretised linear segments
#=============================================================================#
def make_M1(c, d):
    ceval = geometric_field.makeGeometricFieldEvaluatorSparse(c, [d,])
    p = np.array(c.field_parameters)

    def M1():
        x = ceval(p).T
        # return np.sum(((x[1:,:] - x[:-1,:])**2.0).sum(1))
        return np.sum(np.sqrt(((x[1:,:] - x[:-1,:])**2.0).sum(1)))

    return M1

m1 = make_M1(curve, 100)
m1_res = m1()
m1_time = timeit.timeit(m1, number=100)
print('Method 1 - l: {}, time: {}'.format(m1_res, m1_time))

#=============================================================================#
# Common things for numeric integration methods
#=============================================================================#

def make_quad_integrand(x, y, z):
    
    def integrand(t):
        return np.sqrt(
            ((4.0*t-3.0)*z[0]+(4.0-8.0*t)*z[1]+(4.0*t-1.0)*z[2])**2.0+\
            ((4.0*t-3.0)*y[0]+(4.0-8.0*t)*y[1]+(4.0*t-1.0)*y[2])**2.0+\
            ((4.0*t-3.0)*x[0]+(4.0-8.0*t)*x[1]+(4.0*t-1.0)*x[2])**2.0
            )

    return integrand

#=============================================================================#
# METHOD 2: SciPy Quad
#=============================================================================#
def make_M2(c):
    integ = make_quad_integrand(
                c.field_parameters[0].squeeze(),
                c.field_parameters[1].squeeze(),
                c.field_parameters[2].squeeze(),
                )

    def M2():
        return integrate.quad(integ, 0.0, 1.0)

    return M2

m2 = make_M2(curve)
m2_res = m2()
m2_time = timeit.timeit(m2, number=100)
print('Method 2 - l: {}, time: {}'.format(m2_res[0], m2_time))

#=============================================================================#
# METHOD 3: SciPy simps. n must be odd
#=============================================================================#
def make_M3(c, n):
    x = np.linspace(0,1,n)
    integ = make_quad_integrand(
                c.field_parameters[0].squeeze(),
                c.field_parameters[1].squeeze(),
                c.field_parameters[2].squeeze(),
                )
    y = integ(x)

    def M3():
        return integrate.simps(y, x)

    return M3

m3 = make_M3(curve, 101)
m3_res = m3()
m3_time = timeit.timeit(m3, number=100)
print('Method 3 - l: {}, time: {}'.format(m3_res, m3_time))

#=============================================================================#
# METHOD 4: SciPy romb. n must be 2**n+1, e.g. 3, 5, 9, 17, 33
#=============================================================================#
def make_M4(c, n):
    x = np.linspace(0,1,n)
    integ = make_quad_integrand(
                c.field_parameters[0].squeeze(),
                c.field_parameters[1].squeeze(),
                c.field_parameters[2].squeeze(),
                )
    y = integ(x)
    dx = 1.0/n
    def M4():
        return integrate.romb(y, dx=dx)

    return M4

m4 = make_M4(curve, 129)
m4_res = m4()
m4_time = timeit.timeit(m4, number=100)
print('Method 4 - l: {}, time: {}'.format(m4_res, m4_time))

#=============================================================================#
# METHOD 5: SciPy quad with ctypes
#=============================================================================#
import ctypes

lib = ctypes.CDLL('arclengthlib.so')   # Use absolute path to testlib
func = lib.quad   # Assign specific function to name func (for simplicity)
func.restype = ctypes.c_double
func.argtypes = (ctypes.c_int, ctypes.c_double)

def make_M5(c):
    # p = c.field_parameters.squeeze()
    # args = (
    #         p[0,0], p[0,1], p[0,2],
    #         p[1,0], p[1,1], p[1,2],
    #         p[2,0], p[2,1], p[2,2],
    #         )
    args = tuple(c.field_parameters.ravel())

    def M5():
        return integrate.quad(func, 0.0, 1.0, args=args)

    return M5

m5 = make_M5(curve)
m5_res = m5()
m5_time = timeit.timeit(m5, number=100)
print('Method 5 - l: {}, time: {}'.format(m5_res[0], m5_time))