"""
FILE: basis.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: Basis functions for piecewise parametric meshes.
    
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
"""
import logging

from numpy import array, dot, kron, zeros, ones, sqrt, where, newaxis
from scipy.linalg import det

log = logging.getLogger(__name__)

# basis matrices of implemented bases for tensor product basis 
basis_matrices_map = {'linear_lagrange': array([[1.0, -1.0],
                                                [0.0, 1.0]]),
                      'quadratic_lagrange': array([[1.0, -3.0, 2.0],
                                                   [0.0, 4.0, -4.0],
                                                   [0.0, -1.0, 2.0]]),
                      'cubic_lagrange': 0.5 * array([[2.0, -11.0, 18.0, -9.0],
                                                     [0.0, 18.0, -45.0, 27.0],
                                                     [0.0, -9.0, 36.0, -27.0],
                                                     [0.0, 2.0, -9.0, 9.0]]),
                      'cubic_hermite': array([[1.0, 0.0, -3.0, 2.0],
                                              [0.0, 1.0, -2.0, 1.0],
                                              [0.0, 0.0, 3.0, -2.0],
                                              [0.0, 0.0, -1.0, 1.0]])}


# ======================================================================#

# 1D basis functions
def L1(x):
    return array([1.0 - x, x])


def L1d(x):
    return array([[-1.0] * len(x), [1.0] * len(x)])


def L1dd(x):
    return array([[0.0] * len(x), [0.0] * len(x)])


def L2(x):
    x2 = x * x
    return array([1.0 - 3.0 * x + 2.0 * x2,
                  + 4.0 * x - 4.0 * x2,
                  - x + 2.0 * x2])


def L2d(x):
    return array([- 3.0 + 4.0 * x,
                  + 4.0 - 8.0 * x,
                  - 1.0 + 4.0 * x])


def L2dd(x):
    return array([[4.0] * len(x),
                  [-8.0] * len(x),
                  [4.0] * len(x)])


def L3(x):
    x2 = x * x
    x3 = x2 * x

    b0 = 0.5 * (-9.0 * x3 + 18.0 * x2 - 11.0 * x + 2.0)
    b1 = 0.5 * (27.0 * x3 - 45.0 * x2 + 18.0 * x)
    b2 = 0.5 * (-27.0 * x3 + 36.0 * x2 - 9.0 * x)
    b3 = 0.5 * (9.0 * x3 - 9.0 * x2 + 2.0 * x)

    # ~ b0 = ne.evaluate( '0.5*(-9.0*x3 + 18.0*x2 - 11.0*x + 2.0)' )
    # ~ b1 = ne.evaluate( '0.5*(27.0*x3 - 45.0*x2 + 18.0*x)' )
    # ~ b2 = ne.evaluate( '0.5*(-27.0*x3 + 36.0*x2 -  9.0*x)' )
    # ~ b3 = ne.evaluate( '0.5*(9.0*x3 -  9.0*x2 +  2.0*x)' )
    # ~ try:
    return array([b0, b1, b2, b3], dtype=float)
    # ~ except ValueError:
    # ~ pdb.set_trace()

    # ~ return 0.5 * array([  -9.0*x3 + 18.0*x2 - 11.0*x + 2.0,\
    # ~ 27.0*x3 - 45.0*x2 + 18.0*x,\
    # ~ -27.0*x3 + 36.0*x2 -  9.0*x,\
    # ~ 9.0*x3 -  9.0*x2 +  2.0*x ])


def L3d(x):
    x2 = x * x

    # ~ b0 = ne.evaluate( '0.5*(-27.0*x2 + 36.0*x - 11.0)' )
    # ~ b1 = ne.evaluate( '0.5*(81.0*x2 - 90.0*x + 18.0)' )
    # ~ b2 = ne.evaluate( '0.5*(-81.0*x2 + 72.0*x -  9.0)' )
    # ~ b3 = ne.evaluate( '0.5*(27.0*x2 - 18.0*x +  2.0)' )
    # ~
    # ~ return array([b0,b1,b2,b3])

    return 0.5 * array([-27.0 * x2 + 36.0 * x - 11.0,
                        81.0 * x2 - 90.0 * x + 18.0,
                        -81.0 * x2 + 72.0 * x - 9.0,
                        27.0 * x2 - 18.0 * x + 2.0])


def L3dd(x):
    # ~ b0 = ne.evaluate( '-27.0*x + 18.0' )
    # ~ b1 = ne.evaluate( '81.0*x - 45.0' )
    # ~ b2 = ne.evaluate( '-81.0*x + 36.0' )
    # ~ b3 = ne.evaluate( '27.0*x - 9.0' )
    # ~
    # ~ return array([b0,b1,b2,b3])
    # ~
    return array([-27.0 * x + 18.0,
                  81.0 * x - 45.0,
                  -81.0 * x + 36.0,
                  27.0 * x - 9.0])


def L4(x):
    sc = 1 / 3.
    x2 = x * x
    x3 = x2 * x
    x4 = x3 * x

    b0 = sc * (32 * x4 - 80 * x3 + 70 * x2 - 25 * x + 3)
    b1 = sc * (-128 * x4 + 288 * x3 - 208 * x2 + 48 * x)
    b2 = sc * (192 * x4 - 384 * x3 + 228 * x2 - 36 * x)
    b3 = sc * (-128 * x4 + 224 * x3 - 112 * x2 + 16 * x)
    b4 = sc * (32 * x4 - 48 * x3 + 22 * x2 - 3 * x)

    # ~ b0 = ne.evaluate( 'sc*(32*x4-80*x3+70*x2-25*x+3)' )
    # ~ b1 = ne.evaluate( 'sc*(-128*x4+288*x3-208*x2+48*x)' )
    # ~ b2 = ne.evaluate( 'sc*(192*x4-384*x3+228*x2-36*x)' )
    # ~ b3 = ne.evaluate( 'sc*(-128*x4+224*x3-112*x2+16*x)' )
    # ~ b4 = ne.evaluate( 'sc*(32*x4-48*x3+22*x2-3*x)' )

    return array([b0, b1, b2, b3, b4])

    # ~ return array([\
    # ~ sc*(32*x4-80*x3+70*x2-25*x+3),
    # ~ sc*(-128*x4+288*x3-208*x2+48*x),
    # ~ sc*(192*x4-384*x3+228*x2-36*x),
    # ~ sc*(-128*x4+224*x3-112*x2+16*x),
    # ~ sc*(32*x4-48*x3+22*x2-3*x)])


def L4d(x):
    sc = 1 / 3.
    x2 = x * x
    x3 = x2 * x

    b0 = sc * (128 * x3 - 240 * x2 + 140 * x - 25)
    b1 = sc * (-512 * x3 + 864 * x2 - 416 * x + 48)
    b2 = sc * (768 * x3 - 1152 * x2 + 456 * x - 36)
    b3 = sc * (-512 * x3 + 672 * x2 - 224 * x + 16)
    b4 = sc * (128 * x3 - 144 * x2 + 44 * x - 3)

    # ~ b0 = ne.evaluate( 'sc*(128*x3-240*x2+140*x-25)')
    # ~ b1 = ne.evaluate( 'sc*(-512*x3+864*x2-416*x+48)')
    # ~ b2 = ne.evaluate( 'sc*(768*x3-1152*x2+456*x-36)')
    # ~ b3 = ne.evaluate( 'sc*(-512*x3+672*x2-224*x+16)')
    # ~ b4 = ne.evaluate( 'sc*(128*x3-144*x2+44*x-3)')

    return array([b0, b1, b2, b3, b4])

    # ~ return array([
    # ~ sc*(128*x3-240*x2+140*x-25),
    # ~ sc*(-512*x3+864*x2-416*x+48),
    # ~ sc*(768*x3-1152*x2+456*x-36),
    # ~ sc*(-512*x3+672*x2-224*x+16),
    # ~ sc*(128*x3-144*x2+44*x-3)])


def L4dd(x):
    sc = 1 / 3.
    x2 = x * x

    b0 = sc * (384 * x2 - 480 * x + 140)
    b1 = sc * (-1536 * x2 + 1728 * x - 416)
    b2 = sc * (2304 * x2 - 2304 * x + 456)
    b3 = sc * (-1536 * x2 + 1344 * x - 224)
    b4 = sc * (384 * x2 - 288 * x + 44)

    # ~ b0 = ne.evaluate( 'sc*(384*x2-480*x+140)')
    # ~ b1 = ne.evaluate( 'sc*(-1536*x2+1728*x-416)')
    # ~ b2 = ne.evaluate( 'sc*(2304*x2-2304*x+456)')
    # ~ b3 = ne.evaluate( 'sc*(-1536*x2+1344*x-224)')
    # ~ b4 = ne.evaluate( 'sc*(384*x2-288*x+44)')

    return array([b0, b1, b2, b3, b4])

    # ~ return array([
    # ~ sc*(384*x2-480*x+140),
    # ~ sc*(-1536*x2+1728*x-416),
    # ~ sc*(2304*x2-2304*x+456),
    # ~ sc*(-1536*x2+1344*x-224),
    # ~ sc*(384*x2-288*x+44)])


# ======================================================================#
class tensor_product_basis:
    """ Tensor product basis object for quad elements
    dimension: number of dimensions
    basis_matrices: a list of tensor product basis matrices, one for each dimension
    type: a string containing the names of basis in each dimension
    """

    def __init__(self, dimensions, basis_matrices, type):

        self.dimensions = dimensions  # integer of number of dimensions
        self.basis_matrices = basis_matrices  # list of rank 2 arrays
        self.type = type  # string of basis type?
        self.basis_orders = [(b.shape[0] - 1) for b in self.basis_matrices]

    def eval(self, x):
        # evaluates the basis function at the given element coordinates
        # returns a list of values to be dot multiplied with nodal values
        # x should be a list of lists of xi coordinates
        # e.g. [ xi1 , xi2 ,...]

        # check x to be of right dimensionality
        if len(x) != self.dimensions:
            log.debug('ERROR: wrong number of xi, need', self.dimensions)
            return None
        else:
            # get phis for each dimension
            phid = []
            for i in range(self.dimensions):
                x_vector = [x[i] ** p for p in range(self.basis_orders[i] + 1)]
                phid.append(dot(self.basis_matrices[i], x_vector))

            # do tensor product using kron
            # reverse dimensions to get same implementation as cmiss
            phi = phid[0]
            if self.dimensions > 1:
                for d in range(1, self.dimensions):
                    phi = kron(phid[d], phi)

            return phi


# 1D quadratic lagrange
class line_L2(object):
    def __init__(self):
        self.dimensions = 1
        self.basis_order = [2, ]
        self.type = 'line_L2'
        self.tol = 1.0e-12

    def eval(self, x):
        return L2(x)

    def eval_derivatives(self, x, deriv=None):
        """ if deriv==None, evaluate all derivatives
        """
        get_derivatives = {(1,): self.eval_dx0,
                           (2,): self.eval_dx0x0}

        if deriv:
            try:
                d = get_derivatives[deriv](x)
            except KeyError:
                raise Warning('derivative ' + str(deriv) + ' not supported')
        else:
            d = array([get_derivatives[(1,)](x), get_derivatives[(2,)](x)])

        return d

    def eval_dx0(self, x):
        return L2d(x)

    def eval_dx0x0(self, x):
        return L2dd(x)


# 1D cubic lagrange
class line_L3(object):
    def __init__(self):
        self.dimensions = 1
        self.basis_order = [3, ]
        self.type = 'line_L3'
        self.tol = 1.0e-12

    def eval(self, x):
        return L3(x)

    def eval_derivatives(self, x, deriv=None):
        """ if deriv==None, evaluate all derivatives
        """
        get_derivatives = {(1,): self.eval_dx0,
                           (2,): self.eval_dx0x0}

        if deriv:
            try:
                d = get_derivatives[deriv](x)
            except KeyError:
                raise Warning('derivative ' + str(deriv) + ' not supported')
        else:
            d = array([get_derivatives[(1,)](x), get_derivatives[(2,)](x)])

        return d

    def eval_dx0(self, x):
        return L3d(x)

    def eval_dx0x0(self, x):
        return L3dd(x)


# 1D quartic lagrange
class line_L4(object):
    def __init__(self):
        self.dimensions = 1
        self.basis_order = [4, ]
        self.type = 'line_L4'
        self.tol = 1.0e-12

    def eval(self, x):
        return L4(x)

    def eval_derivatives(self, x, deriv=None):
        """ if deriv==None, evaluate all derivatives
        """
        get_derivatives = {(1,): self.eval_dx0,
                           (2,): self.eval_dx0x0}

        if deriv:
            try:
                d = get_derivatives[deriv](x)
            except KeyError:
                raise Warning('derivative ' + str(deriv) + ' not supported')
        else:
            d = array([get_derivatives[(1,)](x), get_derivatives[(2,)](x)])

        return d

    def eval_dx0(self, x):
        return L4d(x)

    def eval_dx0x0(self, x):
        return L4dd(x)


# 2D quad quadratic
class quad_L2_L2(object):
    def __init__(self):
        self.dimensions = 2
        self.basis_order = [2, 2]
        self.type = 'quad_L2_L2'
        self.tol = 1.0e-12

    def tensor(self, Phi1, Phi2):
        Phi = array([
            Phi1[0] * Phi2[0], Phi1[1] * Phi2[0], Phi1[2] * Phi2[0],
            Phi1[0] * Phi2[1], Phi1[1] * Phi2[1], Phi1[2] * Phi2[1],
            Phi1[0] * Phi2[2], Phi1[1] * Phi2[2], Phi1[2] * Phi2[2]])

        return where(abs(Phi) < self.tol, 0.0, Phi)

    def eval(self, x):
        return self.tensor(L2(x[0]), L2(x[1]))

    def eval_derivatives(self, x, deriv=None):
        """ if deriv==None, evaluate all derivatives
        """
        get_derivatives = {(1, 0): self.eval_dx0,
                           (0, 1): self.eval_dx1,
                           (2, 0): self.eval_dx0x0,
                           (0, 2): self.eval_dx1x1,
                           (1, 1): self.eval_dx0x1}

        if deriv:
            try:
                d = get_derivatives[deriv](x)
            except KeyError:
                raise Warning('derivative ' + str(deriv) + ' not supported')
        else:
            d = array([get_derivatives[(1, 0)](x), get_derivatives[(0, 1)](x), get_derivatives[(2, 0)](x),
                       get_derivatives[(0, 2)](x), get_derivatives[(1, 1)](x)])

        return d

    def eval_dx0(self, x):
        return self.tensor(L2d(x[0]), L2(x[1]))

    def eval_dx1(self, x):
        return self.tensor(L2(x[0]), L2d(x[1]))

    def eval_dx0x0(self, x):
        return self.tensor(L2dd(x[0]), L2(x[1]))

    def eval_dx1x1(self, x):
        return self.tensor(L2(x[0]), L2dd(x[1]))

    def eval_dx0x1(self, x):
        return self.tensor(L2d(x[0]), L2d(x[1]))


# 2D quad cubic
class quad_L3_L3(object):
    def __init__(self):
        self.dimensions = 2
        self.basis_order = [3, 3]
        self.type = 'quad_L3_L3'
        self.tol = 1.0e-12

    def tensor(self, Phi1, Phi2):
        Phi = array([
            Phi1[0] * Phi2[0], Phi1[1] * Phi2[0], Phi1[2] * Phi2[0], Phi1[3] * Phi2[0],
            Phi1[0] * Phi2[1], Phi1[1] * Phi2[1], Phi1[2] * Phi2[1], Phi1[3] * Phi2[1],
            Phi1[0] * Phi2[2], Phi1[1] * Phi2[2], Phi1[2] * Phi2[2], Phi1[3] * Phi2[2],
            Phi1[0] * Phi2[3], Phi1[1] * Phi2[3], Phi1[2] * Phi2[3], Phi1[3] * Phi2[3]])

        return where(abs(Phi) < self.tol, 0.0, Phi)

    def eval(self, x):
        return self.tensor(L3(x[0]), L3(x[1]))

    def eval_derivatives(self, x, deriv=None):
        """ if deriv==None, evaluate all derivatives
        """
        get_derivatives = {(1, 0): self.eval_dx0,
                           (0, 1): self.eval_dx1,
                           (2, 0): self.eval_dx0x0,
                           (0, 2): self.eval_dx1x1,
                           (1, 1): self.eval_dx0x1}

        if deriv:
            try:
                d = get_derivatives[deriv](x)
            except KeyError:
                raise Warning('derivative ' + str(deriv) + ' not supported')
        else:
            d = array([get_derivatives[(1, 0)](x), get_derivatives[(0, 1)](x), get_derivatives[(2, 0)](x),
                       get_derivatives[(0, 2)](x), get_derivatives[(1, 1)](x)])

        return d

    def eval_dx0(self, x):
        return self.tensor(L3d(x[0]), L3(x[1]))

    def eval_dx1(self, x):
        return self.tensor(L3(x[0]), L3d(x[1]))

    def eval_dx0x0(self, x):
        return self.tensor(L3dd(x[0]), L3(x[1]))

    def eval_dx1x1(self, x):
        return self.tensor(L3(x[0]), L3dd(x[1]))

    def eval_dx0x1(self, x):
        return self.tensor(L3d(x[0]), L3d(x[1]))


# 3D quad quadratic
class quad_L2_L2_L2(object):
    def __init__(self):
        self.dimensions = 3
        self.basis_order = [2, 2, 2]
        self.type = 'quad_L2_L2_L2'
        self.tol = 1.0e-12

    def tensor(self, p1, p2, p3):
        p = array([
            p1[0] * p2[0] * p3[0], p1[1] * p2[0] * p3[0], p1[2] * p2[0] * p3[0],
            p1[0] * p2[1] * p3[0], p1[1] * p2[1] * p3[0], p1[2] * p2[1] * p3[0],
            p1[0] * p2[2] * p3[0], p1[1] * p2[2] * p3[0], p1[2] * p2[2] * p3[0],

            p1[0] * p2[0] * p3[1], p1[1] * p2[0] * p3[1], p1[2] * p2[0] * p3[1],
            p1[0] * p2[1] * p3[1], p1[1] * p2[1] * p3[1], p1[2] * p2[1] * p3[1],
            p1[0] * p2[2] * p3[1], p1[1] * p2[2] * p3[1], p1[2] * p2[2] * p3[1],

            p1[0] * p2[0] * p3[2], p1[1] * p2[0] * p3[2], p1[2] * p2[0] * p3[2],
            p1[0] * p2[1] * p3[2], p1[1] * p2[1] * p3[2], p1[2] * p2[1] * p3[2],
            p1[0] * p2[2] * p3[2], p1[1] * p2[2] * p3[2], p1[2] * p2[2] * p3[2],
            ])

        return where(abs(p) < self.tol, 0.0, p)

    def eval(self, x):
        return self.tensor(L2(x[0]), L2(x[1]), L2(x[2]))

    def eval_derivatives(self, x, deriv=None):
        """ if deriv==None, evaluate all derivatives
        """
        ### no cross derivatives yet
        get_derivatives = {(1, 0, 0): self.eval_dx0,
                           (0, 1, 0): self.eval_dx1,
                           (0, 0, 1): self.eval_dx2,
                           (2, 0, 0): self.eval_dx0x0,
                           (0, 2, 0): self.eval_dx1x1,
                           (0, 0, 2): self.eval_dx2x2,
                           (1, 1, 0): self.eval_dx0x1,
                           (1, 0, 1): self.eval_dx0x2,
                           (0, 1, 1): self.eval_dx1x2,
                           (1, 1, 1): self.eval_dx0x1x2
                           }

        if deriv:
            try:
                d = get_derivatives[deriv](x)
            except KeyError:
                raise Warning('derivative ' + str(deriv) + ' not supported')
        else:
            d = array([get_derivatives[(1, 0, 0)](x), get_derivatives[(0, 1, 0)](x), get_derivatives[(0, 0, 1)](x),
                       get_derivatives[(2, 0, 0)](x), get_derivatives[(0, 2, 0)](x), get_derivatives[(0, 0, 2)](x),
                       get_derivatives[(1, 1, 0)](x), get_derivatives[(1, 0, 1)](x), get_derivatives[(0, 1, 1)](x),
                       get_derivatives[(1, 1, 1)](x),
                       ])

        return d

    def eval_dx0(self, x):
        return self.tensor(L2d(x[0]), L2(x[1]), L2(x[2]))

    def eval_dx1(self, x):
        return self.tensor(L2(x[0]), L2d(x[1]), L2(x[2]))

    def eval_dx2(self, x):
        return self.tensor(L2(x[0]), L2(x[1]), L2d(x[2]))

    def eval_dx0x0(self, x):
        return self.tensor(L2dd(x[0]), L2(x[1]), L2(x[2]))

    def eval_dx1x1(self, x):
        return self.tensor(L2(x[0]), L2dd(x[1]), L2(x[2]))

    def eval_dx2x2(self, x):
        return self.tensor(L2(x[0]), L2(x[1]), L2dd(x[2]))

    def eval_dx0x1(self, x):
        return self.tensor(L2d(x[0]), L2d(x[1]), L2(x[2]))

    def eval_dx0x2(self, x):
        return self.tensor(L2d(x[0]), L2(x[1]), L2d(x[2]))

    def eval_dx1x2(self, x):
        return self.tensor(L2(x[0]), L2d(x[1]), L2d(x[2]))

    def eval_dx0x1x2(self, x):
        return self.tensor(L2d(x[0]), L2d(x[1]), L2d(x[2]))


# 3D cube cubic element
class quad_L3_L3_L3(object):
    def __init__(self):
        self.dimensions = 3
        self.basis_order = [3, 3, 3]
        self.type = 'quad_L3_L3_L3'
        self.tol = 1.0e-12

    def tensor(self, p1, p2, p3):
        p = array([
            p1[0] * p2[0] * p3[0], p1[1] * p2[0] * p3[0], p1[2] * p2[0] * p3[0], p1[3] * p2[0] * p3[0],
            p1[0] * p2[1] * p3[0], p1[1] * p2[1] * p3[0], p1[2] * p2[1] * p3[0], p1[3] * p2[1] * p3[0],
            p1[0] * p2[2] * p3[0], p1[1] * p2[2] * p3[0], p1[2] * p2[2] * p3[0], p1[3] * p2[2] * p3[0],
            p1[0] * p2[3] * p3[0], p1[1] * p2[3] * p3[0], p1[2] * p2[3] * p3[0], p1[3] * p2[3] * p3[0],

            p1[0] * p2[0] * p3[1], p1[1] * p2[0] * p3[1], p1[2] * p2[0] * p3[1], p1[3] * p2[0] * p3[1],
            p1[0] * p2[1] * p3[1], p1[1] * p2[1] * p3[1], p1[2] * p2[1] * p3[1], p1[3] * p2[1] * p3[1],
            p1[0] * p2[2] * p3[1], p1[1] * p2[2] * p3[1], p1[2] * p2[2] * p3[1], p1[3] * p2[2] * p3[1],
            p1[0] * p2[3] * p3[1], p1[1] * p2[3] * p3[1], p1[2] * p2[3] * p3[1], p1[3] * p2[3] * p3[1],

            p1[0] * p2[0] * p3[2], p1[1] * p2[0] * p3[2], p1[2] * p2[0] * p3[2], p1[3] * p2[0] * p3[2],
            p1[0] * p2[1] * p3[2], p1[1] * p2[1] * p3[2], p1[2] * p2[1] * p3[2], p1[3] * p2[1] * p3[2],
            p1[0] * p2[2] * p3[2], p1[1] * p2[2] * p3[2], p1[2] * p2[2] * p3[2], p1[3] * p2[2] * p3[2],
            p1[0] * p2[3] * p3[2], p1[1] * p2[3] * p3[2], p1[2] * p2[3] * p3[2], p1[3] * p2[3] * p3[2],

            p1[0] * p2[0] * p3[3], p1[1] * p2[0] * p3[3], p1[2] * p2[0] * p3[3], p1[3] * p2[0] * p3[3],
            p1[0] * p2[1] * p3[3], p1[1] * p2[1] * p3[3], p1[2] * p2[1] * p3[3], p1[3] * p2[1] * p3[3],
            p1[0] * p2[2] * p3[3], p1[1] * p2[2] * p3[3], p1[2] * p2[2] * p3[3], p1[3] * p2[2] * p3[3],
            p1[0] * p2[3] * p3[3], p1[1] * p2[3] * p3[3], p1[2] * p2[3] * p3[3], p1[3] * p2[3] * p3[3],
            ])

        # ~ try:
        return where(abs(p) < self.tol, 0.0, p)
        # ~ except ValueError:
        # ~ pdb.set_trace()

    def eval(self, x):
        return self.tensor(L3(x[0]), L3(x[1]), L3(x[2]))

    def eval_derivatives(self, x, deriv=None):
        """ if deriv==None, evaluate all derivatives
        """
        ### no cross derivatives yet
        get_derivatives = {(1, 0, 0): self.eval_dx0,
                           (0, 1, 0): self.eval_dx1,
                           (0, 0, 1): self.eval_dx2,
                           (2, 0, 0): self.eval_dx0x0,
                           (0, 2, 0): self.eval_dx1x1,
                           (0, 0, 2): self.eval_dx2x2,
                           (1, 1, 0): self.eval_dx0x1,
                           (1, 0, 1): self.eval_dx0x2,
                           (0, 1, 1): self.eval_dx1x2,
                           (1, 1, 1): self.eval_dx0x1x2
                           }

        if deriv:
            try:
                d = get_derivatives[deriv](x)
            except KeyError:
                raise Warning('derivative ' + str(deriv) + ' not supported')
        else:
            d = array([get_derivatives[(1, 0, 0)](x), get_derivatives[(0, 1, 0)](x), get_derivatives[(0, 0, 1)](x),
                       get_derivatives[(2, 0, 0)](x), get_derivatives[(0, 2, 0)](x), get_derivatives[(0, 0, 2)](x),
                       get_derivatives[(1, 1, 0)](x), get_derivatives[(1, 0, 1)](x), get_derivatives[(0, 1, 1)](x),
                       get_derivatives[(1, 1, 1)](x),
                       ])

        return d

    def eval_dx0(self, x):
        return self.tensor(L3d(x[0]), L3(x[1]), L3(x[2]))

    def eval_dx1(self, x):
        return self.tensor(L3(x[0]), L3d(x[1]), L3(x[2]))

    def eval_dx2(self, x):
        return self.tensor(L3(x[0]), L3(x[1]), L3d(x[2]))

    def eval_dx0x0(self, x):
        return self.tensor(L3dd(x[0]), L3(x[1]), L3(x[2]))

    def eval_dx1x1(self, x):
        return self.tensor(L3(x[0]), L3dd(x[1]), L3(x[2]))

    def eval_dx2x2(self, x):
        return self.tensor(L3(x[0]), L3d(x[1]), L3dd(x[2]))

    def eval_dx0x1(self, x):
        return self.tensor(L3d(x[0]), L3d(x[1]), L3(x[2]))

    def eval_dx0x2(self, x):
        return self.tensor(L3d(x[0]), L3(x[1]), L3d(x[2]))

    def eval_dx1x2(self, x):
        return self.tensor(L3(x[0]), L3d(x[1]), L3d(x[2]))

    def eval_dx0x1x2(self, x):
        return self.tensor(L3d(x[0]), L3d(x[1]), L3d(x[2]))


# 2D quad quartic element
class quad_L4_L4(object):

    def __init__(self):
        self.dimensions = 2
        self.basis_order = [4, 4]
        self.type = 'quad_L4_L4'
        self.tol = 1.0e-12

    def tensor(self, p1, p2):
        p = array([
            p1[0] * p2[0], p1[1] * p2[0], p1[2] * p2[0], p1[3] * p2[0], p1[4] * p2[0],
            p1[0] * p2[1], p1[1] * p2[1], p1[2] * p2[1], p1[3] * p2[1], p1[4] * p2[1],
            p1[0] * p2[2], p1[1] * p2[2], p1[2] * p2[2], p1[3] * p2[2], p1[4] * p2[2],
            p1[0] * p2[3], p1[1] * p2[3], p1[2] * p2[3], p1[3] * p2[3], p1[4] * p2[3],
            p1[0] * p2[4], p1[1] * p2[4], p1[2] * p2[4], p1[3] * p2[4], p1[4] * p2[4],
            ])
        return where(abs(p) < self.tol, 0.0, p)

    def eval(self, x):
        return self.tensor(L4(x[0]), L4(x[1]))

    def eval_derivatives(self, x, deriv=None):
        """ if deriv==None, evaluate all derivatives
        """
        ### no cross derivatives yet
        get_derivatives = {(1, 0): self.eval_dx0,
                           (0, 1): self.eval_dx1,
                           (2, 0): self.eval_dx0x0,
                           (0, 2): self.eval_dx1x1,
                           (1, 1): self.eval_dx0x1,
                           }

        if deriv:
            try:
                d = get_derivatives[deriv](x)
            except KeyError:
                raise Warning('derivative ' + str(deriv) + ' not supported')
        else:
            x00 = L4(x[0])
            x01 = L4d(x[0])
            x02 = L4dd(x[0])
            x10 = L4(x[1])
            x11 = L4d(x[1])
            x12 = L4dd(x[1])
            d = array([self.tensor(x01, x10), self.tensor(x00, x11),
                       self.tensor(x02, x10), self.tensor(x00, x12),
                       self.tensor(x01, x11)])

            # ~ d = array( [ get_derivatives[(1,0)](x),get_derivatives[(0,1)](x),\
            # ~ get_derivatives[(2,0)](x),get_derivatives[(0,2)](x),get_derivatives[(1,1)](x) ] )

        return d

    def eval_dx0(self, x):
        return self.tensor(L4d(x[0]), L4(x[1]))

    def eval_dx1(self, x):
        return self.tensor(L4(x[0]), L4d(x[1]))

    def eval_dx0x0(self, x):
        return self.tensor(L4dd(x[0]), L4(x[1]))

    def eval_dx1x1(self, x):
        return self.tensor(L4(x[0]), L4dd(x[1]))

    def eval_dx0x1(self, x):
        return self.tensor(L4d(x[0]), L4d(x[1]))


# 2D quad quartic cubic element
class quad_L4_L3(object):

    def __init__(self):
        self.dimensions = 2
        self.basis_order = [4, 3]
        self.type = 'quad_L4_L3'
        self.tol = 1.0e-12

    def tensor(self, p1, p2):
        p = array([
            p1[0] * p2[0], p1[1] * p2[0], p1[2] * p2[0], p1[3] * p2[0], p1[4] * p2[0],
            p1[0] * p2[1], p1[1] * p2[1], p1[2] * p2[1], p1[3] * p2[1], p1[4] * p2[1],
            p1[0] * p2[2], p1[1] * p2[2], p1[2] * p2[2], p1[3] * p2[2], p1[4] * p2[2],
            p1[0] * p2[3], p1[1] * p2[3], p1[2] * p2[3], p1[3] * p2[3], p1[4] * p2[3],
            ])
        return where(abs(p) < self.tol, 0.0, p)

    def eval(self, x):
        return self.tensor(L4(x[0]), L3(x[1]))

    def eval_derivatives(self, x, deriv=None):
        """ if deriv==None, evaluate all derivatives
        """
        ### no cross derivatives yet
        get_derivatives = {(1, 0): self.eval_dx0,
                           (0, 1): self.eval_dx1,
                           (2, 0): self.eval_dx0x0,
                           (0, 2): self.eval_dx1x1,
                           (1, 1): self.eval_dx0x1,
                           }

        if deriv:
            try:
                d = get_derivatives[deriv](x)
            except KeyError:
                raise Warning('derivative ' + str(deriv) + ' not supported')
        else:
            x00 = L4(x[0])
            x01 = L4d(x[0])
            x02 = L4dd(x[0])
            x10 = L3(x[1])
            x11 = L3d(x[1])
            x12 = L3dd(x[1])
            d = array([self.tensor(x01, x10), self.tensor(x00, x11),
                       self.tensor(x02, x10), self.tensor(x00, x12),
                       self.tensor(x01, x11)])

            # ~ d = array( [ get_derivatives[(1,0)](x),get_derivatives[(0,1)](x),\
            # ~ get_derivatives[(2,0)](x),get_derivatives[(0,2)](x),get_derivatives[(1,1)](x) ] )

        return d

    def eval_dx0(self, x):
        return self.tensor(L4d(x[0]), L3(x[1]))

    def eval_dx1(self, x):
        return self.tensor(L4(x[0]), L3d(x[1]))

    def eval_dx0x0(self, x):
        return self.tensor(L4dd(x[0]), L3(x[1]))

    def eval_dx1x1(self, x):
        return self.tensor(L4(x[0]), L3dd(x[1]))

    def eval_dx0x1(self, x):
        return self.tensor(L4d(x[0]), L3d(x[1]))


# 3D cube quartic element
class quad_L4_L4_L4(object):

    def __init__(self):
        self.dimensions = 3
        self.basis_order = [4, 4, 4]
        self.type = 'quad_L4_L4_L4'
        self.tol = 1.0e-12

    def tensor(self, p1, p2, p3):
        p = array([
            p1[0] * p2[0] * p3[0], p1[1] * p2[0] * p3[0], p1[2] * p2[0] * p3[0], p1[3] * p2[0] * p3[0],
            p1[4] * p2[0] * p3[0],
            p1[0] * p2[1] * p3[0], p1[1] * p2[1] * p3[0], p1[2] * p2[1] * p3[0], p1[3] * p2[1] * p3[0],
            p1[4] * p2[1] * p3[0],
            p1[0] * p2[2] * p3[0], p1[1] * p2[2] * p3[0], p1[2] * p2[2] * p3[0], p1[3] * p2[2] * p3[0],
            p1[4] * p2[2] * p3[0],
            p1[0] * p2[3] * p3[0], p1[1] * p2[3] * p3[0], p1[2] * p2[3] * p3[0], p1[3] * p2[3] * p3[0],
            p1[4] * p2[3] * p3[0],
            p1[0] * p2[4] * p3[0], p1[1] * p2[4] * p3[0], p1[2] * p2[4] * p3[0], p1[3] * p2[4] * p3[0],
            p1[4] * p2[4] * p3[0],

            p1[0] * p2[0] * p3[1], p1[1] * p2[0] * p3[1], p1[2] * p2[0] * p3[1], p1[3] * p2[0] * p3[1],
            p1[4] * p2[0] * p3[1],
            p1[0] * p2[1] * p3[1], p1[1] * p2[1] * p3[1], p1[2] * p2[1] * p3[1], p1[3] * p2[1] * p3[1],
            p1[4] * p2[1] * p3[1],
            p1[0] * p2[2] * p3[1], p1[1] * p2[2] * p3[1], p1[2] * p2[2] * p3[1], p1[3] * p2[2] * p3[1],
            p1[4] * p2[2] * p3[1],
            p1[0] * p2[3] * p3[1], p1[1] * p2[3] * p3[1], p1[2] * p2[3] * p3[1], p1[3] * p2[3] * p3[1],
            p1[4] * p2[3] * p3[1],
            p1[0] * p2[4] * p3[1], p1[1] * p2[4] * p3[1], p1[2] * p2[4] * p3[1], p1[3] * p2[4] * p3[1],
            p1[4] * p2[4] * p3[1],

            p1[0] * p2[0] * p3[2], p1[1] * p2[0] * p3[2], p1[2] * p2[0] * p3[2], p1[3] * p2[0] * p3[2],
            p1[4] * p2[0] * p3[2],
            p1[0] * p2[1] * p3[2], p1[1] * p2[1] * p3[2], p1[2] * p2[1] * p3[2], p1[3] * p2[1] * p3[2],
            p1[4] * p2[1] * p3[2],
            p1[0] * p2[2] * p3[2], p1[1] * p2[2] * p3[2], p1[2] * p2[2] * p3[2], p1[3] * p2[2] * p3[2],
            p1[4] * p2[2] * p3[2],
            p1[0] * p2[3] * p3[2], p1[1] * p2[3] * p3[2], p1[2] * p2[3] * p3[2], p1[3] * p2[3] * p3[2],
            p1[4] * p2[3] * p3[2],
            p1[0] * p2[4] * p3[2], p1[1] * p2[4] * p3[2], p1[2] * p2[4] * p3[2], p1[3] * p2[4] * p3[2],
            p1[4] * p2[4] * p3[2],

            p1[0] * p2[0] * p3[3], p1[1] * p2[0] * p3[3], p1[2] * p2[0] * p3[3], p1[3] * p2[0] * p3[3],
            p1[4] * p2[0] * p3[3],
            p1[0] * p2[1] * p3[3], p1[1] * p2[1] * p3[3], p1[2] * p2[1] * p3[3], p1[3] * p2[1] * p3[3],
            p1[4] * p2[1] * p3[3],
            p1[0] * p2[2] * p3[3], p1[1] * p2[2] * p3[3], p1[2] * p2[2] * p3[3], p1[3] * p2[2] * p3[3],
            p1[4] * p2[2] * p3[3],
            p1[0] * p2[3] * p3[3], p1[1] * p2[3] * p3[3], p1[2] * p2[3] * p3[3], p1[3] * p2[3] * p3[3],
            p1[4] * p2[3] * p3[3],
            p1[0] * p2[4] * p3[3], p1[1] * p2[4] * p3[3], p1[2] * p2[4] * p3[3], p1[3] * p2[4] * p3[3],
            p1[4] * p2[4] * p3[3],

            p1[0] * p2[0] * p3[4], p1[1] * p2[0] * p3[4], p1[2] * p2[0] * p3[4], p1[3] * p2[0] * p3[4],
            p1[4] * p2[0] * p3[4],
            p1[0] * p2[1] * p3[4], p1[1] * p2[1] * p3[4], p1[2] * p2[1] * p3[4], p1[3] * p2[1] * p3[4],
            p1[4] * p2[1] * p3[4],
            p1[0] * p2[2] * p3[4], p1[1] * p2[2] * p3[4], p1[2] * p2[2] * p3[4], p1[3] * p2[2] * p3[4],
            p1[4] * p2[2] * p3[4],
            p1[0] * p2[3] * p3[4], p1[1] * p2[3] * p3[4], p1[2] * p2[3] * p3[4], p1[3] * p2[3] * p3[4],
            p1[4] * p2[3] * p3[4],
            p1[0] * p2[4] * p3[4], p1[1] * p2[4] * p3[4], p1[2] * p2[4] * p3[4], p1[3] * p2[4] * p3[4],
            p1[4] * p2[4] * p3[4],
            ])

        return where(abs(p) < self.tol, 0.0, p)

    def eval(self, x):
        return self.tensor(L4(x[0]), L4(x[1]), L4(x[2]))

    def eval_derivatives(self, x, deriv=None):
        """ if deriv==None, evaluate all derivatives
        """
        ### no cross derivatives yet
        get_derivatives = {(1, 0, 0): self.eval_dx0,
                           (0, 1, 0): self.eval_dx1,
                           (0, 0, 1): self.eval_dx2,
                           (2, 0, 0): self.eval_dx0x0,
                           (0, 2, 0): self.eval_dx1x1,
                           (0, 0, 2): self.eval_dx2x2,
                           (1, 1, 0): self.eval_dx0x1,
                           (1, 0, 1): self.eval_dx0x2,
                           (0, 1, 1): self.eval_dx1x2,
                           (1, 1, 1): self.eval_dx0x1x2
                           }

        if deriv:
            try:
                d = get_derivatives[deriv](x)
            except KeyError:
                raise Warning('derivative ' + str(deriv) + ' not supported')
        else:
            d = array([get_derivatives[(1, 0, 0)](x), get_derivatives[(0, 1, 0)](x), get_derivatives[(0, 0, 1)](x),
                       get_derivatives[(2, 0, 0)](x), get_derivatives[(0, 2, 0)](x), get_derivatives[(0, 0, 2)](x),
                       get_derivatives[(1, 1, 0)](x), get_derivatives[(1, 0, 1)](x), get_derivatives[(0, 1, 1)](x),
                       get_derivatives[(1, 1, 1)](x),
                       ])

        return d

    def eval_dx0(self, x):
        return self.tensor(L4d(x[0]), L4(x[1]), L4(x[2]))

    def eval_dx1(self, x):
        return self.tensor(L4(x[0]), L4d(x[1]), L4(x[2]))

    def eval_dx2(self, x):
        return self.tensor(L4(x[0]), L4(x[1]), L4d(x[2]))

    def eval_dx0x0(self, x):
        return self.tensor(L4dd(x[0]), L4(x[1]), L4(x[2]))

    def eval_dx1x1(self, x):
        return self.tensor(L4(x[0]), L4dd(x[1]), L4(x[2]))

    def eval_dx2x2(self, x):
        return self.tensor(L4(x[0]), L4(x[1]), L4dd(x[2]))

    def eval_dx0x1(self, x):
        return self.tensor(L4d(x[0]), L4d(x[1]), L4(x[2]))

    def eval_dx0x2(self, x):
        return self.tensor(L4d(x[0]), L4(x[1]), L4d(x[2]))

    def eval_dx1x2(self, x):
        return self.tensor(L4(x[0]), L4d(x[1]), L4d(x[2]))

    def eval_dx0x1x2(self, x):
        return self.tensor(L4d(x[0]), L4d(x[1]), L4d(x[2]))


class quad_L4_L3_L1(object):

    def __init__(self):
        self.dimensions = 3
        self.basis_order = [4, 3, 1]
        self.type = 'quad_L4_L3_L1'
        self.tol = 1.0e-12

    def tensor(self, p1, p2, p3):
        p = array([
            p1[0] * p2[0] * p3[0], p1[1] * p2[0] * p3[0], p1[2] * p2[0] * p3[0], p1[3] * p2[0] * p3[0],
            p1[4] * p2[0] * p3[0],
            p1[0] * p2[1] * p3[0], p1[1] * p2[1] * p3[0], p1[2] * p2[1] * p3[0], p1[3] * p2[1] * p3[0],
            p1[4] * p2[1] * p3[0],
            p1[0] * p2[2] * p3[0], p1[1] * p2[2] * p3[0], p1[2] * p2[2] * p3[0], p1[3] * p2[2] * p3[0],
            p1[4] * p2[2] * p3[0],
            p1[0] * p2[3] * p3[0], p1[1] * p2[3] * p3[0], p1[2] * p2[3] * p3[0], p1[3] * p2[3] * p3[0],
            p1[4] * p2[3] * p3[0],

            p1[0] * p2[0] * p3[1], p1[1] * p2[0] * p3[1], p1[2] * p2[0] * p3[1], p1[3] * p2[0] * p3[1],
            p1[4] * p2[0] * p3[1],
            p1[0] * p2[1] * p3[1], p1[1] * p2[1] * p3[1], p1[2] * p2[1] * p3[1], p1[3] * p2[1] * p3[1],
            p1[4] * p2[1] * p3[1],
            p1[0] * p2[2] * p3[1], p1[1] * p2[2] * p3[1], p1[2] * p2[2] * p3[1], p1[3] * p2[2] * p3[1],
            p1[4] * p2[2] * p3[1],
            p1[0] * p2[3] * p3[1], p1[1] * p2[3] * p3[1], p1[2] * p2[3] * p3[1], p1[3] * p2[3] * p3[1],
            p1[4] * p2[3] * p3[1],
            ])

        return where(abs(p) < self.tol, 0.0, p)

    def eval(self, x):
        return self.tensor(L4(x[0]), L3(x[1]), L1(x[2]))

    def eval_derivatives(self, x, deriv=None):
        """ if deriv==None, evaluate all derivatives
        """
        ### no cross derivatives yet
        get_derivatives = {(1, 0, 0): self.eval_dx0,
                           (0, 1, 0): self.eval_dx1,
                           (0, 0, 1): self.eval_dx2,
                           (2, 0, 0): self.eval_dx0x0,
                           (0, 2, 0): self.eval_dx1x1,
                           (0, 0, 2): self.eval_dx2x2,
                           (1, 1, 0): self.eval_dx0x1,
                           (1, 0, 1): self.eval_dx0x2,
                           (0, 1, 1): self.eval_dx1x2,
                           (1, 1, 1): self.eval_dx0x1x2
                           }

        if deriv:
            try:
                d = get_derivatives[deriv](x)
            except KeyError:
                raise Warning('derivative ' + str(deriv) + ' not supported')
        else:
            d = array([get_derivatives[(1, 0, 0)](x), get_derivatives[(0, 1, 0)](x), get_derivatives[(0, 0, 1)](x),
                       get_derivatives[(2, 0, 0)](x), get_derivatives[(0, 2, 0)](x), get_derivatives[(0, 0, 2)](x),
                       get_derivatives[(1, 1, 0)](x), get_derivatives[(1, 0, 1)](x), get_derivatives[(0, 1, 1)](x),
                       get_derivatives[(1, 1, 1)](x),
                       ])

    def eval_dx0(self, x):
        return self.tensor(L4d(x[0]), L3(x[1]), L1(x[2]))

    def eval_dx1(self, x):
        return self.tensor(L4(x[0]), L3d(x[1]), L1(x[2]))

    def eval_dx2(self, x):
        return self.tensor(L4(x[0]), L3(x[1]), L1d(x[2]))

    def eval_dx0x0(self, x):
        return self.tensor(L4dd(x[0]), L3(x[1]), L1(x[2]))

    def eval_dx1x2(self, x):
        return self.tensor(L4(x[0]), L3dd(x[1]), L1(x[2]))

    def eval_dx2x2(self, x):
        return self.tensor(L4(x[0]), L3d(x[1]), L1dd(x[2]))

    def eval_dx0x1(self, x):
        return self.tensor(L4d(x[0]), L3d(x[1]), L1(x[2]))

    def eval_dx0x2(self, x):
        return self.tensor(L4d(x[0]), L3(x[1]), L1d(x[2]))

    def eval_dx1x2(self, x):
        return self.tensor(L4(x[0]), L3d(x[1]), L1d(x[2]))

    def eval_dx0x1x2(self, x):
        return self.tensor(L4d(x[0]), L3d(x[1]), L1d(x[2]))


class quad_L4_L4_L1(object):

    def __init__(self):
        self.dimensions = 3
        self.basis_order = [4, 4, 1]
        self.type = 'quad_L4_L4_L1'
        self.tol = 1.0e-12

    def tensor(self, p1, p2, p3):
        p = array([
            p1[0] * p2[0] * p3[0], p1[1] * p2[0] * p3[0], p1[2] * p2[0] * p3[0], p1[3] * p2[0] * p3[0],
            p1[4] * p2[0] * p3[0],
            p1[0] * p2[1] * p3[0], p1[1] * p2[1] * p3[0], p1[2] * p2[1] * p3[0], p1[3] * p2[1] * p3[0],
            p1[4] * p2[1] * p3[0],
            p1[0] * p2[2] * p3[0], p1[1] * p2[2] * p3[0], p1[2] * p2[2] * p3[0], p1[3] * p2[2] * p3[0],
            p1[4] * p2[2] * p3[0],
            p1[0] * p2[3] * p3[0], p1[1] * p2[3] * p3[0], p1[2] * p2[3] * p3[0], p1[3] * p2[3] * p3[0],
            p1[4] * p2[3] * p3[0],
            p1[0] * p2[4] * p3[0], p1[1] * p2[4] * p3[0], p1[2] * p2[4] * p3[0], p1[3] * p2[4] * p3[0],
            p1[4] * p2[4] * p3[0],

            p1[0] * p2[0] * p3[1], p1[1] * p2[0] * p3[1], p1[2] * p2[0] * p3[1], p1[3] * p2[0] * p3[1],
            p1[4] * p2[0] * p3[1],
            p1[0] * p2[1] * p3[1], p1[1] * p2[1] * p3[1], p1[2] * p2[1] * p3[1], p1[3] * p2[1] * p3[1],
            p1[4] * p2[1] * p3[1],
            p1[0] * p2[2] * p3[1], p1[1] * p2[2] * p3[1], p1[2] * p2[2] * p3[1], p1[3] * p2[2] * p3[1],
            p1[4] * p2[2] * p3[1],
            p1[0] * p2[3] * p3[1], p1[1] * p2[3] * p3[1], p1[2] * p2[3] * p3[1], p1[3] * p2[3] * p3[1],
            p1[4] * p2[3] * p3[1],
            p1[0] * p2[4] * p3[1], p1[1] * p2[4] * p3[1], p1[2] * p2[4] * p3[1], p1[3] * p2[4] * p3[1],
            p1[4] * p2[4] * p3[1],
            ])

        return where(abs(p) < self.tol, 0.0, p)

    def eval(self, x):
        return self.tensor(L4(x[0]), L4(x[1]), L1(x[2]))

    def eval_derivatives(self, x, deriv=None):
        """ if deriv==None, evaluate all derivatives
        """
        ### no cross derivatives yet
        get_derivatives = {(1, 0, 0): self.eval_dx0,
                           (0, 1, 0): self.eval_dx1,
                           (0, 0, 1): self.eval_dx2,
                           (2, 0, 0): self.eval_dx0x0,
                           (0, 2, 0): self.eval_dx1x1,
                           (0, 0, 2): self.eval_dx2x2,
                           (1, 1, 0): self.eval_dx0x1,
                           (1, 0, 1): self.eval_dx0x2,
                           (0, 1, 1): self.eval_dx1x2,
                           (1, 1, 1): self.eval_dx0x1x2
                           }

        if deriv:
            try:
                d = get_derivatives[deriv](x)
            except KeyError:
                raise Warning('derivative ' + str(deriv) + ' not supported')
        else:
            d = array([get_derivatives[(1, 0, 0)](x), get_derivatives[(0, 1, 0)](x), get_derivatives[(0, 0, 1)](x),
                       get_derivatives[(2, 0, 0)](x), get_derivatives[(0, 2, 0)](x), get_derivatives[(0, 0, 2)](x),
                       get_derivatives[(1, 1, 0)](x), get_derivatives[(1, 0, 1)](x), get_derivatives[(0, 1, 1)](x),
                       get_derivatives[(1, 1, 1)](x),
                       ])

    def eval_dx0(self, x):
        return self.tensor(L4d(x[0]), L4(x[1]), L1(x[2]))

    def eval_dx1(self, x):
        return self.tensor(L4(x[0]), L4d(x[1]), L1(x[2]))

    def eval_dx2(self, x):
        return self.tensor(L4(x[0]), L4(x[1]), L1d(x[2]))

    def eval_dx0x0(self, x):
        return self.tensor(L4dd(x[0]), L4(x[1]), L1(x[2]))

    def eval_dx1x2(self, x):
        return self.tensor(L4(x[0]), L4dd(x[1]), L1(x[2]))

    def eval_dx2x2(self, x):
        return self.tensor(L4(x[0]), L4d(x[1]), L1dd(x[2]))

    def eval_dx0x1(self, x):
        return self.tensor(L4d(x[0]), L4d(x[1]), L1(x[2]))

    def eval_dx0x2(self, x):
        return self.tensor(L4d(x[0]), L4(x[1]), L1d(x[2]))

    def eval_dx1x2(self, x):
        return self.tensor(L4(x[0]), L4d(x[1]), L1d(x[2]))

    def eval_dx0x1x2(self, x):
        return self.tensor(L4d(x[0]), L4d(x[1]), L1d(x[2]))


# ======================================================================#
# simplex elements                                                     #
# ======================================================================#
class simplex_2d(object):
    """ parent class for simplex 2D basis
    """
    # element vertices
    v0 = [0.0, 0.0]
    v1 = [1.0, 0.0]
    v2 = [0.0, 1.0]
    D = 0.5 * det([[1.0, v0[0], v0[1]],
                   [1.0, v1[0], v1[1]],
                   [1.0, v2[0], v2[1]]])
    D2 = 2.0 * D
    # cartesian -> area coordinate coefficients
    a = array([v1[0] * v2[1] - v2[0] * v1[1],
               v2[0] * v0[1] - v0[0] * v2[1],
               v0[0] * v1[1] - v1[0] * v0[1]])
    b = array([v1[1] - v2[1],
               v2[1] - v0[1],
               v0[1] - v1[1]])
    c = array([v2[0] - v1[0],
               v0[0] - v2[0],
               v1[0] - v0[0]])

    def _cart2areaOld(self, x):
        l0 = (self.a[0] + self.b[0] * x[0] + self.c[0] * x[1]) / (self.D2)
        l1 = (self.a[1] + self.b[1] * x[0] + self.c[1] * x[1]) / (self.D2)
        l2 = 1.0 - l0 - l1

        return (l0, l1, l2)

    def _cart2area(self, x):
        return 1 - x[0] - x[1], x[0], x[1]


class simplex_3d(object):
    """ parent class for simplex 2D basis
    """

    # ~ # element vertices
    # ~ v0 = [ 0.0, 0.0, 0.0 ]
    # ~ v1 = [ 1.0, 0.0, 0.0 ]
    # ~ v2 = [ 0.0, 1.0, 0.0 ]
    # ~ v3 = [ 0.0, 0.0, 1.0 ]
    # ~ D = 0.5 * det( [ [1.0, v0[0], v0[1]],\
    # ~ [1.0, v1[0], v1[1]],\
    # ~ [1.0, v2[0], v2[1]] ] )
    # ~ D2 = 2.0 * D
    # ~ # cartesian -> area coordinate coefficients
    # ~ a = array( [ v1[0]*v2[1] - v2[0]*v1[1],\
    # ~ v2[0]*v0[1] - v0[0]*v2[1],\
    # ~ v0[0]*v1[1] - v1[0]*v0[1] ] )
    # ~ b = array( [ v1[1] - v2[1],\
    # ~ v2[1] - v0[1],\
    # ~ v0[1] - v1[1] ] )
    # ~ c = array( [ v2[0] - v1[0],\
    # ~ v0[0] - v2[0],\
    # ~ v1[0] - v0[0] ] )

    def _cart2area(self, x):
        return 1 - x[0] - x[1] - x[2], x[0], x[1], x[2]


class simplex_L1_L1(simplex_2d):
    """ Triangular element bases are hard coded for now
    """

    def __init__(self):
        self.dimensions = 2
        self.type = 'simplex_L1_L1'
        self.basis_order = 1

    def eval(self, x):
        l0, l1, l2 = self._cart2area(x)
        phi = zeros(3)
        phi[0] = 1 - x[0] - x[1]
        phi[1] = x[0]
        phi[2] = x[1]

        return phi


class simplex_L2_L2(simplex_2d):
    """ Triangular element bases are hard coded for now

    See simplex_L2_L2_derivation.wxm for maxima derivation.
    """

    def __init__(self):
        self.dimensions = 2
        self.type = 'simplex_L2_L2'
        self.basis_order = 2
        self.tol = 1e-12

    def eval(self, x):

        phi = array([
            (x[1] + x[0] - 1) * (2 * x[1] + 2 * x[0] - 1),
            -4 * x[0] * (-1 + x[0] + x[1]),
            # 2*x[0]-1,
            x[0] * (2 * x[0] - 1),
            4 * x[0] * x[1],
            x[1] * (2 * x[1] - 1),
            -4 * x[1] * (-1 + x[0] + x[1]),
        ])

        return phi

    def eval_derivatives(self, x, deriv=None):
        """ if deriv==None, evaluate all derivatives
        """

        get_derivatives = {(1, 0): self.eval_dx0,
                           (0, 1): self.eval_dx1,
                           (2, 0): self.eval_dx0x0,
                           (0, 2): self.eval_dx1x1,
                           (1, 1): self.eval_dx0x1}

        if deriv:
            try:
                d = get_derivatives[deriv](x)
            except KeyError:
                raise Warning('derivative ' + str(deriv) + ' not supported')
        else:
            d = array([
                get_derivatives[(1, 0)](x),
                get_derivatives[(0, 1)](x),
                get_derivatives[(2, 0)](x),
                get_derivatives[(0, 2)](x),
                get_derivatives[(1, 1)](x)
            ])

        return d

    def eval_dx0(self, x):

        phi = array([
            2 * (-1 + x[0] + x[1]) + 2 * x[1] + 2 * x[0] - 1,
            -4 * (-1 + x[0] + x[1]) - 4 * x[0],
            4 * x[0] - 1,
            4 * x[1],
            [0, ] * x.shape[1],
            -4 * x[1],
        ])

        return where(abs(phi) < self.tol, 0.0, phi)

    def eval_dx1(self, x):

        phi = array([
            2 * (-1 + x[0] + x[1]) + 2 * x[1] + 2 * x[0] - 1,
            -4 * x[0],
            [0, ] * x.shape[1],
            4 * x[0],
            4 * x[1] - 1,
            -4 * (-1 + x[0] + x[1]) - 4 * x[1],
        ])

        return where(abs(phi) < self.tol, 0.0, phi)

    def eval_dx0x0(self, x):
        phi = ones([6, x.shape[1]])
        phi = phi * array([4.0, -8.0, 4.0, 0.0, 0.0, 0.0])[:, newaxis]
        return phi

    def eval_dx1x1(self, x):
        phi = ones([6, x.shape[1]])
        phi = phi * array([4.0, 0.0, 0.0, 0.0, 4.0, -8.0])[:, newaxis]
        return phi

    def eval_dx0x1(self, x):
        phi = ones([6, x.shape[1]])
        phi = phi * array([4.0, -4.0, 0.0, 4.0, 0.0, -4.0])[:, newaxis]
        return phi


class simplex_L3_L3(simplex_2d):

    def __init__(self):
        self.dimensions = 2
        self.type = 'simplex_L3_L3'
        self.basis_order = 3
        self.tol = 1e-12

        self.D16 = 16.0 * self.D ** 3.0
        self.D916 = -9.0 / self.D16
        self.D27 = 27.0 / (8.0 * self.D ** 3.0)
        self.D4 = self.D * 4.0
        self.D38 = -3.0 / (8.0 * self.D ** 3.0)
        self.D98 = -9.0 / (8.0 * self.D ** 3.0)

    def _T(self, x):
        T0 = self.a[0] + self.b[0] * x[0] + self.c[0] * x[1]
        T1 = self.a[1] + self.b[1] * x[0] + self.c[1] * x[1]
        T2 = self.a[2] + self.b[2] * x[0] + self.c[2] * x[1]

        return T0, T1, T2

    def eval(self, x):
        T0, T1, T2 = self._cart2area(x)
        phi = array([T0 * (self.D2 - 3.0 * T0) * (self.D4 - 3.0 * T0) / self.D16,
                     self.D916 * T0 * T1 * (self.D2 - 3 * T0),
                     self.D916 * T0 * T1 * (self.D2 - 3 * T1),
                     T1 * (self.D2 - 3.0 * T1) * (self.D4 - 3.0 * T1) / self.D16,
                     self.D916 * T0 * T2 * (self.D2 - 3 * T0),
                     self.D27 * T0 * T1 * T2,
                     self.D916 * T1 * T2 * (self.D2 - 3 * T1),
                     self.D916 * T0 * T2 * (self.D2 - 3 * T2),
                     self.D916 * T1 * T2 * (self.D2 - 3 * T2),
                     T2 * (self.D2 - 3.0 * T2) * (self.D4 - 3.0 * T2) / self.D16])

        return where(abs(phi) < self.tol, 0.0, phi)
        # ~ return phi

    def eval_derivatives(self, x, deriv=None):
        """ if deriv==None, evaluate all derivatives
        """

        get_derivatives = {(1, 0): self.eval_dx0,
                           (0, 1): self.eval_dx1,
                           (2, 0): self.eval_dx0x0,
                           (0, 2): self.eval_dx1x1,
                           (1, 1): self.eval_dx0x1}

        if deriv:
            try:
                d = get_derivatives[deriv](x)
            except KeyError:
                raise Warning('derivative ' + str(deriv) + ' not supported')
        else:
            d = array([get_derivatives[(1, 0)](x), get_derivatives[(0, 1)](x), get_derivatives[(2, 0)](x),
                       get_derivatives[(0, 2)](x), get_derivatives[(1, 1)](x)])

        return d

    def eval_dx0(self, x):
        T0, T1, T2 = self._cart2area(x)

        phi = array([(self.b[0] * (self.D2 - 3.0 * T0) * (self.D4 - 3.0 * T0) - 3.0 * self.b[0] * T0 * (
                    self.D4 - 3.0 * T0) - 3.0 * self.b[0] * T0 * (self.D2 - 3.0 * T0)) / self.D16,
                     self.D916 * self.b[0] * T1 * (self.D2 - 3.0 * T0) + self.D916 * self.b[1] * T0 * (
                                 self.D2 - 3.0 * T0) + 0.5 * self.D27 * self.b[0] * T0 * T1,
                     self.D916 * self.b[0] * T1 * (self.D2 - 3.0 * T1) + self.D916 * self.b[1] * T0 * (
                                 self.D2 - 3.0 * T1) + 0.5 * self.D27 * self.b[1] * T0 * T1,
                     (self.b[1] * (self.D2 - 3.0 * T1) * (self.D4 - 3.0 * T1) - 3.0 * self.b[1] * T1 * (
                                 self.D4 - 3.0 * T1) - 3.0 * self.b[1] * T1 * (self.D2 - 3.0 * T1)) / self.D16,
                     self.D916 * self.b[0] * T2 * (self.D2 - 3.0 * T0) + self.D916 * self.b[2] * T0 * (
                                 self.D2 - 3.0 * T0) + 0.5 * self.D27 * self.b[0] * T0 * T2,
                     self.D27 * (self.b[0] * T1 * T2 + self.b[1] * T0 * T2 + self.b[2] * T0 * T1),
                     self.D916 * self.b[1] * T2 * (self.D2 - 3.0 * T1) + self.D916 * self.b[2] * T1 * (
                                 self.D2 - 3.0 * T1) + 0.5 * self.D27 * self.b[1] * T1 * T2,
                     self.D916 * self.b[0] * T2 * (self.D2 - 3.0 * T2) + self.D916 * self.b[2] * T0 * (
                                 self.D2 - 3.0 * T2) + 0.5 * self.D27 * self.b[2] * T0 * T2,
                     self.D916 * self.b[1] * T2 * (self.D2 - 3.0 * T2) + self.D916 * self.b[2] * T1 * (
                                 self.D2 - 3.0 * T2) + 0.5 * self.D27 * self.b[2] * T1 * T2,
                     (self.b[2] * (self.D2 - 3.0 * T2) * (self.D4 - 3.0 * T2) - 3.0 * self.b[2] * T2 * (
                                 self.D4 - 3.0 * T2) - 3.0 * self.b[2] * T2 * (self.D2 - 3.0 * T2)) / self.D16])

        return where(abs(phi) < self.tol, 0.0, phi)
        # ~ return phi

    def eval_dx1(self, x):
        T0, T1, T2 = self._cart2area(x)

        phi = array([(self.c[0] * (self.D2 - 3.0 * T0) * (self.D4 - 3.0 * T0) - 3.0 * self.c[0] * T0 * (
                    self.D4 - 3.0 * T0) - 3.0 * self.c[0] * T0 * (self.D2 - 3.0 * T0)) / self.D16,
                     self.D916 * self.c[0] * T1 * (self.D2 - 3.0 * T0) + self.D916 * self.c[1] * T0 * (
                                 self.D2 - 3.0 * T0) + 0.5 * self.D27 * self.c[0] * T0 * T1,
                     self.D916 * self.c[0] * T1 * (self.D2 - 3.0 * T1) + self.D916 * self.c[1] * T0 * (
                                 self.D2 - 3.0 * T1) + 0.5 * self.D27 * self.c[1] * T0 * T1,
                     (self.c[1] * (self.D2 - 3.0 * T1) * (self.D4 - 3.0 * T1) - 3.0 * self.c[1] * T1 * (
                                 self.D4 - 3.0 * T1) - 3.0 * self.c[1] * T1 * (self.D2 - 3.0 * T1)) / self.D16,
                     self.D916 * self.c[0] * T2 * (self.D2 - 3.0 * T0) + self.D916 * self.c[2] * T0 * (
                                 self.D2 - 3.0 * T0) + 0.5 * self.D27 * self.c[0] * T0 * T2,
                     self.D27 * (self.c[0] * T1 * T2 + self.c[1] * T0 * T2 + self.c[2] * T0 * T1),
                     self.D916 * self.c[1] * T2 * (self.D2 - 3.0 * T1) + self.D916 * self.c[2] * T1 * (
                                 self.D2 - 3.0 * T1) + 0.5 * self.D27 * self.c[1] * T1 * T2,
                     self.D916 * self.c[0] * T2 * (self.D2 - 3.0 * T2) + self.D916 * self.c[2] * T0 * (
                                 self.D2 - 3.0 * T2) + 0.5 * self.D27 * self.c[2] * T0 * T2,
                     self.D916 * self.c[1] * T2 * (self.D2 - 3.0 * T2) + self.D916 * self.c[2] * T1 * (
                                 self.D2 - 3.0 * T2) + 0.5 * self.D27 * self.c[2] * T1 * T2,
                     (self.c[2] * (self.D2 - 3.0 * T2) * (self.D4 - 3.0 * T2) - 3.0 * self.c[2] * T2 * (
                                 self.D4 - 3.0 * T2) - 3.0 * self.c[2] * T2 * (self.D2 - 3.0 * T2)) / self.D16])

        return where(abs(phi) < self.tol, 0.0, phi)
        # ~ return phi

    def eval_dx0x0(self, x):
        T0, T1, T2 = self._cart2area(x)

        phi = array([self.D38 * self.b[0] * self.b[0] * (self.D4 - 3.0 * T0) + self.D38 * self.b[0] * self.b[0] * (
                    self.D2 - 3.0 * T0) - self.D98 * self.b[0] * self.b[0] * T0,
                     self.D98 * self.b[0] * self.b[1] * (self.D2 - 3.0 * T0) + self.D27 * self.b[0] * self.b[
                         0] * T1 + self.D27 * self.b[0] * self.b[1] * T0,
                     self.D98 * self.b[0] * self.b[1] * (self.D2 - 3.0 * T1) + self.D27 * self.b[0] * self.b[
                         1] * T1 + self.D27 * self.b[1] * self.b[1] * T0,
                     self.D38 * self.b[1] * self.b[1] * (self.D4 - 3.0 * T1) + self.D38 * self.b[1] * self.b[1] * (
                                 self.D2 - 3.0 * T1) - self.D98 * self.b[1] * self.b[1] * T1,
                     self.D98 * self.b[0] * self.b[2] * (self.D2 - 3.0 * T0) + self.D27 * self.b[0] * self.b[
                         0] * T2 + self.D27 * self.b[0] * self.b[2] * T0,
                     2.0 * self.D27 * self.b[0] * self.b[1] * T2 + 2.0 * self.D27 * self.b[0] * self.b[
                         2] * T1 + 2.0 * self.D27 * self.b[1] * self.b[2] * T0,
                     self.D98 * self.b[1] * self.b[2] * (self.D2 - 3.0 * T1) + self.D27 * self.b[1] * self.b[
                         1] * T2 + self.D27 * self.b[1] * self.b[2] * T1,
                     self.D98 * self.b[0] * self.b[2] * (self.D2 - 3.0 * T2) + self.D27 * self.b[0] * self.b[
                         2] * T2 + self.D27 * self.b[2] * self.b[2] * T0,
                     self.D98 * self.b[1] * self.b[2] * (self.D2 - 3.0 * T2) + self.D27 * self.b[1] * self.b[
                         2] * T2 + self.D27 * self.b[2] * self.b[2] * T1,
                     self.D38 * self.b[2] * self.b[2] * (self.D4 - 3.0 * T2) + self.D38 * self.b[2] * self.b[2] * (
                                 self.D2 - 3.0 * T2) - self.D98 * self.b[2] * self.b[2] * T2])

        return where(abs(phi) < self.tol, 0.0, phi)
        # ~ return phi

    def eval_dx1x1(self, x):
        T0, T1, T2 = self._cart2area(x)

        phi = array([self.D38 * self.c[0] * self.c[0] * (self.D4 - 3.0 * T0) + self.D38 * self.c[0] * self.c[0] * (
                    self.D2 - 3.0 * T0) - self.D98 * self.c[0] * self.c[0] * T0,
                     self.D98 * self.c[0] * self.c[1] * (self.D2 - 3.0 * T0) + self.D27 * self.c[0] * self.c[
                         0] * T1 + self.D27 * self.c[0] * self.c[1] * T0,
                     self.D98 * self.c[0] * self.c[1] * (self.D2 - 3.0 * T1) + self.D27 * self.c[0] * self.c[
                         1] * T1 + self.D27 * self.c[1] * self.c[1] * T0,
                     self.D38 * self.c[1] * self.c[1] * (self.D4 - 3.0 * T1) + self.D38 * self.c[1] * self.c[1] * (
                                 self.D2 - 3.0 * T1) - self.D98 * self.c[1] * self.c[1] * T1,
                     self.D98 * self.c[0] * self.c[2] * (self.D2 - 3.0 * T0) + self.D27 * self.c[0] * self.c[
                         0] * T2 + self.D27 * self.c[0] * self.c[2] * T0,
                     2.0 * self.D27 * self.c[0] * self.c[1] * T2 + 2.0 * self.D27 * self.c[0] * self.c[
                         2] * T1 + 2.0 * self.D27 * self.c[1] * self.c[2] * T0,
                     self.D98 * self.c[1] * self.c[2] * (self.D2 - 3.0 * T1) + self.D27 * self.c[1] * self.c[
                         1] * T2 + self.D27 * self.c[1] * self.c[2] * T1,
                     self.D98 * self.c[0] * self.c[2] * (self.D2 - 3.0 * T2) + self.D27 * self.c[0] * self.c[
                         2] * T2 + self.D27 * self.c[2] * self.c[2] * T0,
                     self.D98 * self.c[1] * self.c[2] * (self.D2 - 3.0 * T2) + self.D27 * self.c[1] * self.c[
                         2] * T2 + self.D27 * self.c[2] * self.c[2] * T1,
                     self.D38 * self.c[2] * self.c[2] * (self.D4 - 3.0 * T2) + self.D38 * self.c[2] * self.c[2] * (
                                 self.D2 - 3.0 * T2) - self.D98 * self.c[2] * self.c[2] * T2])

        return where(abs(phi) < self.tol, 0.0, phi)
        # ~ return phi

    def eval_dx0x1(self, x):
        T0, T1, T2 = self._cart2area(x)

        phi = array([self.D38 * self.b[0] * self.c[0] * (self.D4 - 3.0 * T0) + self.D38 * self.b[0] * self.c[0] * (
                    self.D2 - 3.0 * T0) - self.D98 * self.b[0] * self.c[0] * T0,
                     0.5 * self.D98 * self.b[0] * self.c[1] * (self.D2 - 3.0 * T0) + 0.5 * self.D98 * self.b[1] *
                     self.c[0] * (self.D2 - 3.0 * T0) + self.D27 * self.b[0] * self.c[0] * T1 + 0.5 * self.D27 * self.b[
                         0] * self.c[1] * T0 + 0.5 * self.D27 * self.b[1] * self.c[0] * T0,
                     0.5 * self.D98 * self.b[0] * self.c[1] * (self.D2 - 3.0 * T1) + 0.5 * self.D98 * self.b[1] *
                     self.c[0] * (self.D2 - 3.0 * T1) + self.D27 * self.b[0] * self.c[1] * T1 + 0.5 * self.D27 * self.b[
                         1] * self.c[0] * T1 + 0.5 * self.D27 * self.b[1] * self.c[1] * T0,
                     self.D38 * self.b[1] * self.c[1] * (self.D4 - 3.0 * T1) + self.D38 * self.b[1] * self.c[1] * (
                                 self.D2 - 3.0 * T1) - self.D98 * self.b[1] * self.c[1] * T1,
                     0.5 * self.D98 * self.b[0] * self.c[2] * (self.D2 - 3.0 * T0) + 0.5 * self.D98 * self.b[2] *
                     self.c[0] * (self.D2 - 3.0 * T0) + self.D27 * self.b[0] * self.c[0] * T2 + 0.5 * self.D27 * self.b[
                         0] * self.c[2] * T0 + 0.5 * self.D27 * self.b[2] * self.c[0] * T0,
                     self.D27 * self.b[0] * self.c[1] * T2 + self.D27 * self.b[1] * self.c[0] * T2 + self.D27 * self.b[
                         0] * self.c[2] * T1 + self.D27 * self.b[2] * self.c[0] * T1 + self.D27 * self.b[1] * self.c[
                         2] * T0 + self.D27 * self.b[2] * self.c[1] * T0,
                     0.5 * self.D98 * self.b[1] * self.c[2] * (self.D2 - 3.0 * T1) + 0.5 * self.D98 * self.b[2] *
                     self.c[1] * (self.D2 - 3.0 * T1) + self.D27 * self.b[1] * self.c[1] * T2 + 0.5 * self.D27 * self.b[
                         1] * self.c[2] * T1 + 0.5 * self.D27 * self.b[2] * self.c[1] * T1,
                     0.5 * self.D98 * self.b[0] * self.c[2] * (self.D2 - 3.0 * T2) + 0.5 * self.D98 * self.b[2] *
                     self.c[0] * (self.D2 - 3.0 * T2) + self.D27 * self.b[0] * self.c[2] * T2 + 0.5 * self.D27 * self.b[
                         2] * self.c[0] * T2 + 0.5 * self.D27 * self.b[2] * self.c[2] * T0,
                     0.5 * self.D98 * self.b[1] * self.c[2] * (self.D2 - 3.0 * T2) + 0.5 * self.D98 * self.b[2] *
                     self.c[1] * (self.D2 - 3.0 * T2) + self.D27 * self.b[1] * self.c[2] * T2 + 0.5 * self.D27 * self.b[
                         2] * self.c[1] * T2 + 0.5 * self.D27 * self.b[2] * self.c[2] * T1,
                     self.D38 * self.b[2] * self.c[2] * (self.D4 - 3.0 * T2) + self.D38 * self.b[2] * self.c[2] * (
                                 self.D2 - 3.0 * T2) - self.D98 * self.b[2] * self.c[2] * T2])

        # ~ self.b[0]*self.c[0]*( self.D38*( self.D4 + self.D2 - 6.0*T0 ) + self.D98*T0 )
        # ~ 0.5*self.D98*self.b[0]*self.c[1]*(self.D2 - 3.0*T0) + 0.5*self.D98*self.b[1]*self.c[0]*(self.D2 - 3.0*T0) + self.D27*self.b[0]*self.c[0]*T1 + 0.5*self.D27*self.b[0]*self.c[1]*T0 + 0.5*self.D27*self.b[1]*self.c[0]*T0,\
        # ~ 0.5*self.D98*self.b[0]*self.c[1]*(self.D2 - 3.0*T1) + 0.5*self.D98*self.b[1]*self.c[0]*(self.D2 - 3.0*T1) + self.D27*self.b[0]*self.c[1]*T1 + 0.5*self.D27*self.b[1]*self.c[0]*T1 + 0.5*self.D27*self.b[1]*self.c[1]*T0,\
        # ~ self.D38*self.b[1]*self.c[1]*(self.D4 - 3.0*T1) + self.D38*self.b[1]*self.c[1]*(self.D2 - 3.0*T1) - self.D98*self.b[1]*self.c[1]*T1,\
        # ~ 0.5*self.D98*self.b[0]*self.c[2]*(self.D2 - 3.0*T0) + 0.5*self.D98*self.b[2]*self.c[0]*(self.D2 - 3.0*T0) + self.D27*self.b[0]*self.c[0]*T2 + 0.5*self.D27*self.b[0]*self.c[2]*T0 + 0.5*self.D27*self.b[2]*self.c[0]*T0,\
        # ~ self.D27*self.b[0]*self.c[1]*T2 + self.D27*self.b[1]*self.c[0]*T2 + self.D27*self.b[0]*self.c[2]*T1 + self.D27*self.b[2]*self.c[0]*T1 + self.D27*self.b[1]*self.c[2]*T0 + self.D27*self.b[2]*self.c[1]*T0,\
        # ~ 0.5*self.D98*self.b[1]*self.c[2]*(self.D2 - 3.0*T1) + 0.5*self.D98*self.b[2]*self.c[1]*(self.D2 - 3.0*T1) + self.D27*self.b[1]*self.c[1]*T2 + 0.5*self.D27*self.b[1]*self.c[2]*T1 + 0.5*self.D27*self.b[2]*self.c[1]*T1,\
        # ~ 0.5*self.D98*self.b[0]*self.c[2]*(self.D2 - 3.0*T2) + 0.5*self.D98*self.b[2]*self.c[0]*(self.D2 - 3.0*T2) + self.D27*self.b[0]*self.c[2]*T2 + 0.5*self.D27*self.b[2]*self.c[0]*T2 + 0.5*self.D27*self.b[2]*self.c[2]*T0,\
        # ~ 0.5*self.D98*self.b[1]*self.c[2]*(self.D2 - 3.0*T2) + 0.5*self.D98*self.b[2]*self.c[1]*(self.D2 - 3.0*T2) + self.D27*self.b[1]*self.c[2]*T2 + 0.5*self.D27*self.b[2]*self.c[1]*T2 + 0.5*self.D27*self.b[2]*self.c[2]*T1,\
        # ~ self.D38*self.b[2]*self.c[2]*(self.D4 - 3.0*T2) + self.D38*self.b[2]*self.c[2]*(self.D2 - 3.0*T2) - self.D98*self.b[2]*self.c[2]*T2 ])

        return where(abs(phi) < self.tol, 0.0, phi)
        # ~ return phi


class simplex_L4_L1(simplex_2d):

    def __init__(self):
        self.dimensions = 2
        self.type = 'simplex_L4_L1'
        self.basis_order = 4
        self.tol = 1e-12

    def eval(self, x):
        L0, L1, L2 = self._cart2area(x)
        phi = array([
            (32 / 3.) * (L0 - 0.75) * (L0 - 0.5) * (L0 - 0.25) * L0,
            (128 / 3.) * (L0 - 0.5) * (L0 - 0.25) * L0 * L1,
            64 * (L0 - 0.25) * L0 * (L1 - 0.25) * L1,
            (128 / 3.) * L0 * (L1 - 0.5) * (L1 - 0.25) * L1,
            (32 / 3.) * (L1 - 0.75) * (L1 - 0.5) * (L1 - 0.25) * L1,
            L2 ** 4.0,
            ])
        return where(abs(phi) < self.tol, 0.0, phi)

    def eval_derivatives(self, x, deriv=None):
        """ if deriv==None, evaluate all derivatives
        """

        get_derivatives = {(1, 0): self.eval_dx0,
                           (0, 1): self.eval_dx1,
                           (2, 0): self.eval_dx0x0,
                           (0, 2): self.eval_dx1x1,
                           (1, 1): self.eval_dx0x1}

        if deriv:
            try:
                d = get_derivatives[deriv](x)
            except KeyError:
                raise Warning('derivative ' + str(deriv) + ' not supported')
        else:
            d = array([get_derivatives[(1, 0)](x),
                       get_derivatives[(0, 1)](x),
                       get_derivatives[(2, 0)](x),
                       get_derivatives[(0, 2)](x),
                       get_derivatives[(1, 1)](x)])

        return d

    def eval_dx0(self, x):
        raise NotImplementedError

    def eval_dx1(self, x):
        raise NotImplementedError

    def eval_dx0x0(self, x):
        raise NotImplementedError

    def eval_dx1x1(self, x):
        raise NotImplementedError

    def eval_dx0x1(self, x):
        raise NotImplementedError


class simplex_L4_L4(simplex_2d):

    def __init__(self):
        self.dimensions = 2
        self.type = 'simplex_L4_L4'
        self.basis_order = 4
        self.tol = 1e-12
        self.D4 = self.D ** 4.0

    def _T(self, x):
        T0 = self.a[0] + self.b[0] * x[0] + self.c[0] * x[1]
        T1 = self.a[1] + self.b[1] * x[0] + self.c[1] * x[1]
        T2 = self.a[2] + self.b[2] * x[0] + self.c[2] * x[1]
        return T0, T1, T2

    def _cart2area(self, x):
        return 1 - x[0] - x[1], x[0], x[1]

    def eval(self, x):
        L1, L2, L3 = self._cart2area(x)
        phi = array([
            (32 / 3.) * (L1 - 0.75) * (L1 - 0.5) * (L1 - 0.25) * L1,
            (128 / 3.) * (L1 - 0.5) * (L1 - 0.25) * L1 * L2,
            64 * (L1 - 0.25) * L1 * (L2 - 0.25) * L2,
            (128 / 3.) * L1 * (L2 - 0.5) * (L2 - 0.25) * L2,
            (32 / 3.) * (L2 - 0.75) * (L2 - 0.5) * (L2 - 0.25) * L2,
            (128 / 3.) * (L1 - 0.5) * (L1 - 0.25) * L1 * L3,
            128 * (L1 - 0.25) * L1 * L2 * L3,
            128 * L1 * (L2 - 0.25) * L2 * L3,
            (128 / 3.) * (L2 - 0.5) * (L2 - 0.25) * L2 * L3,
            64 * (L1 - 0.25) * L1 * (L3 - 0.25) * L3,
            128 * L1 * L2 * (L3 - 0.25) * L3,
            64 * (L2 - 0.25) * L2 * (L3 - 0.25) * L3,
            (128 / 3.) * L1 * (L3 - 0.5) * (L3 - 0.25) * L3,
            (128 / 3.) * L2 * (L3 - 0.5) * (L3 - 0.25) * L3,
            (32 / 3.) * (L3 - 0.75) * (L3 - 0.5) * (L3 - 0.25) * L3
            ])
        return where(abs(phi) < self.tol, 0.0, phi)

    def eval_derivatives(self, x, deriv=None):
        """ if deriv==None, evaluate all derivatives
        """

        get_derivatives = {(1, 0): self.eval_dx0,
                           (0, 1): self.eval_dx1,
                           (2, 0): self.eval_dx0x0,
                           (0, 2): self.eval_dx1x1,
                           (1, 1): self.eval_dx0x1}

        if deriv:
            try:
                d = get_derivatives[deriv](x)
            except KeyError:
                raise Warning('derivative ' + str(deriv) + ' not supported')
        else:
            d = array([get_derivatives[(1, 0)](x),
                       get_derivatives[(0, 1)](x),
                       get_derivatives[(2, 0)](x),
                       get_derivatives[(0, 2)](x),
                       get_derivatives[(1, 1)](x)])

        return d

    def eval_dx0(self, x):
        a0, a1, a2 = self.a
        b0, b1, b2 = self.b
        c0, c1, c2 = self.c
        D = self.D
        D4 = self.D4
        T0, T1, T2 = self._cart2area(x)
        phi = array([-(b0 * (D - 2.0 * T0) * (D - T0) * (3 * D - 2.0 * T0)) / (6 * D4) + (
                    b0 * (T0) * (D - T0) * (3 * D - 2.0 * T0)) / (3 * D4) + (
                                 b0 * (T0) * (D - 2.0 * T0) * (3 * D - 2.0 * T0)) / (6 * D4) + (
                                 b0 * (T0) * (D - 2.0 * T0) * (D - T0)) / (3 * D4),
                     (4 * b0 * (T1) * (D - 2.0 * T0) * (D - T0)) / (3 * D4) + (
                                 4 * b1 * (T0) * (D - 2.0 * T0) * (D - T0)) / (3 * D4) - (
                                 8 * b0 * (T0) * (T1) * (D - T0)) / (3 * D4) - (
                                 4 * b0 * (T0) * (T1) * (D - 2.0 * T0)) / (3 * D4),
                     (b0 * (T1) * (D - 2.0 * T0) * (D - 2.0 * T1)) / D4 + (
                                 b1 * (T0) * (D - 2.0 * T0) * (D - 2.0 * T1)) / D4 - (
                                 2 * b0 * (T0) * (T1) * (D - 2.0 * T1)) / D4 - (
                                 2 * b1 * (T0) * (T1) * (D - 2.0 * T0)) / D4,
                     (4 * b0 * (T1) * (D - 2.0 * T1) * (D - T1)) / (3 * D4) + (
                                 4 * b1 * (T0) * (D - 2.0 * T1) * (D - T1)) / (3 * D4) - (
                                 8 * b1 * (T0) * (T1) * (D - T1)) / (3 * D4) - (
                                 4 * b1 * (T0) * (T1) * (D - 2.0 * T1)) / (3 * D4),
                     -(b1 * (D - 2.0 * T1) * (D - T1) * (3 * D - 2.0 * T1)) / (6 * D4) + (
                                 b1 * (T1) * (D - T1) * (3 * D - 2.0 * T1)) / (3 * D4) + (
                                 b1 * (T1) * (D - 2.0 * T1) * (3 * D - 2.0 * T1)) / (6 * D4) + (
                                 b1 * (T1) * (D - 2.0 * T1) * (D - T1)) / (3 * D4),
                     (4 * b0 * (T2) * (D - 2.0 * T0) * (D - T0)) / (3 * D4) + (
                                 4 * b2 * (T0) * (D - 2.0 * T0) * (D - T0)) / (3 * D4) - (
                                 8 * b0 * (T0) * (T2) * (D - T0)) / (3 * D4) - (
                                 4 * b0 * (T0) * (T2) * (D - 2.0 * T0)) / (3 * D4),
                     -(4 * b0 * (T1) * (T2) * (D - 2.0 * T0)) / D4 - (4 * b1 * (T0) * (T2) * (D - 2.0 * T0)) / D4 - (
                                 4 * b2 * (T0) * (T1) * (D - 2.0 * T0)) / D4 + (8 * b0 * (T0) * (T1) * (T2)) / D4,
                     -(4 * b0 * (T1) * (T2) * (D - 2.0 * T1)) / D4 - (4 * b1 * (T0) * (T2) * (D - 2.0 * T1)) / D4 - (
                                 4 * b2 * (T0) * (T1) * (D - 2.0 * T1)) / D4 + (8 * b1 * (T0) * (T1) * (T2)) / D4,
                     (4 * b1 * (T2) * (D - 2.0 * T1) * (D - T1)) / (3 * D4) + (
                                 4 * b2 * (T1) * (D - 2.0 * T1) * (D - T1)) / (3 * D4) - (
                                 8 * b1 * (T1) * (T2) * (D - T1)) / (3 * D4) - (
                                 4 * b1 * (T1) * (T2) * (D - 2.0 * T1)) / (3 * D4),
                     (b0 * (T2) * (D - 2.0 * T0) * (D - 2.0 * T2)) / D4 + (
                                 b2 * (T0) * (D - 2.0 * T0) * (D - 2.0 * T2)) / D4 - (
                                 2 * b0 * (T0) * (T2) * (D - 2.0 * T2)) / D4 - (
                                 2 * b2 * (T0) * (T2) * (D - 2.0 * T0)) / D4,
                     -(4 * b0 * (T1) * (T2) * (D - 2.0 * T2)) / D4 - (4 * b1 * (T0) * (T2) * (D - 2.0 * T2)) / D4 - (
                                 4 * b2 * (T0) * (T1) * (D - 2.0 * T2)) / D4 + (8 * b2 * (T0) * (T1) * (T2)) / D4,
                     (b1 * (T2) * (D - 2.0 * T1) * (D - 2.0 * T2)) / D4 + (
                                 b2 * (T1) * (D - 2.0 * T1) * (D - 2.0 * T2)) / D4 - (
                                 2 * b1 * (T1) * (T2) * (D - 2.0 * T2)) / D4 - (
                                 2 * b2 * (T1) * (T2) * (D - 2.0 * T1)) / D4,
                     (4 * b0 * (T2) * (D - 2.0 * T2) * (D - T2)) / (3 * D4) + (
                                 4 * b2 * (T0) * (D - 2.0 * T2) * (D - T2)) / (3 * D4) - (
                                 8 * b2 * (T0) * (T2) * (D - T2)) / (3 * D4) - (
                                 4 * b2 * (T0) * (T2) * (D - 2.0 * T2)) / (3 * D4),
                     (4 * b1 * (T2) * (D - 2.0 * T2) * (D - T2)) / (3 * D4) + (
                                 4 * b2 * (T1) * (D - 2.0 * T2) * (D - T2)) / (3 * D4) - (
                                 8 * b2 * (T1) * (T2) * (D - T2)) / (3 * D4) - (
                                 4 * b2 * (T1) * (T2) * (D - 2.0 * T2)) / (3 * D4),
                     -(b2 * (D - 2.0 * T2) * (D - T2) * (3 * D - 2.0 * T2)) / (6 * D4) + (
                                 b2 * (T2) * (D - T2) * (3 * D - 2.0 * T2)) / (3 * D4) + (
                                 b2 * (T2) * (D - 2.0 * T2) * (3 * D - 2.0 * T2)) / (6 * D4) + (
                                 b2 * (T2) * (D - 2.0 * T2) * (D - T2)) / (3 * D4)
                     ])

        return where(abs(phi) < self.tol, 0.0, phi)

    def eval_dx1(self, x):
        a0, a1, a2 = self.a
        b0, b1, b2 = self.b
        c0, c1, c2 = self.c
        D = self.D
        D4 = self.D4
        T0, T1, T2 = self._cart2area(x)
        phi = array([-(c0 * (D - 2.0 * T0) * (D - T0) * (3 * D - 2.0 * T0)) / (6 * D4) + (
                    c0 * (T0) * (D - T0) * (3 * D - 2.0 * T0)) / (3 * D4) + (
                                 c0 * (T0) * (D - 2.0 * T0) * (3 * D - 2.0 * T0)) / (6 * D4) + (
                                 c0 * (T0) * (D - 2.0 * T0) * (D - T0)) / (3 * D4),
                     (4 * c0 * (T1) * (D - 2.0 * T0) * (D - T0)) / (3 * D4) + (
                                 4 * c1 * (T0) * (D - 2.0 * T0) * (D - T0)) / (3 * D4) - (
                                 8 * c0 * (T0) * (T1) * (D - T0)) / (3 * D4) - (
                                 4 * c0 * (T0) * (T1) * (D - 2.0 * T0)) / (3 * D4),
                     (c0 * (T1) * (D - 2.0 * T0) * (D - 2.0 * T1)) / D4 + (
                                 c1 * (T0) * (D - 2.0 * T0) * (D - 2.0 * T1)) / D4 - (
                                 2 * c0 * (T0) * (T1) * (D - 2.0 * T1)) / D4 - (
                                 2 * c1 * (T0) * (T1) * (D - 2.0 * T0)) / D4,
                     (4 * c0 * (T1) * (D - 2.0 * T1) * (D - T1)) / (3 * D4) + (
                                 4 * c1 * (T0) * (D - 2.0 * T1) * (D - T1)) / (3 * D4) - (
                                 8 * c1 * (T0) * (T1) * (D - T1)) / (3 * D4) - (
                                 4 * c1 * (T0) * (T1) * (D - 2.0 * T1)) / (3 * D4),
                     -(c1 * (D - 2.0 * T1) * (D - T1) * (3 * D - 2.0 * T1)) / (6 * D4) + (
                                 c1 * (T1) * (D - T1) * (3 * D - 2.0 * T1)) / (3 * D4) + (
                                 c1 * (T1) * (D - 2.0 * T1) * (3 * D - 2.0 * T1)) / (6 * D4) + (
                                 c1 * (T1) * (D - 2.0 * T1) * (D - T1)) / (3 * D4),
                     (4 * c0 * (T2) * (D - 2.0 * T0) * (D - T0)) / (3 * D4) + (
                                 4 * c2 * (T0) * (D - 2.0 * T0) * (D - T0)) / (3 * D4) - (
                                 8 * c0 * (T0) * (T2) * (D - T0)) / (3 * D4) - (
                                 4 * c0 * (T0) * (T2) * (D - 2.0 * T0)) / (3 * D4),
                     -(4 * c0 * (T1) * (T2) * (D - 2.0 * T0)) / D4 - (4 * c1 * (T0) * (T2) * (D - 2.0 * T0)) / D4 - (
                                 4 * c2 * (T0) * (T1) * (D - 2.0 * T0)) / D4 + (8 * c0 * (T0) * (T1) * (T2)) / D4,
                     -(4 * c0 * (T1) * (T2) * (D - 2.0 * T1)) / D4 - (4 * c1 * (T0) * (T2) * (D - 2.0 * T1)) / D4 - (
                                 4 * c2 * (T0) * (T1) * (D - 2.0 * T1)) / D4 + (8 * c1 * (T0) * (T1) * (T2)) / D4,
                     (4 * c1 * (T2) * (D - 2.0 * T1) * (D - T1)) / (3 * D4) + (
                                 4 * c2 * (T1) * (D - 2.0 * T1) * (D - T1)) / (3 * D4) - (
                                 8 * c1 * (T1) * (T2) * (D - T1)) / (3 * D4) - (
                                 4 * c1 * (T1) * (T2) * (D - 2.0 * T1)) / (3 * D4),
                     (c0 * (T2) * (D - 2.0 * T0) * (D - 2.0 * T2)) / D4 + (
                                 c2 * (T0) * (D - 2.0 * T0) * (D - 2.0 * T2)) / D4 - (
                                 2 * c0 * (T0) * (T2) * (D - 2.0 * T2)) / D4 - (
                                 2 * c2 * (T0) * (T2) * (D - 2.0 * T0)) / D4,
                     -(4 * c0 * (T1) * (T2) * (D - 2.0 * T2)) / D4 - (4 * c1 * (T0) * (T2) * (D - 2.0 * T2)) / D4 - (
                                 4 * c2 * (T0) * (T1) * (D - 2.0 * T2)) / D4 + (8 * c2 * (T0) * (T1) * (T2)) / D4,
                     (c1 * (T2) * (D - 2.0 * T1) * (D - 2.0 * T2)) / D4 + (
                                 c2 * (T1) * (D - 2.0 * T1) * (D - 2.0 * T2)) / D4 - (
                                 2 * c1 * (T1) * (T2) * (D - 2.0 * T2)) / D4 - (
                                 2 * c2 * (T1) * (T2) * (D - 2.0 * T1)) / D4,
                     (4 * c0 * (T2) * (D - 2.0 * T2) * (D - T2)) / (3 * D4) + (
                                 4 * c2 * (T0) * (D - 2.0 * T2) * (D - T2)) / (3 * D4) - (
                                 8 * c2 * (T0) * (T2) * (D - T2)) / (3 * D4) - (
                                 4 * c2 * (T0) * (T2) * (D - 2.0 * T2)) / (3 * D4),
                     (4 * c1 * (T2) * (D - 2.0 * T2) * (D - T2)) / (3 * D4) + (
                                 4 * c2 * (T1) * (D - 2.0 * T2) * (D - T2)) / (3 * D4) - (
                                 8 * c2 * (T1) * (T2) * (D - T2)) / (3 * D4) - (
                                 4 * c2 * (T1) * (T2) * (D - 2.0 * T2)) / (3 * D4),
                     -(c2 * (D - 2.0 * T2) * (D - T2) * (3 * D - 2.0 * T2)) / (6 * D4) + (
                                 c2 * (T2) * (D - T2) * (3 * D - 2.0 * T2)) / (3 * D4) + (
                                 c2 * (T2) * (D - 2.0 * T2) * (3 * D - 2.0 * T2)) / (6 * D4) + (
                                 c2 * (T2) * (D - 2.0 * T2) * (D - T2)) / (3 * D4)
                     ])
        return where(abs(phi) < self.tol, 0.0, phi)

    def eval_dx0x0(self, x):
        a0, a1, a2 = self.a
        b0, b1, b2 = self.b
        c0, c1, c2 = self.c
        D = self.D
        D4 = self.D4
        T0, T1, T2 = self._cart2area(x)
        phi = array([(2 * b0 ** 2 * (D - T0) * (3 * D - 2.0 * T0)) / (3 * D4) + (
                    b0 ** 2 * (D - 2.0 * T0) * (3 * D - 2.0 * T0)) / (3 * D4) - (
                                 2 * b0 ** 2 * (T0) * (3 * D - 2.0 * T0)) / (3 * D4) + (
                                 2 * b0 ** 2 * (D - 2.0 * T0) * (D - T0)) / (3 * D4) - (
                                 4 * b0 ** 2 * (T0) * (D - T0)) / (3 * D4) - (2 * b0 ** 2 * (T0) * (D - 2.0 * T0)) / (
                                 3 * D4),
                     (8 * b0 * b1 * (D - 2.0 * T0) * (D - T0)) / (3 * D4) - (16 * b0 ** 2 * (T1) * (D - T0)) / (
                                 3 * D4) - (16 * b0 * b1 * (T0) * (D - T0)) / (3 * D4) - (
                                 8 * b0 ** 2 * (T1) * (D - 2.0 * T0)) / (3 * D4) - (
                                 8 * b0 * b1 * (T0) * (D - 2.0 * T0)) / (3 * D4) + (16 * b0 ** 2 * (T0) * (T1)) / (
                                 3 * D4),
                     (2 * b0 * b1 * (D - 2.0 * T0) * (D - 2.0 * T1)) / D4 - (
                                 4 * b0 ** 2 * (T1) * (D - 2.0 * T1)) / D4 - (
                                 4 * b0 * b1 * (T0) * (D - 2.0 * T1)) / D4 - (
                                 4 * b0 * b1 * (T1) * (D - 2.0 * T0)) / D4 - (
                                 4 * b1 ** 2 * (T0) * (D - 2.0 * T0)) / D4 + (8 * b0 * b1 * (T0) * (T1)) / D4,
                     (8 * b0 * b1 * (D - 2.0 * T1) * (D - T1)) / (3 * D4) - (16 * b0 * b1 * (T1) * (D - T1)) / (
                                 3 * D4) - (16 * b1 ** 2 * (T0) * (D - T1)) / (3 * D4) - (
                                 8 * b0 * b1 * (T1) * (D - 2.0 * T1)) / (3 * D4) - (
                                 8 * b1 ** 2 * (T0) * (D - 2.0 * T1)) / (3 * D4) + (16 * b1 ** 2 * (T0) * (T1)) / (
                                 3 * D4),
                     (2 * b1 ** 2 * (D - T1) * (3 * D - 2.0 * T1)) / (3 * D4) + (
                                 b1 ** 2 * (D - 2.0 * T1) * (3 * D - 2.0 * T1)) / (3 * D4) - (
                                 2 * b1 ** 2 * (T1) * (3 * D - 2.0 * T1)) / (3 * D4) + (
                                 2 * b1 ** 2 * (D - 2.0 * T1) * (D - T1)) / (3 * D4) - (
                                 4 * b1 ** 2 * (T1) * (D - T1)) / (3 * D4) - (2 * b1 ** 2 * (T1) * (D - 2.0 * T1)) / (
                                 3 * D4),
                     (8 * b0 * b2 * (D - 2.0 * T0) * (D - T0)) / (3 * D4) - (16 * b0 ** 2 * (T2) * (D - T0)) / (
                                 3 * D4) - (16 * b0 * b2 * (T0) * (D - T0)) / (3 * D4) - (
                                 8 * b0 ** 2 * (T2) * (D - 2.0 * T0)) / (3 * D4) - (
                                 8 * b0 * b2 * (T0) * (D - 2.0 * T0)) / (3 * D4) + (16 * b0 ** 2 * (T0) * (T2)) / (
                                 3 * D4),
                     -(8 * b0 * b1 * (T2) * (D - 2.0 * T0)) / D4 - (8 * b0 * b2 * (T1) * (D - 2.0 * T0)) / D4 - (
                                 8 * b1 * b2 * (T0) * (D - 2.0 * T0)) / D4 + (16 * b0 ** 2 * (T1) * (T2)) / D4 + (
                                 16 * b0 * b1 * (T0) * (T2)) / D4 + (16 * b0 * b2 * (T0) * (T1)) / D4,
                     -(8 * b0 * b1 * (T2) * (D - 2.0 * T1)) / D4 - (8 * b0 * b2 * (T1) * (D - 2.0 * T1)) / D4 - (
                                 8 * b1 * b2 * (T0) * (D - 2.0 * T1)) / D4 + (16 * b0 * b1 * (T1) * (T2)) / D4 + (
                                 16 * b1 ** 2 * (T0) * (T2)) / D4 + (16 * b1 * b2 * (T0) * (T1)) / D4,
                     (8 * b1 * b2 * (D - 2.0 * T1) * (D - T1)) / (3 * D4) - (16 * b1 ** 2 * (T2) * (D - T1)) / (
                                 3 * D4) - (16 * b1 * b2 * (T1) * (D - T1)) / (3 * D4) - (
                                 8 * b1 ** 2 * (T2) * (D - 2.0 * T1)) / (3 * D4) - (
                                 8 * b1 * b2 * (T1) * (D - 2.0 * T1)) / (3 * D4) + (16 * b1 ** 2 * (T1) * (T2)) / (
                                 3 * D4),
                     (2 * b0 * b2 * (D - 2.0 * T0) * (D - 2.0 * T2)) / D4 - (
                                 4 * b0 ** 2 * (T2) * (D - 2.0 * T2)) / D4 - (
                                 4 * b0 * b2 * (T0) * (D - 2.0 * T2)) / D4 - (
                                 4 * b0 * b2 * (T2) * (D - 2.0 * T0)) / D4 - (
                                 4 * b2 ** 2 * (T0) * (D - 2.0 * T0)) / D4 + (8 * b0 * b2 * (T0) * (T2)) / D4,
                     -(8 * b0 * b1 * (T2) * (D - 2.0 * T2)) / D4 - (8 * b0 * b2 * (T1) * (D - 2.0 * T2)) / D4 - (
                                 8 * b1 * b2 * (T0) * (D - 2.0 * T2)) / D4 + (16 * b0 * b2 * (T1) * (T2)) / D4 + (
                                 16 * b1 * b2 * (T0) * (T2)) / D4 + (16 * b2 ** 2 * (T0) * (T1)) / D4,
                     (2 * b1 * b2 * (D - 2.0 * T1) * (D - 2.0 * T2)) / D4 - (
                                 4 * b1 ** 2 * (T2) * (D - 2.0 * T2)) / D4 - (
                                 4 * b1 * b2 * (T1) * (D - 2.0 * T2)) / D4 - (
                                 4 * b1 * b2 * (T2) * (D - 2.0 * T1)) / D4 - (
                                 4 * b2 ** 2 * (T1) * (D - 2.0 * T1)) / D4 + (8 * b1 * b2 * (T1) * (T2)) / D4,
                     (8 * b0 * b2 * (D - 2.0 * T2) * (D - T2)) / (3 * D4) - (16 * b0 * b2 * (T2) * (D - T2)) / (
                                 3 * D4) - (16 * b2 ** 2 * (T0) * (D - T2)) / (3 * D4) - (
                                 8 * b0 * b2 * (T2) * (D - 2.0 * T2)) / (3 * D4) - (
                                 8 * b2 ** 2 * (T0) * (D - 2.0 * T2)) / (3 * D4) + (16 * b2 ** 2 * (T0) * (T2)) / (
                                 3 * D4),
                     (8 * b1 * b2 * (D - 2.0 * T2) * (D - T2)) / (3 * D4) - (16 * b1 * b2 * (T2) * (D - T2)) / (
                                 3 * D4) - (16 * b2 ** 2 * (T1) * (D - T2)) / (3 * D4) - (
                                 8 * b1 * b2 * (T2) * (D - 2.0 * T2)) / (3 * D4) - (
                                 8 * b2 ** 2 * (T1) * (D - 2.0 * T2)) / (3 * D4) + (16 * b2 ** 2 * (T1) * (T2)) / (
                                 3 * D4),
                     (2 * b2 ** 2 * (D - T2) * (3 * D - 2.0 * T2)) / (3 * D4) + (
                                 b2 ** 2 * (D - 2.0 * T2) * (3 * D - 2.0 * T2)) / (3 * D4) - (
                                 2 * b2 ** 2 * (T2) * (3 * D - 2.0 * T2)) / (3 * D4) + (
                                 2 * b2 ** 2 * (D - 2.0 * T2) * (D - T2)) / (3 * D4) - (
                                 4 * b2 ** 2 * (T2) * (D - T2)) / (3 * D4) - (2 * b2 ** 2 * (T2) * (D - 2.0 * T2)) / (
                                 3 * D4)
                     ])
        return where(abs(phi) < self.tol, 0.0, phi)

    def eval_dx1x1(self, x):
        a0, a1, a2 = self.a
        b0, b1, b2 = self.b
        c0, c1, c2 = self.c
        D = self.D
        D4 = self.D4
        T0, T1, T2 = self._cart2area(x)
        phi = array([(2 * c0 ** 2 * (D - T0) * (3 * D - 2.0 * T0)) / (3 * D4) + (
                    c0 ** 2 * (D - 2.0 * T0) * (3 * D - 2.0 * T0)) / (3 * D4) - (
                                 2 * c0 ** 2 * (T0) * (3 * D - 2.0 * T0)) / (3 * D4) + (
                                 2 * c0 ** 2 * (D - 2.0 * T0) * (D - T0)) / (3 * D4) - (
                                 4 * c0 ** 2 * (T0) * (D - T0)) / (3 * D4) - (2 * c0 ** 2 * (T0) * (D - 2.0 * T0)) / (
                                 3 * D4),
                     (8 * c0 * c1 * (D - 2.0 * T0) * (D - T0)) / (3 * D4) - (16 * c0 ** 2 * (T1) * (D - T0)) / (
                                 3 * D4) - (16 * c0 * c1 * (T0) * (D - T0)) / (3 * D4) - (
                                 8 * c0 ** 2 * (T1) * (D - 2.0 * T0)) / (3 * D4) - (
                                 8 * c0 * c1 * (T0) * (D - 2.0 * T0)) / (3 * D4) + (16 * c0 ** 2 * (T0) * (T1)) / (
                                 3 * D4),
                     (2 * c0 * c1 * (D - 2.0 * T0) * (D - 2.0 * T1)) / D4 - (
                                 4 * c0 ** 2 * (T1) * (D - 2.0 * T1)) / D4 - (
                                 4 * c0 * c1 * (T0) * (D - 2.0 * T1)) / D4 - (
                                 4 * c0 * c1 * (T1) * (D - 2.0 * T0)) / D4 - (
                                 4 * c1 ** 2 * (T0) * (D - 2.0 * T0)) / D4 + (8 * c0 * c1 * (T0) * (T1)) / D4,
                     (8 * c0 * c1 * (D - 2.0 * T1) * (D - T1)) / (3 * D4) - (16 * c0 * c1 * (T1) * (D - T1)) / (
                                 3 * D4) - (16 * c1 ** 2 * (T0) * (D - T1)) / (3 * D4) - (
                                 8 * c0 * c1 * (T1) * (D - 2.0 * T1)) / (3 * D4) - (
                                 8 * c1 ** 2 * (T0) * (D - 2.0 * T1)) / (3 * D4) + (16 * c1 ** 2 * (T0) * (T1)) / (
                                 3 * D4),
                     (2 * c1 ** 2 * (D - T1) * (3 * D - 2.0 * T1)) / (3 * D4) + (
                                 c1 ** 2 * (D - 2.0 * T1) * (3 * D - 2.0 * T1)) / (3 * D4) - (
                                 2 * c1 ** 2 * (T1) * (3 * D - 2.0 * T1)) / (3 * D4) + (
                                 2 * c1 ** 2 * (D - 2.0 * T1) * (D - T1)) / (3 * D4) - (
                                 4 * c1 ** 2 * (T1) * (D - T1)) / (3 * D4) - (2 * c1 ** 2 * (T1) * (D - 2.0 * T1)) / (
                                 3 * D4),
                     (8 * c0 * c2 * (D - 2.0 * T0) * (D - T0)) / (3 * D4) - (16 * c0 ** 2 * (T2) * (D - T0)) / (
                                 3 * D4) - (16 * c0 * c2 * (T0) * (D - T0)) / (3 * D4) - (
                                 8 * c0 ** 2 * (T2) * (D - 2.0 * T0)) / (3 * D4) - (
                                 8 * c0 * c2 * (T0) * (D - 2.0 * T0)) / (3 * D4) + (16 * c0 ** 2 * (T0) * (T2)) / (
                                 3 * D4),
                     -(8 * c0 * c1 * (T2) * (D - 2.0 * T0)) / D4 - (8 * c0 * c2 * (T1) * (D - 2.0 * T0)) / D4 - (
                                 8 * c1 * c2 * (T0) * (D - 2.0 * T0)) / D4 + (16 * c0 ** 2 * (T1) * (T2)) / D4 + (
                                 16 * c0 * c1 * (T0) * (T2)) / D4 + (16 * c0 * c2 * (T0) * (T1)) / D4,
                     -(8 * c0 * c1 * (T2) * (D - 2.0 * T1)) / D4 - (8 * c0 * c2 * (T1) * (D - 2.0 * T1)) / D4 - (
                                 8 * c1 * c2 * (T0) * (D - 2.0 * T1)) / D4 + (16 * c0 * c1 * (T1) * (T2)) / D4 + (
                                 16 * c1 ** 2 * (T0) * (T2)) / D4 + (16 * c1 * c2 * (T0) * (T1)) / D4,
                     (8 * c1 * c2 * (D - 2.0 * T1) * (D - T1)) / (3 * D4) - (16 * c1 ** 2 * (T2) * (D - T1)) / (
                                 3 * D4) - (16 * c1 * c2 * (T1) * (D - T1)) / (3 * D4) - (
                                 8 * c1 ** 2 * (T2) * (D - 2.0 * T1)) / (3 * D4) - (
                                 8 * c1 * c2 * (T1) * (D - 2.0 * T1)) / (3 * D4) + (16 * c1 ** 2 * (T1) * (T2)) / (
                                 3 * D4),
                     (2 * c0 * c2 * (D - 2.0 * T0) * (D - 2.0 * T2)) / D4 - (
                                 4 * c0 ** 2 * (T2) * (D - 2.0 * T2)) / D4 - (
                                 4 * c0 * c2 * (T0) * (D - 2.0 * T2)) / D4 - (
                                 4 * c0 * c2 * (T2) * (D - 2.0 * T0)) / D4 - (
                                 4 * c2 ** 2 * (T0) * (D - 2.0 * T0)) / D4 + (8 * c0 * c2 * (T0) * (T2)) / D4,
                     -(8 * c0 * c1 * (T2) * (D - 2.0 * T2)) / D4 - (8 * c0 * c2 * (T1) * (D - 2.0 * T2)) / D4 - (
                                 8 * c1 * c2 * (T0) * (D - 2.0 * T2)) / D4 + (16 * c0 * c2 * (T1) * (T2)) / D4 + (
                                 16 * c1 * c2 * (T0) * (T2)) / D4 + (16 * c2 ** 2 * (T0) * (T1)) / D4,
                     (2 * c1 * c2 * (D - 2.0 * T1) * (D - 2.0 * T2)) / D4 - (
                                 4 * c1 ** 2 * (T2) * (D - 2.0 * T2)) / D4 - (
                                 4 * c1 * c2 * (T1) * (D - 2.0 * T2)) / D4 - (
                                 4 * c1 * c2 * (T2) * (D - 2.0 * T1)) / D4 - (
                                 4 * c2 ** 2 * (T1) * (D - 2.0 * T1)) / D4 + (8 * c1 * c2 * (T1) * (T2)) / D4,
                     (8 * c0 * c2 * (D - 2.0 * T2) * (D - T2)) / (3 * D4) - (16 * c0 * c2 * (T2) * (D - T2)) / (
                                 3 * D4) - (16 * c2 ** 2 * (T0) * (D - T2)) / (3 * D4) - (
                                 8 * c0 * c2 * (T2) * (D - 2.0 * T2)) / (3 * D4) - (
                                 8 * c2 ** 2 * (T0) * (D - 2.0 * T2)) / (3 * D4) + (16 * c2 ** 2 * (T0) * (T2)) / (
                                 3 * D4),
                     (8 * c1 * c2 * (D - 2.0 * T2) * (D - T2)) / (3 * D4) - (16 * c1 * c2 * (T2) * (D - T2)) / (
                                 3 * D4) - (16 * c2 ** 2 * (T1) * (D - T2)) / (3 * D4) - (
                                 8 * c1 * c2 * (T2) * (D - 2.0 * T2)) / (3 * D4) - (
                                 8 * c2 ** 2 * (T1) * (D - 2.0 * T2)) / (3 * D4) + (16 * c2 ** 2 * (T1) * (T2)) / (
                                 3 * D4),
                     (2 * c2 ** 2 * (D - T2) * (3 * D - 2.0 * T2)) / (3 * D4) + (
                                 c2 ** 2 * (D - 2.0 * T2) * (3 * D - 2.0 * T2)) / (3 * D4) - (
                                 2 * c2 ** 2 * (T2) * (3 * D - 2.0 * T2)) / (3 * D4) + (
                                 2 * c2 ** 2 * (D - 2.0 * T2) * (D - T2)) / (3 * D4) - (
                                 4 * c2 ** 2 * (T2) * (D - T2)) / (3 * D4) - (2 * c2 ** 2 * (T2) * (D - 2.0 * T2)) / (
                                 3 * D4),
                     ])
        return where(abs(phi) < self.tol, 0.0, phi)

    def eval_dx0x1(self, x):
        a0, a1, a2 = self.a
        b0, b1, b2 = self.b
        c0, c1, c2 = self.c
        D = self.D
        D4 = self.D4
        T0, T1, T2 = self._cart2area(x)
        phi = array([(2 * b0 * c0 * (D - T0) * (3 * D - 2.0 * T0)) / (3 * D4) + (
                    b0 * c0 * (D - 2.0 * T0) * (3 * D - 2.0 * T0)) / (3 * D4) - (
                                 2 * b0 * c0 * (T0) * (3 * D - 2.0 * T0)) / (3 * D4) + (
                                 2 * b0 * c0 * (D - 2.0 * T0) * (D - T0)) / (3 * D4) - (
                                 4 * b0 * c0 * (T0) * (D - T0)) / (3 * D4) - (2 * b0 * c0 * (T0) * (D - 2.0 * T0)) / (
                                 3 * D4),
                     (4 * b0 * c1 * (D - 2.0 * T0) * (D - T0)) / (3 * D4) + (
                                 4 * b1 * c0 * (D - 2.0 * T0) * (D - T0)) / (3 * D4) - (
                                 16 * b0 * c0 * (T1) * (D - T0)) / (3 * D4) - (8 * b0 * c1 * (T0) * (D - T0)) / (
                                 3 * D4) - (8 * b1 * c0 * (T0) * (D - T0)) / (3 * D4) - (
                                 8 * b0 * c0 * (T1) * (D - 2.0 * T0)) / (3 * D4) - (
                                 4 * b0 * c1 * (T0) * (D - 2.0 * T0)) / (3 * D4) - (
                                 4 * b1 * c0 * (T0) * (D - 2.0 * T0)) / (3 * D4) + (16 * b0 * c0 * (T0) * (T1)) / (
                                 3 * D4),
                     (b0 * c1 * (D - 2.0 * T0) * (D - 2.0 * T1)) / D4 + (
                                 b1 * c0 * (D - 2.0 * T0) * (D - 2.0 * T1)) / D4 - (
                                 4 * b0 * c0 * (T1) * (D - 2.0 * T1)) / D4 - (
                                 2 * b0 * c1 * (T0) * (D - 2.0 * T1)) / D4 - (
                                 2 * b1 * c0 * (T0) * (D - 2.0 * T1)) / D4 - (
                                 2 * b0 * c1 * (T1) * (D - 2.0 * T0)) / D4 - (
                                 2 * b1 * c0 * (T1) * (D - 2.0 * T0)) / D4 - (
                                 4 * b1 * c1 * (T0) * (D - 2.0 * T0)) / D4 + (4 * b0 * c1 * (T0) * (T1)) / D4 + (
                                 4 * b1 * c0 * (T0) * (T1)) / D4,
                     (4 * b0 * c1 * (D - 2.0 * T1) * (D - T1)) / (3 * D4) + (
                                 4 * b1 * c0 * (D - 2.0 * T1) * (D - T1)) / (3 * D4) - (
                                 8 * b0 * c1 * (T1) * (D - T1)) / (3 * D4) - (8 * b1 * c0 * (T1) * (D - T1)) / (
                                 3 * D4) - (16 * b1 * c1 * (T0) * (D - T1)) / (3 * D4) - (
                                 4 * b0 * c1 * (T1) * (D - 2.0 * T1)) / (3 * D4) - (
                                 4 * b1 * c0 * (T1) * (D - 2.0 * T1)) / (3 * D4) - (
                                 8 * b1 * c1 * (T0) * (D - 2.0 * T1)) / (3 * D4) + (16 * b1 * c1 * (T0) * (T1)) / (
                                 3 * D4),
                     (2 * b1 * c1 * (D - T1) * (3 * D - 2.0 * T1)) / (3 * D4) + (
                                 b1 * c1 * (D - 2.0 * T1) * (3 * D - 2.0 * T1)) / (3 * D4) - (
                                 2 * b1 * c1 * (T1) * (3 * D - 2.0 * T1)) / (3 * D4) + (
                                 2 * b1 * c1 * (D - 2.0 * T1) * (D - T1)) / (3 * D4) - (
                                 4 * b1 * c1 * (T1) * (D - T1)) / (3 * D4) - (2 * b1 * c1 * (T1) * (D - 2.0 * T1)) / (
                                 3 * D4),
                     (4 * b0 * c2 * (D - 2.0 * T0) * (D - T0)) / (3 * D4) + (
                                 4 * b2 * c0 * (D - 2.0 * T0) * (D - T0)) / (3 * D4) - (
                                 16 * b0 * c0 * (T2) * (D - T0)) / (3 * D4) - (8 * b0 * c2 * (T0) * (D - T0)) / (
                                 3 * D4) - (8 * b2 * c0 * (T0) * (D - T0)) / (3 * D4) - (
                                 8 * b0 * c0 * (T2) * (D - 2.0 * T0)) / (3 * D4) - (
                                 4 * b0 * c2 * (T0) * (D - 2.0 * T0)) / (3 * D4) - (
                                 4 * b2 * c0 * (T0) * (D - 2.0 * T0)) / (3 * D4) + (16 * b0 * c0 * (T0) * (T2)) / (
                                 3 * D4),
                     -(4 * b0 * c1 * (T2) * (D - 2.0 * T0)) / D4 - (4 * b1 * c0 * (T2) * (D - 2.0 * T0)) / D4 - (
                                 4 * b0 * c2 * (T1) * (D - 2.0 * T0)) / D4 - (
                                 4 * b2 * c0 * (T1) * (D - 2.0 * T0)) / D4 - (
                                 4 * b1 * c2 * (T0) * (D - 2.0 * T0)) / D4 - (
                                 4 * b2 * c1 * (T0) * (D - 2.0 * T0)) / D4 + (16 * b0 * c0 * (T1) * (T2)) / D4 + (
                                 8 * b0 * c1 * (T0) * (T2)) / D4 + (8 * b1 * c0 * (T0) * (T2)) / D4 + (
                                 8 * b0 * c2 * (T0) * (T1)) / D4 + (8 * b2 * c0 * (T0) * (T1)) / D4,
                     -(4 * b0 * c1 * (T2) * (D - 2.0 * T1)) / D4 - (4 * b1 * c0 * (T2) * (D - 2.0 * T1)) / D4 - (
                                 4 * b0 * c2 * (T1) * (D - 2.0 * T1)) / D4 - (
                                 4 * b2 * c0 * (T1) * (D - 2.0 * T1)) / D4 - (
                                 4 * b1 * c2 * (T0) * (D - 2.0 * T1)) / D4 - (
                                 4 * b2 * c1 * (T0) * (D - 2.0 * T1)) / D4 + (8 * b0 * c1 * (T1) * (T2)) / D4 + (
                                 8 * b1 * c0 * (T1) * (T2)) / D4 + (16 * b1 * c1 * (T0) * (T2)) / D4 + (
                                 8 * b1 * c2 * (T0) * (T1)) / D4 + (8 * b2 * c1 * (T0) * (T1)) / D4,
                     (4 * b1 * c2 * (D - 2.0 * T1) * (D - T1)) / (3 * D4) + (
                                 4 * b2 * c1 * (D - 2.0 * T1) * (D - T1)) / (3 * D4) - (
                                 16 * b1 * c1 * (T2) * (D - T1)) / (3 * D4) - (8 * b1 * c2 * (T1) * (D - T1)) / (
                                 3 * D4) - (8 * b2 * c1 * (T1) * (D - T1)) / (3 * D4) - (
                                 8 * b1 * c1 * (T2) * (D - 2.0 * T1)) / (3 * D4) - (
                                 4 * b1 * c2 * (T1) * (D - 2.0 * T1)) / (3 * D4) - (
                                 4 * b2 * c1 * (T1) * (D - 2.0 * T1)) / (3 * D4) + (16 * b1 * c1 * (T1) * (T2)) / (
                                 3 * D4),
                     (b0 * c2 * (D - 2.0 * T0) * (D - 2.0 * T2)) / D4 + (
                                 b2 * c0 * (D - 2.0 * T0) * (D - 2.0 * T2)) / D4 - (
                                 4 * b0 * c0 * (T2) * (D - 2.0 * T2)) / D4 - (
                                 2 * b0 * c2 * (T0) * (D - 2.0 * T2)) / D4 - (
                                 2 * b2 * c0 * (T0) * (D - 2.0 * T2)) / D4 - (
                                 2 * b0 * c2 * (T2) * (D - 2.0 * T0)) / D4 - (
                                 2 * b2 * c0 * (T2) * (D - 2.0 * T0)) / D4 - (
                                 4 * b2 * c2 * (T0) * (D - 2.0 * T0)) / D4 + (4 * b0 * c2 * (T0) * (T2)) / D4 + (
                                 4 * b2 * c0 * (T0) * (T2)) / D4,
                     -(4 * b0 * c1 * (T2) * (D - 2.0 * T2)) / D4 - (4 * b1 * c0 * (T2) * (D - 2.0 * T2)) / D4 - (
                                 4 * b0 * c2 * (T1) * (D - 2.0 * T2)) / D4 - (
                                 4 * b2 * c0 * (T1) * (D - 2.0 * T2)) / D4 - (
                                 4 * b1 * c2 * (T0) * (D - 2.0 * T2)) / D4 - (
                                 4 * b2 * c1 * (T0) * (D - 2.0 * T2)) / D4 + (8 * b0 * c2 * (T1) * (T2)) / D4 + (
                                 8 * b2 * c0 * (T1) * (T2)) / D4 + (8 * b1 * c2 * (T0) * (T2)) / D4 + (
                                 8 * b2 * c1 * (T0) * (T2)) / D4 + (16 * b2 * c2 * (T0) * (T1)) / D4,
                     (b1 * c2 * (D - 2.0 * T1) * (D - 2.0 * T2)) / D4 + (
                                 b2 * c1 * (D - 2.0 * T1) * (D - 2.0 * T2)) / D4 - (
                                 4 * b1 * c1 * (T2) * (D - 2.0 * T2)) / D4 - (
                                 2 * b1 * c2 * (T1) * (D - 2.0 * T2)) / D4 - (
                                 2 * b2 * c1 * (T1) * (D - 2.0 * T2)) / D4 - (
                                 2 * b1 * c2 * (T2) * (D - 2.0 * T1)) / D4 - (
                                 2 * b2 * c1 * (T2) * (D - 2.0 * T1)) / D4 - (
                                 4 * b2 * c2 * (T1) * (D - 2.0 * T1)) / D4 + (4 * b1 * c2 * (T1) * (T2)) / D4 + (
                                 4 * b2 * c1 * (T1) * (T2)) / D4,
                     (4 * b0 * c2 * (D - 2.0 * T2) * (D - T2)) / (3 * D4) + (
                                 4 * b2 * c0 * (D - 2.0 * T2) * (D - T2)) / (3 * D4) - (
                                 8 * b0 * c2 * (T2) * (D - T2)) / (3 * D4) - (8 * b2 * c0 * (T2) * (D - T2)) / (
                                 3 * D4) - (16 * b2 * c2 * (T0) * (D - T2)) / (3 * D4) - (
                                 4 * b0 * c2 * (T2) * (D - 2.0 * T2)) / (3 * D4) - (
                                 4 * b2 * c0 * (T2) * (D - 2.0 * T2)) / (3 * D4) - (
                                 8 * b2 * c2 * (T0) * (D - 2.0 * T2)) / (3 * D4) + (16 * b2 * c2 * (T0) * (T2)) / (
                                 3 * D4),
                     (4 * b1 * c2 * (D - 2.0 * T2) * (D - T2)) / (3 * D4) + (
                                 4 * b2 * c1 * (D - 2.0 * T2) * (D - T2)) / (3 * D4) - (
                                 8 * b1 * c2 * (T2) * (D - T2)) / (3 * D4) - (8 * b2 * c1 * (T2) * (D - T2)) / (
                                 3 * D4) - (16 * b2 * c2 * (T1) * (D - T2)) / (3 * D4) - (
                                 4 * b1 * c2 * (T2) * (D - 2.0 * T2)) / (3 * D4) - (
                                 4 * b2 * c1 * (T2) * (D - 2.0 * T2)) / (3 * D4) - (
                                 8 * b2 * c2 * (T1) * (D - 2.0 * T2)) / (3 * D4) + (16 * b2 * c2 * (T1) * (T2)) / (
                                 3 * D4),
                     (2 * b2 * c2 * (D - T2) * (3 * D - 2.0 * T2)) / (3 * D4) + (
                                 b2 * c2 * (D - 2.0 * T2) * (3 * D - 2.0 * T2)) / (3 * D4) - (
                                 2 * b2 * c2 * (T2) * (3 * D - 2.0 * T2)) / (3 * D4) + (
                                 2 * b2 * c2 * (D - 2.0 * T2) * (D - T2)) / (3 * D4) - (
                                 4 * b2 * c2 * (T2) * (D - T2)) / (3 * D4) - (2 * b2 * c2 * (T2) * (D - 2.0 * T2)) / (
                                 3 * D4),
                     ])
        return where(abs(phi) < self.tol, 0.0, phi)


class simplex_L4_L4_L4(simplex_3d):

    def __init__(self):
        self.dimensions = 3
        self.type = 'simplex_L4_L4_L4'
        self.basis_order = 4
        self.tol = 1e-12

    def eval(self, x):
        L0, L1, L2 = self._cart2area(x)

        phi = array([
            (32 * (L0 - 0.75) * (L0 - 0.5) * (L0 - 0.25) * L0) / 3,
            (128 * (L0 - 0.5) * (L0 - 0.25) * L0 * L1) / 3,
            64 * (L0 - 0.25) * L0 * (L1 - 0.25) * L1,
            (128 * L0 * (L1 - 0.5) * (L1 - 0.25) * L1) / 3,
            (32 * (L1 - 0.75) * (L1 - 0.5) * (L1 - 0.25) * L1) / 3,
            (128 * (L0 - 0.5) * (L0 - 0.25) * L0 * L2) / 3,
            128 * (L0 - 0.25) * L0 * L1 * L2,
            128 * L0 * (L1 - 0.25) * L1 * L2,
            (128 * (L1 - 0.5) * (L1 - 0.25) * L1 * L2) / 3,
            64 * (L0 - 0.25) * L0 * (L2 - 0.25) * L2,
            128 * L0 * L1 * (L2 - 0.25) * L2,
            64 * (L1 - 0.25) * L1 * (L2 - 0.25) * L2,
            (128 * L0 * (L2 - 0.5) * (L2 - 0.25) * L2) / 3,
            (128 * L1 * (L2 - 0.5) * (L2 - 0.25) * L2) / 3,
            (32 * (L2 - 0.75) * (L2 - 0.5) * (L2 - 0.25) * L2) / 3,
            (128 * (L0 - 0.5) * (L0 - 0.25) * L0 * L3) / 3,
            128 * (L0 - 0.25) * L0 * L1 * L3,
            128 * L0 * (L1 - 0.25) * L1 * L3,
            (128 * (L1 - 0.25) * L1 * (L3 - 0.5) * L3) / 3,
            128 * (L0 - 0.25) * L0 * L2 * L3,
            256 * L0 * L1 * L2 * L3,
            128 * (L1 - 0.25) * L1 * L2 * L3,
            128 * L0 * (L2 - 0.25) * L2 * L3,
            128 * L1 * (L2 - 0.25) * L2 * L3,
            (128 * (L2 - 0.5) * (L2 - 0.25) * L2 * L3) / 3,
            64 * (L0 - 0.25) * L0 * (L3 - 0.25) * L3,
            128 * L0 * L1 * (L3 - 0.25) * L3,
            64 * (L1 - 0.25) * L1 * (L3 - 0.25) * L3,
            128 * L0 * L2 * (L3 - 0.25) * L3,
            128 * L1 * L2 * (L3 - 0.25) * L3,
            64 * (L2 - 0.25) * L2 * (L3 - 0.25) * L3,
            (128 * L0 * (L3 - 0.5) * (L3 - 0.25) * L3) / 3,
            (128 * L1 * (L3 - 0.5) * (L3 - 0.25) * L3) / 3,
            (128 * L2 * (L3 - 0.5) * (L3 - 0.25) * L3) / 3,
            (32 * (L3 - 0.75) * (L3 - 0.5) * (L3 - 0.25) * L3) / 3
        ])

        return where(abs(phi) < self.tol, 0.0, phi)

    def eval_derivatives(self, x, deriv=None):
        """ if deriv==None, evaluate all derivatives
        """

        ### no cross derivatives yet
        get_derivatives = {(1, 0, 0): self.eval_dx0,
                           (0, 1, 0): self.eval_dx1,
                           (0, 0, 1): self.eval_dx2,
                           (2, 0, 0): self.eval_dx0x0,
                           (0, 2, 0): self.eval_dx1x1,
                           (0, 0, 2): self.eval_dx2x2,
                           (1, 1, 0): self.eval_dx0x1,
                           (1, 0, 1): self.eval_dx0x2,
                           (0, 1, 1): self.eval_dx1x2,
                           (1, 1, 1): self.eval_dx0x1x2
                           }

        if deriv:
            try:
                d = get_derivatives[deriv](x)
            except KeyError:
                raise Warning('derivative ' + str(deriv) + ' not supported')
        else:
            d = array([get_derivatives[(1, 0, 0)](x), get_derivatives[(0, 1, 0)](x), get_derivatives[(0, 0, 1)](x),
                       get_derivatives[(2, 0, 0)](x), get_derivatives[(0, 2, 0)](x), get_derivatives[(0, 0, 2)](x),
                       get_derivatives[(1, 1, 0)](x), get_derivatives[(1, 0, 1)](x), get_derivatives[(0, 1, 1)](x),
                       get_derivatives[(1, 1, 1)](x),
                       ])

    def eval_dx0(self, x):
        raise NotImplementedError

    def eval_dx1(self, x):
        raise NotImplementedError

    def eval_dx2(self, x):
        raise NotImplementedError

    def eval_dx0x0(self, x):
        raise NotImplementedError

    def eval_dx1x1(self, x):
        raise NotImplementedError

    def eval_dx2x2(self, x):
        raise NotImplementedError

    def eval_dx0x1(self, x):
        raise NotImplementedError

    def eval_dx0x2(self, x):
        raise NotImplementedError

    def eval_dx1x2(self, x):
        raise NotImplementedError

    def eval_dx0x1x2(self, x):
        raise NotImplementedError

    # ======================================================================#


# experimental                                                         #
# ======================================================================#
class simplex_H3_H3:
    dimensions = 2
    type = 'simplex_H3_H3'
    basis_order = 3
    # 10-term cubic basis polynomial coefficients
    A = array([[1.0, 0.0, 0.0, -3.0, 0.0, -3.0, 2.0, 0.0, 0.0, 2.0],
               [0.0, 1.0, 0.0, -2.0, -1.0 / 3.0, 0.0, 1.0, 1.0 / 3.0, -2.0 / 3.0, 0.0],
               [0.0, 0.0, 1.0, 0.0, -1.0 / 3.0, -2.0, 0.0, -2.0 / 3.0, 1.0 / 3.0, 1.0],
               [0.0, 0.0, 0.0, 3.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0],
               [0.0, 0.0, 0.0, 0.0, 1.0 / 3.0, 0.0, 0.0, 2.0 / 3.0, -1.0 / 3.0, 0.0],
               [0.0, 0.0, 0.0, 1.0, -1.0 / 3.0, 0.0, -1.0, -2.0 / 3.0, 1.0 / 3.0, 0.0],
               [0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, -2.0],
               [0.0, 0.0, 0.0, 0.0, -1.0 / 3.0, 1.0, 0.0, 1.0 / 3.0, -2.0 / 3.0, -1.0],
               [0.0, 0.0, 0.0, 0.0, 1.0 / 3.0, 0.0, 0.0, -1.0 / 3.0, 2.0 / 3.0, 0.0]])

    # 16-term cubic with 0 derivative across edge constraint
    # ~ A = array([[ 1.        ,  0.        ,  0.        , -3.        ,  0.        ,
    # ~ -3.        ,  2.        ,  0.        ,  0.        ,  2.        ,
    # ~ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
    # ~ 0.        ],
    # ~ [ 0.        ,  1.        ,  0.        , -2.        ,  0.        ,
    # ~ 0.        ,  1.        ,  0.        , -0.5       ,  0.        ,
    # ~ 0.        ,  0.        , -0.5       ,  0.        ,  0.        ,
    # ~ 0.        ],
    # ~ [ 0.        ,  0.        ,  1.        ,  0.        ,  0.        ,
    # ~ -2.        ,  0.        , -0.5       ,  0.        ,  1.        ,
    # ~ -0.5       ,  0.        ,  0.        ,  0.        ,  0.        ,
    # ~ 0.        ],
    # ~ [ 0.        ,  0.        ,  0.        ,  3.        ,  0.        ,
    # ~ 0.        , -2.        ,  0.        ,  0.        ,  0.        ,
    # ~ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
    # ~ 0.        ],
    # ~ [ 0.        ,  0.        ,  0.        ,  0.        ,  0.25      ,
    # ~ 0.        ,  0.        ,  0.375     , -0.125     ,  0.        ,
    # ~ 0.375     ,  0.        , -0.125     ,  0.        ,  0.        ,
    # ~ 0.        ],
    # ~ [ 0.        ,  0.        , -0.09090909,  0.81818182, -0.18181818,
    # ~ 0.18181818, -0.81818182, -0.18181818,  0.09090909, -0.09090909,
    # ~ -0.18181818,  0.        ,  0.09090909,  0.        ,  0.        ,
    # ~ 0.        ],
    # ~ [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
    # ~ 3.        ,  0.        ,  0.        ,  0.        , -2.        ,
    # ~ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
    # ~ 0.        ],
    # ~ [ 0.        , -0.09090909,  0.        ,  0.18181818, -0.18181818,
    # ~ 0.81818182, -0.09090909,  0.09090909, -0.18181818, -0.81818182,
    # ~ 0.09090909,  0.        , -0.18181818,  0.        ,  0.        ,
    # ~ 0.        ],
    # ~ [ 0.        ,  0.        ,  0.        ,  0.        ,  0.25      ,
    # ~ 0.        ,  0.        , -0.125     ,  0.375     ,  0.        ,
    # ~ -0.125     ,  0.        ,  0.375     ,  0.        ,  0.        ,
    # ~ 0.        ]])

    def __init__(self):
        pass

    def _calc_poly(self, X):
        y = X[1]
        x = X[0]
        x2 = x * x
        x3 = x * x * x
        y2 = y * y
        y3 = y * y * y
        return array([1.0, x, y, x2, x * y, y2, x3, x2 * y, x * y2, y3])
        # ~ return array([ 1.0, x[0], x[1], x[0]**2.0, x[0]*x[1], x[1]**2.0, x[0]**3.0, x[1]*x[0]**2.0, x[0]*x[1]**2.0, x[1]**3.0] )

    # ~ def _calc_poly(self, X):
    # ~ y = X[1]
    # ~ x = X[0]
    # ~ x2 = x*x
    # ~ x3 = x*x*x
    # ~ y2 = y*y
    # ~ y3 = y*y*y
    # ~ return array([ 1.0, x, y, x2, x*y, y2, x3, x2*y, x*y2, y3, x3*y, x2*y2, x*y3, x3*y2, x2*y3, x3*y3 ])

    def eval(self, x):
        X = self._calc_poly(x)
        phi = dot(self.A, X)
        # ~ phi = dot( self.cart_map, phi )

        return phi


class simplex_H3_H3_equi:
    """ tensor product cubic hermite basis function over an equilateral
    simplex element
    """

    dimensions = 2
    type = 'simplex_H3_H3_equi'
    basis_order = 3

    # 10 term cubic over equilateral element
    A = array([[0.5, -1.5, -0.27670243, 0., 0.,
                -1.3609831, 2., 1.10680973, 2., 1.17066414],
               [0.125, -0.25, -0.26073883, -0.5, 0.57735027,
                0.10215054, 1., -0.11174521, -0.33333333, 0.0372484],
               [0., 0., 0.38312645, 0., -1.15470054,
                -0.88479263, 0., 0.77689529, 1.33333333, 0.51083526],
               [0.5, 1.5, -0.27670243, 0., 0.,
                -1.3609831, -2., 1.10680973, -2., 1.17066414],
               [0., 0., 0.38312645, 0., 1.15470054,
                -0.88479263, 0., 0.77689529, -1.33333333, 0.51083526],
               [0.125, 0.25, -0.26073883, -0.5, -0.57735027,
                0.10215054, -1., -0.11174521, 0.33333333, 0.0372484],
               [0., 0., 0.55340487, 0., 0.,
                2.72196621, 0., -2.21361947, 0., -2.34132828],
               [0., 0., 0.10642401, 0., 0.,
                0.42089094, 0., -0.42569605, -1.33333333, -0.62790168],
               [0., 0., 0.10642401, 0., 0.,
                0.42089094, 0., -0.42569605, 1.33333333, -0.62790168]])

    # 16 term with zero deriv across 1-2 edge
    # ~ A = array([[ 0.5       , -1.5       , -0.27670243,  0.        ,  0.07429709,
    # ~ -1.3609831 ,  2.        ,  1.10680973,  1.09383378,  1.17066414,
    # ~ -0.29718834,  0.        ,  0.94728784,  0.        ,  0.        ,
    # ~ 0.        ],
    # ~ [ 0.11401099, -0.21208531, -0.1522682 , -0.45604396,  0.13134035,
    # ~ -0.1043956 ,  0.84834123, -0.03806705,  0.07492665,  0.14803853,
    # ~ 0.03283509,  0.        ,  0.06488839,  0.        ,  0.        ,
    # ~ 0.        ],
    # ~ [ 0.        ,  0.        ,  0.38312645,  0.        , -1.0401592 ,
    # ~ -0.88479263,  0.        ,  0.77689529,  0.68632708,  0.51083526,
    # ~ -0.45816536,  0.        ,  0.59437668,  0.        ,  0.        ,
    # ~ 0.        ],
    # ~ [ 0.5       ,  1.5       , -0.27670243,  0.        , -0.07429709,
    # ~ -1.3609831 , -2.        ,  1.10680973, -1.09383378,  1.17066414,
    # ~ 0.29718834,  0.        , -0.94728784,  0.        ,  0.        ,
    # ~ 0.        ],
    # ~ [ 0.        ,  0.        ,  0.38312645,  0.        ,  1.0401592 ,
    # ~ -0.88479263,  0.        ,  0.77689529, -0.68632708,  0.51083526,
    # ~ 0.45816536,  0.        , -0.59437668,  0.        ,  0.        ,
    # ~ 0.        ],
    # ~ [ 0.11401099,  0.21208531, -0.1522682 , -0.45604396, -0.13134035,
    # ~ -0.1043956 , -0.84834123, -0.03806705, -0.07492665,  0.14803853,
    # ~ -0.03283509,  0.        , -0.06488839,  0.        ,  0.        ,
    # ~ 0.        ],
    # ~ [ 0.        ,  0.        ,  0.55340487,  0.        ,  0.        ,
    # ~ 2.72196621,  0.        , -2.21361947,  0.        , -2.34132828,
    # ~ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
    # ~ 0.        ],
    # ~ [ 0.        ,  0.        ,  0.10642401,  0.        , -0.04953139,
    # ~ 0.42089094,  0.        , -0.42569605, -0.72922252, -0.62790168,
    # ~ 0.19812556,  0.        , -0.63152523,  0.        ,  0.        ,
    # ~ 0.        ],
    # ~ [ 0.        ,  0.        ,  0.10642401,  0.        ,  0.04953139,
    # ~ 0.42089094,  0.        , -0.42569605,  0.72922252, -0.62790168,
    # ~ -0.19812556,  0.        ,  0.63152523,  0.        ,  0.        ,
    # ~ 0.        ]])

    # 16-term with zero derive across edge enforced at 2 points per edge
    # ~ A = array([[ 0.5       , -1.5       , -0.27670243,  0.        ,  0.07429709,
    # ~ -1.3609831 ,  2.        ,  1.10680973,  1.09383378,  1.17066414,
    # ~ -0.29718834,  0.        ,  0.94728784,  0.        ,  0.        ,
    # ~ 0.        ],
    # ~ [ 0.125     , -0.25      ,  0.03608439, -0.5       , -0.07216878,
    # ~ -0.58333333,  1.        , -1.29903811,  0.23809524,  0.4330127 ,
    # ~ 2.59807621,  0.        ,  0.20619652,  0.        ,  0.        ,
    # ~ 0.        ],
    # ~ [ 0.        ,  0.        ,  0.89184325,  0.        , -1.64087265,
    # ~ -2.05962377,  0.        , -1.25797194,  1.651754  ,  1.18912434,
    # ~ 1.94468844,  0.79259251,  0.28054897,  0.95425871,  2.08369382,
    # ~ -0.11440087],
    # ~ [ 0.5       ,  1.5       , -0.27670243,  0.        , -0.07429709,
    # ~ -1.3609831 , -2.        ,  1.10680973, -1.09383378,  1.17066414,
    # ~ 0.29718834,  0.        , -0.94728784,  0.        ,  0.        ,
    # ~ 0.        ],
    # ~ [ 0.        ,  0.        ,  0.89184325,  0.        ,  1.64087265,
    # ~ -2.05962377,  0.        , -1.25797194, -1.651754  ,  1.18912434,
    # ~ -1.94468844,  0.79259251, -0.28054897, -0.95425871,  2.08369382,
    # ~ 0.11440087],
    # ~ [ 0.125     ,  0.25      ,  0.03608439, -0.5       ,  0.07216878,
    # ~ -0.58333333, -1.        , -1.29903811, -0.23809524,  0.4330127 ,
    # ~ -2.59807621,  0.        , -0.20619652,  0.        ,  0.        ,
    # ~ 0.        ],
    # ~ [ 0.        ,  0.        ,  0.55340487,  0.        ,  0.        ,
    # ~ 2.72196621,  0.        , -2.21361947,  0.        , -2.34132828,
    # ~ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
    # ~ 0.        ],
    # ~ [ 0.        ,  0.        ,  0.14667294,  0.        ,  0.34455725,
    # ~ 0.32794003,  0.        , -0.58669175, -0.94423049, -0.57423644,
    # ~ -1.37822898, -1.63555513, -0.90870693, -0.34906643, -2.04125071,
    # ~ 0.23607205],
    # ~ [ 0.        ,  0.        ,  0.14667294,  0.        , -0.34455725,
    # ~ 0.32794003,  0.        , -0.58669175,  0.94423049, -0.57423644,
    # ~ 1.37822898, -1.63555513,  0.90870693,  0.34906643, -2.04125071,
    # ~ -0.23607205]])

    def __init__(self):
        pass

    def _calc_poly(self, X):
        y = X[1]
        x = X[0]
        x2 = x * x
        x3 = x2 * x
        x4 = x3 * x
        y2 = y * y
        y3 = y2 * y
        y4 = y3 * y

        # 10-term cubic
        return array([1.0, x, y, x2, x * y, y2, x3, x2 * y, x * y2, y3])
        # 16-term cubic
        # ~ return array([ 1.0, x, y, x2, x*y, y2, x3, x2*y, x*y2, y3, x3*y, x2*y2, x*y3, x3*y2, x2*y3, x3*y3 ] )

    def eval(self, x):
        X = self._calc_poly(x)
        phi = dot(self.A, X)

        return phi


class simplex_H4_H4_equi:
    """ tensor product cubic hermite basis function over an equilateral
    simplex element
    """

    dimensions = 2
    type = 'simplex_H4_H4_equi'
    basis_order = 3

    # 15-term quadric with zero derive across edge enforced at 2 points per edge
    # ~ A = array([[  3.12500000e-01,  -1.50000000e+00,   0.00000000e+00,
    # ~ 1.50000000e+00,   0.00000000e+00,  -2.50000000e+00,
    # ~ 2.00000000e+00,   0.00000000e+00,   6.00000000e+00,
    # ~ 3.84900179e+00,  -3.00000000e+00,   0.00000000e+00,
    # ~ -2.00000000e+00,  -4.61880215e+00,  -1.66666667e+00],
    # ~ [  9.37500000e-02,  -2.50000000e-01,   3.60843918e-02,
    # ~ -2.50000000e-01,  -7.21687836e-02,  -8.75000000e-01,
    # ~ 1.00000000e+00,  -1.29903811e+00,   1.16666667e+00,
    # ~ 1.29903811e+00,  -5.00000000e-01,   2.59807621e+00,
    # ~ 1.83333333e+00,  -8.66025404e-01,  -5.55555556e-01],
    # ~ [  0.00000000e+00,   0.00000000e+00,   6.13434661e-01,
    # ~ 0.00000000e+00,  -1.80421959e+00,  -1.95833333e+00,
    # ~ 0.00000000e+00,  -1.44337567e-01,   3.50000000e+00,
    # ~ 2.06883846e+00,   0.00000000e+00,   2.59807621e+00,
    # ~ 8.33333333e-01,  -1.63582576e+00,  -7.22222222e-01],
    # ~ [  3.12500000e-01,   1.50000000e+00,   0.00000000e+00,
    # ~ 1.50000000e+00,   0.00000000e+00,  -2.50000000e+00,
    # ~ -2.00000000e+00,   0.00000000e+00,  -6.00000000e+00,
    # ~ 3.84900179e+00,  -3.00000000e+00,   1.29925367e-14,
    # ~ -2.00000000e+00,   4.61880215e+00,  -1.66666667e+00],
    # ~ [  0.00000000e+00,   0.00000000e+00,   6.13434661e-01,
    # ~ 0.00000000e+00,   1.80421959e+00,  -1.95833333e+00,
    # ~ 0.00000000e+00,  -1.44337567e-01,  -3.50000000e+00,
    # ~ 2.06883846e+00,   0.00000000e+00,  -2.59807621e+00,
    # ~ 8.33333333e-01,   1.63582576e+00,  -7.22222222e-01],
    # ~ [  9.37500000e-02,   2.50000000e-01,   3.60843918e-02,
    # ~ -2.50000000e-01,   7.21687836e-02,  -8.75000000e-01,
    # ~ -1.00000000e+00,  -1.29903811e+00,  -1.16666667e+00,
    # ~ 1.29903811e+00,  -5.00000000e-01,  -2.59807621e+00,
    # ~ 1.83333333e+00,   8.66025404e-01,  -5.55555556e-01],
    # ~ [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
    # ~ 0.00000000e+00,   0.00000000e+00,   2.00000000e+00,
    # ~ 0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
    # ~ 1.53960072e+00,   0.00000000e+00,   0.00000000e+00,
    # ~ -8.00000000e+00,  -1.13943198e-14,  -2.66666667e+00],
    # ~ [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
    # ~ 0.00000000e+00,   0.00000000e+00,   1.25000000e+00,
    # ~ 0.00000000e+00,   0.00000000e+00,  -6.66666667e-01,
    # ~ -2.11695099e+00,   0.00000000e+00,   0.00000000e+00,
    # ~ -3.66666667e+00,  -7.69800359e-01,   7.77777778e-01],
    # ~ [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
    # ~ 0.00000000e+00,   0.00000000e+00,   1.25000000e+00,
    # ~ 0.00000000e+00,   0.00000000e+00,   6.66666667e-01,
    # ~ -2.11695099e+00,   0.00000000e+00,   0.00000000e+00,
    # ~ -3.66666667e+00,   7.69800359e-01,   7.77777778e-01]])

    # 15-term quadric with zero derive across edge enforced at 2 points per edge, barycentre = 0,0
    # ~ A = array([[  1.85185185e-01,  -1.11111111e+00,  -6.41500299e-01,
    # ~ 1.33333333e+00,   2.30940108e+00,   0.00000000e+00,
    # ~ 2.00000000e+00,  -1.15470054e+00,   2.00000000e+00,
    # ~ 1.92450090e+00,  -3.00000000e+00,   0.00000000e+00,
    # ~ -2.00000000e+00,  -4.61880215e+00,  -1.66666667e+00],
    # ~ [  5.86419753e-02,  -1.94444444e-01,  -1.97795926e-01,
    # ~ -4.72222222e-01,   3.84900179e-01,  -2.77777778e-02,
    # ~ 1.75000000e+00,  -2.40562612e-01,   4.16666667e-01,
    # ~ 6.57537807e-01,  -5.00000000e-01,   2.59807621e+00,
    # ~ 1.83333333e+00,  -8.66025404e-01,  -5.55555556e-01],
    # ~ [  5.86419753e-02,  -2.68518519e-01,  -6.94958657e-02,
    # ~ 2.77777778e-02,  -1.92450090e-01,  -5.27777778e-01,
    # ~ 7.50000000e-01,   3.36787657e-01,   2.08333333e+00,
    # ~ 1.23488808e+00,   0.00000000e+00,   2.59807621e+00,
    # ~ 8.33333333e-01,  -1.63582576e+00,  -7.22222222e-01],
    # ~ [  1.85185185e-01,   1.11111111e+00,  -6.41500299e-01,
    # ~ 1.33333333e+00,  -2.30940108e+00,   0.00000000e+00,
    # ~ -2.00000000e+00,  -1.15470054e+00,  -2.00000000e+00,
    # ~ 1.92450090e+00,  -3.00000000e+00,   0.00000000e+00,
    # ~ -2.00000000e+00,   4.61880215e+00,  -1.66666667e+00],
    # ~ [  5.86419753e-02,   2.68518519e-01,  -6.94958657e-02,
    # ~ 2.77777778e-02,   1.92450090e-01,  -5.27777778e-01,
    # ~ -7.50000000e-01,   3.36787657e-01,  -2.08333333e+00,
    # ~ 1.23488808e+00,   0.00000000e+00,  -2.59807621e+00,
    # ~ 8.33333333e-01,   1.63582576e+00,  -7.22222222e-01],
    # ~ [  5.86419753e-02,   1.94444444e-01,  -1.97795926e-01,
    # ~ -4.72222222e-01,  -3.84900179e-01,  -2.77777778e-02,
    # ~ -1.75000000e+00,  -2.40562612e-01,  -4.16666667e-01,
    # ~ 6.57537807e-01,  -5.00000000e-01,  -2.59807621e+00,
    # ~ 1.83333333e+00,   8.66025404e-01,  -5.55555556e-01],
    # ~ [  1.85185185e-01,   0.00000000e+00,   1.28300060e+00,
    # ~ -6.66666667e-01,   0.00000000e+00,   2.00000000e+00,
    # ~ 0.00000000e+00,  -4.61880215e+00,   0.00000000e+00,
    # ~ -1.53960072e+00,   0.00000000e+00,   1.93304963e-14,
    # ~ -8.00000000e+00,   0.00000000e+00,  -2.66666667e+00],
    # ~ [  5.86419753e-02,  -7.40740741e-02,   2.67291791e-01,
    # ~ -3.05555556e-01,  -5.77350269e-01,  -1.94444444e-01,
    # ~ 0.00000000e+00,  -2.11695099e+00,  -1.33333333e+00,
    # ~ -1.21885057e+00,   0.00000000e+00,   0.00000000e+00,
    # ~ -3.66666667e+00,  -7.69800359e-01,   7.77777778e-01],
    # ~ [  5.86419753e-02,   7.40740741e-02,   2.67291791e-01,
    # ~ -3.05555556e-01,   5.77350269e-01,  -1.94444444e-01,
    # ~ 0.00000000e+00,  -2.11695099e+00,   1.33333333e+00,
    # ~ -1.21885057e+00,   0.00000000e+00,   0.00000000e+00,
    # ~ -3.66666667e+00,   7.69800359e-01,   7.77777778e-01]])

    # 15-term quadric with zero derive across edge enforced at 2 points per edge
    # value bases replaced with 10-term cubic ones
    A = array([[0.5, -1.5, -0.27670243, 0.0, 0.0, -1.3609831,
                2.0, 1.10680973, 2.0, 1.17066414, 0.0, 0.0, 0.0, 0.0, 0.0],
               [9.37500000e-02, -2.50000000e-01, 3.60843918e-02,
                -2.50000000e-01, -7.21687836e-02, -8.75000000e-01,
                1.00000000e+00, -1.29903811e+00, 1.16666667e+00,
                1.29903811e+00, -5.00000000e-01, 2.59807621e+00,
                1.83333333e+00, -8.66025404e-01, -5.55555556e-01],
               [0.00000000e+00, 0.00000000e+00, 6.13434661e-01,
                0.00000000e+00, -1.80421959e+00, -1.95833333e+00,
                0.00000000e+00, -1.44337567e-01, 3.50000000e+00,
                2.06883846e+00, 0.00000000e+00, 2.59807621e+00,
                8.33333333e-01, -1.63582576e+00, -7.22222222e-01],
               [0.5, 1.5, -0.27670243, 0., 0.,
                -1.3609831, -2., 1.10680973, -2., 1.17066414,
                0.0, 0.0, 0.0, 0.0, 0.0],
               [0.00000000e+00, 0.00000000e+00, 6.13434661e-01,
                0.00000000e+00, 1.80421959e+00, -1.95833333e+00,
                0.00000000e+00, -1.44337567e-01, -3.50000000e+00,
                2.06883846e+00, 0.00000000e+00, -2.59807621e+00,
                8.33333333e-01, 1.63582576e+00, -7.22222222e-01],
               [9.37500000e-02, 2.50000000e-01, 3.60843918e-02,
                -2.50000000e-01, 7.21687836e-02, -8.75000000e-01,
                -1.00000000e+00, -1.29903811e+00, -1.16666667e+00,
                1.29903811e+00, -5.00000000e-01, -2.59807621e+00,
                1.83333333e+00, 8.66025404e-01, -5.55555556e-01],
               [0., 0., 0.55340487, 0., 0.,
                2.72196621, 0., -2.21361947, 0., -2.34132828,
                0.0, 0.0, 0.0, 0.0, 0.0],
               [0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                0.00000000e+00, 0.00000000e+00, 1.25000000e+00,
                0.00000000e+00, 0.00000000e+00, -6.66666667e-01,
                -2.11695099e+00, 0.00000000e+00, 0.00000000e+00,
                -3.66666667e+00, -7.69800359e-01, 7.77777778e-01],
               [0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                0.00000000e+00, 0.00000000e+00, 1.25000000e+00,
                0.00000000e+00, 0.00000000e+00, 6.66666667e-01,
                -2.11695099e+00, 0.00000000e+00, 0.00000000e+00,
                -3.66666667e+00, 7.69800359e-01, 7.77777778e-01]])

    def __init__(self):
        pass

    def _calc_poly(self, X):
        y = X[1]
        x = X[0]
        x2 = x * x
        x3 = x2 * x
        x4 = x3 * x
        y2 = y * y
        y3 = y2 * y
        y4 = y3 * y

        # 15 term quadric
        return array([1.0, x, y, x2, x * y, y2, x3, x2 * y, x * y2, y3, x4, x3 * y, x2 * y2, x * y3, y4])

    def eval(self, x):
        X = self._calc_poly(x)
        phi = dot(self.A, X)

        return phi


class simplex_H2_H2_area:
    """ biquadratic hermite basis for 2D simplex elements
    """
    dimensions = 2
    type = 'simplex_H2_H2_area'
    basis_order = 2

    A = array([[0.32142857, 0.17857143, -0.07142857, 0.32142857, 0.17857143,
                -0.07142857, 0.32142857, 0.17857143, -0.07142857],
               [1.35714286, 0.64285714, 0.14285714, -0.64285714, -0.35714286,
                0.14285714, -0.64285714, -0.35714286, 0.14285714],
               [-0.64285714, -0.35714286, 0.14285714, 1.35714286, 0.64285714,
                0.14285714, -0.64285714, -0.35714286, 0.14285714],
               [-0.64285714, -0.35714286, 0.14285714, -0.64285714, -0.35714286,
                0.14285714, 1.35714286, 0.64285714, 0.14285714],
               [-0.67857143, -0.82142857, -0.07142857, 0.32142857, 0.17857143,
                -0.07142857, 0.32142857, 0.17857143, -0.07142857],
               [0.32142857, 0.17857143, -0.07142857, -0.67857143, -0.82142857,
                -0.07142857, 0.32142857, 0.17857143, -0.07142857],
               [0.32142857, 0.17857143, -0.07142857, 0.32142857, 0.17857143,
                -0.07142857, -0.67857143, -0.82142857, -0.07142857],
               [-1.35714286, -0.64285714, -0.14285714, 0.64285714, 0.35714286,
                0.85714286, 0.64285714, 0.35714286, -0.14285714],
               [0.64285714, 0.35714286, -0.14285714, -1.35714286, -0.64285714,
                -0.14285714, 0.64285714, 0.35714286, 0.85714286],
               [0.64285714, 0.35714286, 0.85714286, 0.64285714, 0.35714286,
                -0.14285714, -1.35714286, -0.64285714, -0.14285714]]).T

    # element vertices
    v1 = [-0.5, 0.0]
    v2 = [0.5, 0.0]
    v3 = [0.0, sqrt(3.0) / 2.0]
    area = 0.5 * det([[1.0, v1[0], v1[1]],
                      [1.0, v2[0], v2[1]],
                      [1.0, v3[0], v3[1]]])

    # cartesian -> area coordinate coefficients
    a = array([v2[0] * v3[1] - v3[0] * v2[1],
               v3[0] * v1[1] - v1[0] * v3[1],
               v1[0] * v2[1] - v2[0] * v1[1]])
    b = array([v2[1] - v3[1],
               v3[1] - v1[1],
               v1[1] - v2[1]])
    c = array([v3[0] - v2[0],
               v1[0] - v3[0],
               v2[0] - v1[0]])

    def __init__(self):
        # ~ self.phi = zeros(9)
        pass

    def _calc_poly(self, L):
        x = L[0]
        y = L[1]
        z = L[2]
        x2 = x * x
        y2 = y * y
        z2 = z * z
        return array([1.0, x, y, z, x2, y2, z2, x * y, y * z, x * z])

    def _cart2area(self, x):
        l0 = (self.a[0] + self.b[0] * x[0] + self.c[0] * x[1]) / (2.0 * self.area)
        l1 = (self.a[1] + self.b[1] * x[0] + self.c[1] * x[1]) / (2.0 * self.area)
        l2 = 1.0 - l0 - l1

        return (l0, l1, l2)

    # ~ def set_element( self, e ):
    # ~ self.delta = e.area
    # ~ self.a = e.a
    # ~ self.b = e.b
    # ~ self.c = e.c

    def eval(self, x, debug=0):
        """ evaluates the basis function at cartesian coordinates
        (x[0], x[1])
        """
        # ~ l0, l1 = x
        # ~ l2 = 1.0 - l1 - l0
        # convert cartesian into area
        L = self._cart2area(x)

        if debug:
            log.debug(L)

        X = self._calc_poly(L)
        phi = dot(self.A, X)

        # ~ self.phi[0] = 2.0*l0 - l0*l0 - 2.0*l0*l1 # u0
        # ~ self.phi[1] = -1*(l0 - l0*l0 - l0*l1)            # du0/dL0
        # ~ self.phi[2] = l0*l2                          # du0/dL2
        # ~
        # ~ self.phi[3] = 2.0*l1 - l1*l1 - 2.0*l1*l2 # u1
        # ~ self.phi[4] = -1*(l1 - l1*l1 - l1*l2)            # du1/dL1
        # ~ self.phi[5] = l0*l1                          # du1/dL0
        # ~
        # ~ self.phi[6] = 2.0*l2 - l2*l2 - 2.0*l0*l2 # u2
        # ~ self.phi[7] = -1*(l2 - l2*l2 - l0*l2)            # du2/dL2
        # ~ self.phi[8] = l1*l2                          # du2/dL1

        return phi


class simplex_H3_H3_area:
    """ biquadratic hermite basis for 2D simplex elements
    """
    dimensions = 2
    type = 'simplex_H3_H3_area'
    basis_order = 3

    def __init__(self):
        self.phi = zeros(9)

    def _cart2area(self, x):
        l0 = (self.a[0] + self.b[0] * x[0] + self.c[0] * x[1]) / (2.0 * self.delta)
        l1 = (self.a[1] + self.b[1] * x[0] + self.c[1] * x[1]) / (2.0 * self.delta)
        l2 = 1.0 - l0 - l1

        return (l0, l1, l2)

    def set_element(self, e):
        self.delta = e.area
        self.a = e.a
        self.b = e.b
        self.c = e.c

    def eval(self, x, debug=0):
        """ evaluates the basis function at cartesian coordinates
        (x[0], x[1])
        """
        # ~ l0, l1 = x
        # ~ l2 = 1-l0-l1

        # convert cartesian into area
        l0, l1, l2 = self._cart2area(x)

        L = 0.5 * l0 * l1 * l2
        l02 = l0 ** 2.0
        l12 = l1 ** 2.0
        l22 = l2 ** 2.0
        self.phi[0] = 3.0 * l02 - 2 * l0 ** 3.0  # u0
        self.phi[2] = -self.c[2] * (l02 * l1 + L) + self.c[1] * (l2 * l02 + L)  # du0/dx2
        self.phi[1] = -self.b[2] * (l02 * l1 + L) + self.b[1] * (l2 * l02 + L)  # du0/dx1

        self.phi[3] = 3.0 * l12 - 2 * l1 ** 3.0  # u1
        self.phi[5] = -self.c[0] * (l12 * l2 + L) + self.c[2] * (l0 * l12 + L)  # du1/dx2
        self.phi[4] = -self.b[0] * (l12 * l2 + L) + self.b[2] * (l0 * l12 + L)  # du1/dx1

        self.phi[6] = 3.0 * l22 - 2 * l2 ** 3.0  # u2
        self.phi[8] = -self.c[1] * (l22 * l0 + L) + self.c[0] * (l1 * l22 + L)  # du2/dx2
        self.phi[7] = -self.b[1] * (l22 * l0 + L) + self.b[0] * (l1 * l22 + L)  # du2/dx1

        if debug:
            log.debug(self.phi)

        return self.phi.copy()


class simplex_B3_B3:
    """ cubic bezier basis for 2D simplex elements
    """
    dimensions = 2
    type = 'simplex_B3_B3'
    basis_order = 3

    def __init__(self):
        self.phi = zeros(10)

    def set_element(self, e):
        self.delta = e.area
        self.a = e.a
        self.b = e.b
        self.c = e.c

    def _cart2area(self, x):
        l0 = (self.a[0] + self.b[0] * x[0] + self.c[0] * x[1]) / (2.0 * self.delta)
        l1 = (self.a[1] + self.b[1] * x[0] + self.c[1] * x[1]) / (2.0 * self.delta)
        l2 = 1.0 - l0 - l1

        return (l0, l1, l2)

    def eval(self, x, debug=0):
        l = self._cart2area(x)

        self.phi[0] = l[0] ** 3.0
        self.phi[1] = 3.0 * l[0] ** 2.0 * l[1]
        self.phi[2] = 3.0 * l[0] * l[1] ** 2.0
        self.phi[3] = l[1] ** 3.0
        self.phi[4] = 3.0 * l[1] ** 2.0 * l[2]
        self.phi[5] = 3.0 * l[1] * l[2] ** 2.0
        self.phi[6] = l[2] ** 3.0
        self.phi[7] = 3.0 * l[2] ** 2.0 * l[0]
        self.phi[8] = 3.0 * l[2] * l[0] ** 2
        self.phi[9] = 6.0 * l[0] * l[1] * l[2]

        return self.phi.copy()


class prism_sL4_sL4_qL4:
    """ prism shaped element basis with quartic lagrange
    """

    def __init__(self):
        self.dimensions = 3
        self.type = 'prism_sL4_sL4_qL4'
        self.basis_order = 4
        self.tol = 1e-12
        self.sL4 = simplex_L4_L4()

    def tensor(self, sp, p3):
        p = array([
            sp[0] * p3[0], sp[1] * p3[0], sp[2] * p3[0], sp[3] * p3[0], sp[4] * p3[0],
            sp[5] * p3[0], sp[6] * p3[0], sp[7] * p3[0], sp[8] * p3[0],
            sp[9] * p3[0], sp[10] * p3[0], sp[11] * p3[0],
            sp[12] * p3[0], sp[13] * p3[0],
            sp[14] * p3[0],

            sp[0] * p3[1], sp[1] * p3[1], sp[2] * p3[1], sp[3] * p3[1], sp[4] * p3[1],
            sp[5] * p3[1], sp[6] * p3[1], sp[7] * p3[1], sp[8] * p3[1],
            sp[9] * p3[1], sp[10] * p3[1], sp[11] * p3[1],
            sp[12] * p3[1], sp[13] * p3[1],
            sp[14] * p3[1],

            sp[0] * p3[2], sp[1] * p3[2], sp[2] * p3[2], sp[3] * p3[2], sp[4] * p3[2],
            sp[5] * p3[2], sp[6] * p3[2], sp[7] * p3[2], sp[8] * p3[2],
            sp[9] * p3[2], sp[10] * p3[2], sp[11] * p3[2],
            sp[12] * p3[2], sp[13] * p3[2],
            sp[14] * p3[2],

            sp[0] * p3[3], sp[1] * p3[3], sp[2] * p3[3], sp[3] * p3[3], sp[4] * p3[3],
            sp[5] * p3[3], sp[6] * p3[3], sp[7] * p3[3], sp[8] * p3[3],
            sp[9] * p3[3], sp[10] * p3[3], sp[11] * p3[3],
            sp[12] * p3[3], sp[13] * p3[3],
            sp[14] * p3[3],

            sp[0] * p3[4], sp[1] * p3[4], sp[2] * p3[4], sp[3] * p3[4], sp[4] * p3[4],
            sp[5] * p3[4], sp[6] * p3[4], sp[7] * p3[4], sp[8] * p3[4],
            sp[9] * p3[4], sp[10] * p3[4], sp[11] * p3[4],
            sp[12] * p3[4], sp[13] * p3[4],
            sp[14] * p3[4],
            ])

        return where(abs(p) < self.tol, 0.0, p)

    def eval(self, x):
        return self.tensor(self.sL4.eval(x[:2]), L4(x[2]))

    def eval_derivatives(self, x, deriv=None):
        """ if deriv==None, evaluate all derivatives
        """
        ### no cross derivatives yet
        get_derivatives = {(1, 0, 0): self.eval_dx0,
                           (0, 1, 0): self.eval_dx1,
                           (0, 0, 1): self.eval_dx2,
                           (2, 0, 0): self.eval_dx0x0,
                           (0, 2, 0): self.eval_dx1x1,
                           (0, 0, 2): self.eval_dx2x2,
                           (1, 1, 0): self.eval_dx0x1,
                           (1, 0, 1): self.eval_dx0x2,
                           (0, 1, 1): self.eval_dx1x2,
                           (1, 1, 1): self.eval_dx0x1x2,
                           }

        if deriv:
            try:
                d = get_derivatives[deriv](x)
            except KeyError:
                raise Warning('derivative ' + str(deriv) + ' not supported')
        else:
            d = array([get_derivatives[(1, 0, 0)](x),
                       get_derivatives[(0, 1, 0)](x),
                       get_derivatives[(0, 0, 1)](x),
                       get_derivatives[(2, 0, 0)](x),
                       get_derivatives[(0, 2, 0)](x),
                       get_derivatives[(0, 0, 2)](x),
                       get_derivatives[(1, 1, 0)](x),
                       get_derivatives[(1, 0, 1)](x),
                       get_derivatives[(0, 1, 1)](x),
                       get_derivatives[(1, 1, 1)](x),
                       ])

        return d

    def eval_dx0(self, x):
        self.tensor(self.sL4.eval_dx0(x[:2]), L4(x[2]))

    def eval_dx1(self, x):
        self.tensor(self.sL4.eval_dx1(x[:2]), L4(x[2]))

    def eval_dx2(self, x):
        self.tensor(self.sL4.eval(x[:2]), L4d(x[2]))

    def eval_dx0x0(self, x):
        self.tensor(self.sL4.eval_dx0x0(x[:2]), L4(x[2]))

    def eval_dx1x1(self, x):
        self.tensor(self.sL4.eval_dx1x1(x[:2]), L4(x[2]))

    def eval_dx2x2(self, x):
        self.tensor(self.sL4.eval(x[:2]), L4dd(x[2]))

    def eval_dx0x1(self, x):
        self.tensor(self.sL4.eval_dx0x1(x[:2]), L4(x[2]))

    def eval_dx0x2(self, x):
        self.tensor(self.sL4.eval_dx0(x[:2]), L4d(x[2]))

    def eval_dx1x2(self, x):
        self.tensor(self.sL4.eval_dx1(x[:2]), L4d(x[2]))

    def eval_dx0x1x2(self, x):
        self.tensor(self.sL4.eval_dx0x1(x[:2]), L4d(x[2]))


class prism_sL4_sL4_qL1:
    """ prism shaped element basis with quartic lagrange in triangle face
    and linear in z
    """

    def __init__(self):
        self.dimensions = 3
        self.type = 'prism_sL4_sL4_qL1'
        self.basis_order = 4
        self.tol = 1e-12
        self.sL4 = simplex_L4_L4()

    def tensor(self, sp, p3):
        p = array([
            sp[0] * p3[0], sp[1] * p3[0], sp[2] * p3[0], sp[3] * p3[0], sp[4] * p3[0],
            sp[5] * p3[0], sp[6] * p3[0], sp[7] * p3[0], sp[8] * p3[0],
            sp[9] * p3[0], sp[10] * p3[0], sp[11] * p3[0],
            sp[12] * p3[0], sp[13] * p3[0],
            sp[14] * p3[0],

            sp[0] * p3[1], sp[1] * p3[1], sp[2] * p3[1], sp[3] * p3[1], sp[4] * p3[1],
            sp[5] * p3[1], sp[6] * p3[1], sp[7] * p3[1], sp[8] * p3[1],
            sp[9] * p3[1], sp[10] * p3[1], sp[11] * p3[1],
            sp[12] * p3[1], sp[13] * p3[1],
            sp[14] * p3[1],
            ])

        return where(abs(p) < self.tol, 0.0, p)

    def eval(self, x):
        return self.tensor(self.sL4.eval(x[:2]), L1(x[2]))

    def eval_derivatives(self, x, deriv=None):
        """ if deriv==None, evaluate all derivatives
        """
        ### no cross derivatives yet
        get_derivatives = {(1, 0, 0): self.eval_dx0,
                           (0, 1, 0): self.eval_dx1,
                           (0, 0, 1): self.eval_dx2,
                           (2, 0, 0): self.eval_dx0x0,
                           (0, 2, 0): self.eval_dx1x1,
                           (0, 0, 2): self.eval_dx2x2,
                           (1, 1, 0): self.eval_dx0x1,
                           (1, 0, 1): self.eval_dx0x2,
                           (0, 1, 1): self.eval_dx1x2,
                           (1, 1, 1): self.eval_dx0x1x2,
                           }

        if deriv:
            try:
                d = get_derivatives[deriv](x)
            except KeyError:
                raise Warning('derivative ' + str(deriv) + ' not supported')
        else:
            d = array([get_derivatives[(1, 0, 0)](x),
                       get_derivatives[(0, 1, 0)](x),
                       get_derivatives[(0, 0, 1)](x),
                       get_derivatives[(2, 0, 0)](x),
                       get_derivatives[(0, 2, 0)](x),
                       get_derivatives[(0, 0, 2)](x),
                       get_derivatives[(1, 1, 0)](x),
                       get_derivatives[(1, 0, 1)](x),
                       get_derivatives[(0, 1, 1)](x),
                       get_derivatives[(1, 1, 1)](x),
                       ])

        return d

    def eval_dx0(self, x):
        self.tensor(self.sL4.eval_dx0(x[:2]), L1(x[2]))

    def eval_dx1(self, x):
        self.tensor(self.sL4.eval_dx1(x[:2]), L1(x[2]))

    def eval_dx2(self, x):
        self.tensor(self.sL4.eval(x[:2]), L1d(x[2]))

    def eval_dx0x0(self, x):
        self.tensor(self.sL4.eval_dx0x0(x[:2]), L1(x[2]))

    def eval_dx1x1(self, x):
        self.tensor(self.sL4.eval_dx1x1(x[:2]), L1(x[2]))

    def eval_dx2x2(self, x):
        self.tensor(self.sL4.eval(x[:2]), L1dd(x[2]))

    def eval_dx0x1(self, x):
        self.tensor(self.sL4.eval_dx0x1(x[:2]), L1(x[2]))

    def eval_dx0x2(self, x):
        self.tensor(self.sL4.eval_dx0(x[:2]), L1d(x[2]))

    def eval_dx1x2(self, x):
        self.tensor(self.sL4.eval_dx1(x[:2]), L1d(x[2]))

    def eval_dx0x1x2(self, x):
        self.tensor(self.sL4.eval_dx0x1(x[:2]), L1d(x[2]))


class prism_sL4_sL1_qL4:
    """ prism shaped element basis with quartic lagrange in triangle face
    and linear in z
    """

    def __init__(self):
        self.dimensions = 3
        self.type = 'prism_sL4_sL4_qL1'
        self.basis_order = 4
        self.tol = 1e-12
        self.sL4 = simplex_L4_L1()

    def tensor(self, sp, p3):
        p = array([
            sp[0] * p3[0], sp[1] * p3[0], sp[2] * p3[0], sp[3] * p3[0], sp[4] * p3[0],
            sp[5] * p3[0],

            sp[0] * p3[1], sp[1] * p3[1], sp[2] * p3[1], sp[3] * p3[1], sp[4] * p3[1],
            sp[5] * p3[1],

            sp[0] * p3[2], sp[1] * p3[2], sp[2] * p3[2], sp[3] * p3[2], sp[4] * p3[2],
            sp[5] * p3[2],

            sp[0] * p3[3], sp[1] * p3[3], sp[2] * p3[3], sp[3] * p3[3], sp[4] * p3[3],
            sp[5] * p3[3],

            sp[0] * p3[4], sp[1] * p3[4], sp[2] * p3[4], sp[3] * p3[4], sp[4] * p3[4],
            sp[5] * p3[4],
            ])

        return where(abs(p) < self.tol, 0.0, p)

    def eval(self, x):
        return self.tensor(self.sL4.eval(x[:2]), L4(x[2]))

    def eval_derivatives(self, x, deriv=None):
        """ if deriv==None, evaluate all derivatives
        """
        ### no cross derivatives yet
        get_derivatives = {(1, 0, 0): self.eval_dx0,
                           (0, 1, 0): self.eval_dx1,
                           (0, 0, 1): self.eval_dx2,
                           (2, 0, 0): self.eval_dx0x0,
                           (0, 2, 0): self.eval_dx1x1,
                           (0, 0, 2): self.eval_dx2x2,
                           (1, 1, 0): self.eval_dx0x1,
                           (1, 0, 1): self.eval_dx0x2,
                           (0, 1, 1): self.eval_dx1x2,
                           (1, 1, 1): self.eval_dx0x1x2,
                           }

        if deriv:
            try:
                d = get_derivatives[deriv](x)
            except KeyError:
                raise Warning('derivative ' + str(deriv) + ' not supported')
        else:
            d = array([get_derivatives[(1, 0, 0)](x),
                       get_derivatives[(0, 1, 0)](x),
                       get_derivatives[(0, 0, 1)](x),
                       get_derivatives[(2, 0, 0)](x),
                       get_derivatives[(0, 2, 0)](x),
                       get_derivatives[(0, 0, 2)](x),
                       get_derivatives[(1, 1, 0)](x),
                       get_derivatives[(1, 0, 1)](x),
                       get_derivatives[(0, 1, 1)](x),
                       get_derivatives[(1, 1, 1)](x),
                       ])

        return d

    def eval_dx0(self, x):
        self.tensor(self.sL4.eval_dx0(x[:2]), L1(x[2]))

    def eval_dx1(self, x):
        self.tensor(self.sL4.eval_dx1(x[:2]), L1(x[2]))

    def eval_dx2(self, x):
        self.tensor(self.sL4.eval(x[:2]), L1d(x[2]))

    def eval_dx0x0(self, x):
        self.tensor(self.sL4.eval_dx0x0(x[:2]), L1(x[2]))

    def eval_dx1x1(self, x):
        self.tensor(self.sL4.eval_dx1x1(x[:2]), L1(x[2]))

    def eval_dx2x2(self, x):
        self.tensor(self.sL4.eval(x[:2]), L1dd(x[2]))

    def eval_dx0x1(self, x):
        self.tensor(self.sL4.eval_dx0x1(x[:2]), L1(x[2]))

    def eval_dx0x2(self, x):
        self.tensor(self.sL4.eval_dx0(x[:2]), L1d(x[2]))

    def eval_dx1x2(self, x):
        self.tensor(self.sL4.eval_dx1(x[:2]), L1d(x[2]))

    def eval_dx0x1x2(self, x):
        self.tensor(self.sL4.eval_dx0x1(x[:2]), L1d(x[2]))


# ======================================================================#
# ======================================================================#
basis_types = {'line_L2': line_L2,
               'line_L3': line_L3,
               'line_L4': line_L4,
               'simplex_linear_lagrange_2d': simplex_L1_L1,
               'simplex_L1_L1': simplex_L1_L1,
               'simplex_quadratic_lagrange_2d': simplex_L2_L2,
               'simplex_L2_L2': simplex_L2_L2,
               'simplex_cubic_lagrange_2d': simplex_L3_L3,
               'simplex_L3_L3': simplex_L3_L3,
               'simplex_L4_L4': simplex_L4_L4,
               'simplex_L4_L1': simplex_L4_L1,
               'simplex_L4_L4_L4': simplex_L4_L4_L4,
               'simplex_quadratic_hermite_2d_area': simplex_H2_H2_area,
               'simplex_H2_H2_area': simplex_H2_H2_area,
               'simplex_cubic_hermite_2d_area': simplex_H3_H3_area,
               'simplex_H3_H3_area': simplex_H3_H3_area,
               'simplex_cubic_hermite_2d': simplex_H3_H3,
               'simplex_H3_H3': simplex_H3_H3,
               'simplex_cubic_hermite_2d_equi': simplex_H3_H3_equi,
               'simplex_H3_H3_equi': simplex_H3_H3_equi,
               'simplex_bezier_2d': simplex_B3_B3,
               'simplex_B3_B3': simplex_B3_B3,
               'simplex_quadric_hermite_2d_equi': simplex_H4_H4_equi,
               'simplex_H4_H4_equi': simplex_H4_H4_equi,
               'quad_cubic_lagrange': quad_L3_L3,
               'quad_L2_L2': quad_L2_L2,
               'quad_L3_L3': quad_L3_L3,
               'quad_L4_L3': quad_L4_L3,
               'quad_L4_L4': quad_L4_L4,
               'quad_L2_L2_L2': quad_L2_L2_L2,
               'quad_cubic_lagrange_3d': quad_L3_L3_L3,
               'quad_L3_L3_L3': quad_L3_L3_L3,
               'quad_quartic_lagrange_3d': quad_L4_L4_L4,
               'quad_L4_L4_L1': quad_L4_L4_L1,
               'quad_L4_L3_L1': quad_L4_L3_L1,
               'quad_L4_L4_L4': quad_L4_L4_L4,
               'prism_sL4_sL1_qL4': prism_sL4_sL1_qL4,
               'prism_sL4_sL4_qL1': prism_sL4_sL4_qL1,
               'prism_sL4_sL4_qL4': prism_sL4_sL4_qL4,
               }


def make_basis(basis_type):
    return basis_types[basis_type]()
