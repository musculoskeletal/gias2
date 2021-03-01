import unittest

from generator import generator, generate
from numpy import array as array
from numpy.linalg import norm

from gias2.common.geoprimitives import Line3D, NonInterceptError, CollinearError, PRECISION

tol = PRECISION


@generator
class TestLine3D(unittest.TestCase):
    """
        Test class to test modules within the Line3D class of Gias2.common.geoprimitives.
        Where array comparison are required, they have been done by evaluating the norm of
        the differences between the arrays.
        Method taken from ENGSCI233 Quality Control lectures, lab material.
    """

    # The below cases should have intersects (hand calculated)
    @generate(
        ('Standard case1', [1, 0, 0], [0, 0, 0], [1, 0, -1], [0, 0, 1], [1, 0, 0]),
        ('Case1 in Z axis', [0, 0, 1], [0, 0, 0], [1, 0, -1], [0, 0, 1], [0, 0, 1]),
        ('Standard case2', [1, 0, 1], [0, 0, 0], [1, 0, 0], [0, 0, 1], [1, 0, 1]),
        ('Case2, off-origin', [1, 0, 1], [-1, -1, -1], [1, 0, 0], [-1, -1, 0], [0, -1, 0]),
        ('Standard case3', [1, 1, 0], [0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]),
        ('Standard case4', [-1, -1, 0], [0, 0, 0], [-1, 0, 0], [0, -1, 0], [-1, -1, 0]),
        ('Case 4, off-origin', [-1, -1, 0], [1, 1, 1], [-1, 0, 0], [1, 0, 1], [0, 0, 1]),
        ('Standard case5', [1, 1, 1], [0, 0, 0], [1, 1, 0], [0, 0, 7], [7, 7, 7]),
        ('Case5, off-origin', [1, 1, 1], [1, 1, 1], [1, 1, 0], [1, 1, 8], [8, 8, 8])
    )
    def test_calcIntercept_intercept_exists(
            self,
            name: str,
            ray1: list,
            origin1: list,
            ray2: list,
            origin2: list,
            expected: list):
        line_a = Line3D(array(ray1), array(origin1))
        line_b = Line3D(array(ray2), array(origin2))
        self.assertTrue(norm(line_a.calcIntercept(line_b)[0] - array(expected)) < tol)

    # These lines are collinear and should result in the CollinearError exception being raised
    @generate(
        ('Same line (collinear) 1D (x)', [1, 0, 0], [0, 0, 0], [1, 0, 0], [0, 0, 0]),
        ('Same line (collinear) 1D (y)', [0, 1, 0], [0, 0, 0], [0, 1, 0], [0, 0, 0]),
        ('Same line (collinear) 1D (z)', [0, 0, 1], [0, 0, 0], [0, 0, 1], [0, 0, 0]),
        ('Same direction, different origin (collinear) 1D (x)', [1, 0, 0], [0, 0, 0], [1, 0, 0], [1, 0, 0]),
        ('Same direction, different origin (collinear) 1D (y)', [0, 1, 0], [0, 0, 0], [0, 1, 0], [0, 1, 0]),
        ('Same direction, different origin (collinear) 1D (z)', [0, 0, 1], [0, 0, 0], [0, 0, 1], [0, 0, 1]),
        ('Same direction, different origin (collinear) 2D (xy)', [1, 1, 0], [0, 0, 0], [-1, -1, 0], [1, 1, 0]),
        ('Same direction, different origin (collinear) 2D (xz)', [1, 0, 1], [0, 0, 0], [-1, 0, -1], [2, 0, 2]),
        ('Same direction, different origin (collinear) 2D (yz)', [0, 1, 1], [0, 0, 0], [0, -1, -1], [0, 1, 1])
    )
    def test_calcIntercept_errors_trivial_cases_collinear(
            self,
            name: str,
            ray1: list,
            origin1: list,
            ray2: list,
            origin2: list):
        line_a = Line3D(array(ray1), array(origin1))
        line_b = Line3D(array(ray2), array(origin2))
        self.assertRaises(CollinearError, lambda: line_a.calcIntercept(line_b))

    # These tests do not have intersects, so should pass when given a NonInterceptError.
    @generate(
        ('Parallel lines, different origin', [1, 0, 0], [0, 0, 0], [1, 0, 0], [1, 1, 1]),
        ('Anti-Parallel lines, different origin', [1, 0, 0], [0, 0, 0], [-1, 0, 0], [1, 1, 1]),
        ('Orthogonal rays, different origins', [1, 0, 0], [0, 0, 0], [0, 1, 0], [1, 1, 1])
    )
    def test_calcIntercept_no_intercept(
            self,
            name: str,
            ray1: list,
            origin1: list,
            ray2: list,
            origin2: list):
        line_a = Line3D(array(ray1), array(origin1))
        line_b = Line3D(array(ray2), array(origin2))
        self.assertRaises(NonInterceptError, lambda: line_a.calcIntercept(line_b))
