import unittest

import numpy as np

from gias2.common import geoprimitives


class GeoPrimitivesTestSuite(unittest.TestCase):

    @staticmethod
    def _make_sphere_point_cloud(n, x, y, z, r):
        # generate n random points at distance r from origin
        _x = np.random.uniform(-r, r, n)
        _y_limits = np.sqrt(r * r - _x * _x)
        _y = np.array([np.random.uniform(-yl, yl) for yl in _y_limits])
        _z = np.sign(np.random.uniform(-1, 1, n)) * np.sqrt(r * r - _x * _x - _y * _y)

        # translate points to x,y,z
        return np.array([_x + x, _y + y, _z + z]).T

    @staticmethod
    def _make_hemisphere_point_cloud(n, x, y, z, r):
        # generate n random points at distance r from origin
        _x = np.random.uniform(0, r, n)
        _y_limits = np.sqrt(r * r - _x * _x)
        _y = np.array([np.random.uniform(-yl, yl) for yl in _y_limits])
        _z = np.sign(np.random.uniform(-1, 1, n)) * np.sqrt(r * r - _x * _x - _y * _y)

        # translate points to x,y,z
        return np.array([_x + x, _y + y, _z + z]).T

    def test_fit_sphere_analytic_full_sphere(self):
        pts = self._make_sphere_point_cloud(100, 1, 2, 3, 10.0)
        centre, r = geoprimitives.fitSphereAnalytic(pts)
        self.assertAlmostEqual(10.0, r, delta=1e-3)
        self.assertAlmostEqual(1, centre[0], delta=1e-3)
        self.assertAlmostEqual(2, centre[1], delta=1e-3)
        self.assertAlmostEqual(3, centre[2], delta=1e-3)

    def test_fit_sphere_analytic_hemisphere(self):
        pts = self._make_hemisphere_point_cloud(100, 1, 2, 3, 10.0)
        centre, r = geoprimitives.fitSphereAnalytic(pts)
        self.assertAlmostEqual(10.0, r, delta=1e-3)
        self.assertAlmostEqual(1, centre[0], delta=1e-3)
        self.assertAlmostEqual(2, centre[1], delta=1e-3)
        self.assertAlmostEqual(3, centre[2], delta=1e-3)
