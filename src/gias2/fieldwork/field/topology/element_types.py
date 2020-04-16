"""
FILE: element_types.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: 
Different types of elements and a constructor function to return
specified types of elements

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

# ======================================================================#
# element classes
# dimension: dimension of the element
# interior: tuple of tuples of the range of value define in the interior 
#           of the elements for each dimension
import logging

from numpy import array, cos, sin, eye, pi, sqrt, linspace, newaxis, all, any, bitwise_and, bitwise_or, hstack
from scipy.linalg import det

log = logging.getLogger(__name__)


class edge(object):
    """ element edge object.
    """

    def __init__(self, edge_num, start, end):
        """ edge number, start point tuple, end point tuple
        """
        self.edge_number = edge_num
        self.start = array(start)
        self.end = array(end)

    def get_elem_coord(self, x):
        """ x is the parameter along the edge 0 to 1
        """
        try:
            # array case
            return self.start * (1.0 - x)[:, newaxis] + self.end * x.T[:, newaxis]
        except TypeError:
            # float case
            return self.start * (1.0 - x) + self.end * x


class face(object):

    def __init__(self, face_num):
        pass


# ======================================================================#
class Element:
    """ Universal element class
    dimension: dimension of the element
    interior: tuple of tuples of the range of value define in the interior
              of the elements for each dimension

    number_of_points: need to be removed. This should be inferred from
    the field basis?
    """

    def __init__(self, dimensions, interior, number_of_points, type):
        self.dimensions = dimensions  # int
        self.interior = array(interior)  # ( (xi1.min, xi1.max), (xi2.min, xi2.max), ...)
        self.type = type  # string
        self.number_of_points = number_of_points  # int
        self.is_element = True
        self.edges = None
        self.edge_points = None

    def get_number_of_ensemble_points(self):
        return self.number_of_points

    def get_edge_points(self, edge):
        return self.edge_points[edge]

    def get_edge_coord(self, edge, x):
        return self.edges[edge].get_elem_coord(x)

    def generate_eval_grid(self, density):
        eval_coords = [linspace(self.interior[i][0], self.interior[i][1], density[i]) for i in range(self.dimensions)]
        grid = []

        if self.dimensions == 1:
            for i in eval_coords[0]:
                if self.is_interior([i]):
                    grid.append(i)
        elif self.dimensions == 2:
            for i in eval_coords[0]:
                for j in eval_coords[1]:
                    if self.is_interior([i, j]):
                        grid.append([i, j])
        elif self.dimensions == 3:
            for i in eval_coords[0]:
                for j in eval_coords[1]:
                    for k in eval_coords[2]:
                        if self.is_interior([i, j, k]):
                            grid.append([i, j, k])

        return array(grid)

    def get_point_edge(self, p):
        """ return the (edge number, point index, edge object) element
        point p is on
        """
        ret = []
        for ei, ep in enumerate(self.edge_points):
            try:
                ret.append((ei, ep.index(p), self.edges[ei]))
            except ValueError:
                pass
        return ret

    def get_edge_by_edge_points(self, p):
        """Given a sequence of element node numbers, return the edge that the
        nodes are on. Also return an integer i denoting the direction of the
        edge that matches the sequence (1 if along, -1 if opposite). Returns
        None, None if no match.
        """
        _p = tuple(p)
        for edge_i, edge_points in enumerate(self.edge_points):
            if _p == edge_points:
                return self.edges[edge_i], 1
            elif _p[::-1] == edge_points:
                return self.edges[edge_i], -1

        return None, None


# ======================================================================#
class line(Element):
    def is_interior(self, coords):

        if len(coords) != 1:
            raise ValueError(
                'wrong number of input coordinates. Need ' + str(self.dimensions) + ' , got ' + str(len(coords)))
        else:
            return all(bitwise_and((self.interior[:, 0] <= coords), coords <= self.interior[:, 1]))

    def is_boundary(self, coords):
        if len(coords) != 1:
            raise ValueError(
                'wrong number of input coordinates. Need ' + str(self.dimensions) + ' , got ' + str(len(coords)))
        else:
            return any(bitwise_or((self.interior[:, 0] == coords), coords == self.interior[:, 1]))


# ======================================================================#
class quad(Element):
    """ quadralateral elements
    """

    def is_interior(self, coords):

        if len(coords) != self.dimensions:
            raise ValueError(
                'wrong number of input coordinates. Need ' + str(self.dimensions) + ' , got ' + str(len(coords)))
        else:
            return all(bitwise_and((self.interior[:, 0] <= coords), coords <= self.interior[:, 1]))

    def is_boundary(self, coords):

        if len(coords) != self.dimensions:
            raise ValueError(
                'wrong number of input coordinates. Need ' + str(self.dimensions) + ' , got ' + str(len(coords)))
        else:
            return any(bitwise_or((self.interior[:, 0] == coords), coords == self.interior[:, 1]))


# ======================================================================#
class prism(Element):
    """ quadralateral elements
    """

    def is_interior(self, coords):

        if len(coords) != self.dimensions:
            raise ValueError(
                'wrong number of input coordinates. Need ' + str(self.dimensions) + ' , got ' + str(len(coords)))
        else:
            return (sum(coords[:2]) <= 1.0) and all(
                bitwise_and((self.interior[:, 0] <= coords), self.interior[:, 1] >= coords))

    def is_boundary(self, coords):

        if len(coords) != self.dimensions:
            raise ValueError(
                'wrong number of input coordinates. Need ' + str(self.dimensions) + ' , got ' + str(len(coords)))
        else:
            return self.is_interior(coords) and ((sum(coords[:2]) == 1.0) or any(
                bitwise_or((self.interior[:, 0] == coords), coords == self.interior[:, 1])))


# ======================================================================#
class tri(Element):
    """ simplex elements
    """
    O = eye(9, 9)

    def is_interior(self, coords):

        if len(coords) != self.dimensions:
            raise ValueError(
                'wrong number of input coordinates. Need ' + str(self.dimensions) + ' , got ' + str(len(coords)))
        else:
            return (sum(coords) <= 1.0) and all(
                bitwise_and((self.interior[:, 0] <= coords), coords <= self.interior[:, 1]))

    def is_boundary(self, coords):

        if len(coords) != self.dimensions:
            raise ValueError(
                'wrong number of input coordinates. Need ' + str(self.dimensions) + ' , got ' + str(len(coords)))
        else:
            return self.is_interior(coords) and ((sum(coords[:2]) == 1.0) or any(
                bitwise_or((self.interior[:, 0] == coords), coords == self.interior[:, 1])))

    def set_angle(self, a):
        """ calculates matrix for transforming derivatives to element
        orientation. Angles in degrees
        """
        # ~ a = a/180.0*pi
        a = a / 180.0 * pi
        self.O = eye(9, 9)
        self.O[1:3, 1:3] = array([[cos(a), sin(a)], [-sin(a), cos(a)]])
        self.O[4:6, 4:6] = array([[-sin(a + pi / 4.0), cos(a + pi / 4.0)], [-cos(a), -sin(a)]])
        self.O[7:9, 7:9] = array([[sin(a), -cos(a)], [sin(a + pi / 4.0), -cos(a + pi / 4.0)]])

        # ~ self.O[1:3, 1:3] = array( [[cos(a), sin(a)], [-sin(a), cos(a)]] )
        # ~ self.O[4:6, 4:6] = array( [[-sin(a), cos(a)], [-cos(a), -sin(a)]] )
        # ~ self.O[7:9, 7:9] = array( [[sin(a), -cos(a)], [sin(a), -cos(a)]] )

    def set_vertices(self, v1, v2, v3):

        # triangle area
        self.area = 0.5 * det([[1.0, v1[0], v1[1]],
                               [1.0, v2[0], v2[1]],
                               [1.0, v3[0], v3[1]]])

        # cartesian -> area coordinate coefficients
        self.a = array([v2[0] * v3[1] - v3[0] * v2[1],
                        v3[0] * v1[1] - v1[0] * v3[1],
                        v1[0] * v2[1] - v2[0] * v1[1]])
        self.b = array([v2[1] - v3[1],
                        v3[1] - v1[1],
                        v1[1] - v2[1]])
        self.c = array([v3[0] - v2[0],
                        v1[0] - v3[0],
                        v2[0] - v1[0]])

        self.v1 = v1
        self.v2 = v2
        self.v3 = v3

        return 1


# ======================================================================#
class triBezier(Element):
    """ simplex elements
    """
    O = eye(9, 9)
    v1 = (0.0, 0.0)
    v2 = (1.0, 0.0)
    v3 = (0.0, 1.0)
    # triangle area
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

    def __init__(self, dimensions, interior, number_of_points, type):
        self.dimensions = dimensions  # int
        self.interior = interior  # ( (xi1.min, xi1.max), (xi2.min, xi2.max), ...)
        self.type = type  # string
        self.number_of_points = number_of_points  # int
        self.is_element = True
        return

    def get_number_of_ensemble_points(self):
        return self.number_of_points

    def is_interior(self, coords):

        if len(coords) != self.dimensions:
            log.debug('ERROR: tri.is_interior: wrong number of input coordinates. Need', self.dimensions, ', got',
                  len(coords))
            return 0
        else:
            return sum(coords) <= 1.0

    def is_boundary(self, coords):

        if len(coords) != self.dimensions:
            log.debug('ERROR: tri.is_interior: wrong number of input coordinates. Need', self.dimensions, ', got',
                  len(coords))
            return 0
        else:
            return sum(coords) == 1.0

        # ======================================================================#


class tri_equilateral(Element):
    """ equilateral 2D simplex element
    """
    O = eye(9, 9)
    a = sqrt(3.0)

    def is_interior(self, coords):

        if len(coords) != 2:
            log.debug('ERROR: quad.is_interior: wrong number of input coordinates. Need', self.dimensions, ', got',
                  len(coords))
            return 0
        else:
            if coords[0] < 0.0:
                return coords[1] <= (coords[0] * self.a + self.a / 2.0)
            else:
                return coords[1] <= (-coords[0] * self.a + self.a / 2.0)

    def is_boundary(self, coords):
        t = 1.0e-10  # tolerance
        if len(coords) != 2:
            log.debug('ERROR: quad.is_interior: wrong number of input coordinates. Need', self.dimensions, ', got',
                  len(coords))
            return 0
        else:
            if coords[0] < 0.0:
                return abs(coords[1] - (coords[0] * self.a + self.a / 2.0)) < t
            else:
                return abs(coords[1] - (-coords[0] * self.a + self.a / 2.0)) < t

    def set_angle(self, a):
        """ calculates matrix for transforming derivatives to element
        orientation. Angles in degrees
        """
        a = a / 180.0 * pi
        self.O = eye(9, 9)
        self.O[1:3, 1:3] = array([[cos(a), sin(a)], [cos(a + pi / 3.0), sin(a + pi / 3.0)]])
        self.O[4:6, 4:6] = array([[-cos(pi / 3.0 - a), sin(pi / 3.0 - a)], [-cos(a), -sin(a)]])
        self.O[7:9, 7:9] = array([[-cos(pi / 3.0 + a), -sin(pi / 3.0 + a)], [cos(pi / 3.0 - a), -sin(pi / 3.0 - a)]])

        # ==================================================================#

    def _equi_evalx(self, y, b):
        """ calculates the x coord on the edge of an equilateral triangle
        give y
        """
        a = sqrt(3.0)
        return ((y - b) / a)

    # ==================================================================#
    def generate_eval_grid(self, d):
        """ generate element points for a equilateral triangle.
        Returns array( [[xcoords], [ycoords]] ).
        d: list of evaluation density in each direction
        int: interior bounds of the element ( (xmin, xmax), (ymin, ymax) )
        """

        e_points = array([[], []])
        y_divs = linspace(self.interior[1][0], self.interior[1][1], d[1])
        xn = d[0]
        for y_row in y_divs:
            x_edge = self._equi_evalx(y_row, self.interior[1][1])
            x = linspace(x_edge, -x_edge, xn)
            e_points = hstack([e_points, array([x, [y_row] * xn])])
            xn -= 1
        return e_points


# ======================================================================#
# ======================================================================#
element_types = {'line2l': line(1, ((0, 1),), 2, 'line2l'),
                 'line2h': line(1, ((0, 1),), 2, 'line2h'),
                 'line3l': line(1, ((0, 1),), 3, 'line3l'),
                 'line4l': line(1, ((0, 1),), 4, 'line4l'),
                 'line5l': line(1, ((0, 1),), 5, 'line5l'),
                 'tri3l': tri(2, ((0, 1), (0, 1)), 3, 'tri3l'),
                 'tri3h': tri(2, ((0, 1), (0, 1)), 3, 'tri3h'),
                 'tri6': tri(2, ((0, 1), (0, 1)), 6, 'tri6'),
                 'tri6_ql': tri(2, ((0, 1), (0, 1)), 6, 'tri6_ql'),
                 'tri10': tri(2, ((0, 1), (0, 1)), 10, 'tri10'),
                 'tri15': tri(2, ((0, 1), (0, 1)), 15, 'tri15'),
                 'tri16': tri(3, ((0, 1), (0, 1), (0, 1)), 16, 'tri16'),
                 'tri35': tri(3, ((0, 1), (0, 1), (0, 1)), 35, 'tri35'),
                 'quad22l': quad(2, ((0, 1), (0, 1)), 4, 'quad22l'),
                 'quad22h': quad(2, ((0, 1), (0, 1)), 4, 'quad22h'),
                 # 'quad23l':
                 # 'quad23h':
                 # 'quad32l':
                 # 'quad32h':\
                 'quad33': quad(2, ((0, 1), (0, 1)), 9, 'quad33'),
                 'quad44': quad(2, ((0, 1), (0, 1)), 16, 'quad44'),
                 'quad54': quad(2, ((0, 1), (0, 1)), 20, 'quad54'),
                 'quad55': quad(2, ((0, 1), (0, 1)), 25, 'quad55'),
                 'quad333': quad(3, ((0, 1), (0, 1), (0, 1)), 27, 'quad333'),
                 'quad444': quad(3, ((0, 1), (0, 1), (0, 1)), 64, 'quad444'),
                 'quad555': quad(3, ((0, 1), (0, 1), (0, 1)), 125, 'quad555'),
                 'quad542': quad(3, ((0, 1), (0, 1), (0, 1)), 40, 'quad542'),
                 'quad552': quad(3, ((0, 1), (0, 1), (0, 1)), 50, 'quad552'),
                 'prism15-5': prism(3, ((0, 1), (0, 1), (0, 1)), 75, 'prism15-5'),
                 'prism15-2': prism(3, ((0, 1), (0, 1), (0, 1)), 30, 'prism15-2'),
                 'prism6-5': prism(3, ((0, 1), (0, 1), (0, 1)), 30, 'prism6-5'),
                 'tri3e': tri_equilateral(2, ((-0.5, 0.5), (0, sqrt(3.0) / 2)), 3, 'tri3e'),
                 'tri10Bezier': triBezier(2, ((0.0, 1.0), (0.0, 1.0)), 9, 'tri10Bezier')}

element_node_coords = {'line3l': array([0.0, 0.5, 1.0]),

                       'line4l': array([0.0, 1. / 3, 2. / 3, 1.0]),

                       'line5l': array([0.0, 0.25, 0.5, 0.75, 1.0]),

                       'tri6': array([(0.0, 0.0), (0.5, 0.0), (1.0, 0.0), (0.5, 0.5), (0.0, 1.0),
                                      (0.0, 0.5)]),

                       'tri6_ql': array([(0.0, 0.0), (0.25, 0.0), (0.5, 0.0), (0.75, 0.0), (1.0, 0.0),
                                         (0.0, 1.0)]),

                       'tri10': array([(0.0, 0.0), (1 / 3.0, 0.0), (2 / 3.0, 0.0), (1.0, 0.0),
                                       (0.0, 1 / 3.0), (1 / 3.0, 1 / 3.0), (2 / 3.0, 1 / 3.0),
                                       (0.0, 2 / 3.0), (1 / 3.0, 2 / 3.0),
                                       (0.0, 1.0)]),

                       'tri15': array([(0.0, 0.0), (0.25, 0.0), (0.5, 0.0), (0.75, 0.0), (1.0, 0.0),
                                       (0.0, 0.25), (0.25, 0.25), (0.5, 0.25), (0.75, 0.25),
                                       (0.0, 0.5), (0.25, 0.5), (0.5, 0.5),
                                       (0.0, 0.75), (0.25, 0.75),
                                       (0.0, 1.0)]),

                       'tri16': array(
                           [(0.0, 0.0, 0.0), (0.25, 0.0, 0.0), (0.5, 0.0, 0.0), (0.75, 0.0, 0.0), (1.0, 0.0, 0.0),
                            (0.0, 0.25, 0.0), (0.25, 0.25, 0.0), (0.5, 0.25, 0.0), (0.75, 0.25, 0.0),
                            (0.0, 0.5, 0.0), (0.25, 0.5, 0.0), (0.5, 0.5, 0.0),
                            (0.0, 0.75, 0.0), (0.25, 0.75, 0.0),
                            (0.0, 1.0, 0.0),
                            (0.0, 0.0, 1.0)]),

                       'tri35': array(
                           [(0.0, 0.0, 0.0), (0.25, 0.0, 0.0), (0.5, 0.0, 0.0), (0.75, 0.0, 0.0), (1.0, 0.0, 0.0),
                            (0.0, 0.25, 0.0), (0.25, 0.25, 0.0), (0.5, 0.25, 0.0), (0.75, 0.25, 0.0),
                            (0.0, 0.5, 0.0), (0.25, 0.5, 0.0), (0.5, 0.5, 0.0),
                            (0.0, 0.75, 0.0), (0.25, 0.75, 0.0),
                            (0.0, 1.0, 0.0),
                            (0.0, 0.0, 0.25), (0.25, 0.0, 0.25), (0.5, 0.0, 0.25), (0.75, 0.0, 0.25),
                            (0.0, 0.25, 0.25), (0.25, 0.25, 0.25), (0.5, 0.25, 0.25),
                            (0.0, 0.5, 0.25), (0.25, 0.5, 0.25),
                            (0.0, 0.75, 0.25),
                            (0.0, 0.0, 0.5), (0.25, 0.0, 0.5), (0.5, 0.0, 0.5),
                            (0.0, 0.25, 0.5), (0.25, 0.25, 0.5),
                            (0.0, 0.5, 0.5),
                            (0.0, 0.0, 0.75), (0.25, 0.0, 0.75),
                            (0.0, 0.25, 0.75),
                            (0.0, 0.0, 1.0)]),

                       'quad33': array([(0.0, 0.0), (0.5, 0.0), (1.0, 0.0),
                                        (0.0, 0.5), (0.5, 0.5), (1.0, 0.5),
                                        (0.0, 1.0), (0.5, 1.0), (1.0, 1.0)]),

                       'quad44': array([(0.0, 0.0), (1 / 3.0, 0.0), (2 / 3.0, 0.0), (1.0, 0.0),
                                        (0.0, 1 / 3.0), (1 / 3.0, 1 / 3.0), (2 / 3.0, 1 / 3.0), (1.0, 1 / 3.0),
                                        (0.0, 2 / 3.0), (1 / 3.0, 2 / 3.0), (2 / 3.0, 2 / 3.0), (1.0, 2 / 3.0),
                                        (0.0, 1.0), (1 / 3.0, 1.0), (2 / 3.0, 1.0), (1.0, 1.0)]),

                       'quad54': array([(0.0, 0.0), (0.25, 0.0), (0.5, 0.0), (0.75, 0.0), (1.0, 0.0),
                                        (0.0, 1. / 3), (0.25, 1. / 3), (0.5, 1. / 3), (0.75, 1. / 3), (1.0, 1. / 3),
                                        (0.0, 2. / 3), (0.25, 2. / 3), (0.5, 2. / 3), (0.75, 2. / 3), (1.0, 2. / 3),
                                        (0.0, 1.0), (0.25, 1.0), (0.5, 1.0), (0.75, 1.0), (1.0, 1.0)]),

                       'quad55': array([(0.0, 0.0), (0.25, 0.0), (0.5, 0.0), (0.75, 0.0), (1.0, 0.0),
                                        (0.0, 0.25), (0.25, 0.25), (0.5, 0.25), (0.75, 0.25), (1.0, 0.25),
                                        (0.0, 0.5), (0.25, 0.5), (0.5, 0.5), (0.75, 0.5), (1.0, 0.5),
                                        (0.0, 0.75), (0.25, 0.75), (0.5, 0.75), (0.75, 0.75), (1.0, 0.75),
                                        (0.0, 1.0), (0.25, 1.0), (0.5, 1.0), (0.75, 1.0), (1.0, 1.0)]),

                       'quad333': array([(0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0),
                                         (0.0, 0.5, 0.0), (0.5, 0.5, 0.0), (1.0, 0.5, 0.0),
                                         (0.0, 1.0, 0.0), (0.5, 1.0, 0.0), (1.0, 1.0, 0.0),

                                         (0.0, 0.0, 0.5), (0.5, 0.0, 0.5), (1.0, 0.0, 0.5),
                                         (0.0, 0.5, 0.5), (0.5, 0.5, 0.5), (1.0, 0.5, 0.5),
                                         (0.0, 1.0, 0.5), (0.5, 1.0, 0.5), (1.0, 1.0, 0.5),

                                         (0.0, 0.0, 1.0), (0.5, 0.0, 1.0), (1.0, 0.0, 1.0),
                                         (0.0, 0.5, 1.0), (0.5, 0.5, 1.0), (1.0, 0.5, 1.0),
                                         (0.0, 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 1.0, 1.0)]),

                       'quad444': array([(0.0, 0.0, 0.0), (1 / 3.0, 0.0, 0.0), (2 / 3.0, 0.0, 0.0), (1.0, 0.0, 0.0),
                                         (0.0, 1 / 3.0, 0.0), (1 / 3.0, 1 / 3.0, 0.0), (2 / 3.0, 1 / 3.0, 0.0),
                                         (1.0, 1 / 3.0, 0.0),
                                         (0.0, 2 / 3.0, 0.0), (1 / 3.0, 2 / 3.0, 0.0), (2 / 3.0, 2 / 3.0, 0.0),
                                         (1.0, 2 / 3.0, 0.0),
                                         (0.0, 1.0, 0.0), (1 / 3.0, 1.0, 0.0), (2 / 3.0, 1.0, 0.0), (1.0, 1.0, 0.0),

                                         (0.0, 0.0, 1 / 3.0), (1 / 3.0, 0.0, 1 / 3.0), (2 / 3.0, 0.0, 1 / 3.0),
                                         (1.0, 0.0, 1 / 3.0),
                                         (0.0, 1 / 3.0, 1 / 3.0), (1 / 3.0, 1 / 3.0, 1 / 3.0),
                                         (2 / 3.0, 1 / 3.0, 1 / 3.0), (1.0, 1 / 3.0, 1 / 3.0),
                                         (0.0, 2 / 3.0, 1 / 3.0), (1 / 3.0, 2 / 3.0, 1 / 3.0),
                                         (2 / 3.0, 2 / 3.0, 1 / 3.0), (1.0, 2 / 3.0, 1 / 3.0),
                                         (0.0, 1.0, 1 / 3.0), (1 / 3.0, 1.0, 1 / 3.0), (2 / 3.0, 1.0, 1 / 3.0),
                                         (1.0, 1.0, 1 / 3.0),

                                         (0.0, 0.0, 2 / 3.0), (1 / 3.0, 0.0, 2 / 3.0), (2 / 3.0, 0.0, 2 / 3.0),
                                         (1.0, 0.0, 2 / 3.0),
                                         (0.0, 1 / 3.0, 2 / 3.0), (1 / 3.0, 1 / 3.0, 2 / 3.0),
                                         (2 / 3.0, 1 / 3.0, 2 / 3.0), (1.0, 1 / 3.0, 2 / 3.0),
                                         (0.0, 2 / 3.0, 2 / 3.0), (1 / 3.0, 2 / 3.0, 2 / 3.0),
                                         (2 / 3.0, 2 / 3.0, 2 / 3.0), (1.0, 2 / 3.0, 2 / 3.0),
                                         (0.0, 1.0, 2 / 3.0), (1 / 3.0, 1.0, 2 / 3.0), (2 / 3.0, 1.0, 2 / 3.0),
                                         (1.0, 1.0, 2 / 3.0),

                                         (0.0, 0.0, 1.0), (1 / 3.0, 0.0, 1.0), (2 / 3.0, 0.0, 1.0), (1.0, 0.0, 1.0),
                                         (0.0, 1 / 3.0, 1.0), (1 / 3.0, 1 / 3.0, 1.0), (2 / 3.0, 1 / 3.0, 1.0),
                                         (1.0, 1 / 3.0, 1.0),
                                         (0.0, 2 / 3.0, 1.0), (1 / 3.0, 2 / 3.0, 1.0), (2 / 3.0, 2 / 3.0, 1.0),
                                         (1.0, 2 / 3.0, 1.0),
                                         (0.0, 1.0, 1.0), (1 / 3.0, 1.0, 1.0), (2 / 3.0, 1.0, 1.0), (1.0, 1.0, 1.0)]),

                       'quad555': array(
                           [(0.0, 0.0, 0.0), (0.25, 0.0, 0.0), (0.5, 0.0, 0.0), (0.75, 0.0, 0.0), (1.0, 0.0, 0.0),
                            (0.0, 0.25, 0.0), (0.25, 0.25, 0.0), (0.5, 0.25, 0.0), (0.75, 0.25, 0.0), (1.0, 0.25, 0.0),
                            (0.0, 0.5, 0.0), (0.25, 0.5, 0.0), (0.5, 0.5, 0.0), (0.75, 0.5, 0.0), (1.0, 0.5, 0.0),
                            (0.0, 0.75, 0.0), (0.25, 0.75, 0.0), (0.5, 0.75, 0.0), (0.75, 0.75, 0.0), (1.0, 0.75, 0.0),
                            (0.0, 1.0, 0.0), (0.25, 1.0, 0.0), (0.5, 1.0, 0.0), (0.75, 1.0, 0.0), (1.0, 1.0, 0.0),

                            (0.0, 0.0, 0.25), (0.25, 0.0, 0.25), (0.5, 0.0, 0.25), (0.75, 0.0, 0.25), (1.0, 0.0, 0.25),
                            (0.0, 0.25, 0.25), (0.25, 0.25, 0.25), (0.5, 0.25, 0.25), (0.75, 0.25, 0.25),
                            (1.0, 0.25, 0.25),
                            (0.0, 0.5, 0.25), (0.25, 0.5, 0.25), (0.5, 0.5, 0.25), (0.75, 0.5, 0.25), (1.0, 0.5, 0.25),
                            (0.0, 0.75, 0.25), (0.25, 0.75, 0.25), (0.5, 0.75, 0.25), (0.75, 0.75, 0.25),
                            (1.0, 0.75, 0.25),
                            (0.0, 1.0, 0.25), (0.25, 1.0, 0.25), (0.5, 1.0, 0.25), (0.75, 1.0, 0.25), (1.0, 1.0, 0.25),

                            (0.0, 0.0, 0.5), (0.25, 0.0, 0.5), (0.5, 0.0, 0.5), (0.75, 0.0, 0.5), (1.0, 0.0, 0.5),
                            (0.0, 0.25, 0.5), (0.25, 0.25, 0.5), (0.5, 0.25, 0.5), (0.75, 0.25, 0.5), (1.0, 0.25, 0.5),
                            (0.0, 0.5, 0.5), (0.25, 0.5, 0.5), (0.5, 0.5, 0.5), (0.75, 0.5, 0.5), (1.0, 0.5, 0.5),
                            (0.0, 0.75, 0.5), (0.25, 0.75, 0.5), (0.5, 0.75, 0.5), (0.75, 0.75, 0.5), (1.0, 0.75, 0.5),
                            (0.0, 1.0, 0.5), (0.25, 1.0, 0.5), (0.5, 1.0, 0.5), (0.75, 1.0, 0.5), (1.0, 1.0, 0.5),

                            (0.0, 0.0, 0.75), (0.25, 0.0, 0.75), (0.5, 0.0, 0.75), (0.75, 0.0, 0.75), (1.0, 0.0, 0.75),
                            (0.0, 0.25, 0.75), (0.25, 0.25, 0.75), (0.5, 0.25, 0.75), (0.75, 0.25, 0.75),
                            (1.0, 0.25, 0.75),
                            (0.0, 0.5, 0.75), (0.25, 0.5, 0.75), (0.5, 0.5, 0.75), (0.75, 0.5, 0.75), (1.0, 0.5, 0.75),
                            (0.0, 0.75, 0.75), (0.25, 0.75, 0.75), (0.5, 0.75, 0.75), (0.75, 0.75, 0.75),
                            (1.0, 0.75, 0.75),
                            (0.0, 1.0, 0.75), (0.25, 1.0, 0.75), (0.5, 1.0, 0.75), (0.75, 1.0, 0.75), (1.0, 1.0, 0.75),

                            (0.0, 0.0, 1.0), (0.25, 0.0, 1.0), (0.5, 0.0, 1.0), (0.75, 0.0, 1.0), (1.0, 0.0, 1.0),
                            (0.0, 0.25, 1.0), (0.25, 0.25, 1.0), (0.5, 0.25, 1.0), (0.75, 0.25, 1.0), (1.0, 0.25, 1.0),
                            (0.0, 0.5, 1.0), (0.25, 0.5, 1.0), (0.5, 0.5, 1.0), (0.75, 0.5, 1.0), (1.0, 0.5, 1.0),
                            (0.0, 0.75, 1.0), (0.25, 0.75, 1.0), (0.5, 0.75, 1.0), (0.75, 0.75, 1.0), (1.0, 0.75, 1.0),
                            (0.0, 1.0, 1.0), (0.25, 1.0, 1.0), (0.5, 1.0, 1.0), (0.75, 1.0, 1.0), (1.0, 1.0, 1.0),
                            ]),

                       'quad542': array(
                           [(0.0, 0.0, 0.0), (0.25, 0.0, 0.0), (0.5, 0.0, 0.0), (0.75, 0.0, 0.0), (1.0, 0.0, 0.0),
                            (0.0, 1 / 3.0, 0.0), (0.25, 1 / 3.0, 0.0), (0.5, 1 / 3.0, 0.0), (0.75, 1 / 3.0, 0.0),
                            (1.0, 1 / 3.0, 0.0),
                            (0.0, 2 / 3.0, 0.0), (0.25, 2 / 3.0, 0.0), (0.5, 2 / 3.0, 0.0), (0.75, 2 / 3.0, 0.0),
                            (1.0, 2 / 3.0, 0.0),
                            (0.0, 1.0, 0.0), (0.25, 1.0, 0.0), (0.5, 1.0, 0.0), (0.75, 1.0, 0.0), (1.0, 1.0, 0.0),

                            (0.0, 0.0, 1.0), (0.25, 0.0, 1.0), (0.5, 0.0, 1.0), (0.75, 0.0, 1.0), (1.0, 0.0, 1.0),
                            (0.0, 1 / 3.0, 1.0), (0.25, 1 / 3.0, 1.0), (0.5, 1 / 3.0, 1.0), (0.75, 1 / 3.0, 1.0),
                            (1.0, 1 / 3.0, 1.0),
                            (0.0, 2 / 3.0, 1.0), (0.25, 2 / 3.0, 1.0), (0.5, 2 / 3.0, 1.0), (0.75, 2 / 3.0, 1.0),
                            (1.0, 2 / 3.0, 1.0),
                            (0.0, 1.0, 1.0), (0.25, 1.0, 1.0), (0.5, 1.0, 1.0), (0.75, 1.0, 1.0), (1.0, 1.0, 1.0),
                            ]),

                       'quad552': array(
                           [(0.0, 0.0, 0.0), (0.25, 0.0, 0.0), (0.5, 0.0, 0.0), (0.75, 0.0, 0.0), (1.0, 0.0, 0.0),
                            (0.0, 0.25, 0.0), (0.25, 0.25, 0.0), (0.5, 0.25, 0.0), (0.75, 0.25, 0.0), (1.0, 0.25, 0.0),
                            (0.0, 0.5, 0.0), (0.25, 0.5, 0.0), (0.5, 0.5, 0.0), (0.75, 0.5, 0.0), (1.0, 0.5, 0.0),
                            (0.0, 0.75, 0.0), (0.25, 0.75, 0.0), (0.5, 0.75, 0.0), (0.75, 0.75, 0.0), (1.0, 0.75, 0.0),
                            (0.0, 1.0, 0.0), (0.25, 1.0, 0.0), (0.5, 1.0, 0.0), (0.75, 1.0, 0.0), (1.0, 1.0, 0.0),

                            (0.0, 0.0, 1.0), (0.25, 0.0, 1.0), (0.5, 0.0, 1.0), (0.75, 0.0, 1.0), (1.0, 0.0, 1.0),
                            (0.0, 0.25, 1.0), (0.25, 0.25, 1.0), (0.5, 0.25, 1.0), (0.75, 0.25, 1.0), (1.0, 0.25, 1.0),
                            (0.0, 0.5, 1.0), (0.25, 0.5, 1.0), (0.5, 0.5, 1.0), (0.75, 0.5, 1.0), (1.0, 0.5, 1.0),
                            (0.0, 0.75, 1.0), (0.25, 0.75, 1.0), (0.5, 0.75, 1.0), (0.75, 0.75, 1.0), (1.0, 0.75, 1.0),
                            (0.0, 1.0, 1.0), (0.25, 1.0, 1.0), (0.5, 1.0, 1.0), (0.75, 1.0, 1.0), (1.0, 1.0, 1.0),
                            ]),

                       'prism15-5': array(
                           [(0.0, 0.0, 0.0), (0.25, 0.0, 0.0), (0.5, 0.0, 0.0), (0.75, 0.0, 0.0), (1.0, 0.0, 0.0),
                            (0.0, 0.25, 0.0), (0.25, 0.25, 0.0), (0.5, 0.25, 0.0), (0.75, 0.25, 0.0),
                            (0.0, 0.5, 0.0), (0.25, 0.5, 0.0), (0.5, 0.5, 0.0),
                            (0.0, 0.75, 0.0), (0.25, 0.75, 0.0),
                            (0.0, 1.0, 0.0),
                            (0.0, 0.0, 0.25), (0.25, 0.0, 0.25), (0.5, 0.0, 0.25), (0.75, 0.0, 0.25), (1.0, 0.0, 0.25),
                            (0.0, 0.25, 0.25), (0.25, 0.25, 0.25), (0.5, 0.25, 0.25), (0.75, 0.25, 0.25),
                            (0.0, 0.5, 0.25), (0.25, 0.5, 0.25), (0.5, 0.5, 0.25),
                            (0.0, 0.75, 0.25), (0.25, 0.75, 0.25),
                            (0.0, 1.0, 0.25),
                            (0.0, 0.0, 0.5), (0.25, 0.0, 0.5), (0.5, 0.0, 0.5), (0.75, 0.0, 0.5), (1.0, 0.0, 0.5),
                            (0.0, 0.25, 0.5), (0.25, 0.25, 0.5), (0.5, 0.25, 0.5), (0.75, 0.25, 0.5),
                            (0.0, 0.5, 0.5), (0.25, 0.5, 0.5), (0.5, 0.5, 0.5),
                            (0.0, 0.75, 0.5), (0.25, 0.75, 0.5),
                            (0.0, 1.0, 0.5),
                            (0.0, 0.0, 0.75), (0.25, 0.0, 0.75), (0.5, 0.0, 0.75), (0.75, 0.0, 0.75), (1.0, 0.0, 0.75),
                            (0.0, 0.25, 0.75), (0.25, 0.25, 0.75), (0.5, 0.25, 0.75), (0.75, 0.25, 0.75),
                            (0.0, 0.5, 0.75), (0.25, 0.5, 0.75), (0.5, 0.5, 0.75),
                            (0.0, 0.75, 0.75), (0.25, 0.75, 0.75),
                            (0.0, 1.0, 0.75),
                            (0.0, 0.0, 1.0), (0.25, 0.0, 1.0), (0.5, 0.0, 1.0), (0.75, 0.0, 1.0), (1.0, 0.0, 1.0),
                            (0.0, 0.25, 1.0), (0.25, 0.25, 1.0), (0.5, 0.25, 1.0), (0.75, 0.25, 1.0),
                            (0.0, 0.5, 1.0), (0.25, 0.5, 1.0), (0.5, 0.5, 1.0),
                            (0.0, 0.75, 1.0), (0.25, 0.75, 1.0),
                            (0.0, 1.0, 1.0),
                            ]),

                       'prism15-2': array(
                           [(0.0, 0.0, 0.0), (0.25, 0.0, 0.0), (0.5, 0.0, 0.0), (0.75, 0.0, 0.0), (1.0, 0.0, 0.0),
                            (0.0, 0.25, 0.0), (0.25, 0.25, 0.0), (0.5, 0.25, 0.0), (0.75, 0.25, 0.0),
                            (0.0, 0.5, 0.0), (0.25, 0.5, 0.0), (0.5, 0.5, 0.0),
                            (0.0, 0.75, 0.0), (0.25, 0.75, 0.0),
                            (0.0, 1.0, 0.0),
                            (0.0, 0.0, 1.0), (0.25, 0.0, 1.0), (0.5, 0.0, 1.0), (0.75, 0.0, 1.0), (1.0, 0.0, 1.0),
                            (0.0, 0.25, 1.0), (0.25, 0.25, 1.0), (0.5, 0.25, 1.0), (0.75, 0.25, 1.0),
                            (0.0, 0.5, 1.0), (0.25, 0.5, 1.0), (0.5, 0.5, 1.0),
                            (0.0, 0.75, 1.0), (0.25, 0.75, 1.0),
                            (0.0, 1.0, 1.0),
                            ]),

                       'prism6-5': array(
                           [(0.0, 0.0, 0.0), (0.25, 0.0, 0.0), (0.5, 0.0, 0.0), (0.75, 0.0, 0.0), (1.0, 0.0, 0.0),
                            (0.0, 1.0, 0.0),
                            (0.0, 0.0, 0.25), (0.25, 0.0, 0.25), (0.5, 0.0, 0.25), (0.75, 0.0, 0.25), (1.0, 0.0, 0.25),
                            (0.0, 1.0, 0.25),
                            (0.0, 0.0, 0.5), (0.25, 0.0, 0.5), (0.5, 0.0, 0.5), (0.75, 0.0, 0.5), (1.0, 0.0, 0.5),
                            (0.0, 1.0, 0.5),
                            (0.0, 0.0, 0.75), (0.25, 0.0, 0.75), (0.5, 0.0, 0.75), (0.75, 0.0, 0.75), (1.0, 0.0, 0.75),
                            (0.0, 1.0, 0.75),
                            (0.0, 0.0, 1.0), (0.25, 0.0, 1.0), (0.5, 0.0, 1.0), (0.75, 0.0, 1.0), (1.0, 0.0, 1.0),
                            (0.0, 1.0, 1.0),
                            ])
                       }

# ~ element_edge_points = { 'tri6': ((0,1,2),(2,3,4),(4,5,0)),
element_edge_points = {'line3l': ((0,), (2,)),
                       'line4l': ((0,), (3,)),
                       'line5l': ((0,), (4,)),
                       'tri6': ((0, 1, 2), (2, 3, 4), (4, 5, 0)),
                       'tri6_ql': ((0, 1, 2, 3, 4), (4, 5), (5, 0)),
                       'tri10': ((0, 1, 2, 3), (3, 6, 8, 9), (9, 7, 4, 0)),
                       'tri15': ((0, 1, 2, 3, 4), (4, 8, 11, 13, 14), (14, 12, 9, 5, 0)),
                       'tri16': ((0, 1, 2, 3, 4), (4, 8, 11, 13, 14), (14, 12, 9, 5, 0), (0, 15), (4, 15), (14, 15)),
                       'tri35': ((0, 1, 2, 3, 4), (4, 8, 11, 13, 14), (14, 12, 9, 5, 0),
                                 (0, 15, 25, 31, 34), (4, 18, 27, 32, 34), (14, 24, 30, 33, 34)),
                       'quad33': ((0, 1, 2), (2, 5, 8), (8, 7, 6), (6, 3, 0)),  # loop
                       # ~ 'quad44': ((0,1,2,3),(3,7,11,15),(12,13,14,15),(0,4,8,12)),
                       'quad44': ((0, 1, 2, 3), (3, 7, 11, 15), (15, 14, 13, 12), (12, 8, 4, 0)),  # loop
                       # ~ 'quad54': ((0,1,2,3,4),(4,9,14,19),(15,16,17,18,19),(0,5,10,15)),
                       'quad54': ((0, 1, 2, 3, 4), (4, 9, 14, 19), (19, 18, 17, 16, 15), (15, 10, 5, 0)),  # loop
                       # ~ 'quad55': ((0,1,2,3,4),(4,9,14,19,24),(20,21,22,23,24),(0,5,10,15,20)),
                       'quad55': ((0, 1, 2, 3, 4), (4, 9, 14, 19, 24), (24, 23, 22, 21, 20), (20, 15, 10, 5, 0)),
                       # loop
                       'quad333': ((0, 1, 2), (2, 5, 8), (6, 7, 8), (0, 3, 6),
                                   (0, 9, 18), (2, 11, 20), (8, 17, 26), (6, 15, 24),
                                   (18, 19, 20), (20, 23, 26), (24, 25, 26), (18, 21, 24)),
                       'quad444': ((0, 1, 2, 3), (3, 7, 11, 15), (12, 13, 14, 15), (0, 4, 8, 12),
                                   (0, 16, 32, 48), (3, 19, 35, 51), (15, 31, 47, 63), (12, 28, 44, 60),
                                   (48, 49, 50, 51), (51, 55, 59, 63), (60, 61, 62, 63), (48, 52, 56, 60)),
                       'quad555': ((0, 1, 2, 3, 4), (4, 9, 14, 19, 24), (20, 21, 22, 23, 24), (0, 5, 10, 15, 20),
                                   (0, 25, 50, 75, 100), (4, 29, 54, 79, 104), (24, 49, 74, 99, 124),
                                   (20, 45, 70, 95, 120),
                                   (100, 101, 102, 103, 104), (104, 109, 114, 119, 124), (120, 121, 122, 123, 124),
                                   (100, 105, 110, 115, 120)),
                       'quad542': ((0, 1, 2, 3, 4), (4, 9, 14, 19), (15, 16, 17, 18, 19), (0, 5, 10, 15),
                                   (20, 21, 22, 23, 24), (24, 29, 34, 39), (35, 36, 37, 38, 39), (20, 25, 30, 35)),
                       'quad552': ((0, 1, 2, 3, 4), (4, 9, 14, 19, 24), (20, 21, 22, 23, 24), (0, 5, 10, 15, 20),
                                   (25, 26, 27, 28, 29), (29, 34, 39, 44, 49), (45, 46, 47, 48, 49),
                                   (25, 30, 35, 40, 45)),
                       'prism15-5': ((0, 1, 2, 3, 4), (4, 8, 11, 13, 14), (14, 12, 9, 5, 0),
                                     (0, 15, 30, 45, 60), (4, 19, 34, 49, 64), (14, 29, 44, 59, 74),
                                     (60, 61, 62, 63, 64), (64, 68, 71, 73, 74), (74, 72, 69, 65, 60)),
                       'prism15-2': ((0, 1, 2, 3, 4), (4, 8, 11, 13, 14), (14, 12, 9, 5, 0),
                                     (0, 15), (4, 19), (14, 29),
                                     (15, 16, 17, 18, 19), (19, 23, 26, 28, 29), (15, 20, 24, 27, 29)),
                       'prism6-5': ((0, 1, 2, 3, 4), (4, 5), (5, 0),
                                    (0, 6, 12, 18, 24), (4, 10, 16, 22, 28), (5, 11, 17, 23, 29),
                                    (24, 25, 26, 27, 28), (28, 29), (29, 0)),

                       }
# ~ 'tri10': ((3,6,8,9),(9,7,4,0),(0,1,2,3)) }

element_edges = {'line3l': None,
                 'line4l': None,
                 'line5l': None,
                 # 'tri6':  ( edge( 2, (0.0,0.0), (1.0,0.0) ), edge( 0, (1.0,0.0), (0.0,1.0) ), edge( 1, (0.0,1.0), (0.0,0.0) ) ),
                 'tri6': (
                 edge(0, (0.0, 0.0), (1.0, 0.0)), edge(1, (1.0, 0.0), (0.0, 1.0)), edge(2, (0.0, 1.0), (0.0, 0.0))),
                 'tri6_ql': (
                 edge(0, (0.0, 0.0), (1.0, 0.0)), edge(1, (1.0, 0.0), (0.0, 1.0)), edge(2, (0.0, 1.0), (0.0, 0.0))),
                 'tri10': (
                 edge(0, (0.0, 0.0), (1.0, 0.0)), edge(1, (1.0, 0.0), (0.0, 1.0)), edge(2, (0.0, 1.0), (0.0, 0.0))),
                 'tri15': (
                 edge(0, (0.0, 0.0), (1.0, 0.0)), edge(1, (1.0, 0.0), (0.0, 1.0)), edge(2, (0.0, 1.0), (0.0, 0.0))),
                 'tri16': (
                 edge(0, (0.0, 0.0), (1.0, 0.0)), edge(1, (1.0, 0.0), (0.0, 1.0)), edge(2, (0.0, 1.0), (0.0, 0.0))),
                 'tri35': (edge(0, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)), edge(1, (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
                           edge(2, (0.0, 1.0, 0.0), (0.0, 0.0, 0.0)),
                           edge(3, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)), edge(4, (1.0, 0.0, 0.0), (0.0, 0.0, 1.0)),
                           edge(5, (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))),
                 'quad33': (
                 edge(0, (0.0, 0.0), (1.0, 0.0)), edge(1, (1.0, 0.0), (1.0, 1.0)), edge(2, (1.0, 1.0), (0.0, 1.0)),
                 edge(3, (0.0, 1.0), (0.0, 0.0))),  # loop
                 # ~ 'quad44': ( edge( 0, (0.0,0.0), (1.0,0.0) ), edge( 1, (1.0,0.0), (1.0,1.0) ), edge( 2, (0.0,1.0), (1.0,1.0) ), edge( 3, (0.0,0.0), (0.0,1.0) ) ),
                 'quad44': (
                 edge(0, (0.0, 0.0), (1.0, 0.0)), edge(1, (1.0, 0.0), (1.0, 1.0)), edge(2, (1.0, 1.0), (0.0, 1.0)),
                 edge(3, (0.0, 1.0), (0.0, 0.0))),  # loop
                 # ~ 'quad54': ( edge( 0, (0.0,0.0), (1.0,0.0) ), edge( 1, (1.0,0.0), (1.0,1.0) ), edge( 2, (0.0,1.0), (1.0,1.0) ), edge( 3, (0.0,0.0), (0.0,1.0) ) ),
                 'quad54': (
                 edge(0, (0.0, 0.0), (1.0, 0.0)), edge(1, (1.0, 0.0), (1.0, 1.0)), edge(2, (1.0, 1.0), (0.0, 1.0)),
                 edge(3, (0.0, 1.0), (0.0, 0.0))),  # loop
                 # ~ 'quad55': ( edge( 0, (0.0,0.0), (1.0,0.0) ), edge( 1, (1.0,0.0), (1.0,1.0) ), edge( 2, (0.0,1.0), (1.0,1.0) ), edge( 3, (0.0,0.0), (0.0,1.0) ) ),
                 'quad55': (
                 edge(0, (0.0, 0.0), (1.0, 0.0)), edge(1, (1.0, 0.0), (1.0, 1.0)), edge(2, (1.0, 1.0), (0.0, 1.0)),
                 edge(3, (0.0, 1.0), (0.0, 0.0))),  # loop
                 'quad333': (edge(0, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)), edge(1, (1.0, 0.0, 0.0), (1.0, 1.0, 0.0)),
                             edge(2, (0.0, 1.0, 0.0), (1.0, 1.0, 0.0)), edge(3, (0.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
                             edge(4, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)), edge(5, (1.0, 0.0, 0.0), (1.0, 0.0, 1.0)),
                             edge(6, (1.0, 1.0, 0.0), (1.0, 1.0, 1.0)), edge(7, (0.0, 1.0, 0.0), (0.0, 1.0, 1.0)),
                             edge(8, (0.0, 0.0, 1.0), (1.0, 0.0, 1.0)), edge(9, (1.0, 0.0, 1.0), (1.0, 1.0, 1.0)),
                             edge(10, (0.0, 1.0, 1.0), (1.0, 1.0, 1.0)), edge(11, (0.0, 0.0, 1.0), (0.0, 1.0, 1.0)),
                             ),
                 'quad444': (edge(0, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)), edge(1, (1.0, 0.0, 0.0), (1.0, 1.0, 0.0)),
                             edge(2, (0.0, 1.0, 0.0), (1.0, 1.0, 0.0)), edge(3, (0.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
                             edge(4, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)), edge(5, (1.0, 0.0, 0.0), (1.0, 0.0, 1.0)),
                             edge(6, (1.0, 1.0, 0.0), (1.0, 1.0, 1.0)), edge(7, (0.0, 1.0, 0.0), (0.0, 1.0, 1.0)),
                             edge(8, (0.0, 0.0, 1.0), (1.0, 0.0, 1.0)), edge(9, (1.0, 0.0, 1.0), (1.0, 1.0, 1.0)),
                             edge(10, (0.0, 1.0, 1.0), (1.0, 1.0, 1.0)), edge(11, (0.0, 0.0, 1.0), (0.0, 1.0, 1.0)),
                             ),
                 'quad555': (edge(0, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)), edge(1, (1.0, 0.0, 0.0), (1.0, 1.0, 0.0)),
                             edge(2, (0.0, 1.0, 0.0), (1.0, 1.0, 0.0)), edge(3, (0.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
                             edge(4, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)), edge(5, (1.0, 0.0, 0.0), (1.0, 0.0, 1.0)),
                             edge(6, (1.0, 1.0, 0.0), (1.0, 1.0, 1.0)), edge(7, (0.0, 1.0, 0.0), (0.0, 1.0, 1.0)),
                             edge(8, (0.0, 0.0, 1.0), (1.0, 0.0, 1.0)), edge(9, (1.0, 0.0, 1.0), (1.0, 1.0, 1.0)),
                             edge(10, (0.0, 1.0, 1.0), (1.0, 1.0, 1.0)), edge(11, (0.0, 0.0, 1.0), (0.0, 1.0, 1.0)),
                             ),
                 'quad542': (edge(0, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)), edge(1, (1.0, 0.0, 0.0), (1.0, 1.0, 0.0)),
                             edge(2, (0.0, 1.0, 0.0), (1.0, 1.0, 0.0)), edge(3, (0.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
                             edge(4, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)), edge(5, (1.0, 0.0, 0.0), (1.0, 0.0, 1.0)),
                             edge(6, (1.0, 1.0, 0.0), (1.0, 1.0, 1.0)), edge(7, (0.0, 1.0, 0.0), (0.0, 1.0, 1.0)),
                             edge(8, (0.0, 0.0, 1.0), (1.0, 0.0, 1.0)), edge(9, (1.0, 0.0, 1.0), (1.0, 1.0, 1.0)),
                             edge(10, (0.0, 1.0, 1.0), (1.0, 1.0, 1.0)), edge(11, (0.0, 0.0, 1.0), (0.0, 1.0, 1.0)),
                             ),
                 'quad552': (edge(0, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)), edge(1, (1.0, 0.0, 0.0), (1.0, 1.0, 0.0)),
                             edge(2, (0.0, 1.0, 0.0), (1.0, 1.0, 0.0)), edge(3, (0.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
                             edge(4, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)), edge(5, (1.0, 0.0, 0.0), (1.0, 0.0, 1.0)),
                             edge(6, (1.0, 1.0, 0.0), (1.0, 1.0, 1.0)), edge(7, (0.0, 1.0, 0.0), (0.0, 1.0, 1.0)),
                             edge(8, (0.0, 0.0, 1.0), (1.0, 0.0, 1.0)), edge(9, (1.0, 0.0, 1.0), (1.0, 1.0, 1.0)),
                             edge(10, (0.0, 1.0, 1.0), (1.0, 1.0, 1.0)), edge(11, (0.0, 0.0, 1.0), (0.0, 1.0, 1.0)),
                             ),
                 'prism15-5': (edge(0, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)), edge(1, (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
                               edge(2, (0.0, 1.0, 0.0), (0.0, 0.0, 0.0)),
                               edge(3, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)), edge(4, (1.0, 0.0, 0.0), (1.0, 0.0, 1.0)),
                               edge(5, (0.0, 1.0, 0.0), (0.0, 1.0, 1.0)),
                               edge(6, (0.0, 0.0, 1.0), (1.0, 0.0, 1.0)), edge(7, (1.0, 0.0, 1.0), (0.0, 1.0, 1.0)),
                               edge(8, (0.0, 1.0, 1.0), (0.0, 0.0, 1.0)),
                               ),
                 'prism15-2': (edge(0, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)), edge(1, (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
                               edge(2, (0.0, 1.0, 0.0), (0.0, 0.0, 0.0)),
                               edge(3, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)), edge(4, (1.0, 0.0, 0.0), (1.0, 0.0, 1.0)),
                               edge(5, (0.0, 1.0, 0.0), (0.0, 1.0, 1.0)),
                               edge(6, (0.0, 0.0, 1.0), (1.0, 0.0, 1.0)), edge(7, (1.0, 0.0, 1.0), (0.0, 1.0, 1.0)),
                               edge(8, (0.0, 1.0, 1.0), (0.0, 0.0, 1.0)),
                               ),
                 'prism6-5': (edge(0, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)), edge(1, (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
                              edge(2, (0.0, 1.0, 0.0), (0.0, 0.0, 0.0)),
                              edge(3, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)), edge(4, (1.0, 0.0, 0.0), (1.0, 0.0, 1.0)),
                              edge(5, (0.0, 1.0, 0.0), (0.0, 1.0, 1.0)),
                              edge(6, (0.0, 0.0, 1.0), (1.0, 0.0, 1.0)), edge(7, (1.0, 0.0, 1.0), (0.0, 1.0, 1.0)),
                              edge(8, (0.0, 1.0, 1.0), (0.0, 0.0, 1.0)),
                              ),

                 }


# ~ 'tri10': [ edge( 0, (1.0,0.0), (0.0,1.0) ), edge( 1, (0.0,1.0), (0.0,0.0) ), edge( 2, (0.0,0.0), (1.0,0.0) ) ] }

# ======================================================================#
def create_element(type):
    """
    Constructor function to return an element of the required type
    with the appropriate number and type of points.
    """
    try:
        element = element_types[type]
    except KeyError:
        raise NotImplementedError('Element type ' + type + ' not implemented')

    try:
        element.node_coordinates = element_node_coords[type]
    except KeyError:
        raise NotImplementedError('Element node coordinate type ' + type + ' not implemented')

    try:
        element.edges = element_edges[type]
    except KeyError:
        raise NotImplementedError('Element edge type ' + type + ' not implemented')

    try:
        element.edge_points = element_edge_points[type]
    except KeyError:
        raise NotImplementedError('Element edge point type ' + type + ' not implemented')

    return element

# ======================================================================#
