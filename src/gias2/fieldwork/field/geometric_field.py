"""
FILE: geometric_field.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: class for a ensemble_field_function representing a mesh geometry.

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import copy
import json
import logging
import os
import shelve
import sys

import numpy
from scipy import sparse
from scipy.interpolate import splprep, splev, splrep
from scipy.optimize import fmin
from scipy.spatial import cKDTree

import gias2.fieldwork.field.ensemble_field_function as EFF
from gias2.common import math
from gias2.common import transform3D
from gias2.fieldwork.field import template_fields
from gias2.fieldwork.field.tools import curvature_tools as CT
from gias2.fieldwork.field.tools import discretisation
from gias2.fieldwork.field.tools import misc
from gias2.fieldwork.field.tools import triangulate
from gias2.fieldwork.field.topology import element_types

log = logging.getLogger(__name__)

try:
    from mayavi import mlab
except ImportError:
    log.debug('Mayavi not imported, 3D visualisation will be disabled')


# !!!!breaks ASM!!!!
# def normaliseVectors(x):
#     return x/(numpy.sqrt((x**2.0).sum(1))[:,numpy.newaxis])

# def normaliseVector(x):
#     return x/numpy.sqrt((x**2.0).sum())

class geometric_field:
    """ Class for an ensemble_field representing the geometry of an 
    object. Either instantiated with an existing ensemble_field, or if 
    none, a new ensemble_field will be instantiated from the 
    template_fields module with the given field_dimension and basis_type.
    """

    def __init__(self, name, dimensions, ensemble_field_function=None, field_dimensions=None, field_basis=None):
        self.name = name
        self.dimensions = dimensions
        self.field_parameters = None
        self.points = []  # list of point objects
        self.points_to_ensemble_map = {}  # { points index: ensemble_point_number }
        self.ensemble_to_points_map = {}  # { ensemble_point_number: points index }
        self.named_points_map = {}  # { point_name: points index }
        self.points_counter = 0
        self.ensemble_point_counter = 0
        self.splines = {}  # { spline name: geometric_spline object }
        # self.smoothCurvature = True

        # ~ for i in range( self.dimensions ):
        # ~ self.field_parameters.append( [] )

        if ensemble_field_function:
            # create the appropriate number of geometric points for the ensemble
            # points in the mesh
            self.ensemble_field_function = ensemble_field_function
            self._create_ensemble_points()
            try:
                self.ensemble_field_function_filename = self.ensemble_field_function.filename
            except AttributeError:
                self.ensemble_field_function_filename = self.ensemble_field_function.name + '.ens'
        else:
            # instantiate a new ensemble_field_function
            self.ensemble_field_function = template_fields.empty_field(self.name, field_dimensions, field_basis)
            # ~ self.ensemble_field_function = EFF.ensemble_field_function( self.name+'field', field_dimensions )
            # ~ self.ensemble_field_function.set_basis( field_basis    )

        # triangulator for drawing the interpolated surface
        self.triangulator = triangulate.mesh_triangulator(self.ensemble_field_function, self.field_parameters, False)

        self.ensemble_field_function_filename = None

    # ==================================================================#
    def save_geometric_field(self, filename=None, field_filename=None, mesh_filename=None, path=''):
        """
        Serialise the geometric_field instance.

        Inputs:
        filename: [str] geometric_field filename. .geof suffix will be added.
        field_filename: [str] ensemble field function filename. .ens suffix will be added.
        mesh_filename: [str] mesh filename. .mesh suffix will be added.
        
        field_filename and mesh_filename are optional. If they are not defined,
        .ens and .mesh files will not be produced. You might do this when you
        are saving a deformed version of a mesh for which you have already have
        .ens and .mesh files.
        """

        if not filename:
            filename = self.name

        if os.path.splitext(filename)[1].lower() != '.geof':
            filename = filename + '.geof'

        return save_gf_json(
            filename, self, ensfn=field_filename, meshfn=mesh_filename, filedir=path,
        )

    def save_geometric_field_shelve(self, filename, field_filename=None, mesh_filename=None, path=''):
        """
        Serialise the geometric_field instance.

        Inputs:
        filename: [str] geometric_field filename. .geof suffix will be added.
        field_filename: [str] ensemble field function filename. .ens suffix will be added.
        mesh_filename: [str] mesh filename. .mesh suffix will be added.
        
        field_filename and mesh_filename are optional. If they are not defined,
        .ens and .mesh files will not be produced. You might do this when you
        are saving a deformed version of a mesh for which you have already have
        .ens and .mesh files.
        """

        if path:
            if path[-1] != '/':
                path += '/'

        try:
            S = shelve.open(os.path.join(path, filename + '.geof'))
        except ImportError:
            import bsddb3
            _db = bsddb3.hashopen(os.path.join(path, filename + '.geof'))
            S = shelve.Shelf(_db)

        S['name'] = self.name
        if field_filename:
            field_loc = self.ensemble_field_function.save_ensemble(
                field_filename, mesh_filename=mesh_filename, path=path
            )
            S['ensemble_field'] = os.path.split(field_loc)[1]
        # ~ else:
        # ~ S['ensemble_field'] = self.ensemble_field_function.save_ensemble( self.ensemble_field_function.name, mesh_filename=mesh_filename )
        else:
            S['ensemble_field'] = None

        S['dimensions'] = self.dimensions
        S['field_parameters'] = self.field_parameters
        S['ensemble_point_counter'] = self.ensemble_point_counter

        S.close()
        return path + filename + '.geof'

    def write_nodes(self, filename, header=None):
        """Write nodal parameters to a text file. Format per line is
        [node number] [node param1] [node param2] ...
        """

        node_numbers = numpy.arange(self.field_parameters.shape[1])
        params = self.field_parameters.squeeze()
        line_pattern = '{:6d} '

        with open(filename, 'w') as f:
            f.write('# Geometric field parameters\n')
            f.write('# name: {}\n'.format(self.name))
            if header is not None:
                f.write('#' + header + '\n')
            for n in node_numbers:
                line = '{:6d} '.format(n) + ' '.join(
                    ['{:12.6f}'.format(x) for x in self.field_parameters[:, n, :].ravel()]) + '\n'
                f.write(line)

    # ==================================================================#
    def _create_ensemble_points(self):
        """ Creates a geometric point for each ensemble point in the
        ensemble field function upon instantiation with an existing 
        ensemble_field_function
        """

        # get total number of points need to make
        n_points = self.ensemble_field_function.get_number_of_ensemble_points()

        for i in range(n_points):
            # make geometric_point object
            point = geometric_point(self.dimensions, ensemble_point_number=i)
            self.points.append(point)
            self.ensemble_to_points_map[i] = self.points_counter
            self.points_to_ensemble_map[self.points_counter] = i
            self.points_counter += 1

        # update ensemble points attribute of points
        self._set_point_ensemble_numbers()

        return 1

    # ==================================================================#
    def get_number_of_points(self):
        """ returns total number of geometric_points in the
        geometric_field
        """
        return len(self.points)

    # ==================================================================#
    def set_field_parameters(self, field_parameters):
        """ Sets the field parameters for the ensemble_field. A 
        geometric point is create at each location.
        
        field_parameters should be a list with self.dimension number of 
        sub lists containing the parameters from each dimension
        """
        # ~ pdb.set_trace()
        # checks
        if len(field_parameters) != self.dimensions:
            raise ValueError('ERROR: geometric_field.set_parameters: parameters inconsistent with number of dimensions')

        p_lengths = [len(p) for p in field_parameters]
        n_ensemble_points = self.ensemble_field_function.get_number_of_ensemble_points()
        if (sum(p_lengths) / self.dimensions) != p_lengths[0]:
            raise ValueError('ERROR: geometric_field.set_field_parameters: parameter vector lengths not equal')

        if p_lengths[0] != n_ensemble_points:
            raise ValueError('ERROR: geometric_field.set_field_parameters: number of given parameters ( ' + str(
                p_lengths[0]) + ') does not match number of field ensemble points ( ' + str(n_ensemble_points) + ')')

        # if nothings wrong so far
        self.field_parameters = numpy.array(field_parameters, dtype=float)

        # set geometric points associated with ensemble point with these parameters
        # ~ for i in self.ensemble_to_points_map.keys():
        # ~ print i
        # ~ print self.ensemble_to_points_map[i]
        # ~ print self.field_parameters[:,i]
        # ~
        # ~ self.points[ self.ensemble_to_points_map[i] ].set_field_parameters( self.field_parameters[:,i] )
        # ~
        for i, ens_i in enumerate(self.ensemble_to_points_map.keys()):
            self.points[self.ensemble_to_points_map[ens_i]].set_field_parameters(self.field_parameters[:, i])

        return 1

    # ==================================================================#
    def transformAffine(self, T):
        x = self.get_field_parameters()[:, :, 0].T
        xT = transform3D.transformAffine(x, T)
        # self.set_field_parameters( xT.T[:,:,numpy.newaxis] )
        self.field_parameters = xT.T[:, :, numpy.newaxis]

    def transformTranslate(self, T):
        x = self.get_field_parameters()[:, :, 0]
        xT = x + T[:, numpy.newaxis]
        # self.set_field_parameters( xT[:,:,numpy.newaxis] )
        self.field_parameters = xT[:, :, numpy.newaxis]

    def transformRigid(self, T):
        x = self.get_field_parameters()[:, :, 0].T
        xT = transform3D.transformRigid3D(x, T)
        # self.set_field_parameters( xT.T[:,:,numpy.newaxis] )
        self.field_parameters = xT.T[:, :, numpy.newaxis]

    def transformRigidRotateAboutCoM(self, T):
        x = self.get_field_parameters()[:, :, 0].T
        xT = transform3D.transformRigid3DAboutCoM(x, T)
        # self.set_field_parameters( xT.T[:,:,numpy.newaxis] )
        self.field_parameters = xT.T[:, :, numpy.newaxis]

    def transformRigidScale(self, T):
        x = self.get_field_parameters()[:, :, 0].T
        xT = transform3D.transformRigidScale3D(x, T)
        # self.set_field_parameters( xT.T[:,:,numpy.newaxis] )
        self.field_parameters = xT.T[:, :, numpy.newaxis]

    def transformRigidScaleRotateAboutCoM(self, T):
        x = self.get_field_parameters()[:, :, 0].T
        xT = transform3D.transformRigidScale3DAboutCoM(x, T)
        # self.set_field_parameters( xT.T[:,:,numpy.newaxis] )
        self.field_parameters = xT.T[:, :, numpy.newaxis]

    def transformRigidRotateAboutP(self, T, P):
        x = self.get_field_parameters()[:, :, 0].T
        xT = transform3D.transformRigid3DAboutP(x, T, P)
        # self.set_field_parameters( xT.T[:,:,numpy.newaxis] )
        self.field_parameters = xT.T[:, :, numpy.newaxis]

    def transformRigidScaleRotateAboutP(self, T, P):
        x = self.get_field_parameters()[:, :, 0].T
        xT = transform3D.transformRigidScale3DAboutP(x, T, P)
        # self.set_field_parameters( xT.T[:,:,numpy.newaxis] )
        self.field_parameters = xT.T[:, :, numpy.newaxis]

    def transformRotateAboutP(self, r, P):
        x = self.get_field_parameters()[:, :, 0].T
        xT = transform3D.transformRotateAboutP(x, r, P)
        # self.set_field_parameters( xT.T[:,:,numpy.newaxis] )
        self.field_parameters = xT.T[:, :, numpy.newaxis]

    def transformRotateAboutAxis(self, theta, p0, p1):
        x = self.get_field_parameters()[:, :, 0].T
        xT = transform3D.transformRotateAboutAxis(x, theta, p0, p1)
        # self.set_field_parameters( xT.T[:,:,numpy.newaxis] )
        self.field_parameters = xT.T[:, :, numpy.newaxis]

    def transformRotateAboutCartCS(self, r, o, v1, v2, v3):
        x = self.get_field_parameters()[:, :, 0].T
        xT = transform3D.transformRotateAboutCartCS(x, r, o, v1, v2, v3)
        # self.set_field_parameters( xT.T[:,:,numpy.newaxis] )
        self.field_parameters = xT.T[:, :, numpy.newaxis]

    # ==================================================================#
    def get_field_parameters(self):

        if self.field_parameters is not None:
            return self.field_parameters.copy()
        else:
            log.debug("no field parameters set")
            return None

            # ==================================================================#

    def add_geometric_point(self, field_parameters, name=None):
        """ Add a point to the geometric field. 
        
        field parameters length must = field dimensions
        returns new point number
        """

        if len(field_parameters) != self.dimensions:
            raise ValueError('ERROR: geometric_field.add_point: parameters length/dimension mismatch')

        self.points.append(geometric_point(self.dimensions, numpy.array(field_parameters)))
        if name:
            self.named_points_map[name] = self.points_counter

        self.points_counter += 1
        return (self.points_counter - 1)

    # ==================================================================#
    def _set_point_ensemble_numbers(self):
        for pi, ei in list(self.points_to_ensemble_map.items()):
            self.points[pi].set_ensemble_point_number(ei)

    # ==================================================================#
    def remove_geometric_points(self, point_numbers):

        old_en2pt = dict(self.ensemble_to_points_map)
        old_pt2en = dict(self.points_to_ensemble_map)

        new_points = []
        new_points_counter = 0
        pt2en = {}
        en2pt = {}
        for pi, p in enumerate(self.points):
            if pi not in point_numbers:
                new_points.append(p)
                pt2en[new_points_counter] = p.get_ensemble_point_number()
                en2pt[p.get_ensemble_point_number()] = new_points_counter
                new_points_counter += 1

        self.points_counter = new_points_counter
        self.points = new_points
        self.ensemble_to_points_map = en2pt
        self.points_to_ensemble_map = pt2en

    # ==================================================================#
    def remove_geometric_point(self, point_number):
        """
        remove a geometric point
        """
        self.points.remove(self.points[point_number])

        # remove from named_points_map
        for k, v in list(self.named_points_map.items()):
            if v == point_number:
                del self.named_points_map[k]

        # remove from ensemble_to_points_map
        for ep, p in list(self.ensemble_to_points_map.items()):
            if p == point_number:
                del self.ensemble_to_points_map[ep]
            if p > point_number:
                self.ensemble_to_points_map[ep] -= 1

        # remove from points_to_ensemble_map
        del self.points_to_ensemble_map[point_number]

        # ~ self.points_counter -= 1

    # ==================================================================#
    def modify_geometric_point(self, point_number, new_parameters):

        try:
            self.points[point_number].set_field_parameters(new_parameters)
        except KeyError:
            log.debug('ERROR: geometric_field.modify_geometric_point: invalid point number {}'.format(point_number))
            return
        else:
            self.field_parameters[:, self.points_to_ensemble_map[point_number], :] = new_parameters

        return

    # ==================================================================#
    def add_geometric_spline(self, points, name, smoothing=3.0, order=3):
        """ add a geometric_spline object given the geometric points
        it will be fitted to.
        
        points is a list of geometric point numbers
        name will be the key to the spline from the self.splines dict
        Returns the spline instance
        """
        # get coordinates of points
        coords = []
        for p in points:
            coords.append(self.get_point_position(p))
        coords = numpy.array(coords).T

        # ~ pdb.set_trace()

        # ~ self.splines[name] = spline_3d_interp_x( coords[0], coords[1], coords[2], order, smoothing )
        self.splines[name] = spline_3d_parametric(coords, order, smoothing)

        return self.splines[name]

    # ==================================================================#
    def add_element(self, element, points, debug=0):
        """ Add an element object to the field, connecting its element
        points to the geometric point numbers given.
         
        element is an element or subfield object.
        points is a list of geometric point numbers in the order of the 
        element points of the element being added.
        """
        if debug:
            log.debug(self.points_to_ensemble_map)
            log.debug(points)

        # check for right number of points
        element_points = element.get_number_of_ensemble_points()
        if element_points != len(points):
            raise ValueError('ERROR: geometric_field.add_element: number of points ' + str(
                len(points)) + ' does not match number of element points ' + str(element_points))

        # element point counter
        ep = 0
        # add element to mesh
        e_number = self.ensemble_field_function.add_element(element)

        # connect element and update mapper
        for p in points:
            # if point is an existing ensemble point
            if p in list(self.points_to_ensemble_map.keys()):
                ensemble_point = self.points_to_ensemble_map[p]
                # get element point tuple for current connection point
                connected_eps = self.ensemble_field_function.mapper.get_ensemble_point_element_points(ensemble_point)
                if debug:
                    log.debug('existing point {}'.format(p))
                    log.debug('connecting {} to {}'.format((e_number, ep), connected_eps))

                # connect
                connected_eps.append((e_number, ep))
                self.ensemble_field_function.connect_element_points(connected_eps)
            else:

                if debug:
                    log.debug('new point {}'.format(p))

                # otherwise point maps to a new ensemble point    
                # ~ new_ensemble_point = max( self.ensemble_to_points_map.keys() ) + 1
                self.points_to_ensemble_map[p] = self.ensemble_point_counter
                self.ensemble_to_points_map[self.ensemble_point_counter] = p
                self.ensemble_point_counter += 1

                # add new field parameters to field_parameters
                fp = self.points[p].get_field_parameters()

                # ~ pdb.set_trace()
                try:
                    self.field_parameters = numpy.hstack([self.field_parameters, numpy.array(fp)[:, numpy.newaxis]])
                except ValueError:
                    self.field_parameters = numpy.array(fp)[:, numpy.newaxis].copy()
                # ~ for i in range( self.dimensions ):
                # ~ self.field_parameters[i].append( fp[i] )

            ep += 1

        # update mapper
        self.ensemble_field_function.map_parameters()

        # update ensemble points attribute of points
        self._set_point_ensemble_numbers()

        # return new element number
        ret = self.ensemble_field_function.mesh.element_counter - 1

        return ret

    # ==================================================================#
    def add_element_with_parameters(self, element, parameters, tol=1.0):
        """ Add an element to the ensemble_field with predefined parameters
        new geometric points will be created for new unique parameter
        coordinates. Non-uniques will be mapped to existing ensemble
        points.
        
        parameters should be [ [d0], [d1], [d2] ]
        [d0] = [ [p0], [p1], [p2] ] where [p0] is a list of the
        parameters at point 0
        
        if distance between a new node and old node is less than tol,
        they will be connected
        """

        # check parameters for right dimensionality and number of points
        element_points = element.get_number_of_ensemble_points()
        parameters = numpy.array(parameters)
        if parameters.shape[0:2] != (self.dimensions, element_points):
            log.debug(
                'ERROR: geometric_field.add_element_with_parameters: parameters length/dimension mismatch. need {}/{}'.format(
                    self.dimensions, element_points
                )
            )
            return None

        # get geometric point numbers for connecting the element to
        existing_positions = self.get_all_point_positions()
        if len(existing_positions) == 0:
            noCompare = 1
        else:
            noCompare = 0

        connect_points = []
        for i in range(element_points):

            if noCompare:
                connect_points.append(self.add_geometric_point(parameters[:, i]))
            else:
                # ith parameter coordinates
                coord = parameters[:, i, 0]
                # calculate distance of coord to existing points
                D = numpy.sqrt(((existing_positions - coord) ** 2.0).sum(1))
                if D.min() < tol:
                    connect_points.append(numpy.argmin(D))
                else:
                    # create new geometric point at current parameter 
                    # coordinates, and record new geometric point number
                    # for connection
                    connect_points.append(self.add_geometric_point(parameters[:, i]))

        # add element
        # ~ print 'new element points:'
        # ~ print connect_points
        ret = self.add_element(element, connect_points)

        return ret

    # ==================================================================#
    def remove_element(self, element_number):
        """ removes element defined by element number. Also removes associated
        geometric points. Geometric point numbers are reset, as are 
        ensemble point numbers.
        """

        F = self.ensemble_field_function

        old_el2ens_map = dict(F.mapper._element_to_ensemble_map)

        del_points = []
        # get non-connected points
        for ep in F.mesh.element_points[element_number]:
            if len(F.mesh.connectivity[ep]) == 0:
                # get ensemble point number of this element point
                del_points.append(F.mapper._element_to_ensemble_map[ep[0]][ep[1]][0][0])

        # remove non-connected geometric points
        self.remove_geometric_points([self.ensemble_to_points_map[p] for p in del_points])
        old_en2pt_map = dict(self.ensemble_to_points_map)

        # remove element
        F.remove_element(element_number)
        F.map_parameters()
        self.ensemble_point_counter = max(F.mapper._ensemble_to_element_map.keys()) + 1

        # make new field_parameters array and ens2pt map
        new_params = numpy.zeros([self.field_parameters.shape[0],
                                  F.get_number_of_ensemble_points(),
                                  self.field_parameters.shape[2]],
                                 dtype=float)

        en2el_map = F.mapper._ensemble_to_element_map
        self.ensemble_to_points_map = {}
        for i, ens in enumerate(en2el_map.keys()):
            elem = list(en2el_map[ens].keys())[0]
            ep = list(en2el_map[ens][elem].keys())[0]
            old_ens = old_el2ens_map[elem][ep][0][0]
            new_params[:, i, :] = self.field_parameters[:, old_ens, :]
            self.ensemble_to_points_map[ens] = old_en2pt_map[old_ens]

        # rebuild points2ens maps    
        self.points_to_ensemble_map = {p: e for e, p in list(self.ensemble_to_points_map.items())}
        # update ensemble points attribute of points
        self._set_point_ensemble_numbers()

        # ~ print 'new params shape:', new_params.shape
        self.set_field_parameters(new_params)
        self.ensemble_field_function.element_param_cache.clear()
        return

        # ==================================================================#

    # ~ def interactive_add_element( self, element, debug = 0):
    # ~ # add an existing element object to the geometric field using an
    # ~ # interactive mayavi interface
    # ~
    # ~ # get number of points needed
    # ~ element_points = element.get_number_of_ensemble_points()
    # ~
    # ~ # generate the scene in which the point picker will operate
    # ~ field_eval_density = 10
    # ~ f = mlab.figure(1)
    # ~ self._plot_points( label = True )
    # ~ self._plot_field( field_eval_density )
    # ~ mlab.show
    # ~
    # ~ # get the points to which the element will connect to using the
    # ~ # mayavi interactive points picker
    # ~

    # ==================================================================#
    def hang_point(self, point, elem, mode, xi):
        """ hang point(geometric point number) off elem (element number)
        using mode='edge' or 'interior' at the edge, or elemnt xi
        """

        pointMap = self.get_point_mesh_info(point)
        ep = []
        for e in pointMap:
            ep.append((e, list(pointMap[e].keys())[0]))

        log.debug('element points to hang: {}'.format(ep))
        self.ensemble_field_function.connect_to_hanging_point(elem, ep, mode, xi)
        self.ensemble_field_function.map_parameters()

        # need to remove hanging point's field parameters
        ensPointNumber = self.points_to_ensemble_map[point]
        fp = self.field_parameters.copy()
        newFp = numpy.hstack((fp[:, :ensPointNumber, :], fp[:, (ensPointNumber + 1):, :]))
        self.set_field_parameters(newFp)

    # ==================================================================#
    def flatten_ensemble_field_function(self):
        self.ensemble_field_function_old = copy.deepcopy(self.ensemble_field_function)
        flatEFF = self.ensemble_field_function.flatten()[0]
        self.ensemble_field_function = flatEFF

    # ==================================================================#
    def invertSurfaceMesh(self, nodeDistTol=1e-3):
        """
        remeshes a 2D mesh so that surface normals are flipped
        """

        invertMap = {
            'quad33': [2, 1, 0, 5, 4, 3, 8, 7, 6],
            'quad44': [3, 2, 1, 0, 7, 6, 5, 4, 11, 10, 9, 8, 15, 14, 13, 12],
            'quad55': [4, 3, 2, 1, 0, 9, 8, 7, 6, 5, 14, 13, 12, 11, 10, 19, 18, 17, 16, 15, 24, 23, 22, 21, 20],
            'quad54': [4, 3, 2, 1, 0, 9, 8, 7, 6, 5, 14, 13, 12, 11, 10, 19, 18, 17, 16, 15],
            'tri6': [2, 1, 0, 5, 4, 3],
            'tri10': [3, 2, 1, 0, 6, 5, 4, 8, 7, 9],
            'tri15': [4, 3, 2, 1, 0, 8, 7, 6, 5, 11, 10, 9, 13, 12, 14]
        }

        oldEFF = copy.deepcopy(self.ensemble_field_function)
        basisTypes = dict([(e, oldEFF.basis[e].type) for e in list(oldEFF.basis.keys())])
        oldParams = self.get_field_parameters()
        M = oldEFF.mapper._element_to_ensemble_map

        self.ensemble_field_function = template_fields.empty_field(self.ensemble_field_function.name, \
                                                                   self.ensemble_field_function.dimensions, \
                                                                   basisTypes)

        # refresh attributes
        self.field_parameters = None
        self.points = []  # list of point objects
        self.points_to_ensemble_map = {}  # { points index: ensemble_point_number }
        self.ensemble_to_points_map = {}  # { ensemble_point_number: points index }
        self.named_points_map = {}  # { point_name: points index }
        self.points_counter = 0
        self.ensemble_point_counter = 0
        self.triangulator = triangulate.mesh_triangulator(self.ensemble_field_function, self.field_parameters, False)

        # re-add elements with reordered parameters
        for elemNumber in list(M.keys()):
            elem = oldEFF.mesh.elements[elemNumber]
            oldNodeNumbers = numpy.array([M[elemNumber][v][0][0] for v in M[elemNumber]])[
                invertMap[elem.type]]  # look up mapping
            elemParams = oldParams[:, oldNodeNumbers, :]
            self.add_element_with_parameters(elem, elemParams, tol=nodeDistTol)

        self.ensemble_field_function.map_parameters()

    # ==================================================================#
    def get_all_point_positions(self):
        """ Returns a list containing a self.dimensions long list of
        coordinates for each geometric point
        """
        # positions = []
        # for p in self.points:
        #     positions.append( p.get_position() )

        # return numpy.array( positions )
        if self.field_parameters is not None:
            return numpy.array(self.field_parameters[:, :, 0].T)
        else:
            return numpy.array([])

    # ==================================================================#
    def get_point_position(self, point):
        """ Returns the coordinates of a point
        """
        # try:
        #     point = self.points[point]
        # except IndexError:
        #     print 'ERROR: geometric_field.get_point_position: invalid point number'
        #     return None
        # else:
        #     return point.get_position()
        if self.field_parameters is not None:
            try:
                return numpy.array(self.field_parameters[:, point, 0])
            except IndexError:
                raise ValueError('Invalid point number')
        else:
            raise RuntimeError('No field parameters set')

    # ==================================================================#
    def get_point_mesh_info(self, point):
        """ returns entry in the ensemble field function mapper
        ensemble to element map entry for the given geometric
        point number
        """
        return self.ensemble_field_function.mapper._ensemble_to_element_map[self.points_to_ensemble_map[point]]

    # ==================================================================#
    def evaluate_geometric_field(self, density, derivs=None):
        """ evaluates the field for all parameter components.
        Returns a list of self.dimension lists
        """
        V = []
        D = []

        # evaluate each coordinate in field
        if not derivs:

            return numpy.array([self.ensemble_field_function.evaluate_field_in_mesh(density, d, unpack=True) for d in
                                self.field_parameters])

            # ~ for d in self.field_parameters:
            # ~ V.append( numpy.array(self.ensemble_field_function.evaluate_field_in_mesh( density, d )).flatten() )
            # ~ return numpy.array(V)
        else:
            for d in self.field_parameters:
                v, d = self.ensemble_field_function.evaluate_field_in_mesh(density, d, derivs)
                V.append(numpy.array(v).flatten())
                D.append(d)

            return numpy.array(V), numpy.array(D)

    # ==================================================================#
    def get_element_numbers(self, coordinates=True):
        """Return the numbers of the mesh elements. If coordinates is True, a set of coordinates
        corresponding to each element number will also be returned. The coordinates are the
        centres of each elements. Useful for displaying element numbers
        """

        elem_numbers = list(self.ensemble_field_function.mesh.elements.keys())

        if coordinates:
            coords = []
            all_node_coords = self.get_all_point_positions()
            for e_i in elem_numbers:
                if self.ensemble_field_function.mesh.elements[e_i].is_element:
                    # evaluate element at 0.5, 0.5 in element coords
                    # coords.append( self.evaluate_geometric_field_at_element_points( e_i, numpy.array([0.5,0.5]) ) )
                    # evaluate centre of mass of element nodes
                    elem_nodes = [v[0][0] for v in
                                  list(self.ensemble_field_function.mapper._element_to_ensemble_map[e_i].values())]
                    coords.append(all_node_coords[elem_nodes, :].mean(0))
                else:
                    # find centre of mass of region nodes
                    region_ens_points = [v[0][0] for v in list(
                        self.ensemble_field_function.mapper._element_to_ensemble_map[e_i].values())]
                    region_nodes = numpy.array(
                        [self.get_point_position(self.ensemble_to_points_map[i]) for i in region_ens_points])
                    coords.append(region_nodes.mean(0))

            return elem_numbers, coords
        else:
            return elem_numbers

    # ==================================================================#
    def evaluate_geometric_field_in_elements(self, density, elements, derivs=None):
        """ evaluates the field for all parameter components in defined
        elements.
        Returns a list of self.dimension lists
        """
        V = []
        D = []

        # evaluate each coordinate in field
        if not derivs:
            for p in self.field_parameters:
                V.append(
                    numpy.hstack([
                        self.ensemble_field_function.evaluate_field_in_element(
                            e, density, p, unpack=True
                        ) for e in elements
                    ])
                )

            return numpy.array(V)
        else:
            for p in self.field_parameters:
                vv = []
                dd = []
                for e in elements:
                    v, d = self.ensemble_field_function.evaluate_field_in_element(
                        e, density, p, derivs
                    )
                    vv.append(v)
                    dd.append(d)

                V.append(numpy.hstack(vv))
                D.append(numpy.hstack(dd))

            return numpy.array(V), numpy.array(D)

    # ==================================================================#
    def evaluate_geometric_field_at_element_points(self, element, XI, derivs=None):
        """
        inputs
        ------
        element : (int) the number of the element to evaluate
        XI : (n*m array) a list of x XI coordinates to evaluate. m is the dimension of the element
        derivs : (tuple|-1) field derivative to evaluate. If a tuple (1,0),
                 evaluate the specified derivative. If -1, evaluate all derivatives.

        returns
        -------
        X : (p*n) n evaluated coordinates, p is the dimension of the geometric field.
            if derivs is defined, X is of shape (p,q,n) where q is the number of
            derivatives + 1. (:,0,:) is always the field value.
        """
        C = []
        # evaluate derivatives of coordinate fields
        for P in self.field_parameters:
            C.append(
                self.ensemble_field_function.evaluate_field_at_element_point(
                    element, XI, parameters=P, derivs=derivs
                )
            )

        return numpy.array(C)

    def evaluate_geometric_field_at_element_points_2(self, EP):
        """
        evaluate in multiple elements. EP is a dictionary 
        {elem_number:[xi]}. Elements are evaluated in order of element
        number.
        """

        elements = list(EP.keys())
        elements.sort()

        if self.ensemble_field_function.is_flat():
            # check for empty lists
            try:
                X = numpy.hstack([
                    self.evaluate_geometric_field_at_element_points(
                        e, numpy.array(EP[e])
                    ) for e in elements if len(numpy.array(EP[e])) > 0
                ])
            except TypeError:
                raise TypeError('bad xi type')
        else:
            # check for empty lists
            # ~ X = numpy.hstack([self.evaluate_geometric_field_at_element_points( e, EP[e] ) for e in elements])
            X = numpy.hstack([
                self.evaluate_geometric_field_at_element_points(
                    e, EP[e]
                ) for e in elements if len(EP[e]) > 0
            ])

        return X

    def evaluate_geometric_field_at_element_points_3(self, EP):
        """
        evaluate in multiple elements. EP is list of element points expressed
        as a tuple (element_number, xi_coordinates). Example of input EP:
        [(0, [0.5, 0.1]), (1, [0.1, 0.1]),]
        """
        evaluator = makeGeometricFieldEvaluatorSparse(self, [1, 1], matPoints=EP)
        X = evaluator(self.field_parameters)
        return X

    # ==================================================================#
    def evaluate_curvature_in_mesh(self, density, smooth=False):
        """ evaluates the guassian K and mean H curvature over the mesh
        at the desired density
        """
        D = []
        V = []

        # evaluate derivatives of coordinate fields
        for P in self.field_parameters:
            v, d = self.ensemble_field_function.evaluate_field_in_mesh(
                density, parameters=P, derivs=-1
            )
            D.append(d)
            V.append(v)
            # ~ D.append( self.ensemble_field_function.evaluate_field_in_mesh( density, parameters=P, derivs=-1 )[1] )

        K, H, k1, k2 = self._calculate_curvature(D)
        V = numpy.array(V).T
        # smooth curvature field
        if smooth:
            H = smoothCurvField(V, H)
            K = smoothCurvField(V, K)
            k1 = smoothCurvField(V, k1)
            k2 = smoothCurvField(V, k2)

        return (K, H, k1, k2)

    # ==================================================================#
    def evaluate_curvature_at_element_points(self, element, XI):
        """ evaluates the guassian K and mean H curvature in an element
        at a list of xi positions
        """

        D = []
        # evaluate derivatives of coordinate fields
        for P in self.field_parameters:
            D.append(
                self.ensemble_field_function.evaluate_field_at_element_point(
                    element, XI, parameters=P, derivs=-1
                )[1]
            )

        return self._calculate_curvature(D)

    # ==================================================================#
    def evaluate_normal_in_mesh(self, d, elemXi=None):
        """ evaluates the normal vector at points on elements for all
        elements in the mesh if d is tuple for discretisation.
        
        If elemXi is not none, it should be a dictionary of element
        number: xi coordinates that define where the normals should be
        calculated
        """

        if elemXi is None:
            d10 = []
            d01 = []
            for p in self.field_parameters:
                d10.append(
                    self.ensemble_field_function.evaluate_derivatives_in_mesh(
                        d, parameters=p, derivs=(1, 0), unpack=True
                    )
                )
                d01.append(
                    self.ensemble_field_function.evaluate_derivatives_in_mesh(
                        d, parameters=p, derivs=(0, 1), unpack=True
                    )
                )
        else:
            d10 = []
            d01 = []
            for p in self.field_parameters:
                self.ensemble_field_function.set_parameters(p)
                d10i = []
                d01i = []
                for elemNumber in list(elemXi.keys()):
                    d10i.append(
                        self.ensemble_field_function.evaluate_field_at_element_point(
                            elemNumber, elemXi[elemNumber], derivs=(1, 0)
                        )[1]
                    )
                    d01i.append(
                        self.ensemble_field_function.evaluate_field_at_element_point(
                            elemNumber, elemXi[elemNumber], derivs=(0, 1)
                        )[1]
                    )

                d10.append(numpy.hstack(d10i))
                d01.append(numpy.hstack(d01i))

        # ~ pdb.set_trace()
        # d10Norm = normaliseVectors( numpy.array(d10).T )
        # d01Norm = normaliseVectors( numpy.array(d01).T )
        # return numpy.cross( d10Norm, d01Norm ).T

        # !!! normaliseVectors BREAKS ASM !!!
        # CTM does its own normalistion, and works regardless of this.
        # ASM training depends on this for normalisation. However, 
        # it doesn't seem to matter whether this is on or not for training
        # as long as its off for segmentation

        N = math.norms(
            numpy.cross(
                numpy.array(d10).T, numpy.array(d01).T
            )
        ).T
        return N

    def evaluate_normal_at_element_point(self, elem, xi):
        d10 = numpy.zeros(3, dtype=float)
        d01 = numpy.zeros(3, dtype=float)
        for i, p in enumerate(self.field_parameters):
            self.ensemble_field_function.set_parameters(p)
            d10[i] = self.ensemble_field_function.evaluate_field_at_element_point(
                elem, xi, derivs=(1, 0)
            )[1]
            d01[i] = self.ensemble_field_function.evaluate_field_at_element_point(
                elem, xi, derivs=(0, 1)
            )[1]

        # d10Norm = normaliseVector(d10)
        # d01Norm = normaliseVector(d01)
        # return numpy.cross( d10Norm, d01Norm )
        return math.norm(numpy.cross(d10, d01))

    # ==================================================================#
    def _calculate_curvature(self, D):

        H = None
        K = None
        k1 = None
        k2 = None
        tol = 1e-12
        # concatenate vectors
        D = numpy.array(D).T
        D = numpy.where(abs(D) < tol, 0.0, D)
        # ~ print 'D ',D.shape
        # ~ A = numpy.zeros( (6,len(D)) )

        # ~ L = numpy.zeros( len(D) )
        # ~ M = numpy.zeros( len(D) )
        # ~ N = numpy.zeros( len(D) )
        # ~ E = numpy.zeros( len(D) )
        # ~ F = numpy.zeros( len(D) )
        # ~ G = numpy.zeros( len(D) )

        i = 0

        n = numpy.cross(D[:, 0, :], D[:, 1, :])
        n = n / numpy.sqrt((n ** 2.0).sum(1))[:, numpy.newaxis]

        E = D[:, 0, 0] * D[:, 0, 0] + D[:, 0, 1] * D[:, 0, 1] + D[:, 0, 2] * D[:, 0, 2]
        F = D[:, 0, 0] * D[:, 1, 0] + D[:, 0, 1] * D[:, 1, 1] + D[:, 0, 2] * D[:, 1, 2]
        G = D[:, 1, 0] * D[:, 1, 0] + D[:, 1, 1] * D[:, 1, 1] + D[:, 1, 2] * D[:, 1, 2]

        L = D[:, 2, 0] * n[:, 0] + D[:, 2, 1] * n[:, 1] + D[:, 2, 2] * n[:, 2]
        M = D[:, 4, 0] * n[:, 0] + D[:, 4, 1] * n[:, 1] + D[:, 4, 2] * n[:, 2]
        N = D[:, 3, 0] * n[:, 0] + D[:, 3, 1] * n[:, 1] + D[:, 3, 2] * n[:, 2]

        # do curvature calculations    
        # ~ for i, d in enumerate(D):
        # ~
        # ~ L[i] = numpy.dot( d[2], n[i] )
        # ~ M[i] = numpy.dot( d[4], n[i] )
        # ~ N[i] = numpy.dot( d[3], n[i] )
        # ~
        # ~ E[i] = numpy.dot( d[0], d[0] )
        # ~ F[i] = numpy.dot( d[0], d[1] )
        # ~ G[i] = numpy.dot( d[1], d[1] )
        # ~
        # ~ A = numpy.where(abs(A)<tol, 0.0, A)

        E = numpy.where(abs(E) < tol, 0.0, E)
        F = numpy.where(abs(F) < tol, 0.0, F)
        G = numpy.where(abs(G) < tol, 0.0, G)
        L = numpy.where(abs(L) < tol, 0.0, L)
        M = numpy.where(abs(M) < tol, 0.0, M)
        N = numpy.where(abs(N) < tol, 0.0, N)

        # ~ [L,M,N,E,F,G] = A

        K = (L * N - M * M) / (E * G - F * F)
        # ~ H = -(L*G - 2.0*M*F + N*E) / (2.0*(E*G - F**2.0)**1.5)
        H = -(L * G - 2.0 * M * F + N * E) / (2.0 * (E * G - F * F))

        K = numpy.where(abs(K) < tol, 0.0, K)
        H = numpy.where(abs(H) < tol, 0.0, H)

        # calculate principal curvatures
        k1 = H + numpy.sqrt(H * H - K)
        k2 = H - numpy.sqrt(H * H - K)

        return K, H, k1, k2

    # ==================================================================#
    # def findElementsContaining(self, points):
    #     """
    #     Given a list of points, find the elements in which they belong. Works by
    #     discretising the GF then using a kdtree search to find the n 
    #     closest material points to each input point. Input point is in the element
    #     where most of the n closest material points come from
    #     """
    #     xid = [10,10]
    #     k = 8

    #     # discretise gf
    #     ep = self.evaluate_geometric_field(xid, unpack=False)

    #     # create array of element numbers corresponding to material points
    #     elemArray = []
    #     for ei, e in enumerate(self.ensemble_field_function.mesh.elements.keys()):
    #         elemArray.append([e]*ep[ei].shape[1])

    #     elemArray = numpy.hstack(elemArray)

    #     # create material point kdtree
    #     ep = numpy.hstack(ep).T
    #     epTree = cKDTree(ep)

    #     # search tree
    #     d, i = epTree.query(points, k=k)

    def _makeXiObj(self, elemNum):
        # initialise
        E = self.ensemble_field_function.mesh.elements[elemNum]
        evaluator = self.ensemble_field_function.evaluators[E.type]
        basis = self.ensemble_field_function.basis[E.type]
        P = []
        for p in self.field_parameters:
            self.ensemble_field_function.set_parameters(p)
            P.append(
                self.ensemble_field_function._get_element_parameters(
                    elemNum
                )
            )

        # objective function
        def findXiObj(xi, target):
            B = basis.eval(xi)
            v = numpy.dot(P, B)
            d = ((v - target) * (v - target)).sum()
            return d

        findXiObj.P = P
        findXiObj.basis = basis

        return findXiObj

    def findXi(self, elem, target, initXi=None, fullOutput=False, findXiObj=None):
        """ returns the xi coordinates of a target point inside the
        specified element
        """

        if initXi is None:
            initXi = numpy.ones(self.ensemble_field_function.dimensions, dtype=float) * 0.5

        if findXiObj is None:
            findXiObj = self._makeXiObj(elem)

        # # initialise
        # E = self.ensemble_field_function.mesh.elements[elem]
        # evaluator = self.ensemble_field_function.evaluators[E.type]
        # basis = self.ensemble_field_function.basis[E.type]
        # P = []
        # for p in self.field_parameters:
        #     self.ensemble_field_function.set_parameters( p )
        #     P.append( self.ensemble_field_function._get_element_parameters( elem ) )

        # # objective function
        # def findXiObj( xi ):
        #     B = basis.eval( xi )
        #     v = numpy.dot(P, B)
        #     #~ v = [ evaluator(B,p) for p in P ]
        #     d = ((v - target)**2.0 ).sum()
        #     return d    

        # initXi = numpy.array( initXi )
        # target = numpy.array( target )
        xiOpt = fmin(findXiObj, initXi, args=(target,), disp=False)
        d = numpy.sqrt(findXiObj(xiOpt, target))
        coord = numpy.dot(findXiObj.P, findXiObj.basis.eval(xiOpt))
        return xiOpt, coord, d

    # ==================================================================#
    def find_closest_material_points(self, dataPoints, initGD=None, verbose=False):
        """
        returns 
        closestMPs = [ (elem, xi),... ]
        closestPoints = [ (x,y,z), (x,y,z),...]
        distances = [ d1, d2,....]
        for each data point
        """

        if verbose:
            log.debug('searching for closest material point for %i points' % (len(dataPoints)))

        if initGD is None:
            initGD = numpy.ones(self.ensemble_field_function.dimensions, dtype=int) * 40

        # initial scattering of EPs
        elemNumbers = list(self.ensemble_field_function.mesh.elements.keys())
        elements = [self.ensemble_field_function.mesh.elements[k] for k in elemNumbers]

        # ~ initEP = self.evaluate_geometric_field_in_elements(initGD, elemNumbers)
        # ~ initEP = self.evaluate_geometric_field(initGD, unpack=False)

        # ~ initEP = [e.T for e in initEP]
        initXi = [e.generate_eval_grid(initGD) for e in elements]

        initEP = []
        elemXiObjs = {}
        for ei, e in enumerate(elemNumbers):
            initEP.append(self.evaluate_geometric_field_at_element_points(e, initXi[ei]).T)
            elemXiObjs[e] = self._makeXiObj(e)

        # if self.ensemble_field_function.is_flat():
        #     initXi, initEP = misc._removeDuplicatesFlat( initXi, initEP )
        # else:
        #     initXi, initEP = misc._removeDuplicates( initXi, initEP )

        initEP = numpy.vstack(initEP)

        # map of initEP to elements
        epElems = numpy.zeros(initEP.shape[0], dtype=int)
        counter = 0
        for elemI, elemXi in enumerate(initXi):
            epElems[counter:counter + len(elemXi)] = elemNumbers[elemI]
            counter += len(elemXi)

        initXi = numpy.vstack(initXi)

        # for each datapoint, find its closest ep, and ep element
        initEPTree = cKDTree(initEP)
        initClosestInd = initEPTree.query(list(dataPoints))[1]
        closestMPs = []
        for epInd in initClosestInd:
            closestMPs.append([epElems[epInd], initXi[epInd]])

        # for each datapoint, do findXi in its closest element
        closestPoints = []
        distances = []
        for pi, p in enumerate(dataPoints):
            closestXi, coord, d = self.findXi(
                closestMPs[pi][0], p,
                initXi=closestMPs[pi][1],
                fullOutput=True,
                findXiObj=elemXiObjs[closestMPs[pi][0]],
            )
            closestXi = closestXi.clip(0.0, 1.0)
            # closestXi = numpy.where(closestXi<0.0, 0.0, closestXi)
            # closestXi = numpy.where(closestXi>1.0, 1.0, closestXi)
            closestMPs[pi][1] = closestXi
            closestPoints.append(coord)
            distances.append(d)

            if verbose:
                sys.stdout.write('closest distance for point %i: %5.3f\r' % (pi, d))
                sys.stdout.flush()

        closestPoints = numpy.array(closestPoints)
        distances = numpy.array(distances)

        return closestMPs, closestPoints, distances

    # ==================================================================#
    def calc_CoM(self):
        x = self.get_all_point_positions()
        return x.mean(0)

    def calc_CoM_2D(self, d, elem=None):
        """ calc the centre of mass of the field using a triangulation 
        of the 2D surface and the area of triangle
        """

        # calc areas of each triangle
        aT, a, T, ep = self.calc_surface_area(d, element=elem, fullOutput=True)

        # calc barycenter of each triangle
        b = ep[T].sum(1) / 3.0

        # ~ pdb.set_trace()
        # centre of mass
        CoM = (a[:, numpy.newaxis] * b).sum(0) / aT
        return CoM

    # ==================================================================#
    def calc_surface_area(self, d, element=None, fullOutput=False):
        """ calculates the approximate area of the mesh, or element if
        defined, as the sum of the area of triangles descritising the
        mesh at a factor of d
        """

        if element is not None:
            ep = numpy.array([self.ensemble_field_function.evaluate_field_in_element(element, d, p) for p in
                              self.field_parameters]).T
            elem = self.ensemble_field_function.mesh.elements[element]
            if elem.is_element:
                T = triangulate.triangulate([elem], d)
            else:
                T = triangulate.triangulate(elem.get_true_elements(), d)
        else:
            ep = numpy.array(
                [self.ensemble_field_function.evaluate_field_in_mesh(d, p) for p in self.field_parameters]).T
            T = triangulate.triangulate(self.ensemble_field_function.mesh.get_true_elements(), d)

        # calculate the side lengths of each triangle
        l1 = numpy.sqrt(((ep[T[:, 0]] - ep[T[:, 1]]) ** 2.0).sum(1))
        l2 = numpy.sqrt(((ep[T[:, 1]] - ep[T[:, 2]]) ** 2.0).sum(1))
        l3 = numpy.sqrt(((ep[T[:, 2]] - ep[T[:, 0]]) ** 2.0).sum(1))

        # heron's formula for calculating area of triangle give lengths of 3 sides
        # warning: unstable for very narrow triangles
        s = (l1 + l2 + l3) * 0.5
        a = numpy.sqrt(s * (s - l1) * (s - l2) * (s - l3))

        if fullOutput:
            return a.sum(), a, T, ep
        else:
            return a.sum()

    # ==================================================================#
    def display_geometric_field(self, field_eval_density, point_glyph='sphere', point_label=None, point_scale=1.0,
                                field_glyph='point', field_scale=1.0, data=None, curvature=None, figure=None,
                                scalar=None, name=None):
        # ~ def display_geometric_field( self, field_eval_density, point_glyph='sphere', point_label=None, point_scale=1.0, field_glyph='point', field_scale=1.0, data=None, curvature=None, **mLabArgs ):
        """ Plots all geometric points and point evaluated over the field
        in a mayavi scene.
        """
        spline_density = 200
        if not figure:
            f = mlab.figure()
        else:
            f = figure

        if not name:
            name = self.name

        # draw points
        if point_glyph:
            self._plot_points(point_glyph, label=point_label, scale=point_scale, figure=f)
        # ~ self._plot_field( field_eval_density, field_glyph, field_scale )

        # draw field
        if self.ensemble_field_function.dimensions == 1:
            # draw curve
            if not self.ensemble_field_function.mesh.get_number_of_elements():
                log.debug('No field found')
            else:
                self._draw_curve(field_eval_density, name=name, figure=f, tube_radius=0.5)
        elif self.ensemble_field_function.dimensions == 2:
            # draw surface
            if not self.ensemble_field_function.mesh.get_number_of_elements():
                log.debug('No field found')
            else:
                if not curvature:
                    self._draw_surface(field_eval_density, name=name, figure=f, scalar=scalar)
                else:
                    self._draw_surface_curvature(field_eval_density, curvature, name=name, figure=f)

        # draw splines
        for s in list(self.splines.values()):
            self._draw_spline(s, spline_density)

        # draw data cloud (if provided )
        if data is not None:
            if data.shape[0] == 3:
                mlab.points3d(data[0], data[1], data[2], mode='point', color=(0, 1, 0), figure=f)
            elif data.shape[0] == 4:
                mlab.points3d(data[0], data[1], data[2], data[3], mode='point', colormap='jet', vmax=0.2, vmin=-0.2,
                              figure=f)

        mlab.show()
        return f

    # ==================================================================#
    def _draw_curve_old(self, density, scene=None, **kwargs):
        E = self.evaluate_geometric_field(density)

        # remove every density points to avoid overlap
        allPoints = numpy.array(E.T)
        points = []
        for i, p in enumerate(allPoints):
            if numpy.mod(i, density[0]) != 0 or i == 0:
                points.append(p)

        # points = E.T

        points = numpy.transpose(points)
        if scene is None:
            line = mlab.plot3d(points[0], points[1], points[2], **kwargs)
        else:
            line = scene.mlab.plot3d(points[0], points[1], points[2], **kwargs)
        # ~ line = mlab.points3d( e[0], e[1], e[2], scale_factor=0.4, **kwargs )
        return line

    def _draw_curve(self, density, scene=None, **kwargs):
        n_elems = self.ensemble_field_function.mesh.get_number_of_true_elements()
        x = self.evaluate_geometric_field(density).T.reshape((n_elems, -1, 3))

        lines = []
        if scene is None:
            for epts in x:
                lines.append(
                    mlab.plot3d(epts[:, 0], epts[:, 1], epts[:, 2], **kwargs)
                )
        else:
            for epts in x:
                lines.append(
                    scene.mlab.plot3d(epts[:, 0], epts[:, 1], epts[:, 2], **kwargs)
                )
        return lines

    # ==================================================================#
    def _draw_surface(self, density, scalar=None, figure=None, name=None, lim=[None, None]):
        """ using triangulate to draw the interpolated surface
        """
        if self.field_parameters is None:
            log.debug('no field parameters')
            return
        else:
            self.triangulator.params = self.field_parameters
            return self.triangulator.draw_surface_simple(density, scalar=scalar, figure=figure, name=name, limits=lim)

    # ==================================================================#
    def _draw_surface_curvature(self, field_eval_density, curvature, name, figure=None):

        if self.field_parameters is None:
            log.debug('no field parameters')
            return
        else:
            K, H, k1, k2 = self.evaluate_curvature_in_mesh(field_eval_density)
            if curvature == 'mean':
                return self._draw_surface(field_eval_density, scalar=H, name=name, figure=figure)
            elif curvature == 'gaussian':
                return self._draw_surface(field_eval_density, scalar=K, name=name, figure=figure)

    # ==================================================================#
    def _draw_surface_curvature_binned(self, field_eval_density, curvature, name, bins, cMin, cMax, fig=None):

        if self.field_parameters is None:
            log.debug('no field parameters')
            return
        else:
            K, H, k1, k2 = self.evaluate_curvature_in_mesh(field_eval_density)
            if curvature == 'mean':
                H = CT.normalise(CT.filterCurv(H, cMin, cMax))
                binInd = numpy.digitize(H, bins)
                return self._draw_surface(field_eval_density, scalar=binInd, name=name, figure=fig,
                                          lim=(1.0, len(bins) - 1))
            elif curvature == 'gaussian':
                K = CT.normalise(CT.filterCurv(K, cMin, cMax))
                binInd = numpy.digitize(K, bins)
                return self._draw_surface(field_eval_density, scalar=binInd, name=name, figure=fig,
                                          lim=(1.0, len(bins) - 1))

    # ==================================================================#
    def _plot_points(self, glyph='sphere', label=None, scale=0.5, figure=None):
        """ uses mayavi points3d to show the positions of all points 
        (with labels if label is true)
        
        label can be 'all', or 'landmarks'
        """

        # get point positions
        p = numpy.array(self.get_all_point_positions())

        if len(p) > 0:
            # ~ f = mlab.figure()
            s = numpy.arange(len(self.points))
            if figure:
                points_plot = mlab.points3d(p[:, 0], p[:, 1], p[:, 2], s, mode='sphere', scale_mode='none',
                                            scale_factor=scale, color=(1, 0, 0), figure=figure)
            else:
                points_plot = mlab.points3d(p[:, 0], p[:, 1], p[:, 2], s, mode='sphere', scale_mode='none',
                                            scale_factor=scale, color=(1, 0, 0))

            # label all ensemble points with their index number
            if label == 'all':
                labels = list(range(len(self.points)))

                labelSceneObjs = [
                    mlab.text3d(
                        p[i, 0], p[i, 1], p[i, 2], str(labels[i]),
                        scale=5.0, color=(1, 1, 1), figure=figure,
                    ) for i in range(len(labels))
                ]

                # for i in range( len(labels ) ):
                #     l = mlab.text( p[i,0], p[i,1], str(labels[i]), z = p[i,2], line_width = 0.01, width = 0.005*len(str(labels[i]))**1.1, figure=figure )

            elif label == 'landmarks':
                m = self.named_points_map
                labels = list(m.keys())

                labelSceneObjs = [
                    mlab.text3d(
                        p[m[label]][0], p[m[label]][1], p[m[label]][2], label,
                        scale=5.0, color=(1, 1, 1), figure=figure,
                    ) for label in labels
                ]

                # for label in labels:
                #     l = mlab.text( p[m[label]][0], p[m[label]][1], label, z = p[m[label]][2], line_width = 0.01, width = 0.005*len( label )**1.1, figure=figure )

            return points_plot
        else:
            return None

    # ==================================================================#
    def _plot_field(self, density, glyph='point', scale_factor=1.0, figure=None):
        """ uses mayavi points3d to show the evaluated field
        """

        if not self.ensemble_field_function.mesh.get_number_of_elements():
            log.debug('No field found')
            return
        else:
            evaluation = []

            # evaluate each coordinate in field
            for d in self.field_parameters:
                evaluation.append(
                    numpy.array(self.ensemble_field_function.evaluate_field_in_mesh(density, d)).flatten())

            if self.dimensions == 3:
                field_plot = mlab.points3d(evaluation[0], evaluation[1], evaluation[2], mode=glyph, scale_factor=0.8,
                                           color=(0.0, 0.5, 1.0), figure=figure)

            return field_plot

    # ==================================================================#
    def _draw_spline(self, spline, density):
        """ evaluates spline and draw using plot3d
        """
        if spline.type == 'parametric':
            e = spline.eval(numpy.linspace(0, 1, density))
            line = mlab.plot3d(e[0], e[1], e[2], tube_radius=0.2)
        elif spline.type == 'interp_x':
            x = numpy.linspace(spline.xb, spline.xe, density)
            e = spline.eval(x)
            line = mlab.plot3d(x, e[0], e[1], tube_radius=0.2)
            return line

    # ==================================================================#
    def makeLineElementsFromPoints(self, points, nodesPerElem, elemBasisMap):
        """ returns a geometric field of line element(s) through the 
        global points points. elemBasisMap = {elemtype:basistype}
        """

        lineGF = geometric_field('line', 3, field_dimensions=1, field_basis=elemBasisMap)
        params = self.field_parameters[:, points, :].copy()
        nElems = (len(points) - 1) / (nodesPerElem - 1)

        # assume only one type of elements and basis
        elemTypes = [list(elemBasisMap.keys())[0]] * nElems
        x = 0
        for nE in range(nElems):
            e = element_types.create_element(elemTypes[nE])
            lineGF.add_element_with_parameters(e, params[:, x:x + nodesPerElem, :], tol=1e-3)
            x += nodesPerElem - 1

        return lineGF

    # ==================================================================#
    def makeLineElementsFromPointSets(self, pointSets, elemTypes, elemBasisMap):
        """ returns a geometric field of line element(s) through the sets of
        global points.

        inputs
        ------
        pointSets : a list of lists of node numbers. Each nested list contains
            the nodes of a line element.
        elemTypes: a dictionary of the number of nodes per line segment mapping
            to the element type name.
        elemBasisMap : a dictionary of the line element types mapping to line 
            basis functions, i.e. {elemtype:basistype}
        """

        lineGF = geometric_field('line', 3, field_dimensions=1, field_basis=elemBasisMap)
        params = self.field_parameters.copy()

        # assume only one type of elements and basis   
        for elemPoints in pointSets:
            e = element_types.create_element(elemTypes[len(elemPoints)])
            lineGF.add_element_with_parameters(
                e, params[:, elemPoints, :], tol=1e-3
            )

        return lineGF

    # ==================================================================#
    def makeElementBoundaryCurve(self, elemNumber, nNodesElemTypeMap, elemBasisMap):

        bNodes = self.ensemble_field_function.get_element_boundary_nodes(elemNumber)
        lineGF = geometric_field('line', 3, field_dimensions=1, field_basis=elemBasisMap)

        for edgeNodes in bNodes:
            e = element_types.create_element(nNodesElemTypeMap[len(edgeNodes)])
            params = self.field_parameters[:, edgeNodes, :].copy()
            lineGF.add_element_with_parameters(e, params, tol=1e-3)

        return lineGF

    # ==================================================================#
    # API functions for Duane's Fitter module                           #
    # ==================================================================#
    def setFitterElementD(self, epD):

        f = self.ensemble_field_function
        self.basisWeights = {}  # {elemNumber: {[xi]:[weights]}}
        self.elementXis = {}  # {elemNumber: [xis]}
        for ei in numpy.sort(list(f.mesh.elements.keys())):
            element = f.mesh.elements[ei]
            elementXis = element.generate_eval_grid(epD).squeeze()
            basisValues = f.basis[element.type].eval(elementXis.T).T

            # ~ pdb.set_trace()

            self.elementXis[ei] = elementXis

            basisValuesMap = {}
            for i, Xi in enumerate(list(elementXis)):
                basisValuesMap[tuple(Xi)] = basisValues[i]

            self.basisWeights[ei] = basisValuesMap

    def addElementPointsToFitter(self, fitter):

        for ei in list(self.ensemble_field_function.mesh.elements.keys()):
            fitter.add_element_points(ei, self._getFitterElementXis(ei))

    def getElementNodeNumbers(self, element):
        elemMap = self.ensemble_field_function.mapper._element_to_ensemble_map[element]
        elemNodeNumbers = [elemMap[eN][0][0] for eN in list(elemMap.keys())]
        return elemNodeNumbers

    def getElementWeights(self, element, xi):
        return self._getElementWeights(element, xi)

    def evaluateElementNumber(self, GD, unpack=False):

        if unpack:
            E = self.ensemble_field_function.flatten()[0]

        x = E.evaluate_field_in_mesh(GD, parameters=self.field_parameters[0], unpack=0)
        cc = []
        for i, e in enumerate(x):
            cc.append(numpy.ones(e.shape) * i)
        cc = numpy.hstack(cc)

        return cc

    def getNDoF(self):
        return self.ensemble_field_function.get_number_of_ensemble_points() * self.dimensions

    def getFitterParameters(self):
        P = self.get_field_parameters()
        return P.squeeze().T.copy()

    def _getElementWeights(self, elementNumber, xi):
        return self.basisWeights[elementNumber][tuple(xi)]

    def _getFitterElementXis(self, elementNumber):
        return self.elementXis[elementNumber]

    # ==================================================================#
    # functions for getting element point indices                       #
    # ==================================================================#
    def getElementPointI(self, GD, elems):
        """
        get the ep indices grouped top-level elements defined in elems
        """
        if isinstance(GD, float):
            ep = self.discretiseAllElementsRegularGeoD(GD, unpack=False)
        else:
            self.ensemble_field_function.set_parameters(self.field_parameters[0])
            ep = self.ensemble_field_function.evaluate_field_in_mesh(GD, unpack=False)

        if elems == 'all':
            elems = sorted(self.ensemble_field_function.mesh.elements.keys())

        epI = []
        for e in elems:
            epCount = 0
            for i in range(e):
                epCount += ep[i].shape[0]

            epI.append(numpy.arange(epCount, epCount + ep[e].shape[0]))

        return epI

    def getElementPointIPerTrueElement(self, GD, subMeshes):
        """ 
        gets the ep indices grouped by true elements within the top level 
        elements defined in subMeshes if mesh is not flat. If flat, 
        returns the ep indices of eps in elements defined in subMeshes
        """

        # evaluate whole mesh
        self.ensemble_field_function.set_parameters(self.field_parameters[0])
        if subMeshes is None:
            subMeshes = list(self.ensemble_field_function.mesh.elements.keys())
            subMeshes.sort()

        if self.ensemble_field_function.is_flat():
            if isinstance(GD, float):
                epPacked = self.discretiseAllElementsRegularGeoD(GD, unpack=False)
            else:
                epPacked = self.ensemble_field_function.evaluate_field_in_mesh(GD, unpack=False)
            epI = []
            I = 0

            # for each element
            for i, e in enumerate(epPacked):
                if i in subMeshes:
                    epI.append(numpy.arange(I, I + e.shape[0]))
                I += len(e)
        else:
            if isinstance(GD, float):
                epSubUnpacked = self.discretiseAllElementsRegularGeoD(GD, unpack=False)
                epSubPacked = None
            else:
                epSubUnpacked = self.ensemble_field_function.evaluate_field_in_mesh(GD, unpack=False)
                epSubPacked = self.ensemble_field_function.evaluate_field_in_mesh(GD, unpack=False, subUnpack=False)
            epI = []

            # for each subMesh of interest
            for s in subMeshes:
                subMeshI = 0
                # find starting index of the current subMesh
                for i in range(s):
                    subMeshI += epSubUnpacked[i].shape[0]

                # for each element in the current subMesh
                for e in epSubPacked[s]:
                    epI.append(numpy.arange(subMeshI, subMeshI + e.shape[0]))
                    subMeshI += e.shape[0]

        return epI

    def getElementPointINested(self, GD, subMeshes):
        """
        get the ep indices grouped by true elements then by top level
        elements defined in subMeshes
        """

        # evaluate whole mesh
        if isinstance(GD, float):
            epSubUnpacked = self.discretiseAllElementsRegularGeoD(GD, unpack=False)
            epSubPacked = None
        else:
            self.ensemble_field_function.set_parameters(self.field_parameters[0])
            epSubUnpacked = self.ensemble_field_function.evaluate_field_in_mesh(GD, unpack=False)
            epSubPacked = self.ensemble_field_function.evaluate_field_in_mesh(GD, unpack=False, subUnpack=False)

        epI = []

        # for each subMesh of interest
        for si, s in enumerate(subMeshes):
            epI.append([])
            subMeshI = 0
            # find starting index of the current subMesh
            for i in range(s):
                subMeshI += epSubUnpacked[i].shape[0]

            # for each element in the current subMesh
            for e in epSubPacked[s]:
                epI[si].append(numpy.arange(subMeshI, subMeshI + e.shape[0]))
                subMeshI += e.shape[0]

        return epI

    # ==================================================================#
    # sub mesh creation methods                                           #
    # ==================================================================#
    def makeSubfieldGF(self, elemI, name=None):
        """
        creates a new GF using subfield elemI and its corresponding
        field params. Does not handle hanging nodes.
        """
        if self.ensemble_field_function.is_flat():
            raise ValueError('not subfields exist')

        try:
            subfield = copy.deepcopy(self.ensemble_field_function.subfields[elemI])
        except KeyError:
            raise ValueError('subfield ' + str(elemI) + ' not found')
        else:
            if name is None:
                name = self.name + '_subfield_' + str(elemI)
            newGF = geometric_field(name, self.dimensions, ensemble_field_function=subfield)
            newParamsI = []
            paramMap = self.ensemble_field_function.mapper._element_to_ensemble_map[elemI]
            for k in list(paramMap.keys()):
                newParamsI.append(paramMap[k][0][0])

            newParams = self.field_parameters[:, newParamsI, :]
            newGF.set_field_parameters(newParams)

            return newGF

    def makeGFFromElements(self, name, elements, basisTypes):
        """
        make a new GF from specified elements
        """
        newGF = geometric_field(name, self.dimensions, field_dimensions=self.ensemble_field_function.dimensions,
                                field_basis=basisTypes)
        for e in elements:
            if not self.ensemble_field_function.is_flat():
                elem = self.ensemble_field_function.subfields[e]
            else:
                elem = self.ensemble_field_function.mesh.elements[e]

            elemParams = numpy.array(
                [self.ensemble_field_function.mapper.get_element_parameters(e, p, do_hack=0) for p in
                 self.get_field_parameters()])[:, :, numpy.newaxis]

            newGF.add_element_with_parameters(elem, elemParams, tol=1e-3)


        return newGF

    # ==================================================================#
    # geometric discretisation methods                                   #        
    # ==================================================================#
    def makeElementEvaluator(self, elementNumber):

        class elementGeometryEvaluator(object):

            def __init__(self, element, basis, evaluator, parameters):
                self.element = element
                self.basis = basis
                self.parameters = parameters
                self.evaluator = evaluator

            def eval(self, xi):
                basisWeights = self.basis.eval(numpy.transpose(xi))
                coords = numpy.array([self.evaluator(basisWeights, p) for p in self.parameters]).T
                return coords

        element = self.ensemble_field_function.mesh.elements[elementNumber]
        elemBasis = self.ensemble_field_function.basis[element.type]
        elemEvaluator = self.ensemble_field_function.evaluators[element.type]
        elemParameters = []
        for p in self.field_parameters:
            self.ensemble_field_function.set_parameters(p)
            elemParameters.append(self.ensemble_field_function._get_element_parameters(elementNumber))
        elemParameters = numpy.array(elemParameters)

        return elementGeometryEvaluator(element, elemBasis, elemEvaluator, elemParameters)

    def discretiseElementRegularGeoD(self, elementNumber, maxDistance, geoCoords=False):
        """
        calculates the xi coordinates that discretises an element by
        regular geometric spacing. If geoCoord==True, also returns the
        geometric coordinates of the discretisation
        """
        elementEvaluator = self.makeElementEvaluator(elementNumber)
        xi = discretisation.discretiseRegularGeoD(maxDistance, elementEvaluator)
        if geoCoords:
            coord = elementEvaluator.eval(xi)
            return xi, coord
        else:
            return xi, None

    def discretiseAllElementsRegularGeoD(self, maxDistance, geoCoords=False, unpack=True):

        elementOutput = []
        for elementNumber in numpy.sort(list(self.ensemble_field_function.mesh.elements.keys())):
            if self.ensemble_field_function.mesh.elements[elementNumber].is_element == True:
                elementOutput.append(self.discretiseElementRegularGeoD(elementNumber, maxDistance, geoCoords=geoCoords))
            else:
                g = self.makeSubfieldGF(elementNumber)
                elementOutput.append(
                    g.discretiseAllElementsRegularGeoD(maxDistance, geoCoords=geoCoords, unpack=unpack))

        xi = [e[0] for e in elementOutput]
        geo = [e[1] for e in elementOutput]

        if geoCoords and unpack:
            geo = numpy.vstack(geo)
            xi = numpy.vstack(xi)
        elif not geoCoords:
            geo = None
            if unpack:
                xi = numpy.vstack(xi)

        return xi, geo

    def generateSurfacePointGrid(self, spacing):
        # ~ self.flatten_ensemble_field_function()
        dataXi, dataCoords = self.discretiseAllElementsRegularGeoD(spacing, geoCoords=True, unpack=False)
        # dataXi, dataCoords = misc._removeDuplicates( dataXi, dataCoords )

        if self.ensemble_field_function.is_flat():
            dataXi, dataCoords = misc._removeDuplicatesFlat(dataXi, dataCoords)
        else:
            dataXi, dataCoords = misc._removeDuplicates(dataXi, dataCoords)

        dataCoords = numpy.vstack([numpy.vstack(c) for c in dataCoords])

        dataXi = misc.dictXiMap(dataXi)
        if self.ensemble_field_function.is_flat():
            regionDataMap = None
        else:
            regionDataMap = misc.makeRegionEPMap(dataXi)

        # ~ self.ensemble_field_function = self.ensemble_field_function_old

        return dataCoords, regionDataMap, dataXi

    def triangulate(self, GD, merge=True, retVertMap=False):
        """Create a triangulated discretisation of the geometric field.
        Inputs:
        GD: 2-tuple of the element dicretisation in each xi direction
        merge: Boolean, whether to connect triangles on element boundaries
        retVertMap: Boolean, return uniqueVertexIndices and vertMap

        Returns:
        P: (nx3) array of vertex coordinates
        T: (mx3) array of face indices
        uniqueVertexIndices: [optional] a list of vertex indices from the 
            original discretisation that is in the merged triangulation
        vertMap: [optional] a dictionary mapping original vertex indices
            to merged vertex indices
        """
        self.flatten_ensemble_field_function()
        P = self.evaluate_geometric_field(GD).T
        T = self.triangulator._triangulate(GD)
        if merge:
            P, T, uniqueVertexIndices, vertMap = self.triangulator._mergePoints2(P)

        T = T[:, ::-1]  # reverse ordering so normal points out

        self.ensemble_field_function = self.ensemble_field_function_old

        if retVertMap:
            return P, T, uniqueVertexIndices, vertMap
        else:
            return P, T

    # ==================================================================#
    # volume functions                                                   #
    # ==================================================================#
    def isInteriorToSurface(self, points, GD=2.0, maxOutDist=None, exactSearch=False):
        """
        classify n points as either inside or outside the GF surface. A 
        point P is interior if dot(PX, N_X)<0, where X is the closest 
        point to P on the GF, and N_X is the normal at X.
        
        if maxOutDist is not None, it should be a float. if maxOutDist 
        is +ve, points outside within distance maxOutDist will be 
        included. If maxOutDist is less than zero, points interior less
        than -maxOutDist from the surface will be excluded.
        
        return an n long binary array. 1 = interior, 0 = exterior
        """

        self.flatten_ensemble_field_function()

        if exactSearch:
            GFXi, GFX, closestDist = self.find_closest_material_points(points, initGD=[80, 80])
            NX = [self.evaluate_normal_in_mesh(None, elemXi={Xi[0]: Xi[1]}).squeeze() for Xi in GFXi]
            NX = numpy.array(NX)
            PX = points - GFX
        else:
            # evaluate points on surface
            GFXi, GFX = self.discretiseAllElementsRegularGeoD(GD, geoCoords=True, unpack=False)
            GFX = numpy.vstack(GFX)
            GFElemXi = dict(list(zip(numpy.sort(list(self.ensemble_field_function.mesh.elements.keys())), GFXi)))

            # evaluate normals of these points
            GFNormals = self.evaluate_normal_in_mesh(None, elemXi=GFElemXi).T

            # kdtree search to find closest GFX for each P
            GFXTree = cKDTree(GFX)
            closestDist, Xi = GFXTree.query(list(points), k=1)

            # calculate P->GFX vectors
            PX = points - GFX[Xi]
            NX = GFNormals[Xi]

        # dot product normals with P->GFX
        dot_ = (NX * PX).sum(1)

        # return 1 for interior (-ve), 0 for exterior (+ve)
        mask = (dot_ < 0.0)

        # outside points within tolerance
        if maxOutDist is not None:
            if maxOutDist > 0.0:
                # find all points within maxOutDist of closest GFX
                mask = mask | (closestDist < maxOutDist)
            else:
                mask = mask & (closestDist > abs(maxOutDist))

        self.ensemble_field_function = self.ensemble_field_function_old

        return mask

    def generateInternalPointsGrid(self, spacing, maxOutDist=None, exactClosestSearch=False):
        """
        generate a grid of points internal to surface mesh with spacing 
        specified by tuple spacing.
        """

        # sample surface to get bounding box
        s = self.evaluate_geometric_field([15, 15]).T
        if maxOutDist is None:
            bboxMin = s.min(0)
            bboxMax = s.max(0)
        else:
            bboxMin = s.min(0) - maxOutDist
            bboxMax = s.max(0) + maxOutDist

        N = (bboxMax - bboxMin) / spacing

        # generate grid of points in bounding box
        PAll = numpy.array([[x, y, z] for z in numpy.linspace(bboxMin[2], bboxMax[2], N[2]) \
                            for y in numpy.linspace(bboxMin[1], bboxMax[1], N[1]) \
                            for x in numpy.linspace(bboxMin[0], bboxMax[0], N[0])])

        # filter out exterior points
        isInterior = self.isInteriorToSurface(PAll, maxOutDist=maxOutDist, exactSearch=exactClosestSearch)

        return PAll[isInterior, :]


# ======================================================================#
class geometric_point(object):
    """ Class of point objects used by geometric_field. Assumes field parameters
    are given as a 2D array with shape = (dimensions, number of parameters per dimensions)
    """

    def __init__(self, dimensions, field_parameters=None, ensemble_point_number=None):

        self.field_parameters = None
        self.ensemble_point_number = None

        self.dimensions = dimensions
        if field_parameters is not None:
            if len(field_parameters) != dimensions:
                raise ValueError('ERROR: geometric_point.__init__: dimension and number of coordinates mismatch')
            else:
                self.field_parameters = field_parameters.copy()
        if ensemble_point_number:
            self.ensemble_point_number = ensemble_point_number

        return

    # ==================================================================#
    def get_dimensions(self):
        return self.dimensions

    # ==================================================================#
    def set_field_parameters(self, field_parameters):
        """ field_parameters should be a 2D array with
        shape = (dimensions, number of parameters per dimensions)
        """

        if len(field_parameters) != self.dimensions:
            raise ValueError(
                'ERROR: geometric_point.set_field_parameters: dimension and number of coordinates mismatch')

        self.field_parameters = field_parameters.copy()
        return 1

    # ==================================================================#
    def get_field_parameters(self):
        return self.field_parameters.copy()

    # ==================================================================#
    def set_ensemble_point_number(self, n):
        self.ensemble_point_number = n

    # ==================================================================#
    def get_ensemble_point_number(self):
        return self.ensemble_point_number

    # ==================================================================#
    def get_position(self):
        return numpy.array(self.field_parameters[:, 0])


# ======================================================================#
class spline_3d_parametric(object):
    """ Class for a spline curve fitted to a number of points in 3D for
    defining landmarks. Interpolates x,y,z given parametric variable u
    """

    def __init__(self, point_coords, smoothing=3.0, order=3, nknots=-1):
        """ point_coords is a list of N lists n long for n points in N-D 
        space.
        """
        self.point_coords = point_coords
        self.smoothing = smoothing
        self.order = order
        self.nknots = nknots
        self.type = 'parametric'

        self.tck, self.u = splprep(point_coords, s=smoothing, k=order, nest=nknots)

    # ==================================================================#
    def eval(self, u):
        """ evaluate the spline g(u) at u. u can be a sequence of values
        """
        u = numpy.array(u)
        x = splev(u, self.tck)

        return x

    # ==================================================================#
    def findClosest(self, p):
        """ returns coords and parameters u of point on the line closest
        to p
        """
        self.p = tuple(p)
        u0 = [0.5]
        uMin = fmin(self._closestObj, u0, disp=0)
        pClosest = self.eval(uMin)
        self.p = None

        return pClosest, uMin

    def _closestObj(self, u):
        pLine = self.eval(u[0])
        d = numpy.sqrt((numpy.subtract(self.p, pLine) ** 2.0).sum())
        return d


# ======================================================================#
class spline_3d_interp_x(object):
    """ Class for a spline curve fitted to a number of points for
    defining landmarks. Interpolates y and z for given x. 
    x must be non-decreasing
    """

    def __init__(self, x, y, z, order=3, smooth=3.0):
        """ x,y,z are lists n long for n points in 3-D 
        space.
        """
        self.xb = x[0]
        self.xe = x[-1]
        self.tcky = splrep(x, y, k=order, s=smooth)
        self.tckz = splrep(x, z, k=order, s=smooth)
        self.type = 'interp_x'

    # ==================================================================#
    def eval(self, x):
        """ returns interpolated y and z values for a list of x values
        """
        newy = splev(x, self.tcky)
        newz = splev(x, self.tckz)

        return (newy, newz)


# ======================================================================#


def smoothCurvField(points, curvature):
    # ~ return curvature
    # ~ return CT.smoothCurvField1( points, curvature )
    return CT.smoothCurvField2(points, curvature)


# =Optimised evaluators=================================================#

def makeGeometricFieldEvaluator(G, evalD):
    # G should be flat

    f = G.ensemble_field_function

    # calculate static basis values for the required evalD
    basisValues = {}
    basisMatrices = {}
    for elementNumber in numpy.sort(list(f.mesh.elements.keys())):
        # get element
        element = f.mesh.elements[elementNumber]
        # calculate basis values
        if basisValues.get(element.type) is None:
            if isinstance(evalD, float):
                evalGrid = G.discretiseElementRegularGeoD(elementNumber, evalD)[0]
            else:
                evalGrid = element.generate_eval_grid(evalD).squeeze()
            basisValues[element.type] = f.basis[element.type].eval(evalGrid.T).T

        basisMatrices[elementNumber] = basisValues[element.type]

    def evaluator(P):
        P3 = P.reshape((3, -1, 1))
        V = []
        for elementNumber in numpy.sort(list(f.mesh.elements.keys())):
            # map global parameters to element parameters
            p1 = f.mapper.get_element_parameters(elementNumber, P3[0], do_hack=1)
            p2 = f.mapper.get_element_parameters(elementNumber, P3[1], do_hack=1)
            p3 = f.mapper.get_element_parameters(elementNumber, P3[2], do_hack=1)
            # combine with static basis values 
            V.append(numpy.dot(basisMatrices[elementNumber], numpy.vstack((p1, p2, p3)).T))

        return numpy.vstack(V).T

    return evaluator


def makeGeometricFieldEvaluatorSparse(G, evalD, epIndex=None, epXi=None, matPoints=None):
    """ create a function for evaluation the geometric field values,
    taking advantage of a precomputed sparse matrix of basis function
    values at fixed element coordinates
    
    This is about 10x faster than the same setup with a dense matrix
    implementation. csc is a little bit faster than csr sparse matrix.
    
    if evalD is a float, regular geometric discretisation is assume, else
    evalD should be a list of number of discretised points in each xi
    direction
    
    epIndex, if defined, is a list of the row numbers (ep numbers) that
    should be in the final system. All other rows are droppped
    
    epXi, if defined, with evalD=None, is a list of lists of
    element xi coordinates at which the field is to be evaluated
    """

    f = G.ensemble_field_function
    if not f.is_flat():
        f = f.flatten()[0]

    if matPoints is not None:
        epMode = 3
        nEPs = len(matPoints)
    elif evalD is None:
        ep = epXi
        nEPs = numpy.sum([e.shape[0] for e in ep])
        # ~ nEPs = numpy.sum( [e.shape[0] for e in ep.values()] )
        epMode = 1
    elif isinstance(evalD, float):
        ep = G.discretiseAllElementsRegularGeoD(evalD, unpack=False)[0]
        nEPs = numpy.sum([e.shape[0] for e in ep])
        epMode = 1
    else:
        ep = G.evaluate_geometric_field(evalD)
        nEPs = ep.shape[1]
        epMode = 2

    d = G.dimensions
    A = numpy.zeros((nEPs, f.get_number_of_ensemble_points()), dtype=float)

    if epMode == 3:
        # ~ pdb.set_trace()
        elemEnsNodes = {}
        for mpI, (elem, xi) in enumerate(matPoints):
            element = f.mesh.elements[elem]
            b = f.basis[element.type].eval(xi)
            ensNodes = elemEnsNodes.get(elem)
            if ensNodes is None:
                emap = f.mapper._element_to_ensemble_map[elem]
                ensNodes = [emap[k][0][0] for k in list(emap.keys())]
                elemEnsNodes[elem] = ensNodes

            A[mpI, ensNodes] = b

    else:
        # calculate static basis values for the required evalD and assemble
        # matrix
        basisValues = {}
        row = 0

        for ei, elementNumber in enumerate(numpy.sort(list(f.mesh.elements.keys()))):
            # get element
            element = f.mesh.elements[elementNumber]

            if epMode == 1:
                # regular geom discretisation
                evalGrid = ep[ei]
                b = f.basis[element.type].eval(evalGrid.T).T  # basis values

                # dictionary epXi
                # ~ evalGrid = ep.get(ei)
                # ~ if evalGrid is None:
                # ~ continue
                # ~ else:
                # ~ b = f.basis[element.type].eval( evalGrid.T ).T # basis values
            elif epMode == 2:
                # regular xi discretisation
                if basisValues.get(element.type) is None:
                    evalGrid = element.generate_eval_grid(evalD).squeeze()
                    basisValues[element.type] = f.basis[element.type].eval(evalGrid.T).T

                b = basisValues[element.type]  # basis values

            # fill in A matrix                                
            emap = f.mapper._element_to_ensemble_map[elementNumber]  # element to ensemble map
            for n in range(b.shape[1]):
                A[row:row + b.shape[0], emap[n][0][0]] = b[:, n]

            row += b.shape[0]

        # if epIndex is defined, remove rows not defined in epIndex
        # ~ pdb.set_trace()
        if epIndex is not None:
            A = A[epIndex, :]

    # sparsify A
    As = sparse.csc_matrix(A)

    def evaluator(P):
        E = As * P.reshape((d, -1)).T
        return E.T

    return evaluator


def makeGeometricFieldElementsEvaluatorSparse(gf, elems, evalD):
    """Create a function to efficiently evaluate given elements
    in a geometric_field at fixed material points defined by 
    a discretisation.

    Inputs:
    gf [geometric_field]: a geometric_field instance
    elems [list]: list of element numbers
    evalD [float or list]: discretisation of the elements
    """
    elemXis = []
    for elem in elems:
        elemXi = gf.discretiseElementRegularGeoD(
            elem, evalD, geoCoords=False
        )[0]
        elemXis += [(elem, xi) for xi in elemXi]
    elemsEvaluator = makeGeometricFieldEvaluatorSparse(
        gf, None, matPoints=elemXis
    )
    return elemsEvaluator


def buildEvaluatorSparse(A, entries, nodeMap, dimensions):
    """
    A = matrix to be filled
    entries = [ [elemN, basisValues], ...]
    
    not used?
    """
    row = 0
    for eN, b in entries:
        ensemblePoints = [nodeMap[eN][n][0][0] for n in range(b.shape[0])]  # element to ensemble map
        A[row, ensemblePoints] = b
        row += 1

    # sparsify A
    As = sparse.csc_matrix(A)

    def evaluator(P):
        E = As * P.reshape((dimensions, -1)).T
        return E.T

    return evaluator


def makeGeometricFieldDerivativesEvaluatorSparse(G, evalD, dim=3, epIndex=None, epXi=None):
    """ create a function for evaluating the geometric field derivatives,
    taking advantage of a precomputed sparse matrix of basis function
    values at fixed element coordinates
    
    This is about 10x faster than the same setup with a dense matrix
    implementation. csc is a little bit faster than csr sparse matrix.
    
    if evalD is a float, regular geometric discretisation is assume, else
    evalD should be a list of number of discretised points in each xi
    direction
    
    epIndex, if defined, is a list of the row numbers (ep numbers) that
    should be in the final system. All other rows are droppped
    """
    f0 = G.ensemble_field_function
    f = f0.flatten()[0]
    G.ensemble_field_function = f

    if evalD is None:
        ep = epXi
        nEPs = numpy.sum([e.shape[0] for e in ep])
        epMode = 1
    elif isinstance(evalD, float):
        ep = G.discretiseAllElementsRegularGeoD(evalD, unpack=False)[0]
        nEPs = numpy.sum([e.shape[0] for e in ep])
        epMode = 1
    else:
        ep = G.evaluate_geometric_field(evalD)
        nEPs = ep.shape[1]
        epMode = 2

    nDerivs = int(G.ensemble_field_function.dimensions ** 2 + 1)

    # calculate static basis values for the required evalD and assemble
    # matrices
    basisValues = {}
    A = [numpy.zeros((nEPs, f.get_number_of_ensemble_points()), dtype=float) for i in range(nDerivs)]
    row = 0
    for ei, elementNumber in enumerate(numpy.sort(list(f.mesh.elements.keys()))):
        # get element
        element = f.mesh.elements[elementNumber]

        # calculate basis values
        if epMode == 1:
            # ~ evalGrid = G.discretiseElementRegularGeoD( elementNumber, evalD )[0]
            evalGrid = ep[ei]
            b = f.basis[element.type].eval_derivatives(evalGrid.T, None)
        elif epMode == 2:
            if basisValues.get(element.type) is None:
                evalGrid = element.generate_eval_grid(evalD)
                basisValues[element.type] = f.basis[element.type].eval_derivatives(evalGrid.T, None)

            b = basisValues[element.type]

        # fill in A matrix
        emap = f.mapper._element_to_ensemble_map[elementNumber]  # element to ensemble map
        for d in range(b.shape[0]):
            for n in range(b.shape[1]):
                A[d][row:row + b.shape[2], emap[n][0][0]] = b[d, n, :]


        row += b.shape[2]

    if epIndex is not None:
        AStacked = []
        for a in A:
            AStacked.append(a[:, epIndex])
        AStacked = numpy.vstack(AStacked)
    else:
        AStacked = numpy.vstack(A)

    # stack all derivative A matrices
    AStackedSparse = sparse.csc_matrix(AStacked)
    nDerivs = len(A)

    def evaluator(P):
        """ uses a A matrix that is the vstack of all derivative A matrices
        """

        Pd = P.reshape((dim, -1)).T
        D = AStackedSparse * Pd
        return D.T.reshape((dim, nDerivs, -1))

    return evaluator


# =============================================================================#
# arc length evaluation

def makeArclengthEvalDisc(c, d):
    """
    Return a function for evaluating all element arclengths in geometric_field
    c by discretisation of each element into d linear segments.
    """

    ceval = makeGeometricFieldEvaluatorSparse(c, [d, ])
    p = numpy.array(c.field_parameters)
    n_elems = c.ensemble_field_function.mesh.get_number_of_true_elements()

    def f(p):
        x = ceval(p).T

        # separate x into element points per element
        _x = x.reshape((n_elems, -1, 3))
        return numpy.sum(
            numpy.sqrt(
                ((_x[:, 1:, :] - _x[:, :-1, :]) ** 2.0).sum(2)
            ), 1)

    return f


# =============================================================================#
# serialisation

def load_gf_shelve(filename, G, ensfn=None, meshfn=None, filedir=None, force=False):
    """
    Deserialise a geometric_field using the shelve package.

    Inputs:
    filename: [str] geometric_field filename (.geof).
    G: a geometric_field instance
    ensfn: [str] ensemble field function filename (.ens).
    meshfn: [str] mesh filename (.mesh).
    """

    if filedir is not None:
        filename = os.path.join(filedir, filename)

    try:
        try:
            S = shelve.open(filename, 'r')
        except ImportError:
            import bsddb3
            _db = bsddb3.hashopen(filename)
            S = shelve.Shelf(_db)
    except:
        raise IOError(filename + ' could not be read')
    else:
        if ensfn:
            F = EFF.load_ensemble(ensfn, meshFilename=meshfn, path=filedir)
        elif S.get('ensemble_field'):
            try:
                F = EFF.load_ensemble(S['ensemble_field'], meshFilename=meshfn, path=filedir)
            except IOError:
                if force:
                    log.debug('WARNING: no ensemble field function loaded, ignoring.')
                    F = None
        else:
            log.debug('WARNING: no ensemble field function loaded')
            F = None

        G.name = S['name']
        G.dimensions = S['dimensions']
        G.ensemble_point_counter = S['ensemble_point_counter']

        if F is not None:
            G.ensemble_field_function = F
            G._create_ensemble_points()
            G.triangulator.f = F
            G.set_field_parameters(S['field_parameters'])
        else:
            G.field_parameters = S['field_parameters']

        return G


class GeometricFieldJSONWriter(object):
    node_ptr = 'node {:08d}'
    dim_ptr = 'dim {:1d}'
    num_pattern = '{:20.16E}'

    def __init__(self, gf):
        """
        Writer class for serialising a geometric field to a
        JSON format file or string.
        """
        self.gf = gf
        self._file_dir = ''

    def write(self, filename, ensfn=None, meshfn=None, filedir=None):
        """
        Write to file a json format file of the eff properties.
        """
        if filedir is not None:
            self._file_dir = filedir

        d = self.serialise(ensfn, meshfn)
        with open(filename, 'w') as f:
            json.dump(d, f, indent=4, sort_keys=True)

        return filename

    def serialise(self, ensfn, meshfn):
        """
        Return a json-compatible dict of the eff properties.
        """
        d = {}
        self._serialise_meta(d)
        self._serialise_ens(d, ensfn, meshfn)
        self._serialise_field_parameters(d)
        return d

    def _serialise_meta(self, gf_dict):
        gf_dict['name'] = self.gf.name
        gf_dict['dimensions'] = self.gf.dimensions
        gf_dict['ensemble_point_counter'] = self.gf.ensemble_point_counter

    def _serialise_ens(self, gf_dict, ensfn, meshfn):
        if ensfn is not None:
            self.gf.ensemble_field_function.save_ensemble(
                ensfn, mesh_filename=meshfn, path=self._file_dir
            )
            gf_dict['ensemble_field'] = os.path.split(ensfn)[1]
        else:
            pass
            # gf_dict['ensemble_field'] = self.gf.ensemble_field_function.save_ensemble(
            #     self.gf.ensemble_field_function.name,
            #     mesh_filename=meshfn,
            #     path=self._file_dir
            #     )

    def _serialise_field_parameters(self, gf_dict):
        d = {}
        node_numbers = numpy.arange(self.gf.field_parameters.shape[1])
        dims = numpy.arange(self.gf.field_parameters.shape[0])
        params = self.gf.field_parameters.squeeze()

        for n in node_numbers:
            node_line = self.node_ptr.format(n)
            d[node_line] = {}
            for dim in dims:
                d[node_line][self.dim_ptr.format(dim)] = ' '.join(
                    [self.num_pattern.format(x) for x in self.gf.field_parameters[dim, n, :]])

        gf_dict['field_parameters'] = d


class GeometricFieldJSONReader(object):
    node_ptr = 'node {:08d}'
    dim_ptr = 'dim {:1d}'
    num_pattern = '{:20.16E}'

    def __init__(self, gf):
        """
        Reader class for setting a geometric field with properties
        from a JSON format file or string.
        """
        self.gf = gf
        self._file_dir = ''
        self.force = False

    def read(self, filename, ensfn=None, meshfn=None, filedir=None, force=False):
        """
        Load gf properties from a gf file.
        """
        self.force = force
        if filedir is not None:
            filename = os.path.join(filedir, filename)
            self._file_dir = filedir
        # else:
        #     self._file_dir = os.path.split(filename)[0]
        with open(filename, 'r') as f:
            gf_dict = json.load(f)

        self.deserialise(gf_dict, ensfn, meshfn)

    def deserialise(self, jsonstr, ensfn, meshfn):
        """
        Set self.gf with properties from the json string.
        file_dir is the directory of any embedded gf files.
        """
        self._parse_meta(jsonstr)
        self._parse_ens(jsonstr, ensfn, meshfn)
        self._parse_field_parameters(jsonstr)

    def _parse_meta(self, gf_dict):
        self.gf.name = gf_dict['name']
        self.gf.dimensions = int(gf_dict['dimensions'])
        self.gf.ensemble_point_counter = int(gf_dict['ensemble_point_counter'])

    def _parse_ens(self, gf_dict, ensfn, meshfn):
        if ensfn is not None:
            self.gf.ensemble_field_function = EFF.load_ensemble(ensfn, meshfn, path=self._file_dir)
        elif gf_dict.get['ensemble_field']:
            try:
                self.gf.ensemble_field_function = EFF.load_ensemble(gf_dict['ensemble_field'], meshfn,
                                                                    path=self._file_dir)
            except IOError:
                if self.force:
                    log.debug('Cannot open {}, ignoring'.format(gf_dict['ensemble_field']))
                    pass
                else:
                    raise IOError('Cannot open {}'.format(gf_dict['ensemble_field']))
        else:
            log.debug('WARNING: no ensemble field function loaded')

        if self.gf.ensemble_field_function is not None:
            self.gf.triangulator.f = self.gf.ensemble_field_function
            self.gf._create_ensemble_points()

    def _parse_field_parameters(self, gf_dict):
        node_numbers = [int(k.split(' ')[1]) for k in gf_dict['field_parameters'].keys()]
        node_keys = sorted(gf_dict['field_parameters'].keys())
        dim_keys = sorted(gf_dict['field_parameters'][node_keys[0]].keys())
        val_len = gf_dict['field_parameters'][node_keys[0]][dim_keys[0]].split(' ').__len__()
        p = numpy.zeros([len(dim_keys), len(node_keys), val_len], dtype=float)
        for ni, nk in enumerate(node_keys):
            node_dict = gf_dict['field_parameters'][nk]
            for di, dk in enumerate(dim_keys):
                p[di, ni, :] = [float(x) for x in node_dict[dk].split(' ')]

        if self.gf.ensemble_field_function is not None:
            self.gf.set_field_parameters(p)
        else:
            self.gf.field_parameters = p


def load_gf_json(filename, gf, ensfn=None, meshfn=None, filedir=None, force=False):
    reader = GeometricFieldJSONReader(gf)
    reader.read(filename, ensfn, meshfn, filedir, force=force)
    return gf


def save_gf_json(filename, gf, ensfn=None, meshfn=None, filedir=None):
    writer = GeometricFieldJSONWriter(gf)
    writer.write(filename, ensfn, meshfn, filedir)


def load_geometric_field(filename, ensFilename=None, meshFilename=None, path=None, force=False):
    """
    Deserialise a geometric_field from either a shelve file or a json file.

    Inputs:
    filename: [str] geometric_field filename (.geof).
    ensFilename: [str] ensemble field function filename (.ens).
    meshFilename: [str] mesh filename (.mesh).
    path: [str] path of the directory of the above files. If defined, the filenames
        above should not include the path to the file.
    """
    gf = geometric_field('none', 1)
    with open(filename, 'r') as f:
        head = f.read(1)
        if head == '{':
            load_gf_json(
                filename, gf, ensfn=ensFilename, meshfn=meshFilename,
                filedir=path, force=force,
            )
        else:
            load_gf_shelve(
                filename, gf, ensfn=ensFilename, meshfn=meshFilename,
                filedir=path, force=force,
            )

    return gf
