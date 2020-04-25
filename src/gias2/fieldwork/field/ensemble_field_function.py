"""
FILE: ensemble_field_function.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION:
Classes for combining information and functions of mesh topology
(topology.mesh.mesh_ensemble, element.element_types), basis 
(basis.basis), and mapping (mapper.mapper)

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import copy
import json
import os
import shelve
import logging
import numpy

from gias2.fieldwork.field import mapper
from gias2.fieldwork.field.basis import basis
from gias2.fieldwork.field.topology import element_types
from gias2.fieldwork.field.topology import mesh

log = logging.getLogger(__name__)
# ~ from field_ctypes import ensemble_evaluators

def eval_line_L3(B, P):
    E = B[0] * P[0] + B[1] * P[1] + B[2] * P[2] + B[3] * P[3]
    return E


def eval_line_L4(B, P):
    E = B[0] * P[0] + B[1] * P[1] + B[2] * P[2] + B[3] * P[3] + B[4] * P[4]
    return E


def eval_simplex_L3_L3(B, P):
    E = B[0] * P[0] + B[1] * P[1] + B[2] * P[2] + B[3] * P[3] + B[4] * P[4] + \
        B[5] * P[5] + B[6] * P[6] + B[7] * P[7] + B[8] * P[8] + B[9] * P[9]

    return E


def eval_simplex_L4_L4_old(B, P):
    return numpy.dot(P, B)


def eval_simplex_L4_L4(B, P):
    E = B[0] * P[0] + B[1] * P[1] + B[2] * P[2] + B[3] * P[3] + \
        B[4] * P[4] + B[5] * P[5] + B[6] * P[6] + B[7] * P[7] + \
        B[8] * P[8] + B[9] * P[9] + B[10] * P[10] + B[11] * P[11] + \
        B[12] * P[12] + B[13] * P[13] + B[14] * P[14]

    return E


# ~ eval_simplex_L4_L4 = ensemble_evaluators.eval

def eval_quad_L3_L3(B, P):
    E = B[0] * P[0] + B[1] * P[1] + B[2] * P[2] + B[3] * P[3] + \
        B[4] * P[4] + B[5] * P[5] + B[6] * P[6] + B[7] * P[7] + \
        B[8] * P[8] + B[9] * P[9] + B[10] * P[10] + B[11] * P[11] + \
        B[12] * P[12] + B[13] * P[13] + B[14] * P[14] + B[15] * P[15]

    return E


def eval_quad_L4_L4_old(B, P):
    E = B[0] * P[0] + B[1] * P[1] + B[2] * P[2] + B[3] * P[3] + B[4] * P[4] + \
        B[5] * P[5] + B[6] * P[6] + B[7] * P[7] + B[8] * P[8] + B[9] * P[9] + \
        B[10] * P[10] + B[11] * P[11] + B[12] * P[12] + B[13] * P[13] + B[14] * P[14] + \
        B[15] * P[15] + B[16] * P[16] + B[17] * P[17] + B[18] * P[18] + B[19] * P[19] + \
        B[20] * P[20] + B[21] * P[21] + B[22] * P[22] + B[23] * P[23] + B[24] * P[24]

    return E


# ~ eval_quad_L4_L4 = ensemble_evaluators.eval

def eval_quad_L2_L2_L2(B, P):
    E = B[0] * P[0] + B[1] * P[1] + B[2] * P[2] + B[3] * P[3] + B[4] * P[4] + \
        B[5] * P[5] + B[6] * P[6] + B[7] * P[7] + B[8] * P[8] + B[9] * P[9] + \
        B[10] * P[10] + B[11] * P[11] + B[12] * P[12] + B[13] * P[13] + B[14] * P[14] + \
        B[15] * P[15] + B[16] * P[16] + B[17] * P[17] + B[18] * P[18] + B[19] * P[19] + \
        B[20] * P[20] + B[21] * P[21] + B[22] * P[22] + B[23] * P[23] + B[24] * P[24] + \
        B[25] * P[25] + B[26] * P[26]

    return E


def eval_quad_L3_L3_L3(B, P):
    E = B[0] * P[0] + B[1] * P[1] + B[2] * P[2] + B[3] * P[3] + B[4] * P[4] + \
        B[5] * P[5] + B[6] * P[6] + B[7] * P[7] + B[8] * P[8] + B[9] * P[9] + \
        B[10] * P[10] + B[11] * P[11] + B[12] * P[12] + B[13] * P[13] + B[14] * P[14] + \
        B[15] * P[15] + B[16] * P[16] + B[17] * P[17] + B[18] * P[18] + B[19] * P[19] + \
        B[20] * P[20] + B[21] * P[21] + B[22] * P[22] + B[23] * P[23] + B[24] * P[24] + \
        B[25] * P[25] + B[26] * P[26] + B[27] * P[27] + B[28] * P[28] + B[29] * P[29] + \
        B[30] * P[30] + B[31] * P[31] + B[32] * P[32] + B[33] * P[33] + B[34] * P[34] + \
        B[35] * P[35] + B[36] * P[36] + B[37] * P[37] + B[38] * P[38] + B[39] * P[39] + \
        B[40] * P[40] + B[41] * P[41] + B[42] * P[42] + B[43] * P[43] + B[44] * P[44] + \
        B[45] * P[45] + B[46] * P[46] + B[47] * P[47] + B[48] * P[48] + B[49] * P[49] + \
        B[50] * P[50] + B[51] * P[51] + B[52] * P[52] + B[53] * P[53] + B[54] * P[54] + \
        B[55] * P[55] + B[56] * P[56] + B[57] * P[57] + B[58] * P[58] + B[59] * P[59] + \
        B[60] * P[60] + B[61] * P[61] + B[62] * P[62] + B[63] * P[63]

    return E


def eval_quad_L4_L4_L4(B, P):
    E = B[0] * P[0] + B[1] * P[1] + B[2] * P[2] + B[3] * P[3] + B[4] * P[4] + \
        B[5] * P[5] + B[6] * P[6] + B[7] * P[7] + B[8] * P[8] + B[9] * P[9] + \
        B[10] * P[10] + B[11] * P[11] + B[12] * P[12] + B[13] * P[13] + B[14] * P[14] + \
        B[15] * P[15] + B[16] * P[16] + B[17] * P[17] + B[18] * P[18] + B[19] * P[19] + \
        B[20] * P[20] + B[21] * P[21] + B[22] * P[22] + B[23] * P[23] + B[24] * P[24] + \
        B[25] * P[25] + B[26] * P[26] + B[27] * P[27] + B[28] * P[28] + B[29] * P[29] + \
        B[30] * P[30] + B[31] * P[31] + B[32] * P[32] + B[33] * P[33] + B[34] * P[34] + \
        B[35] * P[35] + B[36] * P[36] + B[37] * P[37] + B[38] * P[38] + B[39] * P[39] + \
        B[40] * P[40] + B[41] * P[41] + B[42] * P[42] + B[43] * P[43] + B[44] * P[44] + \
        B[45] * P[45] + B[46] * P[46] + B[47] * P[47] + B[48] * P[48] + B[49] * P[49] + \
        B[50] * P[50] + B[51] * P[51] + B[52] * P[52] + B[53] * P[53] + B[54] * P[54] + \
        B[55] * P[55] + B[56] * P[56] + B[57] * P[57] + B[58] * P[58] + B[59] * P[59] + \
        B[60] * P[60] + B[61] * P[61] + B[62] * P[62] + B[63] * P[63] + B[64] * P[64] + \
        B[65] * P[65] + B[66] * P[66] + B[67] * P[67] + B[68] * P[68] + B[69] * P[69] + \
        B[70] * P[70] + B[71] * P[71] + B[72] * P[72] + B[73] * P[73] + B[74] * P[74] + \
        B[75] * P[75] + B[76] * P[76] + B[77] * P[77] + B[78] * P[78] + B[79] * P[79] + \
        B[80] * P[80] + B[81] * P[81] + B[82] * P[82] + B[83] * P[83] + B[84] * P[84] + \
        B[85] * P[85] + B[86] * P[86] + B[87] * P[87] + B[88] * P[88] + B[89] * P[89] + \
        B[90] * P[90] + B[91] * P[91] + B[92] * P[92] + B[93] * P[93] + B[94] * P[94] + \
        B[95] * P[95] + B[96] * P[96] + B[97] * P[97] + B[98] * P[98] + B[99] * P[99] + \
        B[100] * P[100] + B[101] * P[101] + B[102] * P[102] + B[103] * P[103] + B[104] * P[104] + \
        B[105] * P[105] + B[106] * P[106] + B[107] * P[107] + B[108] * P[108] + B[109] * P[109] + \
        B[110] * P[110] + B[111] * P[111] + B[112] * P[112] + B[113] * P[113] + B[114] * P[114] + \
        B[115] * P[115] + B[116] * P[116] + B[117] * P[117] + B[118] * P[118] + B[119] * P[119] + \
        B[120] * P[120] + B[121] * P[121] + B[122] * P[122] + B[123] * P[123] + B[124] * P[124]

    return E


def dot_Evaluator(B, P):
    return numpy.dot(P, B)


class ensemble_field_function:
    """ Class that combines basis, mesh and mapper to create a field.
    Acts as a wrapper for methods in mesh and mapper and handles field
    evaluation.
    """

    evaluators = {'line3l': dot_Evaluator, \
                  'line4l': eval_line_L3, \
                  'line5l': eval_line_L4, \
                  'tri10': eval_simplex_L3_L3, \
                  'tri6': dot_Evaluator, \
                  'tri6_ql': dot_Evaluator, \
                  'tri15': dot_Evaluator, \
                  'tri16': dot_Evaluator, \
                  'tri35': dot_Evaluator, \
                  'quad33': dot_Evaluator, \
                  'quad44': eval_quad_L3_L3, \
                  'quad54': dot_Evaluator, \
                  'quad55': dot_Evaluator, \
                  'quad333': eval_quad_L2_L2_L2, \
                  'quad444': eval_quad_L3_L3_L3, \
                  'quad555': eval_quad_L4_L4_L4, \
                  'quad542': dot_Evaluator, \
                  'quad552': dot_Evaluator, \
                  'prism15-5': dot_Evaluator, \
                  'prism15-2': dot_Evaluator, \
                  'prism6-5': dot_Evaluator, \
                  }

    def __init__(self, name, dimensions, debug=0):
        self.name = name  # string
        self.dimensions = dimensions  # integer
        self.basis = {}  # basis object dictionary { element.type: <basis_obj> }
        self.basis_types = []  # string
        self.mesh = None  # mesh_ensemble object
        self.mapper = mapper.mapper()  # parameter mapper object
        self.parameters = None
        self.subfields = {}  # {element_number: subfield}
        # ~ self.submesh_map = {}            # {submesh: corresponding subfield}
        self.subfield_counter = 0
        self.is_element = False
        self.element_param_cache = {}  # caches parameters for elements

        self.debug = debug

    # ==================================================================#
    def save_ensemble(self, filename=None, mesh_filename=None, path='', subfieldFilenames=None):

        if not filename:
            filename = self.name

        if os.path.splitext(filename)[1].lower() != '.ens':
            filename = filename + '.ens'

        return save_eff_json(
            filename, self, meshfn=mesh_filename, filedir=path, subfieldfns=subfieldFilenames
        )

    def save_ensemble_shelve(self, filename=None, mesh_filename=None, path='', subfieldFilenames=None):
        raise DeprecationWarning('Shelve serialisation will be deprecated')

        if path:
            if path[-1] != '/':
                path += '/'

        if not filename:
            filename = self.name

        try:
            S = shelve.open(path + filename + '.ens')
        except ImportError:
            import bsddb3
            _db = bsddb3.hashopen(path + filename + '.ens')
            S = shelve.Shelf(_db)

        S['name'] = self.name
        S['dimensions'] = self.dimensions
        S['basis'] = dict([[i[0], i[1].type] for i in list(self.basis.items())])
        if mesh_filename:
            mesh_loc = self.mesh.save_mesh(mesh_filename, path=path)
        else:
            pass
            # mesh_loc = self.mesh.save_mesh( self.mesh.name, path=path )

        S['mesh'] = os.path.split(mesh_loc)[1]

        # dont save parameters
        # ~ S['parameters'] = self.parameters
        if subfieldFilenames:
            S['subfields'] = [(i[0], i[1]) for i in list(subfieldFilenames.items())]
        else:
            S['subfields'] = [(s[0], os.path.split(s[1].save_ensemble(path=path))[1]) for s in
                              list(self.subfields.items())]

        S['subfield_counter'] = self.subfield_counter
        S['debug'] = self.debug
        S['custom_map'] = self.mapper._custom_ensemble_order

        return os.path.relpath(path + filename + '.ens')

    def write_elements(self, filename, header=None):
        """Write elements to a text file. Format per line is:
        [element number] [element type] [node_1] [node_2] ...
        where node_x is the ensemble number of element node x
        """

        elem_numbers = sorted(self.mesh.elements.keys())

        with open(filename, 'w') as f:
            f.write('# Ensemble elements\n')
            f.write('# name: {}\n'.format(self.name))
            if header is not None:
                f.write('#' + header + '\n')
            for n in elem_numbers:
                elem_map = self.mapper._element_to_ensemble_map[n]
                elem_ens_points = [elem_map[k][0][0] for k in sorted(elem_map.keys())]
                line = '{:6d} {:s} '.format(n, self.mesh.elements[n].type) + ' '.join(
                    ['{:6d}'.format(p) for p in elem_ens_points]) + '\n'
                f.write(line)

    # ==================================================================#
    def set_basis(self, types):
        """ Instantiate basis function(s) for the specified element types. 
        
        Type is dictionary of {element_type: basis_type}. basis_type 
        is either a string or list of strings for multiple basis 
        functions. If a string, assume same basis type for all dimensions.
        """
        if not isinstance(types, dict):
            raise TypeError('input must be dictionary of {element_type: basis_type}')

        for e_type, b_type in list(types.items()):
            self.basis[e_type] = basis.make_basis(b_type)
            self.basis_types.append(self.basis[e_type].type)

    # ==================================================================#
    def set_new_mesh(self, name):
        """ Instantiate a new mesh_ensemble object.
        """

        self.mesh = mesh.mesh_ensemble(name, self.dimensions)
        if self.mesh:
            try:
                self.mesh_filename = self.mesh.filename
            except AttributeError:
                self.mesh_filename = self.mesh.name + '.mesh'

            return 1
        else:
            return

    # ==================================================================#
    def set_mesh(self, mesh):
        """ Set an existing mesh object to be the field topology
        automatically updates global_map. Mesh must be flat
        """
        if mesh.is_flat():
            self.mesh = mesh
            self.map_parameters()
            return 1
        else:
            log.debug('ERROR: ensemble_field_function.set_mesh: mesh must be flat')
            return

    # ==================================================================#
    def is_flat(self):
        if len(self.subfields) > 0:
            return False
        else:
            return True

    # ==================================================================#
    def add_element(self, element):
        """ Add an existing element (element or ensemble field function 
        object) to the mesh. returns new element number for convenience.
        """

        if element.is_element:
            n = self.mesh.add_element(element)
            if n is None:
                log.debug('ERROR: ensemble_field_function.add_element: unable to add element to mesh')
                return
        else:
            n = self._add_subfield(element)
            if n is None:
                log.debug('ERROR: ensemble_field_function.add_element: unable to add subfield to mesh')
                return

        return n

    # ==================================================================#
    def create_elements(self, type, number_of_elements):
        """ Create new elements of the specified type and add to mesh.
        Returns a tuple of new element numbers for convenience.
        """

        n_all = []
        for i in range(number_of_elements):
            # create element
            element = element_types.create_element(type)

            if element:
                n = self.mesh.add_element(element)
                if n is None:
                    log.debug('ERROR: ensemble_field_function.create_element: unable to add element to mesh')
                    return
                n_all.append(n)
            else:
                log.debug('ERROR ensemble_field_function.create_element: unable to create element')
                return

        return tuple(n_all)

    # ==================================================================#
    def _add_subfield(self, subfield):
        """ Add an existing ensemble field function object as a subfield
        update self.mesh connectivity.
        """

        # ~ self.submesh_map[ subfield.mesh ] = subfield
        # ~ self.subfields[self.subfield_counter] = subfield
        n = self.mesh.add_element(subfield.mesh)
        self.subfields[n] = subfield
        self.subfield_counter += 1
        return n

    # ==================================================================#
    def remove_element(self, element_number):
        """ Remove element element_number from the mesh.
        """

        self.mesh.remove_element(element_number)

        # removed element from the mapper
        self.mapper.set_number_of_ensemble_points(self.mesh.get_number_of_ensemble_points())
        self.mapper.remove_element(element_number)

        # remove from subfields dict in is a subfield
        try:
            del self.subfields[element_number]
        except KeyError:
            pass

        # self.subfield_counter -= 1

        return 1

    # ==================================================================#
    def connect_elements_by_edges(self, elem1, edge1, elem2, edge2, xi_extent=None, invert_order=False):
        self.mesh.connect_elements(elem1, edge1, elem2, edge2, xi_extent=xi_extent, invert_order=invert_order)

        # ==================================================================#

    def connect_element_points(self, points):
        """ Wrapper for mesh_ensemble's connect method.
        """

        if self.mesh.connect(points):
            return 1
        else:
            return

    # ==================================================================#s
    def connect_to_hanging_point(self, host_element, connected_points, mode, arg):
        """ Wrapper for mesh_ensemble's connect_to_hanging_point method.
        """

        if self.mesh.connect_to_hanging_point(host_element, connected_points, mode, arg):
            return 1
        else:
            return

    # ==================================================================#
    def map_parameters(self):
        """ Use mapper object to generates global node and elements 
        numbers using connectivity information.
        """

        self.mapper.set_parent_field(self)
        if self.mapper.do_mapping():
            return 1
        else:
            log.debug('ERROR: ensemble_field_function.map_parameters failed')
            return

    # ==================================================================#
    def get_mapping(self):
        """ Returns the mapping dictionaries generated by the mapper.
        """

        en_to_el = copy.deepcopy(self.mapper._ensemble_to_element_map)
        el_to_en = copy.deepcopy(self.mapper._element_to_ensemble_map)

        return (en_to_el, el_to_en)

    # ==================================================================#
    def get_number_of_ensemble_points(self):
        """ Returns the number of ensemble points in the field by
        querying the mapper.
        """

        return self.mesh.get_number_of_ensemble_points()

    # ==================================================================#
    def get_element_boundary_nodes(self, elemNumber):
        E = self.mesh.elements[elemNumber]
        if not E.is_element:
            raise ValueError('element ' + str(elemNumber) + ' is not a real element')

        # ~ boundaryNodes = []
        # ~ for edge in E.edge_points:
        # ~ boundaryNodes += list(edge)

        boundaryNodesEnsemble = []
        for b in E.edge_points:
            boundaryNodesEnsemble.append([self.mapper._element_to_ensemble_map[elemNumber][bi][0][0] for bi in b])

        return boundaryNodesEnsemble

    # ==================================================================#
    def flatten(self, debug=0):
        """ returns a new ensemble_field_function with a 1 level mesh
        made by flattening the current field mesh. Assumes that all sub-
        fields are of same dimension, basis and element shape.
        
        If mesh is already flat, do nothing
        """

        # don't do anything if already flat
        if self.is_flat():
            log.debug('Mesh is already flat, flattening aborted')
            return self, None, None

        # ==============================================================#
        # create new ESF with same dimensions and basis as self
        f_new = ensemble_field_function(self.name + '_flat', self.dimensions)
        f_new.set_basis(dict([[i[0], i[1].type] for i in list(self.basis.items())]))
        f_new.set_new_mesh(self.mesh.name + '_flat')

        # initialise maintaing a mapping of current global parameter
        # number to new element numbers
        paraMap = {}
        ELEMENT_MAP = {}  # maps current element number to new element number(s) { old_element_number: [ new_element_number,...] }
        ELEMENT_POINT_MAP = {}  # {old element point (e,p): new element point (e,p) }
        flattened_subfields = {}
        # ==============================================================#
        # add elements and connectivity to f_new
        # ~ for e_i in xrange( len(self.mesh.elements) ):
        for e_i in list(self.mesh.elements.keys()):

            E = self.mesh.elements[e_i]
            # ==========================================================#
            # for each current true element, add an element of the same type
            # to f_new, and connect to existing elements in f_new using
            # self.mapper information.
            if E.is_element:
                e_i_new = list(f_new.create_elements(E.type, 1))
                ELEMENT_MAP[e_i] = e_i_new
                for i in range(E.number_of_points):
                    ELEMENT_POINT_MAP[(e_i, i)] = [(e_i_new[0], i)]

                    # ==========================================================#
            # if element is another ESF, check if it is flat. If not,
            # call flatten on it.
            else:
                if E.is_flat():
                    sub = self.subfields[e_i]
                else:
                    sub = self.subfields[e_i].flatten()[0]

                flattened_subfields[e_i] = sub  # save flattened submesh for later use in connecting
                e_i_new = []  # store the element numbers of newly created elements (index is element number in sub)
                # create elements in f_new corresponding to those in sub
                sub_2_new_elem_map = {}  # records the mapping between subfield element numbers and new element numbers

                for e_i_sub in sub.mesh.elements:
                    e_i_new_i = f_new.create_elements(sub.mesh.elements[e_i_sub].type, 1)[0]
                    e_i_new.append(e_i_new_i)
                    sub_2_new_elem_map[e_i_sub] = e_i_new_i

                # a list of new element numbers of elements from submesh e_i
                ELEMENT_MAP[e_i] = e_i_new
                if debug:
                    log.debug(ELEMENT_MAP)
                    # ======================================================#
                # use sub.mesh.connectivity to connect up subfield elements             
                # connect non-hanging points first
                ignore = []  # list of element points already connected so can be ignored
                for k in list(sub.mesh.connectivity.keys()):
                    connected_points_new = []
                    if k[0] == -1:
                        pass
                    else:
                        # get sub-ensemble number (element point number for parent field) for this element point
                        sub_ensemble_i = sub.mapper._element_to_ensemble_map[k[0]][k[1]][0][0]
                        if sub.mapper.has_custom_map:
                            sub_ensemble_i = sub.mapper._custom_ensemble_order[sub_ensemble_i]

                        try:
                            ELEMENT_POINT_MAP[(e_i, sub_ensemble_i)].append((sub_2_new_elem_map[k[0]], k[1]))
                        except KeyError:
                            if debug:
                                log.debug('e_i:' + str(e_i) + ' e_i_new:' + str(e_i_new) + ' k:' + str(k))
                                log.debug(sub_2_new_elem_map[k[0]])
                                log.debug(k[1])
                            ELEMENT_POINT_MAP[(e_i, sub_ensemble_i)] = [(sub_2_new_elem_map[k[0]], k[1])]

                        # if element point was connected
                        # ~ if (sub.mesh.connectivity[k]) and (k not in ignore):
                        if (sub.mesh.connectivity[k]):
                            # gather connected element points and get their new element number from sub_elements
                            connected_points_new.append((sub_2_new_elem_map[k[0]], k[1]))
                            for elem_point in sub.mesh.connectivity[k]:
                                connected_points_new.append((sub_2_new_elem_map[elem_point[0]], elem_point[1]))

                            # connect
                            f_new.connect_element_points(connected_points_new)

                            ignore += connected_points_new

                # now deal with hanging points
                for hp_i in sub.mesh.hanging_points:
                    host_elem = elem_i_new[sub.mesh.hanging_points[hp_i].host_element]
                    host_elem_coords = sub.mesh.hanging_points[hp_i].element_coordinates
                    connected_points = []
                    for elem_point in sub.mesh.connectivity[(-1, hp_i)]:
                        connected_points.append((elem_i_new[elem_point[0]], elem_point[1]))
                    f_new.connect_to_hanging_point(host_elem, host_elem_coords, connected_points)

        # ==============================================================#
        # finalise connectivity between 0-level elements and flattened subfields
        # use self.mesh.connectivity to connect up subfield elements                
        # connect non-hanging points first
        ignore = []  # list of element points already connected so can be ignored
        for k in list(self.mesh.connectivity.keys()):

            connected_points = []
            if k[0] == -1:
                pass
            else:
                # if element point was connected
                # ~ if (self.mesh.connectivity[k]) and (k not in ignore):
                if (self.mesh.connectivity[k]):
                    # gather connected element points

                    # ~ pdb.set_trace()

                    # connected points in current mesh
                    connected_points_old = [k, ] + self.mesh.connectivity[k]

                    connected_points_new = []
                    for elem_point_old in connected_points_old:
                        # if this element point was in a subfield
                        if len(ELEMENT_MAP[elem_point_old[0]]) > 1:
                            # query a map to get the new element numbers in f_new
                            # that correspond to subfield ensemble numbers
                            subfield = flattened_subfields[elem_point_old[0]]

                            # if ensemble point numbers have been custom mapped, need to get original number
                            if subfield.mapper.has_custom_map:
                                elem_point_old = (
                                elem_point_old[0], subfield.mapper._custom_ensemble_order_inverse[elem_point_old[1]])

                            # query this subfield's _ensemble_to_element map to get the elements at the ensemble point (elem_point_old)
                            elem_point_map = subfield.mapper._ensemble_to_element_map[elem_point_old[1]]
                            for elem in elem_point_map:
                                connected_points_new.append(
                                    (ELEMENT_MAP[elem_point_old[0]][elem], list(elem_point_map[elem].keys())[0]))
                        else:
                            # otherwise just get new element number
                            connected_points_new.append((ELEMENT_MAP[elem_point_old[0]][0], elem_point_old[1]))

                    # connect
                    f_new.connect_element_points(connected_points_new)

                    ignore += connected_points_old

        # ~ pdb.set_trace()
        # now deal with hanging points. 
        #### Assume there are no points hanging onto subfields ####
        for hp_i in self.mesh.hanging_points:
            host_elem = element_map[sub.mesh.hanging_points[hp_i].host_element]
            host_elem_coords = self.mesh.hanging_points[hp_i].element_coordinates

            connected_points_new = []
            connected_points_old = self.mesh.connectivity[(-1, hp_i)]
            for elem_point_old in connected_points_old:
                connected_points_new.append((element_map[elem_point_old[0]], elem_point_old[1]))

            f_new.connect_to_hanging_point(host_elem, host_elem_coords, connected_points_new)

        # ==============================================================#
        # Map parameters for f_new
        f_new.map_parameters()

        if f_new.get_number_of_ensemble_points() != self.get_number_of_ensemble_points():
            raise ValueError('number of ensemble points not equal in current and flattened mesh: ' + str(
                self.get_number_of_ensemble_points()) + ' ' + str(f_new.get_number_of_ensemble_points()))

        # set custom ensemble point mapping to get same order as in self
        # for each new ensemble point, get new element points, to get old element
        # points, to get old ensemble point and set { new: old }
        r_EPM = reverse_dict(ELEMENT_POINT_MAP)
        custom_map = {}
        for ens_i_new in list(f_new.mapper._ensemble_to_element_map.keys()):
            # get the 1st element point associated with ensemble point
            ens_i_map = f_new.mapper._ensemble_to_element_map[ens_i_new]
            e_i_new = (list(ens_i_map.keys())[0], list(ens_i_map[list(ens_i_map.keys())[0]].keys())[0])
            # get corresponding old e_i
            e_i_old = r_EPM[e_i_new][0]
            # find old ensemble point
            ens_i_old = self.mapper._element_to_ensemble_map[e_i_old[0]][e_i_old[1]][0][0]
            if self.mapper.has_custom_map:
                ens_i_old = self.mapper._custom_ensemble_order[ens_i_old]
            custom_map[ens_i_new] = ens_i_old

        f_new.mapper.set_custom_ensemble_ordering(custom_map)

        if self.parameters is not None:
            # new_params = numpy.zeros((f_new.mapper.get_number_of_ensemble_points(), self.parameters.shape[1]), dtype=float)
            new_params = numpy.zeros_like(self.parameters, dtype=float)
            for i_new, i_old in list(custom_map.items()):
                new_params[i_new] = self.parameters[i_old]

            f_new.set_parameters(new_params)

        return f_new, ELEMENT_MAP, ELEMENT_POINT_MAP

    # ==================================================================#
    def set_parameters(self, parameters):
        """ Set the parameters for the field. Not sure if this is
        necessary
        """

        # check parameter length
        if len(parameters) != self.mapper.get_number_of_ensemble_points():
            raise ValueError(
                'wrong number of parameter sets, there are {} node, given {} set of nodal parameters'.format(
                    self.mapper.get_number_of_ensemble_points(),
                    len(parameters),
                )
            )
        else:
            self.parameters = parameters
            self.element_param_cache.clear()
            return 1

    # ==================================================================#
    def evaluate_field_at_element_node(self, element, point):
        # gets the values at element int:element point int:point
        # by querying mapper

        element_parameters = self.mapper.get_element_parameters(element, self.parameters)

        if element_parameters is None:
            log.debug(
                'ERROR: ensemble_field_function.evaluate_field_at_element_point: unable to get element parameters from mapper')
            return None
        else:
            try:
                point_values = element_parameters[point]
            except IndexError:
                log.debug('ERROR: ensemble_field_function.evaluate_field_at_element_point: invalid point in element')
                return None
            else:
                return point_values[0]

    # ==================================================================#
    def _get_element_parameters(self, e_i):

        # first look in cache
        try:
            params = self.element_param_cache[e_i]
        except KeyError:
            params = self.mapper.get_element_parameters(e_i, self.parameters)
            self.element_param_cache[e_i] = params

        return params

    # ==================================================================#
    def evaluate_field_in_mesh(self, density, parameters=None, derivs=None, unpack=True, subUnpack=True):
        """ Evaluates the field over the whole mesh, i.e. in all 
        elements. Returns a list of field values.
        
        density is the sampling density per dimension in each element.
        parameters is a list of field parameters, with a list for each 
        ensemble point. If parameters is not passed, will look to see if
        parameters have been set. If not, raises error.
        """
        if parameters is not None:
            self.set_parameters(parameters)

        # check that there are parameters
        if self.parameters is None:
            raise RuntimeError('no parameters passed or set')

        # check for right number of density values
        if isinstance(density, int):
            density = [density] * self.dimensions
        elif len(density) != self.dimensions:
            raise ValueError(
                'ERROR: evaluate_element: needed ' + str(self.dimensions) + ' density values. Got ' + str(len(density)))

        field_values = []
        field_derivatives = []
        element_field_values = None
        basis_values = {}
        basis_derivatives = {}
        eval_grid = {}

        for element_number in numpy.sort(list(self.mesh.elements.keys())):
            # get element
            element = self.mesh.elements[element_number]
            # for each element, get element_parameters
            element_parameters = self._get_element_parameters(element_number)

            if element.is_element:
                # if element is local, evaluate using evaluate_field_in_element
                # assume all local true elements are of the same type (since they use the same basis)
                # so generate grid points once, and calculate basis weights once
                e_type = element.type
                evaluator = self.evaluators[e_type]
                try:
                    # ~ element_field_values = numpy.dot( basis_values[e_type], element_parameters )
                    element_field_values = evaluator(basis_values[e_type], element_parameters)

                    if derivs:
                        # ~ element_field_derivatives = numpy.array( [ numpy.dot( b.T, element_parameters ) for b in basis_derivatives[e_type] ] )
                        element_field_derivatives = numpy.array(
                            [evaluator(b, element_parameters) for b in basis_derivatives[e_type]])
                    # ~ pdb.set_trace()
                except KeyError:
                    eval_grid[e_type] = element.generate_eval_grid(density)
                    # ~ pdb.set_trace()
                    basis_values[e_type] = self.basis[e_type].eval(eval_grid[e_type].T)
                    element_field_values = evaluator(basis_values[e_type], element_parameters)

                    if derivs:
                        if derivs == -1:
                            # ~ pdb.set_trace()
                            basis_derivatives[e_type] = self.basis[element.type].eval_derivatives(eval_grid[e_type].T,
                                                                                                  None)
                        elif type(derivs) == tuple:
                            basis_derivatives[e_type] = self.basis[element.type].eval_derivatives(eval_grid[e_type].T,
                                                                                                  derivs)
                        else:
                            raise NotImplementedError('derivs must be a tuple or -1')

                        element_field_derivatives = numpy.array(
                            [evaluator(b, element_parameters) for b in basis_derivatives[e_type]])

            else:
                # else element is a sub mesh, evaluate the corresponding subfield
                if derivs:
                    element_field_values, element_field_derivatives = self.subfields[
                        element_number].evaluate_field_in_mesh(density, element_parameters, derivs=derivs)
                else:
                    element_field_values = self.subfields[element_number].evaluate_field_in_mesh(density,
                                                                                                 element_parameters,
                                                                                                 derivs=derivs,
                                                                                                 unpack=subUnpack)

            field_values.append(element_field_values)
            if derivs:
                field_derivatives.append(element_field_derivatives)

        if derivs:
            if unpack:
                return numpy.hstack(field_values), numpy.hstack(field_derivatives)
            else:
                return field_values, field_derivatives
        else:
            if unpack:
                return numpy.hstack(field_values)
            else:
                return field_values

    # ==================================================================#
    def evaluate_derivatives_in_mesh(self, density, parameters=None, derivs=None, unpack=True):
        """ Evaluates the field derivatives over the whole mesh, i.e. 
        in all elements. Returns a list of field derivatives.
        
        density is the sampling density per dimension in each element.
        parameters is a list of field parameters, with a list for each 
        ensemble point. If parameters is not passed, will look to see if
        parameters have been set. If not, raises error.
        """
        if parameters is not None:
            self.set_parameters(parameters)

        # check that there are parameters
        if self.parameters is None:
            raise RuntimeError('no parameters passed or set')

        # check for right number of density values
        if isinstance(density, int):
            density = [density] * self.dimensions
        elif len(density) != self.dimensions:
            raise ValueError(
                'ERROR: evaluate_element: needed ' + str(self.dimensions) + ' density values. Got ' + str(len(density)))

        field_derivatives = []
        basis_derivatives = {}
        eval_grid = {}

        for element_number in numpy.sort(list(self.mesh.elements.keys())):
            # get element
            element = self.mesh.elements[element_number]
            # for each element, get element_parameters
            element_parameters = self._get_element_parameters(element_number)

            if element.is_element:
                # if element is local, evaluate using evaluate_field_in_element
                # assume all local true elements are of the same type (since they use the same basis)
                # so generate grid points once, and calculate basis weights once
                e_type = element.type
                evaluator = self.evaluators[e_type]
                try:
                    if derivs == -1:
                        element_field_derivatives = numpy.array(
                            [evaluator(b, element_parameters) for b in basis_derivatives[e_type]])
                    else:
                        element_field_derivatives = evaluator(basis_derivatives[e_type], element_parameters)
                    # ~ pdb.set_trace()
                except KeyError:
                    eval_grid[e_type] = element.generate_eval_grid(density)

                    if derivs == -1:
                        # ~ pdb.set_trace()
                        basis_derivatives[e_type] = self.basis[element.type].eval_derivatives(eval_grid[e_type].T, None)
                        element_field_derivatives = numpy.array(
                            [evaluator(b, element_parameters) for b in basis_derivatives[e_type]])
                    else:
                        basis_derivatives[e_type] = self.basis[element.type].eval_derivatives(eval_grid[e_type].T,
                                                                                              derivs)
                        element_field_derivatives = evaluator(basis_derivatives[e_type], element_parameters)

            else:
                # else element is a sub mesh, evaluate the corresponding subfield
                element_field_derivatives = self.subfields[element_number].evaluate_derivatives_in_mesh(density,
                                                                                                        element_parameters,
                                                                                                        derivs=derivs)

            field_derivatives.append(element_field_derivatives)

        # ~ pdb.set_trace()
        if unpack:
            return numpy.hstack(field_derivatives)
        else:
            return field_derivatives

    # ==================================================================#
    def evaluate_field_in_element(self, element_number, density, parameters=None, derivs=None, unpack=True):
        """ evaluates the field in an element defined by element_number
        """
        if parameters is not None:
            self.set_parameters(parameters)

        # check that there are parameters
        if self.parameters is None:
            raise RuntimeError('no parameters passed or set')

        # check for right number of density values
        if isinstance(density, int):
            density = [density] * self.dimensions
        elif len(density) != self.dimensions:
            log.debug(
                'ERROR: evaluate_field_in_element: needed %s density values. Got %s',
                self.dimensions,
                len(density))
            return 0

        try:
            # get element
            element = self.mesh.elements[element_number]
        except KeyError:
            raise ValueError('invalid element number ' + str(element_number))
        else:
            # get element_parameters
            element_parameters = self._get_element_parameters(element_number)

            if element.is_element:
                # if element is local, evaluate using evaluate_field_in_element
                eval_grid = element.generate_eval_grid(density)
                # ~ element_field_values = self._evaluate_element( element.type, element_parameters, eval_grid )
                element_field_values = self.evaluators[element.type](self.basis[element.type].eval(eval_grid.T),
                                                                     element_parameters)

                if derivs:
                    if derivs == -1:
                        element_deriv_values = self._evaluate_element_derivatives(element.type, element_parameters,
                                                                                  eval_grid, None)
                    elif isinstance(derivs, (list, tuple)):
                        element_deriv_values = self._evaluate_element_derivatives(element.type, element_parameters,
                                                                                  eval_grid, derivs)
                    else:
                        raise NotImplementedError('unrecognised derivative input')

                    return element_field_values, element_deriv_values
                else:
                    return element_field_values
            else:
                # else element is a sub mesh, evaluate the corresponding subfield
                return self.subfields[element_number].evaluate_field_in_mesh(density, element_parameters, derivs=derivs,
                                                                             unpack=unpack)

    # ==================================================================#
    def evaluate_field_at_element_point(self, element_number, xi, parameters=None, derivs=None):
        """ evaluate the specified element at a list of xi locations.
        If derivs = a list of derivative tuples, evaluates those derivatives,
        if derivs = -1, evaluates all derivatives
        derivative tuples: ( order, direction ), direction=3 for cross

        returns
        -------
        element_field_values : field values at the specified xi positions
        element_deriv_values : field derivatives at the specified xi positions
        """
        if parameters is not None:
            self.set_parameters(parameters)

        # check that there are parameters
        if self.parameters is None:
            raise RuntimeError('no parameters passed or set')

        try:
            # get element
            element = self.mesh.elements[element_number]
        except KeyError:
            raise ValueError('invalid element number')
        else:
            # get element_parameters
            element_parameters = self._get_element_parameters(element_number)

            if element.is_element:
                # if element is local, evaluate using evaluate_field_in_element
                # ~ element_field_values = self._evaluate_element( element.type, element_parameters, xi )
                element_field_values = self.evaluators[element.type](
                    self.basis[element.type].eval(xi.T), element_parameters
                )
                if derivs is not None:
                    if derivs == -1:
                        element_deriv_basis_values = self.basis[element.type].eval_derivatives(xi.T, None)
                        element_deriv_values = numpy.array([
                            self.evaluators[element.type](b, element_parameters) for b in element_deriv_basis_values
                        ])
                        # ~ element_deriv_values = self._evaluate_element_derivatives( element.type, element_parameters, numpy.array(xi), None )
                    else:
                        element_deriv_values = self.evaluators[element.type](
                            self.basis[element.type].eval_derivatives(xi.T, derivs), element_parameters
                        )
                        # ~ element_deriv_values = self._evaluate_element_derivatives( element.type, element_parameters, numpy.array(xi), derivs )

                    return element_field_values, element_deriv_values

                else:
                    return element_field_values
            else:
                # ~ raise RuntimeError, 'element is not true element'
                element_field_values = numpy.hstack([
                    self.subfields[element_number].evaluate_field_at_element_point(
                        elemNumber, xi[elemNumber], element_parameters, derivs
                    ) for elemNumber in list(xi.keys())
                ])

                return element_field_values

    # ==================================================================#
    def _evaluate_element(self, element_type, element_parameters, XI):
        """ evaluates element using element_parameters as a list of xi
        positions XI
        """
        if self.debug:
            log.debug('element_parameters:', element_parameters)
            log.debug('XIs:', XI)

        ## cubic hermite simplex testing ##
        # ~ if element_type == 'tri3e':
        # ~ self.basis[element_type].set_element( element )
        ## cubic hermite simplex testing ##

        ## bezier testing ##
        # ~ if self.basis[element_type].type == 'simplex_bezier_2D':
        # ~ self.basis[element_type].set_element( element )
        # ~ #calculates centre point parameters based on other control points (Farin, 1983)
        # ~ b111 = -1.0/6.0 * ( element_parameters[0] + element_parameters[3] + element_parameters[6] ) + 0.25 * (element_parameters[1]+element_parameters[2]+element_parameters[4]+element_parameters[5]+element_parameters[7]+element_parameters[8])
        # ~ element_parameters = numpy.hstack( [element_parameters, b111] )
        ## bezier testing ##

        basis_values = self.basis[element_type].eval(XI.T).T
        eval_values = numpy.dot(basis_values, element_parameters)

        return eval_values

    # ==================================================================#
    def _evaluate_element_derivatives(self, element_type, element_parameters, XI, derivatives):

        # a particular derivative
        if derivatives:
            basis_values = numpy.array(
                [self.basis[element_type].eval_derivatives(numpy.array(XI).T, d).T for d in derivatives])
            eval_values = numpy.dot(basis_values, element_parameters)
        # all derivatives
        else:
            basis_values = self.basis[element_type].eval_derivatives(numpy.array(XI).T, None)
            eval_values = numpy.array([numpy.dot(b.T, element_parameters) for b in basis_values])

        return eval_values

    # ==================================================================#
    def get_element_point_evaluator(self, element_number):
        """ returns an object for evaluating points within element
        element_number
        """

        # check for valid element_number
        try:
            element = self.mesh.elements[element_number]
        except KeyError:
            log.debug('ERROR ensemble_field_function.get_element_point_evaluator: invalid element')
            return None
        else:
            # if element is an element....
            if element.is_element:
                # get element parameters
                parameters = self.mapper.get_element_parameters(element_number)
                # create evaluation object
                evaluator = element_point_evaluator(element, parameters, self.basis[element.type])
            # else element is a subfield
            else:
                # ~ evaluator = self.submesh_map[ element ].get_element_point_evaluator( element_number )
                log.debug('ERROR: ensemble_field_function.get_element_point_evaluator: element', element_number,
                      'is a subfield')
                return None

            return evaluator

    # ~ #==================================================================#
    # ~ def _equi_evalx( self, y, b ):
    # ~ """ calculates the x coord on the edge of an equilateral triangle
    # ~ give y
    # ~ """
    # ~ a = numpy.sqrt(3.0)
    # ~ return ((y - b)/a)
    # ~
    # ~ #==================================================================#
    # ~ def _equi_gen_points( self, d, int ):
    # ~ """ generate element points for a equilateral triangle.
    # ~ Returns array( [[xcoords], [ycoords]] ).
    # ~ d: list of evaluation density in each direction
    # ~ int: interior bounds of the element ( (xmin, xmax), (ymin, ymax) )
    # ~ """
    # ~
    # ~ e_points = numpy.array( [[],[]] )
    # ~ y_divs = numpy.linspace( int[1][0], int[1][1], d[1] )
    # ~ xn = d[0]
    # ~ for y_row in y_divs:
    # ~ x_edge = self._equi_evalx( y_row, int[1][1] )
    # ~ x = numpy.linspace( x_edge,-x_edge, xn )
    # ~ e_points = numpy.hstack( [e_points, numpy.array( [ x, [y_row]*xn] ) ] )
    # ~ xn -= 1
    # ~ return e_points

    # ==================================================================#
    # ~ def evaluate_field_at_element_point( self, element_number, element_coordinates, element_parameters ):
    # ~ # element_number is an integer
    # ~ # element_coordinates should be a list or tuple of xi values
    # ~ #    e.g. [xi1, xi2]
    # ~ # element_parameters should be a list of parameters for each
    # ~ #    element point
    # ~ #    e.g. [ [p1 parameters], [p2 parameters]...]


# ~
# ~ element = self.mesh.elements[ element_number ]
# ~
# ~
# ~
# ~ # check coordinates are in element domain
# ~
# ~ if element.is_interior( element_coordinates ):
# ~
# ~ u_vector = []
# ~
# ~ if ( self.basis_type == 'cubic_hermite' ) and self.dimensions == 2:
# ~ # special mapping of element parameters to the tensor product
# ~ # of basis functions for 2D cubic hermite
# ~ # should be improved!
# ~ # goes [1 2 5 6 3 4 7 8 9 10 13 14 11 12 15 16]
# ~ u_vector = []
# ~ for i in range( 0, self.dimensions+1, 2 ):
# ~ for j in range( 0, self.dimensions+1, 2 ):
# ~ u_vector += list( element_point_values[i:i+self.dimensions, j:j+self.dimensions].ravel() )
# ~
# ~ else:
# ~ u_vector = element_point_values.ravel()
# ~
# ~ if self.debug:
# ~ print 'element_point_values:', element_point_values
# ~ print 'u_vector:', u_vector
# ~
# ~ #
# ~ for i in range( len( element_coordinates ) ):
# ~ # evaluate basis function at element_coordinates
# ~ phi = self.basis.eval( element_coordinates[i] )
# ~
# ~ if self.debug:
# ~ print 'phi:', phi
# ~
# ~ if not isinstance( phi, int ):
# ~ # inner-product with basis
# ~ field_value = numpy.dot( phi, u_vector )
# ~ if self.debug:
# ~ print 'field_values:', field_values
# ~ else:
# ~ print 'ERROR: basis.eval'
# ~ return 0
# ~
# ~ field_values.append( field_value )
# ~
# ~
# ~ return field_values

# ======================================================================#

# derivative mapping from cartesian to triangle edge directions
cos45 = numpy.sqrt(2.0) / 2.0
cart_map = numpy.array([[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], \
                        [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], \
                        [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], \
                        [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0], \
                        [0.0, 0.0, 0.0, 0.0, -cos45, cos45, 0.0, 0.0, 0.0], \
                        [0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0], \
                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0], \
                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0], \
                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, cos45, -cos45]])


def map_derivs(p, wdx1, wdx2):
    m = p.copy()
    m[0, 1] = p[0, 1] * wdx1[0] + p[0, 2] * wdx2[0]  # du0/dL0
    m[0, 2] = p[0, 1] * wdx1[2] + p[0, 2] * wdx2[2]  # du0/dL2

    m[1, 1] = p[1, 1] * wdx1[1] + p[1, 2] * wdx2[1]  # du1/dL1
    m[1, 2] = p[1, 1] * wdx1[0] + p[1, 2] * wdx2[0]  # du1/dL0

    m[2, 1] = p[2, 1] * wdx1[2] + p[2, 2] * wdx2[2]  # du2/dL2
    m[2, 2] = p[2, 1] * wdx1[1] + p[2, 2] * wdx2[1]  # du2/dL1

    return m


class element_point_evaluator:
    """ Class for objects evaluating the field within an element
    """

    def __init__(self, element, parameters, basis_function):
        self.element = element
        self.parameters = parameters
        self.dimensions = self.element.dimensions
        self.basis = basis_function

    def evaluate_at_element_point(self, element_coordinates):
        """
        Evaluates the objects element at the given element coordinates.

        element_coordinates is a list.
        """

        # evaluate basis function
        basis_coeff = self.basis.eval(element_coordinates)
        if basis_coeff is not None:
            try:
                field_value = numpy.dot(basis_coeff, self.parameters)
                return field_value
            except ValueError:
                log.warning('dot product failed, basis_coeff and parameters misaligned')
                return None

        else:
            log.warning('unable to evaluate basis function at element coordinates: %s', element_coordinates)
            return None


def reverse_dict(d):
    d_r = {}
    for k in list(d.keys()):
        if type(d[k]) == list:
            for i in d[k]:
                try:
                    d_r[i].append(k)
                except KeyError:
                    d_r[i] = [k]
        else:
            try:
                d_r[i].append(k)
            except KeyError:
                d_r[i] = [k]

    return d_r


# =============================================================================#
def load_eff_shelve(filename, E, meshfn=None, filedir=None, force=False):
    """
    If force==True, will ignore if mesh cannot be loaded
    """

    if filedir is not None:
        filename = os.path.join(filedir, filename)
    else:
        filedir = ''

    try:
        try:
            S = shelve.open(filename, 'r')
        except ImportError:
            import bsddb3
            _db = bsddb3.hashopen(filename)
            S = shelve.Shelf(_db)
    except:
        raise IOError(filename + ' not found')
    else:
        E.name = S['name']
        E.dimensions = S['dimensions']
        E.debug = S['debug']

        try:
            E.set_basis(S['basis'])
        except KeyError:
            E.set_basis({'tri10': S['basis_type']})  #### HACK

        if force:
            if meshfn:
                E.mesh = mesh.load_mesh(meshfn, path=filedir)
            else:
                if S['mesh'][:len(filedir)] == filedir:
                    try:
                        E.mesh = mesh.load_mesh(S['mesh'][len(filedir):], path=filedir)
                    except IOError:
                        pass
                elif S.get('mesh'):
                    try:
                        E.mesh = mesh.load_mesh(S['mesh'], path=filedir)
                    except IOError:
                        pass
                else:
                    log.debug('WARNING: no mesh loaded')

            for s in S['subfields']:
                if s[1][:len(filedir)] == filedir:
                    try:
                        E.subfields[s[0]] = load_eff_shelve(s[1][len(filedir):], path=filedir)
                    except IOError:
                        pass
                else:
                    try:
                        E.subfields[s[0]] = load_eff_shelve(s[1], path=filedir)
                    except IOError:
                        pass
        else:
            if meshfn:
                E.mesh = mesh.load_mesh(meshfn, path=filedir)
            else:
                if S['mesh'][:len(filedir)] == filedir:
                    E.mesh = mesh.load_mesh(S['mesh'][len(filedir):], path=filedir)
                elif S.get('mesh'):
                    E.mesh = mesh.load_mesh(S['mesh'], path=filedir)
                else:
                    log.debug('WARNING: no mesh loaded')
            for s in S['subfields']:
                sub_eff = ensemble_field_function(None, None)
                if s[1][:len(filedir)] == filedir:
                    E.subfields[s[0]] = load_eff_shelve(s[1][len(filedir):], sub_eff, filedir=filedir)
                else:
                    E.subfields[s[0]] = load_eff_shelve(s[1], sub_eff, filedir=filedir)

        E.subfield_counter = S['subfield_counter']
        if E.mesh is not None:
            E.map_parameters()
            if S['custom_map']:
                E.mapper.set_custom_ensemble_ordering(S['custom_map'])
        else:
            if S['custom_map']:
                E.mapper._custom_ensemble_order = S['custom_map']

        return E


class EFFJSONWriter(object):

    def __init__(self, eff):
        """
        Writer class for serialising an ensemble file function to a
        JSON format file or string.
        """
        self.eff = eff
        self._file_dir = ''

    def write(self, filename, meshfn=None, filedir=None, subfieldfns=None):
        """
        Write to file a json format file of the eff properties.
        """
        if filedir is not None:
            self._file_dir = filedir

        d = self.serialise(meshfn, subfieldfns)
        with open(filename, 'w') as f:
            json.dump(d, f, indent=4, sort_keys=True)

        return filename

    def serialise(self, meshfn, subfieldfns):
        """
        Return a json-compatible dict of the eff properties.
        """
        d = {}
        self._serialise_meta(d)
        self._serialise_basis(d)
        self._serialise_mesh(d, meshfn)
        self._serialise_subfields(d, subfieldfns)
        self._serialise_custom_map(d)
        return d

    def _serialise_meta(self, eff_dict):
        eff_dict['name'] = self.eff.name
        eff_dict['dimensions'] = self.eff.dimensions
        eff_dict['subfield_counter'] = self.eff.subfield_counter

    def _serialise_basis(self, eff_dict):
        eff_dict['basis'] = dict([[i[0], i[1].type] for i in list(self.eff.basis.items())])

    def _serialise_mesh(self, eff_dict, meshfn):
        if meshfn is not None:
            self.eff.mesh.save_mesh(
                meshfn, path=self._file_dir
            )
            eff_dict['mesh'] = os.path.split(meshfn)[1]
        else:
            pass
            # eff_dict['mesh'] = self.eff.mesh.save_mesh(
            #     self.eff.mesh.name, path=self._file_dir
            #     )

    def _serialise_subfields(self, eff_dict, subfieldfns):
        subf_d = {}
        subf_nums = sorted(list(self.eff.subfields.keys()))
        for sn in subf_nums:
            s = self.subfields[en]
            if subfieldfns:
                subf_d[sn] = save_eff_json(
                    subfieldfns[sn],
                    s,
                    self._file_dir,
                )
            else:
                subf_d[sn] = save_eff_json(
                    s.name,
                    s,
                    self._file_dir,
                )

        eff_dict['subfields'] = subf_d

    def _serialise_custom_map(self, eff_dict):
        eff_dict['custom_map'] = self.eff.mapper._custom_ensemble_order


class EFFJSONReader(object):

    def __init__(self, eff):
        """
        Reader class for setting an ensemble field function with properties
        from a JSON format file or string.
        """
        self.eff = eff
        self._file_dir = ''
        self.force = False

    def read(self, filename, meshfn=None, filedir=None, force=False):
        """
        Load eff properties from a eff file. If force==True, will ignore if
        mesh cannot be loaded.
        """
        self.force = force
        if filedir is not None:
            filename = os.path.join(filedir, filename)
            self._file_dir = filedir
        # else:
        #     self._file_dir = os.path.split(filename)[0]
        with open(filename, 'r') as f:
            eff_dict = json.load(f)

        self.deserialise(eff_dict, meshfn)

    def deserialise(self, jsonstr, meshfn):
        """
        Set self.eff with properties from the json string.
        file_dir is the directory of any embedded eff files.
        """
        self._parse_meta(jsonstr)
        self._parse_basis(jsonstr)
        self._parse_mesh(jsonstr, meshfn)
        self._parse_subfields(jsonstr)
        if self.eff.mesh is not None:
            self.eff.map_parameters()
        self._parse_custom_map(jsonstr)

    def _parse_meta(self, eff_dict):
        self.eff.name = eff_dict['name']
        self.eff.dimensions = int(eff_dict['dimensions'])
        self.eff.subfield_counter = int(eff_dict['subfield_counter'])

    def _parse_basis(self, eff_dict):
        basis_dict = dict([(str(i[0]), i[1]) for i in eff_dict['basis'].items()])
        self.eff.set_basis(basis_dict)

    def _parse_mesh(self, eff_dict, meshfn):
        if self.force:
            if meshfn is not None:
                self.eff.mesh = mesh.load_mesh(meshfn, path=self._file_dir)
            elif eff_dict.get('mesh'):
                try:
                    self.eff.mesh = mesh.load_mesh(eff_dict['mesh'], path=self._file_dir)
                except IOError:
                    log.debug('Mesh failed to load, ignoring.')
                    pass
        else:
            if meshfn is not None:
                self.eff.mesh = mesh.load_mesh(meshfn, path=self._file_dir)
            elif eff_dict.get('mesh'):
                self.eff.mesh = mesh.load_mesh(eff_dict['mesh'], path=self._file_dir)
            else:
                log.debug('WARNING: no mesh loaded')

    def _parse_subfields(self, eff_dict):
        for _sn in eff_dict['subfields']:
            sn = int(sn)
            self.eff.subfields[sn] = load_ensemble(
                eff_dict['subfields'][_sn], path=self._file_dir
            )

    def _parse_custom_map(self, eff_dict):
        if self.eff.mesh is not None:
            if eff_dict['custom_map'] is not None:
                cust_map_dict = _intdict(eff_dict['custom_map'])
                self.eff.mapper.set_custom_ensemble_ordering(cust_map_dict)
        else:
            if eff_dict['custom_map'] is not None:
                cust_map_dict = _intdict(eff_dict['custom_map'])
                self.eff.mapper._custom_ensemble_order = cust_map_dict


def load_eff_json(filename, eff, meshfn=None, filedir=None, force=False):
    reader = EFFJSONReader(eff)
    reader.read(filename, meshfn, filedir, force=force)
    return eff


def save_eff_json(filename, eff, meshfn=None, filedir=None, subfieldfns=None):
    writer = EFFJSONWriter(eff)
    writer.write(filename, meshfn, filedir, subfieldfns)


def load_ensemble(filename, meshFilename=None, path=None, force=False):
    mesh = ensemble_field_function(None, None)
    with open(filename, 'r') as f:
        head = f.read(1)
        if head == '{':
            load_eff_json(filename, mesh, meshfn=meshFilename, filedir=path, force=force)
        else:
            load_eff_shelve(filename, mesh, meshfn=meshFilename, filedir=path, force=force)

    return mesh


def _intdict(d):
    """
    convert keys and values of dictionary into ints
    """
    newd = {}
    for k, v in d.items():
        newd[int(k)] = int(v)

    return newd
