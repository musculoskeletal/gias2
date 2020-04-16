"""
FILE: mesh.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION:
class for holding information about the mesh topology.

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import json
import logging
import os
import shelve

from numpy import linspace, sort

from gias2.fieldwork.field.topology import element_types

log = logging.getLogger(__name__)


# =============================================================================#
class mesh_ensemble:
    """ Class for holding information about field topology
    """

    def __init__(self, name, dimensions, debug=0):
        self.name = name
        self.dimensions = dimensions  # number of dimensions of the mesh topology
        self.element_counter = 0  # counter of number of elements
        self.hanging_point_counter = 0  # counter for number of hanging points
        self.number_of_points = None  # number of unique ensemble points

        self.elements = {}  # { element_number: element, element_number: sub mesh_ensemble }

        self.connectivity = {}  # { ( element_num, point_num, ?? ) : [(element_num, point_num, ??), ...]] } !! needs to incorporate topological information: xi orientation and element CS density
        self.hanging_points = {}  # { hanging_point number: hanging_point, ... }
        self.element_points = {}  # { element_number: [ element_points ] } maps element to its element points

        self.is_element = False
        self.submesh_counter = 0
        self.debug = debug

    # ==================================================================#
    def save_mesh(self, filename=None, path=None, submeshFilenames=None):
        if not filename:
            filename = self.name

        if os.path.splitext(filename)[1].lower() != '.mesh':
            filename = filename + '.mesh'

        return save_mesh_json(
            filename, self, filedir=path, subfieldfns=submeshFilenames
        )

    def save_mesh_shelve(self, filename=None, path='', submeshFilenames=None):
        raise DeprecationWarning('Shelve serialisation will be deprecated')

        if path:
            if path[-1] != '/':
                path += '/'

        if not filename:
            filename = self.name

        try:
            S = shelve.open(path + filename + '.mesh')
        except ImportError:
            import bsddb3
            _db = bsddb3.hashopen(path + filename + '.mesh')
            S = shelve.Shelf(_db)

        S['name'] = self.name
        S['dimensions'] = self.dimensions
        S['element_counter'] = self.element_counter
        S['hanging_point_counter'] = self.hanging_point_counter
        S['number_of_points'] = self.number_of_points

        # [(element_number, 0, element_type), (element_number, 1, mesh filename), ...]
        element_list = []
        for e in list(self.elements.items()):
            if e[1].is_element:
                element_list.append((e[0], 0, e[1].type))
            else:
                if submeshFilenames:
                    element_list.append((e[0], 1, os.path.split(path + submeshFilenames[e[0]])[1]))
                else:
                    element_list.append((e[0], 1, os.path.split(e[1].save_mesh(path=path))[1]))

        S['elements'] = element_list

        # (hanging_point_number, host_element, element_coord )
        S['hanging_points'] = [(h[0], h[1].host_element, h[1].element_coordinates) for h in
                               list(self.hanging_points.items())]

        S['connectivity'] = self.connectivity
        S['element_points'] = self.element_points

        S['is_element'] = self.is_element
        S['submesh_counter'] = self.submesh_counter
        S['debug'] = self.debug

        S.close()

        return os.path.relpath(path + filename + '.mesh')

    def load_mesh(self, filename):
        S = shelve.open(filename)

        self.name = S['name'] = self.name
        S['dimensions'] = self.dimensions
        S['element_counter'] = self.element_counter
        S['hanging_point_counter'] = self.hanging_point_counter
        S['number_of_points'] = self.number_of_points

        S['element_list'] = [e.type for e in self.elements]

        S['connectivity'] = self.connectivity
        S['element_points'] = self.element_points

        S['is_element'] = self.is_element
        S['submesh_counter'] = self.submesh_counter
        S['debug'] = self.debug

        S.close()

    # ==================================================================#
    def add_element(self, element):
        """ adds an element to the mesh. Element can either be an element
        object or another mesh_ensemble. Element is added to the elements
        dict as a value with an element number. get_number_of_points is called
        to get the number of entries to add to connectivity and the 
        element_points map.
        
        returns element number for convenience.
        """

        self.elements[self.element_counter] = element
        if not element.is_element:
            self.submesh_counter += 1

        self.element_points[self.element_counter] = []
        for p in range(0, element.get_number_of_ensemble_points()):
            self.connectivity.update({(self.element_counter, p): []})
            self.element_points[self.element_counter].append((self.element_counter, p))

        self.element_counter += 1

        return (self.element_counter - 1)

    # ==================================================================#
    def remove_element(self, element_number):
        """ removes element element_number from the elements dict, and all
        entries the element from connectivity. Decreases element
        counter by 1.
        """

        # check for valid element number
        if element_number not in list(self.elements.keys()):
            raise ValueError('element' + str(element_number) + 'does not exist')

        # removed from elements dict
        # if not self.elements[element_number].is_element:
        #   self.submesh_counter -= 1
        del self.elements[element_number]

        # removed element points from connectivity
        element_points = self.element_points[element_number]
        # for each element point
        for point in element_points:
            # get points its connected to
            connected_points = self.connectivity[point]
            # remove this element point for its connected points
            for connected_point in connected_points:
                self.connectivity[connected_point].remove(point)
            # removed current element point entry
            del self.connectivity[point]

        # removed element points map entry
        del self.element_points[element_number]

        # decrease element counter
        ### shouldnt decrease since this is used when adding new elements - if a non-last element is removed, new element will overwrite last element ###
        # self.element_counter -= 1

        return

    # ==================================================================#
    def get_number_of_elements(self):
        return len(list(self.elements.keys()))

    # ==================================================================#
    def is_flat(self):

        for e in list(self.elements.values()):
            if not e.is_element:
                return False

        return True

    # ==================================================================#
    def get_number_of_true_elements(self):
        n = 0
        for e in list(self.elements.values()):
            if e.is_element:
                n += 1
            else:
                n += e.get_number_of_true_elements()

        return n

    # ==================================================================#
    def get_true_elements(self):
        """ gets a list of pointers to all true elements in a mesh, in
        the order that they are evaluated
        """
        elem_list = []
        for element_number in sort(list(self.elements.keys())):
            # get element
            element = self.elements[element_number]
            if element.is_element:
                elem_list.append(element)
            else:
                # ~ elem_list += _get_field_true_elements( element )
                elem_list += element.get_true_elements()
        return elem_list

    # ==================================================================#
    def get_number_of_ensemble_points(self):
        """ calculates the number of unique ensemble points in the mesh
        by going through the connectivity dict.
        """

        gp = 0  # global points number
        assigned = []  # list of element points already assigned a global number

        # account for points connected to hanging nodes
        # ~ if self.hanging_points:
        # ~ for points in self.hanging_point_connectivity.values():
        # ~ assigned += points

        # loop through each element point in the connectivity dict
        element_points = list(self.connectivity.keys())
        element_points.sort()
        if self.debug:
            log.debug('sorted element points:', element_points)
        for element_point in element_points:
            # check if point has already been assigned a global number
            if (element_point not in assigned):
                # record element points assigned
                assigned += [element_point] + self.connectivity[element_point]
                if self.debug:
                    log.debug('assigned:', assigned)

                # if not a hanging point, count as an ensemble point
                if element_point[0] != -1:
                    gp += 1

        self.number_of_points = gp

        return self.number_of_points

    # ==================================================================#
    def connect_elements_by_edges(self, elem1, edge1, elem2, edge2, xi_extent=None, invert_order=False):
        """ simple way for connect elements along edges. xi_extent allows
        edge-end hanging points to be define, and intermediate hanging
        points are automatically define. 
        """

        elem1_points = self.elements[elem1].get_edge_points(edge1)
        elem2_points = self.elements[elem2].get_edge_points(edge2)
        if invert_order:
            elem2_points = elem2_points[::-1]

        if xi_extent:
            xi = linspace(xi_extent[0], xi_extent[1], len(elem2_points))
            p_i = 0
            for p in elem2_points:
                self.connect_to_hanging_point(elem1, [(elem2, p)], 'edge', (edge1, xi[p_i]))
                p_i += 1
        else:
            if len(elem1_points) == len(elem2_points):
                for i in range(len(elem1_points)):
                    self.connect([(elem1, elem1_points[i]), (elem2, elem2_points[i])])
            else:
                self.connect([(elem1, elem1_points[0]), (elem2, elem2_points[0])])
                self.connect([(elem1, elem1_points[-1]), (elem2, elem2_points[-1])])

                xi = linspace(0.0, 1.0, len(elem2_points))
                for i in range(len(xi) - 2):
                    self.connect_to_hanging_point(elem1, [(elem2, elem2_points[i + 1])], 'edge', xi[i + 1])

    # ==================================================================#
    def connect(self, common_points):
        """ define how the elements are connected
        common_points is a list of tuples of element nodes, or coordinates at the 
        same ensemble mesh position: 
        [ (element_num, node_num1, ??), (element_num2, (element_xi), ??), ...]
        
        !! should incorporate information about element coord orientation and coord density here !!
        !! should be able to connect to any position in an element, not only to other element points => hanging nodes handled naturally !!
        !! => connectivity contains extra topology info!!
        """
        if len(common_points) == 2 and common_points[0] == common_points[1]:
            return

        for point in common_points:

            if point not in list(self.connectivity.keys()):
                log.debug('ERROR: toplogy.connect: point', point, 'does not exist')
                avail_keys = list(self.connectivity.keys())
                avail_keys.sort()
                log.debug('available points to connect:\n', avail_keys)
                return None

            other = list(common_points[:])
            other.remove(point)

            # ~ # update connectivity for other points already connected to point
            # ~ for p in self.connectivity[ point ]:
            # ~ self._do_connect( p, other )
            # ~
            # ~ # update connectivity for point
            # ~ self._do_connect( point, other )

            # for each point in other, check if point is in its entry
            for p in other:
                if point not in self.connectivity[p]:
                    # if not, add point to points in the entry
                    self.connectivity[p].append(point)
                    [self.connect((point, pp)) for pp in self.connectivity[p]]

        return 1

    # ==================================================================#
    def _do_connect(self, point, other_points):
        if len(other_points) == 1 and point == other_points[0]:
            return
        else:
            for other_point in other_points:
                if other_point not in self.connectivity[point]:
                    self.connectivity[point].append(other_point)
                    if point not in self.connectivity[other_point]:
                        self.connectivity[other_point].append(point)

        return

    # ~ def connect_to_hanging_point( self, host_element, element_coordinates, connected_points ):
    # ~ """ host_element: int
    # ~ element_coordinates: [xi2, xi2,...]
    # ~ connected_points: [ (element_num, point_num), ...]
    # ~
    # ~ creates a hanging point in host_element at element_coordinates,
    # ~ and connects connected_points to the hanging point
    # ~ """
    # ~
    # ~ # create hanging point
    # ~ if self._add_hanging_point( host_element, element_coordinates ):
    # ~ # add entry to connectivity
    # ~ connected_points.append( (-1, self.hanging_point_counter - 1 ) )
    # ~ # connect
    # ~ self.connect( connected_points )
    # ~ return 1
    # ~ else:
    # ~ print 'ERROR: topology.connect_to_hanging_point: failed to add hanging_point'
    # ~ return None

    # ==================================================================#
    def connect_to_hanging_point(self, host_element, connected_points, mode, arg):
        """ host_element: int
            connected_points: [ (element_num, point_num), ...]
            mode: 'interior', or 'edge'
            arg: [xi1, xi2,...] for interior, [edge_number, xi] for edge
        
        creates a hanging point in host_element at element_coordinates,
        and connects connected_points to the hanging point
        """

        # check args
        if mode == 'interior':
            if len(arg) != self.dimensions:
                raise ValueError('incorrect number of xi coordinates for interior hanging point')
            else:
                element_coord = arg
        elif mode == 'edge':
            if len(arg) != 2:
                raise ValueError('incorrect arg length, expected 2: edge_number, edge xi')
            else:
                element_coord = self.elements[host_element].get_edge_coord(arg[0], arg[1])

        # create hanging point  
        if self._add_hanging_point(host_element, element_coord):
            # add entry to connectivity
            connected_points.append((-1, self.hanging_point_counter - 1))
            # connect
            self.connect(connected_points)
            return 1
        else:
            raise RuntimeError('failed to add hanging_point')

    # ==================================================================#
    def _add_hanging_point(self, host_element, element_coordinates):
        """ adds a hanging point to the mesh in host_element at 
        element_coordinates
        """

        # check for valid element and number of coordinates
        # will have trouble with child elements with different dimensionality
        if host_element in list(self.elements.keys()):
            if len(element_coordinates) == self.elements[host_element].dimensions:
                # instantiate hanging_point
                hp = hanging_point(host_element, element_coordinates)

                # add hanging point to hanging_points dict and connectivity under element -1
                self.hanging_points[self.hanging_point_counter] = hp
                self.connectivity[(-1, self.hanging_point_counter)] = []
                self.hanging_point_counter += 1

                return 1
            else:
                log.debug('ERROR: topology.add_hanging_point: incorrect number of element coordinates. Need',
                      self.elements[host_element].dimensions, ', got', len(element_coordinates))
                return None
        else:
            log.debug('ERROR: topology.add_hanging_point: element', host_element, 'does not exist')
            return None

    # ==================================================================#
    # ~ def _get_element_points( self, element_number ):
    # ~ # returns a list of tuples containing all elements of element
    # ~ # element_number
    # ~
    # ~ # check for valid element
    # ~ if element_number not in self.elements.keys():
    # ~ print 'ERROR: topology._get_element_points: element', element_number, 'does not exist'
    # ~ return 0
    # ~
    # ~ element_points = []
    # ~ ckeys = self.connectivity.keys()
    # ~ ckeys.sort
    # ~ ckeys_iter = iter(ckeys)
    # ~ done = 0
    # ~
    # ~ # look through sorted connectivity keys until 1st element point
    # ~ # of element
    # ~ while not done:
    # ~ key = ckeys_iter.next()
    # ~ if key[0] = element_number:
    # ~ element_points.append( key )

    # ==================================================================#
    # ~ def add_local_element( self, type ):
    # ~ """ adds a element to the mesh. Type is a string of the of the
    # ~ type of element. The element is appended to the
    # ~ self.local_elements dict with the current element_counter value.
    # ~ List of nodes in that element are appended to the connectivity
    # ~ dict. self.element_counter is increased by 1
    # ~ """
    # ~ # instantiate element
    # ~ element = element_types.create_element( type )
    # ~
    # ~ if element:
    # ~ # add to local_elements dict
    # ~ self.elements[ self.element_counter ] = element
    # ~
    # ~ # add nodes as keys to connectivity
    # ~ for p in range(0, element.number_of_points ):
    # ~ self.connectivity.update( { [element_counter, p]: []} )
    # ~
    # ~ element_number += 1
    # ~ return 1
    # ~ else:
    # ~ print 'ERROR: mesh_ensemble.add_element failed'
    # ~ return 0

    # ==================================================================#
    # ~ def add_sub_mesh( self, sub_mesh ):
    # ~ """ sub_mesh is an ensemble_mesh instance. it is given an element
    # ~ number and appended to the self.elements dict. Queries the sub_mesh's
    # ~ get_number_of_points() to get number of points exposed to the parent.
    # ~ """
    # ~
    # ~
    # ~ # add sub mesh elements/element_points to connectivity dict using
    # ~ # the connectivity in the sub_mesh
    # ~ e_start = len( self.elements )   # element number of elements in the sub mesh are padded by this number
    # ~
    # ~ sub_mesh_elements = sub_mesh.connectivity.keys()
    # ~ sub_mesh_elements.sort()
    # ~
    # ~ for key in sub_mesh_elements:
    # ~
    # ~ sub_element = key[0]
    # ~ sub_point = key [1]
    # ~
    # ~ # add new connectivity key
    # ~ self.connectivity.update( { ( sub_element + e_start, sub_point ): [] }
    # ~ # add values to new connectivity key
    # ~ for ( e, p ) in sub_mesh.connectivity[key]:
    # ~ self.connectivity[ sub_element + e_start, sub_point ].append( ( e + e_start, p ) )
    # ~
    # ~ # map element number to its sub_mesh, and sub mesh element number
    # ~ self.sub_mesh_element_map.update( { ( sub_element + e_start, sub_point ): (sub_mesh_number, sub_element) } )
    # ~
    # ~ return 1


# hanging point object. Stored information about its host_element and
# host element coordinates  
class hanging_point:
    """ Class for holding information about a hanging_point (node)
    """

    def __init__(self, host_element, element_coordinates):
        self.host_element = host_element
        self.element_coordinates = tuple(element_coordinates)

    def get_host_element(self):
        return self.host_element

    def get_element_coordinates(self):
        return self.element_coordinates


# =============================================================================#
def load_mesh_shelve(filename, mesh, filedir=None):
    if filedir is not None:
        filename = os.path.join(filedir, filename)
    else:
        filedir = ''

    if filename[-5:] != '.mesh':
        filename += '.mesh'
    try:
        try:
            S = shelve.open(filename, 'r')
        except ImportError:
            import bsddb3
            _db = bsddb3.hashopen(filename)
            S = shelve.Shelf(_db)
    except:
        raise IOError(filename + ' not found')

    mesh.name = S['name']
    mesh.dimensions = S['dimensions']
    mesh.debug = S['debug']

    # put in elements
    for e in S['elements']:
        if e[1] == 0:
            # true element
            mesh.elements[e[0]] = element_types.create_element(e[2])
        else:
            # submesh
            if e[2][:len(filedir)] == filedir:
                mesh.elements[e[0]] = load_mesh_shelve(
                    e[2][len(filedir):],
                    mesh_ensemble(None, None),
                    filedir=filedir
                )
            else:
                mesh.elements[e[0]] = load_mesh_shelve(
                    e[2],
                    mesh_ensemble(None, None),
                    filedir=filedir
                )

    # create hanging_points
    for h in S['hanging_points']:
        mesh.hanging_points[h[0]] = hanging_point(h[1], h[2])

    # set connectivity
    mesh.connectivity = S['connectivity']
    # set element_points
    mesh.element_points = S['element_points']
    # set others
    mesh.element_counter = S['element_counter']
    mesh.hanging_point_counter = S['hanging_point_counter']
    mesh.number_of_points = S['number_of_points']
    mesh.submesh_counter = S['submesh_counter']

    S.close()
    return mesh


class MeshJSONWriter(object):

    def __init__(self, mesh):
        """
        Writer class for serialising a mesh to a
        JSON format file or string.

        Hanging points are not currently handled.
        """
        self.mesh = mesh
        self._file_dir = ''

    def write(self, filename, filedir=None, submeshfns=None):
        """
        Write to file a json format file of the mesh properties.
        """
        if filedir is not None:
            self._file_dir = filedir

        d = self.serialise(submeshfns)
        with open(os.path.join(self._file_dir, filename), 'w') as f:
            json.dump(d, f, indent=4, sort_keys=True)

        return filename

    def serialise(self, submeshfns):
        """
        Return a json-compatible dict of the mesh properties.
        """
        d = {}
        self._serialise_meta(d)
        self._serialise_elements(d, submeshfns)
        self._serialise_connectivity(d)
        return d

    def _serialise_meta(self, mesh_dict):
        mesh_dict['name'] = self.mesh.name
        mesh_dict['dimensions'] = self.mesh.dimensions
        mesh_dict['number_of_points'] = self.mesh.number_of_points
        mesh_dict['element_counter'] = self.mesh.element_counter
        mesh_dict['submesh_counter'] = self.mesh.submesh_counter

    def _serialise_elements(self, mesh_dict, submeshfn):
        """
        Create dictionary of elements. Each entry is of the format
        [elem_number] : [is_submesh] [elem_name]
        where is_submesh is 1 if element is another mesh else 0.
        If is_submesh is 1, the elem_name is the is name of the mesh
        or submesh[submesh.name].
        If is_submesh is 0, the elem_name is the element type.
        """
        elems_d = {}
        elem_nums = sorted(list(self.mesh.elements.keys()))
        for en in elem_nums:
            e = self.mesh.elements[en]
            if e.is_element:
                elems_d[en] = '{} {}'.format(0, e.type)
            else:
                if submeshfn:
                    if self._file_dir is not None:
                        elems_d[en] = '{} {}'.format(
                            1,
                            save_mesh_json(
                                submeshfn[en],
                                e,
                                self._file_dir,
                            )
                        )
                    else:
                        elems_d[en] = '{} {}'.format(
                            1,
                            save_mesh_json(
                                submeshfn[en],
                                e
                            )
                        )
                else:
                    if self._file_dir is not None:
                        elems_d[en] = '{} {}'.format(
                            1,
                            save_mesh_json(
                                e.name, e, self._file_dir
                            )
                        )
                    else:
                        elems_d[en] = '{} {}'.format(
                            1, save_mesh_json(e.name, e)
                        )

        mesh_dict['elements'] = elems_d

    def _serialise_connectivity(self, mesh_dict):
        """
        Create json-compatible dict of the connectivity dict.
        Each entry is of the format
        [elem_num]_[point_num] : [elem_num]_[point_num] [elem_num]_[point_num] ...
        """

        conn_d = {}
        keys = sorted(list(self.mesh.connectivity.keys()))
        for k in keys:
            k_str = '{}_{}'.format(*k)
            v_str = ' '.join(['{}_{}'.format(v[0], v[1]) for v in self.mesh.connectivity[k]])
            conn_d[k_str] = v_str

        mesh_dict['connectivity'] = conn_d


class MeshJSONReader(object):

    def __init__(self, mesh):
        """
        Reader class for setting a mesh with properties from a
        JSON format file or string.

        Hanging points are not currently handled.
        """
        self.mesh = mesh
        self._file_dir = ''

    def read(self, filename, filedir=None):
        """
        Load mesh properties from a mesh file.
        """
        if filedir is not None:
            filename = os.path.join(filedir, filename)
            self._file_dir = filedir
        # else:
        #     self._file_dir = os.path.split(filename)[0]
        with open(filename, 'r') as f:
            mesh_dict = json.load(f)

        self.deserialise(mesh_dict)

    def deserialise(self, jsonstr):
        """
        Set self.mesh with properties from the json string.
        file_dir is the directory of any embedded mesh files.
        """
        self._parse_meta(jsonstr)
        self._parse_elements(jsonstr)
        self._parse_connectivity(jsonstr)

    def _parse_meta(self, mesh_dict):
        self.mesh.name = mesh_dict['name']
        self.mesh.dimensions = int(mesh_dict['dimensions'])
        self.mesh.number_of_points = int(mesh_dict['number_of_points'])
        self.mesh.element_counter = int(mesh_dict['element_counter'])
        self.mesh.submesh_counter = int(mesh_dict['submesh_counter'])

    def _parse_elements(self, mesh_dict):
        elems_json = mesh_dict['elements']
        for _en in elems_json.keys():
            e = elems_json[_en].split(' ')
            en = int(_en)
            # set element
            if int(e[0]) == 0:
                # true element
                elem = element_types.create_element(e[1])
            else:
                # submesh
                elem = load_mesh(e[1], self._file_dir)

            self.mesh.elements[en] = elem
            self.mesh.element_points[en] = [(en, i) for i in range(elem.get_number_of_ensemble_points())]

        self.mesh.element_counter = max(list(self.mesh.elements.keys())) + 1

    def _parse_connectivity(self, mesh_dict):
        for k, v in mesh_dict['connectivity'].items():
            _k = tuple(int(ki) for ki in k.split('_'))
            _v = [tuple(int(vii) for vii in vi.split('_')) for vi in v.split()]
            self.mesh.connectivity[_k] = _v


def load_mesh_json(filename, mesh, filedir=None):
    reader = MeshJSONReader(mesh)
    reader.read(filename, filedir)
    return mesh


def save_mesh_json(filename, mesh, filedir=None, subfieldfns=None):
    writer = MeshJSONWriter(mesh)
    return writer.write(filename, filedir, subfieldfns)


def load_mesh(filename, path=None):
    mesh = mesh_ensemble(None, None)
    with open(filename, 'r') as f:
        head = f.read(1)
        if head == '{':
            load_mesh_json(filename, mesh, filedir=path)
        else:
            load_mesh_shelve(filename, mesh, filedir=path)

    return mesh
