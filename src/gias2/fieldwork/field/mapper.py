"""
FILE: mapper.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: class to create maps mapping ensemble parameters to element point 
parameters and vice-versa.

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""
import logging

import scipy

log = logging.getLogger(__name__)


class mapper:
    """ Class to create maps mapping ensemble parameters to element point 
    parameters and vice-versa.
    """

    def __init__(self, debug=0):
        self.mesh = None
        self._ensemble_to_element_map = {}  # { global_point_number: { element_number: {point_number: weighting}, ... } }
        # ~ self._element_to_ensemble_map = {}    # { element_number: { point_number: [[ global_point_number, [weightings] ], [...], ...] } }
        self._element_to_ensemble_map = {}  # { element_number: { point_number: [[ global_point_numbers], [weightings] ], ...] } }        ###
        self.number_of_ensemble_points = 0
        self.number_of_dof = 0
        self._custom_ensemble_order = None  # { self generated number: user assigned number }
        self.debug = debug
        self.has_custom_map = False

        return

    # ==================================================================#
    def set_parent_field(self, parent_field):
        """ set new parent field on which mapping is based.
        parent field should contain topological, and basis information 
        about itself and child fields.
        """

        self._ensemble_to_element_map = {}
        self._element_to_ensemble_map = {}

        self.field = parent_field
        self.number_of_ensemble_points = self.field.mesh.get_number_of_ensemble_points()

        if not self.number_of_ensemble_points:
            log.debug('ERROR: mapper.set_parent_field: empty parent mesh')
            return
        else:
            # initialise ensemble to element map
            for i in range(self.number_of_ensemble_points):
                self._ensemble_to_element_map[i] = {}

            # initialise element to ensemble map
            for i in sorted(self.field.mesh.elements.keys()):
                # initialise element entries
                self._element_to_ensemble_map[i] = {}
                for j in range(self.field.mesh.elements[i].get_number_of_ensemble_points()):
                    # initialise element point entries for each element
                    # ~ self._element_to_ensemble_map[i][j] = []
                    self._element_to_ensemble_map[i][j] = ([], [])

            return 1

    # ==================================================================#
    def do_mapping(self):
        """ map the current mesh topology
        """

        gp = 0  # global points number
        assigned = []  # list of element points already assigned a global number
        ep = list(self.field.mesh.connectivity.keys())  # list of all element points and hanging points
        ep.sort()
        i = 0

        if self.debug:
            log.debug('sorted element points:', ep)

        # account for points connected to hanging nodes
        while ep[i][0] == -1:
            assigned += self.field.mesh.connectivity[ep[i]]
            i += 1

        # map conventional points
        for j in range(i, len(ep)):

            if self.debug:
                log.debug('mapping ep:', ep[j])

            # loop through each element point in the connectivity dict
            (e, p) = ep[j]
            # check if point has already been assigned a global number
            if (e, p) not in assigned:

                if self.debug:
                    log.debug('assigning global', gp, 'to', [(e, p)] + self.field.mesh.connectivity[(e, p)])

                # update ensemble to element map
                try:
                    self._ensemble_to_element_map[gp][e][p] = 1.0
                except KeyError:
                    self._ensemble_to_element_map[gp][e] = {p: 1.0}

                # update element to ensemble map
                # ~ self._element_to_ensemble_map[e][p].append( [ gp, [1.0] ] )
                self._element_to_ensemble_map[e][p][0].append(gp)  ###
                self._element_to_ensemble_map[e][p][1].append(1.0)  ###
                assigned.append((e, p))

                # map points connected to the current element point
                for [ec, pc] in self.field.mesh.connectivity[(e, p)]:
                    # update ensemble to element map
                    try:
                        self._ensemble_to_element_map[gp][ec][pc] = 1.0
                    except KeyError:
                        self._ensemble_to_element_map[gp][ec] = {pc: 1.0}

                    # update element to ensemble map
                    # ~ self._element_to_ensemble_map[ec][pc].append( [ gp, [1.0] ] )
                    self._element_to_ensemble_map[ec][pc][0].append(gp)
                    self._element_to_ensemble_map[ec][pc][1].append(1.0)
                    assigned.append((ec, pc))

                if self.debug:
                    log.debug('assigned:', assigned)
                gp += 1

        if self.debug:
            log.debug('element to ensemble map:', self._element_to_ensemble_map)

        # map hanging points
        if self.debug:
            log.debug('i:', i)
        # for each hanging point
        for h in range(i):
            hp = self.field.mesh.hanging_points[h]

            # get its host_element and element_coord
            host_element = hp.get_host_element()
            host_element_c = hp.get_element_coordinates()

            # find the ensemble points associated with the host element
            host_gp = []
            for element_point in list(self._element_to_ensemble_map[host_element].keys()):
                host_gp += self._element_to_ensemble_map[host_element][element_point][0]  ###
                # ~ for [g, w] in self._element_to_ensemble_map[host_element][element_point]:
                # ~ host_gp.append(g)

            # calculate weightings using basis
            # if local
            weights = self.field.basis[self.field.mesh.elements[host_element].type].eval(host_element_c)
            if self.debug:
                log.debug('host_element_c:', host_element_c, 'weights:', weights)
                log.debug('host_gp:', host_gp)

            if isinstance(weights, int):
                log.debug('ERROR: mapper._map_ensemble_to_element: unable to evaluate weights')
                return
            else:
                # update maps for connected element points
                for [e, p] in self.field.mesh.connectivity[ep[h]]:
                    if e == -1:
                        pass
                    else:
                        for host_gp_i in range(len(host_gp)):
                            # update ensemble to element map
                            try:
                                self._ensemble_to_element_map[host_gp[host_gp_i]][e][p] = weights[host_gp_i]
                            except KeyError:
                                self._ensemble_to_element_map[host_gp[host_gp_i]][e] = {p: weights[host_gp_i]}

                            # update element to ensemble map

                            self._element_to_ensemble_map[e][p][0].append(host_gp[host_gp_i])  ###
                            self._element_to_ensemble_map[e][p][1].append(weights[host_gp_i])  ###

        return 1

    # ==================================================================#
    def remove_element(self, element_number):
        """ removes element from the _element_to_ensemble_map
        """

        if element_number in list(self._element_to_ensemble_map.keys()):

            # get element point to ensemble points mapping
            ep_map = self._element_to_ensemble_map[element_number]
            ensemble_points = []
            for element_point in list(ep_map.values()):  ###
                ensemble_points += element_point[0]  ###
            # ~ for element_point in ep_map.values():
            # ~ ensemble_points += [ p[0] for p in element_point ]

            # remove entries in ensemble_to_element_map    
            for ensemble_point in ensemble_points:
                del self._ensemble_to_element_map[ensemble_point][element_number]
                if self._ensemble_to_element_map[ensemble_point] == {}:
                    del self._ensemble_to_element_map[ensemble_point]

            # remove mapping in _element_to_ensemble_map
            del self._element_to_ensemble_map[element_number]

            return
        else:
            raise ValueError('element ' + str(element_number) + 'does not exist')

    # ==================================================================#
    def get_ensemble_point_element_points(self, ensemble_point_number):
        """ returns a list of ( element n, element point n ) mapped to the
        given ensemble_point_number
        """
        element_points = []
        if self.get_number_of_ensemble_points() == 0:
            return element_points

        try:
            element_map = self._ensemble_to_element_map[ensemble_point_number]
        except KeyError:
            raise ValueError
        else:
            for element in element_map:
                for point in element_map[element]:
                    element_points.append((element, point))

            return element_points

    # ==================================================================#
    def get_element_parameters(self, element_number, parameters, do_hack=0):
        """ Uses mapping to return the element parameters for the
        required element. 
        """

        # get mapping for specified element
        try:
            element_map = self._element_to_ensemble_map[element_number]
        except IndexError:
            raise ValueError('invalid element_number {}'.format(element_number))

        # do this check in the ensemble field function
        # ~ if len( parameters ) != self.get_number_of_ensemble_points():
        # ~ print 'ERROR: mapper.get_element_parameters: incorrect number of parameters. Require', self.get_number_of_ensemble_points()
        # ~ return None

        # get list of element points
        points = list(element_map.keys())
        points.sort()

        element_parameters = []
        if do_hack:
            # use only when absolutely no hanging points!!!! Assumes all weights are 1, and 1 to 1 mapping
            element_parameters = scipy.array([parameters[i] for i in [element_map[ep][0] for ep in points]])
            element_parameters = element_parameters.squeeze()
        else:
            # for each element point
            for element_point in points:

                ensemble_i = element_map[element_point][0]
                if self.has_custom_map:
                    ensemble_i = [self._custom_ensemble_order[i] for i in ensemble_i]

                ensemble_parameters = [parameters[i] for i in ensemble_i]  ###
                weights = element_map[element_point][1]  ###
                element_parameters.append(scipy.dot(weights, ensemble_parameters))  ###

            element_parameters = scipy.array(element_parameters).squeeze()

            # ~ # for each global point mapped to this element point
            # ~ for [g, w] in element_map[ element_point ]:
            # ~
            # ~ # if custom ensemble point ordering, map g to the right ensemble point
            # ~ if self.has_custom_map:
            # ~ g = self._custom_ensemble_order[g]
            # ~
            # ~ # calculated the weighted contribution of this global point
            # ~ mapped_values.append( scipy.multiply( parameters[g], w ) )
        # ~
        # ~ # sum contributions and append element_point parameter to element_parameters
        # ~ element_parameters.append( list( scipy.sum( mapped_values, axis = 0 ) ) )

        return element_parameters

    # ==================================================================#
    def get_number_of_ensemble_points(self):
        return self.number_of_ensemble_points

    def set_number_of_ensemble_points(self, x):
        self.number_of_ensemble_points = x

    # ==================================================================#
    def get_number_of_dof(self):
        return self.number_of_dof

    # ==================================================================#
    def set_custom_ensemble_ordering(self, ordering):
        """ defines custom numbering order to ensemble points {original: custom }
        """
        self.has_custom_map = True
        # check for correct length
        if len(ordering) != self.number_of_ensemble_points:
            raise ValueError('Wrong number of entries in new ordering. Requires {}, got {}'.format(
                self.number_of_ensemble_points, len(ordering)
            ))
        else:
            self._custom_ensemble_order = dict(ordering)
            self._custom_ensemble_order_inverse = {}  # {custom number: original number}
            for k in list(self._custom_ensemble_order.keys()):
                self._custom_ensemble_order_inverse[self._custom_ensemble_order[k]] = k
            return 1
