"""
FILE: template_fields_cubic_simplex.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION:
functions to instantiate ensemble_field_functions with common cubic 
triangular mesh topologies.

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import scipy

import gias2.fieldwork.field.ensemble_field_function as E


# =============#
# empty field #
# =============#
def empty_field(name, dimensions, basis_type):
    """ empty dimensions-D field with empty mesh and basis of type 
    basis_type
    """

    f = E.ensemble_field_function(name, dimensions, debug=0)
    f.set_basis(basis_type)
    f.set_new_mesh(name + '_mesh')

    return f


# ===================#
# 4 -triangle patch #
# ===================#
#       14
#       /\
#     12  13 
#    9/_10_\11
#    /\    /\
#   5  6  7  8
#  /____\/____\
# 0  1  2   3  4

# ~ def four_tri_patch():
# ~
# ~ dimensions = 2
# ~ elements = 4
# ~ element_type = 'tri10'
# ~ basis_type = 'triangular_cubic_2d'
# ~
# ~ f = E.ensemble_field_function('4_tri_patch', dimensions, debug = 0)
# ~ f.set_basis( basis_type )
# ~ f.set_new_mesh('4_tri_patch')
# ~ f.create_elements(element_type, elements)
# ~
# ~ # connect elements
# ~ f.connect_element_points( [(0,2),(1,0),(2,0)] )
# ~ f.connect_element_points( [(0,4),(1,4),(3,0)] )
# ~ f.connect_element_points( [(1,2),(2,4),(3,2)] )
# ~ f.connect_element_points( [(0,3),(1,5)] )
# ~ f.connect_element_points( [(1,1),(2,5)] )
# ~ f.connect_element_points( [(1,3),(3,1)] )
# ~
# ~ f.map_parameters()
# ~ f.mapper.set_custom_ensemble_ordering( { 0:0, 1:1, 2:2, 3:6, 4:9, 5:5 ,6:7 ,7:11, 8:10,\
# ~ 9:3, 10:4, 11:8, 12:13, 13:14, 14:12 } )
# ~
# ~ return f

# ===================#
# 7-triangle patch #
# ===================#
#       22_23_ 24
#       /\    /\
#     18  19 20 \21 
#   13/_14_\15_16\17
#    /\    /\    /\
#   7  8  9  10 11 12
#  /____\/____\/____\ 
# 0  1  2   3  4  5  6

# ~ def eight_tri_patch():
# ~
# ~ dimensions = 2
# ~ elements = 8
# ~ element_type = 'tri6'
# ~ basis_type = 'triangular_quadratic_2d'
# ~
# ~ f = E.ensemble_field_function('8_tri_patch', dimensions, debug = 0)
# ~ f.set_basis( basis_type )
# ~ f.set_new_mesh('8_tri_patch')
# ~ f.create_elements(element_type, elements)
# ~
# ~ # connect elements
# ~ f.connect_element_points( [(0,2),(1,0),(2,0)] )
# ~ f.connect_element_points( [(2,2),(3,0),(4,0)] )
# ~
# ~ f.connect_element_points( [(0,3),(1,5)] )
# ~ f.connect_element_points( [(1,1),(2,5)] )
# ~ f.connect_element_points( [(2,3),(3,5)] )
# ~ f.connect_element_points( [(3,1),(4,5)] )
# ~
# ~ f.connect_element_points( [(0,4),(1,4),(5,0)] )
# ~ f.connect_element_points( [(1,3),(5,1)] )
# ~ f.connect_element_points( [(1,2),(2,4),(3,4),(5,2),(6,0),(7,0)] )
# ~ f.connect_element_points( [(3,3),(7,1)] )
# ~ f.connect_element_points( [(3,2),(4,4),(7,2)] )
# ~
# ~ f.connect_element_points( [(5,3),(6,5)] )
# ~ f.connect_element_points( [(6,1),(7,5)] )
# ~
# ~ f.connect_element_points( [(5,4),(6,4)] )
# ~ f.connect_element_points( [(6,2),(7,4)] )
# ~
# ~ f.map_parameters()
# ~ f.mapper.set_custom_ensemble_ordering( { 0:0, 1:1, 2:2, 3:8, 4:13, 5:7 ,6:9,\
# ~ 7:15, 8:14, 9:3, 10:4, 11:10, 12:11,\
# ~ 13:17, 14:16, 15:5, 16:6, 17:12,\
# ~ 18:19, 19:22, 20:18, 21:20,\
# ~ 22:24, 23:23, 24:21 } )
# ~
# ~ return f


# ==================#
# 2 -triangle quad #
# ==================#
# 6 _7_ 8       ______   
#  |\  |  or   /\    /
# 3| 4 |5     /  \  /
#  |__\|     /____\/
# 0  1  2

# ~ def two_quad_patch():
# ~
# ~ dimensions = 2
# ~ elements = 2
# ~ element_type = 'tri6'
# ~ basis_type = 'triangular_quadratic_2d'
# ~
# ~ f = E.ensemble_field_function('2_quad_patch', dimensions, debug = 0)
# ~ f.set_basis( basis_type )
# ~ f.set_new_mesh('2_quad_patch')
# ~ f.create_elements(element_type, elements)
# ~
# ~ # connect elements
# ~ f.connect_element_points( [(0,2),(1,0)] )
# ~ f.connect_element_points( [(0,3),(1,5)] )
# ~ f.connect_element_points( [(0,4),(1,4)] )
# ~
# ~ f.map_parameters()
# ~ f.mapper.set_custom_ensemble_ordering({ 0:0, 1:1, 2:2, 3:4, 4:6, 5:3, 6:5, 7:8, 8:7 })
# ~
# ~ return f

# ==================#
# triangular strip #
# ==================#

def triangular_strip(n):
    dimensions = 2
    elements = n
    element_type = 'tri10'
    basis_type = {'tri10': 'simplex_cubic_lagrange_2d'}

    f = E.ensemble_field_function('triangular_strip_' + str(n), dimensions, debug=0)
    f.set_basis(basis_type)
    f.set_new_mesh('triangular_strip_' + str(n))
    e_i = f.create_elements(element_type, elements)

    # connections
    for i in e_i:
        if i == 0:
            pass
        else:
            f.connect_elements(i, 2, i - 1, 1)

    f.map_parameters()

    # remap
    if scipy.mod(n, 2):
        t1Max = (n / 2 + 1) * 3
        t2Max = (n / 2 + 1) * 6
        t3Max = (n / 2) * 9 + 8
    else:
        t1Max = (n / 2) * 3
        t2Max = n * 3 + 1
        t3Max = (n / 2) * 9 + 2

    # element 0 node map
    t1New = 4
    t2New = t1Max + 4
    t3New = t2Max + 3
    t4New = t3Max + 2
    nodeMap = {0: 0, 1: 1, 2: 2, 3: 3, 4: t1Max + 1, 5: t1Max + 2, 6: t1Max + 3, 7: t2Max + 1, 8: t2Max + 2,
               9: t3Max + 1}
    # remap each element
    for i in range(1, n):
        ePoints = f.mapper._element_to_ensemble_map[i]
        if scipy.mod(i, 2):
            nodeMap[ePoints[1][0][0]] = t4New
            nodeMap[ePoints[2][0][0]] = t4New + 1
            nodeMap[ePoints[3][0][0]] = t4New + 2
            nodeMap[ePoints[5][0][0]] = t3New
            nodeMap[ePoints[6][0][0]] = t3New + 1
            nodeMap[ePoints[8][0][0]] = t2New
            t2New += 1
            t3New += 2
            t4New += 3
        else:
            nodeMap[ePoints[1][0][0]] = t1New
            nodeMap[ePoints[2][0][0]] = t1New + 1
            nodeMap[ePoints[3][0][0]] = t1New + 2
            nodeMap[ePoints[5][0][0]] = t2New
            nodeMap[ePoints[6][0][0]] = t2New + 1
            nodeMap[ePoints[8][0][0]] = t3New
            t1New += 3
            t2New += 2
            t3New += 1

    f.mapper.set_custom_ensemble_ordering(nodeMap)
    return f


def triangular_strip2(n):
    """
                e2 
         21_______24_______ 27
          /\      /\      /
       14/  \   17  \    / 20
    e3 7/    \ 10    \  / 13   e1
       /______\/______\/
      0       3        6
                e0
    
    adjacent triangle elements are flipped so that xi1 and xi2 are
    always along the orthogonal edges of right angle triangles, this
    ensures continuous derivatives and curvature
    """

    dimensions = 2
    elements = n
    element_type = 'tri10'
    basis_type = {'tri10': 'simplex_cubic_lagrange_2d'}

    f = E.ensemble_field_function('triangular_strip_' + str(n), dimensions, debug=0)
    f.set_basis(basis_type)
    f.set_new_mesh('triangular_strip_' + str(n))
    e_i = f.create_elements(element_type, elements)

    # connections
    for i in range(n):
        if i == (n - 1):
            pass
        elif scipy.mod(i, 2):
            f.connect_elements(i, 2, i + 1, 2, invert_order=True)
        else:
            f.connect_elements(i, 1, i + 1, 1, invert_order=True)

    f.map_parameters()

    # remap
    if scipy.mod(n, 2):
        t1Max = (n / 2 + 1) * 3
        t2Max = (n / 2 + 1) * 6
        t3Max = (n / 2) * 9 + 8
    else:
        t1Max = (n / 2) * 3
        t2Max = n * 3 + 1
        t3Max = (n / 2) * 9 + 2

    # element 0 node map
    nodeMap = {0: 0, 1: 1, 2: 2, 3: 3, 4: t1Max + 1, 5: t1Max + 2, 6: t1Max + 3, 7: t2Max + 1, 8: t2Max + 2,
               9: t3Max + 1}
    t1New = 4
    t2New = t1Max + 4
    t3New = t2Max + 3
    t4New = t3Max + 2
    # remap each element
    for i in range(1, n):
        ePoints = f.mapper._element_to_ensemble_map[i]
        if scipy.mod(i, 2):
            nodeMap[ePoints[0][0][0]] = t4New + 2
            nodeMap[ePoints[1][0][0]] = t4New + 1
            nodeMap[ePoints[2][0][0]] = t4New
            nodeMap[ePoints[4][0][0]] = t3New + 1
            nodeMap[ePoints[5][0][0]] = t3New
            nodeMap[ePoints[7][0][0]] = t2New
            t2New += 1
            t3New += 2
            t4New += 3
        else:
            nodeMap[ePoints[1][0][0]] = t1New
            nodeMap[ePoints[2][0][0]] = t1New + 1
            nodeMap[ePoints[3][0][0]] = t1New + 2
            nodeMap[ePoints[5][0][0]] = t2New
            nodeMap[ePoints[6][0][0]] = t2New + 1
            nodeMap[ePoints[8][0][0]] = t3New
            t1New += 3
            t2New += 2
            t3New += 1

    f.mapper.set_custom_ensemble_ordering(nodeMap)
    nTot = f.get_number_of_ensemble_points()
    f.edge_points = (tuple(range(0, t1Max + 1)), (t1Max, t2Max, t3Max, nTot - 1), tuple(range(t3Max + 1, nTot)),
                     (0, t1Max + 1, t2Max + 1, t3Max + 1))
    return f


def triangular_strip2_ring(n):
    """ 
                e1 
         18_______21_______ 18
          /\      /\      /
       12/  \   15  \    / 12
       6/    \  9    \  / 6 
       /______\/______\/
      0       3        0
                e0
    
    adjacent triangle elements are flipped so that xi1 and xi2 are
    always along the orthogonal edges of right angle triangles, this
    ensures continuous derivatives and curvature.
    
    n must be even
    """

    if scipy.mod(n, 2):
        raise ValueError('n must be even')

    # make strip
    r = triangular_strip2(n)

    # connect ends - calculate end node numbers
    r.connect_elements(n - 1, 2, 0, 2, invert_order=True)
    r.map_parameters()

    # remap
    nTot = r.get_number_of_ensemble_points()
    t1Max = nTot / 4 - 1
    t2Max = nTot / 2 - 1
    t3Max = nTot / 4 * 3 - 1

    # element 0 node map
    nodeMap = {0: 0, 1: 1, 2: 2, 3: 3, 4: t1Max + 1, 5: t1Max + 2, 6: t1Max + 3, 7: t2Max + 1, 8: t2Max + 2,
               9: t3Max + 1}
    t1New = 4
    t2New = t1Max + 4
    t3New = t2Max + 3
    t4New = t3Max + 2
    # remap each element
    for i in range(1, n - 1):
        ePoints = r.mapper._element_to_ensemble_map[i]
        if scipy.mod(i, 2):
            nodeMap[ePoints[0][0][0]] = t4New + 2
            nodeMap[ePoints[1][0][0]] = t4New + 1
            nodeMap[ePoints[2][0][0]] = t4New
            nodeMap[ePoints[4][0][0]] = t3New + 1
            nodeMap[ePoints[5][0][0]] = t3New
            nodeMap[ePoints[7][0][0]] = t2New
            t2New += 1
            t3New += 2
            t4New += 3
        else:
            nodeMap[ePoints[1][0][0]] = t1New
            nodeMap[ePoints[2][0][0]] = t1New + 1
            nodeMap[ePoints[3][0][0]] = t1New + 2
            nodeMap[ePoints[5][0][0]] = t2New
            nodeMap[ePoints[6][0][0]] = t2New + 1
            nodeMap[ePoints[8][0][0]] = t3New
            t1New += 3
            t2New += 2
            t3New += 1

    # last element node map
    ePoints = r.mapper._element_to_ensemble_map[n - 1]
    nodeMap[ePoints[1][0][0]] = t4New + 1
    nodeMap[ePoints[2][0][0]] = t4New
    nodeMap[ePoints[5][0][0]] = t3New
    # correct 2nd to last element node 3 mapping to ensemble point 0
    nodeMap[0] = 0

    r.mapper.set_custom_ensemble_ordering(nodeMap)
    r.edge_points = (tuple(range(0, t1Max + 1)), tuple(range(t3Max + 1, nTot)))
    return r

# ==================#
# 3 -triangle quad #
# ==================#
#      9__10__11   
#      /\    /\
#     5  6  7  8
#    /____\/____\
#   0  1  2   3  4

# ~ def three_quad_patch():
# ~
# ~ dimensions = 2
# ~ elements = 3
# ~ element_type = 'tri6'
# ~ basis_type = 'triangular_quadratic_2d'
# ~
# ~ f = E.ensemble_field_function('3_quad_patch', dimensions, debug = 0)
# ~ f.set_basis( basis_type )
# ~ f.set_new_mesh('3_quad_patch')
# ~ f.create_elements(element_type, elements)
# ~
# ~ # connect elements
# ~ f.connect_element_points( [(0,2),(1,0),(2,0)] )
# ~ f.connect_element_points( [(0,3),(1,5)] )
# ~ f.connect_element_points( [(0,4),(1,4)] )
# ~ f.connect_element_points( [(1,1),(2,5)] )
# ~ f.connect_element_points( [(1,2),(2,4)] )
# ~
# ~ f.map_parameters()
# ~ f.mapper.set_custom_ensemble_ordering({ 0:0, 1:1, 2:2, 3:6, 4:9, 5:5,\
# ~ 6:7, 7:11, 8:10,\
# ~ 9:3, 10:4, 11:8 })
# ~
# ~ return f

# ==========================#
# ring of 3-triangle quads #
# ==========================#
#  28____29___30____31___32____33____34___35____ 28
#   |\        /|\        /|\        /|\        /|
#   | \      / | \      / | \      / | \      / |
# 16|  17  18  19 20  21  22 23  24  25 26  27  |16
#   |   \  /   |   \  /   |   \  /   |   \  /   |
#   |____\/____|____\/____|____\/____|____\/____|
#  0   1  2  3 4  5  6  7 8  9 10 11 12 13 14 15 0

# ~ def three_quad_ring( number_of_quads, height, radius, debug=0 ):
# ~ """
# ~ number_of_quads defines the number of 2-triangle quads in the ring
# ~ ( number of triangles / 2 )
# ~ height and radius determine the x,y,z parameters returned
# ~ with the ensemble_field object. z is height
# ~ """
# ~
# ~ if debug:
# ~ print 'Generating three_quad_ring with:', number_of_quads, height, radius, '\n'
# ~
# ~ dimensions = 2
# ~ elements = number_of_quads
# ~ basis_type = 'triangular_quadratic_2d'
# ~ npoints = elements*6
# ~
# ~ # initialise main field ensemble
# ~ ring = E.ensemble_field_function('3_quad_ring', dimensions, debug = 0)
# ~ ring.set_basis( basis_type )
# ~ ring.set_new_mesh('3_quad_ring')
# ~
# ~ # create 3 triangle quads for the ring and add to main field ensemble
# ~ for i in range( elements ):
# ~ ring.add_element( three_quad_patch() )
# ~
# ~ # connect up quads
# ~ for i in range( elements-1 ):
# ~ ring.connect_element_points( [(i,4), (i+1,0)] )
# ~ ring.connect_element_points( [(i,8), (i+1,5)] )
# ~ ring.connect_element_points( [(i,11), (i+1,9)] )
# ~
# ~ # join end of strip to start
# ~ ring.connect_element_points( [(i+1,4), (0,0)] )
# ~ ring.connect_element_points( [(i+1,8), (0,5)] )
# ~ ring.connect_element_points( [(i+1,11), (0,9)] )
# ~
# ~ ring.map_parameters()
# ~
# ~ #==================================================================#
# ~ # remap ensemble point numbering so that node numbering increase first
# ~ # up each column of nodes
# ~ remap = {}
# ~ e = 0
# ~ old = 0
# ~ new_t1 = 0
# ~ new_t2 = elements * 4
# ~ new_t3 = elements * 7
# ~
# ~ for e in range( 0, elements ):
# ~
# ~ if e == 0:
# ~ for old in range( old, old + 5 ):
# ~ remap[old] = new_t1
# ~ new_t1 += 1
# ~
# ~ for old in range( old + 1, old + 5 ):
# ~ remap[old] = new_t2
# ~ new_t2 += 1
# ~
# ~ for old in range( old + 1, old + 4 ):
# ~ remap[old] = new_t3
# ~ new_t3 += 1
# ~ elif e == elements - 1:
# ~ old = 12 + ( e - 1 ) * 9
# ~ for old in range( old, old + 3 ):
# ~ remap[old] = new_t1
# ~ new_t1 += 1
# ~
# ~ for old in range( old + 1, old + 3 ):
# ~ remap[old] = new_t2
# ~ new_t2 += 1
# ~
# ~ for old in range( old + 1, old + 2 ):
# ~ remap[old] = new_t3
# ~ new_t3 += 1
# ~ else:
# ~ old = 12 + ( e - 1 ) * 9
# ~ for old in range( old, old + 4 ):
# ~ remap[old] = new_t1
# ~ new_t1 += 1
# ~
# ~ for old in range( old + 1, old + 4 ):
# ~ remap[old] = new_t2
# ~ new_t2 += 1
# ~
# ~ for old in range( old + 1, old + 3 ):
# ~ remap[old] = new_t3
# ~ new_t3 += 1
# ~
# ~ ring.mapper.set_custom_ensemble_ordering( remap )
# ~
# ~ #==================================================================#
# ~
# ~ if height > 0.0 and radius > 0.0:
# ~ # generate parameters
# ~ theta_t1 = scipy.linspace(0, scipy.pi*2.0, elements * 4 + 1)
# ~ theta_t2 = scipy.linspace(0, scipy.pi*2.0, elements * 3 + 1)
# ~ theta_t3 = scipy.linspace(0, scipy.pi*2.0, elements * 2 + 1)
# ~ theta = [theta_t1, theta_t2, theta_t3]
# ~ Z = scipy.linspace( 0.0, height, 3 )
# ~
# ~ # generate parameters lists
# ~ xparam = []
# ~ yparam = []
# ~ zparam = []
# ~
# ~ for i in range( 3 ):
# ~ x = radius*scipy.cos(theta[i])
# ~ y = radius*scipy.sin(theta[i])
# ~ z = Z[i]
# ~ for i in range( len( theta[i] ) - 1 ):
# ~ xparam.append( x[i] )
# ~ yparam.append( y[i] )
# ~ zparam.append( z )
# ~
# ~ return ( ring, xparam, yparam, zparam )
# ~ else:
# ~ return ring
# ~
# ~ #==========================#
# ~ # ring of 2-triangle quads #
# ~ #==========================#
# ~ #
# ~ # 16 _17_18_19_20_21_22_23__ 16  radius[2]
# ~ #    |\    |\    |\    |\    |
# ~ #   | \   | \   | \   | \   |
# ~ #  8|  9 10 11 12 13 14 15  |8   radius[1]
# ~ #    |   \ |   \ |   \ |   \ |
# ~ #   |____\|____\|____\|____\|
# ~ #  0   1  2  3  4  5  6  7   0   radius[0]
# ~
# ~ def two_quad_ring( number_of_quads, height = None, radius = None ):
# ~ """
# ~ number_of_quads defines the number of 2-triangle quads in the ring
# ~ ( number of triangles / 2 ).
# ~ height and radius determine the x,y,z parameters returned.
# ~ radius can either be scalar, where the cylinder will have uniform
# ~ radius along its length, or a list of 3 values.
# ~ """
# ~
# ~ # initialise main field ensemble
# ~ dimensions = 2
# ~ elements = number_of_quads
# ~ basis_type = 'triangular_quadratic_2d'
# ~ npoints = elements*6
# ~
# ~ ring = E.ensemble_field_function('2_quad_ring', dimensions, debug = 0)
# ~ ring.set_basis( basis_type )
# ~ ring.set_new_mesh('2_quad_ring')
# ~
# ~ # create 2 triangle quads for the ring and add to main field ensemble
# ~ for i in range( elements ):
# ~ ring.add_element( two_quad_patch() )
# ~
# ~ # connect up quads
# ~ for i in range( elements-1 ):
# ~ ring.connect_element_points( [(i,2), (i+1,0)] )
# ~ ring.connect_element_points( [(i,5), (i+1,3)] )
# ~ ring.connect_element_points( [(i,8), (i+1,6)] )
# ~
# ~ # join end of strip to start
# ~ ring.connect_element_points( [(i+1,2), (0,0)] )
# ~ ring.connect_element_points( [(i+1,5), (0,3)] )
# ~ ring.connect_element_points( [(i+1,8), (0,6)] )
# ~
# ~ ring.map_parameters()
# ~
# ~ #==================================================================#
# ~ # remap ensemble point numbering so that node numbering increase first
# ~ # up each column of nodes
# ~ remap = _two_quad_ring_remapper( number_of_quads )
# ~ ring.mapper.set_custom_ensemble_ordering( remap )
# ~
# ~ #==================================================================#
# ~ if height!=None and radius!=None:
# ~ # for quadratic quads, there are elements*2 nodes around the circle
# ~ # generate coordinates for nodes
# ~
# ~ if isinstance( radius, float ) or isinstance( radius, int ):
# ~ radius = [radius]*3
# ~ elif len( radius ) != 3:
# ~ print 'ERROR: radius must be scalar or 3 long'
# ~ return None
# ~
# ~ theta = scipy.linspace( 0.0, 2.0*scipy.pi, elements*2 + 1 )
# ~ zcoords = scipy.linspace( 0.0, height, 3 )
# ~ xparam = []
# ~ yparam = []
# ~ zparam = []
# ~
# ~ for row in range( 3 ):
# ~ xcoords = radius[row] * scipy.cos(theta)
# ~ ycoords = radius[row] * scipy.sin(theta)
# ~ for i in range( len( xcoords ) - 1 ):
# ~ xparam.append( [xcoords[i]] )
# ~ yparam.append( [ycoords[i]] )
# ~ zparam.append( [float(zcoords[row])] )
# ~
# ~ return ( ring, xparam, yparam, zparam )
# ~ else:
# ~ return ring
# ~
# ~
# ~ #========#
# ~ # sphere #
# ~ #========#
# ~
# ~ def sphere( azimuth_divs, incline_divs, radius, inclination_max ):
# ~ """ Full sphere is inclination_max = pi, else sphere is truncated.
# ~ Sphere mesh is composed of tri_hemispheres for the top and bottom
# ~ caps and three_quad_rings in between.
# ~
# ~ Orientation and ensemble point numbering of the rings are flipped at
# ~ the equator. Ensemble point numbering starts at the bottom and
# ~ increases latitudinally (around the sphere) before moving up
# ~ longitudinally.
# ~
# ~ Note: if inclination_max == pi, incline_divs must be even
# ~ """
# ~
# ~ dimensions = 2
# ~ basis_type = 'triangular_quadratic_2d'
# ~
# ~ # initialise main field ensemble
# ~ sphere = E.ensemble_field_function('sphere', dimensions, debug = 0)
# ~ sphere.set_basis( basis_type )
# ~ sphere.set_new_mesh('sphere')
# ~
# ~ # determine number of tiers in the upper hemisphere
# ~ # assign tiers proportional to 90/inclination_max if inclination > 90
# ~ if inclination_max > scipy.pi/2.0:
# ~ upper_tiers = int( round( incline_divs * (scipy.pi/2.0) / inclination_max ) )
# ~ lower_tiers = incline_divs - upper_tiers
# ~ else:
# ~ upper_tiers = incline_divs
# ~ lower_tiers = 0
# ~
# ~ # rings and their number of 3-tri quad elements
# ~ ring_elements = {}
# ~ ring_counter = 0
# ~ if inclination_max == scipy.pi:
# ~ lower_rings = lower_tiers - 1
# ~ else:
# ~ lower_rings = lower_tiers
# ~
# ~ upper_rings = upper_tiers - 1
# ~
# ~ # if truncated sphere, return the number of points at the base of the sphere for connection
# ~ # their numbers will be 0 to connect_points - 1
# ~ connect_points = None
# ~
# ~ #==================================================================#
# ~ # add caps and rings
# ~ # lower hemisphere, node ordering need to be reversed
# ~ # if a full sphere
# ~ if inclination_max == scipy.pi:
# ~ # lower cap first
# ~ lower_cap = tri_hemisphere( 0.0, azimuth_divs )
# ~ remap = _tri_hemisphere_reverse_mapper( azimuth_divs )
# ~ lower_cap.mapper.set_custom_ensemble_ordering( remap )
# ~ sphere.add_element( lower_cap )
# ~
# ~ # then lower rings
# ~ for i in range( lower_tiers - 1 ):
# ~ ring = three_quad_ring( azimuth_divs * 2**i, 0.0, 0.0 )
# ~ ring.mapper.set_custom_ensemble_ordering( _flip_3tri_ring_mapping( ring.mapper._custom_ensemble_order, azimuth_divs * 2**i ) )
# ~ sphere.add_element( ring )
# ~
# ~ ring_elements[ring_counter] = azimuth_divs * 2**i
# ~ ring_counter += 1
# ~
# ~ # if an incomplete sphere
# ~ else:
# ~ # just lower rings
# ~ for i in range( lower_rings ):
# ~ ring = three_quad_ring( azimuth_divs * 2**( (upper_rings - 1) - (lower_rings - 1 - i) ), 0.0, 0.0 )
# ~ ring.mapper.set_custom_ensemble_ordering( _flip_3tri_ring_mapping( ring.mapper._custom_ensemble_order, azimuth_divs * 2**( (upper_rings - 1) - (lower_rings - 1 - i) ) ) )
# ~ sphere.add_element( ring )
# ~
# ~ ring_elements[ring_counter] = azimuth_divs * 2**( (upper_rings - 1) - (lower_rings - 1 - i) )
# ~ ring_counter += 1
# ~
# ~ connect_points = ring_elements[0] * 2
# ~
# ~ # upper hemisphere
# ~ if upper_tiers > 1:
# ~ for i in range( upper_rings ):
# ~ sphere.add_element( three_quad_ring( azimuth_divs * 2**(upper_rings - 1 - i), 0.0, 0.0 ) )
# ~
# ~ ring_elements[ring_counter] = azimuth_divs * 2**(upper_rings - 1 - i)
# ~ ring_counter += 1
# ~
# ~ # upper cap
# ~ sphere.add_element( tri_hemisphere( 0.0, azimuth_divs ) )
# ~
# ~ #==================================================================#
# ~ # connect rings and caps
# ~ # connect lower hemisphere
# ~ element_i = 0
# ~
# ~ if inclination_max == scipy.pi:
# ~ # connect lower cap to closest ring, or upper cap
# ~ for i in range( azimuth_divs * 2 ):
# ~ if not sphere.connect_element_points( [(0, 1+azimuth_divs+i), (1, i)] ):
# ~ print 'ERROR: tri_shapes.sphere: connecting failed'
# ~ return None
# ~
# ~ element_i += 1
# ~
# ~ # connect lower rings
# ~ for lower_ring in range( lower_rings ):
# ~ start = ring_elements[lower_ring] * 5
# ~ for i in range( ring_elements[lower_ring] * 4 ):
# ~ if not sphere.connect_element_points( [(element_i, start+i), (element_i+1, i)] ):
# ~ print 'ERROR: tri_shapes.sphere: connecting failed:', [(element_i, start+i), (element_i+1, i)]
# ~ print 'start:', start
# ~ print 'i_range:', ring_elements[lower_ring] * 4
# ~ return None
# ~
# ~ element_i += 1
# ~
# ~ # connect upper rings and upper cap
# ~ for upper_ring in range( upper_rings ):
# ~ start = ring_elements[upper_ring + lower_rings] * 7
# ~ for i in range( ring_elements[upper_ring + lower_rings] * 2 ):
# ~ if not sphere.connect_element_points( [(element_i, start+i), (element_i+1, i)] ):
# ~ print 'ERROR: tri_shapes.sphere: connecting failed:', [(element_i, start+i), (element_i+1, i)]
# ~ print 'start:', start
# ~ print 'i_range:', ring_elements[upper_ring + lower_rings + 1] * 2
# ~ return None
# ~
# ~ element_i += 1
# ~
# ~ #==================================================================#
# ~ # map, no remap needed
# ~ sphere.map_parameters()
# ~
# ~ #==================================================================#
# ~ # generate parameters
# ~ xparams = []
# ~ yparams = []
# ~ zparams = []
# ~
# ~ phi = scipy.linspace( inclination_max, 0.0, incline_divs*2+1 )
# ~ zcoords = radius * scipy.cos( phi )
# ~ z_i = 0
# ~
# ~ # if lower cap
# ~ if inclination_max == scipy.pi:
# ~ # apex
# ~ xparams.append( [0.0] )
# ~ yparams.append( [0.0] )
# ~ zparams.append( [zcoords[0]] )
# ~
# ~ # 1st row of nodes
# ~ theta = scipy.linspace( 0, 2.0*scipy.pi, azimuth_divs + 1 )
# ~ xcoords = radius * scipy.cos( theta ) * scipy.sin( phi[1] )
# ~ ycoords = radius * scipy.sin( theta ) * scipy.sin( phi[1] )
# ~
# ~ for i in range( azimuth_divs ):
# ~ xparams.append( [xcoords[i]] )
# ~ yparams.append( [ycoords[i]])
# ~ zparams.append( [zcoords[1]] )
# ~
# ~ z_i = 2
# ~
# ~ # lower rings
# ~ for ring in range(lower_rings):
# ~ # for each row of nodes (lower 2 rows)
# ~ for row in range( 2 ):
# ~ theta = scipy.linspace( 0.0, 2.0*scipy.pi, ring_elements[ring] * (row+2) + 1 )
# ~ xcoords = radius * scipy.cos( theta ) * scipy.sin( phi[z_i] )
# ~ ycoords = radius * scipy.sin( theta ) * scipy.sin( phi[z_i] )
# ~
# ~ for i in range( ring_elements[ring] * (row+2) ):
# ~ xparams.append( [xcoords[i]] )
# ~ yparams.append( [ycoords[i]] )
# ~ zparams.append( [zcoords[z_i]] )
# ~
# ~ z_i += 1
# ~
# ~ # upper rings
# ~ for ring in range(upper_rings):
# ~ # for each row of nodes (lower 2 rows)
# ~ for row in range( 2 ):
# ~ theta = scipy.linspace( 0.0, 2.0*scipy.pi, ring_elements[lower_rings + ring] * (4-row) + 1 )
# ~ xcoords = radius * scipy.cos( theta ) * scipy.sin( phi[z_i] )
# ~ ycoords = radius * scipy.sin( theta ) * scipy.sin( phi[z_i] )
# ~
# ~ for i in range( ring_elements[lower_rings + ring] * (4-row) ):
# ~ xparams.append( [xcoords[i]] )
# ~ yparams.append( [ycoords[i]] )
# ~ zparams.append( [zcoords[z_i]] )
# ~
# ~ z_i += 1
# ~
# ~ # upper cap
# ~ for row in range(2):
# ~ theta = scipy.linspace( 0, 2.0*scipy.pi, azimuth_divs * (2 - row) + 1 )
# ~ xcoords = radius * scipy.cos( theta ) * scipy.sin( phi[z_i] )
# ~ ycoords = radius * scipy.sin( theta ) * scipy.sin( phi[z_i] )
# ~
# ~ for i in range( azimuth_divs * (2 - row) ):
# ~ xparams.append( [xcoords[i]] )
# ~ yparams.append( [ycoords[i]] )
# ~ zparams.append( [zcoords[z_i]] )
# ~
# ~ z_i += 1
# ~
# ~ # apex
# ~ xparams.append( [0.0] )
# ~ yparams.append( [0.0] )
# ~ zparams.append( [zcoords[-1]] )
# ~
# ~ #==================================================================#
# ~
# ~ return (sphere, xparams, yparams, zparams, connect_points)
# ~
# ~
# ~
# ~ #================#
# ~ # tri_hemisphere #
# ~ #================#
# ~ #
# ~ #          2
# ~ #    3    |     1
# ~ #         9
# ~ #         |
# ~ #  4--10--12---8--0
# ~ #          |
# ~ #          11
# ~ #     5    |     7
# ~ #         6
# ~ #
# ~
# ~ def tri_hemisphere(radius, elements):
# ~
# ~ dimensions = 2
# ~ element_type = 'tri6'
# ~ basis_type = 'triangular_quadratic_2d'
# ~
# ~ # initialise main field ensemble
# ~ hemi = E.ensemble_field_function('hemisphere', dimensions, debug = 0)
# ~ hemi.set_basis( basis_type )
# ~ hemi.set_new_mesh('hemisphere')
# ~
# ~ hemi.create_elements( element_type, elements )
# ~
# ~ # connect elements
# ~ hemi.connect_element_points( [ (i,4) for i in range( elements ) ] ) # apex
# ~
# ~ for i in range( elements-1 ):
# ~ hemi.connect_element_points( [(i,3), (i+1,5)] )
# ~ hemi.connect_element_points( [(i,2), (i+1,0)] )
# ~
# ~ hemi.connect_element_points( [(i+1,3), (0,5)] )
# ~ hemi.connect_element_points( [(i+1,2), (0,0)] )
# ~
# ~ hemi.map_parameters()
# ~
# ~ #==================================================================#
# ~ # remap
# ~ remap = {}
# ~
# ~ # 1st elements
# ~ remap[0] = 0
# ~ remap[1] = 1
# ~ remap[2] = 2
# ~ remap[5] = elements * 2
# ~ remap[3] = elements * 2 + 1
# ~
# ~ # other elements
# ~ old_t1 = 6
# ~ old_t2 = 8
# ~ new_t1 = 3
# ~ new_t2 = elements * 2 + 2
# ~
# ~ for i in range( 1, elements ):
# ~ remap[old_t1] = new_t1
# ~ if i < elements - 1:
# ~ remap[old_t1 + 1] = new_t1 + 1
# ~ remap[old_t2] = new_t2
# ~
# ~ old_t1 += 3
# ~ new_t1 += 2
# ~ old_t2 += 3
# ~ new_t2 += 1
# ~
# ~ # apex
# ~ remap[4] = new_t2 - 1
# ~
# ~ hemi.mapper.set_custom_ensemble_ordering( remap )
# ~
# ~ #==================================================================#
# ~
# ~ if radius > 0.0:
# ~ # generate parameters
# ~ theta = scipy.linspace( 0, scipy.pi/2.0, 3 )
# ~ phi = scipy.linspace(0, scipy.pi*2.0, elements * 2 + 1)
# ~ xparam = []
# ~ yparam = []
# ~ zparam = []
# ~
# ~ # tier 1
# ~ z = radius * scipy.sin( theta[0] )
# ~ for i in range( elements * 2 ):
# ~ xparam.append( [ radius * scipy.cos( theta[0] ) * scipy.cos( phi[i] ) ] )
# ~ yparam.append( [ radius * scipy.cos( theta[0] ) * scipy.sin( phi[i] ) ] )
# ~ zparam.append( [ z ] )
# ~
# ~ # tier 2
# ~ z = radius * scipy.sin( theta[1] )
# ~ for i in range( 0, elements * 2, 2 ):
# ~ xparam.append( [ radius * scipy.cos( theta[1] ) * scipy.cos( phi[i] ) ] )
# ~ yparam.append( [ radius * scipy.cos( theta[1] ) * scipy.sin( phi[i] ) ] )
# ~ zparam.append( [ z ] )
# ~
# ~ # apex
# ~ xparam.append( [0.0] )
# ~ yparam.append( [0.0] )
# ~ zparam.append( [ radius * scipy.sin( theta[2] ) ])
# ~
# ~ return (hemi, xparam, yparam, zparam )
# ~ else:
# ~ return hemi
# ~
# ~ def four_tri_patch_hemisphere( radius, elements ):
# ~
# ~ dimensions = 2
# ~ basis_type = 'triangular_quadratic_2d'
# ~
# ~ # initialise main field ensemble
# ~ hemi = E.ensemble_field_function('4_tri_patch_hemisphere', dimensions, debug = 0)
# ~ hemi.set_basis( basis_type )
# ~ hemi.set_new_mesh('4_tri_patch_hemisphere')
# ~
# ~ for i in range( elements ):
# ~ hemi.add_element( four_tri_patch() )
# ~
# ~ # connect elements
# ~ hemi.connect_element_points( [ (i,14) for i in range( elements ) ] ) # apex
# ~
# ~ for i in range( elements-1 ):
# ~ hemi.connect_element_points( [(i,4), (i+1,0)] )
# ~ hemi.connect_element_points( [(i,8), (i+1,5)] )
# ~ hemi.connect_element_points( [(i,11), (i+1,9)] )
# ~ hemi.connect_element_points( [(i,13), (i+1,12)] )
# ~
# ~ hemi.connect_element_points( [(i+1,4), (0,0)] )
# ~ hemi.connect_element_points( [(i+1,8), (0,5)] )
# ~ hemi.connect_element_points( [(i+1,11), (0,9)] )
# ~ hemi.connect_element_points( [(i+1,13), (0,12)] )
# ~
# ~ hemi.map_parameters()
# ~
# ~ #==================================================================#
# ~ # remap
# ~ remap = _four_tri_patch_hemisphere_remapper( elements )
# ~ hemi.mapper.set_custom_ensemble_ordering( remap )
# ~
# ~ #==================================================================#
# ~ if radius > 0.0:
# ~ # generate parameters
# ~ theta = scipy.linspace( 0.0, scipy.pi/2.0, 5 )
# ~ phi_t1 = scipy.linspace(0, scipy.pi*2.0, elements * 4 + 1)
# ~ phi_t2 = scipy.linspace(0, scipy.pi*2.0, elements * 3 + 1)
# ~ phi_t3 = scipy.linspace(0, scipy.pi*2.0, elements * 2 + 1)
# ~ phi_t4 = scipy.linspace(0, scipy.pi*2.0, elements * 1 + 1)
# ~ phi = [ phi_t1, phi_t2, phi_t3, phi_t4 ]
# ~ xparam = []
# ~ yparam = []
# ~ zparam = []
# ~
# ~ for i in range( len(theta) - 1 ):
# ~ z = radius * scipy.sin( theta[i] )
# ~ cos_theta = scipy.cos( theta[i] )
# ~ for p in range( len( phi[i] ) - 1 ):
# ~ xparam.append( [ radius * cos_theta * scipy.cos( phi[i][p] ) ] )
# ~ yparam.append( [ radius * cos_theta * scipy.sin( phi[i][p] ) ] )
# ~ zparam.append( [ z ] )
# ~
# ~ # apex
# ~ xparam.append( [0.0] )
# ~ yparam.append( [0.0] )
# ~ zparam.append( [ radius * scipy.sin( theta[-1] ) ])
# ~
# ~ return (hemi, xparam, yparam, zparam)
# ~ else:
# ~ return hemi
# ~
# ~
# ~ #=================#
# ~ # private mappers #
# ~ #=================#
# ~ def _reverse_custom_mapping( map ):
# ~ # map is a dictionary
# ~ max_key = max(map.keys())
# ~ map_r = {}
# ~ for k in map.keys():
# ~ new_key = max_key - map[k]
# ~
# ~ map_r[k] = new_key
# ~
# ~ return map_r
# ~
# ~ def _flip_3tri_ring_mapping( map, n_elements ):
# ~ # map is a dictionary
# ~ lim1 = n_elements * 2 - 1 # 3
# ~ lim2 = n_elements * 5 - 1 # 9
# ~ lim3 = n_elements * 9 - 1 # 17
# ~
# ~ max_key = max(map.keys())
# ~ map_r = {}
# ~ for k in map.keys():
# ~ new_key = max_key - map[k]
# ~ if new_key == lim1:
# ~ new_key = 0
# ~ elif new_key == lim2:
# ~ new_key = lim1 + 1
# ~ elif new_key == lim3:
# ~ new_key = lim2 + 1
# ~ else:
# ~ new_key += 1
# ~
# ~ map_r[k] = new_key
# ~
# ~ return map_r
# ~
# ~
# ~ def _tri_hemisphere_reverse_mapper( elements ):
# ~ # maps tri_hemispheres in reverse (starting from apex) for use in the
# ~ # sphere generator
# ~ remap = {}
# ~
# ~ # 1st elements
# ~ remap[0] = elements + 1
# ~ remap[1] = elements + 2
# ~ remap[2] = elements + 3
# ~ remap[5] = 1
# ~ remap[3] = 2
# ~ remap[4] = 0
# ~
# ~ # other elements
# ~ old_t1 = 6
# ~ old_t2 = 8
# ~ new_t1 = elements + 4
# ~ new_t2 = 3
# ~
# ~ for i in range( 1, elements ):
# ~ remap[old_t1] = new_t1
# ~ if i < elements - 1:
# ~ remap[old_t1 + 1] = new_t1 + 1
# ~ remap[old_t2] = new_t2
# ~
# ~ old_t1 += 3
# ~ new_t1 += 2
# ~ old_t2 += 3
# ~ new_t2 += 1
# ~
# ~ return remap
# ~
# ~ def _four_tri_patch_hemisphere_remapper( elements ):
# ~ remap = {}
# ~
# ~ old_t1 = 0
# ~ old_t2 = 5
# ~ old_t3 = 9
# ~ old_t4 = 12
# ~
# ~ new_t1 = 0
# ~ new_t2 = elements * 4
# ~ new_t3 = new_t2 + 3 * elements
# ~ new_t4 = new_t3 + 2 * elements
# ~
# ~ # 1st element
# ~ for old_t1 in range( old_t1, old_t1 + 5 ):
# ~ remap[old_t1] = new_t1
# ~ new_t1 += 1
# ~
# ~ for old_t2 in range( old_t2, old_t2 + 4 ):
# ~ remap[old_t2] = new_t2
# ~ new_t2 += 1
# ~
# ~ for old_t3 in range( old_t3, old_t3 + 3 ):
# ~ remap[old_t3] = new_t3
# ~ new_t3 += 1
# ~
# ~ for old_t4 in range( old_t4, old_t4 + 2 ):
# ~ remap[old_t4] = new_t4
# ~ new_t4 += 1
# ~
# ~ # other elements
# ~ for e in range( 1, elements ):
# ~ old_t1 = 15 + ( e - 1 ) * 10
# ~ old_t2 = 19 + ( e - 1 ) * 10
# ~ old_t3 = 22 + ( e - 1 ) * 10
# ~ old_t4 = 24 + ( e - 1 ) * 10
# ~
# ~ if e < elements - 1:
# ~ for old_t1 in range( old_t1, old_t1 + 4 ):
# ~ remap[old_t1] = new_t1
# ~ new_t1 += 1
# ~
# ~ for old_t2 in range( old_t2, old_t2 + 3 ):
# ~ remap[old_t2] = new_t2
# ~ new_t2 += 1
# ~
# ~ for old_t3 in range( old_t3, old_t3 + 2 ):
# ~ remap[old_t3] = new_t3
# ~ new_t3 += 1
# ~
# ~ for old_t4 in range( old_t4, old_t4 + 1 ):
# ~ remap[old_t4] = new_t4
# ~ new_t4 += 1
# ~ else:
# ~ for old_t1 in range( old_t1, old_t1 + 3 ):
# ~ remap[old_t1] = new_t1
# ~ new_t1 += 1
# ~
# ~ for old_t2 in range( old_t2 - 1, old_t2 + 1 ):
# ~ remap[old_t2] = new_t2
# ~ new_t2 += 1
# ~
# ~ for old_t3 in range( old_t3 - 2, old_t3 - 1 ):
# ~ remap[old_t3] = new_t3
# ~ new_t3 += 1
# ~
# ~ remap[14] = elements * 10
# ~
# ~ return remap
# ~
# ~ def _two_quad_ring_remapper( elements ):
# ~ remap = {}
# ~ old = scipy.array( [9, 11, 13] )
# ~ old1 = scipy.array( [0, 3, 6] )
# ~ new = scipy.array( [0, elements * 2, elements * 4] )
# ~
# ~ for e in range( 0, elements ):
# ~
# ~ # 1st element
# ~ if e == 0:
# ~ for i in range( 3 ):
# ~ remap[old1[0]] = new[0]
# ~ remap[old1[1]] = new[1]
# ~ remap[old1[2]] = new[2]
# ~ new += 1
# ~ old1 += 1
# ~
# ~ # last element
# ~ elif e == elements - 1:
# ~ oldn = old + (e-1) * 6
# ~ remap[oldn[0]] = new[0]
# ~ remap[oldn[0]+1] = new[1]
# ~ remap[oldn[0]+2] = new[2]
# ~
# ~ else:
# ~ oldn = old + (e-1) * 6
# ~ remap[oldn[0]] = new[0]
# ~ remap[oldn[1]] = new[1]
# ~ remap[oldn[2]] = new[2]
# ~
# ~ remap[oldn[0] + 1] = new[0] + 1
# ~ remap[oldn[1] + 1] = new[1] + 1
# ~ remap[oldn[2] + 1] = new[2] + 1
# ~
# ~ new += 2
# ~
# ~ return remap
# ~
# ~
# ~ #==================#
# ~ # femur primitives #
# ~ #==================#
# ~
# ~ # head and neck =======================================================#
# ~
# ~ def head_neck( h_r, n_l, n_r = None, h_adiv = 5, h_idiv = 3, h_maxi = 3.0/4.0*scipy.pi, n_ldiv = 1 ):
# ~ """ A truncated 3/4 sphere(head) and a truncated cone (neck)
# ~ origin at centre of sphere
# ~ parameters:  h_r = head radius
# ~ n_l = neck length
# ~ n_r = neck radius at the base of the cone, if
# ~ none, = radius at base of sphere
# ~ h_adiv = head sphere azimuth divisions
# ~ h_idiv = head sphere inclination divisions
# ~ h_maxi = head sphere maximum inclination
# ~ n_ldiv = neck cone longitudinal divisions
# ~ """
# ~
# ~ h_r = float(h_r)
# ~ n_l = float(n_l)
# ~
# ~ # initialise main field ensemble
# ~ dimensions = 2
# ~ basis_type = 'triangular_quadratic_2d'
# ~
# ~ head_neck = E.ensemble_field_function('head_neck', dimensions, debug = 0)
# ~ head_neck.set_basis( basis_type )
# ~ head_neck.set_new_mesh('head_neck_mesh')
# ~
# ~ xparams = []
# ~ yparams = []
# ~ zparams = []
# ~
# ~ # generate head sphere ============================================#
# ~ ( h_sphere, h_x, h_y, h_z, h_connect_points_n ) = sphere( h_adiv, h_idiv, h_r, h_maxi )
# ~
# ~ # ensemble points at the base of the sphere for connection and their radius
# ~ h_connect_points = range( h_connect_points_n )
# ~ h_connect_r = h_r * scipy.sin( h_maxi )
# ~ h_connect_z = h_z[0][0]
# ~
# ~ if not n_r:
# ~ n_r = h_connect_r
# ~ else:
# ~ n_r = float(n_r)
# ~
# ~ # generate truncated cone =========================================#
# ~ # generate radius values
# ~ n_r_range = scipy.linspace( n_r, h_connect_r, n_ldiv * 2 + 1)
# ~
# ~ # generate two_quad_rings from base of cone up
# ~ for i in range( n_ldiv ):
# ~ ring_r = n_r_range[ i*2: i*2 + 3 ]
# ~ (ring, ring_x, ring_y, ring_z) = two_quad_ring( h_connect_points_n/2, n_l/n_ldiv, ring_r )
# ~ head_neck.add_element( ring )
# ~
# ~ # keep 1st 2 rows of parameters values
# ~ xparams += ring_x[ 0:h_connect_points_n * 2 ]
# ~ yparams += ring_y[ 0:h_connect_points_n * 2 ]
# ~ # translate in z
# ~ z = ring_z[ 0:h_connect_points_n * 2 ]
# ~ z_shift = h_connect_z - n_l/n_ldiv * ( n_ldiv - i )
# ~ z = list( scipy.add( z, z_shift ) )
# ~ z = [ list(i) for i in z ]
# ~ zparams += z
# ~
# ~ # add sphere to mesh ==============================================#
# ~ head_neck.add_element( h_sphere )
# ~ xparams += h_x
# ~ yparams += h_y
# ~ zparams += h_z
# ~
# ~ # connect top of truncated cone to base of sphere
# ~ for element in range( n_ldiv ):
# ~ for i in range( h_connect_points_n ):
# ~ head_neck.connect_element_points( [(element, h_connect_points_n*n_ldiv*2 + i), (element + 1, i)] )
# ~
# ~ # map
# ~ head_neck.map_parameters()
# ~
# ~ connect_points = h_connect_points
# ~
# ~ return ( head_neck, xparams, yparams, zparams, connect_points )
# ~
# ~
# ~
# ~
# ~
# ~
# ~
# ~
