"""
FILE: meshnormaltest.py
LAST MODIFIED: 24-01-2017
DESCRIPTION: Tests the evaluation of derivative and normals on quadratic
    simplex meshes
    
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
"""

import numpy as np
from gias2.visualisation import fieldvi
from gias2.fieldwork.field import geometric_field
from gias2.fieldwork.field import geometric_field_fitter as GFF
from gias2.fieldwork.field import template_fields

from mayavi import mlab

#==========================================================================#
# main params
#=============================================================================#
normal_d = 4

# visualisation options
do_visualise = True  # whether to show 3D visualisation
vgfd = (8,8)  # mesh discretisation for rendering

# create mesh
eff, px, py, pz = template_fields.two_quad_ring(3, 10, 10)
gf = geometric_field.geometric_field(
        'test', 3, ensemble_field_function=eff
        )
gf.set_field_parameters(np.array([px, py, pz]))
gf.flatten_ensemble_field_function()

# make model evaluators
gf_eval = geometric_field.makeGeometricFieldEvaluatorSparse(gf, vgfd)

# Element boundary smoothness
smooth_normal_ = GFF.normalSmoother2(gf.ensemble_field_function.flatten()[0])
smooth_normal_obj = smooth_normal_.makeObj(normal_d)

# get EP positions of element boundarys
# the EPs of each shared edge
edge_normal_vecs = []
for edge_eps in smooth_normal_.edgeEvalPoints[:1]:
    elem1, xis1, elem2, xis2 = edge_eps

    # evaluate coordinates of points at which normals are evaluated
    ep1 = gf.evaluate_geometric_field_at_element_points(elem1, xis1).T
    ep2 = gf.evaluate_geometric_field_at_element_points(elem2, xis2).T

    # evaluate the normals of points at which normals are evaluated
    en1 = np.array([gf.evaluate_normal_at_element_point(elem1, xi[np.newaxis,]).T for xi in xis1])
    en2 = np.array([gf.evaluate_normal_at_element_point(elem2, xi[np.newaxis,]).T for xi in xis2])

    edge_normal_vecs.append([ep1, ep2, en1, en2])

# Points to show xi field
test_xis = np.array([
    [0.0,0.0], [0.2,0.0],[0.4,0.0],[0.6,0.0],[0.8,0.0],[1,0],
    [0.0,0.2], [0.2,0.2],[0.4,0.2],[0.6,0.2],[0.8,0.2],
    [0.0,0.4], [0.2,0.4],[0.4,0.4],[0.6,0.4],
    [0.0,0.6], [0.2,0.6],[0.4,0.6],
    [0.0,0.8], [0.2,0.8],
    [0.0,1.0],
    ])

test_xi_scalar = np.arange(test_xis.shape[0])
test_xi_X = gf.evaluate_geometric_field_at_element_points(0, test_xis)
test_xi_Xdxi1 = gf.evaluate_geometric_field_at_element_points(0, test_xis, (1,0))[:,1,:]
test_xi_Xdxi2 = gf.evaluate_geometric_field_at_element_points(0, test_xis, (0,1))[:,1,:]

#=============================================================================#
# visualise
#=============================================================================#

# draw edge normals
draw_edge = 0
E = edge_normal_vecs[draw_edge]
s = np.arange(E[0].shape[0])

gf.display_geometric_field([10,10])

mlab.quiver3d(
    E[0][:,0], E[0][:,1], E[0][:,2], 
    E[2][:,0], E[2][:,1], E[2][:,2],
    scalars=s, color=(1,0,0)
    )
mlab.quiver3d(
    E[1][:,0], E[1][:,1], E[1][:,2], 
    E[3][:,0], E[3][:,1], E[3][:,2],
    scalars=s, color=(0,0,1)
    )

mlab.quiver3d(
    test_xi_X[0], test_xi_X[1], test_xi_X[2], 
    test_xi_Xdxi1[0], test_xi_Xdxi1[1],test_xi_Xdxi1[2],
    scalars=test_xi_scalar, color=(0,0,1), name='dxi1'
    )
mlab.quiver3d(
    test_xi_X[0], test_xi_X[1], test_xi_X[2], 
    test_xi_Xdxi2[0], test_xi_Xdxi2[1],test_xi_Xdxi2[2],
    scalars=test_xi_scalar, color=(1,0,0), name='dxi2'
    )