"""
FILE: triangulate.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: functions for triangulating fieldwork meshes.

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""
import logging

from numpy import array, sort, where, unique, zeros, arange
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist, squareform
from scipy.special import comb

log = logging.getLogger(__name__)

try:
    from mayavi import mlab
except ImportError:
    log.debug('WARNING: Mayavi not imported, no 3-D visualisation')


class mesh_triangulator(object):

    def __init__(self, field, params, equi=True):
        self.f = field
        self.params = params
        self.equi = equi
        self.tri = None
        self.P = None

    def draw_surface_simple(self, d, colour='jet', scalar=None, figure=None, name=None, limits=(-0.2, 0.2)):
        """ just creates and returns the triangular_mesh object
        """
        # ~ self.tri = []            # list of triangulated vertices
        self.P = []  # list of [x,y,z] evaluated at each vertex

        # evaluate mesh at triangle vertices
        # evaluate field at the desired element descritisation
        # P[0] = x, P[1] = y, P[2] = z
        for p in self.params:
            self.P.append(self.f.evaluate_field_in_mesh(d, p))
        self.P = array(self.P)

        # triangulate vertices
        T = self._triangulate(d)

        if figure == None:
            figure = mlab.figure()

        # make mesh
        if figure == None:
            if name:
                return mlab.triangular_mesh(self.P[0], self.P[1], self.P[2], T, colormap=colour, scalars=scalar,
                                            name=name, vmin=limits[0], vmax=limits[1])
            else:
                return mlab.triangular_mesh(self.P[0], self.P[1], self.P[2], T, colormap=colour, scalars=scalar,
                                            vmin=limits[0], vmax=limits[1])
        else:
            if name:
                return mlab.triangular_mesh(self.P[0], self.P[1], self.P[2], T, colormap=colour, scalars=scalar,
                                            figure=figure, name=name, vmin=limits[0], vmax=limits[1])
            else:
                return mlab.triangular_mesh(self.P[0], self.P[1], self.P[2], T, colormap=colour, scalars=scalar,
                                            figure=figure, vmin=limits[0], vmax=limits[1])

    def draw_surface(self, d, plot_nodes=True, node_scale=1.0, plot_ep=False, number_ep=False):
        """ draws triangulated surface at discretisation d
        """
        # ~ self.tri = []            # list of triangulated vertices
        self.P = []  # list of [x,y,z] evaluated at each vertex

        # evaluate mesh at triangle vertices
        # evaluate field at the desired element descritisation
        # P[0] = x, P[1] = y, P[2] = z
        for p in self.params:
            self.P.append(self.f.evaluate_field_in_mesh(d, p))
        self.P = array(self.P)

        # triangulate vertices
        T = self._triangulate(d)

        # display surface
        mlab.figure()
        mlab.triangular_mesh(self.P[0], self.P[1], self.P[2], T)
        mlab.outline()
        mlab.xlabel('x')
        # plot nodes
        if plot_nodes:
            mlab.points3d(array(self.params[0])[:, 0], array(self.params[1])[:, 0], array(self.params[2])[:, 0],
                          mode='sphere', scale_factor=node_scale, color=(1, 0, 0))
            for i in range(len(self.params[0])):
                mlab.text(array(self.params[0])[i, 0], array(self.params[1])[i, 0], str(i),
                          z=array(self.params[2])[i, 0], line_width=0.001, width=0.01)
        # plot discretisation points
        if plot_ep:
            mlab.points3d(self.P[0], self.P[1], self.P[2], scale_factor=0.01)
        if number_ep:
            for i in range(self.P.shape[1]):
                mlab.text(self.P[0, i], self.P[1, i], str(i), z=self.P[2, i], line_width=0.001, width=0.01)

        mlab.show()
        return

    def draw_derivatives(self, d, deriv, plot_nodes=True, plot_ep=False, number_ep=False):
        """ deriv = ( x, i ), x=0, 1, or 2 for global position,
        i=xi1 or xi2, ie dx/di
        """
        pass

    def _triangulate(self, d):
        E = self.f.mesh.get_true_elements()
        self.tri = triangulate(E, d)
        return self.tri

    def _mergePoints(self, P):
        self.tri = filterDuplicateVertices(P, self.tri)
        return self.tri

    def _mergePoints2(self, P):
        P, self.tri, remainingVertexI, vertMap = filterDuplicateVertices2(P, self.tri)
        return P, self.tri, remainingVertexI, vertMap


def triangulate(E, d):
    """ forms 3-tuples defining triangles for each element in E 
    evaluated at a descritization of d
    """

    tri = []
    # triangulate
    # ~ e_n = len(self.f.mesh.elements)  # number of elements
    # ~ e_n = f.mesh.get_number_of_true_elements()
    i = 0  # point counter

    # ~ # for each element
    # ~ for e in range(e_n):
    # ~ i_row = d
    # ~ while i_row > 2:
    # ~ # for each row of element points
    # ~ for j in range( i_row-2 ):
    # ~ self.tri.append( (i, i+1, i+i_row) )
    # ~ self.tri.append( (i+1, i+i_row, i+i_row+1) )
    # ~ i += 1
    # ~ # last triangle of the row
    # ~ self.tri.append( (i, i+1, i+i_row) )
    # ~ i_row -= 1
    # ~ i += 2
    # ~ # last triangle of the element
    # ~ self.tri.append(  (i, i+1, i+i_row) )
    # ~ i += 3

    # ==============================================================#
    # ~ E = _get_field_true_elements( self.f.mesh )
    # ~ E = f.mesh.get_true_elements()
    for e in E:
        if 'quad' in e.type:
            for row in range(d[0] - 1):
                # for each row of element points
                for col in range(d[1] - 1):
                    tri.append((i, i + 1, i + d[1]))  # clockwise
                    # ~ tri.append( (i+1, i+d[1], i+d[1]+1) )    # anti-clockwise
                    tri.append((i + 1, i + d[1] + 1, i + d[1]))  # clockwise
                    i += 1
                i += 1
            i += d[1]

        elif 'tri' in e.type:
            i_col = d[0]
            i_row = d[1]
            while i_col > 2:
                # for each row of element points
                for j in range(i_col - 2):
                    tri.append((i, i + 1, i + i_row))  # clockwise
                    # ~ tri.append( (i+1, i+i_row, i+i_row+1) )  # anti-clockwise
                    tri.append((i + 1, i + i_row + 1, i + i_row))  # clockwise
                    i += 1
                # last triangle of the row
                tri.append((i, i + 1, i + i_row))
                i_row -= 1
                i_col -= 1
                i += 2
            # last triangle of the element
            tri.append((i, i + 1, i + i_row))
            i += 3

    return array(tri)


def findDuplicatePoints(P):
    """
    returns list of indices of points that are duplicates in coordinates
    """

    tree = cKDTree(P)
    d, Pi = tree.query(list(P), 2)
    d = array(d)
    PiDup = where(d.sum(1) < 1e-3)[0]  # duplicates have at least 2 closest points at distance 0.0 (one is itself)
    found = set()
    duplicatedGroups = []

    for i in PiDup:
        if i not in found:
            match = P == P[i]
            d = list(where(match.sum(1) == 3)[0])
            # ~ pdb.set_trace()
            duplicatedGroups.append(d)
            found = found.union(d)

    return duplicatedGroups


def findDuplicatePointsTree(P):
    """
    returns list of indices of points that are duplicates in coordinates
    also returns an array of length P.shape[0] that contains the new
    index for each point in P.
    """
    r = 1e-6
    tree = cKDTree(P)
    closest_points = tree.query_ball_tree(tree, r, 2)
    nPoints = P.shape[0]
    skip = set()
    dupGroups = []
    dupMap = zeros(P.shape[0], dtype=int)
    uniqueCount = 0

    for i in range(nPoints):
        if i not in skip:
            # find points too close to point i
            dupi = closest_points[i]

            skip = skip.union(dupi)

            # put i first in the list
            dupi[dupi.index(i)] = dupi[0]
            dupi[0] = i

            if len(dupi) > 1:
                dupGroups.append(dupi)
                dupMap[dupi] = uniqueCount
            else:
                dupMap[i] = uniqueCount

            uniqueCount += 1

    log.debug(uniqueCount, 'unique points')

    return dupGroups, dupMap


def findDuplicatePoints2(P):
    """
    returns list of indices of points that are duplicates in coordinates
    also returns an array of length P.shape[0] that contains the new
    index for each point in P.
    """
    dist = squareform(pdist(P, 'euclidean'))
    nPoints = P.shape[0]
    skip = set()
    dupGroups = []
    dupMap = zeros(P.shape[0], dtype=int)
    uniqueCount = 0

    for i in range(nPoints):
        if i not in skip:
            # find points too close to point i
            dupi = list(where(dist[i, :] < 1e-6)[0])

            skip = skip.union(dupi)

            # put i first in the list
            dupi[dupi.index(i)] = dupi[0]
            dupi[0] = i

            if len(dupi) > 1:
                dupGroups.append(dupi)
                dupMap[dupi] = uniqueCount
            else:
                dupMap[i] = uniqueCount

            uniqueCount += 1

    log.debug(uniqueCount, 'unique points')

    return dupGroups, dupMap


def _vectorFormIndices(n, i, j):
    return comb(n, 2) - comb(n - i, 2) + (j - i - 1)


def findDuplicatePoints3(P):
    """
    returns list of indices of points that are duplicates in coordinates
    also returns a list of indices mapping the each original to a 
    unique-fied list of points. 
    
    does not use squareform
    """
    dist = pdist(P, 'euclidean')
    nPoints = P.shape[0]
    skip = zeros(nPoints, dtype=bool)
    dupGroups = {}
    dupMap = zeros(nPoints, dtype=int)
    uniqueCount = 0

    for i in range(nPoints - 1):
        if not skip[i]:
            # get distances to point i
            j = arange(i + 1, nPoints)
            vecI = _vectorFormIndices(nPoints, i, j).astype(int)
            d = dist[vecI]
            # find points too close to point i
            dupi = list(j[where(d < 1e-6)])

            if len(dupi) > 0:
                # if duplicates
                skip[dupi] = True
                dupMap[dupi] = uniqueCount
                dupGroups[i] = dupi

            skip[i] = True
            dupMap[i] = uniqueCount
            uniqueCount += 1

    # last row of distance matrix does not exist
    # so if last point is not a duplicate, put in a unique entry
    if not skip[nPoints - 1]:
        dupMap[nPoints - 1] = uniqueCount
        skip[nPoints - 1] = True
        uniqueCount += 1

    log.debug(uniqueCount, 'unique points')

    return dupGroups, dupMap


def filterDuplicateVertices(V, T):
    TNew = array(T)
    # find duplicate vertices
    duplicates = findDuplicatePoints(V)

    # for each duplicate vertex v
    for D in duplicates:
        v = D[0]
        vd = D[1:]

        # find the triangles of its duplicate vd, replace vd by v
        for vdi in vd:
            TNew[where(TNew == vdi)] = v

    return TNew


def filterDuplicateVertices2(V, T):
    TNew = array(T, dtype=int)
    # find duplicate vertices
    # duplicates, dupMap = findDuplicatePoints2( V )
    duplicates, dupMap = findDuplicatePointsTree(V)

    # for each duplicate vertex v
    for D in duplicates:

        v = D[0]
        vd = D[1:]

        # find the triangles of its duplicate vd, replace vd by v
        for vdi in vd:
            TNew[where(TNew == vdi)] = v

    # re order points and triangulation
    usedVertexI = unique(TNew)

    TReorder = zeros(TNew.shape, dtype=int)
    VReorder = V[usedVertexI, :]

    for vertexI in usedVertexI:
        TReorder[where(TNew == vertexI)] = where(usedVertexI == vertexI)[0][0]

    return VReorder, TReorder, usedVertexI, dupMap


def filterDuplicateVertices3(V, T):
    TNew = array(T, dtype=int)
    # find duplicate vertices
    duplicates, dupMap = findDuplicatePoints3(V)

    # for each duplicate vertex v
    for v, vd in list(duplicates.items()):

        # for findDuplicatePoints2
        # ~ v = D[0]
        # ~ vd = D[1:]

        # for findDuplicatePoints3
        # ~ v = D[-1]
        # ~ vd = D[:-1]

        # find the triangles of its duplicate vd, replace vd by v
        for vdi in vd:
            TNew[where(TNew == vdi)] = v

    # re order points and triangulation
    # ~ usedVertexI = unique(TNew)
    usedVertexI = unique(dupMap)

    TReorder = zeros(TNew.shape, dtype=int)
    VReorder = V[usedVertexI, :]

    for vertexI in usedVertexI:
        TReorder[where(TNew == vertexI)] = where(usedVertexI == vertexI)[0][0]

    # alternative method, not sure why this doesnt work
    # ~ usedVertexI = unique(dupMap)
    # ~ VReorder = V[ usedVertexI ]
    # ~ TReorder = array( [dupMap[t] for t in T] )

    return VReorder, TReorder, usedVertexI, dupMap


def _get_field_true_elements(mesh):
    """ gets a list of pointers to all true elements in a mesh, in
    the other that they are evaluated
    """
    elem_list = []
    for element_number in sort(list(mesh.elements.keys())):
        # get element
        element = mesh.elements[element_number]
        if element.is_element:
            elem_list.append(element)
        else:
            elem_list += _get_field_true_elements(element)
    return elem_list
