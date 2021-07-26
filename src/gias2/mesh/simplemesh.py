"""
FILE: simplemesh.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: Classes and tools for working with triagulated meshes

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""
import logging
import shelve
from typing import List, Optional, Union, Tuple, Set, Callable

import numpy
import sys
import vtk
from numpy.linalg import svd, eigh

from gias2.common import transform3D
from gias2.mesh import inp
from gias2.registration import alignment_analytic as alignment

log = logging.getLogger(__name__)

try:
    from mayavi import mlab
except ImportError:
    log.debug('WARNING: Mayavi not installed, simpleMesh.disp will not work')


def _load_simple_mesh(filename: str):
    try:
        s = shelve.open(filename, 'r')
    except:
        raise IOError('unable to open ' + filename)

    attrList = ['vertices', 'faces', 'mean_curvature', 'gaussian_curvature', 'k1', 'k2', 'E', 'data']

    out = []
    for a in attrList:
        out.append(s.get(a, None))

    s.close()

    return out


def vrml_2_simple_mesh(vrml_filename: str) -> List['SimpleMesh']:
    """
    Read the meshes in a VRML file
    :param vrml_filename: filename of the VRML file
    :return: a list of SimpleMesh instances
    """
    vrml_2_vtk = vtk.vtkVRMLImporter()
    vrml_2_vtk.SetFileName(vrml_filename)
    vrml_2_vtk.Read()
    vrml_2_vtk.Update()

    actors = vrml_2_vtk.GetRenderer().GetActors()
    actors.InitTraversal()
    number_of_actors = actors.GetNumberOfItems()

    simple_meshes = []

    for i in range(number_of_actors):
        polydata = actors.GetNextActor().GetMapper().GetInput()
        number_of_points = polydata.GetNumberOfPoints()
        number_of_cells = polydata.GetNumberOfCells()

        points = numpy.array([polydata.GetPoint(pi) for pi in range(number_of_points)])

        polys = polydata.GetPolys().GetData()
        polys_data = numpy.array([polys.GetValue(i) for i in range(number_of_cells * 4)])
        tri = polys_data.reshape((-1, 4))[:, 1:4]

        simple_meshes.append(SimpleMesh(points, tri))

    return simple_meshes


def stl_2_simple_mesh(stl_filename: str) -> 'SimpleMesh':
    """
    Read an STL mesh
    :param stl_filename: filename of the stl file
    :return: a SimpleMesh representation of the stl mesh
    """
    stl_reader = vtk.vtkSTLReader()
    stl_reader.SetFileName(stl_filename)
    stl_reader.MergingOn()
    stl_reader.Update()
    polydata = stl_reader.GetOutput()
    number_of_points = polydata.GetNumberOfPoints()
    number_of_cells = polydata.GetNumberOfCells()

    points = numpy.array([polydata.GetPoint(pi) for pi in range(number_of_points)])

    polys = polydata.GetPolys().GetData()
    polys_data = numpy.array([polys.GetValue(i) for i in range(number_of_cells * 4)])
    tri = polys_data.reshape((-1, 4))[:, 1:4]

    sm = SimpleMesh(points, tri)

    return sm


class SimpleMesh(object):
    def __init__(
            self,
            v: Optional[Union[List[List[float]], numpy.ndarray]] = None,
            f: Optional[Union[List[List[int]], numpy.ndarray]] = None,
            H: Optional[Union[List[float], numpy.ndarray]] = None,
            K: Optional[Union[List[float], numpy.ndarray]] = None,
            k1: Optional[Union[List[List[float]], numpy.ndarray]] = None,
            k2: Optional[Union[List[List[float]], numpy.ndarray]] = None,
            E: Optional[Union[List[List[float]], numpy.ndarray]] = None,
            data: Optional[Union[List, numpy.ndarray]] = None):
        """
        A representation of a surface mesh
        :param v: list of vertex coordinates
        :param f: list of face vertex indices (3-tuples or lists)
        :param H: mean curvature of each vertex
        :param K: gaussian curvature of each vertex
        :param k1:
        :param k2:
        :param E:
        :param data: additional scalar or array values associated with each vertex
        """
        self.v = numpy.array(v)
        self.f = numpy.array(f)
        self.H = numpy.array(H)
        self.K = numpy.array(K)
        self.k1 = numpy.array(k1)
        self.k2 = numpy.array(k2)
        self.E = numpy.array(E)
        self.data = data
        self.faceNormals = None
        self.hasFaceNormals = False
        self.faceAreas = None
        self.faceBarycenters = None
        self.boundingBox = None
        self.normCoM = None
        self.principalMoments = None
        self.principalAxes = None

        self.has1Ring = False
        self.faces1Ring = None
        self.has1RingFaces = False
        self.faces1RingFaces = None
        self.hasNeighbourhoods = False
        self.neighbourhoodSize = None
        self.neighbourFaces = None
        self.neighbourVertices = None

        self.vertices1Ring = None
        self.boundaryVertexInd = None
        self.vertexNormals = None
        self.hasVertexNormals = False

        self.inertial_mat: Optional[numpy.ndarray] = None

    def load(self, filename: str) -> None:
        s = _load_simple_mesh(filename)
        self.__init__(*s)

    def save(self, filename: str) -> None:
        s = shelve.open(filename, protocol=2)
        s['vertices'] = self.v
        s['faces'] = self.f
        s['mean_curvature'] = self.H
        s['gaussian_curvature'] = self.K
        s['k1'] = self.k1
        s['k2'] = self.k2
        s['E'] = self.E
        try:
            s['data'] = self.data
        except AttributeError:
            pass

        s.close()
        return

    def exportINP(self, filename: str, name: Optional[str] = None, preamble: Optional[str] = None) -> None:
        elemType = 'R3D3'
        inpWriter = inp.InpWriter(filename, autoFormat=True, nodeOffset=1)

        if preamble is None:
            preamble = 'Exported from GIAS'

        inpWriter.addPreamble(preamble)

        if name is None:
            name = 'mesh'

        inpWriter.addMesh(name, elemType, self.v, self.f)
        inpWriter.write()

    def disp(self, curvature=None, figure=None, scalar=None, lim=(-0.2, 0.2)):
        if figure is not None:
            fig = mlab.figure()
        else:
            fig = figure

        if scalar is not None:
            return mlab.triangular_mesh(self.v[:, 0], self.v[:, 1], self.v[:, 2], self.f, scalars=scalar, figure=fig,
                                        vmax=lim[1], vmin=lim[0])
        elif curvature == 'H':
            return mlab.triangular_mesh(self.v[:, 0], self.v[:, 1], self.v[:, 2], self.f, scalars=self.H, figure=fig,
                                        vmax=lim[1], vmin=lim[0])
        elif curvature == 'K':
            return mlab.triangular_mesh(self.v[:, 0], self.v[:, 1], self.v[:, 2], self.f, scalars=self.K, figure=fig,
                                        vmax=lim[1], vmin=lim[0])
        else:
            return mlab.triangular_mesh(self.v[:, 0], self.v[:, 1], self.v[:, 2], self.f, figure=fig)

    def dispLabel(self, labels: numpy.ndarray, figure=None):
        if figure is None:
            figure = mlab.figure()

        return mlab.triangular_mesh(self.v[:, 0], self.v[:, 1], self.v[:, 2], self.f, scalars=labels, figure=figure,
                                    vmax=labels.max(), vmin=labels.min())

    def setVerticesNeighbourhoods(self, r: int) -> None:
        """ gets the neighbourhood vertices and faces up to radius r for
        each vertex V. r is the number of vertices away from V.

        Populates the attributes self.neighbourFaces and self.neighbourVertices that are lists that hold the sets of
        faces and vertices in the neighbourhood of each face or vertex, respectively, excluding the face or vertex
        itself.
        """

        log.debug('finding neighbourhoods of size {}'.format(r))
        self.neighbourhoodSize = r
        self.neighbourFaces = []
        self.neighbourVertices = []

        if not self.has1Ring:
            self.set1Ring()
        get_neighbour = self.makeNeighbourhoodGetter(r)

        for vi, V in enumerate(self.v):
            sys.stdout.write('\r' + str(vi))
            sys.stdout.flush()
            neighbour_vertices, neighbour_faces = get_neighbour(vi)
            self.neighbourFaces.append(list(neighbour_faces))
            self.neighbourVertices.append(list(neighbour_vertices))

        sys.stdout.write('\n')
        self.hasNeighbourhoods = 1
        return

    def set1Ring(self) -> None:
        """
        for each vertex, get the set of its neighbouring vertices
        and faces.
        """
        log.debug('setting 1-ring for vertices')
        self.faces1Ring = {}
        self.vertices1Ring = {}
        for fi, f in enumerate(self.f):
            for v in f:
                try:
                    self.faces1Ring[v].add(fi)
                except KeyError:
                    self.faces1Ring[v] = {fi}

                try:
                    self.vertices1Ring[v] = self.vertices1Ring[v].union(f)
                except KeyError:
                    self.vertices1Ring[v] = set(f)

        for vi, f in list(self.vertices1Ring.items()):
            f.remove(vi)

        self.has1Ring = True

    def set1RingFaces(self) -> None:
        """
        Create a dict of the adjacent faces of every face in sm
        """
        if not self.has1Ring:
            self.set1Ring()

        log.debug('setting 1-ring for faces')

        faces_1ring_faces = {}
        # share_edge_sets = [None, None, None]
        for fi, f in enumerate(self.f):
            # find 3 adj faces that share a side with f
            shared_edge_set_0 = set(self.faces1Ring[f[0]]).intersection(set(self.faces1Ring[f[1]])).difference([fi])
            shared_edge_set_1 = set(self.faces1Ring[f[0]]).intersection(set(self.faces1Ring[f[2]])).difference([fi])
            shared_edge_set_2 = set(self.faces1Ring[f[1]]).intersection(set(self.faces1Ring[f[2]])).difference([fi])
            faces_1ring_faces[fi] = shared_edge_set_0.union(shared_edge_set_1).union(shared_edge_set_2)

        self.faces1RingFaces = faces_1ring_faces

    def _getAdjacent(
            self,
            vi: int,
            depth: int,
            vertex_list: Set[int],
            face_list: Set[int]) -> Tuple[Set[int], Set[int]]:
        """ recursive gets the adjacent faces and vertices to vertex V.
        uses sets instead of lists
        """
        if self.faces1Ring.get(vi) is None:
            # vertex vi is unconnected or does not exist
            log.info('vertex %d is unconnected or does not exist in mesh', vi)
        else:
            new_vertices = set()
            # add adjacent faces of current vertex to faceList
            for fid in self.faces1Ring[vi]:
                face_list.add(fid)
                # add adjacent vertices to vertexList
                for vid in self.f[fid]:
                    if vid not in vertex_list:
                        new_vertices.add(vid)

            vertex_list = vertex_list.union(new_vertices)

            # recurse for new vertices
            if depth > 1:
                for vid in new_vertices:
                    vertex_list, face_list = self._getAdjacent(vid, depth - 1, vertex_list, face_list)

        return vertex_list, face_list

    def makeNeighbourhoodGetter(self, n_ring: int) -> Callable:

        if not self.has1Ring:
            self.set1Ring()

        def get_neighbour(v_i):
            vertex_list = {v_i}
            face_list = set()
            neighbour_vertices, neighbour_faces = self._getAdjacent(v_i, n_ring, vertex_list, face_list)
            # first element of vertexList is the current vertex - remove
            neighbour_vertices.remove(v_i)
            return neighbour_vertices, neighbour_faces

        return get_neighbour

    def calcFaceProperties(self) -> None:

        face_vertices = numpy.array([self.v[F] for F in self.f])

        v1 = face_vertices[:, 1, :] - face_vertices[:, 0, :]
        v2 = face_vertices[:, 2, :] - face_vertices[:, 0, :]
        v1v2 = numpy.cross(v1, v2)
        self.faceNormals = normalise2(v1v2)
        self.hasFaceNormals = True
        self.faceAreas = 0.5 * mag2(v1v2)
        self.faceBarycenters = (face_vertices[:, 0, :] + (face_vertices[:, 1, :] + face_vertices[:, 2, :])) / 3.0

    def calcVertexNormals(self, sigma: float, nsize: int = 1, normalsout: bool = True) -> None:
        """ calculate the normal at each vertex using normal voting. Considers
        all neighbouring vertices up to nsize edges away.
        """
        log.debug('calculating normals...')

        self.calcFaceProperties()
        if not self.has1Ring:
            self.set1Ring()
        if nsize == 1:
            all_neigh_faces = self.faces1Ring
        else:
            self.setVerticesNeighbourhoods(nsize)
            all_neigh_faces = self.neighbourFaces

        f_bary = self.faceBarycenters
        f_normal = self.faceNormals
        f_area = self.faceAreas
        a_max = self.faceAreas.max()

        v_mat = numpy.zeros((3, 3), dtype=float)
        self.vertexNormals = numpy.zeros((self.v.shape[0], 3), dtype=float)

        # for each vertex get neighbourhood faces 
        for vi, v in enumerate(self.v):
            # for each face i calculate normal vote Ni and weighting
            neigh_faces = all_neigh_faces[vi]
            n_faces = len(neigh_faces)
            if not n_faces:
                raise RuntimeWarning('no faces: vertex').with_traceback(v.ID)

            f_bary_v = numpy.array([f_bary[f] for f in neigh_faces])
            f_normal_v = numpy.array([f_normal[i] for i in neigh_faces])
            f_area_v = numpy.array([f_area[i] for i in neigh_faces])

            # calc votes
            vc = normalise2(f_bary_v - v)
            cos_theta = f_normal_v[:, 0] * vc[:, 0] + f_normal_v[:, 1] * vc[:, 1] + f_normal_v[:, 2] * vc[:, 2]
            normal_ind = f_normal_v - 2.0 * vc * cos_theta[:, numpy.newaxis]
            normal_ind = numpy.where(numpy.isfinite(normal_ind), normal_ind, 0.0)

            # calc vote weights
            g_v = mag2(f_bary_v - v)
            w_i = (f_area_v / a_max) * numpy.exp(-g_v / sigma)

            # form covariance matrix V, and do eigendecomp 
            v_mat[:, :] = 0.0
            w_i = w_i / w_i.sum()  # normalise weights to sum to 1
            for i, n in enumerate(normal_ind):
                v_mat += w_i[i] * numpy.kron(n, n[:, numpy.newaxis])

            try:
                l, e = eigh(v_mat)
            except ValueError:
                log.debug('WARNING: singular V for vertex', vi)
                e = numpy.eye(3)
            else:
                l, e = _sortEigDesc(l, e)

            self.vertexNormals[vi, :] = e[:, 0]

        self.filterVertexNormals()
        self.hasVertexNormals = 1

        # make sure normals point out
        if normalsout:
            if not normals_is_out(self.v, self.vertexNormals):
                self.vertexNormals *= -1.0

            if not normals_is_out(self.faceBarycenters, self.faceNormals):
                self.faceNormals *= -1.0

        return

    def filterVertexNormals(self) -> None:
        """
        Orient vertex normals to be consistent
        """

        log.debug('filtering normals...')
        aligned = numpy.zeros(len(self.v), dtype=bool)
        front = {self.f.min()}
        aligned[self.f.min()] = True

        while front:
            v = front.pop()
            # get vertices immediately ahead of the front
            nvs = numpy.array([vid for vid in self.vertices1Ring[v] if not aligned[vid]], dtype=int)
            # dot product normals to find inverted neighbours
            d = self.vertexNormals[nvs].dot(self.vertexNormals[v])
            self.vertexNormals[nvs[d < 0.0]] *= -1.0
            aligned[nvs] = True
            front = front.union(nvs)

        return

    def calcBoundingBox(self) -> numpy.ndarray:
        self.boundingBox = numpy.array([self.v.min(0), self.v.max(0)]).T
        return self.boundingBox

    def calcCoM(self) -> numpy.ndarray:

        a = self.faceAreas
        x = self.faceBarycenters
        self.CoM = (x * a[:, numpy.newaxis]).sum(0) / sum(a)
        return self.CoM

    def calcNormCoM(self) -> numpy.ndarray:

        box = self.boundingBox
        com = self.CoM
        x = self.faceBarycenters
        self.normCoM = numpy.array([(com[0] - box[0, 0]) / (box[0, 1] - box[0, 0]),
                                    (com[1] - box[1, 0]) / (box[1, 1] - box[1, 0]),
                                    (com[2] - box[2, 0]) / (box[2, 1] - box[2, 0]),
                                    ])

        return self.normCoM

    def calcPMoments(self) -> Tuple[numpy.ndarray, numpy.ndarray]:
        areas = self.faceAreas
        v = self.faceBarycenters - self.CoM

        I11 = ((v[:, 1] * v[:, 1] + v[:, 2] * v[:, 2]) * areas).sum()
        I22 = ((v[:, 0] * v[:, 0] + v[:, 2] * v[:, 2]) * areas).sum()
        I33 = ((v[:, 1] * v[:, 1] + v[:, 0] * v[:, 0]) * areas).sum()
        I12 = -(v[:, 0] * v[:, 1] * areas).sum()
        I13 = -(v[:, 0] * v[:, 2] * areas).sum()
        I23 = -(v[:, 1] * v[:, 2] * areas).sum()

        inertial_mat = numpy.array([[I11, I12, I13], [I12, I22, I23], [I13, I23, I33]])
        self.inertial_mat = inertial_mat

        u, s, vh = svd(inertial_mat)
        self.principalMoments = s.real[::-1]
        self.principalAxes = numpy.fliplr(u.real)

        if self.principalAxes[2, 0] < 0.0:
            self.principalAxes[:, 0] *= -1.0
        if self.principalAxes[0, 1] < 0.0:
            self.principalAxes[:, 1] *= -1.0
        if self.principalAxes[1, 2] < 0.0:
            self.principalAxes[:, 2] *= -1.0

        return self.principalMoments, self.principalAxes

    def alignPAxes(self) -> numpy.ndarray:
        """ rotate mesh to align pAxes with cartesian axes
        """
        p_axes = self.principalAxes
        com = self.CoM
        target_com = numpy.zeros(3)
        target_p_axes = numpy.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=float)
        mat = alignment.calcAffine((com, p_axes), (target_com, target_p_axes))
        self.transformAffine(numpy.vstack((mat, numpy.ones(4))))
        return mat

    def transformAffine(self, t: numpy.ndarray) -> None:
        """ transform mesh vertices by an affine
        transformation matrix T (shape = (3,4))
        """
        # newV = numpy.dot( t, numpy.vstack( (self.v.T, numpy.ones(self.v.shape[0])) ) )[:3,:].T 
        # self.v = newV

        self.v = transform3D.transformAffine(self.v, t)
        if self.vertexNormals is not None:
            self.vertexNormals = numpy.dot(t[:3, :3], self.vertexNormals.T).T
        if self.faceNormals is not None:
            self.faceNormals = numpy.dot(t[:3, :3], self.faceNormals.T).T
        if self.faceBarycenters is not None:
            self.faceBarycenters = transform3D.transformAffine(self.faceBarycenters, t)

    def getBoundaryVertices(self) -> Tuple[List[int], numpy.ndarray]:
        """ 
        Returns the indices and coordinates of vertices on the
        boundary or boundaries of the mesh. Boundary vertices have 
        and unequal number of 1-ring vertices to 1-ring faces.
        """

        if self.boundaryVertexInd != None:
            return self.boundaryVertexInd, self.v[self.boundaryVertexInd]
        else:
            self.set1Ring()

            log.debug('finding boundary vertices')
            boundaryVertexInd = []
            for vi in range(len(self.v)):
                try:
                    if len(self.vertices1Ring[vi]) != len(self.faces1Ring[vi]):
                        boundaryVertexInd.append(vi)
                except KeyError:
                    log.debug("WARNING: no neighbours for vertex", vi)
                    pass

            self.boundaryVertexInd = boundaryVertexInd
            return self.boundaryVertexInd, self.v[self.boundaryVertexInd]

    def getOrderedBoundaryVertices(self) -> List[List[int]]:
        """
        Returns lists of ordered boundary vertex indices. Each list contains the 
        boundary vertices of a boundary on the mesh.
        """

        bv = self.getBoundaryVertices()[0]
        bv_set = set(bv)
        boundaries = []
        while bv_set:
            ordered_bv = [bv_set.pop(), ]
            search_boundary = 1
            # starting from 1st bv, find next bv in its 1-ring
            while search_boundary:
                # get neighbour vertices
                neighv = self.vertices1Ring[ordered_bv[-1]]
                # look through each neighbour
                search_neigh = 1
                prev_bv = ordered_bv[-1]
                for nvi in neighv:
                    # if neighbour is a boundary vertex and isn't the previous one
                    if (nvi in bv_set) and (nvi != prev_bv):
                        ordered_bv.append(nvi)
                        bv_set.remove(nvi)
                        prev_bv = nvi
                        # stop searching neighbours
                        search_neigh = 0
                        break

                # if one wasn't found, terminate the search for this boundary
                if search_neigh == 1:
                    search_boundary = 0

            boundaries.append(ordered_bv)

        return boundaries


def mag(x: numpy.ndarray) -> float:
    return numpy.sqrt((x * x).sum())


def normalise(x: numpy.ndarray) -> numpy.ndarray:
    return x / numpy.sqrt((x * x).sum())


def mag2(x: numpy.ndarray) -> numpy.ndarray:
    return numpy.sqrt((x * x).sum(1))


def normalise2(x: numpy.ndarray) -> numpy.ndarray:
    return x / numpy.sqrt((x * x).sum(1))[:, numpy.newaxis]


def _sortEigDesc(l: numpy.ndarray, e: numpy.ndarray) -> Tuple[numpy.ndarray, numpy.ndarray]:
    """ Sorts evalues and vectors in descending order.
    l is an array of eigenvalues correponding to the eigenvectors in
    the columns of e
    """
    lSortI = abs(l).argsort()[::-1]
    lSort = numpy.array([l[i] for i in lSortI])
    eSort = numpy.array([e[:, i] for i in lSortI]).T

    return lSort, eSort


def normals_is_out(x: numpy.ndarray, xn: numpy.ndarray) -> numpy.ndarray:
    """
    Given a list of vertices and their normals, determine whether each normal is pointing in or out
    :param x: nx3 array of vertex coordinates
    :param xn: nx3 array of vertex normals
    :return: a boolean array, True is pointing out
    """
    check_point = x.max(0) * 10.0

    # find closest point to check_point
    c_i = numpy.argmin(((x - check_point) ** 2.0).sum(1))

    # vector to check point
    v_check_point = check_point - x[c_i]

    return numpy.dot(v_check_point, xn[c_i]) > 0.0
