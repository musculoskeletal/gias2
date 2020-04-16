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


def _loadSimpleMesh(filename):
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


def vrml2SimpleMesh(VRMLFilename):
    vrml2VTK = vtk.vtkVRMLImporter()
    vrml2VTK.SetFileName(VRMLFilename)
    vrml2VTK.Read()
    vrml2VTK.Update()

    actors = vrml2VTK.GetRenderer().GetActors()
    actors.InitTraversal()
    numberOfActors = actors.GetNumberOfItems()

    simpleMeshes = []

    for i in range(numberOfActors):
        polydata = actors.GetNextActor().GetMapper().GetInput()
        numberOfPoints = polydata.GetNumberOfPoints()
        numberOfCells = polydata.GetNumberOfCells()

        points = numpy.array([polydata.GetPoint(pi) for pi in range(numberOfPoints)])

        polys = polydata.GetPolys().GetData()
        polysSize = polys.GetSize()
        # ~ polysData = numpy.array([polys.GetValue(i) for i in xrange(polysSize)])
        polysData = numpy.array([polys.GetValue(i) for i in range(numberOfCells * 4)])
        tri = polysData.reshape((-1, 4))[:, 1:4]

        simpleMeshes.append(SimpleMesh(points, tri))

    return simpleMeshes


def stl2SimpleMesh(STLFilename):
    STLReader = vtk.vtkSTLReader()
    STLReader.SetFileName(STLFilename)
    STLReader.MergingOn()
    STLReader.Update()
    polydata = STLReader.GetOutput()
    numberOfPoints = polydata.GetNumberOfPoints()
    numberOfCells = polydata.GetNumberOfCells()

    points = numpy.array([polydata.GetPoint(pi) for pi in range(numberOfPoints)])

    polys = polydata.GetPolys().GetData()
    polysSize = polys.GetSize()
    # ~ polysData = numpy.array([polys.GetValue(i) for i in xrange(polysSize)])
    polysData = numpy.array([polys.GetValue(i) for i in range(numberOfCells * 4)])
    tri = polysData.reshape((-1, 4))[:, 1:4]

    sm = SimpleMesh(points, tri)

    return sm


class SimpleMesh(object):
    def __init__(self, v=None, f=None, H=None, K=None, k1=None, k2=None, E=None, data=None):
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

    def load(self, filename):
        s = _loadSimpleMesh(filename)
        self.__init__(*s)

    def save(self, filename):
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

    def exportINP(self, filename, name=None, preamble=None):
        elemType = 'R3D3'
        inpWriter = inp.InpWriter(filename, autoFormat=True, nodeOffset=1)

        if preamble == None:
            preamble = 'Exported from GIAS'

        inpWriter.addPreamble(preamble)

        if name == None:
            name = 'mesh'

        inpWriter.addMesh(name, elemType, self.v, self.f)
        inpWriter.write()

    def disp(self, curvature=None, figure=None, scalar=None, lim=[-0.2, 0.2]):
        if figure == None:
            fig = mlab.figure()
        else:
            fig = figure

        if scalar != None:
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

    def dispLabel(self, labels, figure=None):
        if figure == None:
            figure = mlab.figure()

        return mlab.triangular_mesh(self.v[:, 0], self.v[:, 1], self.v[:, 2], self.f, scalars=labels, figure=figure,
                                    vmax=labels.max(), vmin=labels.min())

    def setVerticesNeighbourhoods(self, r):
        """ gets the neighbourhood vertices and faces up to radius r for
        each vertex V. r is the number of vertices away from V.
        """

        log.debug('finding neighbourhoods of size {}'.format(r))
        self.neighbourhoodSize = r
        self.neighbourFaces = []
        self.neighbourVertices = []

        if not self.has1Ring:
            self.set1Ring()
        getNeighbour = self.makeNeighbourhoodGetter(r)

        for vi, V in enumerate(self.v):
            sys.stdout.write('\r' + str(vi))
            sys.stdout.flush()
            neighbourVertices, neighbourFaces = getNeighbour(vi)
            self.neighbourFaces.append(list(neighbourFaces))
            self.neighbourVertices.append(list(neighbourVertices))

        sys.stdout.write('\n')
        self.hasNeighbourhoods = 1
        return

    def set1Ring(self):
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
                    self.faces1Ring[v] = set([fi, ])

                try:
                    self.vertices1Ring[v] = self.vertices1Ring[v].union(f)
                except KeyError:
                    self.vertices1Ring[v] = set(f)

        for vi, f in list(self.vertices1Ring.items()):
            f.remove(vi)

        self.has1Ring = True

    def set1RingFaces(self):
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

    def _getAdjacent(self, vi, depth, vertexList, faceList):
        """ recursive gets the adjacent faces and vertices to vertex V.
        uses sets instead of lists
        """

        newVertices = set()
        # add adjacent faces of current vertex to faceList
        for fid in self.faces1Ring[vi]:
            faceList.add(fid)
            # add adjacent vertices to vertexList
            for vid in self.f[fid]:
                if vid not in vertexList:
                    newVertices.add(vid)

        vertexList = vertexList.union(newVertices)

        # recurse for new vertices
        if depth > 1:
            # ~ print 'recurse, depth =', depth-1
            for vid in newVertices:
                vertexList, faceList = self._getAdjacent(vid, depth - 1, vertexList, faceList)

        return vertexList, faceList

    def makeNeighbourhoodGetter(self, nRing):

        if not self.has1Ring:
            self.set1Ring()

        def getNeighbour(vI):
            vertexList = set([vI, ])
            faceList = set()
            neighbourVertices, neighbourFaces = self._getAdjacent(vI, nRing, vertexList, faceList)
            # first element of vertexList is the current vertex - remove
            neighbourVertices.remove(vI)
            return neighbourVertices, neighbourFaces

        return getNeighbour

    def calcFaceProperties(self):

        faceVertices = numpy.array([self.v[F] for F in self.f])

        v1 = faceVertices[:, 1, :] - faceVertices[:, 0, :]
        v2 = faceVertices[:, 2, :] - faceVertices[:, 0, :]
        v1v2 = numpy.cross(v1, v2)
        self.faceNormals = normalise2(v1v2)
        self.hasFaceNormals = True
        self.faceAreas = 0.5 * mag2(v1v2)
        self.faceBarycenters = (faceVertices[:, 0, :] + (faceVertices[:, 1, :] + faceVertices[:, 2, :])) / 3.0

    def calcVertexNormals(self, sigma, nsize=1, normalsout=True):
        """ calculate the normal at each vertex using normal voting. Considers
        all neighbouring vertices up to nsize edges away.
        """
        # self.saliencyCoeff = saliencyCoeff
        log.debug('calculating normals...')

        self.calcFaceProperties()
        if not self.has1Ring:
            self.set1Ring()
        if nsize == 1:
            allNeighFaces = self.faces1Ring
            allNeighVerts = self.vertices1Ring
        else:
            self.setVerticesNeighbourhoods(nsize)
            allNeighFaces = self.neighbourFaces
            allNeighVerts = self.neighbourVertices

        fBary = self.faceBarycenters
        fNormal = self.faceNormals
        fArea = self.faceAreas
        AMax = self.faceAreas.max()

        V = numpy.zeros((3, 3), dtype=float)
        self.vertexNormals = numpy.zeros((self.v.shape[0], 3), dtype=float)

        # for each vertex get neighbourhood faces 
        for vi, v in enumerate(self.v):
            # for each face i calculate normal vote Ni and weighting
            neighFaces = allNeighFaces[vi]
            nFaces = len(neighFaces)
            if not nFaces:
                raise RuntimeWarning('no faces: vertex').with_traceback(v.ID)

            fBaryV = numpy.array([fBary[f] for f in neighFaces])
            fNormalV = numpy.array([fNormal[i] for i in neighFaces])
            fAreaV = numpy.array([fArea[i] for i in neighFaces])

            # calc votes
            vc = normalise2(fBaryV - v)
            cosTheta = fNormalV[:, 0] * vc[:, 0] + fNormalV[:, 1] * vc[:, 1] + fNormalV[:, 2] * vc[:, 2]
            # ~ pdb.set_trace()
            NI = fNormalV - 2.0 * vc * cosTheta[:, numpy.newaxis]
            NI = numpy.where(numpy.isfinite(NI), NI, 0.0)

            # calc vote weights
            gV = mag2(fBaryV - v)
            WI = (fAreaV / AMax) * numpy.exp(-gV / sigma)

            # form covariance matrix V, and do eigendecomp 
            V[:, :] = 0.0
            WI = WI / WI.sum()  # normalise weights to sum to 1
            for i, n in enumerate(NI):
                V += WI[i] * numpy.kron(n, n[:, numpy.newaxis])

            try:
                l, e = eigh(V)
            except ValueError:
                log.debug('WARNING: singular V for vertex', vi)
                l = numpy.zeros(3)
                e = numpy.eye(3)
            else:
                l, e = _sortEigDesc(l, e)

            # print vi
            # print e
            # print self.vertexNormals.shape
            self.vertexNormals[vi, :] = e[:, 0]

            # # classify geometry at v
            # if self.saliencyCoeff == 0.0:
            #   v.surfaceType = 'surface'
            #   v.normal = vector( e[:,0] )
            #   self.surfaceVertices.append( v.id )
            # else:
            #   classifyVertex( v, l, e, self.saliencyCoeff )

            #   if v.surfaceType == 'surface':
            #       self.surfaceVertices.append( v.id )
            #   elif v.surfaceType == 'edge':
            #       self.edgeVertices.append( v.id )
            #   else:
            #       self.NPVertices.append( v.id )

        self.filterVertexNormals()
        self.hasVertexNormals = 1

        # make sure normals point out
        if normalsout:
            if not normals_is_out(self.v, self.vertexNormals):
                self.vertexNormals *= -1.0

            if not normals_is_out(self.faceBarycenters, self.faceNormals):
                self.faceNormals *= -1.0

        return

    def filterVertexNormals(self):
        """
        Orient vertex normals to be consistent
        """

        log.debug('filtering normals...')
        aligned = numpy.zeros(len(self.v), dtype=bool)
        # front = set([0])
        front = set([self.f.min(), ])
        aligned[self.f.min()] = True

        while front:
            # sys.stdout.write( '\rfront size: '+str(len(front))+' aligned size: '+str(aligned.sum()) )
            # sys.stdout.flush()

            v = front.pop()
            # get vertices immediately ahead of the front
            nvs = numpy.array([vid for vid in self.vertices1Ring[v] if not aligned[vid]], dtype=int)
            # dot product normals to find inverted neighbours
            d = self.vertexNormals[nvs].dot(self.vertexNormals[v])
            self.vertexNormals[nvs[d < 0.0]] *= -1.0
            aligned[nvs] = True
            front = front.union(nvs)

        return

    def calcBoundingBox(self):
        self.boundingBox = numpy.array([self.v.min(0), self.v.max(0)]).T
        return self.boundingBox

    def calcCoM(self):

        a = self.faceAreas
        x = self.faceBarycenters
        self.CoM = (x * a[:, numpy.newaxis]).sum(0) / sum(a)
        return self.CoM

    def calcNormCoM(self):

        box = self.boundingBox
        CoM = self.CoM
        x = self.faceBarycenters
        self.normCoM = numpy.array([(CoM[0] - box[0, 0]) / (box[0, 1] - box[0, 0]),
                                    (CoM[1] - box[1, 0]) / (box[1, 1] - box[1, 0]),
                                    (CoM[2] - box[2, 0]) / (box[2, 1] - box[2, 0]),
                                    ])

        return self.normCoM

    def calcPMoments(self):
        areas = self.faceAreas
        v = self.faceBarycenters - self.CoM

        I11 = ((v[:, 1] * v[:, 1] + v[:, 2] * v[:, 2]) * areas).sum()
        I22 = ((v[:, 0] * v[:, 0] + v[:, 2] * v[:, 2]) * areas).sum()
        I33 = ((v[:, 1] * v[:, 1] + v[:, 0] * v[:, 0]) * areas).sum()
        I12 = -(v[:, 0] * v[:, 1] * areas).sum()
        I13 = -(v[:, 0] * v[:, 2] * areas).sum()
        I23 = -(v[:, 1] * v[:, 2] * areas).sum()

        I = numpy.array([[I11, I12, I13], [I12, I22, I23], [I13, I23, I33]])
        self.I = I

        u, s, vh = svd(I)
        self.principalMoments = s.real[::-1]
        self.principalAxes = numpy.fliplr(u.real)

        # ~ print ' %(one)8.6f\n %(two)8.6f\n %(three)8.6f\n'\
        # ~ %{'one':self.principalAxes[2,0], 'two':self.principalAxes[0,1], 'three':self.principalAxes[1,2]}

        if self.principalAxes[2, 0] < 0.0:
            self.principalAxes[:, 0] *= -1.0
        if self.principalAxes[0, 1] < 0.0:
            self.principalAxes[:, 1] *= -1.0
        if self.principalAxes[1, 2] < 0.0:
            self.principalAxes[:, 2] *= -1.0

        return self.principalMoments, self.principalAxes

    def alignPAxes(self):
        """ rotate mesh to align pAxes with cartesian axes
        """
        pAxes = self.principalAxes
        CoM = self.CoM
        targetCoM = numpy.zeros(3)
        targetPAxes = numpy.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=float)
        T = alignment.calcAffine((CoM, pAxes), (targetCoM, targetPAxes))
        self.transformAffine(numpy.vstack((T, numpy.ones(4))))
        return T

    def transformAffine(self, t):
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

    def getBoundaryVertices(self):
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

    def getOrderedBoundaryVertices(self):
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


def mag(x):
    return numpy.sqrt((x * x).sum())


def normalise(x):
    return x / numpy.sqrt((x * x).sum())


def mag2(x):
    return numpy.sqrt((x * x).sum(1))


def normalise2(x):
    return x / numpy.sqrt((x * x).sum(1))[:, numpy.newaxis]


def _sortEigDesc(l, e):
    """ Sorts evalues and vectors in descending order.
    l is an array of eigenvalues correponding to the eigenvectors in
    the columns of e
    """
    lSortI = abs(l).argsort()[::-1]
    lSort = numpy.array([l[i] for i in lSortI])
    eSort = numpy.array([e[:, i] for i in lSortI]).T

    return lSort, eSort


def normals_is_out(x, xn):
    # check_point = numpy.array([1e6, 1e6, 1e6])
    check_point = x.max(0) * 10.0

    # find closest point to check_point
    c_i = numpy.argmin(((x - check_point) ** 2.0).sum(1))

    # vector to check point
    v_check_point = check_point - x[c_i]

    return numpy.dot(v_check_point, xn[c_i]) > 0.0
