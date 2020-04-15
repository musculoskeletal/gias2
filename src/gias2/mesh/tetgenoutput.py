"""
FILE: tetgenoutput.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: Class for reading output files from tetgen.

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import numpy as np

from gias2.mesh import simplemesh as smt


def _loadTetgenSMesh(filename):
    with open(filename, 'r') as f:

        # find facet list
        nFaces = None
        while nFaces == None:
            l = f.readline()
            if 'facet list' in l:
                l = f.readline()
                nFaces = int(l.split()[0])

        faces = []
        for i in range(nFaces):
            l = f.readline()
            faces.append([int(i) for i in l.split()[1:]])

    return faces


def _loadTetgenNodeFile(filename):
    data = np.loadtxt(filename, dtype=float, skiprows=1, comments='#')
    with open(filename, 'r') as f:
        l = f.readline()
        nNodes, dim, nAttr, boundary = [int(i) for i in l.split()]

    nodeNumbers = data[:, 0].astype(int)
    nodes = data[:, 1:1 + dim].astype(float)
    if nAttr:
        attr = data[:, 1 + dim:1 + dim + nAttr]
    else:
        attr = None
    if boundary:
        bm = data[:, -1]
    else:
        bm = None

    return nodes, nodeNumbers, attr, bm


def _loadTetgenEleFile(filename):
    data = np.loadtxt(filename, dtype=int, skiprows=1, comments='#')
    with open(filename, 'r') as f:
        l = f.readline()
        nElems, nNodesPerElem, region = [int(i) for i in l.split()]

    elemNumbers = data[:, 0]
    elems = data[:, 1:1 + nNodesPerElem]
    if data.shape[1] > (1 + nNodesPerElem):
        attr = data[:, 1 + nNodesPerElem:]
    else:
        attr = None

    return elems, elemNumbers, attr, region


def _loadTetgenFaceFile(filename):
    data = np.loadtxt(filename, dtype=int, skiprows=1, comments='#')
    with open(filename, 'r') as f:
        l = f.readline()
        nFaces, region = [int(i) for i in l.split()]

    faceNumbers = data[:, 0]
    faces = data[:, 1:4]
    bm = data[:, -1]

    return faces, faceNumbers, bm


class TetgenOutput(object):

    def __init__(self, filename=None):
        self.filename = filename

        self.nodes = None
        self.nodeNumbers = None
        self.nodeAttr = None
        self.nodeBM = None

        self.volElems = None
        self.volElemNumbers = None
        self.volElemAttr = None
        self.volElemRegion = None
        self.volElemCentroids = None

        self.faces = None
        self.faceNumbers = None
        self.faceBM = None

        self.surfElems = None

    def load(self, filename=None):
        if filename != None:
            self.filename = filename

        self.readNode()
        self.readEle()
        self.readFaces()
        self.readSMesh()

    def readNode(self):
        self.nodes, self.nodeNumbers, \
        self.nodeAttr, self.nodeBM = _loadTetgenNodeFile(self.filename + '.1.node')

    def readEle(self):
        self.volElems, self.volElemNumbers, \
        self.volElemAttr, self.volElemRegion = _loadTetgenEleFile(self.filename + '.1.ele')

    def readFaces(self):
        self.faces, self.faceNumbers, \
        self.faceBM = _loadTetgenFaceFile(self.filename + '.1.face')

    def readSMesh(self):
        self.surfElems = _loadTetgenSMesh(self.filename + '.1.smesh')

    def exportSimplemesh(self):
        # S = smt.SimpleMesh(v=self.nodes, f=self.surfElems)
        S = smt.SimpleMesh(v=self.nodes, f=self.faces - 1)
        return S

    def calcVolElemCentroids(self):
        volElemShape = self.volElems.shape
        nodes_exp = np.zeros([self.nodeNumbers.max() + 1, self.nodes.shape[1]], dtype=float)
        nodes_exp[self.nodeNumbers, :] = self.nodes
        volElemNodesFlat = np.ravel(self.volElems)
        volElemNodeCoordsFlat = nodes_exp[volElemNodesFlat]
        volElemNodeCoords = volElemNodeCoordsFlat.reshape([volElemShape[0], volElemShape[1], self.nodes.shape[1]])
        self.volElemCentroids = volElemNodeCoords.mean(1)

        return self.volElemCentroids

    def getSurfaceNodes(self):
        """ returns the node number and node coordinates of surface nodes
        """
        surfNodeInds = np.unique(np.ravel(self.surfElems))
        return surfNodeInds, self.nodes[surfNodeInds]
