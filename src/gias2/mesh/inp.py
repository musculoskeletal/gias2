"""
FILE: inp.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION:
classes and functions for reading and writing .inp files for
ABAQUS and FEBio

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""
import logging
from typing import List, Optional, Dict, Tuple, Union

import numpy as np

log = logging.getLogger(__name__)


COMMENTCHARS = '**'
ELEMNODES = {'C3D8R': 8,
             'R3D3': 3,
             'R3D4': 4,
             'C3D4': 4,
             'T3D2': 2,
             'S3': 3,
             'S4': 4,
             }


class Mesh(object):
    """ ABAQUS INP Mesh object
    """

    def __init__(self, name: str):
        self.name: str = name
        self.nodes: Optional[np.ndarray] = None
        self.nodeNumbers: Optional[np.ndarray] = None
        self.elems: Optional[List[List[int]]] = None
        self.elemNumbers: Optional[List[int]] = None
        self.elemType: Optional[str] = None
        self.elsets: dict = {}
        self.surfaces: dict = {}
        self._nodesDict: Dict[int, np.ndarray] = {}

    def getName(self):
        return self.name

    def setNodes(self, nodes: List, node_numbers: Union[np.ndarray, List[int]]) -> None:
        """ Set nodes of the mesh.
        Arguments:
        nodes : a list containing lists of node coordinates
        nodeNumbers : a list of node numbers corresponding to their
                      coordinate.
        """
        self.nodes = np.array(nodes)
        self.nodeNumbers = np.array(node_numbers)
        self._nodesDict = dict(zip(self.nodeNumbers, self.nodes))

    def getNode(self, node_number: int) -> np.ndarray:
        """Returns the coordinates of the node with node number
        nodeNumber
        """
        return self._nodesDict[node_number]

    def getNodes(self) -> np.ndarray:
        """Returns a list of all node coordinates
        """
        return self.nodes

    def getMeshNodes(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Return the list of nodes and node numbers that are actually
        referenced by the mesh elements
        """

        mesh_node_nums = np.unique(np.hstack(self.elems))
        mesh_node_coords = np.array([self._nodesDict[i] for i in mesh_node_nums])
        return mesh_node_nums, mesh_node_coords

    def getNumberOfNodes(self) -> int:
        """Returns the total number of nodes
        """
        return len(self.nodes)

    def setElems(self, elems: List[List], elem_numbers: List[int], elem_type: str) -> None:
        """Set elements of the mesh.
        Arguments:
        elems: a list containing the lists of the node numbers of each 
                element
        elem_numbers: a list of element numbers corresponding to its node lists
        elem_type: a string of the ABAQUS element type
        """
        self.elems = elems
        self.elemNumbers = elem_numbers
        self.elemType = elem_type

    def getElem(self, elem_number: int) -> List[int]:
        """Returns the node numbers of the element with element number
        elemNumber
        """
        return self.elems[self.elemNumbers.index(elem_number)]

    def getElems(self) -> List[List[int]]:
        """Returns a list of all elements' node numbers
        """
        return self.elems

    def getElemType(self) -> str:
        """Returns the element type string
        """
        return self.elemType

    def getNumberOfElems(self) -> int:
        """Returns the total number of elements
        """
        return len(self.elems)

    def setElset(self, name: str, option: Optional[str], val: List) -> None:
        """
        Define an element set.

        inputs:
        name: [str] name of the elset
        option: [None or str] an optional keyword for the type of elset
        val: [list] list of values for the elset
        """

        valid_options = {None, 'GENERATE', 'INSTANCE', 'INTERNAL', 'UNSORTED'}
        if option not in valid_options:
            raise ValueError('Invalid option value {}'.format(option))

        self.elsets[name] = {
            'option': option,
            'value': val
        }

    def getElset(self, name: str) -> np.ndarray:
        """
        Returns the element numbers of the named elset
        """
        if name in self.elsets:
            elset = self.elsets[name]
            if elset['option'] is None:
                return np.array([int(x) for x in elset['value'].split(',')])
            elif elset['option'] == 'GENERATE':
                e0, e1, estep = [int(x) for x in elset['value'].split(',')]
                return np.arange(e0, e1 + 1, estep)
            else:
                raise NotImplementedError('{} elset not implemented'.format(elset['option']))

    def setSurface(self, elset_name: str, **kwargs) -> None:
        self.surfaces[elset_name] = kwargs

    def calcElemCentroids(self) -> np.ndarray:
        node_mapping = dict(zip(self.nodeNumbers, self.nodes))
        elem_shape = np.array(self.elems).shape
        elem_nodes_flat = np.hstack(self.elems)
        elem_node_coords_flat = np.array([node_mapping[i] for i in elem_nodes_flat])
        elem_node_coords = elem_node_coords_flat.reshape([elem_shape[0], elem_shape[1], 3])
        elem_centroids = elem_node_coords.mean(1)
        return elem_centroids


class InpReader(object):
    """INP reading class
    """

    nodeStartString = '*NODE'
    elementStartString = '*ELEMENT'
    elsetStartString = '*ELSET'

    def __init__(self, filename: str):
        self.filename = filename
        self.meshNames: Optional[List[str]] = None

    def readHeader(self):
        """Reads and returns the file header
        """
        header = []
        with open(self.filename, 'r') as f:
            doScan = 1
            while doScan:
                line = next(f)
                if line == (COMMENTCHARS + '\n'):
                    doScan = 0
                elif line[:2] != COMMENTCHARS:
                    doScan = 0
                else:
                    header.append(line[2:].strip())

        return header

    def readMeshNames(self) -> List[str]:
        """Read and returns a set of all the ELSET names in the file
        """
        meshNames = set()
        with open(self.filename, 'r') as f:
            doScan = 1
            while doScan:
                try:
                    l = next(f)
                except StopIteration:
                    doScan = 0
                else:
                    if 'ELSET' in l:
                        for term in l.split(','):
                            if 'ELSET' in term:
                                if term == '*ELSET':
                                    break
                                else:
                                    meshNames.add(term.split('=')[1].strip())
                                    break

        self.meshNames = list(meshNames)
        return self.meshNames

    def readNodes(self) -> Tuple[List[int], List[List[float]]]:
        nodeNumbers = []
        nodes = []

        # read nodes
        with open(self.filename, 'r') as f:
            doScan = 1
            while doScan:
                try:
                    l = next(f)
                except StopIteration:
                    raise IOError('Cannot find nodes starting with {}'.format(self.nodeStartString))
                else:
                    # if ('NSET='+meshName) in l:
                    if (self.nodeStartString) in l.upper():
                        doScan = 0

            doScan = 1
            while doScan:
                l = next(f).strip()
                if '*' not in l:
                    terms = l.split(',')
                    nodeNumbers.append(int(terms[0]))
                    nodes.append([float(t) for t in terms[1:]])
                else:
                    doScan = 0

            log.debug(('loaded %d nodes' % (len(nodes))))

        return nodeNumbers, nodes

    def readMeshOld(self, mesh_name: str) -> Mesh:
        """Reads and returns the mesh with name meshName.
        Arguments:
        meshName: string matching an NSET/ELSET name in the file
        Returns:
        mesh: a Mesh instance with the read-in mesh parameters
        """
        elem_type = None
        elem_numbers = []
        elems = []

        node_numbers, nodes = self.readNodes()

        # read elements
        with open(self.filename, 'r') as f:
            do_scan = 1
            while do_scan:
                try:
                    line = next(f)
                except StopIteration:
                    raise IOError('No ELSET named ' + mesh_name)
                else:
                    if ('ELSET=' + mesh_name) in line:
                        do_scan = 0
                        for term in line.split(','):
                            if 'TYPE' in term.upper():
                                elem_type = term.split('=')[1].strip()

            try:
                en = ELEMNODES[elem_type]
            except KeyError:
                raise RuntimeError('Unsupported element type: ' + elem_type)

            do_scan = 1
            nCount = -1
            elem = []
            while do_scan:
                try:
                    line = next(f).strip()
                except StopIteration:
                    do_scan = 0
                else:
                    if '*' not in line:
                        terms = [int(i) for i in line.split(',') if i]
                        if len(terms) == 0:
                            terms = [int(line.strip())]
                        for t in terms:
                            if nCount == -1:
                                elem_numbers.append(t)
                                nCount += 1
                            elif nCount < en:
                                elem.append(t)
                                nCount += 1
                                if nCount == en:
                                    elems.append(elem)
                                    elem = []
                                    nCount = -1
                            else:
                                # should be here something bad happened
                                raise RuntimeError('Unexpected error when reading mesh %s', mesh_name)
                    else:
                        do_scan = 0

            log.debug(('loaded %s %s elements' % (len(elems), elem_type)))

        # get only nodes of the mesh
        _nodes_dict = dict(zip(node_numbers, nodes))
        mesh_node_nums = np.unique(np.hstack(elems))
        mesh_node_coords = [_nodes_dict[i] for i in mesh_node_nums]

        mesh = Mesh(mesh_name)
        mesh.setNodes(mesh_node_coords, mesh_node_nums)
        mesh.setElems(elems, elem_numbers, elem_type)

        return mesh

    def readMesh(self, mesh_name: Optional[str] = None) -> Mesh:
        """
        Reads and returns the mesh with name meshName.
        Arguments:
        meshName: string matching an NSET/ELSET name in the file. If none,
            reads Element section with no ELSET name.
        Returns:
        mesh: a Mesh instance with the read-in mesh parameters
        """
        node_numbers, nodes = self.readNodes()
        elem_numbers, elems, elemType = self.readElements(elset=mesh_name)

        # get only nodes of the mesh
        _nodes_dict = dict(zip(node_numbers, nodes))
        mesh_node_nums = np.unique(np.hstack(elems))
        mesh_node_coords = [_nodes_dict[i] for i in mesh_node_nums]

        mesh = Mesh(mesh_name)
        mesh.setNodes(mesh_node_coords, mesh_node_nums)
        mesh.setElems(elems, elem_numbers, elemType)

        return mesh

    def readAllMeshes(self) -> Dict[str, Mesh]:
        """Read in all meshes in the file.
        Returns a dictionary in which keys are mesh names and values are the
        meshes.
        """
        meshNames = self.readMeshNames()
        meshes = {}
        for mn in meshNames:
            meshes[mn] = self.readMesh(mn)

        return meshes

    def readElements(self, elset: Optional[str] = None) -> Tuple[List[int], List[List[int]], str]:
        """
        read elements section
        """
        elem_type = None
        elems = []
        elem_numbers = []

        with open(self.filename, 'r') as f:
            do_scan = 1
            while do_scan:
                try:
                    line = next(f)
                except StopIteration:
                    raise IOError('No Elements')
                else:
                    if (self.elementStartString) in line.upper():
                        if elset is not None:
                            # check if these are the right elset
                            if ('ELSET=' + elset).upper() in line.upper():
                                do_scan = 0
                        else:
                            # dont care about elset, we found elements
                            do_scan = 0

                        for term in line.split(','):
                            if 'TYPE' in term.upper():
                                elem_type = term.split('=')[1].strip()

            try:
                en = ELEMNODES[elem_type]
            except KeyError:
                raise RuntimeError('Unsupported element type: ' + elem_type)

            do_scan = 1
            n_count = -1
            elem = []
            while do_scan:
                try:
                    line = next(f).strip()
                except StopIteration:
                    do_scan = 0
                else:
                    if '*' not in line:
                        terms = [int(i) for i in line.split(',') if i]
                        if len(terms) == 0:
                            terms = [int(line.strip())]
                        for t in terms:
                            if n_count == -1:
                                elem_numbers.append(t)
                                n_count += 1
                            elif n_count < en:
                                elem.append(t)
                                n_count += 1
                                if n_count == en:
                                    elems.append(elem)
                                    elem = []
                                    n_count = -1
                            else:
                                # should be here something bad happened
                                raise RuntimeError

                        # elemNumbers.append(int(terms[0]))
                        # elems.append([int(t) for t in terms[1:]])
                    else:
                        do_scan = 0

            log.debug(('loaded %s %s elements' % (len(elems), elem_type)))

        return elem_numbers, elems, elem_type

    def readElset(self, name: str) -> List[List[int]]:
        # read elset
        elset = []

        with open(self.filename, 'r') as f:
            do_scan = 1
            while do_scan:
                try:
                    l = next(f)
                except StopIteration:
                    raise IOError('No ELSET named ' + name)
                else:
                    if (self.elsetStartString) in l.upper():
                        if ('ELSET=' + name).upper() in l.upper():
                            do_scan = 0

            do_scan = 1
            while do_scan:
                l = next(f).strip()
                if '*' not in l:
                    terms = l.split(',')
                    elset += [int(t) for t in terms]
                else:
                    do_scan = 0

            log.debug('loaded {} elements in elset {}'.format(len(elset), name))

        return elset


class InpWriter(object):
    _commentChars = '**'
    _nodeHeaderLine = '*NODE, NSET={}\n'
    _nodeCounterFormat = '{:6d}'
    _nodeCoordFormat = '{:16.10f}'
    _elemHeaderLine = '*ELEMENT, TYPE={}, ELSET={}\n'
    _elemCounterFormat = '{:6d}'
    _elemNodeFormat = '{:10d}'

    def __init__(self, filename: str, autoFormat: bool=True):
        self._meshes: List[Mesh] = []
        self.filename: str = filename
        self._header: Optional[str] = None
        self.autoFormat: bool = autoFormat

    def addHeader(self, header: str) -> None:
        """
        Add commented text to be written at the top of the 
        file.
        """
        self._header = header
        if header[-1:] != '\n':
            self._header = self._header + '\n'

    def addMesh(self, mesh: Mesh) -> None:
        """Add a mesh to be written.
        Argument:
        mesh : a Mesh object
        """

        self._meshes.append(mesh)

    def addSection(self):
        raise NotImplementedError()

    def addMaterial(self):
        raise NotImplementedError()

    def _autoFormat(self):
        max_nodes = max([max(mesh.nodeNumbers) for mesh in self._meshes])
        node_counter_char_length = len(str(max_nodes)) + 1
        self._nodeCounterFormat = '{:' + str(node_counter_char_length) + '}'

        max_elems = max([max(mesh.elemNumbers) for mesh in self._meshes])
        elem_counter_char_length = len(str(max_elems)) + 1
        self._elemCounterFormat = '{:' + str(elem_counter_char_length) + '}'
        self._elemNodeFormat = '{:' + str(node_counter_char_length) + '}'

    def write(self) -> None:
        """
        Write data to file.
        """

        if len(self._meshes) == 0:
            raise RuntimeError('no meshes defined')

        if self.autoFormat:
            self._autoFormat()

        # open file for writing
        with open(self.filename, 'w') as f:
            # write header comments
            if self._header != None:
                f.write(COMMENTCHARS + self._header)
            else:
                f.write(COMMENTCHARS + '\n')

            f.write(COMMENTCHARS + '\n')

            # write each mesh
            for mesh in self._meshes:

                # nodes
                if mesh.nodes is not None:
                    f.write(self._nodeHeaderLine.format(mesh.name))
                    for ni, n in enumerate(mesh.nodes):
                        f.write(self._getNodeLine(mesh.nodeNumbers[ni], n))

                    f.write(COMMENTCHARS + '\n')

                # elems
                if mesh.elems is not None:
                    f.write(self._elemHeaderLine.format(mesh.elemType, mesh.name))
                    for ei, e in enumerate(mesh.elems):
                        f.write(self._getElemLine(mesh.elemNumbers[ei], e))

                    f.write(COMMENTCHARS + '\n')

                # elsets - TODO

                # sections - TODO

                # materials - TODO

    def _getNodeLine(self, i: int, node: List[float]) -> str:
        words = [self._nodeCounterFormat.format(i), ] + \
                [self._nodeCoordFormat.format(x) for x in node]

        return ', '.join(words) + '\n'

    def _getElemLine(self, i: int, elem: List[int]) -> str:
        words = [self._elemCounterFormat.format(i), ] + \
                [self._elemNodeFormat.format(n) for n in elem]

        return ', '.join(words) + '\n'
