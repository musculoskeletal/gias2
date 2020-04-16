"""
FILE: inptest.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: Test INP reading and writing

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""
import logging

from gias2.mesh import inp

log = logging.getLogger(__name__)

inputFilename = 'data/prox_femur.inp'
outputFilename = 'data/prox_femur_out.inp'


def main():
    reader = inp.InpReader(inputFilename)
    header = reader.readHeader()
    log.debug(('header: ' + ' '.join(header)))
    meshnames = reader.readMeshNames()
    log.debug(('mesh names: ' + ', '.join(meshnames)))
    mesh = reader.readMesh(meshnames[0])
    log.debug((mesh.getNumberOfElems()))
    log.debug((mesh.getNumberOfNodes()))
    log.debug((mesh.getNode(10)))
    log.debug((mesh.getElem(10)))
    log.debug((mesh.getElemType()))

    writer = inp.InpWriter(outputFilename)
    writer.addHeader(header[0])
    writer.addMesh(mesh)
    writer.write()


if __name__ == '__main__':
    main()
