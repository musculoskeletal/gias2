"""
FILE: cmissio.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: functions and classes for reading and writing cmiss files:

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""
import logging

from numpy import array, arange, loadtxt

log = logging.getLogger(__name__)


# ======================================================================#
def readIpdataOld(fileName):
    """ reads ipdata file and returns the x y z coords on data points
    and the header if there is one
    """

    try:
        file = open(fileName, 'r')
    except IOError:
        log.debug('ERROR: readIpdata: unable to open', fileName)
        return

    header = None
    lines = file.readlines()
    if lines[0][0] != ' ':
        header = lines[0]
        del lines[0]

    coords = []
    nodeNumbers = []
    for l in lines:
        words = l.split()
        coords.append(array(words[1:4], dtype=float))
        nodeNumbers.append(int(words[0]))
    coords = array(coords)
    nodeNumbers = array(nodeNumbers)

    return coords, nodeNumbers, header


def readIpdata(fileName, hasHeader=True):
    """ reads ipdata file and returns the field values of data points,
    data point numbers, and the header if there is one
    """

    with open(fileName, 'r') as f:
        if hasHeader:
            header = f.readline()
        else:
            header = None

        f.seek(0)

        d = loadtxt(f, skiprows=1)
        nodeNumbers = d[:, 0]
        values = d[:, 1:]

    return values, nodeNumbers, header


# ======================================================================#
def writeIpdata(d, filename, header=None, nodeNumbers=None):
    """ write the coordinates of points in 3D space to ipdata file. Each
    row in d is [x,y,z] of a datapoint. filename is a string, header is
    a string. if ex!=False, uses cmConvert to conver to ex formates. ex 
    can be 'data', 'node' or 'both'
    """

    outputFile = open(filename, 'w')
    if header:
        outputFile.write(header + '\n')

    if nodeNumbers is None:
        nodeNumbers = arange(len(d))

    for i, di in enumerate(d):
        outputFile.write("%10i %3.10f %3.10f %3.10f 1.0 1.0 1.0\n" % (nodeNumbers[i], di[0], di[1], di[2]))
    outputFile.close()

    return


# ======================================================================#
def readIpnode(fileName, extra=False):
    """ reads ipnode node files and returns a list of node parameters,
    if extra=True, also returns a list of whether each parameter is a 
    value or derivative, and a list of the node number of parameter
    """

    try:
        file = open(fileName, 'r')
    except IOError:
        log.debug('ERROR: ipnode_reader: unable to open', fileName)
        return

    parameters = []
    parameterType = []
    node = []
    currentNode = None

    if extra:
        for line in file.readlines():
            if line.find('Node number') != -1:
                currentNode = int(line.strip().split()[-1])

            if line.find('Xj') != -1:
                parameters.append(float(line.strip().split()[-1]))
                parameterType.append('value')
                node.append(currentNode)
            elif line.find('The derivative') != -1:
                parameters.append(float(line.strip().split()[-1]))
                parameterType.append('derivative')
                node.append(currentNode)
    else:
        for line in file.readlines():
            if (line.find('Xj') != -1) or (line.find('The derivative') != -1):
                parameters.append(float(line.strip().split()[-1]))

    file.close()

    if extra:
        return array(parameters), parameterType, node
    else:
        return array(parameters)


# ======================================================================#
def writeIpnode(templateName, writeName, header, data):
    """ writes mesh parameters in data to an ipnode file based on a 
    template ipnode file
    """

    s = '%21.14f'

    try:
        template = open(templateName, 'r')
    except IOError:
        log.debug('ERROR: writeIpnode: unable to open template file', templateName)
        return

    try:
        writeFile = open(writeName, 'w')
    except IOError:
        log.debug('ERROR: writeIpnode: unable to open write file', writeName)
        return

    dataCount = 0
    for line in template.readlines():
        if line.find('Heading') != -1:
            writeFile.write(line.split(':')[0] + ': ' + header + '\n')
        elif (line.find('Xj') != -1) or (line.find('The derivative') != -1):
            writeFile.write(line.split(':')[0] + ': ' + s % (data[dataCount]) + '\n')
            dataCount += 1
        else:
            writeFile.write(line)

    template.close()
    writeFile.close()

    return


# ======================================================================#
def readExdata(filename, fieldsToRead=None):
    """Reads exdata file and returns the descriptions of fields and the 
    field values themselves for each data point. If the fieldsToRead are not
    defined as input list, all fields are returned.
    """

    try:
        f = open(filename, 'r')
    except IOError:
        log.debug('ERROR: unable to open', filename)
        return

    # read header
    header = f.readline().strip()

    # read number of fields
    nFieldsLine = f.readline()
    nFields = int(nFieldsLine.strip().split('=')[-1])

    # parse fields
    fields = {}
    nCompsPerField = []
    for fi in range(nFields):
        fieldLine = f.readline()
        field = {}
        fieldNumber, fieldInfo = fieldLine.split(')')
        fieldNumber = int(fieldNumber.strip())
        fieldInfo = fieldInfo.split(',')
        field['name'] = fieldInfo[0].strip()
        field['type'] = fieldInfo[1].strip()
        field['cs'] = fieldInfo[2].strip()
        field['n_components'] = int(fieldInfo[3].split('=')[1].strip())
        nCompsPerField.append(field['n_components'])
        field['components'] = {}
        # read each component for current field
        for ci in range(field['n_components']):
            componentLine = f.readline()
            componentName = componentLine.split('.')[0].strip()
            componentInfo = componentLine.split('.')[1].strip()
            component = {'info': componentInfo}
            field['components'][componentName] = component
        fields[fieldNumber] = field

    # read node values (all fields)
    readingValues = True
    values = []
    nodes = []
    while readingValues:
        # read line
        try:
            l = next(f)
        except StopIteration:
            readingValues = False
        else:
            # check if now node line
            if 'Node:' in l:
                nodes.append(int(l.strip().split(':')[-1]))
                values.append([])
            # else is a values line
            else:
                # read values on line
                values[-1] += [float(x) for x in l.strip().split()]

    # organise node values into fields
    if fieldsToRead is None:
        fieldsToRead = set(range(1, nFields + 1))
    else:
        fieldsToRead = set(fieldsToRead)

    values2 = []
    for nodeValues in values:
        values2.append([])
        vi = 0
        for fi, fn in enumerate(nCompsPerField):
            fieldValues = nodeValues[vi:vi + fn]
            if fi + 1 in fieldsToRead:
                values2[-1].append(fieldValues)
            vi += fn

    return header, fields, nodes, values2


# ======================================================================#
def writeExdata(templateName, writeName, header, data, fields):
    """ writes data to an exdata file based on a 
    template exdata file
    """
    s = '%.9e'

    try:
        template = open(templateName, 'r')
    except IOError:
        log.debug('ERROR: writeExdata: unable to open template file', templateName)
        return

    try:
        writeFile = open(writeName, 'w')
    except IOError:
        log.debug('ERROR: writeExdata: unable to open write file', writeName)
        return

    lines = template.readlines()

    dataCount = 0
    writingField = False
    currentField = 0
    di = 0
    ni = -1
    for li, l in enumerate(lines):
        if 'Group name' in l:
            lines[li] = l.split(':')[0] + ': ' + header + '\n'
        elif 'Node:' in l:
            currentField = 1
            di = 0
            writingField = True
            ni += 1
        elif writingField:
            if currentField in fields:
                d = data[ni][di]
                if isinstance(d, float) or (len(d) == 1):
                    lines[li] = '  ' + s % (d) + '\n'
                else:
                    lines[li] = '  ' + ' '.join([s % (dx) for dx in d]) + '\n'

                di += 1
            currentField += 1

    writeFile.writelines(lines)

    template.close()
    writeFile.close()

    return


# ======================================================================#
def writeXYZ(data, filename, header=None):
    """ writes 3D coords of point to file. Each line contains the x, y,
    and z coords of a point. Optional head in the 1st line
    """

    fOut = open(filename, 'w')
    if header:
        fOut.write(header + '\n')

    for d in data:
        fOut.write("%(x)10.6f\t%(y)10.6f\t%(z)10.6f\n" % {'x': d[0], 'y': d[1], 'z': d[2]})

    fOut.close()

    return 1
