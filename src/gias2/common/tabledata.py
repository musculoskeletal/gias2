"""
FILE: tabledata.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: representing a table of data.

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import numpy as np

class Classification(object):

    def __init__(self, name, labelsDict, code):
        self.name = name
        self.labels = labelsDict
        self.code = code

class TableData(object):

    def __init__(self, name):
        self.name = name
        self._classifications = {}
        self._headers = None
        self._rowLabels = None
        self._units = None
        self._dataArray = None

    def addClassification(self, classification):
        self._classifications[classification.name] = classification

    def setData(self, headers, rows, units, dataArray):
        self._headers = list(headers)
        self._rowLabels = list(rows)
        self._units = list(units)
        self._dataArray = np.array(dataArray)

    def addDataColumn(self, header, unit, data):
        self._headers.append(header)
        self._units.append(unit)
        self._dataArray = np.hstack([self._dataArray, data[:,np.newaxis]])

    def addDataRow(self, rowLabel, data):
        self._rowLabels = self._rowLabels.append(rowLabel)
        self._dataArray = np.vstack([self._dataArray, data])

    def getClasses(self):
        return list(self._classifications.keys())

    def getLabelsForClass(self, classificationName):
        return self._classifications[classificationName].labels

    def getHeaders(self):
        return self._headers

    def getRowLabels(self):
        return self._rowLabels

    def getUnits(self):
        return self._units

    def getUnitsForHeader(self, header):
        return self._units[self._headers.index(header)]

    def getData(self, header, classificationName=None, classLabel=None):
        data = self._dataArray[:,self._headers.index(header)]
        if classificationName!=None:
            C = self._classifications[classificationName]
            data = data[C.code==C.labels[classLabel]]

        return data

    def getRowLabels(self, classificationName=None, classLabel=None):
        if classificationName==None:
            return self._rowLabels
        else:
            C = self._classifications[classificationName]
            return [self._rowLabels[i] for i in np.where(C.code==C.labels[classLabel])[0]]
