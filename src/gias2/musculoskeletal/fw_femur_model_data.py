"""
FILE: fw_femur_model_data.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: Data about the fieldwork femur model

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

landmarkLabels = {'body': 0,
                  'head': 1,
                  'greatertrochanter': 2,
                  'medialcondyle': 6,
                  'lateralcondyle': 5,
                  'distalshaft': 4,
                  'proximalshaft': 3,
                  'any': -1
                  }

# boundaries defined by the regions they seperate
boundaries = {1: ['head', 'proximalshaft'],
              2: ['greatertrochanter', 'proximalshaft'],
              3: ['proximalshaft', 'distalshaft'],
              4: ['distalshaft', 'medialcondyle'],
              5: ['distalshaft', 'lateralcondyle'],
              }

boundaryNames = {1: 'head',
                 2: 'greatertrochanter',
                 3: 'midshaft',
                 4: 'medialcondyle',
                 5: 'lateralcondyle',
                 }

boundaryNumbers = {'head': 1,
                   'greatertrochanter': 2,
                   'midshaft': 3,
                   'medialcondyle': 4,
                   'lateralcondyle': 5,
                   }
# number of nodes on boundaries       
boundaryNNodes = {1: 16, 2: 24, 3: 12, 4: 20, 5: 20}

# boundary number of regions
regionBoundaries = {'head': (1,),
                    'greatertrochanter': (2,),
                    'medialcondyle': (4,),
                    'lateralcondyle': (5,),
                    'proximalshaft': (1, 2, 3),
                    'distalshaft': (3, 4, 5),
                    'ensemble': (1, 2, 3, 4, 5)
                    }

# boundary node numbers of each region's mesh
regionBoundaryNodes = {'head': {1: [25, 26, 27, 28, 29, 71, 72, 73, 74, 97, 98, 99, 45, 46, 47, 48, ]},
                       # ~ 'greatertrochanter': {2: [24,19,14,9,4,3,2,1,0,99,98,97,68,67,66,65,48,47,46,45,28,27,26,25]},
                       'greatertrochanter': {
                           2: [4, 34, 29, 24, 19, 18, 17, 16, 15, 58, 57, 56, 55, 71, 74, 76, 35, 40, 45, 50, 14, 13,
                               11, 8]},
                       # ~ 'medialcondyle': {4:[0,61,62,63,64,93,94,95,96,113,114,115,44,43,42,41,14,12,9,5]},
                       'medialcondyle': {
                           4: [20, 15, 10, 5, 0, 25, 26, 27, 28, 35, 36, 37, 38, 70, 72, 73, 60, 55, 50, 45]},
                       # ~ 'lateralcondyle': {5:[15,64,63,62,61,83,82,81,38,37,36,35,4,3,2,1,0,18,17,16]},
                       'lateralcondyle': {
                           5: [0, 5, 9, 12, 14, 15, 19, 22, 24, 41, 42, 43, 44, 66, 65, 63, 48, 47, 46, 45]},
                       'proximalshaft': {1: [4, 77, 78, 79, 61, 62, 63, 64, 44, 39, 34, 29, 24, 19, 14, 9],
                                         2: [0, 89, 87, 84, 80, 160, 159, 158, 136, 137, 138, 139, 113, 117, 121, 125,
                                             40, 35, 30, 25, 20, 15, 10, 5],
                                         3: [161, 162, 163, 164, 165, 176, 177, 178, 179, 188, 189, 190]},
                       'distalshaft': {3: [155, 156, 157, 158, 159, 168, 169, 170, 171, 178, 179, 180],
                                       4: [28, 27, 26, 25, 4, 3, 2, 1, 0, 68, 67, 66, 65, 87, 86, 85, 48, 47, 46, 45],
                                       5: [15, 16, 17, 18, 14, 130, 131, 132, 120, 121, 122, 123, 100, 101, 102, 103,
                                           35, 36, 37, 38]},
                       }

wholeMeshBoundaryNodes = {
    1: (  ),
    2: (  ),
    3: (  ),
    4: (  ),
    5: (  ),
}

assemblyElements = ('head', 'greatertrochanter', 'proximalshaft', 'distalshaft', 'lateralcondyle', 'medialcondyle')
assemblyElementsNumbers = {'head': 0, 'greatertrochanter': 1, 'proximalshaft': 2, 'distalshaft': 3, 'lateralcondyle': 4,
                           'medialcondyle': 5}
proximalFemurAssemblyElements = {'head': 0, 'greatertrochanter': 1, 'proximalshaft': 2}
distalFemurAssemblyElements = {'ICN': 0, 'lateralcondyle': 1, 'medialcondyle': 2}
