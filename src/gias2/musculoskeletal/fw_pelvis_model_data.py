"""
FILE: fw_pelvis_model_data.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: Data about the fieldwork pelvis model

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

pelvisCubicBasisTypes = {'tri10': 'simplex_L3_L3', 'quad44': 'quad_L3_L3'}
combinedGFElemNumber = {'left_hemi': 1, 'right_hemi': 0, 'sacrum': 2}

# left_hemi_sacrum ens and mesh using right_hemi ens and mesh submeshes
# left_hemi_right_hemi ens and mesh are flipped: 0: right, 1: left (actual order is still left, right)
landmarksNodes = {'LASIS': 1004,  # anterior superior iliac spine
                  'RASIS': 466,
                  'LPSIS': 924,  # posterior superior iliac spine
                  'RPSIS': 384,
                  'LPS': 631,  # pubis symphysis
                  'RPS': 91,
                  'LIS': 773,  # ischial spine
                  'RIS': 233,
                  'LIT': 566,  # ichial tuberosity (most inferior point)
                  'RIT': 25,
                  'LOFI': 549,  # obturator foramen inferior point
                  'ROFI': 10,
                  'LPAN': 827,  # posterior acetabular notch
                  'RPAN': 289,
                  'RPT': 102,  # pubic tubercle
                  'LPT': 643,
                  }

# lasis = 466 old 464 new
LHLandmarkNodes = {'LASIS': 464,
                   'LPSIS': 384,
                   'LPS': 92,
                   'LIS': 233,
                   'LIT': 25,
                   'LOFI': 10,
                   'LPAN': 289,
                   }

# LHAcetabulumElements = [38,39,40,41,42]
LHAcetabulumElements = [36, 37, 38, 39, 40, 41, 42]
hemiPelvisAcetabulumElements = [36, 35, 38, 39, 40, 41, 42]
hemiPelvisAcetabulumSurroundElements = [64, 66, 44, 43, 27, 28, 29, 26, 24, 25, 23, 20, 21, 16, 31, 30, 33, 47, 34, 72,
                                        65, 14, 13, 15]
femurHeadElements = [0, ]
hemiPelvisSymphysisPubisNodes = [91, 90, 81]
hemiPelvisSacrumNodes = [409, 388, 386, 363, 364, 367, 368, 366, 365, 362, 406]
hemiPelvisAcetabularCupRimNodes = (184, 306, 317, 318, 324, 325, 272, 273)
