"""
FILE: pelvis_hjc_estimation.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: functions and classes for estimating hip joint centre. Contains
predictive methods of Bell (1989) and Leardini (1999).

All methods require hip models to be in the standard ISB anatomic coordinate
system:
    - x: posterior to anterior
    - y: inferior to superior
    - z: left to right

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import numpy as np

# literature data collection from Bell (1989)
# rp only available for adults class
_literatureDataBell = {'Tylkowski': {
    'all': {'rd': 0.29, 'rm': 0.15, 'rp': 0.22},
    'adults': {'rd': 0.30, 'rm': 0.14, 'rp': 0.22},
    'children': {'rd': 0.29, 'rm': 0.16, 'rp': 0.22},
    'boys': {'rd': 0.30, 'rm': 0.15, 'rp': 0.22},
    'girls': {'rd': 0.28, 'rm': 0.16, 'rp': 0.22},
    'men': {'rd': 0.31, 'rm': 0.14, 'rp': 0.22},
    'women': {'rd': 0.29, 'rm': 0.15, 'rp': 0.22}
},
    'Andriacchi': {
        'adults': {'dd': 41, 'dl': 16},
        'children': {'dd': 23, 'dl': 8.2},
        'boys': {'dd': 25, 'dl': 8.7},
        'girls': {'dd': 22, 'dl': 7.8},
        'men': {'dd': 46, 'dl': 17},
        'women': {'dd': 37, 'dl': 15},
    },
    'Seidel': {
        'adults': {'rd': 0.79, 'rm': 0.14, 'rp': 0.34},
        'men': {'rd': 0.79, 'rm': 0.14, 'rp': 0.34},
        'women': {'rd': 0.79, 'rm': 0.14, 'rp': 0.34},
    },
    'Harrington': {
        'all': {'pdx': -0.24, 'cx': -9.9,
                'pwy': -0.30, 'cy': -10.9,
                'pwz': 0.33, 'cz': 7.3,
                },
        'adults': {'pdx': -0.24, 'cx': -9.9,
                   'pwy': -0.30, 'cy': -10.9,
                   'pwz': 0.33, 'cz': 7.3,
                   },
    }
}

_literatureDataMFC = {'Tylkowski': {
    'adults': {'rd': 0.362, 'rm': 0.139, 'rp': 0.154},
    'men': {'rd': 0.354, 'rm': 0.145, 'rp': 0.144},
    'women': {'rd': 0.370, 'rm': 0.134, 'rp': 0.164}
},
    'Andriacchi': {
        'adults': {'dd': 33.0, 'dl': 29.4},
        'men': {'dd': 34.4, 'dl': 27.6},
        'women': {'dd': 31.5, 'dl': 31.3},
    },
    'Seidel': {
        'adults': {'rd': 0.885, 'rm': 0.124, 'rp': 0.235},
        'men': {'rd': 0.892, 'rm': 0.134, 'rp': 0.225},
        'women': {'rd': 0.878, 'rm': 0.114, 'rp': 0.244},
    },
}

# _literatureData = _literatureDataBell
_literatureData = _literatureDataMFC


def calcEucDist(X, Y):
    return np.sqrt((np.subtract(X, Y) ** 2.0).sum(1))


def HJCTylkowski(LASIS, RASIS, popClass):
    """
    Tylkowski's (1982) method predicts the HJC as a point distal, medial,
    and posterior to the ASIS. The distance distal (rd), medial (rm), and 
    posterior (rp) are proportions of the ASIS-ASIS distance of the pelvis.
    """

    rd = _literatureData['Tylkowski'][popClass]['rd']
    rm = _literatureData['Tylkowski'][popClass]['rm']
    rp = _literatureData['Tylkowski'][popClass]['rp']

    D = np.sqrt(((LASIS - RASIS) ** 2.0).sum())
    LHJC = LASIS + [-rp * D, -rd * D, rm * D]
    RHJC = RASIS + [-rp * D, -rd * D, -rm * D]

    return LHJC, RHJC, D


def HJCTylkowskiInverse(LASIS, RASIS, LHJC, RHJC):
    """
    Calculates the coefficients for Tylkowski's method
    """
    D = np.sqrt(((LASIS - RASIS) ** 2.0).sum())
    LCoeffs = (LHJC - LASIS) / D
    RCoeffs = (RHJC - RASIS) / D
    RCoeffs[:, 2] = RCoeffs[:, 2] * -1.0
    return LCoeffs, RCoeffs


def HJCAndriacchi(LASIS, RASIS, PS, popClass):
    """
    Andriacchi's (1982) method predicts the HJC as a point distal and lateral
    of the midpoint of a line between the ASIS and the pubis symphysis. This 
    distances distal (dd) and medial (dm) are defined in absolute values.
    """

    dd = _literatureData['Andriacchi'][popClass]['dd']
    dl = _literatureData['Andriacchi'][popClass]['dl']

    LO = (LASIS + PS) / 2.0
    RO = (RASIS + PS) / 2.0
    LHJC = LO + [0, -dd, -dl]
    RHJC = RO + [0, -dd, dl]

    return LHJC, RHJC, LO, RO


def HJCAndriacchiInverse(LASIS, RASIS, PS, LHJC, RHJC):
    """
    Calculates the constants for Andriacchi's method
    """
    LO = (LASIS + PS) / 2.0
    RO = (RASIS + PS) / 2.0
    LConst = LHJC - LO
    LConst[:, 2] = LConst[:, 2] * -1.0
    RConst = RHJC - RO
    return LConst, RConst


def HJCBell(LASIS, RASIS, PS, popClass):
    """
    Bell's (1989) method uses Andriacchi's method to predict frontal plane 
    position, and Tylkowski's method to predict antero-posterior position
    """

    LHJC_A, RHJC_A, LO, RO = HJCAndriacchi(LASIS, RASIS, PS, popClass)
    LHJC_T, RHJC_T, D = HJCTylkowski(LASIS, RASIS, popClass)

    LHJC = np.array([LHJC_T[0], LHJC_A[1], LHJC_A[2]])
    RHJC = np.array([RHJC_T[0], RHJC_A[1], RHJC_A[2]])

    return LHJC, RHJC, D, LO, RO


def HJCSeidel(LASIS, RASIS, LPSIS, RPSIS, PS, popClass):
    """
    Seidel's (1995) method predicts HJC as proportions of hip 
    width (ASIS to ASSIS), height (ASIS-ASIS to PS), and 
    depth (ASIS to PSIS)
    """

    W = np.sqrt(((LASIS - RASIS) ** 2.0).sum())
    HL = abs(LASIS[1] - PS[1])
    HR = abs(RASIS[1] - PS[1])
    H = (HL + HR) / 2.0
    DL = abs(LASIS[0] - LPSIS[0])
    DR = abs(RASIS[0] - RPSIS[0])
    D = (DL + DR) / 2.0

    rd = _literatureData['Seidel'][popClass]['rd']
    rm = _literatureData['Seidel'][popClass]['rm']
    rp = _literatureData['Seidel'][popClass]['rp']

    LHJC = LASIS + [-rp * D, -rd * H, rm * W]
    RHJC = RASIS + [-rp * D, -rd * H, -rm * W]

    return LHJC, RHJC, W, H, D


def HJCSeidelInverse(LASIS, RASIS, PS, H, D, LHJC, RHJC):
    W = np.sqrt(((LASIS - RASIS) ** 2.0).sum(1))
    LCoeffs = (LHJC - LASIS) / np.array([D, H, W]).T
    RCoeffs = (RHJC - RASIS) / np.array([D, H, -W]).T

    return LCoeffs, RCoeffs


def HJCHarrington(LASIS, RASIS, LPSIS, RPSIS, popClass):
    """ Harrington (2007) method for predicting right HJC from 
    pelvic depth (PD), and pelvis width (PW) in absolute
    value with the pelvis in the ISB anatomical coordinate system.
    x axis is anterior to posterior
    y axis is inferior to superior
    z axis is left to right
    """
    # calc PW and PD
    PW = np.sqrt(((LASIS - RASIS) ** 2.0).sum())
    AC = (LASIS + RASIS) * 0.5
    PC = (LPSIS + RPSIS) * 0.5
    PD = np.sqrt(((AC - PC) ** 2.0).sum())

    # get regression coeffs and constants
    pdx = _literatureData['Harrington'][popClass]['pdx']
    cx = _literatureData['Harrington'][popClass]['cx']
    pwy = _literatureData['Harrington'][popClass]['pwy']
    cy = _literatureData['Harrington'][popClass]['cy']
    pwz = _literatureData['Harrington'][popClass]['pwz']
    cz = _literatureData['Harrington'][popClass]['cz']

    # calculate HJC
    RHJC = np.array([
        pdx * PD + cx,
        pwy * PW + cy,
        pwz * PW + cz,
    ])

    LHJC = np.array([
        pdx * PD + cx,
        pwy * PW + cy,
        -pwz * PW - cz,
    ])

    return LHJC, RHJC, PW, PD
