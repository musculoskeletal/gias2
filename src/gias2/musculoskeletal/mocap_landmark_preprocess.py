"""
FILE: mocap_landmark_preprocess.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: module for preprocessing mocap landmarks

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

from gias2.common import math


class InsufficientLandmarksError(Exception):
    pass


def preprocess_lower_limb(marker_radius, skin_pad, LASIS, RASIS, sacral, LEC,
                          MEC, LM, MM, ):
    """Returns adjusted coordinates of LASIS, RASIS, Sacral, LEC, MEC, LM, MM
    """
    pelvis_output = preprocess_pelvis(
        marker_radius, skin_pad, LASIS, RASIS, None, None, sacral
    )
    femur_output = preprocess_femur(
        marker_radius, skin_pad, LEC, MEC
    )
    tibiafibula_output = preprocess_tibiafibula(
        marker_radius, skin_pad, LM, MM
    )
    return tuple(list(pelvis_output) + list(femur_output) + list(tibiafibula_output))


def preprocess_pelvis(marker_radius, skin_pad, LASIS, RASIS, LPSIS, RPSIS, sacral):
    if sacral is None:
        if (LPSIS is None) or (RPSIS is None):
            raise (
                InsufficientLandmarksError,
                'Either Sacral or the PSISs must be defined'
            )

    if (LASIS is None) or (RASIS is None):
        raise (InsufficientLandmarksError, 'LASIS and RASIS must be defined')

    # calculate AP axis
    oa = (LASIS + RASIS) / 2.0
    if (LPSIS is not None) and (RPSIS is not None):
        op = (LPSIS + RPSIS) / 2.0
    else:
        op = sacral

    ap = math.norm(op - oa)

    # shift ASIS posteriorly
    LASIS2 = LASIS + ap * (marker_radius + skin_pad)
    RASIS2 = RASIS + ap * (marker_radius + skin_pad)

    # shift PSIS or sacrum anteriorly
    if LPSIS is not None:
        LPSIS2 = LPSIS - ap * (marker_radius + skin_pad)
    else:
        LPSIS2 = LPSIS
    if RPSIS is not None:
        RPSIS2 = RPSIS - ap * (marker_radius + skin_pad)
    else:
        RPSIS2 = RPSIS
    if sacral is not None:
        sacral2 = sacral - ap * (marker_radius + skin_pad)
    else:
        sacral2 = sacral

    return LASIS2, RASIS2, LPSIS2, RPSIS2, sacral2


def preprocess_femur(marker_radius, skin_pad, LEC, MEC):
    if LEC is None:
        raise (InsufficientLandmarksError, 'LEC must be defined')

    if MEC is None:
        raise (InsufficientLandmarksError, 'MEC must be defined')

    # calculate epicondylar axis
    ML = math.norm(LEC - MEC)

    # shift MEC laterally
    MEC2 = MEC + ML * (marker_radius + skin_pad)

    # shift LEC medially
    LEC2 = LEC - ML * (marker_radius + skin_pad)

    return LEC2, MEC2


def preprocess_tibiafibula(marker_radius, skin_pad, LM, MM):
    if LM is None:
        raise (InsufficientLandmarksError, 'LM must be defined')

    if MM is None:
        raise (InsufficientLandmarksError, 'MM must be defined')

    # calculate epimalleolus axis
    ML = math.norm(LM - MM)

    # shift MM laterally
    MM2 = MM + ML * (marker_radius + skin_pad)

    # shift LEC medially
    LM2 = LM - ML * (marker_radius + skin_pad)

    return LM2, MM2


preprocessors = {
    'pelvis': preprocess_pelvis,
    'femur': preprocess_femur,
    'tibiafibula': preprocess_tibiafibula,
}
