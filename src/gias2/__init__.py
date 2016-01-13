"""
===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

from gias2.common import transform3D, geoprimitives
from gias2.fieldwork.field import geometric_field, scalar_field
from gias2.fieldwork.field.tools import fitting_tools, mesh_fitter
from gias2.image_analysis import fw_segmentation_tools, image_tools
from gias2.io import cmissio
from gias2.learning import PCA, PCA_fitting, kernelregression, ols
from gias2.mesh import vtktools, simplemesh, inp, tetgenoutput
from gias2.musculoskeletal import fw_model_landmarks, model_alignment,\
    model_alignment_multi, pelvis_hjc_estimation, fw_femur_measurements,\
    fw_pelvis_measurements
from gias2.musculoskeletal.bonemodels import bonemodels, lowerlimbatlasfit,\
    lowerlimbatlasfitscaling
from gias2.registration import alignment_analytic, alignment_fitting

import sys
# needed to open PCA files pickled before PCA got moved to learning
sys.modules['gias2.common.PCA'] = PCA

