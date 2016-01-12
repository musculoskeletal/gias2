"""
===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import os, sys

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    # Using __file__ will not work if py2exe is used,
    # Possible problem of OSX10.6 also.
    sys.path.insert(0, current_dir)

# import class that derives itself from the step mountpoint.
from mayaviscenewidget import MayaviSceneWidget
from mayaviviewerobjects import colours, MayaviViewerObjectsContainer, MayaviViewerSceneObject, MayaviViewerObject
from mayaviviewerfieldworkmeasurements import MayaviViewerFemurMeasurements
from mayaviviewerfieldworkmodel import MayaviViewerFieldworkModel
from mayaviviewergiasscan import MayaviViewerGiasScan
from mayaviviewerdatapoints import MayaviViewerDataPoints
from mayaviviewerlandmark import MayaviViewerLandmark
from mayaviviewerimageplane import MayaviViewerImagePlane