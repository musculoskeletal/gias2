"""
===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import sys

# needed to open PCA files pickled by gias(1)
import gias2
from gias2 import learning
from gias2.common import transform3D, geoprimitives
from gias2.fieldwork.field import geometric_field
from gias2.learning import PCA, PCA_fitting
from gias2.version import __version__, __version_info__

# for compatibility with gias(1) layout
sys.modules['gias'] = gias2
sys.modules['gias.learning'] = learning
sys.modules['gias.learning.PCA'] = PCA
