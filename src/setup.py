"""
===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

#!/usr/bin/env python
from setuptools import setup

def readme():
     with open('readme.md', 'r') as f:
          return f.read()

version = '0.1'
install_requires = [
     'numpy',
     'scipy',
     'scikit-learn',
     'matplotlib'
]

setup(
     name='gias2',
     version=version,
     description='A library of musculoskeletal modelling tools.',
     long_description=readme(),
     author='MAP Client Developers',
     url='https://bitbucket.org/jangle/gias2',
     install_requires=install_requires,
     keywords='musculoskeletal map mapclient',
     license='mozilla',
     classifiers=[
          'Development Status :: 3 - Alpha',
          'License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)',
          'Programming Language :: Python :: 2.7',
          'Topic :: Scientific/Engineering :: Medical Science Apps.'
     ],
)