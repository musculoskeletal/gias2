"""
===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

#!/usr/bin/env python
from setuptools import setup, find_packages

def readme():
     with open('README.md', 'r') as f:
          return f.read()

#=============================================================================#
name = 'gias2'
version = '0.3.2'
install_requires = [
     'numpy >= 1.6.1',
     'scipy >= 0.9',
     'scikit-learn >= 0.15',
]
package_data = {
     'gias2': [
          'src/gias2/examples/data/*', 
          'src/gias2/examples/data/tetgen_mesh/*',
     ],
}
include_package_data = True
description = 'A library of musculoskeletal modelling tools.'
author = 'MAP Client Developers'
url = 'https://bitbucket.org/jangle/gias2'
keywords = 'musculoskeletal map mapclient'
license = 'mozilla'
classifiers = [
     'Development Status :: 3 - Alpha',
     'License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)',
     'Programming Language :: Python :: 2.7',
     'Topic :: Scientific/Engineering :: Medical Science Apps.'
]

#=============================================================================#
if __name__ == '__main__':
     setup(
          name=name,
          version=version,
          description=description,
          long_description=readme(),
          packages=find_packages(where="src"),
          package_data=package_data,
          include_package_data=include_package_data,
          package_dir={"": "src"},
          classifiers=classifiers,
          author=author,
          url=url,
          install_requires=install_requires,
          keywords=keywords,
          license=license,  
     )