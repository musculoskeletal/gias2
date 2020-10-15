"""
===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import re
import codecs
from os import path

# !/usr/bin/env python
import sys
from setuptools import setup, find_packages
from setuptools.extension import Extension

here = path.abspath(path.dirname(__file__))


def read(*parts):
    with codecs.open(path.join(here, *parts), 'r') as fp:
        return fp.read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


def readme():
    with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
        return f.read()


def requirements():
    requirement_path = path.join(here, 'requirements.txt')
    with open(requirement_path, 'r') as f:
        return list(f.read().splitlines())


# =============================================================================#
name = 'gias2'
version = find_version('src', 'gias2', 'version.py')
setup_requires = [
    'setuptools>=40.0',
    'cython>=0.29.7',
    'numpy>=1.16.2'
]
install_requires = requirements()
if sys.version_info.major == 2:
    install_requires.append('ConfigParser')
else:
    install_requires.append('configparser')

package_data = {
    'gias2': [
        'src/gias2/examples/data/*',
        'src/gias2/examples/outputs/*.md',
        'src/gias2/examples/data/tetgen_mesh/*',
        'src/gias2/examples/fieldwork/data/*',
        'src/gias2/examples/fieldwork/fit_whole_pelvis_data/*',
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
    'Programming Language :: Python :: 3.5',
    'Topic :: Scientific/Engineering :: Medical Science Apps.'
]
entry_points = {
    'console_scripts': [
        'gias-rbfreg=gias2.applications.giasrbfreg:main',
        'gias-rigidreg=gias2.applications.giasrigidreg:main',
        'gias-pcreg=gias2.applications.giaspcreg:main',
        'gias-trainpcashapemodel=gias2.applications.giastrainpcashapemodel:main',
        'gias-surfacedistance=gias2.applications.giassurfacedistance:main',
        'gias-hmfinp2surf=gias2.applications.giashmfinp2surf:main',
        'gias-inpsampledicom=gias2.applications.giasinpsampledicom:main'
    ]
}


# Numerical python is a transitive dependency that is required for build only.
# If it is required for install, then all distributions that use it will likely
# fail as they will have to express the dependency.
#
# see
#  - https://stackoverflow.com/questions/19919905/how-to-bootstrap-numpy-installation-in-setup-py
def any_argv(s):
    for arg in sys.argv:
        if arg == s:
            return True
    return False


if any_argv('bdist_wheel'):
    print('Performing build with numpy')
    import numpy
    from Cython.Build import cythonize

    np_include_dirs = [numpy.get_include()]
    cython_modules = cythonize(
        [
            Extension('gias2.image_analysis.asm_search_c', ['src/gias2/image_analysis/asm_search_c.pyx']),
            Extension('gias2.image_analysis.integralimagec', ['src/gias2/image_analysis/integralimagec.pyx']),
        ])
    print(np_include_dirs)
else:
    print('Not performing a build with numpy')
    np_include_dirs = None
    cython_modules = None

# =============================================================================#
if __name__ == '__main__':
    setup(
        name=name,
        version=version,
        description=description,
        long_description=readme(),
        long_description_content_type='text/x-rst',
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
        entry_points=entry_points,
        ext_modules=cython_modules,
        include_dirs=np_include_dirs,
        setup_requires=setup_requires
    )
