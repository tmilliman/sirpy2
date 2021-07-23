# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import print_function

import io
import re
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext

from setuptools import setup
from setuptools import find_packages


with open("README.rst") as f:
    readme = f.read()

with open("LICENSE") as f:
    license = f.read()

setup(
    name="sirpy2",
    version="0.1.0",
    description="Package to read .sir files with python",
    long_description=readme,
    author="Thomas Milliman",
    author_email="thomas.milliman@unh.edu",
    url="https://github.com/tmilliman/sirpy2",
    license=license,
    packages=find_packages(exclude=("tests", "docs")),
    py_modules=[splitext(basename(path))[0] for path in glob("sir2py/*.py")],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
        "Operating System :: POSIX",
        "Operating System :: Microsoft :: Windows",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    install_requires=["numpy>=1.13.0", "matplotlib>=1.5.3", "pypng>=0.0.20"],
    extras_require={
        # eg:
        "rst": ["docutils>=0.11"],
    },
    entry_points={
        "console_scripts": [
            "printsirhead=sirpy2.printsirhead:main",
            "showsir=sirpy2.showsir:main",
            "sir2png=sirpy2.sir2png:main",
        ]
    },
)
