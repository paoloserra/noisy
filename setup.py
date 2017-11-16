#!/usr/bin/env python

import os

try:
  from setuptools import setup
except ImportError as e:
  from distutils.core import setup

requirements = [
'numpy>=1.11.0',
'python-casacore>=2.1.2',
'matplotlib>=2.1.0'
]

PACKAGE_NAME = 'noisy'

__version__ = 0.1

setup(name = PACKAGE_NAME,
    version = __version__,
    description = "Noise estimation tool",
    author = "Paolo Serra",
    author_email = "bhugo@ska.ac.za",
    url = "https://github.com/paoloserra/noisy",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Astronomy"
    ],
    packages=[PACKAGE_NAME],
    install_requires = requirements,
    include_package_data = True,
    scripts = ["bin/noisy_predictrms.py"])
