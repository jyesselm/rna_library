#!/usr/bin/env python

import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

with open("requirements.txt", "r") as f:
    requirements = f.read().splitlines()

setup(
    name="rna_library",
    version="0.1.0",
    description="TODO",
    author="Chris Jurich",
    author_email="cjurich2@huskers.unl.edu",
    packages=[
        "rna_library",
        "rna_library.core",
        "rna_library.structure",
    ],
    include_package_data=True,
    install_requires=requirements,
)
