#!/usr/bin/env python

import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(name='pyMotif',
        version='0.1.0',
        description='TODO',
        author='Chris Jurich',
        author_email='cjurich2@huskers.unl.edu',
        packages=['pyMotif'],
        )
