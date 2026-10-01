#!/usr/bin/env python
"""
setup.py file for GWforge package
"""

from glob import glob

from setuptools import setup

setup(scripts=glob("bin/gwforge_*"))
