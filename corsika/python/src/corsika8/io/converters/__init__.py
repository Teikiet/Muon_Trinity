"""
The 'io' module provides for reading CORSIKA8 output files.

(c) Copyright 2018 CORSIKA Project, corsika8@kit.edu

This software is distributed under the terms of the 3-clause BSD license.
See file LICENSE for a full version of the license.
"""

from .arrow_to_numpy import convert_to_numpy

__all__ = ["convert_to_numpy"]
