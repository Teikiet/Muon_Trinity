#
# (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
#
# See file AUTHORS for a list of contributors.
#
# This software is distributed under the terms of the 3-clause BSD license.
# See file LICENSE for a full version of the license.
#

# Look for an executable called sphinx-build
find_program (SPHINX_EXECUTABLE
              NAMES sphinx-build
              DOC "Path to sphinx-build executable")

include (FindPackageHandleStandardArgs)

#Handle standard arguments to find_package like REQUIRED and QUIET
find_package_handle_standard_args (Sphinx
                                  "Failed to find sphinx-build executable"
                                  SPHINX_EXECUTABLE)
