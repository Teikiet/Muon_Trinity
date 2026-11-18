#
# (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
#
# See file AUTHORS for a list of contributors.
#
# This software is distributed under the terms of the 3-clause BSD license.
# See file LICENSE for a full version of the license.
#

add_library (PhysUnits INTERFACE)

target_compile_options (PhysUnits
  INTERFACE
  -I${CMAKE_SOURCE_DIR}/externals/phys_units
  )

set (PhysUnits_FOUND True)
