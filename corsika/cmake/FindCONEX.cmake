#
# (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
#
# See file AUTHORS for a list of contributors.
#
# This software is distributed under the terms of the 3-clause BSD license.
# See file LICENSE for a full version of the license.
#

# - find Conex
#
# This module defines
# CONEX_PREFIX
# CONEX_INCLUDE_DIR
# CONEX_LIBRARY

set (SEARCH_conex_
  ${WITH_CONEX}
  ${CONEXROOT}
  ${CONEX_ROOT}
  $ENV{CONEXROOT}
  $ENV{CONEX_ROOT}
  )

find_path (CONEX_PREFIX
  NAMES lib/${CMAKE_SYSTEM_NAME}
  PATHS ${SEARCH_conex_}
  DOC "The CONEX root directory"
  NO_DEFAULT_PATH
  )

find_path (CONEX_INCLUDE_DIR
  NAMES ConexDynamicInterface.h
  PATHS ${CONEX_PREFIX}
  PATH_SUFFIXES src
  DOC "The CONEX include directory"
  )

find_library (CONEX_LIBRARY
  NAMES libCONEXdynamic.a
  PATHS ${CONEX_PREFIX}
  PATH_SUFFIXES lib/${CMAKE_SYSTEM_NAME}
  DOC "The CONEX library"
  )

# standard cmake infrastructure:
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (CONEX
  "Did not find system-level CONEX."
  CONEX_INCLUDE_DIR CONEX_LIBRARY CONEX_PREFIX)
mark_as_advanced (CONEX_INCLUDE_DIR CONEX_LIBRARY CONEX_PREFIX)
