#
# (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
#
# See file AUTHORS for a list of contributors.
#
# This software is distributed under the terms of the 3-clause BSD license.
# See file LICENSE for a full version of the license.
#

#################################################
#
# run pythia8-config and interpret result
#

function (_Pythia8_CONFIG_ option variable type doc)
  execute_process (COMMAND ${Pythia8_CONFIG} ${option}
    OUTPUT_VARIABLE _local_out_
    RESULT_VARIABLE _local_res_)
  string (REGEX REPLACE "\n$" "" _local_out_ "${_local_out_}")
  if (NOT ${_local_res_} EQUAL 0)
    message ("Error in running ${Pythia8_CONFIG} ${option}")
  else ()
    set (${variable} "${_local_out_}" CACHE ${type} ${doc})
  endif ()
endfunction (_Pythia8_CONFIG_)
  
#################################################
#
# take directory and assume standard install layout
#

function (_Pythia8_LAYOUT_ dir variable type doc)
  set (${variable} "${dir}" CACHE ${type} ${doc})
endfunction (_Pythia8_LAYOUT_)


#################################################
# 
# Searched Pythia8 on system. Expect pythia8-config in PATH, or typical installation location
#
# This module defines
# HAVE_Pythia8
# Pythia8_INCLUDE_DIR   where to locate Pythia.h file
# Pythia8_LIBRARY       where to find the libpythia8 library
# Pythia8_LIBRARIES     (not cached) the libraries to link against to use Pythia8
# Pythia8_VERSION       version of Pythia8 if found
#

set (_SEARCH_Pythia8_
  ${PROJECT_BINARY_DIR}/ThirdParty/pythia8-install
  ${PYTHIA8}
  $ENV{PYTHIA8}
  ${PYTHIA8DIR}
  $ENV{PYTHIA8DIR}
  ${PYTHIA8_ROOT}
  $ENV{PYTHIA8_ROOT}
  ${PYTHIA8_DIR}
  $ENV{PYTHIA8_DIR}
  ${Pythia8_DIR}
  $ENV{Pythia8_DIR}
  /opt/pythia8
  )

find_file (Pythia8_Pythia_h_LOC
  NAME Pythia.h
  PATHS ${_SEARCH_Pythia8_}
  PATH_SUFFIXES include/Pythia8
  DOC "The location of the Pythia8/Pythia.h script"
  REQUIRED)

if ("${Pythia8_Pythia_h_LOC}" STREQUAL "Pythia8_Pythia_h_LOC-NOTFOUND")
  message (FATAL_ERROR "Did not find SYSTEM-level Pythia8 in: \"${_SEARCH_Pythia8_}\"")
endif ()

string (REPLACE "/include/Pythia8/Pythia.h" "" Pythia8_DIR ${Pythia8_Pythia_h_LOC})

set (Pythia8_CONFIG ${Pythia8_DIR}/bin/pythia-config)
if (Pythia8_CONFIG)
  set (HAVE_Pythia8 1 CACHE BOOL "presence of pythia8, found via pythia8-config")

  # pythia-config is not relocatable
  #_Pythia8_CONFIG_ ("--prefix" Pythia8_PREFIX PATH "location of pythia8 installation")
  #_Pythia8_CONFIG_ ("--includedir" Pythia8_INCLUDE_DIR PATH "pythia8 include directory")
  #_Pythia8_CONFIG_ ("--libdir" Pythia8_LIBRARY STRING "the pythia8 libs")
  #_Pythia8_CONFIG_ ("--datadir" Pythia8_DATA_DIR PATH "the pythia8 data dir")
  _Pythia8_LAYOUT_ ("${Pythia8_DIR}" Pythia8_PREFIX PATH "location of pythia8 installation")
  _Pythia8_LAYOUT_ ("${Pythia8_DIR}/include" Pythia8_INCLUDE_DIR PATH "pythia8 include directory")
  _Pythia8_LAYOUT_ ("${Pythia8_DIR}/lib" Pythia8_LIBRARY STRING "the pythia8 libs")
  _Pythia8_LAYOUT_ ("${Pythia8_DIR}/share/Pythia8/xmldoc" Pythia8_DATA_DIR PATH "the pythia8 data dir")
  
  # read the config string
  file (READ "${Pythia8_INCLUDE_DIR}/Pythia8/Pythia.h" Pythia8_TMP_PYTHIA_H)
  string (REGEX MATCH "#define PYTHIA_VERSION_INTEGER ([0-9]*)" _ ${Pythia8_TMP_PYTHIA_H})  
  set (Pythia8_VERSION ${CMAKE_MATCH_1})
  message (STATUS "Found Pythia8 version: ${Pythia8_VERSION}")
  
endif ()

# standard cmake infrastructure:
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (Pythia8 REQUIRED_VARS Pythia8_PREFIX Pythia8_INCLUDE_DIR Pythia8_LIBRARY VERSION_VAR Pythia8_VERSION)
mark_as_advanced (Pythia8_DATA_DIR Pythia8_INCLUDE_DIR Pythia8_LIBRARY)
