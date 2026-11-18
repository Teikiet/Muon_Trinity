#
# (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
#
# See file AUTHORS for a list of contributors.
#
# This software is distributed under the terms of the 3-clause BSD license.
# See file LICENSE for a full version of the license.
#


#+++++++++++++++++++++++++++++
# as long as there still are modules using it:
#
enable_language (Fortran)
set (CMAKE_Fortran_FLAGS "-std=legacy -Wfunction-elimination")


#+++++++++++++++++++++++++++++
# Build types settings
#
# setup coverage build type
set (CMAKE_CXX_FLAGS_COVERAGE "-g --coverage")
set (CMAKE_EXE_LINKER_FLAGS_COVERAGE "--coverage")
set (CMAKE_SHARED_LINKER_FLAGS_COVERAGE "--coverage")
# set a flag to inform code that we are in debug mode
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_C8_DEBUG_")


#+++++++++++++++++++++++++++++
# Build type selection
#
# Set the possible values of build type for cmake-gui and command line check
set (ALLOWED_BUILD_TYPES Debug Release MinSizeRel RelWithDebInfo Coverage)
set_property (CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${ALLOWED_BUILD_TYPES})

set (DEFAULT_BUILD_TYPE "Release")
if (EXISTS "${CMAKE_SOURCE_DIR}/.git")
  set (DEFAULT_BUILD_TYPE "Debug")
endif ()

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message (STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as no other was specified.")
  set (CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE
      STRING "Choose the type of build." FORCE)
else (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  # Ignore capitalization when build type is selected manually and check for valid setting
  string (TOLOWER ${CMAKE_BUILD_TYPE} SELECTED_LOWER)
  string (TOLOWER "${ALLOWED_BUILD_TYPES}" BUILD_TYPES_LOWER)
  if (NOT SELECTED_LOWER IN_LIST BUILD_TYPES_LOWER)
    message (FATAL_ERROR "Unknown build type: ${CMAKE_BUILD_TYPE} [allowed: ${ALLOWED_BUILD_TYPES}]")
  endif ()
  message (STATUS "Build type is: ${CMAKE_BUILD_TYPE}")
endif (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)



#
# Floating point exception support - select implementation to use
#
include (CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX ("fenv.h" HAS_FEENABLEEXCEPT)
if (HAS_FEENABLEEXCEPT) # FLOATING_POINT_ENVIRONMENT
    set (CORSIKA_HAS_FEENABLEEXCEPT 1)
    set_property (DIRECTORY ${CMAKE_HOME_DIRECTORY} APPEND PROPERTY COMPILE_DEFINITIONS "CORSIKA_HAS_FEENABLEEXCEPT")
endif ()


#
# General OS Detection
#
if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set (CORSIKA_OS_WINDOWS TRUE)
    set (CORSIKA_OS "Windows")
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set (CORSIKA_OS_LINUX TRUE)
    set (CORSIKA_OS "Linux")
    # check for RedHat/CentOS compiler from the software collections (scl)
    string (FIND ${CMAKE_CXX_COMPILER} "/opt/rh/devtoolset-" index)
    if (${index} EQUAL 0)
      set (CORSIKA_SCL_CXX TRUE)
    endif ()
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set (CORSIKA_OS_MAC TRUE)  
    set (CORSIKA_OS "Mac") 
endif ()
