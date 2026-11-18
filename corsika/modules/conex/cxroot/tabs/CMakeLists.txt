#################################
# Distribution of partly compressed data tables+files for air shower physics
# Ralf Ulrich, 2020
#

#################################
# This CMakeLists.txt may be included in your project
# to access corsika-data by adding
# add_subdirectory (corsika_data)
# to your CMakeLists.txt
# it creates the cmake taget "CorsikaData"
# and it builds the libCorsikaData.a archive.
#

##################################
# for standalong running
#
cmake_minimum_required (VERSION 3.9)
project (CorsikaData)

###################################
# set useful varialbles upstream
#
get_directory_property (hasParent PARENT_DIRECTORY)
if (hasParent)
  set (CorsikaData_DIR ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
  set (CorsikaData_VERSION "1.1.0" PARENT_SCOPE)
endif (hasParent)

###################################
# CorsikaData requires to run with conan, typically inside CORSIKA8
#
if (NOT TARGET BZip2::BZip2)
  message (FATAL_ERROR "Please configure with conan bzip2 package")
endif (NOT TARGET BZip2::BZip2)
#
if (NOT TARGET Catch2::Catch2WithMain)
  message (FATAL_ERROR "Please configure with conan catch2 package")
endif (NOT TARGET Catch2::Catch2WithMain)

enable_testing ()
if (${CMAKE_VERSION} VERSION_LESS "3.12.0")  
  list (APPEND CMAKE_CTEST_ARGUMENTS "--output-on-failure")
else (${CMAKE_VERSION} VERSION_LESS "3.12.0")
  set (CTEST_OUTPUT_ON_FAILURE 1) # this is for new versions of cmake
endif (${CMAKE_VERSION} VERSION_LESS "3.12.0")

add_subdirectory (readLib)
