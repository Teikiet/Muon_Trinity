set (
  CNPY_SOURCES
  cnpy.cpp
  )

set (
  CNPY_HEADERS
  cnpy.hpp
  )

add_library (cnpy STATIC ${CNPY_SOURCES})
find_package(ZLIB REQUIRED)
# target dependencies on other libraries (also the header onlys)
target_link_libraries (
  cnpy
  PUBLIC
  ZLIB::ZLIB
  )

target_include_directories (
  cnpy 
  INTERFACE 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:corsika_external/cnpy/include>
  )

install (
  FILES cnpy.hpp
  DESTINATION corsika_external/cnpy/include
  )

install (
  TARGETS cnpy
  EXPORT CORSIKA8PublicTargets
  LIBRARY DESTINATION corsika_external/cnpy/lib
  ARCHIVE DESTINATION corsika_external/cnpy/lib
  )

