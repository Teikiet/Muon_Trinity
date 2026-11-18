set (CMAKE_BUILT_TYPE Debug)
#set (CMAKE_Fortran_FLAGS "-std=legacy -Wfunction-elimination -O0")

set (
  MODEL_HEADERS
  epos-lhcr-public.hpp
  epos-lhcr-types.hpp
)

enable_language (Fortran)
enable_language (C)

add_library (EPOS_LHCR_static STATIC
  epos-lhcr-private.cpp
  epos-af-lhcr.f   epos-fra-lhcr.f  epos-jps-lhcr.f  epos-tim-lhcr.f
  epos-app-lhcr.f  epos-hg-lhcr.f   eposm.f          eposu_no.f
  epos-bas-lhcr.f  epos-hnb-lhcr.f  epos-omg-lhcr.f  epos-uti-lhcr.f
  epos-con-lhcr.f  epos-ico-lhcr.f  epos-par-lhcr.f  epos-xan-lhcr.f
  epos-dky-lhcr.f  epos-ids-lhcr.f  epos-qsh-lhcr.f  epos-xpr-lhcr.f
  epos-dro-lhcr.f  epos-int-lhcr.f  epos-rsh-lhcr.f
  epos-ems-lhcr.f  epos_j.f         epos-sem-lhcr.f
  al.cpp           memory.c
)

target_link_libraries(EPOS_LHCR_static PUBLIC modules_common)
target_compile_options(EPOS_LHCR_static PRIVATE "-mcmodel=large")

add_library (EPOS_LHCR SHARED "epos-lhcr-public.cpp")
target_link_libraries(EPOS_LHCR PRIVATE EPOS_LHCR_static)

target_include_directories (EPOS_LHCR PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

set_target_properties (
  EPOS_LHCR_static
  PROPERTIES
  POSITION_INDEPENDENT_CODE 1
  CXX_VISIBILITY_PRESET hidden
  # unfortunately ignored by gfortran:
  Fortran_VISIBILITY_PRESET hidden
)

set_target_properties(EPOS_LHCR PROPERTIES
  FRAMEWORK TRUE
  PUBLIC_HEADER "${MODEL_HEADERS}"
)

target_link_options(EPOS_LHCR PRIVATE "LINKER:--exclude-libs,ALL") # hide symbols from static/private lib
target_include_directories (
  EPOS_LHCR
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/corsika_modules/epos-lhcr>
)

install (
  TARGETS EPOS_LHCR
  EXPORT CORSIKA8PublicTargets
#  ARCHIVE DESTINATION lib/corsika
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/corsika"
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/corsika_modules/epos-lhcr
)

# add EPOS to CORSIKA8 build
add_dependencies (CORSIKA8 EPOS_LHCR)
target_link_libraries (CORSIKA8 INTERFACE EPOS_LHCR)
