enable_language (Fortran)
set (CMAKE_Fortran_FLAGS "-I${CMAKE_CURRENT_SOURCE_DIR}/src/ -I${CMAKE_CURRENT_BINARY_DIR}/src/ -I${CMAKE_CURRENT_SOURCE_DIR}/src/urqmd13/")
message (STATUS "***** Configuring CxRoot/CONEX version")


set (MODEL_HEADERS
  src/conexConfig.h
  src/conexHEModels.h
  )


#+++++++++++++++++++++++++++++
# Location, saved in conexConfigAuto.h
#
if (TARGET CorsikaData)
   set (TABDIR "${CorsikaData_DIR}/")   
else (TARGET CorsikaData)
  set (TABDIR_PRE "${CMAKE_CURRENT_SOURCE_DIR}/tabs")
  if (DEFINED ENV{CORSIKA_DATA})
    set (TABDIR_PRE "$ENV{CORSIKA_DATA}/CONEX")
  else ()
    if (DEFINED ENV{CONEXTABROOT})
      set (TABDIR_PRE "$ENV{CONEXTABROOT}")
    endif ()
  endif ()
  set (TABDIR "${TABDIR_PRE}/")
  add_subdirectory (${TABDIR} tabs)
endif (TARGET CorsikaData)

# we need two version of conexConfigAuto.h. One for the build and one for install area.
# build version:
set (CFG_CONEX_CONFIG_PATH "${CMAKE_CURRENT_BINARY_DIR}/cfg/")
configure_file (src/conexConfigAuto.h.in src/conexConfigAuto.h) 
# install version:
set (CFG_CONEX_CONFIG_PATH "${CMAKE_INSTALL_PREFIX}/share/corsika/cfg/")
configure_file (src/conexConfigAuto.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/src/conexConfigAuto_install.h
  )
install (
  FILES ${CMAKE_CURRENT_BINARY_DIR}/src/conexConfigAuto_install.h 
  RENAME conexConfigAuto.h 
  DESTINATION include/corsika_modules/conex/
  )


#+++++++++++++++++++++++++++++
# Config file generation
# Note: if we also copy/relocate TABDIR we have to rething this, too:
#
file (GLOB CONFIG_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.paramin")

foreach (CONFIG_FILE_ABS ${CONFIG_FILES})
  get_filename_component (CONFIG_FILE ${CONFIG_FILE_ABS} NAME)
  string (REPLACE paramin param  CONFIG_OUT ${CONFIG_FILE}) 
  message (STATUS "==[make]==> generating cfg/${CONFIG_OUT} from : src/${CONFIG_FILE}")
  set (CONEXTABROOT ${TABDIR})
  set (BZ2 ".bz2")
  configure_file (src/${CONFIG_FILE} ${CMAKE_CURRENT_BINARY_DIR}/cfg/${CONFIG_OUT})
endforeach ()


#+++++++++++++++++++++++++++++
# Model selection, and also library generation
#
set (ALL_MODELS
  epos
  qgsjet
  qgsjetII
  qgsjetIII
  dpm
  sibyll
  )

# we only build and use SIBYLL now, there is not yet
# a mechanism to dynamically change that
set (MODELS
  sibyll
  )

# Our simplified "dict":
set (DEF_epos __EPOS__)
set (DEF_sibyll __SIBYLL21__)
set (DEF_qgsjetII __QGSJETII__)
set (DEF_qgsjetIII __QGSJETIII__)
set (DEF_qII __QGSJETII__)
set (DEF_qIII __QGSJETIII__)
set (DEF_dpmjet __DPMJET__)
set (DEF_qgsjet __QGSJET__)

# Another simplified "dict":
set (OBJ_sibyll SIBYLL23d_static)

set (CONEX_LIB_SOURCES
  src/conex_sub.F
  src/conex_mod.F
  src/conex2r.F
  src/leadingInteractionsData.cc
  src/conexConfigAuto.h
  )

foreach (MODEL ${MODELS})
  message (STATUS "adding: libCONEX${MODEL}")
  add_library (CONEX${MODEL} SHARED ${CONEX_LIB_SOURCES})
  set (MODEL_DEF DEF_${MODEL})
  set_target_properties (
    CONEX${MODEL}
    PROPERTIES
    COMPILE_DEFINITIONS "${${MODEL_DEF}};CX_NO_ROOT;CORSIKA_8;__CORSIKA8__;__URQMD__;__CXSUB__;__SORT_FOR_ENERGY__;CORSIKA_DATA=${CORSIKA_DATA}"
    POSITION_INDEPENDENT_CODE 1
    )

  target_compile_options (
    CONEX${MODEL}
    PUBLIC
    $<$<COMPILE_LANGUAGE:Fortran>:-std=legacy;-Wfunction-elimination;-w;-fno-second-underscore;-fno-automatic;-Wno-conversion>
    )

  target_include_directories (
    CONEX${MODEL}
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/src>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include/corsika_modules/conex>
    )

  if (NOT TARGET ${OBJ_${MODEL}})
    message(FATAL_ERROR "cxroot link: target \"${OBJ_${MODEL}}\" does not exist")
  endif()

  # make sure to link to actual model library
  target_link_libraries (
    CONEX${MODEL} PRIVATE
    ${OBJ_${MODEL}}
    UrQMD_static      # AND urqmd
    #EPOS_LHC_static
    )

  #
  # conex also read tables from corsika-data
  #
  if (TARGET CorsikaData)
    target_link_libraries (
      CONEX${MODEL} PRIVATE CorsikaData
      ) 
  endif (TARGET CorsikaData)

  #
  # add to combined CORSIKA8 target
  #
  add_dependencies (CORSIKA8 CONEX${MODEL})
  target_link_libraries (CORSIKA8 INTERFACE CONEX${MODEL})

  #
  # install and export
  #  
  install (
    TARGETS CONEX${MODEL}
    EXPORT CORSIKA8PublicTargets
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/corsika"
    PUBLIC_HEADER DESTINATION include/corsika_modules/conex
    )
  
endforeach ()

# Sibyll has a variable called NDEBUG, which clashes with the "-DNDEBUG" compile option
# that is set automatically by CMake for non-debug builds. Let's undefine it again here..
target_compile_options(CONEXsibyll PRIVATE "$<$<COMPILE_LANGUAGE:Fortran>:-UNDEBUG>")


install (
  DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cfg
  DESTINATION share/corsika
  )
