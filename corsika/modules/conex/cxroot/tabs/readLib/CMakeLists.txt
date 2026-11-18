set (
    files
    source/Interface.bz2.cc 
    )

add_library (CorsikaData STATIC ${files})
find_package(BZip2 REQUIRED)
target_link_libraries (CorsikaData
  PUBLIC
  BZip2::BZip2
  )

set_target_properties (
  CorsikaData
  PROPERTIES
  PUBLIC_HEADER corsika_data/Interface.h
  POSITION_INDEPENDENT_CODE 1
  )

target_include_directories (
  CorsikaData
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/externals/compress_stream/include>
  $<INSTALL_INTERFACE:include/corsika_modules/data>
  )

if (TARGET CORSIKA8)
  install (
    TARGETS CorsikaData
    EXPORT CORSIKA8PublicTargets
    LIBRARY DESTINATION lib/corsika
    ARCHIVE DESTINATION lib/corsika
    PUBLIC_HEADER DESTINATION include/corsika_modules/data
    )
else (TARGET CORSIKA8)
  install (
    TARGETS CORSIKA8PublicTargets
    LIBRARY DESTINATION lib/corsika
    ARCHIVE DESTINATION lib/corsika
    PUBLIC_HEADER DESTINATION include/corsika_modules/data
    )
endif (TARGET CORSIKA8)

# add unit test, if run inside CORSIKA
if (COMMAND CORSIKA_ADD_TEST)
  if (CORSIKA_DATA_WITH_TEST)
    message (STATUS "CorsikaData uses CORSIKA_ADD_TEST")
    CORSIKA_ADD_TEST (testData SOURCES source/testData.cc)
    target_link_libraries (testData PRIVATE CorsikaData)
    target_compile_definitions (
      testData
      PRIVATE
      STANDALONE=1
      TESTDATA="${CMAKE_CURRENT_SOURCE_DIR}/source/test.data"
      )
  endif (CORSIKA_DATA_WITH_TEST)
else (COMMAND CORSIKA_ADD_TEST)
  message (STATUS "CorsikaData did not detector CORSIKA8")
  add_executable (testData source/testData.cc)
  target_link_libraries (testData PRIVATE CorsikaData Catch2::Catch2WithMain)
  target_compile_options (testData PRIVATE -g) # do not skip asserts
  target_include_directories (testData PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_compile_definitions (
     testData
     PRIVATE
     STANDALONE=1
     TESTDATA="${CMAKE_CURRENT_SOURCE_DIR}/source/test.data"
     )
  file (MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/test_outputs/)
  add_test (NAME testData COMMAND testData -o ${PROJECT_BINARY_DIR}/test_outputs/junit-testData.xml -s -r junit)
endif (COMMAND CORSIKA_ADD_TEST)
