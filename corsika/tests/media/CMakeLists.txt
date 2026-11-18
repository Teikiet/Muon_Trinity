set (test_media_sources
  TestMain.cpp
  testNuclearComposition.cpp
  testEnvironment.cpp
  testShowerAxis.cpp
  testMedium.cpp
  refractivity/testExponential.cpp
  refractivity/testGladstoneDale.cpp
  refractivity/testLinearTabulated.cpp
  refractivity/testUniform.cpp
  testMagneticField.cpp
  testCORSIKA7Atmospheres.cpp
  atmosphere_table/testAtmosphereTable.cpp
  )

CORSIKA_ADD_TEST (testMedia SOURCES ${test_media_sources})
find_package(Boost REQUIRED)
target_link_libraries (testMedia PRIVATE Boost::boost) # for math/quadrature/gauss_kronrod.hpp

target_compile_definitions (
  testMedia
  PRIVATE
  REFDATADIR="${CMAKE_CURRENT_SOURCE_DIR}"
)
