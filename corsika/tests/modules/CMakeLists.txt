set (test_modules_sources
  TestMain.cpp
  testStackInspector.cpp
  testTracking.cpp
  ## testExecTime.cpp --> needs to be fixed, see #326
  testObservationPlane.cpp
  testObservationVolume.cpp
  testQGSJetII.cpp
  testQGSJetIII.cpp
  testProposal.cpp
  testPythia8.cpp
  testTAUOLA.cpp
  testUrQMD.cpp
  # testCONEX.cpp # TODO: re-enable when Sibyll initialization is fixed
  testParticleCut.cpp
  testProductionProfile.cpp
  testSibyll.cpp
  testEpos.cpp
  testEposLhcr.cpp
  testRadio.cpp
  testEMThinning.cpp
  testSophia.cpp
  testCherenkovRefractiveIndexTable.cpp
  testCherenkovCurvedRayTable.cpp
  #testRandom.cpp
)

if (NOT "$CACHE{C8_FLUKA_ENABLED}" STREQUAL "")
    list(APPEND test_modules_sources "testFluka.cpp")
    message("testFluka enabled")
endif()

CORSIKA_ADD_TEST (testModules SOURCES ${test_modules_sources})

# Separate target: the ON/OFF transmission measurement must not trigger a full
# rebuild of all of testModules (QGSJet, Sibyll, Epos, Pythia, ...).
CORSIKA_ADD_TEST (testCherenkovE2E SOURCES testCherenkovE2E.cpp)

# Chord-route geometry grid against the independent spherical reference (2g-h2g).
CORSIKA_ADD_TEST (testCherenkovChordRoute SOURCES testCherenkovChordRoute.cpp)

# Curved tracer (no table) against the refracted reference (2g-h2h).
CORSIKA_ADD_TEST (testCherenkovCurvedTracerGrid SOURCES testCherenkovCurvedTracerGrid.cpp)

# Step 3: table prototype at 3 km detector, 90 deg pointing, 10 deg full FoV (2g-h2i).
CORSIKA_ADD_TEST (testCherenkovTable3km SOURCES testCherenkovTable3km.cpp)

# Trace fix track (2g-h2k): frozen acceptance set for the legacy curved-ray tracer.
CORSIKA_ADD_TEST (testCherenkovTracerAcceptance SOURCES testCherenkovTracerAcceptance.cpp)

# Telescope-mode track (2g-h2l): frozen acceptance set for the plane tracer at
# 90 deg pointing (the mode had no test before this).
CORSIKA_ADD_TEST (testCherenkovPlaneTracerAcceptance SOURCES testCherenkovPlaneTracerAcceptance.cpp)

# Wave 2b (2g-h2n): the exact chord optical depth of the absorption table, and the
# measured size of the Chapman air-mass bias it replaced.
CORSIKA_ADD_TEST (testCherenkovAbsorptionChord SOURCES testCherenkovAbsorptionChord.cpp)

# Phase 3: RK4 (reference/fallback) vs the optional adaptive DOPRI5 5(4) on a frozen ray
# set, against a converged reference.  Reports cost, asserts accuracy.
CORSIKA_ADD_TEST (testCherenkovTracerBenchmark SOURCES testCherenkovTracerBenchmark.cpp)

# Phase 4, step 1: the regularised (turning-point) tracer against regular DOPRI5.
CORSIKA_ADD_TEST (testCherenkovTangentTracer SOURCES testCherenkovTangentTracer.cpp)

# Phase 4, step 3: the impact-aware accuracy controller, swept over plane orientations
# (the shallow rays are the ones the plain tolerance rule cannot see).
CORSIKA_ADD_TEST (testCherenkovImpactController SOURCES testCherenkovImpactController.cpp)

target_compile_definitions (
  testModules
  PRIVATE
  REFDATADIR="${CMAKE_CURRENT_SOURCE_DIR}"
)
