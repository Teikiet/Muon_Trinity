add_executable (c8_air_shower c8_air_shower.cpp)
target_link_libraries (c8_air_shower PRIVATE CORSIKA8)

set_target_properties(c8_air_shower PROPERTIES
  # use relative path that also works when cmake-install is invoked with --prefix
  INSTALL_RPATH "$ORIGIN/../lib/corsika"
)

if(WITH_FLUKA)
  message("compiling c8_air_shower.cpp with FLUKA")
  target_compile_definitions(c8_air_shower PRIVATE WITH_FLUKA)
else()
  message("compiling c8_air_shower.cpp with UrQMD")
endif()

install (
  TARGETS c8_air_shower DESTINATION bin
)

# --- Add test_cherenkov ---

#add_executable(test_cherenkov test_cherenkov.cpp)

# Link against the main CORSIKA8 target AND the Cherenkov module
#target_link_libraries(test_cherenkov 
#  PRIVATE 
#    CORSIKA8
#)

# Set runtime paths so it finds the libraries
#set_target_properties(test_cherenkov PROPERTIES
#  INSTALL_RPATH "$ORIGIN/../lib/corsika"
#)

# Allow it to be installed to the bin folder
#install(TARGETS test_cherenkov DESTINATION bin)

# --- Add test_muons ---

add_executable(muon_cherenkov muon_cherenkov.cpp)

# Link against the main CORSIKA8 target AND the Cherenkov module
target_link_libraries(muon_cherenkov 
  PRIVATE 
    CORSIKA8
)

# Set runtime paths so it finds the libraries
set_target_properties(muon_cherenkov PROPERTIES
  INSTALL_RPATH "$ORIGIN/../lib/corsika"
)

# Allow it to be installed to the bin folder
install(TARGETS muon_cherenkov DESTINATION bin)