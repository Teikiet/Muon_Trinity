# Cherenkov module - header only

#set(
#  HEADERS
#  Cherenkov.hpp
#  Cherenkov.inl
#)

#add_library(Cherenkov INTERFACE)

#target_link_libraries(
#  Cherenkov
#  INTERFACE
#    CorsikaFramework::Core
#    CorsikaFramework::Process
#)

#target_include_directories(
#  Cherenkov
#  INTERFACE
#    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
#    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
#    $<INSTALL_INTERFACE:include/corsika>
#)

# Optional, but nice to have: official module alias
#add_library(CorsikaModules::Cherenkov ALIAS Cherenkov)

#install(
#  FILES ${HEADERS}
#  DESTINATION include/corsika/modules/cherenkov
#)
