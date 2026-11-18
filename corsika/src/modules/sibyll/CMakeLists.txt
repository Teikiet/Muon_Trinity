set (input_dir ${PROJECT_SOURCE_DIR}/src/modules/sibyll)
set (output_dir ${PROJECT_BINARY_DIR}/corsika/modules/sibyll)

file (MAKE_DIRECTORY ${output_dir})

add_custom_command (
  OUTPUT  ${output_dir}/Generated.inc
  COMMAND ${input_dir}/code_generator.py 
          ${PROJECT_BINARY_DIR}/corsika/framework/core/particle_db.pkl
          ${input_dir}/sibyll_codes.dat
  DEPENDS ${input_dir}/code_generator.py
          ${input_dir}/sibyll_codes.dat
          ${PROJECT_BINARY_DIR}/corsika/framework/core/particle_db.pkl
          GenParticlesHeaders
  WORKING_DIRECTORY ${output_dir}
  COMMENT "Generate conversion tables for particle codes SIBYLL <-> CORSIKA"
  VERBATIM
  )

set_source_files_properties (
  ${output_dir}/Generated.inc
  PROPERTIES GENERATED TRUE
)

add_custom_target (SourceDirLinkSib DEPENDS ${output_dir}/Generated.inc)
add_dependencies (CORSIKA8 SourceDirLinkSib)

install (
  FILES ${output_dir}/Generated.inc
  DESTINATION include/corsika/modules/sibyll
)
