set (input_dir ${PROJECT_SOURCE_DIR}/src/modules/fluka)
set (output_dir ${PROJECT_BINARY_DIR}/corsika/modules/fluka)

file (MAKE_DIRECTORY ${output_dir})

add_custom_command (
  OUTPUT  ${output_dir}/Generated.inc
  COMMAND ${input_dir}/code_generator.py 
          ${PROJECT_BINARY_DIR}/corsika/framework/core/particle_db.pkl
          ${input_dir}/fluka_codes.dat
  DEPENDS ${input_dir}/code_generator.py
          ${input_dir}/fluka_codes.dat
          GenParticlesHeaders # for particle_db.pkl
  WORKING_DIRECTORY
          ${output_dir}/
  COMMENT "Generate conversion tables for particle codes FLUKA <-> CORSIKA"
  VERBATIM
  )

add_custom_target (SourceDirLinkFLUKA DEPENDS ${output_dir}/Generated.inc)
add_dependencies (CORSIKA8 SourceDirLinkFLUKA)

install (
  FILES ${output_dir}/Generated.inc
  DESTINATION include/corsika/modules/fluka
)
