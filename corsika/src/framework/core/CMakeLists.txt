set (input_dir ${PROJECT_SOURCE_DIR}/src/framework/core)
set (output_dir ${PROJECT_BINARY_DIR}/corsika/framework/core)

file (MAKE_DIRECTORY ${output_dir})

add_custom_command (
  OUTPUT  ${output_dir}/GeneratedParticleProperties.inc
          ${output_dir}/GeneratedParticleClasses.inc
          ${output_dir}/particle_db.pkl
  COMMAND ${input_dir}/make_input_file.py
  DEPENDS ${input_dir}/make_input_file.py
          ${input_dir}/PDG_Ids.txt
          ${input_dir}/NuclearData.txt
  WORKING_DIRECTORY ${output_dir}
  COMMENT "Get pdg data and produce C++ source code GeneratedParticle[...].inc"
  VERBATIM
  )

set_source_files_properties (
  ${output_dir}/GeneratedParticleProperties.inc
  ${output_dir}/GeneratedParticleClasses.inc
  ${output_dir}/particle_db.pkl
  PROPERTIES GENERATED TRUE
  )

add_custom_target (GenParticlesHeaders
  DEPENDS ${output_dir}/GeneratedParticleProperties.inc
          ${output_dir}/GeneratedParticleClasses.inc
          ${output_dir}/particle_db.pkl
  )
add_dependencies (CORSIKA8 GenParticlesHeaders)

install (
  FILES
  ${output_dir}/GeneratedParticleProperties.inc
  ${output_dir}/GeneratedParticleClasses.inc
  DESTINATION include/corsika/framework/core
)
