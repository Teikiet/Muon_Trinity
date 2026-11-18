set (input_dir ${PROJECT_SOURCE_DIR}/src/media)
set (output_dir ${PROJECT_BINARY_DIR}/corsika/media)

file (MAKE_DIRECTORY ${output_dir})

add_custom_command (       
  OUTPUT  ${output_dir}/GeneratedMediaProperties.inc     
  COMMAND ${input_dir}/readProperties.py ${input_dir}/properties8.dat                                       
  DEPENDS ${input_dir}/readProperties.py
          ${input_dir}/properties8.dat         
  WORKING_DIRECTORY
          ${output_dir}
  COMMENT "Read NIST properties8 data file and produce C++ source code GeneratedMediaProperties.inc"
  VERBATIM
  )

set_source_files_properties (
  ${output_dir}/GeneratedMediaProperties.inc
  PROPERTIES GENERATED TRUE
)

add_custom_target (GenMediaProperties DEPENDS ${output_dir}/GeneratedMediaProperties.inc)
add_dependencies (CORSIKA8 GenMediaProperties)

install (
  FILES ${output_dir}/GeneratedMediaProperties.inc
  DESTINATION include/corsika/media
  )
