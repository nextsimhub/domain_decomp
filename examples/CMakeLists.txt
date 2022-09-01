function(create_example)
  # Define the supported set of keywords
  set(prefix ARG)
  set(noValueKeywords "")
  set(oneValueKeywords TARGET SOURCES)
  set(multiValueKeywords "")

  # Process the arguments passed in
  cmake_parse_arguments(${prefix}
    "${noValueKeywords}"
    "${oneValueKeywords}"
    "${multiValueKeywords}"
    ${ARGN}
    )

  add_executable(${ARG_TARGET}
    ${ARG_SOURCES}
    )
  target_link_libraries(${ARG_TARGET}
    PRIVATE Boost::program_options ${PROJECT_NAME}::${LIB_NAME}
    )
  target_include_directories(${ARG_TARGET}
    PRIVATE ${CMAKE_SOURCE_DIR}
    )
  install(TARGETS ${ARG_TARGET} DESTINATION ${CMAKE_INSTALL_BINDIR}/examples)
endfunction()

# Add examples
create_example(TARGET zoltan_comm SOURCES zoltan_comm.cpp)
