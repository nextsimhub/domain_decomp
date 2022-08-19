function(create_test)
  # Define the supported set of keywords
  set(prefix ARG)
  set(noValueKeywords "")
  set(oneValueKeywords TARGET MAX_NUM_PROCS)
  set(multiValueKeywords SOURCES)

  # Process the arguments passed in
  cmake_parse_arguments(${prefix}
    "${noValueKeywords}"
    "${oneValueKeywords}"
    "${multiValueKeywords}"
    ${ARGN}
    )

  add_executable(${ARG_TARGET}
    ${ARG_SOURCES}
    main.cpp
    )
  target_link_libraries(${ARG_TARGET}
    PRIVATE Catch2::Catch2 ${PROJECT_NAME}::${LIB_NAME}
    )
  target_include_directories(${ARG_TARGET}
    PRIVATE ${CMAKE_SOURCE_DIR}
    )
  foreach(numProcs RANGE 1 ${ARG_MAX_NUM_PROCS})
    add_test(NAME ${ARG_TARGET}_n${numProcs} COMMAND ${MPIEXEC} --oversubscribe ${MPIEXEC_NUMPROC_FLAG} ${numProcs} ${MPIEXEC_PREFLAGS} ./${ARG_TARGET})
  endforeach()
endfunction()

# Add tests
create_test(TARGET test_grid_0 SOURCES test_grid_0.cpp MAX_NUM_PROCS 2)
create_test(TARGET test_grid_1 SOURCES test_grid_1.cpp MAX_NUM_PROCS 2)
create_test(TARGET test_grid_2 SOURCES test_grid_2.cpp MAX_NUM_PROCS 2)
create_test(TARGET test_zoltan_0 SOURCES test_zoltan_partitioner_0.cpp MAX_NUM_PROCS 2)
create_test(TARGET test_zoltan_1 SOURCES test_zoltan_partitioner_1.cpp MAX_NUM_PROCS 2)
create_test(TARGET test_zoltan_2 SOURCES test_zoltan_partitioner_2.cpp MAX_NUM_PROCS 4)

# Copy grid files used for testing
file(COPY test_0.nc DESTINATION ${CMAKE_BINARY_DIR}/test)
file(COPY test_1.nc DESTINATION ${CMAKE_BINARY_DIR}/test)
file(COPY test_2.nc DESTINATION ${CMAKE_BINARY_DIR}/test)
