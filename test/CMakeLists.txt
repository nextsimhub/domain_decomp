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

# Generate input files from respective cdl files
add_custom_command(
    OUTPUT test_0.nc
    COMMAND
        ncgen -b -o test_0.nc ${CMAKE_CURRENT_SOURCE_DIR}/test_0.cdl
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/test_0.cdl
)
add_custom_command(
    OUTPUT test_1.nc
    COMMAND
        ncgen -b -o test_1.nc ${CMAKE_CURRENT_SOURCE_DIR}/test_1.cdl
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/test_1.cdl
)
add_custom_command(
    OUTPUT test_2.nc
    COMMAND
        ncgen -b -o test_2.nc ${CMAKE_CURRENT_SOURCE_DIR}/test_2.cdl
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/test_2.cdl
)
add_custom_target(
    generate_input_files
    ALL
    DEPENDS test_0.nc test_1.nc test_2.nc
)
