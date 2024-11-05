function(create_test)
  # Define the supported set of keywords
  set(prefix ARG)
  set(noValueKeywords "")
  set(oneValueKeywords TARGET)
  set(multiValueKeywords SOURCES)

  # Process the arguments passed in
  cmake_parse_arguments(${prefix}
    "${noValueKeywords}"
    "${oneValueKeywords}"
    "${multiValueKeywords}"
    ${ARGN}
    )

  add_executable(${ARG_TARGET} ${ARG_SOURCES} MainMPI.cpp)
  target_link_libraries(${ARG_TARGET}
    PRIVATE doctest::doctest ${PROJECT_NAME}::${LIB_NAME}
    )
  target_include_directories(${ARG_TARGET}
    PRIVATE ${CMAKE_SOURCE_DIR} ${DOCTEST_INCLUDE_DIR}
    )

  add_test(NAME ${ARG_TARGET} COMMAND ${MPIEXEC} --oversubscribe ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./${ARG_TARGET})

endfunction()

# Add tests
create_test(TARGET test_grid_0 SOURCES test_grid_0.cpp)
create_test(TARGET test_grid_1 SOURCES test_grid_1.cpp)
create_test(TARGET test_grid_2 SOURCES test_grid_2.cpp)
create_test(TARGET test_zoltan_0 SOURCES test_zoltan_partitioner_0.cpp)
create_test(TARGET test_zoltan_1 SOURCES test_zoltan_partitioner_1.cpp)
create_test(TARGET test_zoltan_2 SOURCES test_zoltan_partitioner_2.cpp)

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
