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

  add_executable(${ARG_TARGET} ${ARG_SOURCES} MainMPI.cpp)
  target_link_libraries(${ARG_TARGET}
    PRIVATE doctest::doctest ${PROJECT_NAME}::${LIB_NAME}
    )
  target_include_directories(${ARG_TARGET}
    PRIVATE ${CMAKE_SOURCE_DIR} ${DOCTEST_INCLUDE_DIR}
    )
endfunction()

# Add tests
create_test(TARGET test_grid_0_1 SOURCES test_grid_0_1.cpp)
create_test(TARGET test_grid_0_2 SOURCES test_grid_0_2.cpp)
create_test(TARGET test_grid_1_1 SOURCES test_grid_1_1.cpp)
create_test(TARGET test_grid_1_2 SOURCES test_grid_1_2.cpp)
create_test(TARGET test_grid_2_1 SOURCES test_grid_2_1.cpp)
create_test(TARGET test_grid_2_2 SOURCES test_grid_2_2.cpp)
create_test(TARGET test_zoltan_0_1 SOURCES test_zoltan_partitioner_0_1.cpp)
create_test(TARGET test_zoltan_0_2 SOURCES test_zoltan_partitioner_0_2.cpp)
create_test(TARGET test_zoltan_1_1 SOURCES test_zoltan_partitioner_1_1.cpp)
create_test(TARGET test_zoltan_1_2 SOURCES test_zoltan_partitioner_1_2.cpp)
create_test(TARGET test_zoltan_2_1 SOURCES test_zoltan_partitioner_2_1.cpp)
create_test(TARGET test_zoltan_2_2 SOURCES test_zoltan_partitioner_2_2.cpp)
create_test(TARGET test_zoltan_2_3 SOURCES test_zoltan_partitioner_2_3.cpp)
create_test(TARGET test_zoltan_2_4 SOURCES test_zoltan_partitioner_2_4.cpp)

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
