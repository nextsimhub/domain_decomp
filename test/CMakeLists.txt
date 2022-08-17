function(add_my_test)
  # Define the supported set of keywords
  set(prefix ARG)
  set(noValueArgs "")
  set(oneValueArgs TARGET MAX_NUM_PROCS)
  set(multiValueArgs SOURCES)

  # Process the arguments passed in
  cmake_parse_arguments(${prefix}
    "${noValueArgs}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN})

  add_executable(${ARG_TARGET}
    ${CMAKE_SOURCE_DIR}/Grid.cpp
    ${CMAKE_SOURCE_DIR}/Partitioner.cpp
    ${CMAKE_SOURCE_DIR}/ZoltanPartitioner.cpp
    ${ARG_SOURCES}
    main.cpp)
  target_include_directories(${ARG_TARGET} PRIVATE ${CMAKE_SOURCE_DIR} ${netCDF_INCLUDE_DIR} ${Zoltan_INCLUDE_DIRS})
  target_link_directories(${ARG_TARGET} PRIVATE ${netCDF_LIB_DIR})
  target_link_libraries(${ARG_TARGET} PRIVATE MPI::MPI_CXX netcdf netcdf-cxx4 ${Zoltan_LIBRARIES} Boost::program_options Catch2::Catch2)
  foreach(numProcs RANGE 1 ${ARG_MAX_NUM_PROCS})
    add_test(NAME ${ARG_TARGET}_n${numProcs} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${numProcs} ${MPIEXEC_PREFLAGS} ./${ARG_TARGET})
  endforeach()
endfunction()

# Add tests
add_my_test(TARGET test_grid_0 SOURCES test_grid_0.cpp MAX_NUM_PROCS 2)
add_my_test(TARGET test_grid_1 SOURCES test_grid_1.cpp MAX_NUM_PROCS 2)
add_my_test(TARGET test_grid_2 SOURCES test_grid_2.cpp MAX_NUM_PROCS 2)
add_my_test(TARGET test_zoltan_0 SOURCES test_zoltan_partitioner_0.cpp MAX_NUM_PROCS 2)
add_my_test(TARGET test_zoltan_1 SOURCES test_zoltan_partitioner_1.cpp MAX_NUM_PROCS 2)
add_my_test(TARGET test_zoltan_2 SOURCES test_zoltan_partitioner_2.cpp MAX_NUM_PROCS 4)

# Copy grid files used for testing
file(COPY test_0.nc DESTINATION ${CMAKE_BINARY_DIR}/test)
file(COPY test_1.nc DESTINATION ${CMAKE_BINARY_DIR}/test)
file(COPY test_2.nc DESTINATION ${CMAKE_BINARY_DIR}/test)
