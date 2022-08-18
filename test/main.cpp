/*!
 * @file main.cpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @date 13 August 2022
 */

#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>
#include <mpi.h>

// Global options
int global_argc;
char** global_argv;

int main(int argc, char* argv[])
{
  // Initialize MPI
  MPI_Init(&argc, &argv);

  // Assign globals, so they can be used within test cases
  global_argc = argc;
  global_argv = argv;

  int result = Catch::Session().run(argc, argv);

  // Finalize MPI
  MPI_Finalize();

  return result;
}
