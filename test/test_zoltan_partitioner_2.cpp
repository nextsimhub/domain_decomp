/*!
 * @file test_zoltan_partitioner_0.cpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @date 17 August 2022
 */

#include "Grid.hpp"
#include "Partitioner.hpp"

#include <catch2/catch.hpp>
#include <mpi.h>

extern int global_argc;
extern char** global_argv;

TEST_CASE("ZoltanPartitioner: non-default dimension naming", "[zoltan]")
{
  // Build grid from netCDF file
  Grid* grid
      = Grid::create(MPI_COMM_WORLD, "./test_2.nc", "m", "n", "land_mask");

  // Create a Zoltan partitioner
  Partitioner* partitioner = Partitioner::Factory::create(
      MPI_COMM_WORLD, global_argc, global_argv, PartitionerType::Zoltan_RCB);

  // Partition grid
  partitioner->partition(*grid);

  int mpi_size, mpi_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  assert(mpi_size <= 4);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  if (mpi_size == 1) {
    REQUIRE(grid->get_local_ext_0() == 6);
    REQUIRE(grid->get_local_ext_1() == 4);
    REQUIRE(grid->get_global_0() == 0);
    REQUIRE(grid->get_global_1() == 0);
  } else if (mpi_size == 2) {
    REQUIRE(grid->get_local_ext_0() == 3);
    REQUIRE(grid->get_local_ext_1() == 4);
    if (mpi_rank == 0) {
      REQUIRE(grid->get_global_0() == 0);
      REQUIRE(grid->get_global_1() == 0);
    } else {
      REQUIRE(grid->get_global_0() == 3);
      REQUIRE(grid->get_global_1() == 0);
    }
  } else if (mpi_size == 3) {
    REQUIRE(grid->get_local_ext_0() == 2);
    REQUIRE(grid->get_local_ext_1() == 4);
    if (mpi_rank == 0) {
      REQUIRE(grid->get_global_0() == 0);
      REQUIRE(grid->get_global_1() == 0);
    } else if (mpi_rank == 1) {
      REQUIRE(grid->get_global_0() == 2);
      REQUIRE(grid->get_global_1() == 0);
    } else {
      REQUIRE(grid->get_global_0() == 4);
      REQUIRE(grid->get_global_1() == 0);
    }
  } else {
    if (mpi_rank == 0) {
      REQUIRE(grid->get_local_ext_0() == 1);
      REQUIRE(grid->get_local_ext_1() == 4);
      REQUIRE(grid->get_global_0() == 0);
      REQUIRE(grid->get_global_1() == 0);
    } else if (mpi_rank == 1) {
      REQUIRE(grid->get_local_ext_0() == 2);
      REQUIRE(grid->get_local_ext_1() == 4);
      REQUIRE(grid->get_global_0() == 1);
      REQUIRE(grid->get_global_1() == 0);
    } else if (mpi_rank == 2) {
      REQUIRE(grid->get_local_ext_0() == 1);
      REQUIRE(grid->get_local_ext_1() == 4);
      REQUIRE(grid->get_global_0() == 3);
      REQUIRE(grid->get_global_1() == 0);
    } else {
      REQUIRE(grid->get_local_ext_0() == 2);
      REQUIRE(grid->get_local_ext_1() == 4);
      REQUIRE(grid->get_global_0() == 4);
      REQUIRE(grid->get_global_1() == 0);
    }
  }

  // Cleanup
  delete grid;
  delete partitioner;
}
