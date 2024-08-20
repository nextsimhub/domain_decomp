/*!
 * @file test_grid_0.cpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @date 20 August 2024
 */

#include "Grid.hpp"

#include <catch2/catch_all.hpp>
#include <mpi.h>

TEST_CASE("Grid: non-default dimension naming", "[grid]")
{
  // Build grid from netCDF file
  Grid* grid
      = Grid::create(MPI_COMM_WORLD, "./test_2.nc", "m", "n", "land_mask");

  int mpi_size, mpi_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  assert(mpi_size <= 2);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  if (grid->get_blk_factor_0() == 1 && grid->get_blk_factor_1() == 1) {
    REQUIRE(grid->get_global_ext_0() == 6);
    REQUIRE(grid->get_global_ext_1() == 4);
    if (mpi_size == 1) {
      REQUIRE(grid->get_num_objects() == 24);
      REQUIRE(grid->get_num_nonzero_objects() == 12);
    } else if (mpi_size == 2) {
      REQUIRE(grid->get_num_objects() == 12);
      REQUIRE(grid->get_num_nonzero_objects() == 6);
    }
  } else if (grid->get_blk_factor_0() == 2 && grid->get_blk_factor_1() == 2) {
    REQUIRE(grid->get_global_ext_0() == 3);
    REQUIRE(grid->get_global_ext_1() == 2);
    if (mpi_size == 1) {
      REQUIRE(grid->get_num_objects() == 6);
      REQUIRE(grid->get_num_nonzero_objects() == 6);
    } else if (mpi_size == 2) {
      if (mpi_rank == 0) {
        REQUIRE(grid->get_num_objects() == 4);
        REQUIRE(grid->get_num_nonzero_objects() == 4);
      } else {
        REQUIRE(grid->get_num_objects() == 2);
        REQUIRE(grid->get_num_nonzero_objects() == 2);
      }
    }
  }

  const int* mask = grid->get_land_mask();
  for (int i = 0; i < grid->get_num_objects(); i++) {
    if (grid->get_blk_factor_0() == 1 && grid->get_blk_factor_1() == 1) {
      if ((i % grid->get_global_ext_1()) % 2 == 0) {
        REQUIRE(mask[i] == 0);
      }
    } else if (grid->get_blk_factor_0() == 2 && grid->get_blk_factor_1() == 2) {
      REQUIRE(mask[i] == 1);
    }
  }

  // Cleanup
  delete grid;
}
