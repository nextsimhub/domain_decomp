/*!
 * @file test_grid_0.cpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @date 23 Aug 2024
 */

#include "Grid.hpp"

#include <catch2/catch_all.hpp>
#include <mpi.h>

TEST_CASE("Grid: no land", "[grid]")
{
    // Build grid from netCDF file
    Grid* grid = Grid::create(MPI_COMM_WORLD, "./test_1.nc");

    int mpi_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    assert(mpi_size <= 2);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    REQUIRE(grid->get_global_ext_0() == 6);
    REQUIRE(grid->get_global_ext_1() == 4);
    if (mpi_size == 1) {
        REQUIRE(grid->get_num_objects() == 24);
        REQUIRE(grid->get_num_nonzero_objects() == 24);
    } else if (mpi_size == 2) {
        REQUIRE(grid->get_num_objects() == 12);
        REQUIRE(grid->get_num_nonzero_objects() == 12);
    }

    const int* mask = grid->get_land_mask();
    for (int i = 0; i < grid->get_num_objects(); i++) {
        REQUIRE(mask[i] == 1);
    }

    // Cleanup
    delete grid;
}
