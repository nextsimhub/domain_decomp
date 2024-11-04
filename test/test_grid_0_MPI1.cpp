/*!
 * @file test_grid_0_1.cpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @date 04 Nov 2024
 */
#include <doctest/extensions/doctest_mpi.h>

#include "Grid.hpp"

MPI_TEST_CASE("Grid: all land, 1 MPI rank", 1)
{
    // Build grid from netCDF file
    Grid* grid = Grid::create(MPI_COMM_WORLD, "./test_0.nc");

    REQUIRE(grid->get_global_ext()[0] == 6);
    REQUIRE(grid->get_global_ext()[1] == 4);
    REQUIRE(grid->get_num_nonzero_objects() == 0);
    REQUIRE(grid->get_num_objects() == 24);

    const int* mask = grid->get_land_mask();
    for (int i = 0; i < grid->get_num_objects(); i++) {
        REQUIRE(mask[i] == 0);
    }

    // Cleanup
    delete grid;
}
