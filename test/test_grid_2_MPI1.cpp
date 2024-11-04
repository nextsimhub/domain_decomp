/*!
 * @file test_grid_2_1.cpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @date 04 Nov 2024
 */

#include "Grid.hpp"

#include <doctest/extensions/doctest_mpi.h>

MPI_TEST_CASE("Grid: non-default dimension naming, 1 MPI rank", 1)
{
    // Build grid from netCDF file
    Grid* grid = Grid::create(MPI_COMM_WORLD, "./test_2.nc", "m", "n", { 1, 0 }, "land_mask");

    REQUIRE(grid->get_global_ext()[0] == 6);
    REQUIRE(grid->get_global_ext()[1] == 4);
    REQUIRE(grid->get_num_objects() == 24);
    REQUIRE(grid->get_num_nonzero_objects() == 12);

    const int* mask = grid->get_land_mask();
    for (int i = 0; i < grid->get_num_objects(); i++) {
        if ((i / grid->get_local_ext()[0]) % 2 == 0) {
            REQUIRE(mask[i] == 0);
        }
    }

    // Cleanup
    delete grid;
}
