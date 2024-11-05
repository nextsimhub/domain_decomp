/*!
 * @file test_grid_2.cpp
 * @author Joe Wallwork <jw2423@cam.ac.uk>
 * @date 05 Nov 2024
 */

#include "Grid.hpp"

#include <doctest/extensions/doctest_mpi.h>

MPI_TEST_CASE("Grid: non-default dimension naming, 1 MPI rank", 1)
{
    // Build grid from netCDF file
    Grid* grid = Grid::create(test_comm, "./test_2.nc", "m", "n", { 1, 0 }, "land_mask");

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

MPI_TEST_CASE("Grid: non-default dimension naming, 2 MPI ranks", 2)
{
    // Build grid from netCDF file
    Grid* grid = Grid::create(test_comm, "./test_2.nc", "m", "n", { 1, 0 }, "land_mask");

    REQUIRE(grid->get_global_ext()[0] == 6);
    REQUIRE(grid->get_global_ext()[1] == 4);
    REQUIRE(grid->get_num_objects() == 12);
    REQUIRE(grid->get_num_nonzero_objects() == 6);

    const int* mask = grid->get_land_mask();
    for (int i = 0; i < grid->get_num_objects(); i++) {
        if ((i / grid->get_local_ext()[0]) % 2 == 0) {
            REQUIRE(mask[i] == 0);
        }
    }

    // Cleanup
    delete grid;
}
