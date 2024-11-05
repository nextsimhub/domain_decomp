/*!
 * @file test_zoltan_partitioner_1.cpp
 * @author Joe Wallwork <jw2423@cam.ac.uk>
 * @date 05 Nov 2024
 */

#include "Grid.hpp"
#include "Partitioner.hpp"

#include <doctest/extensions/doctest_mpi.h>

extern int global_argc;
extern char** global_argv;

MPI_TEST_CASE("ZoltanPartitioner: no land, 1 MPI rank", 1)
{
    // Build grid from netCDF file
    Grid* grid = Grid::create(test_comm, "./test_1.nc");

    // Create a Zoltan partitioner
    Partitioner* partitioner = Partitioner::Factory::create(
        test_comm, global_argc, global_argv, PartitionerType::Zoltan_RCB);

    // Partition grid
    partitioner->partition(*grid);

    int global_0, global_1, local_ext_0, local_ext_1;
    partitioner->get_bounding_box(global_0, global_1, local_ext_0, local_ext_1);
    REQUIRE(local_ext_0 == 6);
    REQUIRE(local_ext_1 == 4);
    REQUIRE(global_0 == 0);
    REQUIRE(global_1 == 0);

    // Cleanup
    delete grid;
    delete partitioner;
}

MPI_TEST_CASE("ZoltanPartitioner: no land, 2 MPI ranks", 2)
{
    // Build grid from netCDF file
    Grid* grid = Grid::create(test_comm, "./test_1.nc");

    // Create a Zoltan partitioner
    Partitioner* partitioner = Partitioner::Factory::create(
        test_comm, global_argc, global_argv, PartitionerType::Zoltan_RCB);

    // Partition grid
    partitioner->partition(*grid);

    int global_0, global_1, local_ext_0, local_ext_1;
    partitioner->get_bounding_box(global_0, global_1, local_ext_0, local_ext_1);
    REQUIRE(local_ext_0 == 3);
    REQUIRE(local_ext_1 == 4);
    if (test_rank == 0) {
        REQUIRE(global_0 == 0);
        REQUIRE(global_1 == 0);
    } else {
        REQUIRE(global_0 == 3);
        REQUIRE(global_1 == 0);
    }

    // Cleanup
    delete grid;
    delete partitioner;
}
