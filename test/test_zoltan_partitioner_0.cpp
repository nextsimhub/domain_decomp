/*!
 * @file test_zoltan_partitioner_0.cpp
 * @author Joe Wallwork <jw2423@cam.ac.uk>
 * @date 05 Nov 2024
 */

#include "Grid.hpp"
#include "Partitioner.hpp"

#include <doctest/extensions/doctest_mpi.h>

extern int global_argc;
extern char** global_argv;

MPI_TEST_CASE("ZoltanPartitioner: all land, 1 MPI rank", 1)
{
    // Build grid from netCDF file
    Grid* grid = Grid::create(test_comm, "./test_0.nc");

    // Create a Zoltan partitioner
    Partitioner* partitioner = Partitioner::Factory::create(
        test_comm, global_argc, global_argv, PartitionerType::Zoltan_RCB);

    // Partition grid
    partitioner->partition(*grid);

    int global0, global1, localExt0, localExt1;
    partitioner->getBoundingBox(global0, global1, localExt0, localExt1);
    REQUIRE(localExt0 == 6);
    REQUIRE(localExt1 == 4);
    REQUIRE(global0 == 0);
    REQUIRE(global1 == 0);

    // Cleanup
    delete grid;
    delete partitioner;
}

MPI_TEST_CASE("ZoltanPartitioner: all land, 2 MPI ranks", 2)
{
    // Build grid from netCDF file
    Grid* grid = Grid::create(test_comm, "./test_0.nc");

    // Create a Zoltan partitioner
    Partitioner* partitioner = Partitioner::Factory::create(
        test_comm, global_argc, global_argv, PartitionerType::Zoltan_RCB);

    // Partition grid
    partitioner->partition(*grid);

    int global0, global1, localExt0, localExt1;
    partitioner->getBoundingBox(global0, global1, localExt0, localExt1);
    REQUIRE(localExt0 == 3);
    REQUIRE(localExt1 == 4);
    if (test_rank == 0) {
        REQUIRE(global0 == 0);
        REQUIRE(global1 == 0);
    } else {
        REQUIRE(global0 == 3);
        REQUIRE(global1 == 0);
    }

    // Cleanup
    delete grid;
    delete partitioner;
}
