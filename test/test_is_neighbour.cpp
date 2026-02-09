/*!
 * @file test_is_neighbour.cpp
 * @author Nirav Shah <nvs31@cam.ac.uk>
 * @date 19 January 2026
 */

#include "Grid.hpp"
#include "Partitioner.hpp"
#include "Utils.hpp"
#include <doctest/extensions/doctest_mpi.h>

#include <iostream>

extern int global_argc;
extern char** global_argv;

MPI_TEST_CASE("Neighbour, 4 MPI ranks", 4)
{
    // Build grid from netCDF file
    Grid* grid = Grid::create(test_comm, "./test_3.nc");

    // Create a Zoltan partitioner
    Partitioner* partitioner = Partitioner::Factory::create(
        test_comm, global_argc, global_argv, PartitionerType::Zoltan_RCB);

    // Partition grid
    partitioner->partition(*grid);

    int global_0, global_1, local_ext_0, local_ext_1;
    partitioner->get_bounding_box(global_0, global_1, local_ext_0, local_ext_1);

    // Gather bounding boxes for all processes
    std::vector<Point> origins(4);
    std::vector<Point> extents(4);
    std::vector<Domain> domains(4);
    std::vector<int> tmp0(4);
    std::vector<int> tmp1(4);

    CHECK_MPI(MPI_Allgather(&global_0, 1, MPI_INT, tmp0.data(), 1, MPI_INT, test_comm));
    CHECK_MPI(MPI_Allgather(&global_1, 1, MPI_INT, tmp1.data(), 1, MPI_INT, test_comm));

    // origin points mark the bottom-left corner of each domain
    for (int p = 0; p < 4; p++) {
        origins[p].x = tmp0[p];
        origins[p].y = tmp1[p];
    }

    CHECK_MPI(MPI_Allgather(&local_ext_0, 1, MPI_INT, tmp0.data(), 1, MPI_INT, test_comm));
    CHECK_MPI(MPI_Allgather(&local_ext_1, 1, MPI_INT, tmp1.data(), 1, MPI_INT, test_comm));

    // extents can be used to find the top-right corner of each domain
    for (int p = 0; p < 4; p++) {
        extents[p].x = tmp0[p];
        extents[p].y = tmp1[p];
    }

    // generate domains
    for (int p = 0; p < 4; p++) {
        domains[p].p1.x = origins[p].x;
        domains[p].p1.y = origins[p].y;
        domains[p].p2.x = origins[p].x + extents[p].x;
        domains[p].p2.y = origins[p].y + extents[p].y;
    }

    // Non-periodic cases
    bool neighbour;

    // RIGHT
    if (test_rank == 0) {
        neighbour = partitioner->is_neighbour(domains[0], domains[2], RIGHT);
        REQUIRE(neighbour);
    }

    // BOTTOM
    if (test_rank == 1) {
        neighbour = partitioner->is_neighbour(domains[1], domains[0], BOTTOM);
        REQUIRE(neighbour);
    }

    // TOP
    if (test_rank == 2) {
        neighbour = partitioner->is_neighbour(domains[2], domains[3], TOP);
        REQUIRE(neighbour);
    }

    // LEFT
    if (test_rank == 3) {
        neighbour = partitioner->is_neighbour(domains[3], domains[1], LEFT);
        REQUIRE(neighbour);
    }

    // Periodic cases
    bool periodic_neighbour;

    // LEFT
    if (test_rank == 0) {
        periodic_neighbour = partitioner->is_neighbour(domains[0], domains[2], LEFT, true, true);
        REQUIRE(periodic_neighbour);
    }

    // TOP
    if (test_rank == 1) {
        periodic_neighbour = partitioner->is_neighbour(domains[1], domains[0], TOP, true, true);
        REQUIRE(periodic_neighbour);
    }

    // BOTTOM
    if (test_rank == 2) {
        periodic_neighbour = partitioner->is_neighbour(domains[2], domains[3], BOTTOM, true, true);
        REQUIRE(periodic_neighbour);
    }

    // RIGHT
    if (test_rank == 3) {
        periodic_neighbour = partitioner->is_neighbour(domains[3], domains[1], RIGHT, true, true);
        REQUIRE(periodic_neighbour);
    }

    // Cleanup
    delete grid;
    delete partitioner;
}

MPI_TEST_CASE("Neighbour, 4 MPI ranks", 4)
{
    // Build grid from netCDF file
    Grid* grid = Grid::create(
        test_comm, "./test_4.nc", "x", "y", { 1, 0 }, "land_mask", false, true, true);

    // Create a Zoltan partitioner
    Partitioner* partitioner = Partitioner::Factory::create(
        test_comm, global_argc, global_argv, PartitionerType::Zoltan_RCB);

    // Partition grid
    partitioner->partition(*grid);

    int global_0, global_1, local_ext_0, local_ext_1;
    partitioner->get_bounding_box(global_0, global_1, local_ext_0, local_ext_1);

    // Gather bounding boxes for all processes
    std::vector<Point> origins(4);
    std::vector<Point> extents(4);
    std::vector<Domain> domains(4);
    std::vector<int> tmp0(4);
    std::vector<int> tmp1(4);

    CHECK_MPI(MPI_Allgather(&global_0, 1, MPI_INT, tmp0.data(), 1, MPI_INT, test_comm));
    CHECK_MPI(MPI_Allgather(&global_1, 1, MPI_INT, tmp1.data(), 1, MPI_INT, test_comm));

    // origin points mark the bottom-left corner of each domain
    for (int p = 0; p < 4; p++) {
        origins[p].x = tmp0[p];
        origins[p].y = tmp1[p];
    }

    CHECK_MPI(MPI_Allgather(&local_ext_0, 1, MPI_INT, tmp0.data(), 1, MPI_INT, test_comm));
    CHECK_MPI(MPI_Allgather(&local_ext_1, 1, MPI_INT, tmp1.data(), 1, MPI_INT, test_comm));

    // extents can be used to find the top-right corner of each domain
    for (int p = 0; p < 4; p++) {
        extents[p].x = tmp0[p];
        extents[p].y = tmp1[p];
    }

    // generate domains
    for (int p = 0; p < 4; p++) {
        domains[p].p1.x = origins[p].x;
        domains[p].p1.y = origins[p].y;
        domains[p].p2.x = origins[p].x + extents[p].x;
        domains[p].p2.y = origins[p].y + extents[p].y;
    }

    // Non-periodic cases
    bool neighbour;

    // RIGHT
    if (test_rank == 0) {
        neighbour = partitioner->is_neighbour(domains[0], domains[2], RIGHT);
        REQUIRE(neighbour);
    }

    // BOTTOM
    if (test_rank == 1) {
        neighbour = partitioner->is_neighbour(domains[1], domains[0], BOTTOM);
        REQUIRE(neighbour);
    }

    // TOP
    if (test_rank == 2) {
        neighbour = partitioner->is_neighbour(domains[2], domains[3], TOP);
        REQUIRE(neighbour);
    }

    // LEFT
    if (test_rank == 3) {
        neighbour = partitioner->is_neighbour(domains[3], domains[0], LEFT);
        REQUIRE(neighbour == false);
    }

    // Periodic cases
    bool periodic_neighbour;

    // LEFT
    if (test_rank == 0) {
        periodic_neighbour = partitioner->is_neighbour(domains[0], domains[2], LEFT, true, true);
        REQUIRE(periodic_neighbour);
    }

    // TOP
    if (test_rank == 1) {
        periodic_neighbour = partitioner->is_neighbour(domains[1], domains[0], TOP, true, true);
        REQUIRE(periodic_neighbour);
    }

    // BOTTOM
    if (test_rank == 2) {
        periodic_neighbour = partitioner->is_neighbour(domains[2], domains[3], BOTTOM, true, true);
        REQUIRE(periodic_neighbour);
    }

    // RIGHT
    if (test_rank == 3) {
        periodic_neighbour = partitioner->is_neighbour(domains[3], domains[0], RIGHT, true, true);
        REQUIRE(periodic_neighbour == false);
    }

    // Cleanup
    delete grid;
    delete partitioner;
}
