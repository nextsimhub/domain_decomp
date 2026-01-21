/*!
 * @file halo_corner_start.cpp
 * @author Nirav Shah <nvs31@cam.ac.uk>
 * @date 25 Nov 2025
 */

#include "Grid.hpp"
#include "Partitioner.hpp"
#include "Utils.hpp"
#include <doctest/extensions/doctest_mpi.h>

extern int global_argc;
extern char** global_argv;

MPI_TEST_CASE("Corner neighbour: Non-periodic, 4 MPI ranks", 4)
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
    bool corner_neighbour;

    if (test_rank == 0) {
        corner_neighbour
            = partitioner->is_corner_neighbour(domains[test_rank], domains[3], TOP_RIGHT);
        REQUIRE(corner_neighbour);
    }

    if (test_rank == 1) {
        corner_neighbour
            = partitioner->is_corner_neighbour(domains[test_rank], domains[2], BOTTOM_RIGHT);
        REQUIRE(corner_neighbour);
    }

    if (test_rank == 2) {
        corner_neighbour
            = partitioner->is_corner_neighbour(domains[test_rank], domains[1], TOP_LEFT);
        REQUIRE(corner_neighbour);
    }

    if (test_rank == 3) {
        corner_neighbour
            = partitioner->is_corner_neighbour(domains[test_rank], domains[0], BOTTOM_LEFT);
        REQUIRE(corner_neighbour);
    }

    // Periodic cases
    bool corner_neighbour_periodic;

    if (test_rank == 0) {
        corner_neighbour_periodic = partitioner->is_corner_neighbour(
            domains[test_rank], domains[3], TOP_LEFT, true, true);
        REQUIRE(corner_neighbour_periodic);
    }

    if (test_rank == 1) {
        corner_neighbour_periodic = partitioner->is_corner_neighbour(
            domains[test_rank], domains[2], BOTTOM_LEFT, true, true);
        REQUIRE(corner_neighbour_periodic);
    }

    if (test_rank == 2) {
        corner_neighbour_periodic = partitioner->is_corner_neighbour(
            domains[test_rank], domains[1], TOP_RIGHT, true, true);
        REQUIRE(corner_neighbour_periodic);
    }

    if (test_rank == 3) {
        corner_neighbour_periodic = partitioner->is_corner_neighbour(
            domains[test_rank], domains[0], BOTTOM_RIGHT, true, true);
        REQUIRE(corner_neighbour_periodic);
    }

    // Cleanup
    delete grid;
    delete partitioner;
}
