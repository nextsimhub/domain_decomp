/*!
 * @file test_haloCornerBufferPositions.cpp
 * @author Nirav Shah <nvs31@cam.ac.uk>
 * @date 20 Nov 2025
 */

#include "Grid.hpp"
#include "Partitioner.hpp"
#include "Utils.hpp"
#include <doctest/extensions/doctest_mpi.h>

extern int global_argc;
extern char** global_argv;

MPI_TEST_CASE("test corner neighbours and metadata for corner buffers", 4)
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
    std::vector<Point> origins(test_nb_procs);
    std::vector<Point> extents(test_nb_procs);
    std::vector<Domain> domains(test_nb_procs);
    std::vector<int> tmp0(test_nb_procs);
    std::vector<int> tmp1(test_nb_procs);

    CHECK_MPI(MPI_Allgather(&global_0, 1, MPI_INT, tmp0.data(), 1, MPI_INT, test_comm));
    CHECK_MPI(MPI_Allgather(&global_1, 1, MPI_INT, tmp1.data(), 1, MPI_INT, test_comm));

    // origin points mark the bottom-left corner of each domain
    for (int p = 0; p < test_nb_procs; p++) {
        origins[p].x = tmp0[p];
        origins[p].y = tmp1[p];
    }

    CHECK_MPI(MPI_Allgather(&local_ext_0, 1, MPI_INT, tmp0.data(), 1, MPI_INT, test_comm));
    CHECK_MPI(MPI_Allgather(&local_ext_1, 1, MPI_INT, tmp1.data(), 1, MPI_INT, test_comm));

    // extents can be used to find the top-right corner of each domain
    for (int p = 0; p < test_nb_procs; p++) {
        extents[p].x = tmp0[p];
        extents[p].y = tmp1[p];
    }

    // generate domains
    for (int p = 0; p < test_nb_procs; p++) {
        domains[p].p1.x = origins[p].x;
        domains[p].p1.y = origins[p].y;
        domains[p].p2.x = origins[p].x + extents[p].x;
        domains[p].p2.y = origins[p].y + extents[p].y;
    }

    int start;
    struct HaloCornerInfo {
        bool isPeriodic;
        int rank;
        int startPos;
    };

    const bool pxOn = true, pyOn = true;

    std::map<Corner, HaloCornerInfo> cornerInfoExpected, cornerInfoActual;

    if (test_rank == 0) {
        cornerInfoExpected[BOTTOM_LEFT] = { true, 3, 5 };
        cornerInfoExpected[BOTTOM_RIGHT] = { true, 3, 11 };
        cornerInfoExpected[TOP_RIGHT] = { false, 2, 15 };
        cornerInfoExpected[TOP_LEFT] = { true, 2, 6 };
    } else if (test_rank == 1) {
        cornerInfoExpected[BOTTOM_LEFT] = { true, 2, 5 };
        cornerInfoExpected[BOTTOM_RIGHT] = { false, 2, 14 };
        cornerInfoExpected[TOP_RIGHT] = { true, 2, 13 };
        cornerInfoExpected[TOP_LEFT] = { true, 2, 4 };
    } else if (test_rank == 2) {
        cornerInfoExpected[BOTTOM_LEFT] = { true, 1, 7 };
        cornerInfoExpected[BOTTOM_RIGHT] = { true, 1, 15 };
        cornerInfoExpected[TOP_RIGHT] = { true, 1, 14 };
        cornerInfoExpected[TOP_LEFT] = { false, 1, 6 };
    } else if (test_rank == 3) {
        cornerInfoExpected[BOTTOM_LEFT] = { false, 1, 5 };
        cornerInfoExpected[BOTTOM_RIGHT] = { true, 1, 13 };
        cornerInfoExpected[TOP_RIGHT] = { true, 0, 8 };
        cornerInfoExpected[TOP_LEFT] = { true, 0, 3 };
    }

    for (auto corner : corners) {
        for (int p = 0; p < test_nb_procs; p++) {
            // periodic neighbours
            if (partitioner->is_corner_neighbour(
                    domains[test_rank], domains[p], corner, pxOn, pyOn)) {
                partitioner->haloCornerBufferPositions(
                    domains[test_rank], domains[p], corner, start);
                cornerInfoActual[corner] = { true, p, start };
            }
            // non-periodic neighbours
            if (p != test_rank) {
                if (partitioner->is_corner_neighbour(domains[test_rank], domains[p], corner)) {
                    partitioner->haloCornerBufferPositions(
                        domains[test_rank], domains[p], corner, start);
                    cornerInfoActual[corner] = { false, p, start };
                }
            }
        }
    }

    for (auto corner : corners) {
        // check computed corner info (actual) against the expected information
        REQUIRE(cornerInfoActual[corner].isPeriodic == cornerInfoExpected[corner].isPeriodic);
        REQUIRE(cornerInfoActual[corner].rank == cornerInfoExpected[corner].rank);
        REQUIRE(cornerInfoActual[corner].startPos == cornerInfoExpected[corner].startPos);
    }

    // Cleanup
    delete grid;
    delete partitioner;
}
