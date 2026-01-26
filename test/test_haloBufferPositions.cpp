/*!
 * @file test_haloBufferPositions.cpp
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

MPI_TEST_CASE("Corner neighbour: Non-periodic, 4 MPI ranks", 4)
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
    int recv;

    struct HaloEdgeInfo {
        bool isPeriodic;
        int startPos;
        int recvPos;
    };

    const bool pxOn = true, pyOn = true;

    // The int in tuple<Edge, int> refers to rank. In the case of edge neighbours,
    // it is possible that one edge may corresspond to multiple neighbour unlike
    // corner neighbour where a given vertice can correspond to only one neighbour.
    std::map<std::tuple<Edge, int>, HaloEdgeInfo> edgeInfoExpected, edgeInfoActual;

    if (test_rank == 0) {
        edgeInfoExpected[{ BOTTOM, 1 }] = { true, 8, 0 };
        edgeInfoExpected[{ RIGHT, 2 }] = { false, 13, 3 };
        edgeInfoExpected[{ TOP, 1 }] = { false, 0, 5 };
        edgeInfoExpected[{ LEFT, 2 }] = { true, 4, 8 };
    } else if (test_rank == 1) {
        edgeInfoExpected[{ BOTTOM, 0 }] = { false, 5, 0 };
        edgeInfoExpected[{ RIGHT, 3 }] = { false, 10, 6 };
        edgeInfoExpected[{ RIGHT, 2 }] = { false, 15, 3 };
        edgeInfoExpected[{ TOP, 0 }] = { true, 0, 8 };
        edgeInfoExpected[{ LEFT, 3 }] = { true, 4, 14 };
        edgeInfoExpected[{ LEFT, 2 }] = { true, 6, 11 };
    } else if (test_rank == 2) {
        edgeInfoExpected[{ BOTTOM, 3 }] = { true, 6, 0 };
        edgeInfoExpected[{ RIGHT, 1 }] = { true, 11, 6 };
        edgeInfoExpected[{ RIGHT, 0 }] = { true, 8, 4 };
        edgeInfoExpected[{ TOP, 3 }] = { false, 0, 9 };
        edgeInfoExpected[{ LEFT, 1 }] = { false, 3, 15 };
        edgeInfoExpected[{ LEFT, 0 }] = { false, 3, 13 };
    } else if (test_rank == 3) {
        edgeInfoExpected[{ BOTTOM, 2 }] = { false, 9, 0 };
        edgeInfoExpected[{ RIGHT, 1 }] = { true, 14, 4 };
        edgeInfoExpected[{ TOP, 2 }] = { true, 0, 6 };
        edgeInfoExpected[{ LEFT, 1 }] = { false, 6, 10 };
    }

    for (auto edge : edges) {
        for (int p = 0; p < test_nb_procs; p++) {
            // periodic neighbours
            if (partitioner->is_neighbour(domains[test_rank], domains[p], edge, pxOn, pyOn)) {
                partitioner->haloBufferPositions(domains[test_rank], domains[p], edge, start, recv);
                edgeInfoActual[{ edge, p }] = { true, start, recv };
            }

            // non-periodic neighbours
            if (p != test_rank) {
                if (partitioner->is_neighbour(domains[test_rank], domains[p], edge)) {
                    partitioner->haloBufferPositions(
                        domains[test_rank], domains[p], edge, start, recv);
                    edgeInfoActual[{ edge, p }] = { false, start, recv };
                }
            }
        }
    }

    for (auto edge : edges) {
        for (int p = 0; p < test_nb_procs; p++) {
            if (p != test_rank && edgeInfoExpected.count({ edge, p }) != 0) {
                // check computed edge info (actual) against the expected information
                REQUIRE(edgeInfoActual[{ edge, p }].isPeriodic
                    == edgeInfoExpected[{ edge, p }].isPeriodic);
                REQUIRE(
                    edgeInfoActual[{ edge, p }].startPos == edgeInfoExpected[{ edge, p }].startPos);
                REQUIRE(
                    edgeInfoActual[{ edge, p }].recvPos == edgeInfoExpected[{ edge, p }].recvPos);
            }
        }
    }

    // Cleanup
    delete grid;
    delete partitioner;
}
