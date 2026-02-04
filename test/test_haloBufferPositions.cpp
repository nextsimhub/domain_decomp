/*!
 * @file test_haloBufferPositions.cpp
 * @author Nirav Shah <nvs31@cam.ac.uk>
 * @date 19 January 2026
 */

#include "Grid.hpp"
#include "Partitioner.hpp"
#include "Utils.hpp"
#include <doctest/extensions/doctest_mpi.h>
#include <list>

#include <iostream>

extern int global_argc;
extern char** global_argv;

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
        int rank;
        int startPos;
        int recvPos;
    };

    const bool pxOn = true, pyOn = true;

    // The int in tuple<Edge, int> refers to rank. In the case of edge neighbours,
    // it is possible that one edge may corresspond to multiple neighbour unlike
    // corner neighbour where a given vertice can correspond to only one neighbour.
    std::map<Edge, std::list<HaloEdgeInfo>> edgeInfoExpected, edgeInfoActual;

    if (test_rank == 0) {
        edgeInfoExpected[LEFT].push_back({ true, 2, 4, 8 });
        edgeInfoExpected[RIGHT].push_back({ false, 2, 13, 3 });
        edgeInfoExpected[BOTTOM].push_back({ true, 1, 8, 0 });
        edgeInfoExpected[TOP].push_back({ false, 1, 0, 5 });
    } else if (test_rank == 1) {
        edgeInfoExpected[LEFT].push_back({ true, 2, 6, 11 });
        edgeInfoExpected[LEFT].push_back({ true, 3, 4, 14 });
        edgeInfoExpected[RIGHT].push_back({ false, 2, 15, 3 });
        edgeInfoExpected[RIGHT].push_back({ false, 3, 10, 6 });
        edgeInfoExpected[BOTTOM].push_back({ false, 0, 5, 0 });
        edgeInfoExpected[TOP].push_back({ true, 0, 0, 8 });
    } else if (test_rank == 2) {
        edgeInfoExpected[LEFT].push_back({ false, 0, 3, 13 });
        edgeInfoExpected[LEFT].push_back({ false, 1, 3, 15 });
        edgeInfoExpected[RIGHT].push_back({ true, 0, 8, 4 });
        edgeInfoExpected[RIGHT].push_back({ true, 1, 11, 6 });
        edgeInfoExpected[BOTTOM].push_back({ true, 3, 6, 0 });
        edgeInfoExpected[TOP].push_back({ false, 3, 0, 9 });
    } else if (test_rank == 3) {
        edgeInfoExpected[LEFT].push_back({ false, 1, 6, 10 });
        edgeInfoExpected[RIGHT].push_back({ true, 1, 14, 4 });
        edgeInfoExpected[BOTTOM].push_back({ false, 2, 9, 0 });
        edgeInfoExpected[TOP].push_back({ true, 2, 0, 6 });
    }

    bool periodic_check_neighbour; 
    periodic_check_neighbour = (partitioner->is_neighbour(domains[3], domains[0], RIGHT, pxOn, pyOn));

    bool check_neighbour;
    check_neighbour = (partitioner->is_neighbour(domains[3], domains[0], LEFT));

    for (auto edge : edges) {
        for (int p = 0; p < test_nb_procs; p++) {
            // periodic neighbours
            if (p != test_rank) {
                if (partitioner->is_neighbour(domains[test_rank], domains[p], edge, pxOn, pyOn)) {
                    std::cout << "Periodic" << test_rank << p << edge << std::endl;
                    partitioner->haloBufferPositions(
                        domains[test_rank], domains[p], edge, start, recv);
                    edgeInfoActual[edge].push_back({ true, p, start, recv });
                }
            }

            // non-periodic neighbours
            if (p != test_rank) {
                if (partitioner->is_neighbour(domains[test_rank], domains[p], edge)) {
                    std::cout << "Non-periodic" << test_rank << p << edge << std::endl;
                    partitioner->haloBufferPositions(
                        domains[test_rank], domains[p], edge, start, recv);
                    edgeInfoActual[edge].push_back({ false, p, start, recv });
                }
            }
        }
    }

    for (auto edge : edges) {
        std::cout << test_rank << edge << edgeInfoActual[edge].size()
                  << edgeInfoExpected[edge].size() << std::endl;
        REQUIRE(edgeInfoActual[edge].size() == edgeInfoExpected[edge].size());
        if (edgeInfoActual[edge].empty()) {
            REQUIRE(edgeInfoExpected[edge].empty());
        } else {
            for (auto edgeInfoA : edgeInfoActual[edge]) {
                for (auto edgeInfoE : edgeInfoExpected[edge]) {
                    if (edgeInfoA.rank == edgeInfoE.rank) {
                        REQUIRE(edgeInfoA.isPeriodic == edgeInfoE.isPeriodic);
                        REQUIRE(edgeInfoA.startPos == edgeInfoE.startPos);
                        REQUIRE(edgeInfoA.recvPos == edgeInfoE.recvPos);
                    }
                }
            }
        }
    }

    // Cleanup
    delete grid;
    delete partitioner;
}
