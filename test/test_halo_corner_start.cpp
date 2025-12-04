/*!
 * @file halo_corner_start.cpp
 * @author Nirav Shah <nvs31@cam.ac.uk>
 * @date 20 Nov 2025
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

    // Non-periodic cases
    // Top right
    if (test_rank == 0){
        for (int p = 0; p < test_nb_procs; p++) {
        if (p != test_rank) {
            if (partitioner->is_corner_neighbour(domains[test_rank], domains[p], TOP_RIGHT)){
                start = partitioner->halo_corner_start(domains[test_rank], domains[p], TOP_RIGHT);
            }
        }
        }
        REQUIRE(start == 0);
    }

    // Bottom right
    if (test_rank == 1){
        for (int p = 0; p < test_nb_procs; p++) {
        if (p != test_rank) {
            if (partitioner->is_corner_neighbour(domains[test_rank], domains[p], BOTTOM_RIGHT)){
                start = partitioner->halo_corner_start(domains[test_rank], domains[p], BOTTOM_RIGHT);
            }
        }
        }
        REQUIRE(start == 6);
    }

    // Top left
    if (test_rank == 2){
        for (int p = 0; p < test_nb_procs; p++) {
        if (p != test_rank) {
            if (partitioner->is_corner_neighbour(domains[test_rank], domains[p], TOP_LEFT)){
                start = partitioner->halo_corner_start(domains[test_rank], domains[p], TOP_LEFT);
            }
        }
        }
        REQUIRE(start == 2);
    }

    // Bottom left
    if (test_rank == 3){
        for (int p = 0; p < test_nb_procs; p++) {
        if (p != test_rank) {
            if (partitioner->is_corner_neighbour(domains[test_rank], domains[p], BOTTOM_LEFT)){
                start = partitioner->halo_corner_start(domains[test_rank], domains[p], BOTTOM_LEFT);
            }
        }
        }
        REQUIRE(start == 8);
    }

    // Periodic cases
    // Top right
    if (test_rank == 2){
        for (int p = 0; p < test_nb_procs; p++) {
        if (p != test_rank) {
            if (partitioner->is_corner_neighbour(domains[test_rank], domains[p], TOP_RIGHT, true, true)){
                start = partitioner->halo_corner_start(domains[test_rank], domains[p], TOP_RIGHT, true, true);
            }
        }
        }
        REQUIRE(start == 0);
    }

    // Bottom right
    if (test_rank == 3){
        for (int p = 0; p < test_nb_procs; p++) {
        if (p != test_rank) {
            if (partitioner->is_corner_neighbour(domains[test_rank], domains[p], BOTTOM_RIGHT, true, true)){
                start = partitioner->halo_corner_start(domains[test_rank], domains[p], BOTTOM_RIGHT, true, true);
            }
        }
        }
        REQUIRE(start == 6);
    }

    // Top left
    if (test_rank == 0){
        for (int p = 0; p < test_nb_procs; p++) {
        if (p != test_rank) {
            if (partitioner->is_corner_neighbour(domains[test_rank], domains[p], TOP_LEFT, true, true)){
                start = partitioner->halo_corner_start(domains[test_rank], domains[p], TOP_LEFT, true, true);
            }
        }
        }
        REQUIRE(start == 2);
    }

    // Bottom left
    if (test_rank == 1){
        for (int p = 0; p < test_nb_procs; p++) {
        if (p != test_rank) {
            if (partitioner->is_corner_neighbour(domains[test_rank], domains[p], BOTTOM_LEFT, true, true)){
                start = partitioner->halo_corner_start(domains[test_rank], domains[p], BOTTOM_LEFT, true, true);
            }
        }
        }
        REQUIRE(start == 8);
    }

    // Cleanup
    delete grid;
    delete partitioner;

}
