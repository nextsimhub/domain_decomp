/*!
 * @file main.cpp
 * @date 1 May 2022
 * @author Athena Elafrou <ae488@cam.ac.uk>
 */

#include <cstdio>

#include "Grid.hpp"
#include "Partitioner.hpp"

#include <mpi.h>

using namespace std;

int main(int argc, char* argv[])
{
  // Initialize MPI
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);

  // Build grid from netCDF file
  Grid* grid = Grid::create(comm, argc, argv, argv[1]);

  // Create a Zoltan partitioner
  Partitioner* partitioner
      = Partitioner::Factory::create(comm, argc, argv, PartitionerType::Zoltan);

  // Partition grid
  partitioner->partition(*grid);

  // Store partitioning results in netCDF file
  int num_procs;
  MPI_Comm_size(comm, &num_procs);
  partitioner->save_mask("partition_mask_" + to_string(num_procs) + ".nc");
  partitioner->save_metadata("partition_metadata_" + to_string(num_procs)
                             + ".nc");

  // Cleanup
  delete grid;
  delete partitioner;

  // Finalize MPI
  MPI_Finalize();

  return 0;
}
