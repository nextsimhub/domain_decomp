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

  // Select partitioner
  Partitioner* partitioner = Partitioner::create(comm, argc, argv);

  // Partition grid
  partitioner->partition(*grid);

  // Store partitioning results in netCDF file
  int num_procs;
  MPI_Comm_size(comm, &num_procs);
  partitioner->save("partition_" + to_string(num_procs) + ".nc");

  // Cleanup
  delete grid;
  delete partitioner;

  // Finalize MPI
  MPI_Finalize();

  return 0;
}
