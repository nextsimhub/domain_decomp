#include <cstdio>

#include <mpi.h>
#include "grid.h"

using namespace std;

int main(int argc, char* argv[])
{
  // Initialize MPI
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);

  // Build distributed grid from netCDF file
  Grid* grid = new Grid(comm, argc, argv, argv[1]);

  // Perform domain decomposition
  grid->partition();

  // Store partitioning results in netCDF file
  int num_procs;
  MPI_Comm_size(comm, &num_procs);
  grid->save("rcb_" + to_string(num_procs) + ".nc");

  // Cleanup
  delete grid;

  MPI_Finalize();

  return 0;
}
