/*!
 * @file Partitioner.cpp
 * @date 25 June 2022
 * @author Athena Elafrou <ae488@cam.ac.uk>
 */

#include "Partitioner.hpp"
#include "ZoltanPartitioner.hpp"

#include <netcdf>
#include <netcdf_par.h>

Partitioner::Partitioner(MPI_Comm comm, int argc, char** argv)
{
  _comm = comm;
  MPI_Comm_size(comm, &_num_procs);
  MPI_Comm_rank(comm, &_rank);
}

void Partitioner::save(const std::string& filename) const
{
  // Use C API for parallel I/O
  const int NDIMS = 2;
  int ncid, vid, dimids[NDIMS];
  size_t start[NDIMS], count[NDIMS];

  nc_create_par(filename.c_str(), NC_MPIIO | NC_NETCDF4, _comm, MPI_INFO_NULL,
                &ncid);
  nc_put_att_int(ncid, NC_GLOBAL, "num_processors", NC_INT, 1, &_num_procs);

  // Create 2 dimensions
  // The values to be written are associated with the netCDF variable by
  // assuming that the last dimension of the netCDF variable varies fastest in
  // the C interface
  nc_def_dim(ncid, "x", _global_dim_x, &dimids[0]);
  nc_def_dim(ncid, "y", _global_dim_y, &dimids[1]);
  // Create variables
  nc_def_var(ncid, "part_id", NC_INT, NDIMS, dimids, &vid);
  // Write metadata to file
  nc_enddef(ncid);
  // Set up slab for this process
  start[0] = _global_top_x;
  start[1] = _global_top_y;
  count[0] = _local_dim_x;
  count[1] = _local_dim_y;

  // Store data
  nc_var_par_access(ncid, vid, NC_COLLECTIVE);
  nc_put_vara_int(ncid, vid, start, count, _proc_id.data());
  nc_close(ncid);
}

Partitioner* Partitioner::Factory::create(MPI_Comm comm, int argc, char** argv,
                                          PartitionerType type)
{
  if (type == PartitionerType::Zoltan)
    return ZoltanPartitioner::create(comm, argc, argv);
  else
    throw std::runtime_error("Invalid partitioner!");
}
