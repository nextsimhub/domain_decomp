/*!
 * @file Partitioner.cpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @date 25 June 2022
 */

#include "Partitioner.hpp"
#include "ZoltanPartitioner.hpp"

#include <cmath>

#include <netcdf>
#include <netcdf_par.h>

Partitioner::Partitioner(MPI_Comm comm, int argc, char** argv)
{
  _comm = comm;
  MPI_Comm_size(comm, &_num_procs);
  MPI_Comm_rank(comm, &_rank);
}

void Partitioner::save_mask(const std::string& filename) const
{
  // Use C API for parallel I/O
  const int NDIMS = 2;
  int ncid, vid, dimids[NDIMS];
  size_t start[NDIMS], count[NDIMS];

  nc_create_par(filename.c_str(), NC_MPIIO | NC_NETCDF4, _comm, MPI_INFO_NULL,
                &ncid);
  nc_put_att_int(ncid, NC_GLOBAL, "num_processes", NC_SHORT, 1, &_num_procs);

  // Create 2 dimensions
  // The values to be written are associated with the netCDF variable by
  // assuming that the last dimension of the netCDF variable varies fastest in
  // the C interface
  nc_def_dim(ncid, "x", _global_ext_0, &dimids[0]);
  nc_def_dim(ncid, "y", _global_ext_1, &dimids[1]);

  // Create variables
  nc_def_var(ncid, "pid", NC_SHORT, NDIMS, dimids, &vid);

  // Write metadata to file
  nc_enddef(ncid);

  // Set up slab for this process
  start[0] = _global_0;
  start[1] = _global_1;
  count[0] = _local_ext_0;
  count[1] = _local_ext_1;

  // Store data
  nc_var_par_access(ncid, vid, NC_COLLECTIVE);
  nc_put_vara_int(ncid, vid, start, count, _proc_id.data());
  nc_close(ncid);
}

void Partitioner::save_metadata(const std::string& filename) const
{
  // Use C API for parallel I/O
  int ncid;
  nc_create_par(filename.c_str(), NC_MPIIO | NC_NETCDF4, _comm, MPI_INFO_NULL,
                &ncid);

  // Define dimensions
  int dimid;
  nc_def_dim(ncid, "P", _num_procs, &dimid);

  // Define variables
  int top_x_vid, top_y_vid;
  int cnt_x_vid, cnt_y_vid;
  nc_def_var(ncid, "global_x", NC_INT, 1, &dimid, &top_x_vid);
  nc_def_var(ncid, "global_y", NC_INT, 1, &dimid, &top_y_vid);
  nc_def_var(ncid, "local_extent_x", NC_INT, 1, &dimid, &cnt_x_vid);
  nc_def_var(ncid, "local_extent_y", NC_INT, 1, &dimid, &cnt_y_vid);

  // Write metadata to file
  nc_enddef(ncid);

  // Set up slab for this process
  const size_t start = _rank;

  // Store data
  nc_var_par_access(ncid, top_x_vid, NC_COLLECTIVE);
  nc_put_var1_int(ncid, top_x_vid, &start, &_global_0_cur);
  nc_var_par_access(ncid, top_y_vid, NC_COLLECTIVE);
  nc_put_var1_int(ncid, top_y_vid, &start, &_global_1_cur);
  nc_var_par_access(ncid, cnt_x_vid, NC_COLLECTIVE);
  nc_put_var1_int(ncid, cnt_x_vid, &start, &_local_ext_0_cur);
  nc_var_par_access(ncid, cnt_y_vid, NC_COLLECTIVE);
  nc_put_var1_int(ncid, cnt_y_vid, &start, &_local_ext_1_cur);

  nc_close(ncid);
}

Partitioner* Partitioner::Factory::create(MPI_Comm comm, int argc, char** argv,
                                          PartitionerType type)
{
  if (type == PartitionerType::Zoltan_RCB)
    return ZoltanPartitioner::create(comm, argc, argv);
  else
    throw std::runtime_error("Invalid partitioner!");
}
