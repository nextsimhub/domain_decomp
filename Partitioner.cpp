/*!
 * @file Partitioner.cpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @date 25 June 2022
 */

#include "Partitioner.hpp"
#include "ZoltanPartitioner.hpp"

#include <cmath>
#include <stdexcept>

#include <netcdf.h>
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
  int nc_id, nc_mode;
  nc_mode = NC_MPIIO | NC_NETCDF4;
  NC_CHECK(
      nc_create_par(filename.c_str(), nc_mode, _comm, MPI_INFO_NULL, &nc_id));
  NC_CHECK(nc_put_att_int(nc_id, NC_GLOBAL, "num_processes", NC_SHORT, 1,
                          &_num_procs));

  // Create 2 dimensions
  // The values to be written are associated with the netCDF variable by
  // assuming that the last dimension of the netCDF variable varies fastest in
  // the C interface
  const int NDIMS = 2;
  int dimid[NDIMS];
  NC_CHECK(nc_def_dim(nc_id, "x", _global_ext_0, &dimid[0]));
  NC_CHECK(nc_def_dim(nc_id, "y", _global_ext_1, &dimid[1]));

  // Create variables
  int mask_nc_id;
  NC_CHECK(nc_def_var(nc_id, "pid", NC_SHORT, NDIMS, dimid, &mask_nc_id));

  // Write metadata to file
  NC_CHECK(nc_enddef(nc_id));

  // Set up slab for this process
  size_t start[NDIMS], count[NDIMS];
  start[0] = _global_0;
  start[1] = _global_1;
  count[0] = _local_ext_0;
  count[1] = _local_ext_1;
  // Store data
  NC_CHECK(nc_var_par_access(nc_id, mask_nc_id, NC_COLLECTIVE));
  NC_CHECK(nc_put_vara_int(nc_id, mask_nc_id, start, count, _proc_id.data()));
  NC_CHECK(nc_close(nc_id));
}

void Partitioner::save_metadata(const std::string& filename) const
{
  // Use C API for parallel I/O
  int nc_id, nc_mode;
  nc_mode = NC_MPIIO | NC_NETCDF4;
  NC_CHECK(
      nc_create_par(filename.c_str(), nc_mode, _comm, MPI_INFO_NULL, &nc_id));

  // Define dimensions
  int dimid;
  NC_CHECK(nc_def_dim(nc_id, "P", _num_procs, &dimid));

  // Define variables
  int top_x_vid, top_y_vid;
  int cnt_x_vid, cnt_y_vid;
  NC_CHECK(nc_def_var(nc_id, "global_x", NC_INT, 1, &dimid, &top_x_vid));
  NC_CHECK(nc_def_var(nc_id, "global_y", NC_INT, 1, &dimid, &top_y_vid));
  NC_CHECK(nc_def_var(nc_id, "local_extent_x", NC_INT, 1, &dimid, &cnt_x_vid));
  NC_CHECK(nc_def_var(nc_id, "local_extent_y", NC_INT, 1, &dimid, &cnt_y_vid));

  // Write metadata to file
  NC_CHECK(nc_enddef(nc_id));

  // Set up slab for this process
  const size_t start = _rank;

  // Store data
  NC_CHECK(nc_var_par_access(nc_id, top_x_vid, NC_COLLECTIVE));
  NC_CHECK(nc_put_var1_int(nc_id, top_x_vid, &start, &_global_0_new));
  NC_CHECK(nc_var_par_access(nc_id, top_y_vid, NC_COLLECTIVE));
  NC_CHECK(nc_put_var1_int(nc_id, top_y_vid, &start, &_global_1_new));
  NC_CHECK(nc_var_par_access(nc_id, cnt_x_vid, NC_COLLECTIVE));
  NC_CHECK(nc_put_var1_int(nc_id, cnt_x_vid, &start, &_local_ext_0_new));
  NC_CHECK(nc_var_par_access(nc_id, cnt_y_vid, NC_COLLECTIVE));
  NC_CHECK(nc_put_var1_int(nc_id, cnt_y_vid, &start, &_local_ext_1_new));

  NC_CHECK(nc_close(nc_id));
}

Partitioner* Partitioner::Factory::create(MPI_Comm comm, int argc, char** argv,
                                          PartitionerType type)
{
  if (type == PartitionerType::Zoltan_RCB)
    return ZoltanPartitioner::create(comm, argc, argv);
  else
    throw std::runtime_error("Invalid partitioner!");
}
