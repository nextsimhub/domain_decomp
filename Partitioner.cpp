/*!
 * @file Partitioner.cpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @date 25 June 2022
 */

#include "Partitioner.hpp"
#include "ZoltanPartitioner.hpp"

#include <algorithm>
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
  nc_def_dim(ncid, "x", _global_dim_x, &dimids[0]);
  nc_def_dim(ncid, "y", _global_dim_y, &dimids[1]);

  // Create variables
  nc_def_var(ncid, "pid", NC_SHORT, NDIMS, dimids, &vid);

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

void Partitioner::save_metadata(const std::string& filename) const
{
  // Find local bounding boxes for each process
  std::vector<int> min_x(_num_procs, _global_dim_x);
  std::vector<int> min_y(_num_procs, _global_dim_y);
  std::vector<int> max_x(_num_procs, -1);
  std::vector<int> max_y(_num_procs, -1);
  for (int i = 0; i < _proc_id.size(); i++) {
    // Find global 2D coordinates of element
    int x = i / _local_dim_y + (_rank / _num_procs_y) * _local_dim_x;
    int y = i % _local_dim_y + (_rank % _num_procs_y) * _local_dim_y;
    if (_proc_id[i] != -1) {
      if (x > max_x[_proc_id[i]])
        max_x[_proc_id[i]] = x;
      if (y > max_y[_proc_id[i]])
        max_y[_proc_id[i]] = y;
      if (x < min_x[_proc_id[i]])
        min_x[_proc_id[i]] = x;
      if (y < min_y[_proc_id[i]])
        min_y[_proc_id[i]] = y;
    }
  }

  // Find global bounding boxes for each process
  std::vector<int> global_min_x(_num_procs);
  std::vector<int> global_min_y(_num_procs);
  std::vector<int> global_max_x(_num_procs);
  std::vector<int> global_max_y(_num_procs);
  MPI_Allreduce(min_x.data(), global_min_x.data(), _num_procs, MPI_INT, MPI_MIN,
                _comm);
  MPI_Allreduce(min_y.data(), global_min_y.data(), _num_procs, MPI_INT, MPI_MIN,
                _comm);
  MPI_Allreduce(max_x.data(), global_max_x.data(), _num_procs, MPI_INT, MPI_MAX,
                _comm);
  MPI_Allreduce(max_y.data(), global_max_y.data(), _num_procs, MPI_INT, MPI_MAX,
                _comm);

  // Expand bounding boxes to account for land points in x dimension
  std::vector<int> coords(_num_procs_x * 2);
  for (int j = 0; j < _num_procs_y; j++) {
    const int step
        = (_num_procs_x == _num_procs || _num_procs_y == _num_procs) ? 1 : 2;
    // Sort points on x dimension for processes in column j
    int cnt = 0;
    for (int p = j; p < _num_procs; p += step) {
      coords[2 * cnt] = global_min_x[p];
      coords[2 * cnt + 1] = global_max_x[p];
      cnt++;
    }
    std::sort(coords.begin(), coords.end());
    for (int i = 1; i < coords.size() - 1; i += 2) {
      if (coords[i + 1] - coords[i] != 1)
        coords[i] += coords[i + 1] - coords[i] - 1;
    }
    // Correct first and last point if necessary
    coords[0] = 0;
    coords[2 * _num_procs_x - 1] = _global_dim_x - 1;

    // Update bounding boxes
    cnt = 0;
    for (int p = j; p < _num_procs; p += step) {
      global_min_x[p] = coords[2 * cnt];
      global_max_x[p] = coords[2 * cnt + 1];
      cnt++;
    }
  }

  // Expand bounding boxes to account for land points in y dimension
  coords.clear();
  coords.resize(_num_procs_y * 2);
  for (int j = 0; j < _num_procs_x; j++) {
    const int step
        = (_num_procs_x == _num_procs || _num_procs_y == _num_procs) ? 1 : 2;
    // Sort points on y dimension for processes in row j
    int cnt = 0;
    for (int p = step * j; p < (j + 1) * _num_procs_y; p++) {
      coords[2 * cnt] = global_min_y[p];
      coords[2 * cnt + 1] = global_max_y[p];
      cnt++;
    }
    std::sort(coords.begin(), coords.end());
    for (int i = 1; i < coords.size() - 1; i += 2) {
      if (coords[i + 1] - coords[i] != 1)
        coords[i] += coords[i + 1] - coords[i] - 1;
    }
    // Correct first and last point if necessary
    coords[0] = 0;
    coords[2 * _num_procs_y - 1] = _global_dim_y - 1;

    cnt = 0;
    for (int p = step * j; p < (j + 1) * _num_procs_y; p++) {
      global_min_y[p] = coords[2 * cnt];
      global_max_y[p] = coords[2 * cnt + 1];
      cnt++;
    }
  }

  int global_top_x = global_min_x[_rank];
  int global_top_y = global_min_y[_rank];
  int local_dim_x = global_max_x[_rank] - global_top_x + 1;
  int local_dim_y = global_max_y[_rank] - global_top_y + 1;

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
  nc_put_var1_int(ncid, top_x_vid, &start, &global_top_x);
  nc_var_par_access(ncid, top_y_vid, NC_COLLECTIVE);
  nc_put_var1_int(ncid, top_y_vid, &start, &global_top_y);
  nc_var_par_access(ncid, cnt_x_vid, NC_COLLECTIVE);
  nc_put_var1_int(ncid, cnt_x_vid, &start, &local_dim_x);
  nc_var_par_access(ncid, cnt_y_vid, NC_COLLECTIVE);
  nc_put_var1_int(ncid, cnt_y_vid, &start, &local_dim_y);

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
