/*!
 * @file Grid.cpp
 * @date 1 May 2022
 * @author Athena Elafrou <ae488@cam.ac.uk>
 */

#include "Grid.hpp"

#include <algorithm>
#include <cmath>

#include <netcdf>
#include <netcdf_par.h>

static void find_factors(int n, int& factor_a, int& factor_b)
{
  for (int i = 2; i * i <= n; i += 2) {
    if (n % i == 0) {
      factor_a = i;
      factor_b = n / factor_a;
    }
  }
}

Grid* Grid::create(MPI_Comm comm, int argc, char** argv,
                   const std::string& filename)
{
  return new Grid(comm, argc, argv, filename);
}

Grid::Grid(MPI_Comm comm, int argc, char** argv, const std::string& filename)
    : _comm(comm)
{
  netCDF::NcFile nc_file(filename, netCDF::NcFile::read);

  MPI_Comm_rank(comm, &_rank);

  // Retrieve the dimension of each axis
  _global_dim_x = nc_file.getDim(x_axis_id).getSize();
  _global_dim_y = nc_file.getDim(y_axis_id).getSize();

  // Initiallly we partition assuming there is no land mask
  // Figure out my subset of objects
  MPI_Comm_size(comm, &_num_procs);

  // Start from a 2D decomposition
  int num_procs_y = 1;
  int num_procs_x = _num_procs;
  find_factors(_num_procs, num_procs_y, num_procs_x);

  _local_dim_x = ceil((float)_global_dim_x / (num_procs_x));
  _local_dim_y = ceil((float)_global_dim_y / (num_procs_y));
  _global_top_y = (_rank / num_procs_x) * _local_dim_y;
  _global_top_x = (_rank % num_procs_x) * _local_dim_x;
  if ((_rank % num_procs_x) == num_procs_x - 1) {
    _local_dim_x = _global_dim_x - (_rank % num_procs_x) * _local_dim_x;
  }
  if ((_rank / num_procs_x) == num_procs_y - 1) {
    _local_dim_y = _global_dim_y - (_rank / num_procs_x) * _local_dim_y;
  }
  _num_objects = _local_dim_x * _local_dim_y;

  // Retrieve the land mask
  netCDF::NcVar mask_var = nc_file.getVar(mask_id); // (y, x)
  std::vector<size_t> start(2);
  std::vector<size_t> count(2);
  std::vector<ptrdiff_t> stride(2);
  // Coordinate of first element
  start[0] = _global_top_y; // y
  start[1] = _global_top_x; // x
  // Number of elements in every dimension
  count[0] = _local_dim_y; // y
  count[1] = _local_dim_x; // x
  // Stride in every dimension
  stride[0] = 1; // y
  stride[1] = 1; // x
  _land_mask.resize(_num_objects);
  mask_var.getVar(start, count, stride, _land_mask.data());

  // Apply land mask
  for (int i = 0; i < _num_objects; i++) {
    // The convention is that sea data points will have a positive value
    // and land points a negative value
    if (_land_mask[i] > 0) {
      _num_nonzero_objects++;
      int _local_y = i / _local_dim_x;
      int _local_x = i % _local_dim_x;
      int _global_y = _local_y + _global_top_y;
      int _global_x = _local_x + _global_top_x;
      _index_map.push_back(_global_y * _global_dim_x + _global_x);
    }
  }

  // Set unique object IDs
  _object_ids.resize(_num_nonzero_objects);
  for (int i = 0; i < _num_nonzero_objects; i++) {
    _object_ids[i] = _index_map[i];
  }

  nc_file.close();
}

int Grid::get_num_nonzero_objects() const { return _num_nonzero_objects; }
const int* Grid::get_nonzero_object_ids() const { return _object_ids.data(); }
