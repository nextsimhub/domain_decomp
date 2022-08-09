/*!
 * @file Grid.cpp
 * @date 1 May 2022
 * @author Athena Elafrou <ae488@cam.ac.uk>
 */

#include "Grid.hpp"

#include <algorithm>
#include <cmath>

#include <netcdf>

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
  // The data group contains the information we need
  netCDF::NcGroup data_group(nc_file.getGroup(data_id));

  MPI_Comm_rank(comm, &_rank);

  // Retrieve the dimension of each dimension
  _global_dim_x = data_group.getDim(x_dim_id).getSize();
  _global_dim_y = data_group.getDim(y_dim_id).getSize();

  // Initiallly we partition assuming there is no land mask
  // Figure out my subset of objects
  MPI_Comm_size(comm, &_num_procs);

  // Start from a 2D decomposition
  _num_procs_x = _num_procs;
  _num_procs_y = 1;
  find_factors(_num_procs, _num_procs_x, _num_procs_y);

  _local_dim_x = ceil((float)_global_dim_x / (_num_procs_x));
  _local_dim_y = ceil((float)_global_dim_y / (_num_procs_y));
  _global_top_x = (_rank / _num_procs_y) * _local_dim_x;
  _global_top_y = (_rank % _num_procs_y) * _local_dim_y;
  if ((_rank % _num_procs_y) == _num_procs_y - 1) {
    _local_dim_y = _global_dim_y - (_rank % _num_procs_y) * _local_dim_y;
  }
  if ((_rank / _num_procs_y) == _num_procs_x - 1) {
    _local_dim_x = _global_dim_x - (_rank / _num_procs_y) * _local_dim_x;
  }
  _num_objects = _local_dim_x * _local_dim_y;

  std::vector<size_t> start(2);
  std::vector<size_t> count(2);
  std::vector<ptrdiff_t> stride(2);
  // Coordinate of first element
  start[0] = _global_top_x;
  start[1] = _global_top_y;
  // Number of elements in every dimension
  count[0] = _local_dim_x;
  count[1] = _local_dim_y;
  // Stride in every dimension
  stride[0] = 1;
  stride[1] = 1;

  // Retrieve the land mask, if available
  netCDF::NcVar mask_var = data_group.getVar(mask_id); // (x, y)
  if (!mask_var.isNull()) {
    _land_mask.resize(_num_objects);
    mask_var.getVar(start, count, stride, _land_mask.data());

    // Apply land mask
    for (int i = 0; i < _num_objects; i++) {
      // The convention is that sea data points will have a positive value
      // and land points a zero value
      if (_land_mask[i] > 0) {
        int _local_x = i / _local_dim_y;
        int _local_y = i % _local_dim_y;
        int _global_x = _local_x + _global_top_x;
        int _global_y = _local_y + _global_top_y;
        _object_id.push_back(_global_x * _global_dim_y + _global_y);
        _sparse_to_dense.push_back(i);
        _num_nonzero_objects++;
      }
    }
  } else {
    _num_nonzero_objects = _num_objects;
    for (int i = 0; i < _num_objects; i++) {
      int _local_x = i / _local_dim_y;
      int _local_y = i % _local_dim_y;
      int _global_x = _local_x + _global_top_x;
      int _global_y = _local_y + _global_top_y;
      _object_id.push_back(_global_x * _global_dim_y + _global_y);
    }
  }

  nc_file.close();
}

int Grid::get_num_objects() const { return _num_objects; }
int Grid::get_num_nonzero_objects() const { return _num_nonzero_objects; }
int Grid::get_global_dim_x() const { return _global_dim_x; }
int Grid::get_global_dim_y() const { return _global_dim_y; }
int Grid::get_local_dim_x() const { return _local_dim_x; }
int Grid::get_local_dim_y() const { return _local_dim_y; }
int Grid::get_global_top_x() const { return _global_top_x; }
int Grid::get_global_top_y() const { return _global_top_y; }
int Grid::get_num_procs_x() const { return _num_procs_x; }
int Grid::get_num_procs_y() const { return _num_procs_y; }
const int* Grid::get_land_mask() const { return _land_mask.data(); }
const int* Grid::get_sparse_to_dense() const { return _sparse_to_dense.data(); }
const int* Grid::get_nonzero_object_ids() const { return _object_id.data(); }
