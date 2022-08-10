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

Grid* Grid::create(MPI_Comm comm, const std::string& filename,
                   const std::string dim0_id, const std::string dim1_id,
                   const std::string mask_id)
{
  return new Grid(comm, filename, dim0_id, dim1_id, mask_id);
}

Grid::Grid(MPI_Comm comm, const std::string& filename,
           const std::string& dim0_id, const std::string& dim1_id,
           const std::string& mask_id)
    : _comm(comm)
{
  netCDF::NcFile nc_file(filename, netCDF::NcFile::read);
  // The data group contains the information we need
  netCDF::NcGroup data_group(nc_file.getGroup(data_id));

  MPI_Comm_rank(comm, &_rank);

  // Retrieve the extent of each dimension of interest. The dimensions of
  // interest are the spatial dimensions of the grid. These are named "x" and
  // "y" by default.
  _global_ext_0 = data_group.getDim(dim0_id).getSize();
  _global_ext_1 = data_group.getDim(dim1_id).getSize();

  // Initially we partition assuming there is no land mask
  // Figure out my subset of objects
  MPI_Comm_size(comm, &_num_procs);

  // Start from a 2D decomposition
  _num_procs_0 = _num_procs;
  _num_procs_1 = 1;
  find_factors(_num_procs, _num_procs_0, _num_procs_1);

  _local_ext_0 = ceil((float)_global_ext_0 / (_num_procs_0));
  _local_ext_1 = ceil((float)_global_ext_1 / (_num_procs_1));
  _global_0 = (_rank / _num_procs_1) * _local_ext_0;
  _global_1 = (_rank % _num_procs_1) * _local_ext_1;
  if ((_rank % _num_procs_1) == _num_procs_1 - 1) {
    _local_ext_1 = _global_ext_1 - (_rank % _num_procs_1) * _local_ext_1;
  }
  if ((_rank / _num_procs_1) == _num_procs_0 - 1) {
    _local_ext_0 = _global_ext_0 - (_rank / _num_procs_1) * _local_ext_0;
  }
  _num_objects = _local_ext_0 * _local_ext_1;

  // Retrieve the land mask, if available
  netCDF::NcVar mask_var = data_group.getVar(mask_id);
  if (!mask_var.isNull()) {
    // Verify the order of dimensions provided is correct by comparing to the
    // dimension order of the mask variable
    std::vector<netCDF::NcDim> dims = mask_var.getDims();
    if (dims[0].getName() != dim0_id || dims[1].getName() != dim1_id) {
      throw std::runtime_error("Dimension ordering provided does not match "
                               "ordering in netCDF grid file");
    }

    _land_mask.resize(_num_objects);
    std::vector<size_t> start(2);
    std::vector<size_t> count(2);
    std::vector<ptrdiff_t> stride(2);
    // Coordinate of first element
    start[0] = _global_0;
    start[1] = _global_1;
    // Number of elements in every extension
    count[0] = _local_ext_0;
    count[1] = _local_ext_1;
    // Stride in every extension
    stride[0] = 1;
    stride[1] = 1;
    mask_var.getVar(start, count, stride, _land_mask.data());

    // Apply land mask
    for (int i = 0; i < _num_objects; i++) {
      // The convention is that sea data points will have a positive value
      // and land points a zero value
      if (_land_mask[i] > 0) {
        int _local_0 = i / _local_ext_1;
        int _local_1 = i % _local_ext_1;
        int _global_0 = _local_0 + _global_0;
        int _global_1 = _local_1 + _global_1;
        _object_id.push_back(_global_0 * _global_ext_1 + _global_1);
        _sparse_to_dense.push_back(i);
        _num_nonzero_objects++;
      }
    }
  } else {
    _num_nonzero_objects = _num_objects;
    for (int i = 0; i < _num_objects; i++) {
      int _local_0 = i / _local_ext_1;
      int _local_1 = i % _local_ext_1;
      int _global_0 = _local_0 + _global_0;
      int _global_1 = _local_1 + _global_1;
      _object_id.push_back(_global_0 * _global_ext_1 + _global_1);
    }
  }

  nc_file.close();
}

int Grid::get_num_objects() const { return _num_objects; }
int Grid::get_num_nonzero_objects() const { return _num_nonzero_objects; }
int Grid::get_global_ext_0() const { return _global_ext_0; }
int Grid::get_global_ext_1() const { return _global_ext_1; }
int Grid::get_local_ext_0() const { return _local_ext_0; }
int Grid::get_local_ext_1() const { return _local_ext_1; }
void Grid::set_local_ext_0(int val) { _local_ext_0 = val; }
void Grid::set_local_ext_1(int val) { _local_ext_1 = val; }
int Grid::get_global_0() const { return _global_0; }
int Grid::get_global_1() const { return _global_1; }
void Grid::set_global_0(int val) { _global_0 = val; }
void Grid::set_global_1(int val) { _global_1 = val; }
int Grid::get_num_procs_0() const { return _num_procs_0; }
int Grid::get_num_procs_1() const { return _num_procs_1; }
const int* Grid::get_land_mask() const { return _land_mask.data(); }
const int* Grid::get_sparse_to_dense() const { return _sparse_to_dense.data(); }
const int* Grid::get_nonzero_object_ids() const { return _object_id.data(); }
