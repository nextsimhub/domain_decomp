/*!
 * @file Grid.cpp
 * @date 23 Aug 2024
 * @author Athena Elafrou <ae488@cam.ac.uk>
 */

#include "Grid.hpp"
#include "Utils.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <netcdf.h>
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

Grid* Grid::create(MPI_Comm comm, const std::string& filename, bool ignore_mask)
{
    return new Grid(comm, filename, "x", "y", "mask", ignore_mask);
}

Grid* Grid::create(MPI_Comm comm, const std::string& filename, const std::string dim0_name,
    const std::string dim1_name, const std::string mask_name, bool ignore_mask)
{
    return new Grid(comm, filename, dim0_name, dim1_name, mask_name, ignore_mask);
}

Grid::Grid(MPI_Comm comm, const std::string& filename, const std::string& dim0_name,
    const std::string& dim1_name, const std::string& mask_name, bool ignore_mask)
    : _comm(comm)
{
    // Use C API for parallel I/O
    int nc_id, nc_o_mode;
    nc_o_mode = NC_NOWRITE;
    NC_CHECK(nc_open_par(filename.c_str(), nc_o_mode, _comm, MPI_INFO_NULL, &nc_id));

    // Extract group ID in case of enhanced data model
    int data_nc_id;
    int ret = nc_inq_ncid(nc_id, "data", &data_nc_id);
    if (ret != NC_NOERR)
        data_nc_id = nc_id;
    int dim0_nc_id, dim1_nc_id;
    // I have switched the order of the dimensions here to reflect the change of
    // indices in nextsim-dg
    NC_CHECK(nc_inq_dimid(data_nc_id, dim1_name.c_str(), &dim0_nc_id));
    NC_CHECK(nc_inq_dimid(data_nc_id, dim0_name.c_str(), &dim1_nc_id));

    CHECK_MPI(MPI_Comm_rank(comm, &_rank));

    // Retrieve the extent of each dimension of interest. The dimensions of
    // interest are the spatial dimensions of the grid. These are named "x" and
    // "y" by default.
    // I have switched the order of the dimensions here to reflect the change of
    // indices in nextsim-dg
    size_t tmp_0, tmp_1;
    NC_CHECK(nc_inq_dimlen(data_nc_id, dim1_nc_id, &tmp_0));
    NC_CHECK(nc_inq_dimlen(data_nc_id, dim0_nc_id, &tmp_1));
    _global_ext_0 = static_cast<int>(tmp_0);
    _global_ext_1 = static_cast<int>(tmp_1);

    // Initially we partition assuming there is no land mask
    // Figure out my subset of objects
    CHECK_MPI(MPI_Comm_size(comm, &_num_procs));

    // Start from a 2D decomposition
    _num_procs_0 = _num_procs;
    _num_procs_1 = 1;
    find_factors(_num_procs, _num_procs_0, _num_procs_1);

    _local_ext_0 = ceil((float)_global_ext_0 / (_num_procs_0));
    _local_ext_1 = ceil((float)_global_ext_1 / (_num_procs_1));
    _global_0 = (_rank / _num_procs_1) * _local_ext_0;
    _global_1 = (_rank % _num_procs_1) * _local_ext_1;

    if ((_rank / _num_procs_1) == _num_procs_0 - 1) {
        _local_ext_0 = _global_ext_0 - (_rank / _num_procs_1) * _local_ext_0;
    }
    if ((_rank % _num_procs_1) == _num_procs_1 - 1) {
        _local_ext_1 = _global_ext_1 - (_rank % _num_procs_1) * _local_ext_1;
    }

    _num_objects = _local_ext_0 * _local_ext_1;

    // Retrieve the land mask, if available and enabled
    int mask_nc_id;
    int nc_err;
    nc_err = nc_inq_varid(data_nc_id, mask_name.c_str(), &mask_nc_id);

    if (!ignore_mask && nc_err == NC_NOERR && nc_err != NC_ENOTVAR) {
        // Data reads are independent by default, so we need to switch to
        // collective for improved parallel I/O performance
        NC_CHECK(nc_var_par_access(data_nc_id, mask_nc_id, NC_COLLECTIVE));

        // Verify the order of dimensions provided is correct by comparing to the
        // dimension order of the mask variable
        const int NDIMS = 2;
        int dim_id[NDIMS];
        char dim_name[NDIMS][128];
        NC_CHECK(nc_inq_vardimid(data_nc_id, mask_nc_id, &dim_id[0]));
        // I have switched the order of the dimensions here to reflect the change of
        // indices in nextsim-dg
        NC_CHECK(nc_inq_dimname(data_nc_id, dim_id[1], &dim_name[0][0]));
        NC_CHECK(nc_inq_dimname(data_nc_id, dim_id[0], &dim_name[1][0]));
        if (dim_name[0] != dim0_name || dim_name[1] != dim1_name) {
            throw std::runtime_error("Dimension ordering provided does not match "
                                     "ordering in netCDF grid file");
        }

        _land_mask.resize(_num_objects);
        size_t start[NDIMS], count[NDIMS];
        // I have switched the order of the dimensions here to reflect the change of
        // indices in nextsim-dg
        // Coordinate of first element
        start[1] = _global_0;
        start[0] = _global_1;
        // Number of elements in every extension
        count[1] = _local_ext_0;
        count[0] = _local_ext_1;
        NC_CHECK(nc_get_vara_int(data_nc_id, mask_nc_id, start, count, _land_mask.data()));

        // create copy of land mask ready to transpose
        std::vector<int> _land_mask_copy(_land_mask);

        // transpose _land_mask to reflect the change of indices in nextsim-dg
        int index = 0;
        for (size_t j = 0; j < count[1]; j++) {
            for (size_t i = 0; i < count[0]; i++) {
                _land_mask[index] = _land_mask_copy[i * count[1] + j];
                index++;
            }
        }

        // Apply land mask
        for (int i = 0; i < _num_objects; i++) {
            // The convention is that sea data points will have a positive value
            // and land points a zero value
            if (_land_mask[i] > 0) {
                int local_0 = i / _local_ext_1;
                int local_1 = i % _local_ext_1;
                int global_0 = local_0 + _global_0;
                int global_1 = local_1 + _global_1;
                _object_id.push_back(global_0 * _global_ext_1 + global_1);
                _sparse_to_dense.push_back(i);
                _num_nonzero_objects++;
            }
        }
    } else {
        _num_nonzero_objects = _num_objects;
        for (int i = 0; i < _num_objects; i++) {
            int local_0 = i / _local_ext_1;
            int local_1 = i % _local_ext_1;
            int global_0 = local_0 + _global_0;
            int global_1 = local_1 + _global_1;
            _object_id.push_back(global_0 * _global_ext_1 + global_1);
        }
    }

    NC_CHECK(nc_close(nc_id));
}

int Grid::get_num_objects() const { return _num_objects; }

int Grid::get_num_nonzero_objects() const { return _num_nonzero_objects; }

int Grid::get_global_ext_0() const { return _global_ext_0; }

int Grid::get_global_ext_1() const { return _global_ext_1; }

int Grid::get_local_ext_0() const { return _local_ext_0; }

int Grid::get_local_ext_1() const { return _local_ext_1; }

int Grid::get_global_0() const { return _global_0; }

int Grid::get_global_1() const { return _global_1; }

int Grid::get_num_procs_0() const { return _num_procs_0; }

int Grid::get_num_procs_1() const { return _num_procs_1; }

const int* Grid::get_land_mask() const { return _land_mask.data(); }

const int* Grid::get_sparse_to_dense() const { return _sparse_to_dense.data(); }

const int* Grid::get_nonzero_object_ids() const { return _object_id.data(); }

void Grid::get_bounding_box(int& global_0, int& global_1, int& local_ext_0, int& local_ext_1) const
{
    global_0 = _global_0;
    global_1 = _global_1;
    local_ext_0 = _local_ext_0;
    local_ext_1 = _local_ext_1;
}
