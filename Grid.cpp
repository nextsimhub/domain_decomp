/*!
 * @file Grid.cpp
 * @date 05 Nov 2024
 * @author Athena Elafrou <ae488@cam.ac.uk>
 */

#include "Grid.hpp"
#include "Utils.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <netcdf.h>
#include <netcdf_par.h>
#include <vector>

static std::vector<int> find_factors(const int n)
{
    int factor_a = -1;
    int factor_b = -1;
    for (int i = 2; i * i <= n; i += 2) {
        if (n % i == 0) {
            factor_a = i;
            factor_b = n / factor_a;
        }
    }
    if (factor_a == -1 || factor_b == -1) {
        // if factor is not found then use 1D splitting
        return { n, 1 };
    } else {
        // if factor is found then use the 2D factor split
        return { factor_a, factor_b };
    }
}

Grid* Grid::create(
    MPI_Comm comm, const std::string& filename, bool ignore_mask, bool px, bool py, bool tripolar)
{
    return new Grid(comm, filename, "x", "y", std::vector<int>({ 1, 0 }), "mask", ignore_mask, px,
        py, tripolar);
}

Grid* Grid::create(MPI_Comm comm, const std::string& filename, const std::string xdim_name,
    const std::string ydim_name, const std::vector<int> dim_order, const std::string mask_name,
    bool ignore_mask, bool px, bool py, bool tripolar)
{
    return new Grid(
        comm, filename, xdim_name, ydim_name, dim_order, mask_name, ignore_mask, px, py, tripolar);
}

void Grid::ReadGridExtents(const std::string& filename)
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

    int dim_ids[NDIMS];
    size_t tmp[NDIMS];
    for (int i = 0; i < NDIMS; i++) {
        // get the dimension id of each dimension
        NC_CHECK(nc_inq_dimid(data_nc_id, _dim_names[_dim_order[i]].c_str(), &dim_ids[i]));
        // get the extent of each dimension
        NC_CHECK(nc_inq_dimlen(data_nc_id, dim_ids[i], &tmp[i]));
        _globalExt[_dim_order[i]] = static_cast<int>(tmp[i]);
    }

    NC_CHECK(nc_close(nc_id));
}

void Grid::ReadGridMask(const std::string& filename, const std::string& mask_name)
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

    // Retrieve the land mask, if available and enabled
    int mask_nc_id;
    int nc_err;
    nc_err = nc_inq_varid(data_nc_id, mask_name.c_str(), &mask_nc_id);

    if (nc_err == NC_NOERR && nc_err != NC_ENOTVAR) {

        // Data reads are independent by default, so we need to switch to
        // collective for improved parallel I/O performance
        NC_CHECK(nc_var_par_access(data_nc_id, mask_nc_id, NC_COLLECTIVE));

        int dim_id[NDIMS];
        char dim_name[NDIMS][128];

        NC_CHECK(nc_inq_vardimid(data_nc_id, mask_nc_id, &dim_id[0]));

        for (int i = 0; i < NDIMS; i++) {
            // Verify the order of dimensions provided is correct by comparing
            // to the dimension order of the mask variable
            NC_CHECK(nc_inq_dimname(data_nc_id, dim_id[i], &dim_name[i][0]));
            if (std::string(dim_name[i]) != _dim_names[_dim_order[i]]) {
                throw std::runtime_error(
                    "Dimension ordering provided does not match ordering in netCDF grid file");
            }
        }

        _land_mask.resize(_num_objects);
        size_t start[NDIMS], count[NDIMS];

        for (int i = 0; i < NDIMS; i++) {
            // Position of first element
            start[_dim_order[i]] = static_cast<size_t>(_global[i]);
            // Number of elements to read
            count[_dim_order[i]] = static_cast<size_t>(_localExt[i]);
        }

        NC_CHECK(nc_get_vara_int(data_nc_id, mask_nc_id, start, count, _land_mask.data()));
    }

    NC_CHECK(nc_close(nc_id));
}

Grid::Grid(MPI_Comm comm, const std::string& filename, const std::string& xdim_name,
    const std::string& ydim_name, const std::vector<int>& dim_order, const std::string& mask_name,
    bool ignore_mask, bool px, bool py, bool tripolar)
    : _comm(comm)
    , _dim_names({ xdim_name, ydim_name })
    , _dim_order(dim_order)
    , _px(px)
    , _py(py)
    , _tripolar(tripolar)
{

    // Tripolar grid cannot be periodic
    if (_tripolar && (_px || _py)) {
        // TODO:
        //  - Improve error message
        //  - Verify conventions our for error handling
        throw std::invalid_argument(
            "Tripolar topology and periodic BCs are toggled at the same time");
    }

    CHECK_MPI(MPI_Comm_rank(comm, &_rank));

    Grid::ReadGridExtents(filename);

    CHECK_MPI(MPI_Comm_size(comm, &_totalNumProcs));

    // Initially we partition assuming there is no land mask
    // Start from a naive 2D decomposition
    _numProcs = find_factors(_totalNumProcs);

    // split into chunks
    for (size_t i = 0; i < NDIMS; i++) {
        _localExt[i] = ceil((float)_globalExt[i] / (_numProcs[i]));
    }

    // need to account for residuals
    _global[0] = (_rank / _numProcs[1]) * _localExt[0];
    _global[1] = (_rank % _numProcs[1]) * _localExt[1];

    if ((_rank / _numProcs[1]) == _numProcs[0] - 1) {
        _localExt[0] = _globalExt[0] - (_rank / _numProcs[1]) * _localExt[0];
    }
    if ((_rank % _numProcs[1]) == _numProcs[1] - 1) {
        _localExt[1] = _globalExt[1] - (_rank % _numProcs[1]) * _localExt[1];
    }

    // number of objects for this process
    _num_objects = _localExt[0] * _localExt[1];

    if (!ignore_mask) {

        Grid::ReadGridMask(filename, mask_name);

        // Apply land mask
        for (int i = 0; i < _num_objects; i++) {
            // The convention is that sea data points will have a positive value
            // and land points a zero value
            if (_land_mask[i] > 0) {
                // find index position in the global grid
                int temp_i = i % _localExt[0] + _global[0];
                int temp_j = i / _localExt[0] + _global[1];
                _global_id.push_back(temp_j * _globalExt[0] + temp_i);
                // store local id for mapping
                _local_id.push_back(i);
                _num_nonzero_objects++;
            }
        }
    } else {
        _num_nonzero_objects = _num_objects;
        for (int i = 0; i < _num_objects; i++) {
            int temp_i = i % _localExt[0];
            int temp_j = i / _localExt[0];
            _global_id.push_back(temp_j * _globalExt[0] + _global[0] + temp_i);
            _local_id.push_back(i);
        }
    }
}

int Grid::get_num_objects() const { return _num_objects; }

int Grid::get_num_nonzero_objects() const { return _num_nonzero_objects; }

bool Grid::get_px() const { return _px; }

bool Grid::get_py() const { return _py; }

bool Grid::get_tripolar() const { return _tripolar; }

std::vector<int> Grid::getNumProcs() const { return _numProcs; }

std::vector<int> Grid::getGlobalExt() const { return _globalExt; }

std::vector<int> Grid::getLocalExt() const { return _localExt; }

std::vector<int> Grid::get_global() const { return _global; }

const int* Grid::get_land_mask() const { return _land_mask.data(); }

const int* Grid::get_sparse_to_dense() const { return _local_id.data(); }

const int* Grid::get_nonzero_object_ids() const { return _global_id.data(); }

void Grid::get_bounding_box(int& global0, int& global1, int& localExt0, int& localExt1) const
{
    global0 = _global[0];
    global1 = _global[1];
    localExt0 = _localExt[0];
    localExt1 = _localExt[1];
}
