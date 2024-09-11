/*!
 * @file Partitioner.cpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @date 11 Sep 2024
 */

#include "Partitioner.hpp"
#include "Utils.hpp"
#include "ZoltanPartitioner.hpp"

#include <cmath>
#include <stdexcept>

#include <netcdf.h>
#include <netcdf_par.h>

Partitioner::Partitioner(MPI_Comm comm)
{
    _comm = comm;
    CHECK_MPI(MPI_Comm_size(comm, &_num_procs));
    CHECK_MPI(MPI_Comm_rank(comm, &_rank));
}

void Partitioner::get_bounding_box(
    int& global_0, int& global_1, int& local_ext_0, int& local_ext_1) const
{
    global_0 = _global_new[0];
    global_1 = _global_new[1];
    local_ext_0 = _local_ext_new[0];
    local_ext_1 = _local_ext_new[1];
}

void Partitioner::get_neighbours(
    std::vector<std::vector<int>>& ids, std::vector<std::vector<int>>& halo_sizes) const
{
    for (int idx = 0; idx < 4; idx++) {
        for (auto it = _neighbours[idx].begin(); it != _neighbours[idx].end(); ++it) {
            ids[idx].push_back(it->first);
            halo_sizes[idx].push_back(it->second);
        }
    }
}

void Partitioner::save_mask(const std::string& filename) const
{
    // Use C API for parallel I/O
    int nc_id, nc_mode;
    nc_mode = NC_CLOBBER | NC_NETCDF4;
    NC_CHECK(nc_create_par(filename.c_str(), nc_mode, _comm, MPI_INFO_NULL, &nc_id));
    NC_CHECK(nc_put_att_int(nc_id, NC_GLOBAL, "num_processes", NC_INT, 1, &_num_procs));

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
    NC_CHECK(nc_def_var(nc_id, "pid", NC_INT, NDIMS, dimid, &mask_nc_id));

    // Write metadata to file
    NC_CHECK(nc_enddef(nc_id));

    // Set up slab for this process
    size_t start[NDIMS], count[NDIMS];
    start[0] = _global[0];
    start[1] = _global[1];
    count[0] = _local_ext[0];
    count[1] = _local_ext[1];

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
    NC_CHECK(nc_create_par(filename.c_str(), nc_mode, _comm, MPI_INFO_NULL, &nc_id));

    // Prepare neighbour data
    std::vector<std::vector<int>> ids = { {}, {}, {}, {} };
    std::vector<std::vector<int>> halos = { {}, {}, {}, {} };
    get_neighbours(ids, halos);
    std::vector<int> num_neighbours
        = { (int)ids[0].size(), (int)ids[1].size(), (int)ids[2].size(), (int)ids[3].size() };

    // Compute global dimensions and offsets
    std::vector<int> dims = { 0, 0, 0, 0 }, offsets = { 0, 0, 0, 0 };
    for (int idx = 0; idx < 4; idx++) {
        CHECK_MPI(MPI_Allreduce(&num_neighbours[idx], &dims[idx], 1, MPI_INT, MPI_SUM, _comm));
        CHECK_MPI(MPI_Exscan(&num_neighbours[idx], &offsets[idx], 1, MPI_INT, MPI_SUM, _comm));
    }

    // Create 2 dimensions
    // The values to be written are associated with the netCDF variable by
    // assuming that the last dimension of the netCDF variable varies fastest in
    // the C interface
    const int NDIMS = 2;
    int dimid_global[NDIMS];
    NC_CHECK(nc_def_dim(nc_id, "NX", _global_ext_0, &dimid_global[0]));
    NC_CHECK(nc_def_dim(nc_id, "NY", _global_ext_1, &dimid_global[1]));

    // Define dimensions in netCDF file
    int dimid;
    std::vector<int> dimids(4);
    NC_CHECK(nc_def_dim(nc_id, "P", _num_procs, &dimid));
    std::vector<std::string> dim_letters = { "L", "R", "B", "T" };
    for (int idx = 0; idx < 4; idx++) {
        NC_CHECK(nc_def_dim(nc_id, dim_letters[idx].c_str(), dims[idx], &dimids[idx]));
    }

    // Define groups in netCDF file
    int bbox_gid, connectivity_gid;
    NC_CHECK(nc_def_grp(nc_id, "bounding_boxes", &bbox_gid));
    NC_CHECK(nc_def_grp(nc_id, "connectivity", &connectivity_gid));

    // Define variables in netCDF file
    int top_x_vid, top_y_vid;
    int cnt_x_vid, cnt_y_vid;
    std::vector<int> num_vid(4);
    std::vector<int> ids_vid(4);
    std::vector<int> halos_vid(4);
    std::vector<std::string> dim_names = { "left", "right", "bottom", "top" };
    // Bounding boxes group
    NC_CHECK(nc_def_var(bbox_gid, "domain_x", NC_INT, 1, &dimid, &top_x_vid));
    NC_CHECK(nc_def_var(bbox_gid, "domain_y", NC_INT, 1, &dimid, &top_y_vid));
    NC_CHECK(nc_def_var(bbox_gid, "domain_extent_x", NC_INT, 1, &dimid, &cnt_x_vid));
    NC_CHECK(nc_def_var(bbox_gid, "domain_extent_y", NC_INT, 1, &dimid, &cnt_y_vid));
    // Connectivity group
    for (int idx = 0; idx < 4; idx++) {
        NC_CHECK(nc_def_var(connectivity_gid, (dim_names[idx] + "_neighbours").c_str(), NC_INT, 1,
            &dimid, &num_vid[idx]));
        NC_CHECK(nc_def_var(connectivity_gid, (dim_names[idx] + "_neighbour_ids").c_str(), NC_INT,
            1, &dimids[idx], &ids_vid[idx]));
        NC_CHECK(nc_def_var(connectivity_gid, (dim_names[idx] + "_neighbour_halos").c_str(), NC_INT,
            1, &dimids[idx], &halos_vid[idx]));
    }

    // Write metadata to file
    NC_CHECK(nc_enddef(nc_id));

    // Set up slab for this process
    size_t start, count;

    // Store data
    start = _rank;
    NC_CHECK(nc_var_par_access(bbox_gid, top_x_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_var1_int(bbox_gid, top_x_vid, &start, &_global_new[0]));
    NC_CHECK(nc_var_par_access(bbox_gid, top_y_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_var1_int(bbox_gid, top_y_vid, &start, &_global_new[1]));
    NC_CHECK(nc_var_par_access(bbox_gid, cnt_x_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_var1_int(bbox_gid, cnt_x_vid, &start, &_local_ext_new[0]));
    NC_CHECK(nc_var_par_access(bbox_gid, cnt_y_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_var1_int(bbox_gid, cnt_y_vid, &start, &_local_ext_new[1]));
    for (int idx = 0; idx < 4; idx++) {
        NC_CHECK(nc_var_par_access(connectivity_gid, num_vid[idx], NC_COLLECTIVE));
        NC_CHECK(nc_put_var1_int(connectivity_gid, num_vid[idx], &start, &num_neighbours[idx]));
        start = offsets[idx];
        count = num_neighbours[idx];
        NC_CHECK(nc_var_par_access(connectivity_gid, ids_vid[idx], NC_COLLECTIVE));
        NC_CHECK(nc_put_vara_int(connectivity_gid, ids_vid[idx], &start, &count, ids[idx].data()));
        NC_CHECK(nc_var_par_access(connectivity_gid, halos_vid[idx], NC_COLLECTIVE));
        NC_CHECK(
            nc_put_vara_int(connectivity_gid, halos_vid[idx], &start, &count, halos[idx].data()));
    }
    NC_CHECK(nc_close(nc_id));
}

Partitioner* Partitioner::Factory::create(
    MPI_Comm comm, int argc, char** argv, PartitionerType type)
{
    if (type == PartitionerType::Zoltan_RCB)
        return ZoltanPartitioner::create(comm, argc, argv);
    else
        throw std::runtime_error("Invalid partitioner!");
}

void Partitioner::discover_neighbours()
{
    // Gather bounding boxes for all processes
    std::vector<std::vector<int>> top_left
        = { std::vector<int>(_num_procs, -1), std::vector<int>(_num_procs, -1) };
    std::vector<std::vector<int>> top_right
        = { std::vector<int>(_num_procs, -1), std::vector<int>(_num_procs, -1) };
    std::vector<std::vector<int>> bottom_left
        = { std::vector<int>(_num_procs, -1), std::vector<int>(_num_procs, -1) };
    std::vector<std::vector<int>> bottom_right
        = { std::vector<int>(_num_procs, -1), std::vector<int>(_num_procs, -1) };
    for (int idx = 0; idx < 1; idx++) {
        CHECK_MPI(
            MPI_Allgather(&_global_new[idx], 1, MPI_INT, top_left[idx].data(), 1, MPI_INT, _comm));
        CHECK_MPI(MPI_Allgather(
            &_local_ext_new[idx], 1, MPI_INT, bottom_right[idx].data(), 1, MPI_INT, _comm));
    }
    for (int i = 0; i < _num_procs; i++)
        top_right[0][i] = top_left[0][i];
    for (int i = 0; i < _num_procs; i++)
        top_right[1][i] = top_left[1][i] + bottom_right[1][i] - 1;
    for (int i = 0; i < _num_procs; i++)
        bottom_left[0][i] = top_left[0][i] + bottom_right[0][i] - 1;
    for (int i = 0; i < _num_procs; i++)
        bottom_left[1][i] = top_left[1][i];
    for (int i = 0; i < _num_procs; i++)
        bottom_right[0][i] += top_left[0][i] - 1;
    for (int i = 0; i < _num_procs; i++)
        bottom_right[1][i] += top_left[1][i] - 1;

    // Find my top neighbours and their halo sizes
    for (int i = 0; i < _num_procs; i++) {
        if (i != _rank) {
            if (top_left[1][_rank] >= bottom_left[1][i] && top_left[1][_rank] <= bottom_right[1][i]
                && bottom_right[1][i] <= top_right[1][_rank]
                && (top_left[0][_rank] - bottom_left[0][i] == 1)) {
                int halo_size = bottom_right[1][i] - top_left[1][_rank] + 1;
                _neighbours[3].insert(std::pair<int, int>(i, halo_size));
            }
            if (top_right[1][_rank] >= bottom_left[1][i]
                && top_right[1][_rank] <= bottom_right[1][i]
                && bottom_left[1][i] >= top_left[1][_rank]
                && (top_right[0][_rank] - bottom_right[0][i] == 1)) {
                int halo_size = top_right[1][_rank] - bottom_left[1][i] + 1;
                _neighbours[3].insert(std::pair<int, int>(i, halo_size));
            }
        }
    }

    // Find my bottom neighbours
    for (int i = 0; i < _num_procs; i++) {
        if (i != _rank) {
            if (bottom_left[1][_rank] >= top_left[1][i] && bottom_left[1][_rank] <= top_right[1][i]
                && top_right[1][i] <= bottom_right[1][_rank]
                && (top_left[0][i] - bottom_left[0][_rank] == 1)) {
                int halo_size = top_right[1][i] - bottom_left[1][_rank] + 1;
                _neighbours[2].insert(std::pair<int, int>(i, halo_size));
            }
            if (bottom_right[1][_rank] >= top_left[1][i]
                && bottom_right[1][_rank] <= top_right[1][i]
                && top_left[1][i] >= bottom_left[1][_rank]
                && (top_right[0][i] - bottom_right[0][_rank] == 1)) {
                int halo_size = bottom_right[1][_rank] - top_left[1][i] + 1;
                _neighbours[2].insert(std::pair<int, int>(i, halo_size));
            }
        }
    }

    // Find my left neighbours
    for (int i = 0; i < _num_procs; i++) {
        if (i != _rank) {
            if (top_left[0][_rank] >= top_right[0][i] && top_left[0][_rank] <= bottom_right[0][i]
                && bottom_left[0][_rank] <= bottom_right[0][i]
                && (top_left[1][_rank] - top_right[1][i] == 1)) {
                int halo_size = bottom_right[0][i] - top_left[0][_rank] + 1;
                _neighbours[0].insert(std::pair<int, int>(i, halo_size));
            }
            if (bottom_left[0][_rank] >= top_right[0][i]
                && bottom_left[0][_rank] <= bottom_right[0][i]
                && top_left[0][_rank] <= top_right[0][i]
                && (bottom_left[1][_rank] - top_right[1][i] == 1)) {
                int halo_size = bottom_left[0][_rank] - top_right[0][i] + 1;
                _neighbours[0].insert(std::pair<int, int>(i, halo_size));
            }
        }
    }

    // Find my right neighbours
    for (int i = 0; i < _num_procs; i++) {
        if (i != _rank) {
            if (top_right[0][_rank] >= top_left[0][i] && top_right[0][_rank] <= bottom_left[0][i]
                && bottom_right[0][_rank] >= bottom_left[0][i]
                && (top_left[1][i] - top_right[1][_rank] == 1)) {
                int halo_size = bottom_left[0][i] - top_right[0][_rank] + 1;
                _neighbours[1].insert(std::pair<int, int>(i, halo_size));
            }
            if (bottom_right[0][_rank] >= top_left[0][i]
                && bottom_right[0][_rank] <= bottom_left[0][i]
                && top_right[0][_rank] <= top_left[0][i]
                && (top_left[1][i] - top_right[1][_rank] == 1)) {
                int halo_size = bottom_right[0][_rank] - top_left[0][i] + 1;
                _neighbours[1].insert(std::pair<int, int>(i, halo_size));
            }
        }
    }
}
