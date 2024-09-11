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
    global_0 = _global_0_new;
    global_1 = _global_1_new;
    local_ext_0 = _local_ext_0_new;
    local_ext_1 = _local_ext_1_new;
}

void Partitioner::get_neighbours(std::vector<int>& ids, std::vector<int>& halo_sizes, int idx) const
{
    for (auto it = _neighbours[idx].begin(); it != _neighbours[idx].end(); ++it) {
        ids.push_back(it->first);
        halo_sizes.push_back(it->second);
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
    NC_CHECK(nc_create_par(filename.c_str(), nc_mode, _comm, MPI_INFO_NULL, &nc_id));

    // Prepare neighbour data
    std::vector<int> top_ids, bottom_ids, left_ids, right_ids;
    std::vector<int> top_halos, bottom_halos, left_halos, right_halos;
    get_neighbours(left_ids, left_halos, 0);
    get_neighbours(right_ids, right_halos, 1);
    get_neighbours(bottom_ids, bottom_halos, 2);
    get_neighbours(top_ids, top_halos, 3);
    std::vector<int> num_neighbours = { (int)left_ids.size(), (int)right_ids.size(),
        (int)bottom_ids.size(), (int)top_ids.size() };

    // Compute global dimensions
    std::vector<int> dims = { 0, 0, 0, 0 };
    CHECK_MPI(MPI_Allreduce(&num_neighbours[0], &dims[0], 1, MPI_INT, MPI_SUM, _comm));
    CHECK_MPI(MPI_Allreduce(&num_neighbours[1], &dims[1], 1, MPI_INT, MPI_SUM, _comm));
    CHECK_MPI(MPI_Allreduce(&num_neighbours[2], &dims[2], 1, MPI_INT, MPI_SUM, _comm));
    CHECK_MPI(MPI_Allreduce(&num_neighbours[3], &dims[3], 1, MPI_INT, MPI_SUM, _comm));

    // Compute global offsets
    std::vector<int> offsets = { 0, 0, 0, 0 };
    CHECK_MPI(MPI_Exscan(&num_neighbours[0], &offsets[0], 1, MPI_INT, MPI_SUM, _comm));
    CHECK_MPI(MPI_Exscan(&num_neighbours[1], &offsets[1], 1, MPI_INT, MPI_SUM, _comm));
    CHECK_MPI(MPI_Exscan(&num_neighbours[2], &offsets[2], 1, MPI_INT, MPI_SUM, _comm));
    CHECK_MPI(MPI_Exscan(&num_neighbours[3], &offsets[3], 1, MPI_INT, MPI_SUM, _comm));

    // Create 2 dimensions
    // The values to be written are associated with the netCDF variable by
    // assuming that the last dimension of the netCDF variable varies fastest in
    // the C interface
    const int NDIMS = 2;
    int dimid_global[NDIMS];
    NC_CHECK(nc_def_dim(nc_id, "NX", _global_ext_0, &dimid_global[0]));
    NC_CHECK(nc_def_dim(nc_id, "NY", _global_ext_1, &dimid_global[1]));

    // Define dimensions in netCDF file
    int dimid, top_dimid, bottom_dimid, left_dimid, right_dimid;
    NC_CHECK(nc_def_dim(nc_id, "P", _num_procs, &dimid));
    NC_CHECK(nc_def_dim(nc_id, "T", dims[3], &top_dimid));
    NC_CHECK(nc_def_dim(nc_id, "B", dims[2], &bottom_dimid));
    NC_CHECK(nc_def_dim(nc_id, "L", dims[0], &left_dimid));
    NC_CHECK(nc_def_dim(nc_id, "R", dims[1], &right_dimid));

    // Define groups in netCDF file
    int bbox_gid, connectivity_gid;
    NC_CHECK(nc_def_grp(nc_id, "bounding_boxes", &bbox_gid));
    NC_CHECK(nc_def_grp(nc_id, "connectivity", &connectivity_gid));

    // Define variables in netCDF file
    int top_x_vid, top_y_vid;
    int cnt_x_vid, cnt_y_vid;
    int top_num_vid, bottom_num_vid, left_num_vid, right_num_vid;
    int top_ids_vid, bottom_ids_vid, left_ids_vid, right_ids_vid;
    int top_halos_vid, bottom_halos_vid, left_halos_vid, right_halos_vid;
    // Bounding boxes group
    NC_CHECK(nc_def_var(bbox_gid, "domain_x", NC_INT, 1, &dimid, &top_x_vid));
    NC_CHECK(nc_def_var(bbox_gid, "domain_y", NC_INT, 1, &dimid, &top_y_vid));
    NC_CHECK(nc_def_var(bbox_gid, "domain_extent_x", NC_INT, 1, &dimid, &cnt_x_vid));
    NC_CHECK(nc_def_var(bbox_gid, "domain_extent_y", NC_INT, 1, &dimid, &cnt_y_vid));
    // Connectivity group
    NC_CHECK(nc_def_var(connectivity_gid, "top_neighbours", NC_INT, 1, &dimid, &top_num_vid));
    NC_CHECK(
        nc_def_var(connectivity_gid, "top_neighbour_ids", NC_INT, 1, &top_dimid, &top_ids_vid));
    NC_CHECK(
        nc_def_var(connectivity_gid, "top_neighbour_halos", NC_INT, 1, &top_dimid, &top_halos_vid));
    NC_CHECK(nc_def_var(connectivity_gid, "bottom_neighbours", NC_INT, 1, &dimid, &bottom_num_vid));
    NC_CHECK(nc_def_var(
        connectivity_gid, "bottom_neighbour_ids", NC_INT, 1, &bottom_dimid, &bottom_ids_vid));
    NC_CHECK(nc_def_var(
        connectivity_gid, "bottom_neighbour_halos", NC_INT, 1, &bottom_dimid, &bottom_halos_vid));
    NC_CHECK(nc_def_var(connectivity_gid, "left_neighbours", NC_INT, 1, &dimid, &left_num_vid));
    NC_CHECK(
        nc_def_var(connectivity_gid, "left_neighbour_ids", NC_INT, 1, &left_dimid, &left_ids_vid));
    NC_CHECK(nc_def_var(
        connectivity_gid, "left_neighbour_halos", NC_INT, 1, &left_dimid, &left_halos_vid));
    NC_CHECK(nc_def_var(connectivity_gid, "right_neighbours", NC_INT, 1, &dimid, &right_num_vid));
    NC_CHECK(nc_def_var(
        connectivity_gid, "right_neighbour_ids", NC_INT, 1, &right_dimid, &right_ids_vid));
    NC_CHECK(nc_def_var(
        connectivity_gid, "right_neighbour_halos", NC_INT, 1, &right_dimid, &right_halos_vid));

    // Write metadata to file
    NC_CHECK(nc_enddef(nc_id));

    // Set up slab for this process
    size_t start, count;

    // Store data
    start = _rank;
    NC_CHECK(nc_var_par_access(bbox_gid, top_x_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_var1_int(bbox_gid, top_x_vid, &start, &_global_0_new));
    NC_CHECK(nc_var_par_access(bbox_gid, top_y_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_var1_int(bbox_gid, top_y_vid, &start, &_global_1_new));
    NC_CHECK(nc_var_par_access(bbox_gid, cnt_x_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_var1_int(bbox_gid, cnt_x_vid, &start, &_local_ext_0_new));
    NC_CHECK(nc_var_par_access(bbox_gid, cnt_y_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_var1_int(bbox_gid, cnt_y_vid, &start, &_local_ext_1_new));

    NC_CHECK(nc_var_par_access(connectivity_gid, top_num_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_var1_int(connectivity_gid, top_num_vid, &start, &num_neighbours[3]));
    NC_CHECK(nc_var_par_access(connectivity_gid, bottom_num_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_var1_int(connectivity_gid, bottom_num_vid, &start, &num_neighbours[2]));
    NC_CHECK(nc_var_par_access(connectivity_gid, left_num_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_var1_int(connectivity_gid, left_num_vid, &start, &num_neighbours[0]));
    NC_CHECK(nc_var_par_access(connectivity_gid, right_num_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_var1_int(connectivity_gid, right_num_vid, &start, &num_neighbours[1]));

    start = offsets[0];
    count = num_neighbours[0];
    NC_CHECK(nc_var_par_access(connectivity_gid, left_ids_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_vara_int(connectivity_gid, left_ids_vid, &start, &count, left_ids.data()));
    NC_CHECK(nc_var_par_access(connectivity_gid, left_halos_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_vara_int(connectivity_gid, left_halos_vid, &start, &count, left_halos.data()));
    start = offsets[1];
    count = num_neighbours[1];
    NC_CHECK(nc_var_par_access(connectivity_gid, right_ids_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_vara_int(connectivity_gid, right_ids_vid, &start, &count, right_ids.data()));
    NC_CHECK(nc_var_par_access(connectivity_gid, right_halos_vid, NC_COLLECTIVE));
    NC_CHECK(
        nc_put_vara_int(connectivity_gid, right_halos_vid, &start, &count, right_halos.data()));
    start = offsets[2];
    count = num_neighbours[2];
    NC_CHECK(nc_var_par_access(connectivity_gid, bottom_ids_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_vara_int(connectivity_gid, bottom_ids_vid, &start, &count, bottom_ids.data()));
    NC_CHECK(nc_var_par_access(connectivity_gid, bottom_halos_vid, NC_COLLECTIVE));
    NC_CHECK(
        nc_put_vara_int(connectivity_gid, bottom_halos_vid, &start, &count, bottom_halos.data()));
    start = offsets[3];
    count = num_neighbours[3];
    NC_CHECK(nc_var_par_access(connectivity_gid, top_ids_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_vara_int(connectivity_gid, top_ids_vid, &start, &count, top_ids.data()));
    NC_CHECK(nc_var_par_access(connectivity_gid, top_halos_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_vara_int(connectivity_gid, top_halos_vid, &start, &count, top_halos.data()));

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
    std::vector<int> top_left_0(_num_procs, -1);
    std::vector<int> top_left_1(_num_procs, -1);
    std::vector<int> top_right_0(_num_procs, -1);
    std::vector<int> top_right_1(_num_procs, -1);
    std::vector<int> bottom_left_0(_num_procs, -1);
    std::vector<int> bottom_left_1(_num_procs, -1);
    std::vector<int> bottom_right_0(_num_procs, -1);
    std::vector<int> bottom_right_1(_num_procs, -1);
    CHECK_MPI(MPI_Allgather(&_global_0_new, 1, MPI_INT, top_left_0.data(), 1, MPI_INT, _comm));
    CHECK_MPI(MPI_Allgather(&_global_1_new, 1, MPI_INT, top_left_1.data(), 1, MPI_INT, _comm));
    CHECK_MPI(
        MPI_Allgather(&_local_ext_0_new, 1, MPI_INT, bottom_right_0.data(), 1, MPI_INT, _comm));
    CHECK_MPI(
        MPI_Allgather(&_local_ext_1_new, 1, MPI_INT, bottom_right_1.data(), 1, MPI_INT, _comm));
    for (int i = 0; i < _num_procs; i++)
        top_right_0[i] = top_left_0[i];
    for (int i = 0; i < _num_procs; i++)
        top_right_1[i] = top_left_1[i] + bottom_right_1[i] - 1;
    for (int i = 0; i < _num_procs; i++)
        bottom_left_0[i] = top_left_0[i] + bottom_right_0[i] - 1;
    for (int i = 0; i < _num_procs; i++)
        bottom_left_1[i] = top_left_1[i];
    for (int i = 0; i < _num_procs; i++)
        bottom_right_0[i] += top_left_0[i] - 1;
    for (int i = 0; i < _num_procs; i++)
        bottom_right_1[i] += top_left_1[i] - 1;

    // Find my top neighbours and their halo sizes
    for (int i = 0; i < _num_procs; i++) {
        if (i != _rank) {
            if (top_left_1[_rank] >= bottom_left_1[i] && top_left_1[_rank] <= bottom_right_1[i]
                && bottom_right_1[i] <= top_right_1[_rank]
                && (top_left_0[_rank] - bottom_left_0[i] == 1)) {
                int halo_size = bottom_right_1[i] - top_left_1[_rank] + 1;
                _neighbours[3].insert(std::pair<int, int>(i, halo_size));
            }
            if (top_right_1[_rank] >= bottom_left_1[i] && top_right_1[_rank] <= bottom_right_1[i]
                && bottom_left_1[i] >= top_left_1[_rank]
                && (top_right_0[_rank] - bottom_right_0[i] == 1)) {
                int halo_size = top_right_1[_rank] - bottom_left_1[i] + 1;
                _neighbours[3].insert(std::pair<int, int>(i, halo_size));
            }
        }
    }

    // Find my bottom neighbours
    for (int i = 0; i < _num_procs; i++) {
        if (i != _rank) {
            if (bottom_left_1[_rank] >= top_left_1[i] && bottom_left_1[_rank] <= top_right_1[i]
                && top_right_1[i] <= bottom_right_1[_rank]
                && (top_left_0[i] - bottom_left_0[_rank] == 1)) {
                int halo_size = top_right_1[i] - bottom_left_1[_rank] + 1;
                _neighbours[2].insert(std::pair<int, int>(i, halo_size));
            }
            if (bottom_right_1[_rank] >= top_left_1[i] && bottom_right_1[_rank] <= top_right_1[i]
                && top_left_1[i] >= bottom_left_1[_rank]
                && (top_right_0[i] - bottom_right_0[_rank] == 1)) {
                int halo_size = bottom_right_1[_rank] - top_left_1[i] + 1;
                _neighbours[2].insert(std::pair<int, int>(i, halo_size));
            }
        }
    }

    // Find my left neighbours
    for (int i = 0; i < _num_procs; i++) {
        if (i != _rank) {
            if (top_left_0[_rank] >= top_right_0[i] && top_left_0[_rank] <= bottom_right_0[i]
                && bottom_left_0[_rank] <= bottom_right_0[i]
                && (top_left_1[_rank] - top_right_1[i] == 1)) {
                int halo_size = bottom_right_0[i] - top_left_0[_rank] + 1;
                _neighbours[0].insert(std::pair<int, int>(i, halo_size));
            }
            if (bottom_left_0[_rank] >= top_right_0[i] && bottom_left_0[_rank] <= bottom_right_0[i]
                && top_left_0[_rank] <= top_right_0[i]
                && (bottom_left_1[_rank] - top_right_1[i] == 1)) {
                int halo_size = bottom_left_0[_rank] - top_right_0[i] + 1;
                _neighbours[0].insert(std::pair<int, int>(i, halo_size));
            }
        }
    }

    // Find my right neighbours
    for (int i = 0; i < _num_procs; i++) {
        if (i != _rank) {
            if (top_right_0[_rank] >= top_left_0[i] && top_right_0[_rank] <= bottom_left_0[i]
                && bottom_right_0[_rank] >= bottom_left_0[i]
                && (top_left_1[i] - top_right_1[_rank] == 1)) {
                int halo_size = bottom_left_0[i] - top_right_0[_rank] + 1;
                _neighbours[1].insert(std::pair<int, int>(i, halo_size));
            }
            if (bottom_right_0[_rank] >= top_left_0[i] && bottom_right_0[_rank] <= bottom_left_0[i]
                && top_right_0[_rank] <= top_left_0[i]
                && (top_left_1[i] - top_right_1[_rank] == 1)) {
                int halo_size = bottom_right_0[_rank] - top_left_0[i] + 1;
                _neighbours[1].insert(std::pair<int, int>(i, halo_size));
            }
        }
    }
}
