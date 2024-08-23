/*!
 * @file Partitioner.cpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @date 23 Aug 2024
 */

#include "Partitioner.hpp"
#include "Utils.hpp"
#include "ZoltanPartitioner.hpp"

#include <cmath>
#include <iostream> // TODO: Remove
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

void Partitioner::get_top_neighbours(std::vector<int>& ids, std::vector<int>& halo_sizes) const
{
    for (auto it = _top_neighbours.begin(); it != _top_neighbours.end(); ++it) {
        ids.push_back(it->first);
        halo_sizes.push_back(it->second);
    }
}

void Partitioner::get_bottom_neighbours(std::vector<int>& ids, std::vector<int>& halo_sizes) const
{
    for (auto it = _bottom_neighbours.begin(); it != _bottom_neighbours.end(); ++it) {
        ids.push_back(it->first);
        halo_sizes.push_back(it->second);
    }
}

void Partitioner::get_left_neighbours(std::vector<int>& ids, std::vector<int>& halo_sizes) const
{
    for (auto it = _left_neighbours.begin(); it != _left_neighbours.end(); ++it) {
        ids.push_back(it->first);
        halo_sizes.push_back(it->second);
    }
}

void Partitioner::get_right_neighbours(std::vector<int>& ids, std::vector<int>& halo_sizes) const
{
    for (auto it = _right_neighbours.begin(); it != _right_neighbours.end(); ++it) {
        ids.push_back(it->first);
        halo_sizes.push_back(it->second);
    }
}

void Partitioner::get_top_neighbours_periodic(
    std::vector<int>& ids, std::vector<int>& halo_sizes) const
{
    if (_p0) {
        // TODO: remove stdout
        std::cout << "_top_neighbour_ids_periodic.data(): " << std::endl;
        for (auto it = _top_neighbours_periodic.begin(); it != _top_neighbours_periodic.end();
             ++it) {
            std::cout << it->first << ", ";
        }
        std::cout << std::endl;
        std::cout << "_top_neighbour_halo_sizes_periodic.data(): " << std::endl;
        for (auto it = _top_neighbours_periodic.begin(); it != _top_neighbours_periodic.end();
             ++it) {
            std::cout << it->second << ", ";
        }
        std::cout << std::endl;
        for (auto it = _top_neighbours_periodic.begin(); it != _top_neighbours_periodic.end();
             ++it) {
            ids.push_back(it->first);
            halo_sizes.push_back(it->second);
        }
    }
}

void Partitioner::get_bottom_neighbours_periodic(
    std::vector<int>& ids, std::vector<int>& halo_sizes) const
{
    if (_p0) {
        // TODO: remove stdout
        std::cout << "_bottom_neighbour_ids_periodic.data(): " << std::endl;
        for (auto it = _bottom_neighbours_periodic.begin(); it != _bottom_neighbours_periodic.end();
             ++it) {
            std::cout << it->first << ", ";
        }
        std::cout << std::endl;
        std::cout << "_bottom_neighbour_halo_sizes_periodic.data(): " << std::endl;
        for (auto it = _bottom_neighbours_periodic.begin(); it != _bottom_neighbours_periodic.end();
             ++it) {
            std::cout << it->second << ", ";
        }
        std::cout << std::endl;
        for (auto it = _bottom_neighbours_periodic.begin(); it != _bottom_neighbours_periodic.end();
             ++it) {
            ids.push_back(it->first);
            halo_sizes.push_back(it->second);
        }
    }
}

void Partitioner::get_left_neighbours_periodic(
    std::vector<int>& ids, std::vector<int>& halo_sizes) const
{
    if (_p1) {
        // TODO: remove stdout
        std::cout << "_left_neighbour_ids_periodic.data(): " << std::endl;
        for (auto it = _left_neighbours_periodic.begin(); it != _left_neighbours_periodic.end();
             ++it) {
            std::cout << it->first << ", ";
        }
        std::cout << std::endl;
        std::cout << "_left_neighbour_halo_sizes_periodic.data(): " << std::endl;
        for (auto it = _left_neighbours_periodic.begin(); it != _left_neighbours_periodic.end();
             ++it) {
            std::cout << it->second << ", ";
        }
        std::cout << std::endl;
        for (auto it = _left_neighbours_periodic.begin(); it != _left_neighbours_periodic.end();
             ++it) {
            ids.push_back(it->first);
            halo_sizes.push_back(it->second);
        }
    }
}

void Partitioner::get_right_neighbours_periodic(
    std::vector<int>& ids, std::vector<int>& halo_sizes) const
{
    if (_p1) {
        // TODO: remove stdout
        std::cout << "_right_neighbour_ids_periodic.data(): " << std::endl;
        for (auto it = _right_neighbours_periodic.begin(); it != _right_neighbours_periodic.end();
             ++it) {
            std::cout << it->first << ", ";
        }
        std::cout << std::endl;
        std::cout << "_right_neighbour_halo_sizes_periodic.data(): " << std::endl;
        for (auto it = _right_neighbours_periodic.begin(); it != _right_neighbours_periodic.end();
             ++it) {
            std::cout << it->second << ", ";
        }
        std::cout << std::endl;
        for (auto it = _right_neighbours_periodic.begin(); it != _right_neighbours_periodic.end();
             ++it) {
            ids.push_back(it->first);
            halo_sizes.push_back(it->second);
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
    start[0] = _global_0;
    start[1] = _global_1;
    count[0] = _local_ext_0;
    count[1] = _local_ext_1;

    // Store data
    NC_CHECK(nc_var_par_access(nc_id, mask_nc_id, NC_COLLECTIVE));
    NC_CHECK(nc_put_vara_int(nc_id, mask_nc_id, start, count, _proc_id.data()));
    NC_CHECK(nc_close(nc_id));
}

void Partitioner::discover_periodic_neighbours()
{
    // FIXME: Dimensions are transposed from what I expected

    // Gather starts and counts for each dimension
    std::vector<int> start0(_num_procs), start1(_num_procs);
    std::vector<int> count0(_num_procs), count1(_num_procs);
    CHECK_MPI(MPI_Allgather(&_global_0, 1, MPI_INT, start0.data(), 1, MPI_INT, _comm));
    CHECK_MPI(MPI_Allgather(&_global_1, 1, MPI_INT, start1.data(), 1, MPI_INT, _comm));
    CHECK_MPI(MPI_Allgather(&_local_ext_0, 1, MPI_INT, count0.data(), 1, MPI_INT, _comm));
    CHECK_MPI(MPI_Allgather(&_local_ext_1, 1, MPI_INT, count1.data(), 1, MPI_INT, _comm));

    // Reshape _proc_id into 2D array
    std::vector<std::vector<int>> proc_id(_global_ext_0);
    for (int p = 0; p < _num_procs; p++) {
        for (int i = start0[p], k = 0; k < count0[p]; i++, k++) {
            proc_id[i] = std::vector<int>(_global_ext_1);
            for (int j = start1[p], l = 0; l < count1[p]; j++, l++) {
                proc_id[i][j] = p;
            }
        }
    }
    // TODO: remove stdout
    if (!_rank) {
        std::cout << "proc_id.data(): " << std::endl;
        for (int i = 0; i < _global_ext_0; i++) {
            for (int j = 0; j < _global_ext_1; j++) {
                std::cout << proc_id[i][j] << ", ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    };

    // Determine top and bottom periodic neighbours
    if (_p0) {
        std::vector<int> halo_sizes_t(_num_procs), halo_sizes_b(_num_procs);
        for (int i = 0; i < _global_ext_0; i++) {
            halo_sizes_t[proc_id[i][_global_ext_1 - 1]]++;
            halo_sizes_b[proc_id[i][0]]++;
        }
        for (int p = 0; p < _num_procs; p++) {
            if (halo_sizes_t[p]) {
                _top_neighbours_periodic.insert(std::pair<int, int>(p, halo_sizes_t[p]));
            }
            if (halo_sizes_b[p]) {
                _bottom_neighbours_periodic.insert(std::pair<int, int>(p, halo_sizes_b[p]));
            }
        }
    }
    // Determine left and right periodic neighbours
    if (_p1) {
        std::vector<int> halo_sizes_l(_num_procs), halo_sizes_r(_num_procs);
        for (int j = 0; j < _global_ext_1; j++) {
            halo_sizes_l[proc_id[_global_ext_0 - 1][j]]++;
            halo_sizes_r[proc_id[0][j]]++;
        }
        for (int p = 0; p < _num_procs; p++) {
            if (halo_sizes_l[p]) {
                _left_neighbours_periodic.insert(std::pair<int, int>(p, halo_sizes_l[p]));
            }
            if (halo_sizes_r[p]) {
                _right_neighbours_periodic.insert(std::pair<int, int>(p, halo_sizes_r[p]));
            }
        }
    }
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
    get_top_neighbours(top_ids, top_halos);
    get_bottom_neighbours(bottom_ids, bottom_halos);
    get_left_neighbours(left_ids, left_halos);
    get_right_neighbours(right_ids, right_halos);
    int top_num_neighbours = top_ids.size();
    int bottom_num_neighbours = bottom_ids.size();
    int left_num_neighbours = left_ids.size();
    int right_num_neighbours = right_ids.size();

    // Prepare periodic neighbour data
    std::vector<int> top_ids_p, bottom_ids_p, left_ids_p, right_ids_p;
    std::vector<int> top_halos_p, bottom_halos_p, left_halos_p, right_halos_p;
    get_top_neighbours_periodic(top_ids_p, top_halos_p);
    get_bottom_neighbours_periodic(bottom_ids_p, bottom_halos_p);
    get_left_neighbours_periodic(left_ids_p, left_halos_p);
    get_right_neighbours_periodic(right_ids_p, right_halos_p);
    int top_num_neighbours_p = top_ids_p.size();
    int bottom_num_neighbours_p = bottom_ids_p.size();
    int left_num_neighbours_p = left_ids_p.size();
    int right_num_neighbours_p = right_ids_p.size();

    // Compute global dimensions
    int top_dim, bottom_dim, left_dim, right_dim;
    CHECK_MPI(MPI_Allreduce(&top_num_neighbours, &top_dim, 1, MPI_INT, MPI_SUM, _comm));
    CHECK_MPI(MPI_Allreduce(&bottom_num_neighbours, &bottom_dim, 1, MPI_INT, MPI_SUM, _comm));
    CHECK_MPI(MPI_Allreduce(&left_num_neighbours, &left_dim, 1, MPI_INT, MPI_SUM, _comm));
    CHECK_MPI(MPI_Allreduce(&right_num_neighbours, &right_dim, 1, MPI_INT, MPI_SUM, _comm));

    // Compute global dimensions for periodic case
    int top_dim_p, bottom_dim_p, left_dim_p, right_dim_p;
    CHECK_MPI(MPI_Allreduce(&top_num_neighbours_p, &top_dim_p, 1, MPI_INT, MPI_SUM, _comm));
    CHECK_MPI(MPI_Allreduce(&bottom_num_neighbours_p, &bottom_dim_p, 1, MPI_INT, MPI_SUM, _comm));
    CHECK_MPI(MPI_Allreduce(&left_num_neighbours_p, &left_dim_p, 1, MPI_INT, MPI_SUM, _comm));
    CHECK_MPI(MPI_Allreduce(&right_num_neighbours_p, &right_dim_p, 1, MPI_INT, MPI_SUM, _comm));

    // Compute global offsets
    int top_offset = 0, bottom_offset = 0, left_offset = 0, right_offset = 0;
    CHECK_MPI(MPI_Exscan(&top_num_neighbours, &top_offset, 1, MPI_INT, MPI_SUM, _comm));
    CHECK_MPI(MPI_Exscan(&bottom_num_neighbours, &bottom_offset, 1, MPI_INT, MPI_SUM, _comm));
    CHECK_MPI(MPI_Exscan(&left_num_neighbours, &left_offset, 1, MPI_INT, MPI_SUM, _comm));
    CHECK_MPI(MPI_Exscan(&right_num_neighbours, &right_offset, 1, MPI_INT, MPI_SUM, _comm));

    // Compute global offsets for periodic case
    int top_offset_p = 0, bottom_offset_p = 0, left_offset_p = 0, right_offset_p = 0;
    CHECK_MPI(MPI_Exscan(&top_num_neighbours_p, &top_offset_p, 1, MPI_INT, MPI_SUM, _comm));
    CHECK_MPI(MPI_Exscan(&bottom_num_neighbours_p, &bottom_offset_p, 1, MPI_INT, MPI_SUM, _comm));
    CHECK_MPI(MPI_Exscan(&left_num_neighbours_p, &left_offset_p, 1, MPI_INT, MPI_SUM, _comm));
    CHECK_MPI(MPI_Exscan(&right_num_neighbours_p, &right_offset_p, 1, MPI_INT, MPI_SUM, _comm));

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
    NC_CHECK(nc_def_dim(nc_id, "T", top_dim, &top_dimid));
    NC_CHECK(nc_def_dim(nc_id, "B", bottom_dim, &bottom_dimid));
    NC_CHECK(nc_def_dim(nc_id, "L", left_dim, &left_dimid));
    NC_CHECK(nc_def_dim(nc_id, "R", right_dim, &right_dimid));

    // Define periodic dimensions in netCDF file
    int top_dimid_p, bottom_dimid_p, left_dimid_p, right_dimid_p;
    NC_CHECK(nc_def_dim(nc_id, "T_periodic", top_dim_p, &top_dimid_p));
    NC_CHECK(nc_def_dim(nc_id, "B_periodic", bottom_dim_p, &bottom_dimid_p));
    NC_CHECK(nc_def_dim(nc_id, "L_periodic", left_dim_p, &left_dimid_p));
    NC_CHECK(nc_def_dim(nc_id, "R_periodic", right_dim_p, &right_dimid_p));

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
    int top_num_vid_p, bottom_num_vid_p, left_num_vid_p, right_num_vid_p;
    int top_ids_vid_p, bottom_ids_vid_p, left_ids_vid_p, right_ids_vid_p;
    int top_halos_vid_p, bottom_halos_vid_p, left_halos_vid_p, right_halos_vid_p;
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
    // Periodic members of connectivity group
    NC_CHECK(
        nc_def_var(connectivity_gid, "top_neighbours_periodic", NC_INT, 1, &dimid, &top_num_vid_p));
    NC_CHECK(nc_def_var(
        connectivity_gid, "top_neighbour_ids_periodic", NC_INT, 1, &top_dimid_p, &top_ids_vid_p));
    NC_CHECK(nc_def_var(connectivity_gid, "top_neighbour_halos_periodic", NC_INT, 1, &top_dimid_p,
        &top_halos_vid_p));
    NC_CHECK(nc_def_var(
        connectivity_gid, "bottom_neighbours_periodic", NC_INT, 1, &dimid, &bottom_num_vid_p));
    NC_CHECK(nc_def_var(connectivity_gid, "bottom_neighbour_ids_periodic", NC_INT, 1,
        &bottom_dimid_p, &bottom_ids_vid_p));
    NC_CHECK(nc_def_var(connectivity_gid, "bottom_neighbour_halos_periodic", NC_INT, 1,
        &bottom_dimid_p, &bottom_halos_vid_p));
    NC_CHECK(nc_def_var(
        connectivity_gid, "left_neighbours_periodic", NC_INT, 1, &dimid, &left_num_vid_p));
    NC_CHECK(nc_def_var(connectivity_gid, "left_neighbour_ids_periodic", NC_INT, 1, &left_dimid_p,
        &left_ids_vid_p));
    NC_CHECK(nc_def_var(connectivity_gid, "left_neighbour_halos_periodic", NC_INT, 1, &left_dimid_p,
        &left_halos_vid_p));
    NC_CHECK(nc_def_var(
        connectivity_gid, "right_neighbours_periodic", NC_INT, 1, &dimid, &right_num_vid_p));
    NC_CHECK(nc_def_var(connectivity_gid, "right_neighbour_ids_periodic", NC_INT, 1, &right_dimid_p,
        &right_ids_vid_p));
    NC_CHECK(nc_def_var(connectivity_gid, "right_neighbour_halos_periodic", NC_INT, 1,
        &right_dimid_p, &right_halos_vid_p));

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
    NC_CHECK(nc_put_var1_int(connectivity_gid, top_num_vid, &start, &top_num_neighbours));
    NC_CHECK(nc_var_par_access(connectivity_gid, bottom_num_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_var1_int(connectivity_gid, bottom_num_vid, &start, &bottom_num_neighbours));
    NC_CHECK(nc_var_par_access(connectivity_gid, left_num_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_var1_int(connectivity_gid, left_num_vid, &start, &left_num_neighbours));
    NC_CHECK(nc_var_par_access(connectivity_gid, right_num_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_var1_int(connectivity_gid, right_num_vid, &start, &right_num_neighbours));
    NC_CHECK(nc_var_par_access(connectivity_gid, top_num_vid_p, NC_COLLECTIVE));
    NC_CHECK(nc_put_var1_int(connectivity_gid, top_num_vid_p, &start, &top_num_neighbours_p));
    NC_CHECK(nc_var_par_access(connectivity_gid, bottom_num_vid_p, NC_COLLECTIVE));
    NC_CHECK(nc_put_var1_int(connectivity_gid, bottom_num_vid_p, &start, &bottom_num_neighbours_p));
    NC_CHECK(nc_var_par_access(connectivity_gid, left_num_vid_p, NC_COLLECTIVE));
    NC_CHECK(nc_put_var1_int(connectivity_gid, left_num_vid_p, &start, &left_num_neighbours_p));
    NC_CHECK(nc_var_par_access(connectivity_gid, right_num_vid_p, NC_COLLECTIVE));
    NC_CHECK(nc_put_var1_int(connectivity_gid, right_num_vid_p, &start, &right_num_neighbours_p));

    start = top_offset;
    count = top_num_neighbours;
    NC_CHECK(nc_var_par_access(connectivity_gid, top_ids_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_vara_int(connectivity_gid, top_ids_vid, &start, &count, top_ids.data()));
    NC_CHECK(nc_var_par_access(connectivity_gid, top_halos_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_vara_int(connectivity_gid, top_halos_vid, &start, &count, top_halos.data()));
    start = bottom_offset;
    count = bottom_num_neighbours;
    NC_CHECK(nc_var_par_access(connectivity_gid, bottom_ids_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_vara_int(connectivity_gid, bottom_ids_vid, &start, &count, bottom_ids.data()));
    NC_CHECK(nc_var_par_access(connectivity_gid, bottom_halos_vid, NC_COLLECTIVE));
    NC_CHECK(
        nc_put_vara_int(connectivity_gid, bottom_halos_vid, &start, &count, bottom_halos.data()));
    start = left_offset;
    count = left_num_neighbours;
    NC_CHECK(nc_var_par_access(connectivity_gid, left_ids_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_vara_int(connectivity_gid, left_ids_vid, &start, &count, left_ids.data()));
    NC_CHECK(nc_var_par_access(connectivity_gid, left_halos_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_vara_int(connectivity_gid, left_halos_vid, &start, &count, left_halos.data()));
    start = right_offset;
    count = right_num_neighbours;
    NC_CHECK(nc_var_par_access(connectivity_gid, right_ids_vid, NC_COLLECTIVE));
    NC_CHECK(nc_put_vara_int(connectivity_gid, right_ids_vid, &start, &count, right_ids.data()));
    NC_CHECK(nc_var_par_access(connectivity_gid, right_halos_vid, NC_COLLECTIVE));
    NC_CHECK(
        nc_put_vara_int(connectivity_gid, right_halos_vid, &start, &count, right_halos.data()));

    start = top_offset_p;
    count = top_num_neighbours_p;
    NC_CHECK(nc_var_par_access(connectivity_gid, top_ids_vid_p, NC_COLLECTIVE));
    NC_CHECK(nc_put_vara_int(connectivity_gid, top_ids_vid_p, &start, &count, top_ids_p.data()));
    NC_CHECK(nc_var_par_access(connectivity_gid, top_halos_vid_p, NC_COLLECTIVE));
    NC_CHECK(
        nc_put_vara_int(connectivity_gid, top_halos_vid_p, &start, &count, top_halos_p.data()));
    start = bottom_offset_p;
    count = bottom_num_neighbours_p;
    NC_CHECK(nc_var_par_access(connectivity_gid, bottom_ids_vid_p, NC_COLLECTIVE));
    NC_CHECK(
        nc_put_vara_int(connectivity_gid, bottom_ids_vid_p, &start, &count, bottom_ids_p.data()));
    NC_CHECK(nc_var_par_access(connectivity_gid, bottom_halos_vid_p, NC_COLLECTIVE));
    NC_CHECK(nc_put_vara_int(
        connectivity_gid, bottom_halos_vid_p, &start, &count, bottom_halos_p.data()));
    start = left_offset_p;
    count = left_num_neighbours_p;
    NC_CHECK(nc_var_par_access(connectivity_gid, left_ids_vid_p, NC_COLLECTIVE));
    NC_CHECK(nc_put_vara_int(connectivity_gid, left_ids_vid_p, &start, &count, left_ids_p.data()));
    NC_CHECK(nc_var_par_access(connectivity_gid, left_halos_vid_p, NC_COLLECTIVE));
    NC_CHECK(
        nc_put_vara_int(connectivity_gid, left_halos_vid_p, &start, &count, left_halos_p.data()));
    start = right_offset_p;
    count = right_num_neighbours_p;
    NC_CHECK(nc_var_par_access(connectivity_gid, right_ids_vid_p, NC_COLLECTIVE));
    NC_CHECK(
        nc_put_vara_int(connectivity_gid, right_ids_vid_p, &start, &count, right_ids_p.data()));
    NC_CHECK(nc_var_par_access(connectivity_gid, right_halos_vid_p, NC_COLLECTIVE));
    NC_CHECK(
        nc_put_vara_int(connectivity_gid, right_halos_vid_p, &start, &count, right_halos_p.data()));

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
                _top_neighbours.insert(std::pair<int, int>(i, halo_size));
            }
            if (top_right_1[_rank] >= bottom_left_1[i] && top_right_1[_rank] <= bottom_right_1[i]
                && bottom_left_1[i] >= top_left_1[_rank]
                && (top_right_0[_rank] - bottom_right_0[i] == 1)) {
                int halo_size = top_right_1[_rank] - bottom_left_1[i] + 1;
                _top_neighbours.insert(std::pair<int, int>(i, halo_size));
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
                _bottom_neighbours.insert(std::pair<int, int>(i, halo_size));
            }
            if (bottom_right_1[_rank] >= top_left_1[i] && bottom_right_1[_rank] <= top_right_1[i]
                && top_left_1[i] >= bottom_left_1[_rank]
                && (top_right_0[i] - bottom_right_0[_rank] == 1)) {
                int halo_size = bottom_right_1[_rank] - top_left_1[i] + 1;
                _bottom_neighbours.insert(std::pair<int, int>(i, halo_size));
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
                _left_neighbours.insert(std::pair<int, int>(i, halo_size));
            }
            if (bottom_left_0[_rank] >= top_right_0[i] && bottom_left_0[_rank] <= bottom_right_0[i]
                && top_left_0[_rank] <= top_right_0[i]
                && (bottom_left_1[_rank] - top_right_1[i] == 1)) {
                int halo_size = bottom_left_0[_rank] - top_right_0[i] + 1;
                _left_neighbours.insert(std::pair<int, int>(i, halo_size));
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
                _right_neighbours.insert(std::pair<int, int>(i, halo_size));
            }
            if (bottom_right_0[_rank] >= top_left_0[i] && bottom_right_0[_rank] <= bottom_left_0[i]
                && top_right_0[_rank] <= top_left_0[i]
                && (top_left_1[i] - top_right_1[_rank] == 1)) {
                int halo_size = bottom_right_0[_rank] - top_left_0[i] + 1;
                _right_neighbours.insert(std::pair<int, int>(i, halo_size));
            }
        }
    }
}
