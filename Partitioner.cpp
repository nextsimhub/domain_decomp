/*!
 * @file Partitioner.cpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @date 05 Nov 2024
 */

#include "Partitioner.hpp"
#include "DomainUtils.hpp"
#include "Utils.hpp"
#include "ZoltanPartitioner.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include <netcdf.h>
#include <netcdf_par.h>

bool Partitioner::is_neighbour(
    const Domain d1, const Domain d2, const Edge edge, const bool is_px, const bool is_py)
{
    if (edge == TOP) {
        // Check if TOP neighbour i.e., the bottom of domain d2 must match the top of domain d1.
        // The logic for the other edges is essentially the same.
        if (is_py) {
            if (d1.p2.y == d2.p1.y + _global_ext[1]) {
                return true;
            }
        } else {
            if (d1.p2.y == d2.p1.y) {
                return true;
            }
        }
    } else if (edge == BOTTOM) {
        if (is_py) {
            if (d1.p1.y == d2.p2.y - _global_ext[1]) {
                return true;
            }
        } else {
            if (d1.p1.y == d2.p2.y) {
                return true;
            }
        }
    } else if (edge == LEFT) {
        if (is_px) {
            if (d1.p1.x == d2.p2.x - _global_ext[0]) {
                return true;
            }
        } else {
            if (d1.p1.x == d2.p2.x) {
                return true;
            }
        }
    } else if (edge == RIGHT) {
        if (is_px) {
            if (d1.p2.x == d2.p1.x + _global_ext[0]) {
                return true;
            }
        } else {

            if (d1.p2.x == d2.p1.x) {
                return true;
            }
        }
    } else {
        std::cerr << "ERROR: edge must be LEFT, RIGHT, BOTTOM, TOP." << std::endl;
        exit(EXIT_FAILURE);
    }
    return false;
}

int Partitioner::halo_start(const Domain d1, const Domain d2, const Edge edge)
{
    int start = 0;
    if (edge == TOP) {
        // offset between domains
        int dx = std::max(d1.p1.x, d2.p1.x) - d2.p1.x;
        start = dx;
    } else if (edge == BOTTOM) {
        int dx = std::max(d1.p1.x, d2.p1.x) - d2.p1.x;
        start = (d2.get_height() - 1) * d2.get_width() + dx;
    } else if (edge == LEFT) {
        int dy = std::max(d1.p1.y, d2.p1.y) - d2.p1.y;
        start = ((dy + 1) * d2.get_width()) - 1;
    } else if (edge == RIGHT) {
        int dy = std::max(d1.p1.y, d2.p1.y) - d2.p1.y;
        start = dy * d2.get_width();
    } else {
        std::cerr << "ERROR: edge must be LEFT, RIGHT, BOTTOM, TOP." << std::endl;
        exit(EXIT_FAILURE);
    }
    return start;
}

Partitioner::Partitioner(MPI_Comm comm)
{
    _comm = comm;
    CHECK_MPI(MPI_Comm_size(comm, &_total_num_procs));
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

void Partitioner::get_neighbour_info(std::vector<std::vector<int>>& ids,
    std::vector<std::vector<int>>& halo_sizes, std::vector<std::vector<int>>& halo_starts) const
{
    for (auto edge : edges) {
        for (auto it = _neighbours[edge].begin(); it != _neighbours[edge].end(); ++it) {
            ids[edge].push_back(it->first);
            halo_sizes[edge].push_back(it->second);
        }
        for (auto it = _halo_starts[edge].begin(); it != _halo_starts[edge].end(); ++it) {
            halo_starts[edge].push_back(it->second);
        }
    }
}

void Partitioner::get_neighbour_info_periodic(std::vector<std::vector<int>>& ids,
    std::vector<std::vector<int>>& halo_sizes, std::vector<std::vector<int>>& halo_starts) const
{
    for (auto edge : edges) {
        if (((edge == LEFT || edge == RIGHT) && _px) || ((edge == TOP || edge == BOTTOM) && _py)) {
            for (auto it = _neighbours_p[edge].begin(); it != _neighbours_p[edge].end(); ++it) {
                ids[edge].push_back(it->first);
                halo_sizes[edge].push_back(it->second);
            }
            for (auto it = _halo_starts_p[edge].begin(); it != _halo_starts_p[edge].end(); ++it) {
                halo_starts[edge].push_back(it->second);
            }
        }
    }
}

void Partitioner::save_mask(const std::string& filename) const
{
    // Use C API for parallel I/O
    int nc_id, nc_mode;
    nc_mode = NC_CLOBBER | NC_NETCDF4;
    NC_CHECK(nc_create_par(filename.c_str(), nc_mode, _comm, MPI_INFO_NULL, &nc_id));
    NC_CHECK(nc_put_att_int(nc_id, NC_GLOBAL, "num_processes", NC_INT, 1, &_total_num_procs));

    // Create 2 dimensions
    // The values to be written are associated with the netCDF variable by
    // assuming that the last dimension of the netCDF variable varies fastest in
    // the C interface
    int dimid[NDIMS];
    // for nextsimdg code we always want the output in the order of yx
    std::vector<std::string> dim_chars = { "y", "x" };
    for (int idx = 0; idx < NDIMS; idx++) {
        NC_CHECK(
            nc_def_dim(nc_id, dim_chars[idx].c_str(), _global_ext[NDIMS - 1 - idx], &dimid[idx]));
    }

    // Create variables
    int mask_nc_id;
    NC_CHECK(nc_def_var(nc_id, "pid", NC_INT, NDIMS, dimid, &mask_nc_id));

    // Write metadata to file
    NC_CHECK(nc_enddef(nc_id));

    // Set up slab for this process
    size_t start[NDIMS], count[NDIMS];
    for (int idx = 0; idx < NDIMS; idx++) {
        start[idx] = _global[NDIMS - 1 - idx];
        count[idx] = _local_ext[NDIMS - 1 - idx];
    }

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

    // Create 2 dimensions
    // The values to be written are associated with the netCDF variable by
    // assuming that the last dimension of the netCDF variable varies fastest in
    // the C interface
    const int NDIMS = 2; // TODO: Why redeclared?
    int dimid_global[NDIMS];
    for (int idx = 0; idx < NDIMS; idx++) {
        NC_CHECK(nc_def_dim(
            nc_id, global_extent_names[idx].c_str(), _global_ext[idx], &dimid_global[idx]));
    }

    // Prepare neighbour data
    std::vector<std::vector<int>> ids(N_EDGE), halos(N_EDGE), halo_starts(N_EDGE);
    get_neighbour_info(ids, halos, halo_starts);
    std::vector<int> num_neighbours(N_EDGE), dims(N_EDGE, 0), offsets(N_EDGE, 0);
    for (auto edge : edges) {
        num_neighbours[edge] = (int)ids[edge].size();
        CHECK_MPI(MPI_Allreduce(&num_neighbours[edge], &dims[edge], 1, MPI_INT, MPI_SUM, _comm));
        CHECK_MPI(MPI_Exscan(&num_neighbours[edge], &offsets[edge], 1, MPI_INT, MPI_SUM, _comm));
    }

    // Prepare periodic neighbour data
    std::vector<std::vector<int>> ids_p(N_EDGE), halos_p(N_EDGE), halo_starts_p(N_EDGE);
    get_neighbour_info_periodic(ids_p, halos_p, halo_starts_p);
    std::vector<int> num_neighbours_p(N_EDGE), dims_p(N_EDGE, 0), offsets_p(N_EDGE, 0);
    for (auto edge : edges) {
        num_neighbours_p[edge] = (int)ids_p[edge].size();
        CHECK_MPI(
            MPI_Allreduce(&num_neighbours_p[edge], &dims_p[edge], 1, MPI_INT, MPI_SUM, _comm));
        CHECK_MPI(
            MPI_Exscan(&num_neighbours_p[edge], &offsets_p[edge], 1, MPI_INT, MPI_SUM, _comm));
    }

    // Define dimensions in netCDF file
    int dimid;
    std::vector<int> dimids(N_EDGE);
    NC_CHECK(nc_def_dim(nc_id, "P", _total_num_procs, &dimid));
    for (auto edge : edges) {
        NC_CHECK(nc_def_dim(nc_id, dir_chars[edge].c_str(), dims[edge], &dimids[edge]));
    }

    // Define periodic dimensions in netCDF file
    std::vector<int> dimids_p(N_EDGE);
    for (auto edge : edges) {
        NC_CHECK(nc_def_dim(
            nc_id, (dir_chars[edge] + "_periodic").c_str(), dims_p[edge], &dimids_p[edge]));
    }

    // Define groups in netCDF file
    int bbox_gid, connectivity_gid;
    NC_CHECK(nc_def_grp(nc_id, "bounding_boxes", &bbox_gid));
    NC_CHECK(nc_def_grp(nc_id, "connectivity", &connectivity_gid));

    // Define variables in netCDF file
    int top_vid[NDIMS];
    int cnt_vid[NDIMS];
    int num_vid[N_EDGE];
    for (int idx = 0; idx < NDIMS; idx++) {
        // Bounding boxes group
        NC_CHECK(nc_def_var(
            bbox_gid, ("domain_" + dim_chars[idx]).c_str(), NC_INT, 1, &dimid, &top_vid[idx]));
        NC_CHECK(nc_def_var(bbox_gid, ("domain_extent_" + dim_chars[idx]).c_str(), NC_INT, 1,
            &dimid, &cnt_vid[idx]));
    }

    int ids_vid[N_EDGE];
    int halos_vid[N_EDGE];
    int halo_starts_vid[N_EDGE];
    for (auto edge : edges) {
        // Connectivity group
        NC_CHECK(nc_def_var(connectivity_gid, (dir_names[edge] + "_neighbours").c_str(), NC_INT, 1,
            &dimid, &num_vid[edge]));
        NC_CHECK(nc_def_var(connectivity_gid, (dir_names[edge] + "_neighbour_ids").c_str(), NC_INT,
            1, &dimids[edge], &ids_vid[edge]));
        NC_CHECK(nc_def_var(connectivity_gid, (dir_names[edge] + "_neighbour_halos").c_str(),
            NC_INT, 1, &dimids[edge], &halos_vid[edge]));
        NC_CHECK(nc_def_var(connectivity_gid, (dir_names[edge] + "_neighbour_halo_starts").c_str(),
            NC_INT, 1, &dimids[edge], &halo_starts_vid[edge]));
    }

    int num_vid_p[N_EDGE];
    int ids_vid_p[N_EDGE];
    int halos_vid_p[N_EDGE];
    int halo_starts_vid_p[N_EDGE];
    for (auto edge : edges) {
        // Periodic members of connectivity group
        NC_CHECK(nc_def_var(connectivity_gid, (dir_names[edge] + "_neighbours_periodic").c_str(),
            NC_INT, 1, &dimid, &num_vid_p[edge]));
        NC_CHECK(nc_def_var(connectivity_gid, (dir_names[edge] + "_neighbour_ids_periodic").c_str(),
            NC_INT, 1, &dimids_p[edge], &ids_vid_p[edge]));
        NC_CHECK(
            nc_def_var(connectivity_gid, (dir_names[edge] + "_neighbour_halos_periodic").c_str(),
                NC_INT, 1, &dimids_p[edge], &halos_vid_p[edge]));
        NC_CHECK(nc_def_var(connectivity_gid,
            (dir_names[edge] + "_neighbour_halo_starts_periodic").c_str(), NC_INT, 1,
            &dimids_p[edge], &halo_starts_vid_p[edge]));
    }

    // Write metadata to file
    NC_CHECK(nc_enddef(nc_id));

    // Store data
    for (int idx = 0; idx < NDIMS; idx++) {
        size_t start = _rank;
        NC_CHECK(nc_var_par_access(bbox_gid, top_vid[idx], NC_COLLECTIVE));
        NC_CHECK(nc_put_var1_int(bbox_gid, top_vid[idx], &start, &_global_new[idx]));
        NC_CHECK(nc_var_par_access(bbox_gid, cnt_vid[idx], NC_COLLECTIVE));
        NC_CHECK(nc_put_var1_int(bbox_gid, cnt_vid[idx], &start, &_local_ext_new[idx]));
    }
    for (auto edge : edges) {
        // Numbers of neighbours
        size_t start = _rank;
        NC_CHECK(nc_var_par_access(connectivity_gid, num_vid[edge], NC_COLLECTIVE));
        NC_CHECK(nc_put_var1_int(connectivity_gid, num_vid[edge], &start, &num_neighbours[edge]));
        // Numbers of neighbours for periodic dimensions
        NC_CHECK(nc_var_par_access(connectivity_gid, num_vid_p[edge], NC_COLLECTIVE));
        NC_CHECK(
            nc_put_var1_int(connectivity_gid, num_vid_p[edge], &start, &num_neighbours_p[edge]));
        // IDs and halos
        start = offsets[edge];
        size_t count = num_neighbours[edge];
        NC_CHECK(nc_var_par_access(connectivity_gid, ids_vid[edge], NC_COLLECTIVE));
        NC_CHECK(
            nc_put_vara_int(connectivity_gid, ids_vid[edge], &start, &count, ids[edge].data()));
        NC_CHECK(nc_var_par_access(connectivity_gid, halos_vid[edge], NC_COLLECTIVE));
        NC_CHECK(
            nc_put_vara_int(connectivity_gid, halos_vid[edge], &start, &count, halos[edge].data()));
        NC_CHECK(nc_var_par_access(connectivity_gid, halo_starts_vid[edge], NC_COLLECTIVE));
        NC_CHECK(nc_put_vara_int(
            connectivity_gid, halo_starts_vid[edge], &start, &count, halo_starts[edge].data()));
        // IDs and halos for periodic dimensions
        start = offsets_p[edge];
        count = num_neighbours_p[edge];
        NC_CHECK(nc_var_par_access(connectivity_gid, ids_vid_p[edge], NC_COLLECTIVE));
        NC_CHECK(
            nc_put_vara_int(connectivity_gid, ids_vid_p[edge], &start, &count, ids_p[edge].data()));
        NC_CHECK(nc_var_par_access(connectivity_gid, halos_vid_p[edge], NC_COLLECTIVE));
        NC_CHECK(nc_put_vara_int(
            connectivity_gid, halos_vid_p[edge], &start, &count, halos_p[edge].data()));
        NC_CHECK(nc_put_vara_int(
            connectivity_gid, halo_starts_vid_p[edge], &start, &count, halo_starts_p[edge].data()));
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
    /*

       In the netcdf file data are stored in this order
     O┌───►X
      │                             20           30
      │  0 ┌─────────────────────────┬────────────┐
     Y▼    │                         │            │
           │                         │            │
           │            0            │            │
           │                         │            │
           │                         │            │
        12 ├─────────────────────────┤      2     │
           │                         │            │
           │                         │            │
           │            1            │            │
           │                         │            │
           │                         │            │
        24 └─────────────────────────┴────────────┘


        But in the calculation of TLBR neighbours we need to flip the Y-axis

                                    20           30
        24 ┌─────────────────────────┬────────────┐
           │                         │            │
           │                         │            │
           │            1            │            │
           │                         │            │
           │                         │            │
        12 ├─────────────────────────┤      2     │
           │                         │            │
           │                         │            │
           │            0            │            │
           │                         │            │
     Y▲    │                         │            │
      │  0 └─────────────────────────┴────────────┘
      │
     0└───►X
     */

    // Gather bounding boxes for all processes
    std::vector<Point> origins(_total_num_procs);
    std::vector<Point> extents(_total_num_procs);
    std::vector<Domain> domains(_total_num_procs);
    std::vector<int> tmp0(_total_num_procs);
    std::vector<int> tmp1(_total_num_procs);

    CHECK_MPI(MPI_Allgather(&_global_new[0], 1, MPI_INT, tmp0.data(), 1, MPI_INT, _comm));
    CHECK_MPI(MPI_Allgather(&_global_new[1], 1, MPI_INT, tmp1.data(), 1, MPI_INT, _comm));

    // origin points mark the bottom-left corner of each domain
    for (int p = 0; p < _total_num_procs; p++) {
        origins[p].x = tmp0[p];
        origins[p].y = tmp1[p];
    }

    CHECK_MPI(MPI_Allgather(&_local_ext_new[0], 1, MPI_INT, tmp0.data(), 1, MPI_INT, _comm));
    CHECK_MPI(MPI_Allgather(&_local_ext_new[1], 1, MPI_INT, tmp1.data(), 1, MPI_INT, _comm));

    // extents can be used to find the top-right corner of each domain
    for (int p = 0; p < _total_num_procs; p++) {
        extents[p].x = tmp0[p];
        extents[p].y = tmp1[p];
    }

    // generate domains
    for (int p = 0; p < _total_num_procs; p++) {
        domains[p].p1.x = origins[p].x;
        domains[p].p1.y = origins[p].y;
        domains[p].p2.x = origins[p].x + extents[p].x;
        domains[p].p2.y = origins[p].y + extents[p].y;
    }

    for (int p = 0; p < _total_num_procs; p++) {

        // When finding neighbours *within* the domain, we don't check against the current rank
        // because a subdomain can't be a neighbour of itself.
        if (p != _rank) {

            for (auto edge : edges) {
                if (is_neighbour(domains[_rank], domains[p], edge)) {
                    int halo_size = domain_overlap(domains[_rank], domains[p], edge);
                    if (halo_size > 0) {
                        _neighbours[edge].insert(std::pair<int, int>(p, halo_size));
                        int start = halo_start(domains[_rank], domains[p], edge);
                        _halo_starts[edge].insert(std::pair<int, int>(p, start));
                    }
                }
            }
        }

        // When finding neighours *across periodic boundaries*, we need to check against the
        // current rank, too, because a subdomain can be a periodic neighbour of itself.
        for (auto edge : edges) {
            if (is_neighbour(domains[_rank], domains[p], edge, _px, _py)) {
                int halo_size = domain_overlap(domains[_rank], domains[p], edge);
                if (halo_size > 0) {
                    _neighbours_p[edge].insert(std::pair<int, int>(p, halo_size));
                    int start = halo_start(domains[_rank], domains[p], edge);
                    _halo_starts_p[edge].insert(std::pair<int, int>(p, start));
                }
            }
        }
    }
}
