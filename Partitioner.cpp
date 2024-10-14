/*!
 * @file Partitioner.cpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @date 12 Sep 2024
 */

#include "DomainUtils.hpp"
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

void Partitioner::get_neighbours(
    std::vector<std::vector<int>>& ids, std::vector<std::vector<int>>& halo_sizes) const
{
    for (int idx = 0; idx < NNBRS; idx++) {
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
    NC_CHECK(nc_put_att_int(nc_id, NC_GLOBAL, "num_processes", NC_INT, 1, &_total_num_procs));

    // Create 2 dimensions
    // The values to be written are associated with the netCDF variable by
    // assuming that the last dimension of the netCDF variable varies fastest in
    // the C interface
    const int NDIMS = 2; // TODO: Why redeclared?
    int dimid[NDIMS];
    std::vector<std::string> dim_chars = { "x", "y" };
    for (int idx = 0; idx < NDIMS; idx++) {
        NC_CHECK(nc_def_dim(nc_id, dim_chars[idx].c_str(), _global_ext[idx], &dimid[idx]));
    }

    // Create variables
    int mask_nc_id;
    NC_CHECK(nc_def_var(nc_id, "pid", NC_INT, NDIMS, dimid, &mask_nc_id));

    // Write metadata to file
    NC_CHECK(nc_enddef(nc_id));

    // Set up slab for this process
    size_t start[NDIMS], count[NDIMS];
    for (int idx = 0; idx < NDIMS; idx++) {
        start[idx] = _global[idx];
        count[idx] = _local_ext[idx];
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

    // There are two neighbours for each dimension
    const int NNBRS = NDIMS * 2; // TODO: Why redeclared?

    // Prepare neighbour data
    std::vector<std::vector<int>> ids(NNBRS), halos(NNBRS);
    get_neighbours(ids, halos);
    std::vector<int> num_neighbours(NNBRS), dims(NNBRS, 0), offsets(NNBRS, 0);
    for (int idx = 0; idx < NNBRS; idx++) {
        num_neighbours[idx] = (int)ids[idx].size();
        CHECK_MPI(MPI_Allreduce(&num_neighbours[idx], &dims[idx], 1, MPI_INT, MPI_SUM, _comm));
        CHECK_MPI(MPI_Exscan(&num_neighbours[idx], &offsets[idx], 1, MPI_INT, MPI_SUM, _comm));
    }

    // Define dimensions in netCDF file
    int dimid;
    std::vector<int> dimids(NNBRS);
    NC_CHECK(nc_def_dim(nc_id, "P", _total_num_procs, &dimid));
    for (int idx = 0; idx < NNBRS; idx++) {
        NC_CHECK(nc_def_dim(nc_id, dir_chars[idx].c_str(), dims[idx], &dimids[idx]));
    }

    // Define groups in netCDF file
    int bbox_gid, connectivity_gid;
    NC_CHECK(nc_def_grp(nc_id, "bounding_boxes", &bbox_gid));
    NC_CHECK(nc_def_grp(nc_id, "connectivity", &connectivity_gid));

    // Define variables in netCDF file
    int top_vid[NDIMS];
    int cnt_vid[NDIMS];
    int num_vid[NNBRS];
    int ids_vid[NNBRS];
    int halos_vid[NNBRS];
    for (int idx = 0; idx < NDIMS; idx++) {
        // Bounding boxes group
        NC_CHECK(nc_def_var(
            bbox_gid, ("domain_" + dim_chars[idx]).c_str(), NC_INT, 1, &dimid, &top_vid[idx]));
        NC_CHECK(nc_def_var(bbox_gid, ("domain_extent_" + dim_chars[idx]).c_str(), NC_INT, 1,
            &dimid, &cnt_vid[idx]));
    }
    for (int idx = 0; idx < NNBRS; idx++) {
        // Connectivity group
        NC_CHECK(nc_def_var(connectivity_gid, (dir_names[idx] + "_neighbours").c_str(), NC_INT, 1,
            &dimid, &num_vid[idx]));
        NC_CHECK(nc_def_var(connectivity_gid, (dir_names[idx] + "_neighbour_ids").c_str(), NC_INT,
            1, &dimids[idx], &ids_vid[idx]));
        NC_CHECK(nc_def_var(connectivity_gid, (dir_names[idx] + "_neighbour_halos").c_str(), NC_INT,
            1, &dimids[idx], &halos_vid[idx]));
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
    for (int idx = 0; idx < NNBRS; idx++) {
        size_t start = _rank;
        NC_CHECK(nc_var_par_access(connectivity_gid, num_vid[idx], NC_COLLECTIVE));
        NC_CHECK(nc_put_var1_int(connectivity_gid, num_vid[idx], &start, &num_neighbours[idx]));
        start = offsets[idx];
        size_t count = num_neighbours[idx];
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
        if (p != _rank) {

            // Find my left neighbours
            // i.e., the left edge of domain[_rank] has to match the right edge of domain[p]
            if (domains[_rank].p1.x == domains[p].p2.x) {
                // next compute overlap (if overlap is zero then they are not neighbours)
                // domain overlap is equivalent to the amount of data transfered in halo exchnage.
                int halo_size = domainOverlap(domains[_rank], domains[p], 'y');
                if (halo_size) {
                    _neighbours[0].insert(std::pair<int, int>(p, halo_size));
                }
            }

            // Find my right neighbours
            if (domains[_rank].p2.x == domains[p].p1.x) {
                int halo_size = domainOverlap(domains[_rank], domains[p], 'y');
                if (halo_size) {
                    _neighbours[1].insert(std::pair<int, int>(p, halo_size));
                }
            }

            // Find my top neighbours
            if (domains[_rank].p2.y == domains[p].p1.y) {
                int halo_size = domainOverlap(domains[_rank], domains[p], 'x');
                if (halo_size) {
                    _neighbours[2].insert(std::pair<int, int>(p, halo_size));
                }
            }

            // Find my bottom neighbours
            if (domains[_rank].p1.y == domains[p].p2.y) {
                int halo_size = domainOverlap(domains[_rank], domains[p], 'x');
                if (halo_size) {
                    _neighbours[3].insert(std::pair<int, int>(p, halo_size));
                }
            }
        }
    }
}
