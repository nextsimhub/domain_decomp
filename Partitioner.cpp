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

/*!
 * @brief Compute modulo operation with non-negative result
 *
 * This function computes a % n (which is the modulo operator) and ensures the result is always
 * non-negative. This is particularly useful for periodic boundary conditions where we need to wrap
 * around coordinates that might be negative.
 *
 * @param a The dividend
 * @param n The divisor (modulus)
 * @return The non-negative remainder of a divided by n
 *
 * @note The standard % operator in C++ can return negative results when a is negative.
 *       This function adjusts the result to always be in the range [0, n-1].
 *       e.g.,
 *       mod(-1, 5);  // 4
 *       mod(6, 5);   // 1
 *       mod(4, 5);   // 4
 *       mod(-6, 5);  // 4
 */
int mod(int a, int n) { return ((a % n) + n) % n; }

bool Partitioner::isNeighbour(
    const Domain d1, const Domain d2, const Edge edge, const bool isPx, const bool isPy)
{
    // For TOP & BOTTOM Edges check that the domains share a y-coordinate AND that the horizontal
    // overlap is non-zero

    // For LEFT & RIGHT Edges check that the domains share a x-coordinate AND that the vertical
    // overlap is non-zero
    if (edge == TOP) {
        if (isPy) {
            return d1.p2.y == d2.p1.y + _globalExt[1] && d2.p1.x < d1.p2.x && d2.p2.x > d1.p1.x;
        } else {
            return d1.p2.y == d2.p1.y && d2.p1.x < d1.p2.x && d2.p2.x > d1.p1.x;
        }
    } else if (edge == BOTTOM) {
        if (isPy) {
            return d1.p1.y == d2.p2.y - _globalExt[1] && d2.p1.x < d1.p2.x && d2.p2.x > d1.p1.x;
        } else {
            return d1.p1.y == d2.p2.y && d2.p1.x < d1.p2.x && d2.p2.x > d1.p1.x;
        }
    } else if (edge == LEFT) {
        if (isPx) {
            return d1.p1.x == d2.p2.x - _globalExt[0] && d2.p1.y < d1.p2.y && d2.p2.y > d1.p1.y;
        } else {
            return d1.p1.x == d2.p2.x && d2.p1.y < d1.p2.y && d2.p2.y > d1.p1.y;
        }
    } else if (edge == RIGHT) {
        if (isPx) {
            return d1.p2.x == d2.p1.x + _globalExt[0] && d2.p1.y < d1.p2.y && d2.p2.y > d1.p1.y;
        } else {
            return d1.p2.x == d2.p1.x && d2.p1.y < d1.p2.y && d2.p2.y > d1.p1.y;
        }
    } else {
        std::cerr << "ERROR: edge must be LEFT, RIGHT, BOTTOM, TOP." << std::endl;
        exit(EXIT_FAILURE);
    }
}

bool Partitioner::isCornerNeighbour(
    const Domain d1, const Domain d2, const Corner corner, const bool isPx, const bool isPy)
{
    // create helper vars for domain 1
    auto left = d1.p1.x;
    auto right = d1.p2.x;
    auto top = d1.p2.y;
    auto bottom = d1.p1.y;

    // adjust for periodic boundaries if domain lies on one of the outer boundaries
    if (isPx) {
        if (left == 0) {
            left = _globalExt[0];
        }
        if (right == _globalExt[0]) {
            right = 0;
        }
    }
    if (isPy) {
        if (bottom == 0) {
            bottom = _globalExt[1];
        }
        if (top == _globalExt[1]) {
            top = 0;
        }
    }
    if (corner == TOP_LEFT) {
        // Check if top and left coordinates of domain 1 fall in domain 2 (d2)
        // Similar logic applies to the other corners
        return top >= d2.p1.y && top < d2.p2.y && left > d2.p1.x && left <= d2.p2.x;
    } else if (corner == TOP_RIGHT) {
        return top >= d2.p1.y && top < d2.p2.y && right < d2.p2.x && right >= d2.p1.x;
    } else if (corner == BOTTOM_RIGHT) {
        return bottom <= d2.p2.y && bottom > d2.p1.y && right >= d2.p1.x && right < d2.p2.x;
    } else if (corner == BOTTOM_LEFT) {
        return bottom <= d2.p2.y && bottom > d2.p1.y && left <= d2.p2.x && left > d2.p1.x;
    } else {
        std::cerr << "ERROR: corner must be TOP_LEFT, TOP_RIGHT, BOTTOM_RIGHT, BOTTOM_LEFT."
                  << std::endl;
        exit(EXIT_FAILURE);
    }
}

void Partitioner::haloEdgeBufferPositions(
    const Domain d1, const Domain d2, const Edge edge, int& sendPos, int& recvPos)
{
    // sendPos is the index where we will read the halo data from processes sending their data
    // recvPos is the index where we will write the halo data into that rank's recv buffer

    // Both indices (sendPos and recvPos) are relative indices in the send and recv buffers used
    // in NextSim. The dimension of each buffer may be different for each rank. Each rank will have
    // it's own send and recv buffer. The size of each buffer will depend on each domain's
    // perimeter. The send buffer for each rank is formed by taking the outer edges of the 2D domain
    // and unrolling them into a 1D array in the order of Bottom, Right, Top and Left edges. The
    // recv buffer is formed from the data gathered during the halo exchange. It is also laid out a
    // similar way in memory.

    // Detailed description of Halo exchange logic is available at:
    // https://nextsim-dg.readthedocs.io/en/latest/halo-exchange.html

    sendPos = 0;
    if (edge == TOP) {
        // dx is the offset between domains
        int dx = std::max(d1.p1.x, d2.p1.x) - d2.p1.x;
        // in this case d2 is the TOP nieghbour of d1, which means it will need to share elements
        // from it's BOTTOM edge, which is first in the 1D perimeter array. Therefore there is no
        // additional offset
        sendPos = dx;
    } else if (edge == BOTTOM) {
        // dx is the offset between domains
        int dx = std::max(d1.p1.x, d2.p1.x) - d2.p1.x;
        // in this case d2 is the BOTTOM nieghbour of d1, which means it will need to share elements
        // from it's TOP edge, which comes third in the 1D perimeter array (i.e., Bottom, Right and
        // then Top). Therefore we have to account for an additional offset.
        sendPos = d2.getHeight() + d2.getWidth() + dx;
    } else if (edge == LEFT) {
        int dy = std::max(d1.p1.y, d2.p1.y) - d2.p1.y;
        sendPos = d2.getWidth() + dy;
    } else if (edge == RIGHT) {
        int dy = std::max(d1.p1.y, d2.p1.y) - d2.p1.y;
        sendPos = 2 * d2.getWidth() + d2.getHeight() + dy;
    } else {
        std::cerr << "ERROR: edge must be LEFT, RIGHT, BOTTOM, TOP." << std::endl;
        exit(EXIT_FAILURE);
    }
    // this logic here is similar to sendPos but it is a mirror reflection between L<->R and T<->B
    // and 1<->2
    recvPos = 0;
    if (edge == TOP) {
        int dx = std::max(d1.p1.x, d2.p1.x) - d1.p1.x;
        recvPos = d1.getHeight() + d1.getWidth() + dx;
    } else if (edge == BOTTOM) {
        int dx = std::max(d1.p1.x, d2.p1.x) - d1.p1.x;
        recvPos = dx;
    } else if (edge == LEFT) {
        int dy = std::max(d1.p1.y, d2.p1.y) - d1.p1.y;
        recvPos = 2 * d1.getWidth() + d1.getHeight() + dy;
    } else if (edge == RIGHT) {
        int dy = std::max(d1.p1.y, d2.p1.y) - d1.p1.y;
        recvPos = d1.getWidth() + dy;
    } else {
        std::cerr << "ERROR: edge must be LEFT, RIGHT, BOTTOM, TOP." << std::endl;
        exit(EXIT_FAILURE);
    }
}

void Partitioner::haloCornerBufferPositions(
    const Domain d1, const Domain d2, const Corner corner, int& sendPos)
{
    int globalX = _globalExt[0];
    int globalY = _globalExt[1];
    int xpos, ypos;

    sendPos = 0;
    if (corner == TOP_RIGHT) {
        // for a TOP_RIGHT corner the corner neighbour can either be along the other domains bottom
        // or left edge. We need to check so we know where to look in the send buffer.
        xpos = mod(d1.p2.x, globalX);
        ypos = mod(d1.p2.y, globalY);
        if (ypos >= d2.p1.y) {
            // this first case is true if the corner neighbour lies on the left edge
            int dy = ypos - d2.p1.y;
            sendPos = 2 * d2.getWidth() + d2.getHeight() + dy;
        } else {
            // this second case works if the corner neighbour lies on the bottom edge (or in the
            // corner of both e.g., the bottom left corner)
            int dx = xpos - d2.p1.x;
            sendPos = dx;
        }
    } else if (corner == TOP_LEFT) {
        xpos = mod(d1.p1.x - 1, globalX);
        ypos = mod(d1.p2.y, globalY);
        if (ypos >= d2.p1.y) {
            int dy = ypos - d2.p1.y;
            sendPos = d2.getWidth() + dy;
        } else {
            int dx = xpos - d2.p1.x;
            sendPos = dx;
        }
    } else if (corner == BOTTOM_LEFT) {
        xpos = mod(d1.p1.x - 1, globalX);
        ypos = mod(d1.p1.y - 1, globalY);
        if (d2.p2.y >= ypos) {
            int dy = ypos - d2.p1.y;
            sendPos = d2.getWidth() + dy;
        } else {
            int dx = xpos - d2.p1.x;
            sendPos = d2.getWidth() + d2.getHeight() + dx;
        }
    } else if (corner == BOTTOM_RIGHT) {
        xpos = mod(d1.p2.x, globalX);
        ypos = mod(d1.p1.y - 1, globalY);
        if (d2.p2.y >= ypos) {
            int dy = ypos - d2.p1.y;
            sendPos = 2 * d2.getWidth() + d2.getHeight() + dy;
        } else {
            int dx = xpos - d2.p1.x;
            sendPos = d2.getWidth() + d2.getHeight() + dx;
        }
    } else {
        std::cerr << "ERROR: corner must be TOP_LEFT, TOP_RIGHT, BOTTOM_LEFT, BOTTOM_RIGHT."
                  << std::endl;
        exit(EXIT_FAILURE);
    }
}

Partitioner::Partitioner(MPI_Comm comm)
{
    _comm = comm;
    CHECK_MPI(MPI_Comm_size(comm, &_totalNumProcs));
    CHECK_MPI(MPI_Comm_rank(comm, &_rank));
}

void Partitioner::getBoundingBox(int& global0, int& global1, int& localExt0, int& localExt1) const
{
    global0 = _globalNew[0];
    global1 = _globalNew[1];
    localExt0 = _localExtNew[0];
    localExt1 = _localExtNew[1];
}

void Partitioner::getNeighbourInfo(std::array<std::vector<int>, N_EDGE>& ids,
    std::array<std::vector<int>, N_EDGE>& haloSizes, std::array<std::vector<int>, N_EDGE>& haloSend,
    std::array<std::vector<int>, N_EDGE>& haloRecv,
    std::array<std::vector<int>, N_CORNER>& cornerIds,
    std::array<std::vector<int>, N_CORNER>& cornerSend) const
{
    for (auto edge : edges) {
        for (auto it = _neighbours[edge].begin(); it != _neighbours[edge].end(); ++it) {
            ids[edge].push_back(it->first);
            haloSizes[edge].push_back(it->second);
            haloSend[edge].push_back(_sendPos[edge].at(it->first));
            haloRecv[edge].push_back(_recvPos[edge].at(it->first));
        }
    }

    for (auto corner : corners) {
        for (auto it = _cornerNeighbours[corner].begin(); it != _cornerNeighbours[corner].end();
             ++it) {
            cornerIds[corner].push_back(it->first);
            cornerSend[corner].push_back(_cornerSendPos[corner].at(it->first));
        }
    }
}

void Partitioner::saveMask(const std::string& filename) const
{
    // Use C API for parallel I/O
    int nc_id, nc_mode;
    nc_mode = NC_CLOBBER | NC_NETCDF4;
    NC_CHECK(nc_create_par(filename.c_str(), nc_mode, _comm, MPI_INFO_NULL, &nc_id));
    NC_CHECK(nc_put_att_int(nc_id, NC_GLOBAL, "num_processes", NC_INT, 1, &_totalNumProcs));

    // Create 2 dimensions
    // The values to be written are associated with the netCDF variable by
    // assuming that the last dimension of the netCDF variable varies fastest in
    // the C interface
    int dimid[NDIMS];
    // for nextsimdg code we always want the output in the order of yx
    std::vector<std::string> dim_chars = { "y", "x" };
    for (int idx = 0; idx < NDIMS; idx++) {
        NC_CHECK(
            nc_def_dim(nc_id, dim_chars[idx].c_str(), _globalExt[NDIMS - 1 - idx], &dimid[idx]));
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
        count[idx] = _localExt[NDIMS - 1 - idx];
    }

    // Store data
    NC_CHECK(nc_var_par_access(nc_id, mask_nc_id, NC_COLLECTIVE));
    NC_CHECK(nc_put_vara_int(nc_id, mask_nc_id, start, count, _procId.data()));
    NC_CHECK(nc_close(nc_id));
}

struct DimInfo {
    std::vector<int> numNeighbours;
    std::vector<int> dims;
    std::vector<int> offsets;
};

template <typename T, std::size_t N>
static DimInfo compute_dims(
    const std::array<T, N>& items, const std::array<std::vector<int>, N>& data, MPI_Comm comm)
{
    DimInfo info;
    info.numNeighbours.resize(items.size());
    info.dims.resize(items.size(), 0);
    info.offsets.resize(items.size(), 0);
    for (std::size_t i = 0; i < items.size(); i++) {
        info.numNeighbours[i] = (int)data[items[i]].size();
        CHECK_MPI(MPI_Allreduce(&info.numNeighbours[i], &info.dims[i], 1, MPI_INT, MPI_SUM, comm));
        CHECK_MPI(MPI_Exscan(&info.numNeighbours[i], &info.offsets[i], 1, MPI_INT, MPI_SUM, comm));
    }
    return info;
}

void Partitioner::saveMetadata(const std::string& filename) const
{
    // Use C API for parallel I/O
    int nc_id, nc_mode;
    nc_mode = NC_MPIIO | NC_NETCDF4;
    NC_CHECK(nc_create_par(filename.c_str(), nc_mode, _comm, MPI_INFO_NULL, &nc_id));

    // Mark tripolar topology by attribute presence (written only when set).
    if (_tripolar) {
        int one = 1;
        NC_CHECK(nc_put_att_int(nc_id, NC_GLOBAL, "tripolar", NC_INT, 1, &one));
    }

    // utility lambdas for netcdf operations
    // define a new dimension
    auto def_dim = [&](const std::string& name, int len, int& dimid) {
        NC_CHECK(nc_def_dim(nc_id, name.c_str(), len, &dimid));
    };

    // write an array to netcdf file
    auto write_array = [&](int gid, int vid, size_t start, size_t count, const int* data) {
        NC_CHECK(nc_var_par_access(gid, vid, NC_COLLECTIVE));
        NC_CHECK(nc_put_vara_int(gid, vid, &start, &count, data));
    };

    // write a scalar to netcdf file
    auto write_scalar = [&](int gid, int vid, size_t start, int val) {
        NC_CHECK(nc_var_par_access(gid, vid, NC_COLLECTIVE));
        NC_CHECK(nc_put_var1_int(gid, vid, &start, &val));
    };

    // ---- Global dimensions (NX, NY) ----
    int dimid_global[Partitioner::NDIMS];
    for (int idx = 0; idx < Partitioner::NDIMS; idx++) {
        def_dim(globalExtentNames[idx], _globalExt[idx], dimid_global[idx]);
    }

    // ---- Prepare neighbour data ----
    std::array<std::vector<int>, N_EDGE> ids, halos, haloSend, haloRecv;
    std::array<std::vector<int>, N_CORNER> corner_ids, cornerSend;
    getNeighbourInfo(ids, halos, haloSend, haloRecv, corner_ids, cornerSend);

    // ---- Compute edge dimensions (MPI Allreduce/Exscan) ----
    DimInfo edge_info = compute_dims(edges, ids, _comm);
    DimInfo corner_info = compute_dims(corners, corner_ids, _comm);

    // ---- Define netCDF dimensions ----
    int dimid;
    std::vector<int> edge_dimids(N_EDGE);
    std::vector<int> corner_dimids(N_CORNER);

    NC_CHECK(nc_def_dim(nc_id, "P", _totalNumProcs, &dimid));
    for (auto edge : edges) {
        def_dim(dir_chars[edge], edge_info.dims[edge], edge_dimids[edge]);
    }
    for (auto corner : corners) {
        def_dim(corner_dir_chars[corner], corner_info.dims[corner], corner_dimids[corner]);
    }

    // ---- Define groups ----
    int bbox_gid, c_grid;
    NC_CHECK(nc_def_grp(nc_id, "bounding_boxes", &bbox_gid));
    NC_CHECK(nc_def_grp(nc_id, "connectivity", &c_grid));

    // ---- Define variables: bounding boxes ----
    int top_vid[Partitioner::NDIMS];
    int cnt_vid[Partitioner::NDIMS];
    auto def_bbox_var = [&](const std::string& prefix, const auto idx, int& vid) {
        NC_CHECK(nc_def_var(bbox_gid, (prefix + dim_chars[idx]).c_str(), NC_INT, 1, &dimid, &vid));
    };
    for (int idx = 0; idx < Partitioner::NDIMS; idx++) {
        def_bbox_var("domain_", idx, top_vid[idx]);
        def_bbox_var("domain_extent_", idx, cnt_vid[idx]);
    }

    // ---- Define variables: edge connectivity ----
    int num_vid[N_EDGE], ids_vid[N_EDGE], halos_vid[N_EDGE];
    int haloSend_vid[N_EDGE], haloRecv_vid[N_EDGE];

    // lambda to define netcdf variables (for edges)
    auto def_edge_var = [&](const auto edge, const std::string& name, const int* dimArr,
                            auto* vid) {
        NC_CHECK(
            nc_def_var(c_grid, (dir_names[edge] + name).c_str(), NC_INT, 1, dimArr, &vid[edge]));
    };
    for (auto edge : edges) {
        def_edge_var(edge, "_neighbours", &dimid, num_vid);
        const auto* edgeDimId = &edge_dimids[edge];
        def_edge_var(edge, "_neighbour_ids", edgeDimId, ids_vid);
        def_edge_var(edge, "_neighbour_halos", edgeDimId, halos_vid);
        def_edge_var(edge, "_neighbour_halo_send", edgeDimId, haloSend_vid);
        def_edge_var(edge, "_neighbour_halo_recv", edgeDimId, haloRecv_vid);
    }

    // ---- Define variables: corner connectivity ----
    int num_corner_vid[N_CORNER], ids_corner_vid[N_CORNER];
    int cornerSend_vid[N_CORNER];
    // lambda to define netcdf variables (for corners)
    auto def_corner_var
        = [&](const auto corner, const std::string& name, const int* dimArr, int* vid) {
              NC_CHECK(nc_def_var(c_grid, (corner_dir_names[corner] + name).c_str(), NC_INT, 1,
                  dimArr, &vid[corner]));
          };
    for (auto corner : corners) {
        def_corner_var(corner, "_neighbours", &dimid, num_corner_vid);
        const auto* cornerDimId = &corner_dimids[corner];
        def_corner_var(corner, "_neighbour_ids", cornerDimId, ids_corner_vid);
        def_corner_var(corner, "_neighbour_send", cornerDimId, cornerSend_vid);
    }

    // ---- Write ----
    NC_CHECK(nc_enddef(nc_id));

    // Bounding boxes: one value per process
    for (int idx = 0; idx < Partitioner::NDIMS; idx++) {
        size_t start = _rank;
        write_scalar(bbox_gid, top_vid[idx], start, _globalNew[idx]);
        write_scalar(bbox_gid, cnt_vid[idx], start, _localExtNew[idx]);
    }

    // Edge connectivity
    for (auto edge : edges) {
        size_t start = _rank;
        size_t count = edge_info.numNeighbours[edge];
        write_scalar(c_grid, num_vid[edge], start, edge_info.numNeighbours[edge]);
        start = edge_info.offsets[edge];
        write_array(c_grid, ids_vid[edge], start, count, ids[edge].data());
        write_array(c_grid, halos_vid[edge], start, count, halos[edge].data());
        write_array(c_grid, haloSend_vid[edge], start, count, haloSend[edge].data());
        write_array(c_grid, haloRecv_vid[edge], start, count, haloRecv[edge].data());
    }

    // Corner connectivity
    for (auto corner : corners) {
        size_t start = _rank;
        size_t count = corner_info.numNeighbours[corner];
        write_scalar(c_grid, num_corner_vid[corner], start, corner_info.numNeighbours[corner]);
        start = corner_info.offsets[corner];
        write_array(c_grid, ids_corner_vid[corner], start, count, corner_ids[corner].data());
        write_array(c_grid, cornerSend_vid[corner], start, count, cornerSend[corner].data());
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
    std::vector<Point> origins(_totalNumProcs);
    std::vector<Point> extents(_totalNumProcs);
    std::vector<Domain> domains(_totalNumProcs);
    std::vector<int> tmp0(_totalNumProcs);
    std::vector<int> tmp1(_totalNumProcs);

    CHECK_MPI(MPI_Allgather(&_globalNew[0], 1, MPI_INT, tmp0.data(), 1, MPI_INT, _comm));
    CHECK_MPI(MPI_Allgather(&_globalNew[1], 1, MPI_INT, tmp1.data(), 1, MPI_INT, _comm));

    // origin points mark the bottom-left corner of each domain
    for (int p = 0; p < _totalNumProcs; p++) {
        origins[p].x = tmp0[p];
        origins[p].y = tmp1[p];
    }

    CHECK_MPI(MPI_Allgather(&_localExtNew[0], 1, MPI_INT, tmp0.data(), 1, MPI_INT, _comm));
    CHECK_MPI(MPI_Allgather(&_localExtNew[1], 1, MPI_INT, tmp1.data(), 1, MPI_INT, _comm));

    // extents can be used to find the top-right corner of each domain
    for (int p = 0; p < _totalNumProcs; p++) {
        extents[p].x = tmp0[p];
        extents[p].y = tmp1[p];
    }

    // generate domains
    for (int p = 0; p < _totalNumProcs; p++) {
        domains[p].p1.x = origins[p].x;
        domains[p].p1.y = origins[p].y;
        domains[p].p2.x = origins[p].x + extents[p].x;
        domains[p].p2.y = origins[p].y + extents[p].y;
    }

    const bool usePeriodicX = _px || _tripolar;
    const bool usePeriodicY = _py;

    for (int p = 0; p < _totalNumProcs; p++) {

        // When finding neighbours *within* the domain, we don't check against the current rank
        // because a subdomain can't be a neighbour of itself.
        if (p != _rank) {

            // check edge neighours
            for (auto edge : edges) {
                if (isNeighbour(domains[_rank], domains[p], edge)) {
                    int haloSize = domainOverlap(domains[_rank], domains[p], edge);
                    if (haloSize > 0) {
                        _neighbours[edge].insert(std::pair<int, int>(p, haloSize));
                        int sendPos = 0;
                        int recvPos = 0;
                        haloEdgeBufferPositions(domains[_rank], domains[p], edge, sendPos, recvPos);
                        _sendPos[edge].insert(std::pair<int, int>(p, sendPos));
                        _recvPos[edge].insert(std::pair<int, int>(p, recvPos));
                    }
                }
            }

            // check corner neighours
            for (auto corner : corners) {
                if (isCornerNeighbour(domains[_rank], domains[p], corner)) {
                    _cornerNeighbours[corner].insert(std::pair<int, int>(p, 1));
                    int sendPos = 0;
                    haloCornerBufferPositions(domains[_rank], domains[p], corner, sendPos);
                    _cornerSendPos[corner].insert(std::pair<int, int>(p, sendPos));
                }
            }
        }

        // When finding neighbours *across periodic boundaries*, we need to check against the
        // current rank, too, because a subdomain can be a periodic neighbour of itself.
        // check periodic edge neighbours
        for (auto edge : edges) {
            if (isNeighbour(domains[_rank], domains[p], edge, usePeriodicX, usePeriodicY)) {
                int haloSize = domainOverlap(domains[_rank], domains[p], edge);
                if (haloSize > 0) {
                    _neighbours[edge].insert(std::pair<int, int>(p, haloSize));
                    int sendPos = 0;
                    int recvPos = 0;
                    haloEdgeBufferPositions(domains[_rank], domains[p], edge, sendPos, recvPos);
                    _sendPos[edge].insert(std::pair<int, int>(p, sendPos));
                    _recvPos[edge].insert(std::pair<int, int>(p, recvPos));
                }
            }
        }

        // check periodic corner neighbours
        for (auto corner : corners) {
            if (isCornerNeighbour(domains[_rank], domains[p], corner, usePeriodicX, usePeriodicY)) {
                _cornerNeighbours[corner].insert(std::pair<int, int>(p, 1));
                int sendPos = 0;
                haloCornerBufferPositions(domains[_rank], domains[p], corner, sendPos);
                _cornerSendPos[corner].insert(std::pair<int, int>(p, sendPos));
            }
        }
    }

    // For tripolar topology apply special treatment to the top edge
    // To correctly resolve the neighbourhood relations, we just need to create
    // an 'image' of the domain on the other side of the tripolar top edge.
    // This is just a point symmetry through the "middle point" of the top edge
    //
    // The only difficulty comes from calculating the send and receive buffer
    // positions.
    //
    // For `send` buffers these are calculated on the 'target' domain.
    // Hence we need to use the image of the current domain to get correct orientation.
    //
    // For the `recv` buffers, these are calculated on the 'current' domain.
    // Hence we use the image of the target
    //
    if (_tripolar) {
        const Domain& thisDomain = domains[_rank];
        const Point symmetryPoint = { _globalExt[0] / 2, _globalExt[1] };

        // Resolve Edge Neighbours
        for (int p = 0; p < _totalNumProcs; p++) {
            const Domain image = pointReflection(symmetryPoint, domains[p]);
            const Domain selfImage = pointReflection(symmetryPoint, thisDomain);

            // Now we can check for the neighbourhood relation with the image
            if (isNeighbour(thisDomain, image, TOP)) {
                const int haloSize = domainOverlap(thisDomain, image, TOP);

                if (haloSize > 0) {
                    _neighbours[TOP].insert(std::pair<int, int>(p, haloSize));

                    int sendPos, recvPos, dontCare;

                    // Calculate the receive position on self
                    haloEdgeBufferPositions(thisDomain, image, TOP, dontCare, recvPos);

                    // Calculate the send position on the image
                    haloEdgeBufferPositions(selfImage, domains[p], BOTTOM, sendPos, dontCare);

                    _sendPos[TOP].insert(std::pair<int, int>(p, sendPos));
                    _recvPos[TOP].insert(std::pair<int, int>(p, recvPos));
                }
            }

            if (isCornerNeighbour(thisDomain, image, TOP_LEFT)) {
                _cornerNeighbours[TOP_LEFT].insert(std::pair<int, int>(p, 1));

                int sendPos;
                haloCornerBufferPositions(selfImage, domains[p], BOTTOM_RIGHT, sendPos);

                _cornerSendPos[TOP_LEFT].insert(std::pair<int, int>(p, sendPos));
            }

            if (isCornerNeighbour(thisDomain, image, TOP_RIGHT)) {
                _cornerNeighbours[TOP_RIGHT].insert(std::pair<int, int>(p, 1));

                int sendPos;
                haloCornerBufferPositions(selfImage, domains[p], BOTTOM_LEFT, sendPos);

                _cornerSendPos[TOP_RIGHT].insert(std::pair<int, int>(p, sendPos));
            }
        }
    }
}
