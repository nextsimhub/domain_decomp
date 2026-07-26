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

void Partitioner::getNeighbourInfoPeriodic(std::array<std::vector<int>, N_EDGE>& ids,
    std::array<std::vector<int>, N_EDGE>& haloSizes, std::array<std::vector<int>, N_EDGE>& haloSend,
    std::array<std::vector<int>, N_EDGE>& haloRecv,
    std::array<std::vector<int>, N_CORNER>& cornerIds,
    std::array<std::vector<int>, N_CORNER>& cornerSend) const
{
    for (auto edge : edges) {
        if (((edge == LEFT || edge == RIGHT) && _px) || ((edge == TOP || edge == BOTTOM) && _py)) {

            for (auto it = _neighbours_p[edge].begin(); it != _neighbours_p[edge].end(); ++it) {
                ids[edge].push_back(it->first);
                haloSizes[edge].push_back(it->second);
                haloSend[edge].push_back(_sendPos_p[edge].at(it->first));
                haloRecv[edge].push_back(_recvPos_p[edge].at(it->first));
            }
        }
    }

    for (auto corner : corners) {
        for (auto it = _cornerNeighbours_p[corner].begin(); it != _cornerNeighbours_p[corner].end();
             ++it) {
            cornerIds[corner].push_back(it->first);
            cornerSend[corner].push_back(_cornerSendPos_p[corner].at(it->first));
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

void Partitioner::saveMetadata(const std::string& filename) const
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
        NC_CHECK(
            nc_def_dim(nc_id, globalExtentNames[idx].c_str(), _globalExt[idx], &dimid_global[idx]));
    }

    // Prepare neighbour data
    std::array<std::vector<int>, N_EDGE> ids, halos, haloSend, haloRecv;
    std::array<std::vector<int>, N_CORNER> corner_ids, cornerSend;
    getNeighbourInfo(ids, halos, haloSend, haloRecv, corner_ids, cornerSend);

    std::vector<int> num_neighbours(N_EDGE), dims(N_EDGE, 0), offsets(N_EDGE, 0);
    for (auto edge : edges) {
        num_neighbours[edge] = (int)ids[edge].size();
        CHECK_MPI(MPI_Allreduce(&num_neighbours[edge], &dims[edge], 1, MPI_INT, MPI_SUM, _comm));
        CHECK_MPI(MPI_Exscan(&num_neighbours[edge], &offsets[edge], 1, MPI_INT, MPI_SUM, _comm));
    }

    std::vector<int> numCornerNeighbours(N_CORNER), corner_dims(N_CORNER, 0),
        corner_offsets(N_CORNER, 0);
    for (auto corner : corners) {
        numCornerNeighbours[corner] = (int)corner_ids[corner].size();
        CHECK_MPI(MPI_Allreduce(
            &numCornerNeighbours[corner], &corner_dims[corner], 1, MPI_INT, MPI_SUM, _comm));
        CHECK_MPI(MPI_Exscan(
            &numCornerNeighbours[corner], &corner_offsets[corner], 1, MPI_INT, MPI_SUM, _comm));
    }

    // Prepare periodic neighbour data
    std::array<std::vector<int>, N_EDGE> ids_p, halos_p, haloSend_p, haloRecv_p;
    std::array<std::vector<int>, N_CORNER> corner_ids_p, cornerSend_p;
    getNeighbourInfoPeriodic(ids_p, halos_p, haloSend_p, haloRecv_p, corner_ids_p, cornerSend_p);
    std::vector<int> num_neighbours_p(N_EDGE), dims_p(N_EDGE, 0), offsets_p(N_EDGE, 0);
    for (auto edge : edges) {
        num_neighbours_p[edge] = (int)ids_p[edge].size();
        CHECK_MPI(
            MPI_Allreduce(&num_neighbours_p[edge], &dims_p[edge], 1, MPI_INT, MPI_SUM, _comm));
        CHECK_MPI(
            MPI_Exscan(&num_neighbours_p[edge], &offsets_p[edge], 1, MPI_INT, MPI_SUM, _comm));
    }

    std::vector<int> numCornerNeighbours_p(N_CORNER), corner_dims_p(N_CORNER, 0),
        corner_offsets_p(N_CORNER, 0);
    for (auto corner : corners) {
        numCornerNeighbours_p[corner] = (int)corner_ids_p[corner].size();
        CHECK_MPI(MPI_Allreduce(
            &numCornerNeighbours_p[corner], &corner_dims_p[corner], 1, MPI_INT, MPI_SUM, _comm));
        CHECK_MPI(MPI_Exscan(
            &numCornerNeighbours_p[corner], &corner_offsets_p[corner], 1, MPI_INT, MPI_SUM, _comm));
    }

    // Define dimensions in netCDF file
    int dimid;
    std::vector<int> dimids(N_EDGE);
    NC_CHECK(nc_def_dim(nc_id, "P", _totalNumProcs, &dimid));
    for (auto edge : edges) {
        NC_CHECK(nc_def_dim(nc_id, dir_chars[edge].c_str(), dims[edge], &dimids[edge]));
    }

    std::vector<int> corner_dimids(N_CORNER);
    for (auto corner : corners) {
        NC_CHECK(nc_def_dim(
            nc_id, corner_dir_chars[corner].c_str(), corner_dims[corner], &corner_dimids[corner]));
    }

    // Define periodic dimensions in netCDF file
    std::vector<int> dimids_p(N_EDGE);
    for (auto edge : edges) {
        NC_CHECK(nc_def_dim(
            nc_id, (dir_chars[edge] + "_periodic").c_str(), dims_p[edge], &dimids_p[edge]));
    }

    std::vector<int> corner_dimids_p(N_CORNER);
    for (auto edge : edges) {
        NC_CHECK(nc_def_dim(nc_id, (corner_dir_chars[edge] + "_periodic").c_str(),
            corner_dims_p[edge], &corner_dimids_p[edge]));
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
    int haloSend_vid[N_EDGE];
    int haloRecv_vid[N_EDGE];
    for (auto edge : edges) {
        // Connectivity group
        NC_CHECK(nc_def_var(connectivity_gid, (dir_names[edge] + "_neighbours").c_str(), NC_INT, 1,
            &dimid, &num_vid[edge]));
        NC_CHECK(nc_def_var(connectivity_gid, (dir_names[edge] + "_neighbour_ids").c_str(), NC_INT,
            1, &dimids[edge], &ids_vid[edge]));
        NC_CHECK(nc_def_var(connectivity_gid, (dir_names[edge] + "_neighbour_halos").c_str(),
            NC_INT, 1, &dimids[edge], &halos_vid[edge]));
        NC_CHECK(nc_def_var(connectivity_gid, (dir_names[edge] + "_neighbour_halo_send").c_str(),
            NC_INT, 1, &dimids[edge], &haloSend_vid[edge]));
        NC_CHECK(nc_def_var(connectivity_gid, (dir_names[edge] + "_neighbour_halo_recv").c_str(),
            NC_INT, 1, &dimids[edge], &haloRecv_vid[edge]));
    }

    int num_corner_vid[N_CORNER];
    int ids_corner_vid[N_CORNER];
    int cornerSend_vid[N_CORNER];
    for (auto corner : corners) {
        // Connectivity group
        NC_CHECK(nc_def_var(connectivity_gid, (corner_dir_names[corner] + "_neighbours").c_str(),
            NC_INT, 1, &dimid, &num_corner_vid[corner]));
        NC_CHECK(nc_def_var(connectivity_gid, (corner_dir_names[corner] + "_neighbour_ids").c_str(),
            NC_INT, 1, &corner_dimids[corner], &ids_corner_vid[corner]));
        NC_CHECK(
            nc_def_var(connectivity_gid, (corner_dir_names[corner] + "_neighbour_send").c_str(),
                NC_INT, 1, &corner_dimids[corner], &cornerSend_vid[corner]));
    }

    int num_vid_p[N_EDGE];
    int ids_vid_p[N_EDGE];
    int halos_vid_p[N_EDGE];
    int haloSend_vid_p[N_EDGE];
    int haloRecv_vid_p[N_EDGE];
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
            (dir_names[edge] + "_neighbour_halo_send_periodic").c_str(), NC_INT, 1, &dimids_p[edge],
            &haloSend_vid_p[edge]));
        NC_CHECK(nc_def_var(connectivity_gid,
            (dir_names[edge] + "_neighbour_halo_recv_periodic").c_str(), NC_INT, 1, &dimids_p[edge],
            &haloRecv_vid_p[edge]));
    }

    int num_corner_vid_p[N_CORNER];
    int ids_corner_vid_p[N_CORNER];
    int cornerSend_vid_p[N_CORNER];
    for (auto corner : corners) {
        // Connectivity group
        NC_CHECK(nc_def_var(connectivity_gid,
            (corner_dir_names[corner] + "_neighbours_periodic").c_str(), NC_INT, 1, &dimid,
            &num_corner_vid_p[corner]));
        NC_CHECK(nc_def_var(connectivity_gid,
            (corner_dir_names[corner] + "_neighbour_ids_periodic").c_str(), NC_INT, 1,
            &corner_dimids_p[corner], &ids_corner_vid_p[corner]));
        NC_CHECK(nc_def_var(connectivity_gid,
            (corner_dir_names[corner] + "_neighbour_send_periodic").c_str(), NC_INT, 1,
            &corner_dimids_p[corner], &cornerSend_vid_p[corner]));
    }

    // Write metadata to file
    NC_CHECK(nc_enddef(nc_id));

    // Store data
    for (int idx = 0; idx < NDIMS; idx++) {
        size_t start = _rank;
        NC_CHECK(nc_var_par_access(bbox_gid, top_vid[idx], NC_COLLECTIVE));
        NC_CHECK(nc_put_var1_int(bbox_gid, top_vid[idx], &start, &_globalNew[idx]));
        NC_CHECK(nc_var_par_access(bbox_gid, cnt_vid[idx], NC_COLLECTIVE));
        NC_CHECK(nc_put_var1_int(bbox_gid, cnt_vid[idx], &start, &_localExtNew[idx]));
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
        NC_CHECK(nc_var_par_access(connectivity_gid, haloSend_vid[edge], NC_COLLECTIVE));
        NC_CHECK(nc_put_vara_int(
            connectivity_gid, haloSend_vid[edge], &start, &count, haloSend[edge].data()));
        NC_CHECK(nc_put_vara_int(
            connectivity_gid, haloRecv_vid[edge], &start, &count, haloRecv[edge].data()));
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
            connectivity_gid, haloSend_vid_p[edge], &start, &count, haloSend_p[edge].data()));
        NC_CHECK(nc_put_vara_int(
            connectivity_gid, haloRecv_vid_p[edge], &start, &count, haloRecv_p[edge].data()));
    }

    for (auto corner : corners) {
        // Numbers of "corner" neighbours
        size_t start = _rank;
        NC_CHECK(nc_var_par_access(connectivity_gid, num_corner_vid[corner], NC_COLLECTIVE));
        NC_CHECK(nc_put_var1_int(
            connectivity_gid, num_corner_vid[corner], &start, &numCornerNeighbours[corner]));
        // Numbers of corner neighbours for periodic dimensions
        NC_CHECK(nc_var_par_access(connectivity_gid, num_corner_vid_p[corner], NC_COLLECTIVE));
        NC_CHECK(nc_put_var1_int(
            connectivity_gid, num_corner_vid_p[corner], &start, &numCornerNeighbours_p[corner]));

        // "Corner" IDs and halos
        start = corner_offsets[corner];
        size_t count = numCornerNeighbours[corner];
        NC_CHECK(nc_var_par_access(connectivity_gid, ids_corner_vid[corner], NC_COLLECTIVE));
        NC_CHECK(nc_put_vara_int(
            connectivity_gid, ids_corner_vid[corner], &start, &count, corner_ids[corner].data()));
        NC_CHECK(nc_var_par_access(connectivity_gid, cornerSend_vid[corner], NC_COLLECTIVE));
        NC_CHECK(nc_put_vara_int(
            connectivity_gid, cornerSend_vid[corner], &start, &count, cornerSend[corner].data()));

        // "Corner" IDs and halos for periodic dimensions
        start = corner_offsets_p[corner];
        count = numCornerNeighbours_p[corner];
        NC_CHECK(nc_var_par_access(connectivity_gid, ids_corner_vid_p[corner], NC_COLLECTIVE));
        NC_CHECK(nc_put_vara_int(connectivity_gid, ids_corner_vid_p[corner], &start, &count,
            corner_ids_p[corner].data()));
        NC_CHECK(nc_var_par_access(connectivity_gid, halos_vid_p[corner], NC_COLLECTIVE));
        NC_CHECK(nc_put_vara_int(connectivity_gid, cornerSend_vid_p[corner], &start, &count,
            cornerSend_p[corner].data()));
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

void Partitioner::setTotalNumProcs(int totalNumProcs) { _totalNumProcs = totalNumProcs; }

void Partitioner::setComm(MPI_Comm comm) { _comm = comm; }

std::vector<int> Partitioner::getLocalExtNew() const { return _localExtNew; }

void Partitioner::setLocalExtNew(const std::vector<int>& localExtNew)
{
    _localExtNew = localExtNew;
}

std::vector<int> Partitioner::getGlobalNew() const { return _globalNew; }

void Partitioner::setGlobalNew(const std::vector<int>& globalNew) { _globalNew = globalNew; }

std::vector<int> Partitioner::getGlobal() const { return _global; }

void Partitioner::setGlobal(const std::vector<int>& global) { _global = global; }

std::vector<int> Partitioner::getGlobalExt() const { return _globalExt; }

void Partitioner::setGlobalExt(const std::vector<int>& globalExt) { _globalExt = globalExt; }

std::vector<int> Partitioner::getProcId() const { return _procId; }

void Partitioner::setProcId(const std::vector<int>& procId) { _procId = procId; }

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

        // When finding neighours *across periodic boundaries*, we need to check against the
        // current rank, too, because a subdomain can be a periodic neighbour of itself.
        // check periodic edge neighours
        for (auto edge : edges) {
            if (isNeighbour(domains[_rank], domains[p], edge, _px, _py)) {
                int haloSize = domainOverlap(domains[_rank], domains[p], edge);
                if (haloSize > 0) {
                    _neighbours_p[edge].insert(std::pair<int, int>(p, haloSize));
                    int sendPos = 0;
                    int recvPos = 0;
                    haloEdgeBufferPositions(domains[_rank], domains[p], edge, sendPos, recvPos);
                    _sendPos_p[edge].insert(std::pair<int, int>(p, sendPos));
                    _recvPos_p[edge].insert(std::pair<int, int>(p, recvPos));
                }
            }
        }

        // check periodic corner neighours
        for (auto corner : corners) {
            if (isCornerNeighbour(domains[_rank], domains[p], corner, _px, _py)) {
                if (_cornerNeighbours[corner].size() > 0) {
                    // skip if we have already counted as a non-periodic neighbour
                    continue;
                }
                _cornerNeighbours_p[corner].insert(std::pair<int, int>(p, 1));
                int sendPos = 0;
                haloCornerBufferPositions(domains[_rank], domains[p], corner, sendPos);
                _cornerSendPos_p[corner].insert(std::pair<int, int>(p, sendPos));
            }
        }
    }
}
