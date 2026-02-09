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

bool Partitioner::is_neighbour(
    const Domain d1, const Domain d2, const Edge edge, const bool is_px, const bool is_py)
{
    // For TOP & BOTTOM Edges check that the domains share a y-coordinate AND that the horizontal
    // overlap is non-zero

    // For LEFT & RIGHT Edges check that the domains share a x-coordinate AND that the vertical
    // overlap is non-zero
    if (edge == TOP) {
        if (is_py) {
            return d1.p2.y == d2.p1.y + _global_ext[1] && d2.p1.x < d1.p2.x && d2.p2.x > d1.p1.x;
        } else {
            return d1.p2.y == d2.p1.y && d2.p1.x < d1.p2.x && d2.p2.x > d1.p1.x;
        }
    } else if (edge == BOTTOM) {
        if (is_py) {
            return d1.p1.y == d2.p2.y - _global_ext[1] && d2.p1.x < d1.p2.x && d2.p2.x > d1.p1.x;
        } else {
            return d1.p1.y == d2.p2.y && d2.p1.x < d1.p2.x && d2.p2.x > d1.p1.x;
        }
    } else if (edge == LEFT) {
        if (is_px) {
            return d1.p1.x == d2.p2.x - _global_ext[0] && d2.p1.y < d1.p2.y && d2.p2.y > d1.p1.y;
        } else {
            return d1.p1.x == d2.p2.x && d2.p1.y < d1.p2.y && d2.p2.y > d1.p1.y;
        }
    } else if (edge == RIGHT) {
        if (is_px) {
            return d1.p2.x == d2.p1.x + _global_ext[0] && d2.p1.y < d1.p2.y && d2.p2.y > d1.p1.y;
        } else {
            return d1.p2.x == d2.p1.x && d2.p1.y < d1.p2.y && d2.p2.y > d1.p1.y;
        }
    } else {
        std::cerr << "ERROR: edge must be LEFT, RIGHT, BOTTOM, TOP." << std::endl;
        exit(EXIT_FAILURE);
    }
}

bool Partitioner::is_corner_neighbour(
    const Domain d1, const Domain d2, const Vertex vertex, const bool is_px, const bool is_py)
{
    // create helper vars for domain 1
    auto left = d1.p1.x;
    auto right = d1.p2.x;
    auto top = d1.p2.y;
    auto bottom = d1.p1.y;

    // adjust for periodic boundaries if domain lies on one of the outer boundaries
    if (is_px) {
        if (left == 0) {
            left = _global_ext[0];
        }
        if (right == _global_ext[0]) {
            right = 0;
        }
    }
    if (is_py) {
        if (bottom == 0) {
            bottom = _global_ext[1];
        }
        if (top == _global_ext[1]) {
            top = 0;
        }
    }
    if (vertex == TOP_LEFT) {
        // Check if top and left coordinates of domain 1 fall in domain 2 (d2)
        // Similar logic applies to the other corners
        return top >= d2.p1.y && top < d2.p2.y && left > d2.p1.x && left <= d2.p2.x;
    } else if (vertex == TOP_RIGHT) {
        return top >= d2.p1.y && top < d2.p2.y && right < d2.p2.x && right >= d2.p1.x;
    } else if (vertex == BOTTOM_RIGHT) {
        return bottom <= d2.p2.y && bottom > d2.p1.y && right >= d2.p1.x && right < d2.p2.x;
    } else if (vertex == BOTTOM_LEFT) {
        return bottom <= d2.p2.y && bottom > d2.p1.y && left <= d2.p2.x && left > d2.p1.x;
    } else {
        std::cerr << "ERROR: vertex must be TOP_LEFT, TOP_RIGHT, BOTTOM_RIGHT, BOTTOM_LEFT."
                  << std::endl;
        exit(EXIT_FAILURE);
    }
}

void Partitioner::haloEdgeBufferPositions(
    const Domain d1, const Domain d2, const Edge edge, int& send_pos, int& recv_pos)
{
    // send_pos is the index where we will read the halo data from processes sending their data
    // recv_pos is the index where we will write the halo data into that rank's recv buffer

    // Both indices (send_pos and recv_pos) are relative indices in the send and recv buffers used
    // in NextSim. The dimension of each buffer may be different for each rank. Each rank will have
    // it's own send and recv buffer. The size of each buffer will depend on each domain's
    // perimeter. The send buffer for each rank is formed by taking the outer edges of the 2D domain
    // and unrolling them into a 1D array in the order of Bottom, Right, Top and Left edges. The
    // recv buffer is formed from the data gathered during the halo exchange. It is also laid out a
    // similar way in memory.

    send_pos = 0;
    if (edge == TOP) {
        // dx is the offset between domains
        int dx = std::max(d1.p1.x, d2.p1.x) - d2.p1.x;
        // in this case d2 is the TOP nieghbour of d1, which means it will need to share elements
        // from it's BOTTOM edge, which is first in the 1D perimeter array. Therefore there is no
        // additional offset
        send_pos = dx;
    } else if (edge == BOTTOM) {
        // dx is the offset between domains
        int dx = std::max(d1.p1.x, d2.p1.x) - d2.p1.x;
        // in this case d2 is the BOTTOM nieghbour of d1, which means it will need to share elements
        // from it's TOP edge, which comes third in the 1D perimeter array (i.e., Bottom, Right and
        // then Top). Therefore we have to account for an additional offset.
        send_pos = d2.get_height() + d2.get_width() + dx;
    } else if (edge == LEFT) {
        int dy = std::max(d1.p1.y, d2.p1.y) - d2.p1.y;
        send_pos = d2.get_width() + dy;
    } else if (edge == RIGHT) {
        int dy = std::max(d1.p1.y, d2.p1.y) - d2.p1.y;
        send_pos = 2 * d2.get_width() + d2.get_height() + dy;
    } else {
        std::cerr << "ERROR: edge must be LEFT, RIGHT, BOTTOM, TOP." << std::endl;
        exit(EXIT_FAILURE);
    }
    // this logic here is similar to send_pos but it is a mirror reflection between L<->R and T<->B
    // and 1<->2
    recv_pos = 0;
    if (edge == TOP) {
        int dx = std::max(d1.p1.x, d2.p1.x) - d1.p1.x;
        recv_pos = d1.get_height() + d1.get_width() + dx;
    } else if (edge == BOTTOM) {
        int dx = std::max(d1.p1.x, d2.p1.x) - d1.p1.x;
        recv_pos = dx;
    } else if (edge == LEFT) {
        int dy = std::max(d1.p1.y, d2.p1.y) - d1.p1.y;
        recv_pos = 2 * d1.get_width() + d1.get_height() + dy;
    } else if (edge == RIGHT) {
        int dy = std::max(d1.p1.y, d2.p1.y) - d1.p1.y;
        recv_pos = d1.get_width() + dy;
    } else {
        std::cerr << "ERROR: edge must be LEFT, RIGHT, BOTTOM, TOP." << std::endl;
        exit(EXIT_FAILURE);
    }
}

void Partitioner::haloCornerBufferPositions(
    const Domain d1, const Domain d2, const Vertex vertex, int& send_pos)
{
    int globalX = _global_ext[0];
    int globalY = _global_ext[1];
    int xpos, ypos;

    send_pos = 0;
    if (vertex == TOP_RIGHT) {
        // for a TOP_RIGHT vertex the corner neighbour can either be along the other domains bottom
        // or left edge. We need to check so we know where to look in the send buffer.
        xpos = mod(d1.p2.x, globalX);
        ypos = mod(d1.p2.y, globalY);
        if (ypos >= d2.p1.y) {
            // this first case is true if the corner neighbour lies on the left edge
            int dy = ypos - d2.p1.y;
            send_pos = 2 * d2.get_width() + d2.get_height() + dy;
        } else {
            // this second case works if the corner neighbour lies on the bottom edge (or in the
            // corner of both e.g., the bottom left corner)
            int dx = xpos - d2.p1.x;
            send_pos = dx;
        }
    } else if (vertex == TOP_LEFT) {
        xpos = mod(d1.p1.x - 1, globalX);
        ypos = mod(d1.p2.y, globalY);
        if (ypos >= d2.p1.y) {
            int dy = ypos - d2.p1.y;
            send_pos = d2.get_width() + dy;
        } else {
            int dx = xpos - d2.p1.x;
            send_pos = dx;
        }
    } else if (vertex == BOTTOM_LEFT) {
        xpos = mod(d1.p1.x - 1, globalX);
        ypos = mod(d1.p1.y - 1, globalY);
        if (d2.p2.y >= ypos) {
            int dy = ypos - d2.p1.y;
            send_pos = d2.get_width() + dy;
        } else {
            int dx = xpos - d2.p1.x;
            send_pos = d2.get_width() + d2.get_height() + dx;
        }
    } else if (vertex == BOTTOM_RIGHT) {
        xpos = mod(d1.p2.x, globalX);
        ypos = mod(d1.p1.y - 1, globalY);
        if (d2.p2.y >= ypos) {
            int dy = ypos - d2.p1.y;
            send_pos = 2 * d2.get_width() + d2.get_height() + dy;
        } else {
            int dx = xpos - d2.p1.x;
            send_pos = d2.get_width() + d2.get_height() + dx;
        }
    } else {
        std::cerr << "ERROR: vertex must be TOP_LEFT, TOP_RIGHT, BOTTOM_LEFT, BOTTOM_RIGHT."
                  << std::endl;
        exit(EXIT_FAILURE);
    }
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

void Partitioner::get_neighbour_info(std::array<std::vector<int>, N_EDGE>& ids,
    std::array<std::vector<int>, N_EDGE>& halo_sizes,
    std::array<std::vector<int>, N_EDGE>& halo_send,
    std::array<std::vector<int>, N_EDGE>& halo_recv,
    std::array<std::vector<int>, N_VERTEX>& corner_ids,
    std::array<std::vector<int>, N_VERTEX>& corner_send) const
{
    for (auto edge : edges) {
        for (auto it = _neighbours[edge].begin(); it != _neighbours[edge].end(); ++it) {
            ids[edge].push_back(it->first);
            halo_sizes[edge].push_back(it->second);
            halo_send[edge].push_back(_send_pos[edge].at(it->first));
            halo_recv[edge].push_back(_recv_pos[edge].at(it->first));
        }
    }

    for (auto vertex : vertices) {
        for (auto it = _corner_neighbours[vertex].begin(); it != _corner_neighbours[vertex].end();
             ++it) {
            corner_ids[vertex].push_back(it->first);
            corner_send[vertex].push_back(_corner_send_pos[vertex].at(it->first));
        }
    }
}

void Partitioner::get_neighbour_info_periodic(std::array<std::vector<int>, N_EDGE>& ids,
    std::array<std::vector<int>, N_EDGE>& halo_sizes,
    std::array<std::vector<int>, N_EDGE>& halo_send,
    std::array<std::vector<int>, N_EDGE>& halo_recv,
    std::array<std::vector<int>, N_VERTEX>& corner_ids,
    std::array<std::vector<int>, N_VERTEX>& corner_send) const
{
    for (auto edge : edges) {
        if (((edge == LEFT || edge == RIGHT) && _px) || ((edge == TOP || edge == BOTTOM) && _py)) {

            for (auto it = _neighbours_p[edge].begin(); it != _neighbours_p[edge].end(); ++it) {
                ids[edge].push_back(it->first);
                halo_sizes[edge].push_back(it->second);
                halo_send[edge].push_back(_send_pos_p[edge].at(it->first));
                halo_recv[edge].push_back(_recv_pos_p[edge].at(it->first));
            }
        }
    }

    for (auto vertex : vertices) {
        for (auto it = _corner_neighbours_p[vertex].begin();
             it != _corner_neighbours_p[vertex].end(); ++it) {
            corner_ids[vertex].push_back(it->first);
            corner_send[vertex].push_back(_corner_send_pos_p[vertex].at(it->first));
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
    std::array<std::vector<int>, N_EDGE> ids, halos, halo_send, halo_recv;
    std::array<std::vector<int>, N_VERTEX> corner_ids, corner_send;
    get_neighbour_info(ids, halos, halo_send, halo_recv, corner_ids, corner_send);

    std::vector<int> num_neighbours(N_EDGE), dims(N_EDGE, 0), offsets(N_EDGE, 0);
    for (auto edge : edges) {
        num_neighbours[edge] = (int)ids[edge].size();
        CHECK_MPI(MPI_Allreduce(&num_neighbours[edge], &dims[edge], 1, MPI_INT, MPI_SUM, _comm));
        CHECK_MPI(MPI_Exscan(&num_neighbours[edge], &offsets[edge], 1, MPI_INT, MPI_SUM, _comm));
    }

    std::vector<int> num_corner_neighbours(N_VERTEX), corner_dims(N_VERTEX, 0),
        corner_offsets(N_VERTEX, 0);
    for (auto vertex : vertices) {
        num_corner_neighbours[vertex] = (int)corner_ids[vertex].size();
        CHECK_MPI(MPI_Allreduce(
            &num_corner_neighbours[vertex], &corner_dims[vertex], 1, MPI_INT, MPI_SUM, _comm));
        CHECK_MPI(MPI_Exscan(
            &num_corner_neighbours[vertex], &corner_offsets[vertex], 1, MPI_INT, MPI_SUM, _comm));
    }

    // Prepare periodic neighbour data
    std::array<std::vector<int>, N_EDGE> ids_p, halos_p, halo_send_p, halo_recv_p;
    std::array<std::vector<int>, N_VERTEX> corner_ids_p, corner_send_p;
    get_neighbour_info_periodic(
        ids_p, halos_p, halo_send_p, halo_recv_p, corner_ids_p, corner_send_p);
    std::vector<int> num_neighbours_p(N_EDGE), dims_p(N_EDGE, 0), offsets_p(N_EDGE, 0);
    for (auto edge : edges) {
        num_neighbours_p[edge] = (int)ids_p[edge].size();
        CHECK_MPI(
            MPI_Allreduce(&num_neighbours_p[edge], &dims_p[edge], 1, MPI_INT, MPI_SUM, _comm));
        CHECK_MPI(
            MPI_Exscan(&num_neighbours_p[edge], &offsets_p[edge], 1, MPI_INT, MPI_SUM, _comm));
    }

    std::vector<int> num_corner_neighbours_p(N_VERTEX), corner_dims_p(N_VERTEX, 0),
        corner_offsets_p(N_VERTEX, 0);
    for (auto vertex : vertices) {
        num_corner_neighbours_p[vertex] = (int)corner_ids_p[vertex].size();
        CHECK_MPI(MPI_Allreduce(
            &num_corner_neighbours_p[vertex], &corner_dims_p[vertex], 1, MPI_INT, MPI_SUM, _comm));
        CHECK_MPI(MPI_Exscan(&num_corner_neighbours_p[vertex], &corner_offsets_p[vertex], 1,
            MPI_INT, MPI_SUM, _comm));
    }

    // Define dimensions in netCDF file
    int dimid;
    std::vector<int> dimids(N_EDGE);
    NC_CHECK(nc_def_dim(nc_id, "P", _total_num_procs, &dimid));
    for (auto edge : edges) {
        NC_CHECK(nc_def_dim(nc_id, dir_chars[edge].c_str(), dims[edge], &dimids[edge]));
    }

    std::vector<int> corner_dimids(N_VERTEX);
    for (auto vertex : vertices) {
        NC_CHECK(nc_def_dim(
            nc_id, corner_dir_chars[vertex].c_str(), corner_dims[vertex], &corner_dimids[vertex]));
    }

    // Define periodic dimensions in netCDF file
    std::vector<int> dimids_p(N_EDGE);
    for (auto edge : edges) {
        NC_CHECK(nc_def_dim(
            nc_id, (dir_chars[edge] + "_periodic").c_str(), dims_p[edge], &dimids_p[edge]));
    }

    std::vector<int> corner_dimids_p(N_VERTEX);
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
    int halo_send_vid[N_EDGE];
    int halo_recv_vid[N_EDGE];
    for (auto edge : edges) {
        // Connectivity group
        NC_CHECK(nc_def_var(connectivity_gid, (dir_names[edge] + "_neighbours").c_str(), NC_INT, 1,
            &dimid, &num_vid[edge]));
        NC_CHECK(nc_def_var(connectivity_gid, (dir_names[edge] + "_neighbour_ids").c_str(), NC_INT,
            1, &dimids[edge], &ids_vid[edge]));
        NC_CHECK(nc_def_var(connectivity_gid, (dir_names[edge] + "_neighbour_halos").c_str(),
            NC_INT, 1, &dimids[edge], &halos_vid[edge]));
        NC_CHECK(nc_def_var(connectivity_gid, (dir_names[edge] + "_neighbour_halo_send").c_str(),
            NC_INT, 1, &dimids[edge], &halo_send_vid[edge]));
        NC_CHECK(nc_def_var(connectivity_gid, (dir_names[edge] + "_neighbour_halo_recv").c_str(),
            NC_INT, 1, &dimids[edge], &halo_recv_vid[edge]));
    }

    int num_corner_vid[N_VERTEX];
    int ids_corner_vid[N_VERTEX];
    int corner_send_vid[N_VERTEX];
    for (auto vertex : vertices) {
        // Connectivity group
        NC_CHECK(nc_def_var(connectivity_gid, (corner_dir_names[vertex] + "_neighbours").c_str(),
            NC_INT, 1, &dimid, &num_corner_vid[vertex]));
        NC_CHECK(nc_def_var(connectivity_gid, (corner_dir_names[vertex] + "_neighbour_ids").c_str(),
            NC_INT, 1, &corner_dimids[vertex], &ids_corner_vid[vertex]));
        NC_CHECK(
            nc_def_var(connectivity_gid, (corner_dir_names[vertex] + "_neighbour_send").c_str(),
                NC_INT, 1, &corner_dimids[vertex], &corner_send_vid[vertex]));
    }

    int num_vid_p[N_EDGE];
    int ids_vid_p[N_EDGE];
    int halos_vid_p[N_EDGE];
    int halo_send_vid_p[N_EDGE];
    int halo_recv_vid_p[N_EDGE];
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
            &halo_send_vid_p[edge]));
        NC_CHECK(nc_def_var(connectivity_gid,
            (dir_names[edge] + "_neighbour_halo_recv_periodic").c_str(), NC_INT, 1, &dimids_p[edge],
            &halo_recv_vid_p[edge]));
    }

    int num_corner_vid_p[N_VERTEX];
    int ids_corner_vid_p[N_VERTEX];
    int corner_send_vid_p[N_VERTEX];
    for (auto vertex : vertices) {
        // Connectivity group
        NC_CHECK(nc_def_var(connectivity_gid,
            (corner_dir_names[vertex] + "_neighbours_periodic").c_str(), NC_INT, 1, &dimid,
            &num_corner_vid_p[vertex]));
        NC_CHECK(nc_def_var(connectivity_gid,
            (corner_dir_names[vertex] + "_neighbour_ids_periodic").c_str(), NC_INT, 1,
            &corner_dimids_p[vertex], &ids_corner_vid_p[vertex]));
        NC_CHECK(nc_def_var(connectivity_gid,
            (corner_dir_names[vertex] + "_neighbour_send_periodic").c_str(), NC_INT, 1,
            &corner_dimids_p[vertex], &corner_send_vid_p[vertex]));
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
        NC_CHECK(nc_var_par_access(connectivity_gid, halo_send_vid[edge], NC_COLLECTIVE));
        NC_CHECK(nc_put_vara_int(
            connectivity_gid, halo_send_vid[edge], &start, &count, halo_send[edge].data()));
        NC_CHECK(nc_put_vara_int(
            connectivity_gid, halo_recv_vid[edge], &start, &count, halo_recv[edge].data()));
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
            connectivity_gid, halo_send_vid_p[edge], &start, &count, halo_send_p[edge].data()));
        NC_CHECK(nc_put_vara_int(
            connectivity_gid, halo_recv_vid_p[edge], &start, &count, halo_recv_p[edge].data()));
    }

    for (auto vertex : vertices) {
        // Numbers of "corner" neighbours
        size_t start = _rank;
        NC_CHECK(nc_var_par_access(connectivity_gid, num_corner_vid[vertex], NC_COLLECTIVE));
        NC_CHECK(nc_put_var1_int(
            connectivity_gid, num_corner_vid[vertex], &start, &num_corner_neighbours[vertex]));
        // Numbers of corner neighbours for periodic dimensions
        NC_CHECK(nc_var_par_access(connectivity_gid, num_corner_vid_p[vertex], NC_COLLECTIVE));
        NC_CHECK(nc_put_var1_int(
            connectivity_gid, num_corner_vid_p[vertex], &start, &num_corner_neighbours_p[vertex]));

        // "Corner" IDs and halos
        start = corner_offsets[vertex];
        size_t count = num_corner_neighbours[vertex];
        NC_CHECK(nc_var_par_access(connectivity_gid, ids_corner_vid[vertex], NC_COLLECTIVE));
        NC_CHECK(nc_put_vara_int(
            connectivity_gid, ids_corner_vid[vertex], &start, &count, corner_ids[vertex].data()));
        NC_CHECK(nc_var_par_access(connectivity_gid, corner_send_vid[vertex], NC_COLLECTIVE));
        NC_CHECK(nc_put_vara_int(
            connectivity_gid, corner_send_vid[vertex], &start, &count, corner_send[vertex].data()));

        // "Corner" IDs and halos for periodic dimensions
        start = corner_offsets_p[vertex];
        count = num_corner_neighbours_p[vertex];
        NC_CHECK(nc_var_par_access(connectivity_gid, ids_corner_vid_p[vertex], NC_COLLECTIVE));
        NC_CHECK(nc_put_vara_int(connectivity_gid, ids_corner_vid_p[vertex], &start, &count,
            corner_ids_p[vertex].data()));
        NC_CHECK(nc_var_par_access(connectivity_gid, halos_vid_p[vertex], NC_COLLECTIVE));
        NC_CHECK(nc_put_vara_int(connectivity_gid, corner_send_vid_p[vertex], &start, &count,
            corner_send_p[vertex].data()));
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

            // check edge neighours
            for (auto edge : edges) {
                if (is_neighbour(domains[_rank], domains[p], edge)) {
                    int halo_size = domain_overlap(domains[_rank], domains[p], edge);
                    if (halo_size > 0) {
                        _neighbours[edge].insert(std::pair<int, int>(p, halo_size));
                        int sendPos = 0;
                        int recvPos = 0;
                        haloEdgeBufferPositions(domains[_rank], domains[p], edge, sendPos, recvPos);
                        _send_pos[edge].insert(std::pair<int, int>(p, sendPos));
                        _recv_pos[edge].insert(std::pair<int, int>(p, recvPos));
                    }
                }
            }

            // check corner neighours
            for (auto vertex : vertices) {
                if (is_corner_neighbour(domains[_rank], domains[p], vertex)) {
                    _corner_neighbours[vertex].insert(std::pair<int, int>(p, 1));
                    int sendPos = 0;
                    haloCornerBufferPositions(domains[_rank], domains[p], vertex, sendPos);
                    _corner_send_pos[vertex].insert(std::pair<int, int>(p, sendPos));
                }
            }
        }

        // When finding neighours *across periodic boundaries*, we need to check against the
        // current rank, too, because a subdomain can be a periodic neighbour of itself.
        // check periodic edge neighours
        for (auto edge : edges) {
            if (is_neighbour(domains[_rank], domains[p], edge, _px, _py)) {
                int halo_size = domain_overlap(domains[_rank], domains[p], edge);
                if (halo_size > 0) {
                    _neighbours_p[edge].insert(std::pair<int, int>(p, halo_size));
                    int sendPos = 0;
                    int recvPos = 0;
                    haloEdgeBufferPositions(domains[_rank], domains[p], edge, sendPos, recvPos);
                    _send_pos_p[edge].insert(std::pair<int, int>(p, sendPos));
                    _recv_pos_p[edge].insert(std::pair<int, int>(p, recvPos));
                }
            }
        }

        // check periodic corner neighours
        for (auto vertex : vertices) {
            if (is_corner_neighbour(domains[_rank], domains[p], vertex, _px, _py)) {
                if (_corner_neighbours[vertex].size() > 0) {
                    // skip if we have already counted as a non-periodic neighbour
                    continue;
                }
                _corner_neighbours_p[vertex].insert(std::pair<int, int>(p, 1));
                int sendPos = 0;
                haloCornerBufferPositions(domains[_rank], domains[p], vertex, sendPos);
                _corner_send_pos_p[vertex].insert(std::pair<int, int>(p, sendPos));
            }
        }
    }
}
