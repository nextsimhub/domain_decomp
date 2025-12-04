/*!
 * @file Partitioner.hpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @date 05 Nov 2024
 */

#pragma once

#include <map>

#include "DomainUtils.hpp"
#include "Grid.hpp"
#include "domain_decomp_export.hpp"

/*!
 * Supported partitioners.
 */
enum class LIB_EXPORT PartitionerType {
    Zoltan_RCB /*!< Recursive Coordinate Bisection (RCB) geometric partitioning
                  algorithm from the Zoltan toolkit */
};

/*!
 * @class Partitioner
 * @brief Abstract polymorphic class that encapsulates a 2D grid partitioner.
 */

class LIB_EXPORT Partitioner {
public:

    // Disallow compiler-generated special functions
    Partitioner(const Partitioner&) = delete;
    Partitioner& operator=(const Partitioner&) = delete;

    /*!
     * @brief Destructor.
     */
    virtual ~Partitioner() {};

    /*!
     * @brief Partitions a 2D grid into rectangular boxes, one per process.
     *
     * Partitions a 2D grid into rectangular boxes, one per process, taking into
     * account a land mask, if provided. The Grid object is updated with the new
     * partitioning information.
     */
    virtual void partition(Grid& grid) = 0;

    /*!
     * @brief Returns the new bounding box for this process after partitioning.
     *
     * @param global_0 Global coordinate in the 1st dimension of the upper left
     * corner.
     * @param global_1 Global coordinate in the 2nd dimension of the upper left
     * corner.
     * @param local_ext_0 Local extent in the 1st dimension of the grid.
     * @param local_ext_1 Local extent in the 2nd dimension of the grid.
     */
    void get_bounding_box(int& global_0, int& global_1, int& local_ext_0, int& local_ext_1) const;

    /*!
     * @brief Returns vectors containing the MPI ranks, halo sizes and halo starting indices of the
     * neighbours of this process in the domain interior after partitioning. The neighbours are
     * ordered left, right, bottom, top.
     *
     * @param ids MPI ranks of the neighbours for each direction
     * @param corner_ids MPI ranks of the "corner" neighbours for each direction
     * @param halo_sizes Halo sizes of the neighbours for each direction
     * @param halo_send index in send buffer to get halo data
     * @param halo_recv index in recv buffer to put halo data
     */
    void get_neighbour_info(std::array<std::vector<int>, N_EDGE>& ids,
        std::array<std::vector<int>, N_EDGE>& halo_sizes,
        std::array<std::vector<int>, N_EDGE>& halo_send,
        std::array<std::vector<int>, N_EDGE>& halo_recv,
        std::array<std::vector<int>, N_VERTEX>& corner_ids,
        std::array<std::vector<int>, N_VERTEX>& corner_send) const;
    /*!
     * @brief Returns vectors containing the MPI ranks, halo sizes and halo starting indices of
     * the neighbours of this process across periodic boundaries after partitioning. The
     * neighbours are ordered left, right, bottom, top.
     *
     * @param ids MPI ranks of the periodic neighbours for each direction
     * @param halo_sizes Halo sizes of the periodic neighbours for each direction
     * @param halo_send index in send buffer to get halo data
     * @param halo_recv index in recv buffer to put halo data
     */
    void get_neighbour_info_periodic(std::array<std::vector<int>, N_EDGE>& ids,
        std::array<std::vector<int>, N_EDGE>& halo_sizes,
        std::array<std::vector<int>, N_EDGE>& halo_send,
        std::array<std::vector<int>, N_EDGE>& halo_recv,
        std::array<std::vector<int>, N_VERTEX>& corner_ids,
        std::array<std::vector<int>, N_VERTEX>& corner_send) const;

    /*!
     * @brief Saves the partition IDs of the latest 2D domain decomposition in a
     * NetCDF file.
     *
     * Saves the partition IDs of the latest 2D domain decomposition in a NetCDF
     * file. The NetCDF file contains dimensions x and y and integer variable
     * pid(x, y) which defines the partition ID of each point in the grid.
     *
     * @param filename Name of the NetCDF file.
     */
    void save_mask(const std::string& filename) const;

    /*!
     * @brief Saves the boxes and connectivity information of the latest domain
     * decomposition in a NetCDF file.
     *
     * Saves the boxes and connectivity information of the latest 2D domain
     * decomposition in a NetCDF file. The NetCDF file contains a dimension P
     * equal to the number of partitions and integer variables domain_x(P),
     * domain_y(P), domain_extent_x(P) and domain_extent_y(P). Variables domain_x
     * and domain_y are defined as the coordinates of the upper left corner of the
     * box for each partition, while the domain_extent_x and domain_extent_y
     * variables define the local extent of the x and y dimensions respectively.
     * The file also defines the variables X_neighbours(P), X_neighbour_ids(X_dim)
     * and X_neighbour_halos(X_dim), where X is top/bottom/left/right, which
     * correspond to the number of neighbours per process, the neighbour IDs and
     * halo sizes of each process sorted from lower to higher MPI rank.
     *
     * @param filename Name of the NetCDF file.
     */
    void save_metadata(const std::string& filename) const;

protected:
    // Construct a partitioner
    // We are using the named constructor idiom so that objects can only be
    // created in the heap to ensure it's dtor is executed before MPI_Finalize()
    Partitioner(MPI_Comm comm);

    // Discover the neighbours and halo sizes of the process after partitioning
    void discover_neighbours();

protected:
    MPI_Comm _comm; // MPI communicator
    int _rank = -1; // Process rank
    int _total_num_procs = -1; // Total number of processes in communicator
    static const int NDIMS = 2; // Number of dimensions
    static const int NNBRS = 2 * NDIMS; // Number of neighbours (two per dimension)
    bool _px = false; // Periodic boundary in the x-direction
    bool _py = false; // Periodic boundary in the y-direction

    // Letters used for each dimension
    std::vector<std::string> dim_chars = { "x", "y" };

    // Letters used for each direction
    std::vector<std::string> dir_chars = { "L", "R", "B", "T" };

    // Letters used for each "corner"
    std::vector<std::string> corner_dir_chars = { "TL", "TR", "BR", "BL" };

    // Names used for each direction
    std::vector<std::string> dir_names = { "left", "right", "bottom", "top" };

    // Names used for each "corner"
    std::vector<std::string> corner_dir_names
        = { "top_left", "top_right", "bottom_right", "bottom_left" };

    // Names used for global dimension extents
    std::vector<std::string> global_extent_names = { "NX", "NY" };

    // Total number of processes in each dimension
    std::vector<int> _num_procs = std::vector<int>(NDIMS, -1);

    // Global extents in each dimension
    std::vector<int> _global_ext = std::vector<int>(NDIMS, 0);

    // Local extents in each dimension
    std::vector<int> _local_ext = std::vector<int>(NDIMS, 0);

    // Global coordinates of upper left corner
    std::vector<int> _global = std::vector<int>(NDIMS, -1);

    // Local extents in each dimension (after partitioning)
    std::vector<int> _local_ext_new = std::vector<int>(NDIMS, 0);

    // Global coordinates of upper left corner (after partitioning)
    std::vector<int> _global_new = std::vector<int>(NDIMS, -1);

    // Process ids of partition (dense form)
    std::vector<int> _proc_id = {};

    // Vector of maps of neighbours to their halo sizes after partitioning
    std::vector<std::map<int, int>> _neighbours = std::vector<std::map<int, int>>(NNBRS);

    // Vector of maps of "corner" neighbours to their halo sizes after partitioning (corners are all
    // of size 1)
    std::vector<std::map<int, int>> _corner_neighbours = std::vector<std::map<int, int>>(NNBRS);

    // Vector of maps of neighbours to their send buffer indices - index of data to fetch from send
    // buffer
    std::vector<std::map<int, int>> _send_pos = std::vector<std::map<int, int>>(NNBRS);

    // Vector of maps of neighbours to their recv (receive) buffer indices - index where data will
    // be stored in the recv buffer
    std::vector<std::map<int, int>> _recv_pos = std::vector<std::map<int, int>>(NNBRS);

    // Vector of maps of "corner" neighbours to their halo start indices after partitioning
    std::vector<std::map<int, int>> _corner_send_pos = std::vector<std::map<int, int>>(NNBRS);

    // Vector of maps of periodic neighbours to their halo sizes after partitioning
    std::vector<std::map<int, int>> _neighbours_p = std::vector<std::map<int, int>>(NNBRS);

    // Vector of maps of periodic neighbours to their send buffer indices - index of data to fetch
    // from send buffer
    std::vector<std::map<int, int>> _send_pos_p = std::vector<std::map<int, int>>(NNBRS);

    // Vector of maps of periodic neighbours to their recv (receive) buffer indices - index where
    // data will be stored in the recv buffer
    std::vector<std::map<int, int>> _recv_pos_p = std::vector<std::map<int, int>>(NNBRS);

    // Vector of maps of "corner" neighbours to their halo sizes after partitioning
    std::vector<std::map<int, int>> _corner_neighbours_p = std::vector<std::map<int, int>>(NNBRS);

    // Vector of maps of "corner" neighbours to their halo start indices after partitioning
    std::vector<std::map<int, int>> _corner_send_pos_p = std::vector<std::map<int, int>>(NNBRS);

public:
    struct LIB_EXPORT Factory {
        /*!
         * @brief Factory function for creating grid partitioners.
         *
         * @param comm MPI communicator.
         * @param argc The number of arguments.
         * @param argv The argument vector.
         * @param type Type of partitioner.
         * @return A Partitioner object.
         */
        static Partitioner* create(MPI_Comm comm, int argc, char** argv, PartitionerType type);
    };

    bool is_neighbour(const Domain d1, const Domain d2, const Edge edge, const bool is_px = false,
        const bool is_py = false);
    bool is_corner_neighbour(const Domain d1, const Domain d2, const Vertex vertex,
        const bool is_px = false, const bool is_py = false);

    /*!
     * @brief Check if two domains are neighbouring. If true, then domain 2 is the [edge] neighbour
     * of domain 1, relative to domain 1. e.g., if domain 2 is to the right of domain 1, the
     * following call will return true:
     * is_neighbour(d1, d2, RIGHT)
     *
     * @param d1 first domain
     * @param d2 second domain
     * @param edge LEFT, RIGHT, BOTTOM or TOP
     * @param is_px are we looking for periodic neighbour in x-direction?
     * @param is_py are we looking for periodic neighbour in y-direction?
     * @return bool
     */
   
    /*!
     * @brief Compute the start location of the halo for a given pair of neighbouring domains.
     *
     * TODO: This needs to be updated to reflect changes to halo start
     *
     * For example, in the diagram below. If we want to compute the start location for the halo of
     * domain 1 which is the LEFT neighbour of domain 2, halo_start should return 4 (see square
     * bracket in diagram below).
     *         memory             domains          halo required
     *   4┌───────────┬─────┐ ┌──────────┬────┐  ┌──────────┬────┐
     *    │ 5 6 7 8 9 │ 6 7 │ │          │    │  │         9│    │
     *    │           │     │ │    1     │    │  │          │    │
     *    │ 0 1 2 3[4]│ 4 5 │ │          │    │  │         4│    │
     *   2├───────────┤     │ ├──────────┤  2 │  ├──────────┤    │
     *    │ 5 6 7 8 9 │ 2 3 │ │          │    │  │          │    │
     *    │           │     │ │    0     │    │  │          │    │
     * Y▲ │ 0 1 2 3 4 │ 0 1 │ │          │    │  │          │    │
     *  │ └───────────┴─────┘ └──────────┴────┘  └──────────┴────┘
     *  │0            5     7
     *  └───►X
     *
     * Another example. If we want the start location of the halo for domain 0 which is domain 1's
     * BOTTOM neighbour. We want halo_start to return 5. (see square bracket in diagram below).
     *         memory             domains          halo required
     *   4┌───────────┬─────┐ ┌──────────┬────┐  ┌──────────┬────┐
     *    │ 5 6 7 8 9 │ 6 7 │ │          │    │  │          │    │
     *    │           │     │ │    1     │    │  │          │    │
     *    │ 0 1 2 3 4 │ 4 5 │ │          │    │  │          │    │
     *   2├───────────┤     │ ├──────────┤  2 │  ├──────────┤    │
     *    │[5]6 7 8 9 │ 2 3 │ │          │    │  │5 6 7 8 9 │    │
     *    │           │     │ │    0     │    │  │          │    │
     * Y▲ │ 0 1 2 3 4 │ 0 1 │ │          │    │  │          │    │
     *  │ └───────────┴─────┘ └──────────┴────┘  └──────────┴────┘
     *  │0            5     7
     *  └───►X
     *
     * @param d1 first domain
     * @param d2 second domain
     * @param edge LEFT, RIGHT, BOTTOM or TOP
     * @return starting index of halo for the flattened domain array
     */
    void haloBufferPositions(
        const Domain d1, const Domain d2, const Edge edge, int& sendPos, int& recvPos);
    int halo_corner_start(const Domain d1, const Domain d2, const Vertex vertex,
        const bool is_px = false, const bool is_py = false);
};
