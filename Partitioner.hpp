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
     * @brief Initializes the partitioner with grid parameters.
     *
     * Sets the method variables (grid extents, bounding box, periodicity) from
     * the provided grid.
     *
     * @param grid Reference to the grid object.
     */
    virtual void initialize(Grid& grid) = 0;

    /*!
     * @brief Returns the new bounding box for this process after partitioning.
     *
     * @param global0 Global coordinate in the 1st dimension of the upper left
     * corner.
     * @param global1 Global coordinate in the 2nd dimension of the upper left
     * corner.
     * @param localExt0 Local extent in the 1st dimension of the grid.
     * @param localExt1 Local extent in the 2nd dimension of the grid.
     */
    void getBoundingBox(int& global0, int& global1, int& localExt0, int& localExt1) const;

    /*!
     * @brief Returns vectors containing the MPI ranks, halo sizes and halo starting indices of the
     * neighbours of this process in the domain interior after partitioning. The neighbours are
     * ordered left, right, bottom, top.
     *
     * @param ids MPI ranks of the neighbours for each direction
     * @param corner_ids MPI ranks of the "corner" neighbours for each direction
     * @param haloSizes Halo sizes of the neighbours for each direction
     * @param haloSend index in send buffer to get halo data
     * @param haloRecv index in recv buffer to put halo data
     */
    void getNeighbourInfo(std::array<std::vector<int>, N_EDGE>& ids,
        std::array<std::vector<int>, N_EDGE>& haloSizes,
        std::array<std::vector<int>, N_EDGE>& haloSend,
        std::array<std::vector<int>, N_EDGE>& haloRecv,
        std::array<std::vector<int>, N_CORNER>& cornerIds,
        std::array<std::vector<int>, N_CORNER>& cornerSend) const;
    /*!
     * @brief Returns vectors containing the MPI ranks, halo sizes and halo starting indices of
     * the neighbours of this process across periodic boundaries after partitioning. The
     * neighbours are ordered left, right, bottom, top.
     *
     * @param ids MPI ranks of the periodic neighbours for each direction
     * @param haloSizes Halo sizes of the periodic neighbours for each direction
     * @param haloSend index in send buffer to get halo data
     * @param haloRecv index in recv buffer to put halo data
     */
    void getNeighbourInfoPeriodic(std::array<std::vector<int>, N_EDGE>& ids,
        std::array<std::vector<int>, N_EDGE>& haloSizes,
        std::array<std::vector<int>, N_EDGE>& haloSend,
        std::array<std::vector<int>, N_EDGE>& haloRecv,
        std::array<std::vector<int>, N_CORNER>& cornerIds,
        std::array<std::vector<int>, N_CORNER>& cornerSend) const;

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
    void saveMask(const std::string& filename) const;

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
    void saveMetadata(const std::string& filename) const;

    /*!
     * @brief Sets the total number of processes.
     *
     * @param totalNumProcs Total number of processes in communicator.
     */
    void setTotalNumProcs(int totalNumProcs);

    /*!
     * @brief Sets the MPI communicator.
     *
     * @param comm MPI communicator.
     */
    void setComm(MPI_Comm comm);

    /*!
     * @brief Returns the local extents in each dimension after partitioning.
     *
     * @return A vector of size NDIMS containing the local extents.
     */
    std::vector<int> getLocalExtNew() const;

    /*!
     * @brief Sets the local extents in each dimension after partitioning.
     *
     * @param localExtNew A vector of size NDIMS containing the local extents.
     */
    void setLocalExtNew(const std::vector<int>& localExtNew);

    /*!
     * @brief Returns the global coordinates of the upper left corner after partitioning.
     *
     * @return A vector of size NDIMS containing the global coordinates.
     */
    std::vector<int> getGlobalNew() const;

    /*!
     * @brief Sets the global coordinates of the upper left corner after partitioning.
     *
     * @param globalNew A vector of size NDIMS containing the global coordinates.
     */
    void setGlobalNew(const std::vector<int>& globalNew);

    /*!
     * @brief Returns the global coordinates of the upper left corner.
     *
     * @return A vector of size NDIMS containing the global coordinates.
     */
    std::vector<int> getGlobal() const;

    /*!
     * @brief Sets the global coordinates of the upper left corner.
     *
     * @param global A vector of size NDIMS containing the global coordinates.
     */
    void setGlobal(const std::vector<int>& global);

    /*!
     * @brief Returns the global extents in each dimension.
     *
     * @return A vector of size NDIMS containing the global extents.
     */
    std::vector<int> getGlobalExt() const;

    /*!
     * @brief Sets the global extents in each dimension.
     *
     * @param globalExt A vector of size NDIMS containing the global extents.
     */
    void setGlobalExt(const std::vector<int>& globalExt);

    /*!     * @brief Returns the process IDs of the latest 2D domain decomposition.
     *
     * @return A vector containing the partition ID of each point in the grid.
     */
    std::vector<int> getProcId() const;

    /*!     * @brief Sets the process IDs of the latest 2D domain decomposition.
     *
     * @param procId A vector containing the partition ID of each point in the grid.
     */
    void setProcId(const std::vector<int>& procId);

    // Discover the neighbours and halo sizes of the process after partitioning
    void discover_neighbours();

protected:
    // Construct a partitioner
    // We are using the named constructor idiom so that objects can only be
    // created in the heap to ensure it's dtor is executed before MPI_Finalize()
    Partitioner(MPI_Comm comm);

protected:
    MPI_Comm _comm; // MPI communicator
    int _rank = -1; // Process rank
    int _totalNumProcs = -1; // Total number of processes in communicator
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
    std::vector<std::string> globalExtentNames = { "NX", "NY" };

    // Total number of processes in each dimension
    std::vector<int> _numProcs = std::vector<int>(NDIMS, -1);

    // Global extents in each dimension
    std::vector<int> _globalExt = std::vector<int>(NDIMS, 0);

    // Local extents in each dimension
    std::vector<int> _localExt = std::vector<int>(NDIMS, 0);

    // Global coordinates of upper left corner
    std::vector<int> _global = std::vector<int>(NDIMS, -1);

    // Local extents in each dimension (after partitioning)
    std::vector<int> _localExtNew = std::vector<int>(NDIMS, 0);

    // Global coordinates of upper left corner (after partitioning)
    std::vector<int> _globalNew = std::vector<int>(NDIMS, -1);

    // Process ids of partition (dense form)
    std::vector<int> _procId = {};

    // Vector of maps of neighbours to their halo sizes after partitioning
    std::vector<std::map<int, int>> _neighbours = std::vector<std::map<int, int>>(NNBRS);

    // Vector of maps of "corner" neighbours to their halo sizes after partitioning (corners are all
    // of size 1)
    std::vector<std::map<int, int>> _cornerNeighbours = std::vector<std::map<int, int>>(NNBRS);

    // Vector of maps of neighbours to their send buffer indices - index of data to fetch from send
    // buffer
    std::vector<std::map<int, int>> _sendPos = std::vector<std::map<int, int>>(NNBRS);

    // Vector of maps of neighbours to their recv (receive) buffer indices - index where data will
    // be stored in the recv buffer
    std::vector<std::map<int, int>> _recvPos = std::vector<std::map<int, int>>(NNBRS);

    // Vector of maps of "corner" neighbours to their halo start indices after partitioning
    std::vector<std::map<int, int>> _cornerSendPos = std::vector<std::map<int, int>>(NNBRS);

    // Vector of maps of periodic neighbours to their halo sizes after partitioning
    std::vector<std::map<int, int>> _neighbours_p = std::vector<std::map<int, int>>(NNBRS);

    // Vector of maps of periodic neighbours to their send buffer indices - index of data to fetch
    // from send buffer
    std::vector<std::map<int, int>> _sendPos_p = std::vector<std::map<int, int>>(NNBRS);

    // Vector of maps of periodic neighbours to their recv (receive) buffer indices - index where
    // data will be stored in the recv buffer
    std::vector<std::map<int, int>> _recvPos_p = std::vector<std::map<int, int>>(NNBRS);

    // Vector of maps of "corner" neighbours to their halo sizes after partitioning
    std::vector<std::map<int, int>> _cornerNeighbours_p = std::vector<std::map<int, int>>(NNBRS);

    // Vector of maps of "corner" neighbours to their halo start indices after partitioning
    std::vector<std::map<int, int>> _cornerSendPos_p = std::vector<std::map<int, int>>(NNBRS);

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

    /*!
     * @brief Check if two domains are neighbouring. If true, then domain 2 is the [edge] neighbour
     * of domain 1, relative to domain 1. e.g., if domain 2 is to the right of domain 1, the
     * following call will return true:
     * isNeighbour(d1, d2, RIGHT)
     *
     * @param d1 first domain
     * @param d2 second domain
     * @param edge LEFT, RIGHT, BOTTOM or TOP
     * @param isPx are we looking for periodic neighbour in x-direction?
     * @param isPy are we looking for periodic neighbour in y-direction?
     * @return bool
     */
    bool isNeighbour(const Domain d1, const Domain d2, const Edge edge, const bool isPx = false,
        const bool isPy = false);
    bool isCornerNeighbour(const Domain d1, const Domain d2, const Corner corner,
        const bool isPx = false, const bool isPy = false);

    /*!
     * @brief Compute the sendPos and recvPos for halo exchange.
     *
     * See https://nextsim-dg.readthedocs.io/en/latest/halo-exchange.html for a diagram explaining
     * how the start locations are calculated.
     *
     * @param d1 first domain
     * @param d2 second domain
     * @param edge LEFT, RIGHT, BOTTOM or TOP
     * @param sendPos position to get data from send buffer
     * @param recvPos position to store data in recv buffer
     */
    void haloEdgeBufferPositions(
        const Domain d1, const Domain d2, const Edge edge, int& sendPos, int& recvPos);

    /*!
     * @brief Compute the sendPos for corner neighbours
     *
     * @param d1 first domain
     * @param d2 second domain
     * @param edge TOP_LEFT, TOP_RIGHT, BOTTOM_RIGHT, BOTTOM_LEFT
     * @param sendPos position to get data from send buffer
     */
    void haloCornerBufferPositions(
        const Domain d1, const Domain d2, const Corner corner, int& sendPos);
};
