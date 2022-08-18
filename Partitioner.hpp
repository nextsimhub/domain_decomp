/*!
 * @file Partitioner.hpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @date 25 June 2022
 */

#pragma once

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
class LIB_EXPORT Partitioner
{
public:
  // Disallow compiler-generated special functions
  Partitioner(const Partitioner&) = delete;
  Partitioner& operator=(const Partitioner&) = delete;

  /*!
   * @brief Destructor.
   */
  virtual ~Partitioner(){};

  /*!
   * @brief Partitions a 2D grid into rectangular boxes, one per process.
   *
   * Partitions a 2D grid into rectangular boxes, one per process, taking into
   * account a land mask, if provided. The Grid object is updated with the new
   * partitioning information.
   */
  virtual void partition(Grid& grid) = 0;

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
   * @brief Saves the boxes of the latest domain decomposition in a NetCDF file.
   *
   * Saves the boxes of the latest 2D domain decomposition in a NetCDF file. The
   * NetCDF file contains a dimension P equal to the number of partitions and
   * integer variables global_x(P), global_y(P), local_extent_x(P) and
   * local_extent_y(P). Variables global_x and global_y are defined as the
   * coordinates of the upper left corner of the box for each partition, while
   * the local_extent_x and local_extent_y variables define the local extent of
   * the x and y dimensions respectively.
   *
   * @param filename Name of the NetCDF file.
   */
  void save_metadata(const std::string& filename) const;

protected:
  // Construct a partitioner
  // We are using the named constructor idiom so that objects can only be
  // created in the heap to ensure it's dtor is executed before MPI_Finalize()
  Partitioner(MPI_Comm comm, int argc, char** argv);

protected:
  MPI_Comm _comm;        // MPI communicator
  int _rank = -1;        // Process rank
  int _num_procs = -1;   // Total number of processes in communicator
  int _num_procs_0 = -1; // Total number of processes in first dimension
  int _num_procs_1 = -1; // Total number of processes in second dimension
  int _global_ext_0 = 0; // Global extent in first extension
  int _global_ext_1 = 0; // Global extent in second extension
  int _local_ext_0 = 0;  // Local extent in first dimension (original)
  int _local_ext_1 = 0;  // Local extent in second dimension (original)
  int _global_0 = -1;    // Global coordinate in first dimension of upper left
                         // corner (original)
  int _global_1 = -1;    // Global coordinate in second dimension of upper left
                         // corner (original)
  int _local_ext_0_new
      = 0; /* Local extent in first dimension (after partitioning) */
  int _local_ext_1_new
      = 0; /* Local extent in second dimension (after partitioning) */
  int _global_0_new = -1; /* Global coordinate in first dimension of upper left
                             corner (after partitioning) */
  int _global_1_new = -1; /* Global coordinate in second dimension of upper left
                             corner (after partitioning) */
  std::vector<int> _proc_id = {}; // Process ids of partition (dense form)

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
    static Partitioner* create(MPI_Comm comm, int argc, char** argv,
                               PartitionerType type);
  };
};
