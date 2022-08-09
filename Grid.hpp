/*!
 * @file Grid.hpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @date 1 May 2022
 */

#pragma once

#include <string>
#include <vector>

#include <mpi.h>

/*!
 * @class Grid
 * @brief A class that encapsulates a distributed 2D grid describing the ocean.
 *
 * A class that encapsulates a distributed 2D grid of dimensions x and y, where
 * x is the first dimension and y is the second, with the last dimension varying
 * the fastest in terms of storage. The grid is partitioned evenly among
 * processes using a 2D decomposition, ignoring any land mask. The grid can be
 * subsequently re-partitioned differently using a Partitioner.
 */
class Grid
{
  // Define grid metadata
  const std::string x_dim_id = "x";
  const std::string y_dim_id = "y";
  const std::string data_id = "data";
  const std::string mask_id = "mask";

public:
  // Disallow compiler-generated special functions
  Grid(const Grid&) = delete;
  Grid& operator=(const Grid&) = delete;

  /*!
   * @brief Destructor.
   */
  ~Grid(){};

  /*!
   * @brief Constructs a distributed grid from a NetCDF file describing the
   * global domain.
   *
   * Constructs a distributed grid from a NetCDF file. We follow the convention
   * that the NetCDF grid file defines x as the first dimension and y as the
   * second and all arrays are dimensioned (nx, ny) where nx and ny are the grid
   * extents in first and second dimension. Data will be read with the last
   * dimension varying the fastest in terms of ordering (NetCDF C/C++ interface
   * convention).
   *
   * @param comm MPI communicator.
   * @param argc The number of arguments.
   * @param argv The argument vector.
   * @param filename Grid file in NetCDF format.
   * @return A distributed 2D grid object partitioned evenly in terms of grid
   * points.
   */
  // We are using the named constructor idiom so that objects can only be
  // created in the heap to ensure it's dtor is executed before MPI_Finalize()
  static Grid* create(MPI_Comm comm, int argc, char** argv,
                      const std::string& filename);

  /*!
   * @brief Returns the total number of objects in the local domain.
   *
   * @return Total number of objects in the local domain.
   */
  int get_num_objects() const;

  /*!
   * @brief Returns the number of non-land objects in the local domain.
   *
   * @return Number of non-land objects in the local domain.
   */
  int get_num_nonzero_objects() const;

  /*!
   * @brief Returns the global extent in x dimension.
   *
   * @return Global extent in x dimension.
   */
  int get_global_dim_x() const;

  /*!
   * @brief Returns the global extent in y dimension.
   *
   * @return Global extent in y dimension.
   */
  int get_global_dim_y() const;

  /*!
   * @brief Returns the local extent in x dimension.
   *
   * @return Local extent in x dimension.
   */
  int get_local_dim_x() const;

  /*!
   * @brief Returns the local extent in y dimension.
   *
   * @return Local extent in y dimension.
   */
  int get_local_dim_y() const;

  /*!
   * @brief Returns the global x coordinate of the top left object in this
   * process's partition.
   *
   * @return Global x coordinate of the top left object in this process's
   * partition.
   */
  int get_global_top_x() const;

  /*!
   * @brief Returns the global y coordinate of the top left object in this
   * process's partition.
   *
   * @return Global y coordinate of the top left object in this process's
   * partition.
   */
  int get_global_top_y() const;

  /*!
   * @brief Returns the number of processes in the x dimension.
   *
   * @return Number of processes in the x dimension.
   */
  int get_num_procs_x() const;

  /*!
   * @brief Returns the number of processes in the x dimension.
   *
   * @return Number of processes in the x dimension.
   */
  int get_num_procs_y() const;

  /*!
   * @brief Returns the global land mask, where x is the first dimension and y
   * is the second, with y varying the fastest in terms of storage.
   *
   * @return Global land mask.
   */
  const int* get_land_mask() const;

  /*!
   * @brief Returns the index mapping of sparse to dense representation, where x
   * is the first dimension and y is the second, with y varying the fastest in
   * terms of storage.
   *
   * @return Index mapping of sparse to dense representation.
   */
  const int* get_sparse_to_dense() const;

  /*!
   * @brief Returns the IDs of the non-land objects in the local domain, where x
   * is the first dimension and y is the second, with y varying the fastest in
   * terms of storage.
   *
   * @return IDs of the non-land objects in the local domain.
   */
  const int* get_nonzero_object_ids() const;

private:
  // Construct a ditributed grid from a NetCDF file describing the global domain
  Grid(MPI_Comm comm, int argc, char** argv, const std::string& filename);

private:
  MPI_Comm _comm;                   // MPI communicator
  int _rank = -1;                   // Process rank
  int _num_procs = -1;              // Total number of processes in communicator
  int _num_procs_x = -1;            // Total number of processes in x dimension
  int _num_procs_y = -1;            // Total number of processes in y dimension
  size_t _global_dim_x = 0;         // global longitude dimension
  size_t _global_dim_y = 0;         // global latitude dimension
  size_t _local_dim_x = 0;          // local longitude dimension
  size_t _local_dim_y = 0;          // local latitude dimension
  int _global_top_x = 0;            // global top left longitude
  int _global_top_y = 0;            // global top left latitude
  int _num_objects = 0;             // number of grid objects ignoring land mask
  int _num_nonzero_objects = 0;     // number of non-land grid objects
  std::vector<int> _land_mask = {}; // land mask values
  std::vector<int> _sparse_to_dense = {}; // map from sparse to dense index
  std::vector<int> _object_id = {};       // unique non-land object IDs
};
