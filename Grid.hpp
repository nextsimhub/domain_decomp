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
 * A class that encapsulates a distributed 2D grid. The dimensions of the grid
 * as well as the ordering of dimensions is extracted from the netCDF grid file,
 * with the assumption that all variables defined in the netCDF file follow the
 * same convention in terms of dimension ordering. The grid is partitioned
 * evenly among processes using a 2D decomposition, ignoring any land mask. The
 * grid can be subsequently re-partitioned differently using a Partitioner.
 */
class Grid
{
  // Default grid metadata naming conventions
  const std::string data_id = "data";

public:
  // Disallow compiler-generated special functions
  Grid(const Grid&) = delete;
  Grid& operator=(const Grid&) = delete;

  /*!
   * @brief Destructor.
   */
  ~Grid(){};

  /*!
   * @brief Constructs a distributed 2D grid from a netCDF file describing the
   * global domain.
   *
   * Constructs a distributed grid of dimensions (dim0, dim1) from a netCDF
   * file. The dimension sizes as well as the ordering of dimensions is
   * extracted from the netCDF grid file and mapped appropriately to the (dim0,
   * dim1) dimensions of the Grid class, where dim1 is defined as the fastest
   * increasing dimension. For example, if the dimensions of interest in the
   * netCDF file are named x and y and the variables are dimensioned as (y, x),
   * then the y dimension will be mapped to dim0 and x to dim1 of the Grid
   * class, since the netCDF C/C++ convention is that the last dimension in the
   * CDL notation is the fastest increasing dimension. The default dimension
   * names for the netCDF file are "x" for dim0 and "y" for dim1, and the
   * default name for the land mask variable is "mask". the We assume that all
   * variables defined in the netCDF file follow the same convention in terms of
   * dimension ordering.
   *
   * @param comm MPI communicator.
   * @param filename Grid file in netCDF format.
   * @return A distributed 2D grid object partitioned evenly in terms of grid
   * points.
   */
  // We are using the named constructor idiom so that objects can only be
  // created in the heap to ensure it's dtor is executed before MPI_Finalize()
  static Grid* create(MPI_Comm comm, const std::string& filename,
                      const std::string dim0_id = "x",
                      const std::string dim1_id = "y",
                      const std::string mask_id = "mask");

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
   * @brief Returns the global extent in the first dimension of the grid.
   *
   * @return Global extent in first dimension.
   */
  int get_global_ext_0() const;

  /*!
   * @brief Returns the global extent in the second dimension of the grid.
   *
   * @return Global extent in second dimension.
   */
  int get_global_ext_1() const;

  /*!
   * @brief Returns the local extent in the first dimension of the grid.
   *
   * @return Local extent in first dimension.
   */
  int get_local_ext_0() const;

  /*!
   * @brief Returns the local extent in the second dimension of the grid.
   *
   * @return Local extent in second dimension.
   */
  int get_local_ext_1() const;

  /*!
   * @brief Sets the local extent in the first dimension of the grid.
   *
   * @param val Value to set.
   */
  void set_local_ext_0(int val);

  /*!
   * @brief Sets the local extent in the second dimension of the grid.
   *
   * @param val Value to set.
   */
  void set_local_ext_1(int val);

  /*!
   * @brief Returns the global coordinate in the first dimension of the upper
   * left corner in this process's partition.
   *
   * @return Global coordinate in the first dimension of the upper left corner
   * in this process's partition.
   */
  int get_global_0() const;

  /*!
   * @brief Returns the global coordinate in the second dimension of the upper
   * left corner in this process's partition.
   *
   * @return Global coordinate in the second dimension of the upper left corner
   * in this process's partition.
   */
  int get_global_1() const;

  /*!
   * @brief Sets the global coordinate in the first dimension of the upper
   * left corner in this process's partition.
   *
   * @param val Coordinate to set.
   */
  void set_global_0(int val);

  /*!
   * @brief Sets the global coordinate in the second dimension of the upper
   * left corner in this process's partition.
   *
   * @param val Coordinate to set.
   */
  void set_global_1(int val);

  /*!
   * @brief Returns the number of processes in the first dimension of the grid.
   *
   * @return Number of processes in first dimension.
   */
  int get_num_procs_0() const;

  /*!
   * @brief Returns the number of processes in the second dimension of the grid.
   *
   * @return Number of processes in second dimension.
   */
  int get_num_procs_1() const;

  /*!
   * @brief Returns the global land mask dimensioned (dim0, dim1), where dim0 is
   * the first dimension and dim1 is the second, with dim1 varying the fastest
   * in terms of storage.
   *
   * @return Global land mask.
   */
  const int* get_land_mask() const;

  /*!
   * @brief Returns the index mapping of sparse to dense representation, where
   * dim0 is the first dimension and dim1 is the second, with dim1 varying the
   * fastest in terms of storage.
   *
   * @return Index mapping of sparse to dense representation.
   */
  const int* get_sparse_to_dense() const;

  /*!
   * @brief Returns the IDs of the non-land objects in the local domain, where
   * dim0 is the first dimension and dim1 is the second, with dim1 varying the
   * fastest in terms of storage.
   *
   * @return IDs of the non-land objects in the local domain.
   */
  const int* get_nonzero_object_ids() const;

private:
  // Construct a ditributed grid from a NetCDF file describing the global domain
  Grid(MPI_Comm comm, const std::string& filename, const std::string& dim0_id,
       const std::string& dim1_id, const std::string& mask_id);

private:
  MPI_Comm _comm;           // MPI communicator
  int _rank = -1;           // Process rank
  int _num_procs = -1;      // Total number of processes in communicator
  int _num_procs_0 = -1;    // Total number of processes in first dimension
  int _num_procs_1 = -1;    // Total number of processes in second dimension
  size_t _global_ext_0 = 0; // Global extent in first dimension
  size_t _global_ext_1 = 0; // Global extent in second dimension
  size_t _local_ext_0 = 0;  // Local extent in first dimension
  size_t _local_ext_1 = 0;  // Local extent in second dimension
  int _global_0 = -1; // Upper left corner global coordinate in first dimension
  int _global_1 = -1; // Upper left corner global coordinate in second dimension
  int _num_objects = 0;             // Number of grid objects ignoring land mask
  int _num_nonzero_objects = 0;     // Number of non-land grid objects
  std::vector<int> _land_mask = {}; // Land mask values
  std::vector<int> _sparse_to_dense = {}; // Map from sparse to dense index
  std::vector<int> _object_id = {};       // Unique non-land object IDs
};
