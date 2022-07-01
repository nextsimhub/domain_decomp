/*!
 * @file Grid.hpp
 * @date 1 May 2022
 * @author Athena Elafrou <ae488@cam.ac.uk>
 */

#pragma once

#include <string>
#include <vector>

#include <mpi.h>

//! A class that encapsulates the grid of the model
class Grid
{
  // Define grid metadata
  const std::string x_axis_id = "x";
  const std::string y_axis_id = "y";
  const std::string data_id = "data";
  const std::string mask_id = "mask";

public:
  // Disallow compiler-generated special functions
  Grid(const Grid&) = delete;
  Grid& operator=(const Grid&) = delete;

  // Destructor
  ~Grid(){};

  // Construct a distributed grid from a netCDF file describing the global
  // domain. We are using the named constructor idiom so that objects can only
  // be created in the heap to ensure it's dtor is executed before
  // MPI_Finalize()
  static Grid* create(MPI_Comm comm, int argc, char** argv,
                      const std::string& filename);

  // Returns the total number of objects in the local domain
  int get_num_objects() const;

  // Returns the number of non-land objects in the local domain
  int get_num_nonzero_objects() const;

  // Returns the global x dimension size
  int get_global_dim_x() const;

  // Returns the global y dimension size
  int get_global_dim_y() const;

  // Returns the local x dimension size
  int get_local_dim_x() const;

  // Returns the local y dimension size
  int get_local_dim_y() const;

  // Returns the global x coordinate of the top left object in this process's
  // partition
  int get_global_top_x() const;

  // Returns the global y coordinate of the top left object in this process's
  // partition
  int get_global_top_y() const;

  // Returns the land mask
  const int* get_land_mask() const;

  // Returns the index mapping of sparse to dense representation
  const int* get_sparse_to_dense() const;

  // Returns the IDs of the non-land objects in the local domain
  const int* get_nonzero_object_ids() const;

private:
  // Construct a ditributed grid from a netCDF file describing the global domain
  Grid(MPI_Comm comm, int argc, char** argv, const std::string& filename);

private:
  MPI_Comm _comm;                   // MPI communicator
  int _rank = -1;                   // Process rank
  int _num_procs = -1;              // Total number of processes in communicator
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
