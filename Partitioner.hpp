/*!
 * @file Partitioner.hpp
 * @date 25 June 2022
 * @author Athena Elafrou <ae488@cam.ac.uk>
 */

#pragma once

#include "Grid.hpp"

//! A class that encapsulates a grid partitioner
class Partitioner
{
public:
  // Disallow compiler-generated special functions
  Partitioner(const Partitioner&) = delete;
  Partitioner& operator=(const Partitioner&) = delete;

  // Destructor
  ~Partitioner() { delete _zoltan; }

  // Construct a Zoltan partitioner
  // We are using the named constructor idiom so that objects can only be
  // created in the heap to ensure it's dtor is executed before MPI_Finalize()
  static Partitioner* create(MPI_Comm comm, int argc, char** argv);

  // Partition a grid
  void partition(Grid& grid);

  // Save the results of the domain decomposition in a netCDF file
  void save(const std::string& filename) const;

private:
  // Construct a Zoltan partitioner
  Partitioner(MPI_Comm comm, int argc, char** argv);

private:
  MPI_Comm _comm;                 // MPI communicator
  int _rank = -1;                 // Process rank
  int _num_procs = -1;            // Total number of processes in communicator
  size_t _global_dim_x = 0;       // Global longitude dimension
  size_t _global_dim_y = 0;       // Global latitude dimension
  size_t _local_dim_x = 0;        // Local longitude dimension
  size_t _local_dim_y = 0;        // Local latitude dimension
  int _global_top_x = 0;          // Global top left longitude
  int _global_top_y = 0;          // Global top left latitude
  Zoltan* _zoltan = nullptr;      // Zoltan object
  std::vector<int> _proc_id = {}; // Process ids of partition (dense form)
};
