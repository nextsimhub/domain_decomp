/*!
 * @file Partitioner.hpp
 * @date 25 June 2022
 * @author Athena Elafrou <ae488@cam.ac.uk>
 */

#pragma once

#include "Grid.hpp"

enum class PartitionerType { Zoltan };

/**
 * @brief Abstract polymorphic class that encapsulates a grid partitioner
 */
class Partitioner
{
public:
  // Disallow compiler-generated special functions
  Partitioner(const Partitioner&) = delete;
  Partitioner& operator=(const Partitioner&) = delete;

  // Destructor
  virtual ~Partitioner(){};

  // Partition a grid
  virtual void partition(Grid& grid) = 0;

  // Save the results of the domain decomposition in a netCDF file
  void save_mask(const std::string& filename) const;

  // Save the results of the domain decomposition in a netCDF file
  void save_metadata(const std::string& filename) const;

protected:
  // Construct a partitioner
  // We are using the named constructor idiom so that objects can only be
  // created in the heap to ensure it's dtor is executed before MPI_Finalize()
  Partitioner(MPI_Comm comm, int argc, char** argv);

protected:
  MPI_Comm _comm;                 // MPI communicator
  int _rank = -1;                 // Process rank
  int _num_procs = -1;            // Total number of processes in communicator
  int _num_procs_x = -1;          // Total number of processes in x axis
  int _num_procs_y = -1;          // Total number of processes in y axis
  int _global_dim_x = 0;          // Global longitude dimension
  int _global_dim_y = 0;          // Global latitude dimension
  int _local_dim_x = 0;           // Local longitude dimension (original)
  int _local_dim_y = 0;           // Local latitude dimension (original)
  int _global_top_x = -1;         // Global top left longitude (original)
  int _global_top_y = -1;         // Global top left latitude (original)
  std::vector<int> _proc_id = {}; // Process ids of partition (dense form)

public:
  struct Factory {
    static Partitioner* create(MPI_Comm comm, int argc, char** argv,
                               PartitionerType type);
  };
};
