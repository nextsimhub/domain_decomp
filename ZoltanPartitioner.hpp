/*!
 * @file ZoltanPartitioner.hpp
 * @date 25 June 2022
 * @author Athena Elafrou <ae488@cam.ac.uk>
 */

#pragma once

#include "Grid.hpp"
#include "Partitioner.hpp"

#include <zoltan_cpp.h>

//! A class that encapsulates a grid partitioner
class ZoltanPartitioner final : public Partitioner
{
public:
  // Disallow compiler-generated special functions
  ZoltanPartitioner(const ZoltanPartitioner&) = delete;
  ZoltanPartitioner& operator=(const ZoltanPartitioner&) = delete;

  // Destructor
  ~ZoltanPartitioner() { delete _zoltan; }

  // Construct a partitioner
  // We are using the named constructor idiom so that objects can only be
  // created in the heap to ensure it's dtor is executed before MPI_Finalize()
  static ZoltanPartitioner* create(MPI_Comm comm, int argc, char** argv);

  // Partition a grid
  void partition(Grid& grid) override;

protected:
  // Construct a Zoltan partitioner
  ZoltanPartitioner(MPI_Comm comm, int argc, char** argv);

private:
  Zoltan* _zoltan = nullptr; // Zoltan object
};
