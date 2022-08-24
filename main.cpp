/*!
 * @file main.cpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @date 1 May 2022
 */

#include <cstdio>
#include <iostream>

#include "Grid.hpp"
#include "Partitioner.hpp"

#include <boost/program_options.hpp>
#include <mpi.h>

using namespace std;
namespace po = boost::program_options;

int main(int argc, char* argv[])
{
  // Initialize MPI
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);

  // Configure command line options for dimension names and mask variable name
  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "Display this help message")
    ("grid,g", po::value<string>()->required(),
     "NetCDF grid file")
    ("dim0", po::value<string>()->default_value("x"),
     "First spatial dimension in netCDF grid file")
    ("dim1", po::value<string>()->default_value("y"),
     "Second spatial dimension in netCDF grid file")
    ("blk0", po::value<int>()->default_value(1),
     "Blocking factor in first dimension")
    ("blk1", po::value<int>()->default_value(1),
     "Blocking factor in second dimension")
    ("mask,m", po::value<string>()->default_value("mask"),
     "Mask variable name in netCDF grid file")
    ("ignore-mask", po::bool_switch()->default_value(false),
     "Ignore mask in netCDF grid file");

  // Parse optional command line options
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  if (vm.count("help")) {
    cout << "Usage: " << argv[0] << " [options]\n" << desc;
    return 0;
  }
  try {
    po::notify(vm);
  } catch (po::error& e) {
    cerr << "ERROR: " << e.what() << endl;
    return 1;
  }

  // Build grid from netCDF file
  Grid* grid
      = Grid::create(comm, vm["grid"].as<string>(), vm["dim0"].as<string>(),
                     vm["dim1"].as<string>(), vm["mask"].as<string>(),
                     vm["blk0"].as<int>(), vm["blk1"].as<int>(),
		     vm["ignore-mask"].as<bool>());

  // Create a Zoltan partitioner
  Partitioner* partitioner = Partitioner::Factory::create(
      comm, argc, argv, PartitionerType::Zoltan_RCB);

  // Partition grid
  partitioner->partition(*grid);

  // Store partitioning results in netCDF file
  int num_procs;
  MPI_Comm_size(comm, &num_procs);
  partitioner->save_mask("partition_mask_" + to_string(num_procs) + ".nc");
  partitioner->save_metadata("partition_metadata_" + to_string(num_procs)
                             + ".nc");

  // Cleanup
  delete grid;
  delete partitioner;

  // Finalize MPI
  MPI_Finalize();

  return 0;
}
