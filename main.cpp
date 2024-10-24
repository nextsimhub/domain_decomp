/*!
 * @file main.cpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @date 23 Aug 2024
 */

#include <cstdio>
#include <iostream>

#include "Grid.hpp"
#include "Partitioner.hpp"

#include <boost/program_options.hpp>
#include <mpi.h>
#include <vector>

using namespace std;
namespace po = boost::program_options;

void validateOrder(const std::string DimOrderStr)
{
    if (DimOrderStr != "xy" && DimOrderStr != "yx") {
        cerr << "ERROR: invalid option. [order] must be either 'xy' or 'yx'." << endl;
    }
    return;
}

std::vector<int> dimOrderFromStr(const std::string& DimOrderStr)
{
    if (DimOrderStr[0] == 'x') {
        return std::vector<int>({ 0, 1 });
    } else {
        return std::vector<int>({ 1, 0 });
    }
}

int main(int argc, char* argv[])
{
    // Initialize MPI
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Init(&argc, &argv);

    // Configure command line options for dimension names and mask variable name
    po::options_description desc("Options");
    // clang-format off
    desc.add_options()
        ("help,h", "Display this help message")
        ("grid,g", po::value<string>()->required(), "NetCDF grid file")
        ("xdim,x", po::value<string>()->default_value("x"), "Name of x dimension in netCDF grid file")
        ("ydim,y", po::value<string>()->default_value("y"), "Name of y dimension in netCDF grid file")
        ("order,o", po::value<string>()->default_value("yx"), "Order of dimensions in netCDF grid file, e.g., 'yx' or 'xy'")
        ("mask,m", po::value<string>()->default_value("mask"), "Mask variable name in netCDF grid file")
        ("ignore-mask,i", po::bool_switch()->default_value(false), "Ignore mask in netCDF grid file");
    // clang-format on

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

    validateOrder(vm["order"].as<string>());
    std::vector<int> order = dimOrderFromStr(vm["order"].as<string>());

    // Build grid from netCDF file
    Grid* grid = Grid::create(comm, vm["grid"].as<string>(), vm["xdim"].as<string>(),
        vm["ydim"].as<string>(), order, vm["mask"].as<string>(), vm["ignore-mask"].as<bool>());

    // Create a Zoltan partitioner
    Partitioner* partitioner
        = Partitioner::Factory::create(comm, argc, argv, PartitionerType::Zoltan_RCB);

    // Partition grid
    partitioner->partition(*grid);

    // Store partitioning results in netCDF file
    int num_procs;
    MPI_Comm_size(comm, &num_procs);
    partitioner->save_mask("partition_mask_" + to_string(num_procs) + ".nc");
    partitioner->save_metadata("partition_metadata_" + to_string(num_procs) + ".nc");

    // Cleanup
    delete grid;
    delete partitioner;

    // Finalize MPI
    MPI_Finalize();

    return 0;
}
