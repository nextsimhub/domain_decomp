/*!
 * @file main.cpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @date 05 Nov 2024
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

bool validateOrder(const std::string DimOrderStr)
{
    if (DimOrderStr != "xy" && DimOrderStr != "yx") {
        cerr << "ERROR: invalid option. [order] must be either 'xy' or 'yx'." << endl;
        return false;
    }
    return true;
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
        ("ignore-mask,i", po::bool_switch()->default_value(false), "Ignore mask in netCDF grid file")
        ("periodic-x,px", po::bool_switch()->default_value(false), "Periodicity in x-direction")
        ("periodic-y,py", po::bool_switch()->default_value(false), "Periodicity in y-direction")
        ("output-prefix,op", po::value<string>()->default_value(""), "Prefix for output filenames")
        ("tripolar,t", po::bool_switch()->default_value(false),
            "Use split-communicator tri-polar decomposition");
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

    if (!validateOrder(vm["order"].as<string>())) {
        return 1;
    }
    std::vector<int> order = dimOrderFromStr(vm["order"].as<string>());
    string prefix = vm["output-prefix"].as<string>();
    if (prefix.find('/') != std::string::npos) {
        throw std::invalid_argument("prefix must not contain '/'");
    }
    if (!prefix.empty()) {
        prefix += "_";
    }

    // Capture tripolar flag
    bool tripolar = vm["tripolar"].as<bool>();

    // Build grid from netCDF file
    Grid* grid = Grid::create(comm, vm["grid"].as<string>(), vm["xdim"].as<string>(),
        vm["ydim"].as<string>(), order, vm["mask"].as<string>(), vm["ignore-mask"].as<bool>(),
        vm["periodic-x"].as<bool>(), vm["periodic-y"].as<bool>());

    // only used for tripolar case
    Grid* subgrid = new Grid(*grid);

    // Split communicator into two sub-communicators by X-coordinate (East/West)
    int g0, g1, le0, le1;
    grid->get_bounding_box(g0, g1, le0, le1);

    MPI_Comm sub_comm = comm; // default: no split
    if (tripolar) {
        int xMid = grid->getGlobalExt()[0] / 2; // X midpoint for East/West split
        int color = (g0 < xMid) ? 0 : 1; // 0 = West, 1 = East
        int world_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, &sub_comm);
    }

    // Create a Zoltan partitioner (on sub_comm if tripolar, otherwise on comm)
    Partitioner* partitioner
        = Partitioner::Factory::create(comm, argc, argv, PartitionerType::Zoltan_RCB);

    Partitioner* subpartitioner;

    // Partition grid
    if (tripolar) {
        // Create a copy of the grid on each rank
        // overwrite global position and extents for subgrids
        subpartitioner
            = Partitioner::Factory::create(sub_comm, argc, argv, PartitionerType::Zoltan_RCB);
        subgrid->set_global({ 0, g1 });
        auto globalExt = grid->getGlobalExt();
        subgrid->set_globalExt({ globalExt[0] / 2, globalExt[1] });
        subgrid->set_comm(sub_comm);
        subgrid->recompute_ids();
        subpartitioner->partition(*subgrid);
    } else {
        partitioner->partition(*grid);
    }

    // Store partitioning results in netCDF file
    int numProcs;
    MPI_Comm_size(comm, &numProcs);

    if (tripolar) {
        auto globalExt = subpartitioner->getGlobalNew();
        partitioner->setGlobalNew({ globalExt[0] + g0, globalExt[1] });
        auto localExt = subpartitioner->getLocalExtNew();
        partitioner->setLocalExtNew(localExt);
        partitioner->setTotalNumProcs(numProcs);
        partitioner->setGlobalExt(grid->getGlobalExt());
    }

    // Find my neighbours
    partitioner->discover_neighbours();

    // Free the sub-communicator
    if (tripolar) {
        MPI_Comm_free(&sub_comm);
    }

    // TODO: gather _procId results across both halves onto MPI_COMM_WORLD and
    // re-run neighbour discovery. Until then, saveMask/saveMetadata use the
    // partitioner's _comm which is the (now-freed) sub_comm in tripolar mode.
    // For now we skip saving in tripolar mode to avoid using a freed communicator.
    if (!tripolar) {
        partitioner->saveMask(prefix + "partition_mask_" + to_string(numProcs) + ".nc");
        partitioner->saveMetadata(prefix + "partition_metadata_" + to_string(numProcs) + ".nc");
    } else {
        partitioner->saveMetadata(prefix + "partition_metadata_" + to_string(numProcs) + ".nc");
    }

    // Cleanup
    delete grid;
    delete subgrid;
    delete partitioner;

    // Finalize MPI
    MPI_Finalize();

    return 0;
}
