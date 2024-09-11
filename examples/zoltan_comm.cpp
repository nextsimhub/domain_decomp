/*!
 * @file zoltan_comm.cpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @date 11 Sep 2024
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
    desc.add_options()("help,h", "Display this help message")(
        "grid,g", po::value<string>()->required(), "NetCDF grid file")("dim0",
        po::value<string>()->default_value("x"),
        "First spatial dimension in netCDF grid file")("dim1",
        po::value<string>()->default_value("y"), "Second spatial dimension in netCDF grid file")(
        "blk0", po::value<int>()->default_value(1), "Blocking factor in first dimension")(
        "blk1", po::value<int>()->default_value(1), "Blocking factor in second dimension")("mask,m",
        po::value<string>()->default_value("mask"), "Mask variable name in netCDF grid file")(
        "ignore-mask", po::bool_switch()->default_value(false), "Ignore mask in netCDF grid file");

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
    Grid* grid = Grid::create(comm, vm["grid"].as<string>(), vm["dim0"].as<string>(),
        vm["dim1"].as<string>(), vm["mask"].as<string>(), vm["blk0"].as<int>(),
        vm["blk1"].as<int>(), vm["ignore-mask"].as<bool>());

    // Create a Zoltan partitioner
    Partitioner* partitioner
        = Partitioner::Factory::create(comm, argc, argv, PartitionerType::Zoltan_RCB);

    // Partition grid
    partitioner->partition(*grid);

    // Retrieve neighbours
    vector<vector<int>> id_vec = { {}, {}, {}, {} };
    vector<int> top_halos, bottom_halos, left_halos, right_halos;
    partitioner->get_neighbours(id_vec[0], left_halos, 0);
    partitioner->get_neighbours(id_vec[1], right_halos, 1);
    partitioner->get_neighbours(id_vec[2], bottom_halos, 2);
    partitioner->get_neighbours(id_vec[3], top_halos, 3);

    // MPI ranks of neighbours in order: top, bottom, left, right
    vector<int> ids(id_vec[3]);
    ids.insert(ids.end(), id_vec[2].begin(), id_vec[2].end());
    ids.insert(ids.end(), id_vec[0].begin(), id_vec[0].end());
    ids.insert(ids.end(), id_vec[1].begin(), id_vec[1].end());

    // Create distributed neighbourhood communicator
    MPI_Comm comm_dist_graph;
    int err = MPI_Dist_graph_create_adjacent(comm, ids.size(), ids.data(), MPI_UNWEIGHTED,
        ids.size(), ids.data(), MPI_UNWEIGHTED, MPI_INFO_NULL, 0, &comm_dist_graph);
    if (err != MPI_SUCCESS)
        std::cerr << "MPI error" << std::endl;

    // Get my rank in the new communicator
    int mpi_rank;
    MPI_Comm_rank(comm_dist_graph, &mpi_rank);

    // Get information on my partition
    int global_0, global_1, local_ext_0, local_ext_1;
    partitioner->get_bounding_box(global_0, global_1, local_ext_0, local_ext_1);

    // Define data array
    vector<double> data((local_ext_0 + 2) * (local_ext_1 + 2), mpi_rank);

    // Define counts and displacements for MPI halo exchange
    vector<int> scounts(ids.size());
    vector<int> rcounts(ids.size());
    vector<MPI_Aint> sdispls(ids.size());
    vector<MPI_Aint> rdispls(ids.size());

    // Define MPI derived datatypes
    const int ndims { 2 };
    int sizes[ndims] = { local_ext_0 + 2, local_ext_1 + 2 };

    // Send subarray types
    vector<MPI_Datatype> subar_top(id_vec[3].size());
    vector<MPI_Datatype> subar_bottom(id_vec[2].size());
    vector<MPI_Datatype> subar_left(id_vec[0].size());
    vector<MPI_Datatype> subar_right(id_vec[1].size());

    // Receive subarray types
    vector<MPI_Datatype> ghost_top(id_vec[3].size());
    vector<MPI_Datatype> ghost_bottom(id_vec[2].size());
    vector<MPI_Datatype> ghost_left(id_vec[0].size());
    vector<MPI_Datatype> ghost_right(id_vec[1].size());

    // Mixed subarray types
    vector<MPI_Datatype> sendtypes(ids.size());
    vector<MPI_Datatype> recvtypes(ids.size());

    int cnt = 0;
    int offset = 0;

    // Top neighbours
    for (int i = 0; i < static_cast<int>(id_vec[3].size()); i++) {
        int subsizes[ndims] = { 1, top_halos[i] };

        int send_top_start[ndims] = { 1, offset + 1 };
        MPI_Type_create_subarray(
            ndims, sizes, subsizes, send_top_start, MPI_ORDER_C, MPI_DOUBLE, &subar_top[i]);
        MPI_Type_commit(&subar_top[i]);
        sendtypes[cnt] = subar_top[i];
        scounts[cnt] = 1;
        sdispls[cnt] = 0;

        int recv_top_start[ndims] = { 0, offset + 1 };
        MPI_Type_create_subarray(
            ndims, sizes, subsizes, recv_top_start, MPI_ORDER_C, MPI_DOUBLE, &ghost_top[i]);
        MPI_Type_commit(&ghost_top[i]);
        recvtypes[cnt] = ghost_top[i];
        rcounts[cnt] = 1;
        rdispls[cnt] = 0;

        offset += top_halos[i];
        cnt++;
    }

    // Bottom neighbours
    offset = 0;
    for (int i = 0; i < static_cast<int>(id_vec[2].size()); i++) {
        int subsizes[ndims] = { 1, bottom_halos[i] };

        int send_bottom_start[ndims] = { local_ext_0, offset + 1 };
        MPI_Type_create_subarray(
            ndims, sizes, subsizes, send_bottom_start, MPI_ORDER_C, MPI_DOUBLE, &subar_bottom[i]);
        MPI_Type_commit(&subar_bottom[i]);
        sendtypes[cnt] = subar_bottom[i];
        scounts[cnt] = 1;
        sdispls[cnt] = 0;

        int recv_bottom_start[ndims] = { local_ext_0 + 1, offset + 1 };
        MPI_Type_create_subarray(
            ndims, sizes, subsizes, recv_bottom_start, MPI_ORDER_C, MPI_DOUBLE, &ghost_bottom[i]);
        MPI_Type_commit(&ghost_bottom[i]);
        recvtypes[cnt] = ghost_bottom[i];
        rcounts[cnt] = 1;
        rdispls[cnt] = 0;

        offset += bottom_halos[i];
        cnt++;
    }

    // Left neighbours
    offset = 0;
    for (int i = 0; i < static_cast<int>(id_vec[0].size()); i++) {
        int subsizes[ndims] = { left_halos[i], 1 };

        int send_left_start[ndims] = { offset + 1, 1 };
        MPI_Type_create_subarray(
            ndims, sizes, subsizes, send_left_start, MPI_ORDER_C, MPI_DOUBLE, &subar_left[i]);
        MPI_Type_commit(&subar_left[i]);
        sendtypes[cnt] = subar_left[i];
        scounts[cnt] = 1;
        sdispls[cnt] = 0;

        int recv_left_start[ndims] = { offset + 1, 0 };
        MPI_Type_create_subarray(
            ndims, sizes, subsizes, recv_left_start, MPI_ORDER_C, MPI_DOUBLE, &ghost_left[i]);
        MPI_Type_commit(&ghost_left[i]);
        recvtypes[cnt] = ghost_left[i];
        rcounts[cnt] = 1;
        rdispls[cnt] = 0;

        offset += left_halos[i];
        cnt++;
    }

    // Right neighbours
    offset = 0;
    for (int i = 0; i < static_cast<int>(id_vec[1].size()); i++) {
        int subsizes[ndims] = { right_halos[i], 1 };

        int send_right_start[ndims] = { offset + 1, local_ext_1 };
        MPI_Type_create_subarray(
            ndims, sizes, subsizes, send_right_start, MPI_ORDER_C, MPI_DOUBLE, &subar_right[i]);
        MPI_Type_commit(&subar_right[i]);
        sendtypes[cnt] = subar_right[i];
        scounts[cnt] = 1;
        sdispls[cnt] = 0;

        int recv_right_start[ndims] = { offset + 1, local_ext_1 + 1 };
        MPI_Type_create_subarray(
            ndims, sizes, subsizes, recv_right_start, MPI_ORDER_C, MPI_DOUBLE, &ghost_right[i]);
        MPI_Type_commit(&ghost_right[i]);
        recvtypes[cnt] = ghost_right[i];
        rcounts[cnt] = 1;
        rdispls[cnt] = 0;

        offset += right_halos[i];
        cnt++;
    }

    // Perform halo exchange
    err = MPI_Neighbor_alltoallw(data.data(), scounts.data(), sdispls.data(), sendtypes.data(),
        data.data(), rcounts.data(), rdispls.data(), recvtypes.data(), comm_dist_graph);
    if (err != MPI_SUCCESS)
        std::cerr << "MPI error" << std::endl;

    // Cleanup
    delete grid;
    delete partitioner;

    for (int i = 0; i < static_cast<int>(sendtypes.size()); i++)
        MPI_Type_free(&sendtypes[i]);
    for (int i = 0; i < static_cast<int>(recvtypes.size()); i++)
        MPI_Type_free(&recvtypes[i]);

    MPI_Finalize();

    return 0;
}
