#pragma once

#include <vector>

#include <mpi.h>
#include <zoltan_cpp.h>

using namespace std;

// Define grid metadata
const string x_axis_id = "x";
const string y_axis_id = "y";
const string mask_id = "top_level";
const string lon_id = "nav_lon";
const string lat_id = "nav_lat";

class Grid
{
public:
  // Disallow compiler-generated functions
  Grid() = delete;
  Grid(const Grid&) = delete;
  Grid& operator=(const Grid&) = delete;

  ~Grid();

  // Construct a ditributed grid from a netCDF file describing the global domain
  Grid(MPI_Comm comm, int argc, char** argv, const string& filename);

  // Balance the domain decomposition using RCB
  void partition();

  // Save the results of the domain decomposition in a netCDF file
  void save(const string& filename);

  // Returns the number of non-land objects in the local domain
  int get_num_nonzero_objects() const;

private:
  // Returns the IDs of the non-land objects in the local domain
  const int* get_nonzero_object_ids() const;

  // Returns the number of objects that are currently assigned to the processor
  static int get_num_objects(void* data, int* ierr);

  // Returns an array of unique global and local IDs for all objects assigned to
  // the processor
  static void get_object_list(void* data, int num_gid_entries,
                              int num_lid_entries, ZOLTAN_ID_PTR global_ids,
                              ZOLTAN_ID_PTR local_ids, int wgt_dim,
                              float* obj_wgts, int* ierr);

  // Returns the number of values needed to express the geometry of an object
  static int get_num_geometry(void* data, int* ierr);

  // Returns a vector of geometry values for a list of given objects
  // For object i (specified by global_ids[i*num_gid_entries] and
  // local_ids[i*num_lid_entries], i=0,1,...,num_obj-1), coordinate values
  // should be stored in geom_vec[i*num_dim:(i+1)*num_dim-1]
  static void get_geometry_list(void* data, int num_gid_entries,
                                int num_lid_entries, int num_obj,
                                ZOLTAN_ID_PTR global_ids,
                                ZOLTAN_ID_PTR local_ids, int num_dim,
                                double* geom_vec, int* ierr);

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
  vector<double> _lon = {};         // longitude values
  vector<double> _lat = {};         // latitude values
  vector<double> _lon_nonzero = {}; // longitude values for non-land objects
  vector<double> _lat_nonzero = {}; // latitude values for non-land objects
  vector<int> _land_mask = {};      // land mask values
  vector<int> map = {}; // index mapping from sparse to dense representation
  vector<int> _object_ids = {};      // unique non-land object IDs
  Zoltan* _zoltan = nullptr;         // Zoltan object
  vector<int> _zoltan_part_ids = {}; // processor ids for each object in the
                                     // local domain after partitioning
};
