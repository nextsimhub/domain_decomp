/*!
 * @file Partitioner.hpp
 * @date 25 June 2022
 * @author Athena Elafrou <ae488@cam.ac.uk>
 */

#include "Partitioner.hpp"

#include <netcdf>
#include <netcdf_par.h>

static int get_num_objects(void* data, int* ierr)
{
  Grid* grid = (Grid*)data;
  *ierr = ZOLTAN_OK;

  return grid->get_num_nonzero_objects();
}

static void get_object_list(void* data, int num_gid_entries,
                            int num_lid_entries, ZOLTAN_ID_PTR global_ids,
                            ZOLTAN_ID_PTR local_ids, int wgt_dim,
                            float* obj_wgts, int* ierr)
{
  Grid* grid = (Grid*)data;
  *ierr = ZOLTAN_OK;

  // In this example, return the IDs of our objects, but no weights
  // Zoltan will assume equally weighted objects
  for (int i = 0; i < grid->get_num_nonzero_objects(); i++) {
    global_ids[i] = grid->get_nonzero_object_ids()[i];
    local_ids[i] = i;
  }
  return;
}

static int get_num_geometry(void* data, int* ierr)
{
  *ierr = ZOLTAN_OK;
  return 2;
}

static void get_geometry_list(void* data, int num_gid_entries,
                              int num_lid_entries, int num_obj,
                              ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
                              int num_dim, double* geom_vec, int* ierr)
{
  Grid* grid = (Grid*)data;

  if ((num_gid_entries != 1) || (num_lid_entries != 1) || (num_dim != 2)) {
    *ierr = ZOLTAN_FATAL;
    return;
  }

  *ierr = ZOLTAN_OK;

  // Use grid coordinates
  for (int i = 0; i < num_obj; i++) {
    geom_vec[2 * i] = grid->get_nonzero_object_ids()[i] / grid->_global_dim_x;
    geom_vec[2 * i + 1]
        = grid->get_nonzero_object_ids()[i] % grid->_global_dim_x;
  }

  return;
}

Partitioner::Partitioner(MPI_Comm comm, int argc, char** argv)
{
  _comm = comm;

  // Initialize Zoltan
  float version;
  int ret = Zoltan_Initialize(argc, argv, &version);
  if (ret != ZOLTAN_OK) {
    std::cout << "Zoltan initialization failed on process " << _rank
              << std::endl;
    exit(EXIT_FAILURE);
  }

  // Create a Zoltan object
  _zoltan = new Zoltan(comm);
  if (ret != ZOLTAN_OK) {
    std::cout << "Creating Zoltan object failed on process " << _rank
              << std::endl;
    exit(EXIT_FAILURE);
  }
}

Partitioner* Partitioner::create(MPI_Comm comm, int argc, char** argv)
{
  return new Partitioner(comm, argc, argv);
}

void Partitioner::partition(Grid& grid)
{
  // Set Zoltan parameters for RCB partitioning
  // General parameters
  // _zoltan->Set_Param("DEBUG_LEVEL", "7");
  _zoltan->Set_Param("NUM_GID_ENTRIES", "1");
  _zoltan->Set_Param("NUM_LID_ENTRIES", "1");
  _zoltan->Set_Param("CHECK_GEOM", "1");
  _zoltan->Set_Param("KEEP_CUTS", "1");
  // RCB parameters
  std::string method = "RCB";
  _zoltan->Set_Param("LB_METHOD", "RCB");
  _zoltan->Set_Param("RCB_OUTPUT_LEVEL", "1");
  _zoltan->Set_Param("RCB_RECTILINEAR_BLOCKS", "1");
  _zoltan->Set_Param("RCB_SET_DIRECTIONS", "3");
  // Query functions
  _zoltan->Set_Num_Obj_Fn(get_num_objects, &grid);
  _zoltan->Set_Obj_List_Fn(get_object_list, &grid);
  _zoltan->Set_Num_Geom_Fn(get_num_geometry, &grid);
  _zoltan->Set_Geom_Multi_Fn(get_geometry_list, &grid);

  // Partition
  int changes;
  int num_gid_entries;
  int num_lid_entries;
  int num_import;
  ZOLTAN_ID_PTR import_global_ids;
  ZOLTAN_ID_PTR import_local_ids;
  int* import_procs;
  int* import_to_part;
  int num_export;
  ZOLTAN_ID_PTR export_global_ids;
  ZOLTAN_ID_PTR export_local_ids;
  int* export_procs;
  int* export_to_part;

  int ret = _zoltan->LB_Partition(
      changes, num_gid_entries, num_lid_entries, num_import, import_global_ids,
      import_local_ids, import_procs, import_to_part, num_export,
      export_global_ids, export_local_ids, export_procs, export_to_part);

  if (ret != ZOLTAN_OK) {
    std::cout << "Partitioning failed on process " << _rank << std::endl;
    delete _zoltan;
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  std::vector<int> sparse_part_ids(grid.get_num_nonzero_objects(), _rank);
  for (int i = 0; i < num_export; i++) {
    sparse_part_ids[export_local_ids[i]] = export_to_part[i];
  }

  _dense_part_ids.resize(grid._num_objects, _rank);
  int cnt = 0;
  for (int i = 0; i < grid._num_objects; i++) {
    if (grid._land_mask[i] == 1) {
      _dense_part_ids[i] = sparse_part_ids[cnt++];
    }
  }
  sparse_part_ids.clear();

  // Free the arrays allocated by Zoltan
  Zoltan::LB_Free_Part(&import_global_ids, &import_local_ids, &import_procs,
                       &import_to_part);
  Zoltan::LB_Free_Part(&export_global_ids, &export_local_ids, &export_procs,
                       &export_to_part);

  _global_dim_y = grid._global_dim_y;
  _global_dim_x = grid._global_dim_x;
  _local_dim_y = grid._local_dim_y;
  _local_dim_x = grid._local_dim_x;
  _global_top_y = grid._global_top_y;
  _global_top_x = grid._global_top_x;
}

void Partitioner::save(const std::string& filename) const
{
  // Use C API for parallel I/O
  const int NDIMS = 2;
  int ncid, vid, vid_lat, vid_lon, dimids[NDIMS];
  size_t start[NDIMS], count[NDIMS];

  nc_create_par(filename.c_str(), NC_MPIIO | NC_NETCDF4, _comm, MPI_INFO_NULL,
                &ncid);
  nc_put_att_int(ncid, NC_GLOBAL, "num_processors", NC_INT, 1, &_num_procs);

  // Create 2 dimensions
  // The values to be written are associated with the netCDF variable by
  // assuming that the last dimension of the netCDF variable varies fastest in
  // the C interface
  nc_def_dim(ncid, "y", _global_dim_y, &dimids[0]);
  nc_def_dim(ncid, "x", _global_dim_x, &dimids[1]);
  // Create variables
  nc_def_var(ncid, "part_id", NC_INT, NDIMS, dimids, &vid);
  nc_def_var(ncid, "lat", NC_DOUBLE, NDIMS, dimids, &vid_lat);
  nc_def_var(ncid, "lon", NC_DOUBLE, NDIMS, dimids, &vid_lon);
  // int init_val = -1;
  // Create a fill value for land data
  // nc_def_var_fill(ncid, vid, NC_FILL, &init_val);
  // Write metadata to file
  nc_enddef(ncid);
  // Set up slab for this process
  start[0] = _global_top_y;
  start[1] = _global_top_x;
  count[0] = _local_dim_y;
  count[1] = _local_dim_x;

  // Store data
  nc_var_par_access(ncid, vid, NC_COLLECTIVE);
  nc_put_vara_int(ncid, vid, start, count, _dense_part_ids.data());
  nc_close(ncid);
}
