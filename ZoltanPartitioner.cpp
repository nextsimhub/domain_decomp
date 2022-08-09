/*!
 * @file ZoltanPartitioner.cpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @date 25 June 2022
 */

#include "ZoltanPartitioner.hpp"

#include <netcdf>
#include <netcdf_par.h>
#include <unordered_map>

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
    geom_vec[2 * i]
        = grid->get_nonzero_object_ids()[i] / grid->get_global_dim_y();
    geom_vec[2 * i + 1]
        = grid->get_nonzero_object_ids()[i] % grid->get_global_dim_y();
  }

  return;
}

ZoltanPartitioner::ZoltanPartitioner(MPI_Comm comm, int argc, char** argv)
    : Partitioner(comm, argc, argv)
{
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

ZoltanPartitioner* ZoltanPartitioner::create(MPI_Comm comm, int argc,
                                             char** argv)
{
  return new ZoltanPartitioner(comm, argc, argv);
}

void ZoltanPartitioner::partition(Grid& grid)
{
  _num_procs_x = grid.get_num_procs_x();
  _num_procs_y = grid.get_num_procs_y();
  _global_dim_x = grid.get_global_dim_x();
  _global_dim_y = grid.get_global_dim_y();
  _local_dim_x = grid.get_local_dim_x();
  _local_dim_y = grid.get_local_dim_y();
  _global_top_x = grid.get_global_top_x();
  _global_top_y = grid.get_global_top_y();

  // Set Zoltan parameters for RCB partitioning
  // General parameters
  //_zoltan->Set_Param("DEBUG_LEVEL", "7");
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

  if (grid.get_num_objects() != grid.get_num_nonzero_objects()) {
    const int* land_mask = grid.get_land_mask();
    _proc_id.resize(grid.get_num_objects(), -1);
    for (int i = 0; i < grid.get_num_objects(); i++) {
      if (land_mask[i] > 0) {
        _proc_id[i] = _rank;
      }
    }

    const int* sparse_to_dense = grid.get_sparse_to_dense();
    for (int i = 0; i < num_export; i++) {
      _proc_id[sparse_to_dense[export_local_ids[i]]] = export_procs[i];
    }
  } else {
    _proc_id.resize(grid.get_num_objects(), _rank);
    for (int i = 0; i < num_export; i++) {
      _proc_id[export_local_ids[i]] = export_procs[i];
    }
  }

  // Free the arrays allocated by Zoltan
  Zoltan::LB_Free_Part(&import_global_ids, &import_local_ids, &import_procs,
                       &import_to_part);
  Zoltan::LB_Free_Part(&export_global_ids, &export_local_ids, &export_procs,
                       &export_to_part);
}
