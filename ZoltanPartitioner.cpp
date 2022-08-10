/*!
 * @file ZoltanPartitioner.cpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @date 25 June 2022
 */

#include "ZoltanPartitioner.hpp"

#include <algorithm>
#include <unordered_map>

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
    geom_vec[2 * i]
        = grid->get_nonzero_object_ids()[i] / grid->get_global_ext_1();
    geom_vec[2 * i + 1]
        = grid->get_nonzero_object_ids()[i] % grid->get_global_ext_1();
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
  // Load initial grid state
  _num_procs_0 = grid.get_num_procs_0();
  _num_procs_1 = grid.get_num_procs_1();
  _global_ext_0 = grid.get_global_ext_0();
  _global_ext_1 = grid.get_global_ext_1();
  _local_ext_0 = grid.get_local_ext_0();
  _local_ext_1 = grid.get_local_ext_1();
  _global_0 = grid.get_global_0();
  _global_1 = grid.get_global_1();

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

  // Find new bounding boxes for each process
  std::vector<int> min_x(_num_procs, _global_ext_0);
  std::vector<int> min_y(_num_procs, _global_ext_1);
  std::vector<int> max_x(_num_procs, -1);
  std::vector<int> max_y(_num_procs, -1);
  for (int i = 0; i < _proc_id.size(); i++) {
    // Find global 2D coordinates of element
    int x = i / _local_ext_1 + (_rank / _num_procs_1) * _local_ext_0;
    int y = i % _local_ext_1 + (_rank % _num_procs_1) * _local_ext_1;
    if (_proc_id[i] != -1) {
      if (x > max_x[_proc_id[i]])
        max_x[_proc_id[i]] = x;
      if (y > max_y[_proc_id[i]])
        max_y[_proc_id[i]] = y;
      if (x < min_x[_proc_id[i]])
        min_x[_proc_id[i]] = x;
      if (y < min_y[_proc_id[i]])
        min_y[_proc_id[i]] = y;
    }
  }

  // Find global bounding boxes for each process and update grid
  std::vector<int> global_min_x(_num_procs);
  std::vector<int> global_min_y(_num_procs);
  std::vector<int> global_max_x(_num_procs);
  std::vector<int> global_max_y(_num_procs);
  MPI_Allreduce(min_x.data(), global_min_x.data(), _num_procs, MPI_INT, MPI_MIN,
                _comm);
  MPI_Allreduce(min_y.data(), global_min_y.data(), _num_procs, MPI_INT, MPI_MIN,
                _comm);
  MPI_Allreduce(max_x.data(), global_max_x.data(), _num_procs, MPI_INT, MPI_MAX,
                _comm);
  MPI_Allreduce(max_y.data(), global_max_y.data(), _num_procs, MPI_INT, MPI_MAX,
                _comm);

  // Expand bounding boxes to account for land points in x dimension
  std::vector<int> coords(_num_procs_0 * 2);
  for (int j = 0; j < _num_procs_1; j++) {
    const int step
        = (_num_procs_0 == _num_procs || _num_procs_1 == _num_procs) ? 1 : 2;
    // Sort points on x dimension for processes in column j
    int cnt = 0;
    for (int p = j; p < _num_procs; p += step) {
      coords[2 * cnt] = global_min_x[p];
      coords[2 * cnt + 1] = global_max_x[p];
      cnt++;
    }
    std::sort(coords.begin(), coords.end());
    for (int i = 1; i < coords.size() - 1; i += 2) {
      if (coords[i + 1] - coords[i] != 1)
        coords[i] += coords[i + 1] - coords[i] - 1;
    }
    // Correct first and last point if necessary
    coords[0] = 0;
    coords[2 * _num_procs_0 - 1] = _global_ext_0 - 1;

    // Update bounding boxes
    cnt = 0;
    for (int p = j; p < _num_procs; p += step) {
      global_min_x[p] = coords[2 * cnt];
      global_max_x[p] = coords[2 * cnt + 1];
      cnt++;
    }
  }

  // Expand bounding boxes to account for land points in y dimension
  coords.clear();
  coords.resize(_num_procs_1 * 2);
  for (int j = 0; j < _num_procs_0; j++) {
    const int step
        = (_num_procs_0 == _num_procs || _num_procs_1 == _num_procs) ? 1 : 2;
    // Sort points on y dimension for processes in row j
    int cnt = 0;
    for (int p = step * j; p < (j + 1) * _num_procs_1; p++) {
      coords[2 * cnt] = global_min_y[p];
      coords[2 * cnt + 1] = global_max_y[p];
      cnt++;
    }
    std::sort(coords.begin(), coords.end());
    for (int i = 1; i < coords.size() - 1; i += 2) {
      if (coords[i + 1] - coords[i] != 1)
        coords[i] += coords[i + 1] - coords[i] - 1;
    }
    // Correct first and last point if necessary
    coords[0] = 0;
    coords[2 * _num_procs_1 - 1] = _global_ext_1 - 1;

    cnt = 0;
    for (int p = step * j; p < (j + 1) * _num_procs_1; p++) {
      global_min_y[p] = coords[2 * cnt];
      global_max_y[p] = coords[2 * cnt + 1];
      cnt++;
    }
  }

  // Update internal state
  _global_0_cur = global_min_x[_rank];
  _global_1_cur = global_min_y[_rank];
  _local_ext_0_cur = global_max_x[_rank] - _global_0_cur + 1;
  _local_ext_1_cur = global_max_y[_rank] - _global_1_cur + 1;

  // Update grid with new boxes
  grid.set_global_0(_global_0_cur);
  grid.set_global_1(_global_1_cur);
  grid.set_local_ext_0(_local_ext_0_cur);
  grid.set_local_ext_1(_local_ext_1_cur);
}
