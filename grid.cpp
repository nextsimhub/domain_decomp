#include "grid.h"

#include <algorithm>
#include <cmath>

#include <netcdf>
#include <netcdf_par.h>

static void find_factors(int n, int& factor_a, int& factor_b)
{
  for (int i = 2; i * i <= n; i += 2) {
    if (n % i == 0) {
      factor_a = i;
      factor_b = n / factor_a;
    }
  }
}

Grid* Grid::create(MPI_Comm comm, int argc, char** argv, const string& filename)
{
  return new Grid(comm, argc, argv, filename);
}

Grid::Grid(MPI_Comm comm, int argc, char** argv, const string& filename)
    : _comm(comm)
{
  netCDF::NcFile nc_file(filename, netCDF::NcFile::read);

  MPI_Comm_rank(comm, &_rank);

  // Retrieve the dimension of each axis
  _global_dim_x = nc_file.getDim(x_axis_id).getSize();
  _global_dim_y = nc_file.getDim(y_axis_id).getSize();

  // Initiallly we partition assuming there is no land mask
  // Figure out my subset of objects
  MPI_Comm_size(comm, &_num_procs);

  // Start from a 2D decomposition
  int num_procs_y = 1;
  int num_procs_x = _num_procs;
  find_factors(_num_procs, num_procs_y, num_procs_x);

  _local_dim_x = ceil((float)_global_dim_x / (num_procs_x));
  _local_dim_y = ceil((float)_global_dim_y / (num_procs_y));
  _global_top_y = (_rank / num_procs_x) * _local_dim_y;
  _global_top_x = (_rank % num_procs_x) * _local_dim_x;
  if ((_rank % num_procs_x) == num_procs_x - 1) {
    _local_dim_x = _global_dim_x - (_rank % num_procs_x) * _local_dim_x;
  }
  if ((_rank / num_procs_x) == num_procs_y - 1) {
    _local_dim_y = _global_dim_y - (_rank / num_procs_x) * _local_dim_y;
  }
  _num_objects = _local_dim_x * _local_dim_y;

  // Retrieve the coordinate data of interest for partitioining
  netCDF::NcVar mask_var = nc_file.getVar(mask_id); // (y, x)
  vector<size_t> start(2);
  vector<size_t> count(2);
  vector<ptrdiff_t> stride(2);
  // Coordinate of first element
  start[0] = _global_top_y; // y
  start[1] = _global_top_x; // x
  // Number of elements in every dimension
  count[0] = _local_dim_y; // y
  count[1] = _local_dim_x; // x
  // Stride in every dimension
  stride[0] = 1; // y
  stride[1] = 1; // x
  _land_mask.resize(_num_objects);
  mask_var.getVar(start, count, stride, _land_mask.data());

  netCDF::NcVar lon_var = nc_file.getVar(lon_id); // (y, x)
  start[0] = _global_top_y;                       // y
  start[1] = _global_top_x;                       // x
  count[0] = _local_dim_y;                        // y
  count[1] = _local_dim_x;                        // x
  _lon.resize(_num_objects);
  lon_var.getVar(start, count, stride, _lon.data());

  netCDF::NcVar lat_var = nc_file.getVar(lat_id); // (y, x)
  _lat.resize(_num_objects);
  lat_var.getVar(start, count, stride, _lat.data());

  // Apply land mask
  for (int i = 0; i < _lat.size(); i++) {
    if (_land_mask[i] == 1) {
      _lat_nonzero.push_back(_lat[i]);
      _lon_nonzero.push_back(_lon[i]);
      int _local_y = i / _local_dim_x;
      int _local_x = i % _local_dim_x;
      int _global_y = _local_y + _global_top_y;
      int _global_x = _local_x + _global_top_x;
      map.push_back(_global_y * _global_dim_x + _global_x);
    }
  }

  _num_nonzero_objects = _lat_nonzero.size();

  // Set unique object IDs
  _object_ids.resize(_num_nonzero_objects);
  for (int i = 0; i < _num_nonzero_objects; i++) {
    _object_ids[i] = map[i];
  }

  nc_file.close();

  // Initialize Zoltan
  float version;
  int ret = Zoltan_Initialize(argc, argv, &version);
  if (ret != ZOLTAN_OK) {
    cout << "Zoltan initialization failed on process " << _rank << endl;
    exit(EXIT_FAILURE);
  }

  // Create a Zoltan object
  _zoltan = new Zoltan(comm);
  if (ret != ZOLTAN_OK) {
    cout << "Creating Zoltan object failed on process " << _rank << endl;
    exit(EXIT_FAILURE);
  }
}

Grid::~Grid() { delete _zoltan; }

void Grid::partition()
{
  // Set Zoltan parameters for RCB partitioning
  // General parameters
  // _zoltan->Set_Param("DEBUG_LEVEL", "7");
  _zoltan->Set_Param("NUM_GID_ENTRIES", "1");
  _zoltan->Set_Param("NUM_LID_ENTRIES", "1");
  _zoltan->Set_Param("CHECK_GEOM", "1");
  _zoltan->Set_Param("KEEP_CUTS", "1");
  // RCB parameters
  string method = "RCB";
  _zoltan->Set_Param("LB_METHOD", "RCB");
  _zoltan->Set_Param("RCB_OUTPUT_LEVEL", "1");
  _zoltan->Set_Param("RCB_RECTILINEAR_BLOCKS", "1");
  _zoltan->Set_Param("RCB_SET_DIRECTIONS", "3");
  // Query functions
  _zoltan->Set_Num_Obj_Fn(Grid::get_num_objects, this);
  _zoltan->Set_Obj_List_Fn(Grid::get_object_list, this);
  _zoltan->Set_Num_Geom_Fn(Grid::get_num_geometry, this);
  _zoltan->Set_Geom_Multi_Fn(Grid::get_geometry_list, this);

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
    cout << "Partitioning failed on process " << _rank << endl;
    delete _zoltan;
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  _zoltan_part_ids.resize(get_num_nonzero_objects());
  for (int i = 0; i < get_num_nonzero_objects(); i++) {
    _zoltan_part_ids[i] = _rank;
  }
  for (int i = 0; i < num_export; i++) {
    _zoltan_part_ids[export_local_ids[i]] = export_to_part[i];
  }

  // Free the arrays allocated by Zoltan
  Zoltan::LB_Free_Part(&import_global_ids, &import_local_ids, &import_procs,
                       &import_to_part);
  Zoltan::LB_Free_Part(&export_global_ids, &export_local_ids, &export_procs,
                       &export_to_part);
}

void Grid::save(const string& filename) const
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

  vector<int> dense_parts(_num_objects, -1);
  int cnt = 0;
  for (int i = 0; i < _num_objects; i++) {
    if (_land_mask[i] == 1) {
      dense_parts[i] = _zoltan_part_ids[cnt++];
    }
  }

  // Store data
  nc_var_par_access(ncid, vid, NC_COLLECTIVE);
  nc_put_vara_int(ncid, vid, start, count, dense_parts.data());
  nc_put_vara_double(ncid, vid_lat, start, count, _lat.data());
  nc_put_vara_double(ncid, vid_lon, start, count, _lon.data());
  nc_close(ncid);
}

int Grid::get_num_nonzero_objects() const { return _num_nonzero_objects; }
const int* Grid::get_nonzero_object_ids() const { return _object_ids.data(); }

int Grid::get_num_objects(void* data, int* ierr)
{
  Grid* grid = (Grid*)data;
  *ierr = ZOLTAN_OK;

  return grid->get_num_nonzero_objects();
}

void Grid::get_object_list(void* data, int num_gid_entries, int num_lid_entries,
                           ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
                           int wgt_dim, float* obj_wgts, int* ierr)
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

int Grid::get_num_geometry(void* data, int* ierr)
{
  *ierr = ZOLTAN_OK;
  return 2;
}

void Grid::get_geometry_list(void* data, int num_gid_entries,
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

  // Use geographic coordinates
  // for (int i = 0; i < num_obj ; i++) {
  //   geom_vec[2*i] = grid->lat_nonzero[i];
  //   geom_vec[2*i + 1] = grid->lon_nonzero[i];
  // }

  // Use grid coordinates
  for (int i = 0; i < num_obj; i++) {
    geom_vec[2 * i] = grid->get_nonzero_object_ids()[i] / grid->_global_dim_x;
    geom_vec[2 * i + 1]
        = grid->get_nonzero_object_ids()[i] % grid->_global_dim_x;
  }

  return;
}
