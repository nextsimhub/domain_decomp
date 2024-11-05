/*!
 * @file ZoltanPartitioner.cpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @date 05 Nov 2024
 */

#include "ZoltanPartitioner.hpp"
#include "Utils.hpp"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <unordered_map>

#include <netcdf.h>
#include <netcdf_par.h>

static int get_num_objects(void* data, int* ierr)
{
    Grid* grid = (Grid*)data;
    *ierr = ZOLTAN_OK;

    return grid->get_num_nonzero_objects();
}

static void get_object_list(
    void* data, int, int, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int, float*, int* ierr)
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

static int get_num_geometry(void*, int* ierr)
{
    *ierr = ZOLTAN_OK;
    return 2;
}

static void get_geometry_list(void* data, int num_gid_entries, int num_lid_entries, int num_obj,
    ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int num_dim, double* geom_vec, int* ierr)
{
    Grid* grid = (Grid*)data;

    if ((num_gid_entries != 1) || (num_lid_entries != 1) || (num_dim != 2)) {
        *ierr = ZOLTAN_FATAL;
        return;
    }

    *ierr = ZOLTAN_OK;

    // Use grid coordinates
    for (int i = 0; i < num_obj; i++) {
        geom_vec[2 * i] = grid->get_nonzero_object_ids()[i] % grid->get_global_ext()[0];
        geom_vec[2 * i + 1] = grid->get_nonzero_object_ids()[i] / grid->get_global_ext()[0];
    }

    return;
}

ZoltanPartitioner::ZoltanPartitioner(MPI_Comm comm, int argc, char** argv)
    : Partitioner(comm)
{
    // Initialize Zoltan
    float version;
    int ret = Zoltan_Initialize(argc, argv, &version);
    if (ret != ZOLTAN_OK) {
        std::cerr << "Zoltan initialization failed on process " << _rank << std::endl;
        exit(EXIT_FAILURE);
    }

    // Create a Zoltan object
    _zoltan = std::make_unique<Zoltan>(comm);
    if (ret != ZOLTAN_OK) {
        std::cerr << "Creating Zoltan object failed on process " << _rank << std::endl;
        exit(EXIT_FAILURE);
    }
}

ZoltanPartitioner* ZoltanPartitioner::create(MPI_Comm comm, int argc, char** argv)
{
    return new ZoltanPartitioner(comm, argc, argv);
}

void ZoltanPartitioner::partition(Grid& grid)
{
    // Load initial grid state
    _num_procs = grid.get_num_procs();
    _global_ext = grid.get_global_ext();
    grid.get_bounding_box(_global[0], _global[1], _local_ext[0], _local_ext[1]);
    _px = grid.get_px();
    _py = grid.get_py();

    if (_total_num_procs == 1) {
        for (int idx = 0; idx < 2; idx++) {
            _global_new[idx] = _global[idx];
            _local_ext_new[idx] = _local_ext[idx];
        }

        if (grid.get_num_objects() != grid.get_num_nonzero_objects()) {
            const int* land_mask = grid.get_land_mask();
            _proc_id.resize(grid.get_num_objects(), -1);
            for (int i = 0; i < grid.get_num_objects(); i++) {
                if (land_mask[i] > 0) {
                    _proc_id[i] = _rank;
                }
            }
        } else {
            _proc_id.resize(grid.get_num_objects(), _rank);
        }

        return;
    }

    // Set Zoltan parameters for RCB partitioning
    // General parameters
    _zoltan->Set_Param("DEBUG_LEVEL", "0");
    _zoltan->Set_Param("LB_METHOD", "RCB");
    _zoltan->Set_Param("LB_APPROACH", "REPARTITION");
    _zoltan->Set_Param("NUM_GID_ENTRIES", "1");
    _zoltan->Set_Param("NUM_LID_ENTRIES", "1");
    _zoltan->Set_Param("CHECK_GEOM", "1");
    _zoltan->Set_Param("KEEP_CUTS", "1");
    _zoltan->Set_Param("REMAP", "0");
    // RCB parameters
    _zoltan->Set_Param("RCB_OUTPUT_LEVEL", "1");
    _zoltan->Set_Param("RCB_RECTILINEAR_BLOCKS", "1");
    _zoltan->Set_Param("RCB_LOCK_DIRECTIONS", "1");
    _zoltan->Set_Param("RCB_SET_DIRECTIONS", "1");
    _zoltan->Set_Param("RCB_RECOMPUTE_BOX", "1");
    _zoltan->Set_Param("AVERAGE_CUTS", "1");
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

    int ret = _zoltan->LB_Partition(changes, num_gid_entries, num_lid_entries, num_import,
        import_global_ids, import_local_ids, import_procs, import_to_part, num_export,
        export_global_ids, export_local_ids, export_procs, export_to_part);

    if (ret != ZOLTAN_OK) {
        std::cerr << "Partitioning failed on process " << _rank << std::endl;
        CHECK_MPI(MPI_Finalize());
        exit(EXIT_FAILURE);
    }

    // Find new bounding boxes for each process
    if (changes == 1) {
        int ndim;
        double mins[3];
        double maxs[3];
        _zoltan->RCB_Box(_rank, ndim, mins[0], mins[1], mins[2], maxs[0], maxs[1], maxs[2]);
        for (int idx = 0; idx < 2; idx++) {
            _global_new[idx] = (mins[idx] == -DBL_MAX) ? 0 : std::ceil(mins[idx]);
            int global_lower = (maxs[idx] == DBL_MAX) ? _global_ext[idx] : std::ceil(maxs[idx]);
            _local_ext_new[idx] = global_lower - _global_new[idx];
        }
    } else {
        for (int idx = 0; idx < 2; idx++) {
            _global_new[idx] = _global[idx];
            _local_ext_new[idx] = _local_ext[idx];
        }
    }

    // Adapt to blocking
    std::vector<int> global_ext_orig = grid.get_global_ext();
    for (int idx = 0; idx < 2; idx++) {
        if (_global_new[idx] + _local_ext_new[idx] > global_ext_orig[idx]) {
            _local_ext_new[idx] = global_ext_orig[idx] - _global_new[idx];
        }
    }

    // Find my neighbours
    discover_neighbours();

    // Find the process IDs of each grid point I own
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
    Zoltan::LB_Free_Part(&import_global_ids, &import_local_ids, &import_procs, &import_to_part);
    Zoltan::LB_Free_Part(&export_global_ids, &export_local_ids, &export_procs, &export_to_part);
}
