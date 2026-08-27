#include "analysis_io.h"

#include <hdf5.h>
#include <dirent.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "config.h"
#include "log.h"
#include "setup_mpi.h"

#ifdef GAMERA_NONORTHOGONAL_BACKEND
#include "nonorthogonal_background.h"
#include "nonorthogonal_driver.h"
#include "nonorthogonal_grid.h"
#include "nonorthogonal_legacy_adapter.h"
#include "nonorthogonal_storage.h"
#endif

enum {
  GAMERA_ANALYSIS_SCHEMA_VERSION = 1,
  GAMERA_ANALYSIS_FIELD_COUNT = 8,
  GAMERA_ANALYSIS_COORDINATE_COUNT = 3,
  GAMERA_ANALYSIS_NAME_BYTES = 16
};

#if defined(GAMERA_NONORTHOGONAL_BACKEND) && defined(H5_HAVE_PARALLEL)

static int patch_consensus(int local_failed) {
  int failed = local_failed;
  return MPI_Allreduce(&local_failed, &failed, 1, MPI_INT, MPI_MAX,
                       comm_cart) == MPI_SUCCESS
             ? failed
             : 1;
}

int initialize_analysis_sequence(void) {
  int sequence = -1;
  int directory_errno = 0;
  if (patch_rank == 0) {
    DIR *directory = opendir(".");
    if (directory == NULL) {
      sequence = -2;
      directory_errno = errno;
    } else {
      char prefix[64];
      snprintf(prefix, sizeof(prefix), "analysis_p%d_", patch_id);
      struct dirent *entry = NULL;
      while ((entry = readdir(directory)) != NULL) {
        int candidate = -1;
        char trailing = '\0';
        if (strncmp(entry->d_name, prefix, strlen(prefix)) == 0 &&
            sscanf(entry->d_name + strlen(prefix), "%6d.h5%c", &candidate,
                   &trailing) == 1 &&
            candidate > sequence) {
          sequence = candidate;
        }
      }
      closedir(directory);
    }
  }
  const int broadcast_status =
      MPI_Bcast(&sequence, 1, MPI_INT, 0, comm_cart);
  if (broadcast_status != MPI_SUCCESS || sequence < -1) {
    log_error("Compact-analysis sequence recovery failed: patch=%d "
              "patch_rank=%d sequence=%d opendir_errno=%d mpi_status=%d",
              patch_id, patch_rank, sequence, directory_errno,
              broadcast_status);
    return -1;
  }
  analysis_seq_num = sequence;
  return 0;
}

static hid_t create_parallel_file(const char *filename) {
  hid_t access = H5Pcreate(H5P_FILE_ACCESS);
  if (patch_consensus(access < 0)) {
    if (access >= 0) {
      H5Pclose(access);
    }
    return -1;
  }
  const int access_failed =
      H5Pset_fapl_mpio(access, comm_cart, MPI_INFO_NULL) < 0;
  if (patch_consensus(access_failed)) {
    H5Pclose(access);
    return -1;
  }
  const hid_t file =
      H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, access);
  H5Pclose(access);
  return file;
}

static int write_root_dataset(hid_t file, const char *name, hid_t file_type,
                              hid_t memory_type, int dimensions,
                              const hsize_t *extent, const void *value) {
  hid_t space = H5Screate_simple(dimensions, extent, NULL);
  hid_t dataset =
      space < 0 ? -1 : H5Dcreate2(file, name, file_type, space, H5P_DEFAULT,
                                  H5P_DEFAULT, H5P_DEFAULT);
  int failed = space < 0 || dataset < 0;
  if (!patch_consensus(failed) && patch_rank == 0 &&
      H5Dwrite(dataset, memory_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, value) <
          0) {
    failed = 1;
  }
  if (dataset >= 0) {
    H5Dclose(dataset);
  }
  if (space >= 0) {
    H5Sclose(space);
  }
  return patch_consensus(failed) ? -1 : 0;
}

static int write_parallel_hyperslab(hid_t file, const char *name,
                                    hid_t file_type, hid_t memory_type,
                                    int dimensions,
                                    const hsize_t *global_extent,
                                    const hsize_t *local_extent,
                                    const hsize_t *offset,
                                    const void *values) {
  hid_t file_space = H5Screate_simple(dimensions, global_extent, NULL);
  hid_t memory_space = H5Screate_simple(dimensions, local_extent, NULL);
  hid_t dataset = file_space < 0
                      ? -1
                      : H5Dcreate2(file, name, file_type, file_space,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t transfer = H5Pcreate(H5P_DATASET_XFER);
  int failed = file_space < 0 || memory_space < 0 || dataset < 0 ||
               transfer < 0;
  failed = patch_consensus(failed);
  if (!failed &&
      (H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset, NULL,
                           local_extent, NULL) < 0 ||
       H5Pset_dxpl_mpio(transfer, H5FD_MPIO_COLLECTIVE) < 0)) {
    failed = 1;
  }
  failed = patch_consensus(failed);
  if (!failed &&
      H5Dwrite(dataset, memory_type, memory_space, file_space, transfer,
               values) < 0) {
    failed = 1;
  }
  if (transfer >= 0) {
    H5Pclose(transfer);
  }
  if (dataset >= 0) {
    H5Dclose(dataset);
  }
  if (memory_space >= 0) {
    H5Sclose(memory_space);
  }
  if (file_space >= 0) {
    H5Sclose(file_space);
  }
  return patch_consensus(failed) ? -1 : 0;
}

static int write_common_metadata(hid_t file, double snapshot_time,
                                 int sequence) {
  const int schema = GAMERA_ANALYSIS_SCHEMA_VERSION;
  const int coordinate_backend = GAMERA_COORDINATE_BACKEND_ID;
  const int problem_id = GAMERA_PROBLEM_ID;
  const int global_cells[3] = {config.ni_global, config.nj_global,
                               config.nk_global};
  const double logical_lower[3] = {x1min_global, x2min_global, x3min_global};
  const double logical_upper[3] = {x1max_global, x2max_global, x3max_global};
  const hsize_t one[1] = {1};
  const hsize_t three[1] = {3};
  const hsize_t normalization_count[1] = {
      sizeof(Norm_t) / sizeof(double)};
  const double time_seconds = snapshot_time * norm_config.Time_Norm;
  if (write_root_dataset(file, "/analysis_schema", H5T_STD_I32LE,
                         H5T_NATIVE_INT, 1, one, &schema) != 0 ||
      write_root_dataset(file, "/patch_id", H5T_STD_I32LE, H5T_NATIVE_INT,
                         1, one, &patch_id) != 0 ||
      write_root_dataset(file, "/sequence", H5T_STD_I32LE, H5T_NATIVE_INT,
                         1, one, &sequence) != 0 ||
      write_root_dataset(file, "/coordinate_backend", H5T_STD_I32LE,
                         H5T_NATIVE_INT, 1, one, &coordinate_backend) != 0 ||
      write_root_dataset(file, "/problem_id", H5T_STD_I32LE,
                         H5T_NATIVE_INT, 1, one, &problem_id) != 0 ||
      write_root_dataset(file, "/global_cells", H5T_STD_I32LE,
                         H5T_NATIVE_INT, 1, three, global_cells) != 0 ||
      write_root_dataset(file, "/logical_lower", H5T_IEEE_F64LE,
                         H5T_NATIVE_DOUBLE, 1, three, logical_lower) != 0 ||
      write_root_dataset(file, "/logical_upper", H5T_IEEE_F64LE,
                         H5T_NATIVE_DOUBLE, 1, three, logical_upper) != 0 ||
      write_root_dataset(file, "/time_code", H5T_IEEE_F64LE,
                         H5T_NATIVE_DOUBLE, 1, one, &snapshot_time) != 0 ||
      write_root_dataset(file, "/time_seconds", H5T_IEEE_F64LE,
                         H5T_NATIVE_DOUBLE, 1, one, &time_seconds) != 0 ||
      write_root_dataset(file, "/normalization", H5T_IEEE_F64LE,
                         H5T_NATIVE_DOUBLE, 1, normalization_count,
                         &norm_config) != 0) {
    return -1;
  }

  static const char field_names[GAMERA_ANALYSIS_FIELD_COUNT]
                               [GAMERA_ANALYSIS_NAME_BYTES] = {
      "rho", "p", "vx", "vy", "vz", "Bx", "By", "Bz"};
  hid_t string_type = H5Tcopy(H5T_C_S1);
  int string_failed =
      string_type < 0 ||
      (string_type >= 0 &&
       (H5Tset_size(string_type, GAMERA_ANALYSIS_NAME_BYTES) < 0 ||
        H5Tset_strpad(string_type, H5T_STR_NULLTERM) < 0));
  if (patch_consensus(string_failed)) {
    if (string_type >= 0) {
      H5Tclose(string_type);
    }
    return -1;
  }
  const hsize_t names[1] = {GAMERA_ANALYSIS_FIELD_COUNT};
  const int result = write_root_dataset(file, "/field_names", string_type,
                                        string_type, 1, names, field_names);
  H5Tclose(string_type);
  return result;
}

static int write_analysis_grid(void) {
  char filename[128];
  snprintf(filename, sizeof(filename), "analysis_grid_p%d.h5", patch_id);
  int exists = 0;
  if (patch_rank == 0) {
    exists = access(filename, F_OK) == 0;
  }
  if (MPI_Bcast(&exists, 1, MPI_INT, 0, comm_cart) != MPI_SUCCESS) {
    return -1;
  }
  if (exists) {
    return 0;
  }

  const size_t local_cells = (size_t)config.ni * (size_t)config.nj *
                             (size_t)config.nk;
  if (local_cells == 0U ||
      local_cells > SIZE_MAX /
                        (GAMERA_ANALYSIS_COORDINATE_COUNT * sizeof(double))) {
    return -1;
  }
  double *coordinates =
      (double *)malloc(GAMERA_ANALYSIS_COORDINATE_COUNT * local_cells *
                       sizeof(*coordinates));
  const int vertex_count_i =
      config.ni + (proc_coords[0] + 1 == config.proc_dims[0]);
  const int vertex_count_j =
      config.nj + (proc_coords[1] + 1 == config.proc_dims[1]);
  const int vertex_count_k =
      config.nk + (proc_coords[2] + 1 == config.proc_dims[2]);
  const size_t local_vertices = (size_t)vertex_count_i *
                                (size_t)vertex_count_j *
                                (size_t)vertex_count_k;
  double *vertices =
      local_vertices >
              SIZE_MAX /
                  (GAMERA_ANALYSIS_COORDINATE_COUNT * sizeof(double))
          ? NULL
          : (double *)malloc(GAMERA_ANALYSIS_COORDINATE_COUNT *
                             local_vertices * sizeof(*vertices));
  if (patch_consensus(coordinates == NULL || vertices == NULL)) {
    free(coordinates);
    free(vertices);
    return -1;
  }
  for (int i = 0; i < config.ni; ++i) {
    for (int j = 0; j < config.nj; ++j) {
      for (int k = 0; k < config.nk; ++k) {
        const size_t local =
            ((size_t)i * (size_t)config.nj + (size_t)j) *
                (size_t)config.nk +
            (size_t)k;
        const int ii = i + NG;
        const int jj = j + NG;
        const int kk = k + NG;
        coordinates[0U * local_cells + local] = x1ctr[ii][jj][kk];
        coordinates[1U * local_cells + local] = x2ctr[ii][jj][kk];
        coordinates[2U * local_cells + local] = x3ctr[ii][jj][kk];
      }
    }
  }
  for (int i = 0; i < vertex_count_i; ++i) {
    for (int j = 0; j < vertex_count_j; ++j) {
      for (int k = 0; k < vertex_count_k; ++k) {
        const size_t local =
            ((size_t)i * (size_t)vertex_count_j + (size_t)j) *
                (size_t)vertex_count_k +
            (size_t)k;
        const int ii = i + NG;
        const int jj = j + NG;
        const int kk = k + NG;
        vertices[0U * local_vertices + local] = x1[ii][jj][kk];
        vertices[1U * local_vertices + local] = x2[ii][jj][kk];
        vertices[2U * local_vertices + local] = x3[ii][jj][kk];
      }
    }
  }

  hid_t file = create_parallel_file(filename);
  int failed = patch_consensus(file < 0);
  const hsize_t global[4] = {GAMERA_ANALYSIS_COORDINATE_COUNT,
                             (hsize_t)config.ni_global,
                             (hsize_t)config.nj_global,
                             (hsize_t)config.nk_global};
  const hsize_t local[4] = {GAMERA_ANALYSIS_COORDINATE_COUNT,
                            (hsize_t)config.ni, (hsize_t)config.nj,
                            (hsize_t)config.nk};
  const hsize_t offset[4] = {
      0, (hsize_t)(proc_coords[0] * config.ni),
      (hsize_t)(proc_coords[1] * config.nj),
      (hsize_t)(proc_coords[2] * config.nk)};
  const hsize_t global_vertices[4] = {
      GAMERA_ANALYSIS_COORDINATE_COUNT, (hsize_t)config.ni_global + 1U,
      (hsize_t)config.nj_global + 1U, (hsize_t)config.nk_global + 1U};
  const hsize_t local_vertex_extent[4] = {
      GAMERA_ANALYSIS_COORDINATE_COUNT, (hsize_t)vertex_count_i,
      (hsize_t)vertex_count_j, (hsize_t)vertex_count_k};
  if (!failed &&
      (write_common_metadata(file, time_sim, -1) != 0 ||
       write_parallel_hyperslab(file, "/coordinates", H5T_IEEE_F64LE,
                                H5T_NATIVE_DOUBLE, 4, global, local, offset,
                                coordinates) != 0 ||
       write_parallel_hyperslab(file, "/vertices", H5T_IEEE_F64LE,
                                H5T_NATIVE_DOUBLE, 4, global_vertices,
                                local_vertex_extent, offset, vertices) != 0)) {
    failed = 1;
  }
  if (file >= 0) {
    H5Fclose(file);
  }
  free(coordinates);
  free(vertices);
  if (patch_consensus(failed)) {
    return -1;
  }
  if (patch_rank == 0) {
    log_info("Wrote static compact analysis grid %s", filename);
  }
  return 0;
}

int dump_analysis_hdf5(void) {
  const gamera_no_grid *grid = gamera_no_legacy_grid();
  const gamera_no_storage *storage = gamera_no_legacy_storage();
#ifdef GAMERA_NONORTHOGONAL_HAS_BACKGROUND_FIELD
  const gamera_no_background_field *background =
      gamera_no_driver_background_field();
#else
  const gamera_no_background_field *background = NULL;
#endif
  int invalid = grid == NULL || storage == NULL ||
                storage->cell_magnetic == NULL ||
                (background != NULL && background->cell_magnetic == NULL);
  if (patch_consensus(invalid) || write_analysis_grid() != 0) {
    return -1;
  }

  const size_t local_cells = (size_t)config.ni * (size_t)config.nj *
                             (size_t)config.nk;
  if (local_cells == 0U ||
      local_cells >
          SIZE_MAX / (GAMERA_ANALYSIS_FIELD_COUNT * sizeof(float))) {
    return -1;
  }
  float *state = (float *)malloc(GAMERA_ANALYSIS_FIELD_COUNT * local_cells *
                                 sizeof(*state));
  if (patch_consensus(state == NULL)) {
    free(state);
    return -1;
  }
  int nonfinite = 0;
  for (int i = 0; i < config.ni; ++i) {
    for (int j = 0; j < config.nj; ++j) {
      for (int k = 0; k < config.nk; ++k) {
        const size_t local =
            ((size_t)i * (size_t)config.nj + (size_t)j) *
                (size_t)config.nk +
            (size_t)k;
        const size_t cell = gamera_no_index3(
            grid->cell_extent, (size_t)(i + NG), (size_t)(j + NG),
            (size_t)(k + NG));
        const double fluid[5] = {
            gas[0][gas_rho][i + NG][j + NG][k + NG],
            gas[0][gas_p][i + NG][j + NG][k + NG],
            gas[0][gas_v1][i + NG][j + NG][k + NG],
            gas[0][gas_v2][i + NG][j + NG][k + NG],
            gas[0][gas_v3][i + NG][j + NG][k + NG]};
        for (int field = 0; field < 5; ++field) {
          state[(size_t)field * local_cells + local] = (float)fluid[field];
          nonfinite |= !isfinite(fluid[field]);
        }
        for (int component = 0; component < 3; ++component) {
          double magnetic = storage->cell_magnetic[cell].value[component];
          if (background != NULL) {
            magnetic += background->cell_magnetic[cell].value[component];
          }
          state[(size_t)(5 + component) * local_cells + local] =
              (float)magnetic;
          nonfinite |= !isfinite(magnetic);
        }
      }
    }
  }
  if (patch_consensus(nonfinite)) {
    free(state);
    log_error("Refusing to write non-finite compact analysis state");
    return -1;
  }

  ++analysis_seq_num;
  char filename[128];
  snprintf(filename, sizeof(filename), "analysis_p%d_%06d.h5", patch_id,
           analysis_seq_num);
  hid_t file = create_parallel_file(filename);
  int failed = patch_consensus(file < 0);
  const hsize_t global[4] = {GAMERA_ANALYSIS_FIELD_COUNT,
                             (hsize_t)config.ni_global,
                             (hsize_t)config.nj_global,
                             (hsize_t)config.nk_global};
  const hsize_t local[4] = {GAMERA_ANALYSIS_FIELD_COUNT,
                            (hsize_t)config.ni, (hsize_t)config.nj,
                            (hsize_t)config.nk};
  const hsize_t offset[4] = {
      0, (hsize_t)(proc_coords[0] * config.ni),
      (hsize_t)(proc_coords[1] * config.nj),
      (hsize_t)(proc_coords[2] * config.nk)};
  if (!failed &&
      (write_common_metadata(file, time_sim, analysis_seq_num) != 0 ||
       write_parallel_hyperslab(file, "/state", H5T_IEEE_F32LE,
                                H5T_NATIVE_FLOAT, 4, global, local, offset,
                                state) != 0)) {
    failed = 1;
  }
  if (file >= 0) {
    H5Fclose(file);
  }
  free(state);
  if (patch_consensus(failed)) {
    return -1;
  }
  if (patch_rank == 0) {
    log_info("Wrote compact analysis snapshot %s at %.9g code (%.9g s)",
             filename, time_sim, time_sim * norm_config.Time_Norm);
  }
  return 0;
}

#else

int initialize_analysis_sequence(void) { return -1; }

int dump_analysis_hdf5(void) {
  log_error("Compact analysis output requires the non-orthogonal backend "
            "and a parallel HDF5 build");
  return -1;
}

#endif
