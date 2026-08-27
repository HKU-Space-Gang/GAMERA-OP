#include <dirent.h>
#include <errno.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "log.h"

#include "config.h"
#ifdef GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC
#include "nonorthogonal_backend/nonorthogonal_checkpoint_manifest_v2.h"
#include "nonorthogonal_backend/nonorthogonal_inner_magnetosphere_online_checkpoint.h"
#endif
#include "problem.h"
#include "solver.h"
#include "utils.h"
#include "curvilinear.h"
#include "setup_mpi.h"

static int allocate_standalone_arrays() {
  // standalone arrays
  x1 = (double ***)alloc_3d_array(config.NI, config.NJ, config.NK,
                                  sizeof(double));
  if (x1 == NULL) {
    return -1;
  }
  x2 = (double ***)alloc_3d_array(config.NI, config.NJ, config.NK,
                                  sizeof(double));
  if (x2 == NULL) {
    return -1;
  }
  x3 = (double ***)alloc_3d_array(config.NI, config.NJ, config.NK,
                                  sizeof(double));
  if (x3 == NULL) {
    return -1;
  }

  x1c = (double ***)alloc_3d_array(config.NI, config.NJ, config.NK,
                                   sizeof(double));
  if (x1c == NULL) {
    return -1;
  }
  x2c = (double ***)alloc_3d_array(config.NI, config.NJ, config.NK,
                                   sizeof(double));
  if (x2c == NULL) {
    return -1;
  }
  x3c = (double ***)alloc_3d_array(config.NI, config.NJ, config.NK,
                                   sizeof(double));
  if (x3c == NULL) {
    return -1;
  }

  x1ctr = (double ***)alloc_3d_array(config.NI, config.NJ, config.NK,
                                    sizeof(double));
  if (x1ctr == NULL) {
      return -1;
  }
  x2ctr = (double ***)alloc_3d_array(config.NI, config.NJ, config.NK,
                                      sizeof(double));
  if (x2ctr == NULL) {
      return -1;
  }
  x3ctr = (double ***)alloc_3d_array(config.NI, config.NJ, config.NK,
                                      sizeof(double));
  if (x3ctr == NULL) {
      return -1;
  }

  dx1 = (double ***)alloc_3d_array(config.NI, config.NJ, config.NK,
                                   sizeof(double));
  if (dx1 == NULL) {
    return -1;
  }
  dx2 = (double ***)alloc_3d_array(config.NI, config.NJ, config.NK,
                                   sizeof(double));
  if (dx2 == NULL) {
    return -1;
  }
  dx3 = (double ***)alloc_3d_array(config.NI, config.NJ, config.NK,
                                   sizeof(double));
  if (dx3 == NULL) {
    return -1;
  }
  return 0;
}

static int allocate_geometry_arrays() {
  // ------- 4D geometry arrays --------------
  geo = (double ****)alloc_4d_array_contiguous(NF_geo, config.NI, config.NJ,
                                               config.NK, sizeof(double));
  if (geo == NULL) {
    return -1;
  }

  if (doRingAverage == 1) {
    geo_ring = (double ****) alloc_4d_array_contiguous(NF_geo_ring, Ring_Ni, Ring_Nj,
                                                       Ring_Nk, sizeof(double));
    if (geo_ring == NULL) {
      return -1;
    }
  }
  return 0;
}

static int check_dims_from_grid_file() {
  log_info("Checking grid dimensions from grid.h5 if it exists");
  int dims[3];
  // chekc if grid.h5 exists
  if (access("grid.h5", F_OK) == 0) {
    log_info("grid.h5 exists. Checking dimensions.");
    if (read_data_dims("grid.h5", "x1grid", dims, 3) != 0) {
      log_info("Failed to read dimensions from grid.h5");
      return -1;
    }
    int ni_global = dims[2] - 2 * NG - 1;
    int nj_global = dims[1] - 2 * NG - 1;
    int nk_global = dims[0] - 2 * NG - 1;
    log_info("grid.h5 exists. ni_global = %d, nj_global = %d, nk_global = %d",
             ni_global, nj_global, nk_global);
    if (config.ni_global != ni_global) {
      log_error("ni_global from grid.h5 != config.ni_global");
      return -1;
    }
    if (config.nj_global != nj_global) {
      log_error("nj_global from grid.h5 != config.nj_global");
      return -1;
    }
    if (config.nk_global != nk_global) {
      log_error("nk_global from grid.h5 != config.nk_global");
      return -1;
    }
    return 0;
  } else {
    log_info("grid.h5 does not exist. ingore.\n");
    return 1;
  }
}

int read_griddata() {
  double ***x1_temp = (double ***)alloc_3d_array(
      config.nk_global + 2 * NG + 1, config.nj_global + 2 * NG + 1,
      config.ni_global + 2 * NG + 1, sizeof(double));
  if (x1_temp == NULL) {
    return -1;
  }
  double ***x2_temp = (double ***)alloc_3d_array(
      config.nk_global + 2 * NG + 1, config.nj_global + 2 * NG + 1,
      config.ni_global + 2 * NG + 1, sizeof(double));
  if (x2_temp == NULL) {
    return -1;
  }
  double ***x3_temp = (double ***)alloc_3d_array(
      config.nk_global + 2 * NG + 1, config.nj_global + 2 * NG + 1,
      config.ni_global + 2 * NG + 1, sizeof(double));
  if (x3_temp == NULL) {
    return -1;
  }
  log_info("Reading grid data from grid.h5");
  int status = read_grid_data_hdf5("grid.h5", x1_temp, x2_temp, x3_temp);
  if (status < 0) {
    return -1;
  }
  for (int i = isg + proc_coords[0] * config.ni;
       i <= (ieg + 1) + proc_coords[0] * config.ni; i++) {
    for (int j = jsg + proc_coords[1] * config.nj;
         j <= (jeg + 1) + proc_coords[1] * config.nj; j++) {
      for (int k = ksg + proc_coords[2] * config.nk;
           k <= (keg + 1) + proc_coords[2] * config.nk; k++) {
        x1[i - proc_coords[0] * config.ni][j - proc_coords[1] * config.nj]
          [k - proc_coords[2] * config.nk] = x1_temp[k][j][i];
        x2[i - proc_coords[0] * config.ni][j - proc_coords[1] * config.nj]
          [k - proc_coords[2] * config.nk] = x2_temp[k][j][i];
        x3[i - proc_coords[0] * config.ni][j - proc_coords[1] * config.nj]
          [k - proc_coords[2] * config.nk] = x3_temp[k][j][i];
      }
    }
  }

  x1min_global = x1_temp[0][0][NG];
  x1max_global = x1_temp[0][0][config.ni_global + NG];
  x2min_global = x2_temp[0][NG][0];
  x2max_global = x2_temp[0][config.nj_global + NG][0];
  x3min_global = x3_temp[NG][0][0];
  x3max_global = x3_temp[config.nk_global + NG][0][0];
  x1min = x1[is][js][ks];
  x1max = x1[ie+1][je][ke];
  x2min = x2[is][js][ks];
  x2max = x2[ie][je+1][ke];
  x3min = x3[is][js][ks];
  x3max = x3[ie][je][ke+1];
  log_info("x1min_global = %f, x1max_global = %f, x2min_global = %f, x2max_global = %f, x3min_global = %f, x3max_global = %f",
           x1min_global, x1max_global, x2min_global, x2max_global,
           x3min_global, x3max_global);
  log_info("x1min = %f, x1max = %f, x2min = %f, x2max = %f, x3min = %f, x3max = %f for rank %d",
           x1min, x1max, x2min, x2max, x3min, x3max, rank);

  free_3d_array((void ***)x1_temp);
  free_3d_array((void ***)x2_temp);
  free_3d_array((void ***)x3_temp);

  return 0;
}
void default_decomposition() {
  double dx1_proc = (x1max_global - x1min_global) / config.proc_dims[0];
  double dx2_proc = (x2max_global - x2min_global) / config.proc_dims[1];
  double dx3_proc = (x3max_global - x3min_global) / config.proc_dims[2];
  x1min = x1min_global + dx1_proc * proc_coords[0];
  x1max = x1min_global + dx1_proc * (proc_coords[0] + 1);
  x2min = x2min_global + dx2_proc * proc_coords[1];
  x2max = x2min_global + dx2_proc * (proc_coords[1] + 1);
  x3min = x3min_global + dx3_proc * proc_coords[2];
  x3max = x3min_global + dx3_proc * (proc_coords[2] + 1);

  // Setting up the grid and initial conditions
  for (int i = isg; i <= (ieg + 1); i++) {
    for (int j = jsg; j <= (jeg + 1); j++) {
      for (int k = ksg; k <= (keg + 1); k++) {
        x1[i][j][k] = x1min - (x1max - x1min) * NO2 / config.ni +
                      (i - isg) * (x1max - x1min) / config.ni;
        x2[i][j][k] = x2min - (x2max - x2min) * NO2 / config.nj +
                      (j - jsg) * (x2max - x2min) / config.nj;
        x3[i][j][k] = x3min - (x3max - x3min) * NO2 / config.nk +
                      (k - ksg) * (x3max - x3min) / config.nk;
      }
    }
  }
}

int set_domain_decomposition() {
  int status = check_dims_from_grid_file();
  if (status == 0) {
    return read_griddata();
  } else if (status < 0) {
    return -1;
  } else {
    log_info(
        "Neither data file nor grid file provided. Use default decomposition.");
    default_decomposition();
  }
  log_info("Successfully decomposed domain");
  return 0;
}


static int allocate_data_arrays() {
  // allocate memory for the Ring average related arrays
  if (doRingAverage == 1) {
    gem_ring = (double ****) alloc_4d_array_contiguous(NF_gem_ring, Ring_Ni, Ring_Nj,
                                                       Ring_Nk, sizeof(double));
    if (gem_ring == NULL) {
      return -1;
    }
    gas_ring = (double *****) alloc_5d_array_with_4d_contiguous(NS1, NF_gas_ring, Ring_Ni, Ring_Nj,
                                                                Ring_Nk, sizeof(double));
    if (gas_ring == NULL) {
      return -1;
    }
  }

  // allocate memory for the grid
  gem = (double ****)alloc_4d_array_contiguous(NF_gem, config.NI, config.NJ,
                                               config.NK, sizeof(double));
  if (gem == NULL) {
    return -1;
  }
  // single gas NS1 = NS =1;
  // for Multi-gas MHD: NS1 = NS+1;
  // s = 0 -> gas 1, s = 1 -> gas 2, etc. s = NF -> bulk values
  gas = (double *****)alloc_5d_array_with_4d_contiguous(
      NS1, NF_gas, config.NI, config.NJ, config.NK, sizeof(double));
  if (gas == NULL) {
    return -1;
  }

  rec = (double ****)alloc_4d_array_contiguous(NF_rec, config.NI,
                                               config.NJ, config.NK, sizeof(double));
  if (rec == NULL) {
    return -1;
  }

  return 0;
}

int allocate_solver() {
  log_info("Initializing grid");
  if (allocate_standalone_arrays() != 0) {
    return -1;
  }
  log_info("Allocating geometry arrays");
  if (allocate_geometry_arrays() != 0) {
    return -1;
  }
  log_info("Allocating data arrays");
  if (allocate_data_arrays() != 0) {
    return -1;
  }
  return 0;
}

static int checkpoint_manifest_complete(const char *directory) {
  char path[768];
  if (snprintf(path, sizeof(path), "%s/manifest.json", directory) < 0) {
    return 0;
  }
  FILE *file = fopen(path, "r");
  if (file == NULL) {
    return 0;
  }
  char text[1024];
  const size_t count = fread(text, 1, sizeof(text) - 1, file);
  fclose(file);
  text[count] = '\0';
  return strstr(text, "\"complete\": true") != NULL;
}

#ifdef GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC
static int parse_checkpoint_sequence(const char *name, uint64_t *sequence) {
  static const char prefix[] = "checkpoint_";
  if (name == NULL || sequence == NULL ||
      strncmp(name, prefix, sizeof(prefix) - 1U) != 0) {
    return 0;
  }
  const char *digits = name + sizeof(prefix) - 1U;
  if (*digits == '\0') {
    return 0;
  }
  for (const char *cursor = digits; *cursor != '\0'; ++cursor) {
    if (*cursor < '0' || *cursor > '9') {
      return 0;
    }
  }
  errno = 0;
  char *end = NULL;
  const unsigned long long parsed = strtoull(digits, &end, 10);
  if (errno != 0 || end == digits || *end != '\0' ||
      parsed > (unsigned long long)INT_MAX) {
    return 0;
  }
  *sequence = (uint64_t)parsed;
  return 1;
}

/* Online H+ restart discovery is one root-selected strict transaction.  A
 * directory without a manifest is merely incomplete and may be skipped when
 * some other complete checkpoint exists.  If checkpoint directories exist
 * but none has a complete manifest, fail stop rather than silently starting a
 * fresh run over interrupted restart state.  An existing malformed or
 * incompatible manifest is likewise fail-stop and is never hidden by an
 * older checkpoint or a legacy root-level rank file. */
static int get_online_checkpoint_filename(const char *rank_prefix,
                                          char *filename,
                                          size_t filename_size) {
  int world_size = 0;
  int status = 0;
  uint64_t selected_sequence = 0U;
  int selected = 0;
  int checkpoint_directory_seen = 0;
  if (rank_prefix == NULL || filename == NULL || filename_size == 0U ||
      MPI_Comm_size(MPI_COMM_WORLD, &world_size) != MPI_SUCCESS ||
      world_size <= 0) {
    status = -1;
  }
  if (rank == 0 && status == 0) {
    DIR *root = opendir("restart");
    if (root == NULL) {
      status = errno == ENOENT ? 1 : -1;
    } else {
      struct dirent *entry = NULL;
      while (status == 0 && (entry = readdir(root)) != NULL) {
        static const char checkpoint_prefix[] = "checkpoint_";
        if (strncmp(entry->d_name, checkpoint_prefix,
                    sizeof(checkpoint_prefix) - 1U) != 0) {
          continue;
        }
        char directory[512];
        const int directory_count =
            snprintf(directory, sizeof(directory), "restart/%s",
                     entry->d_name);
        if (directory_count < 0 ||
            (size_t)directory_count >= sizeof(directory)) {
          status = -1;
          break;
        }
        struct stat directory_information;
        if (stat(directory, &directory_information) != 0) {
          status = -1;
          break;
        }
        if (!S_ISDIR(directory_information.st_mode)) {
          continue;
        }
        checkpoint_directory_seen = 1;
        uint64_t sequence = 0U;
        if (!parse_checkpoint_sequence(entry->d_name, &sequence)) {
          status = -1;
          break;
        }
        char manifest_path[768];
        const int manifest_count =
            snprintf(manifest_path, sizeof(manifest_path), "%s/manifest.json",
                     directory);
        if (manifest_count < 0 ||
            (size_t)manifest_count >= sizeof(manifest_path)) {
          status = -1;
          break;
        }
        struct stat manifest_information;
        if (stat(manifest_path, &manifest_information) != 0) {
          if (errno == ENOENT) {
            continue;
          }
          status = -1;
          break;
        }
        gamera_no_checkpoint_manifest_v2 manifest;
        memset(&manifest, 0, sizeof(manifest));
        if (!S_ISREG(manifest_information.st_mode) ||
            gamera_no_checkpoint_manifest_v2_read_strict(manifest_path,
                                                         &manifest) != 0 ||
            manifest.checkpoint_sequence != sequence ||
            manifest.rank_file_count != (uint64_t)world_size ||
            strcmp(manifest.provider_sidecar_filename,
                   GAMERA_NO_IM_ONLINE_CHECKPOINT_SIDECAR) != 0 ||
            manifest.provider_sidecar_schema !=
                gamera_no_im_online_checkpoint_expected_schema()) {
          status = -1;
          break;
        }
        char sidecar_path[768];
        const int sidecar_count =
            snprintf(sidecar_path, sizeof(sidecar_path), "%s/%s", directory,
                     manifest.provider_sidecar_filename);
        struct stat sidecar_information;
        if (sidecar_count < 0 ||
            (size_t)sidecar_count >= sizeof(sidecar_path) ||
            stat(sidecar_path, &sidecar_information) != 0 ||
            !S_ISREG(sidecar_information.st_mode) ||
            sidecar_information.st_size <= 0 ||
            (uint64_t)sidecar_information.st_size !=
                manifest.provider_sidecar_file_size) {
          status = -1;
          break;
        }
        if (!selected || sequence > selected_sequence) {
          selected = 1;
          selected_sequence = sequence;
        }
      }
      if (closedir(root) != 0 && status == 0) {
        status = -1;
      }
      if (status == 0 && !selected) {
        if (checkpoint_directory_seen) {
          log_error("Online H+ restart contains checkpoint directories but "
                    "no complete compatible online-H+ manifest; refusing "
                    "fresh fallback");
          status = -1;
        } else {
          status = 1;
        }
      }
    }
  }
  if (MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&selected_sequence, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD) !=
          MPI_SUCCESS) {
    return -1;
  }
  if (status != 0) {
    return status;
  }
  const int count = snprintf(
      filename, filename_size, "restart/checkpoint_%06llu/%s%06llu.h5",
      (unsigned long long)selected_sequence, rank_prefix,
      (unsigned long long)selected_sequence);
  const int local_failed = count < 0 || (size_t)count >= filename_size ||
                           access(filename, R_OK) != 0;
  int global_failed = 0;
  if (MPI_Allreduce(&local_failed, &global_failed, 1, MPI_INT, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS ||
      global_failed) {
    return -1;
  }
  return 0;
}
#endif

static int get_last_checkpoint_filename(const char *prefix, char *filename,
                                        size_t filename_size) {
  DIR *root = opendir("restart");
  if (root == NULL) {
    return errno == ENOENT ? 1 : -1;
  }
  struct dirent *entry;
  char selected[256] = {0};
  while ((entry = readdir(root)) != NULL) {
    if (strncmp(entry->d_name, "checkpoint_", 11) != 0 ||
        (selected[0] != '\0' && strcmp(entry->d_name, selected) <= 0)) {
      continue;
    }
    char directory[512];
    char candidate[768];
    if (snprintf(directory, sizeof(directory), "restart/%s", entry->d_name) < 0 ||
        !checkpoint_manifest_complete(directory) ||
        snprintf(candidate, sizeof(candidate), "%s/%s%s.h5", directory,
                 prefix, entry->d_name + 11) < 0 || access(candidate, R_OK) != 0) {
      continue;
    }
    strncpy(selected, entry->d_name, sizeof(selected) - 1);
  }
  closedir(root);
  if (selected[0] == '\0') {
    return 1;
  }
  if (snprintf(filename, filename_size, "restart/%s/%s%s.h5", selected,
               prefix, selected + 11) < 0) {
    return -1;
  }
  return 0;
}

static int get_last_h5_filename(const char *base_name, char *filename,
                                size_t filename_size) {
  DIR *dir;
  struct dirent *entry;
  char prefix[256];

  if (patch_count > 1) {
    snprintf(prefix, sizeof(prefix), "%s_p%d_%02d-%02d-%02d_", base_name,
             patch_id, proc_coords[0], proc_coords[1], proc_coords[2]);
  } else {
    snprintf(prefix, sizeof(prefix), out_prefix_pattern, base_name,
             proc_coords[0], proc_coords[1], proc_coords[2]);
  }

#ifdef GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC
  if (im_online_hplus_config.enabled &&
      strcmp(base_name, BASE_RESTART_NAME) == 0) {
    return get_online_checkpoint_filename(prefix, filename, filename_size);
  }
#endif

  dir = opendir(".");
  if (dir == NULL) {
    log_error("Error opening current directory for restart files");
    return -1;
  }

  char last_filename[512] = {0};

  while ((entry = readdir(dir)) != NULL) {
    if (strstr(entry->d_name, prefix) == entry->d_name &&
        strcmp(entry->d_name + strlen(entry->d_name) - 3, ".h5") == 0) {
      // Compare to find the lexicographically last filename
      if (strcmp(entry->d_name, last_filename) > 0) {
        strncpy(last_filename, entry->d_name, sizeof(last_filename) - 1);
      }
    }
  }

  closedir(dir);

  if (last_filename[0] == '\0') {
    if (strcmp(base_name, BASE_RESTART_NAME) == 0) {
      const int checkpoint_status =
          get_last_checkpoint_filename(prefix, filename, filename_size);
      if (checkpoint_status <= 0) {
        return checkpoint_status;
      }
    }
    log_warn("No restart file found with pattern %s*.h5", prefix);
    return 1;
  }

  strncpy(filename, last_filename, filename_size - 1);
  filename[filename_size - 1] = '\0';
  return 0;
}

static int read_restart_file(const char *base_name) {
  char filename[768];
  int status;
  status = get_last_h5_filename(base_name, filename, sizeof(filename));
  if (status != 0) {
    return status;
  }
  log_info("Reading restart file: %s", filename);
  status = read_data_hdf5(filename);
  if (status < 0) {
    return -1;
  }
  const int restart_path_count = snprintf(
      restart_source_filename, sizeof(restart_source_filename), "%s",
      filename);
  if (restart_path_count < 0 ||
      (size_t)restart_path_count >= sizeof(restart_source_filename)) {
    log_error("Restart source path is too long: %s", filename);
    restart_source_filename[0] = '\0';
    return -1;
  }
  // parse filename to get hdf_seq_num
  const char *sequence = strrchr(filename, '_');
  if (sequence == NULL || sscanf(sequence + 1, "%d.h5", &hdf_seq_num) != 1) {
    log_error("Unable to parse restart sequence from %s", filename);
    return -1;
  }
  log_info("Successfully read restart file %s, hdf_seq_num = %d", filename,
           hdf_seq_num);
  return 0;
}

/*
 * Restart discovery and HDF reads are rank-local, but the Yin-Yang H(div)
 * path contains world collectives.  Reject a mixed checkpoint generation or
 * scheme here, before one subset can skip reconciliation while another
 * subset enters MPI_Allgather.
 */
static int validate_restart_consensus(int local_status) {
  int status_min = local_status;
  int status_max = local_status;
  if (MPI_Allreduce(&local_status, &status_min, 1, MPI_INT, MPI_MIN,
                    MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Allreduce(&local_status, &status_max, 1, MPI_INT, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS) {
    log_error("Unable to establish restart-file consensus");
    return -1;
  }
  if (status_min != status_max) {
    log_error("Inconsistent restart availability/read status across ranks: "
              "min=%d max=%d",
              status_min, status_max);
    return -1;
  }
  if (status_min != 0) {
    return status_min;
  }

  const int local_metadata[3] = {hdf_seq_num, log_seq_num,
                                 restart_has_yinyang_hdiv};
  int metadata_min[3];
  int metadata_max[3];
  double time_min = time_sim;
  double time_max = time_sim;
  double dt0_min = dt0;
  double dt0_max = dt0;
  const int local_dt0_invalid = !isfinite(dt0) || dt0 <= 0.0;
  int global_dt0_invalid = local_dt0_invalid;
  if (MPI_Allreduce(local_metadata, metadata_min, 3, MPI_INT, MPI_MIN,
                    MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Allreduce(local_metadata, metadata_max, 3, MPI_INT, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Allreduce(&time_sim, &time_min, 1, MPI_DOUBLE, MPI_MIN,
                    MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Allreduce(&time_sim, &time_max, 1, MPI_DOUBLE, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Allreduce(&dt0, &dt0_min, 1, MPI_DOUBLE, MPI_MIN,
                    MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Allreduce(&dt0, &dt0_max, 1, MPI_DOUBLE, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Allreduce(&local_dt0_invalid, &global_dt0_invalid, 1, MPI_INT,
                    MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS) {
    log_error("Unable to validate restart metadata across ranks");
    return -1;
  }
  for (int item = 0; item < 3; ++item) {
    if (metadata_min[item] != metadata_max[item]) {
      log_error("Inconsistent restart metadata field %d across ranks: "
                "min=%d max=%d",
                item, metadata_min[item], metadata_max[item]);
      return -1;
    }
  }
  if (!isfinite(time_min) || !isfinite(time_max) || time_min != time_max) {
    log_error("Inconsistent restart time across ranks: min=%.17g max=%.17g",
              time_min, time_max);
    return -1;
  }
  if (global_dt0_invalid != 0 || !isfinite(dt0_min) || !isfinite(dt0_max) ||
      dt0_min <= 0.0 || dt0_min != dt0_max) {
    log_error("Invalid or inconsistent restart dt0 across ranks: min=%.17g "
              "max=%.17g",
              dt0_min, dt0_max);
    return -1;
  }
  return 0;
}


int initialize_solver(const char *base_name) {
  log_info("Allocating solver");
  if (allocate_solver() != 0) {
    return -1;
  }
  log_info("Trying to read restart file");
  /* New compact-analysis runs keep restart checkpoints under a distinct
   * prefix.  Fall back to the historical mhd prefix so existing runs remain
   * restartable without conversion. */
  int status = read_restart_file(BASE_RESTART_NAME);
  if (status > 0 && !im_online_hplus_config.enabled) {
    status = read_restart_file(base_name);
  }
  status = validate_restart_consensus(status);
  if (status < 0) {
    return status;
  } else if (status > 0) {
    restart_source_filename[0] = '\0';
    status = set_domain_decomposition();
    if (status < 0) {
      return -1;
    }
#ifdef GAMERA_NONORTHOGONAL_BACKEND
    status = problem_grid_init();
    if (status != 0) {
      log_error("Error generating non-orthogonal problem grid");
      return -1;
    }
#endif
    // start from scratch
    read_restart = false;
    time_sim = 0.0;
    time_sim_start = 0.0;
  }
//  set_domain_decomposition();
  log_info("Calculating geometry arrays");
  if (set_geometry_arrays() != 0) {
    log_error("Error calculating geometry arrays");
    return -1;
  }


  log_info("Solver initialized");
  return 0;
}

int finalize_solver() {
  log_info("Finalizing solver");
  log_info("Freeing 3D arrays");
  if (x1) {
    free_3d_array((void ***)x1);
  }
  if (x2) {
    free_3d_array((void ***)x2);
  }
  if (x3) {
    free_3d_array((void ***)x3);
  }
  if (x1c) {
    free_3d_array((void ***)x1c);
  }
  if (x2c) {
    free_3d_array((void ***)x2c);
  }
  if (x3c) {
    free_3d_array((void ***)x3c);
  }
  if (x1ctr) {
    free_3d_array((void ***)x1ctr);
  }
  if (x2ctr) {
      free_3d_array((void ***)x2ctr);
  }
  if (x3ctr) {
      free_3d_array((void ***)x3ctr);
  }
  if (dx1) {
    free_3d_array((void ***)dx1);
  }
  if (dx2) {
    free_3d_array((void ***)dx2);
  }
  if (dx3) {
    free_3d_array((void ***)dx3);
  }

  if (doRingAverage == 1) {
    log_info("Freeing Ring average arrays");
    if (geo_ring) {
      free_4d_array_contiguous((void ****) geo_ring, NF_geo_ring);
    }
    if (gem_ring) {
      free_4d_array_contiguous((void ****) gem_ring, NF_gem_ring);
    }
    if (gas_ring) {
      free_5d_array_with_4d_contiguous((void *****) gas_ring, NS1, NF_gas_ring);
    }
  }

  log_info("Freeing 4D array geo");
  if (geo) {
    free_4d_array_contiguous((void ****)geo, NF_geo);
  }
  log_info("Freeing 4D array gem");
  if (gem) {
      free_4d_array_contiguous((void ****)gem, NF_gem);
  }
  log_info("Freeing 5D array gas");
  if (gas) {
      free_5d_array_with_4d_contiguous((void *****)gas, NS1, NF_gas);
  }
  log_info("Solver finalized");
  return 0;
}
