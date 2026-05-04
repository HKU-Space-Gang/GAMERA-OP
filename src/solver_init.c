#include <dirent.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "log.h"

#include "config.h"
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

static int get_last_h5_filename(const char *base_name, char *filename) {
  DIR *dir;
  struct dirent *entry;
  char prefix[256];

  snprintf(prefix, sizeof(prefix), out_prefix_pattern, base_name,
           proc_coords[0], proc_coords[1], proc_coords[2]);

  dir = opendir(".");
  if (dir == NULL) {
    log_error("Error opening current directory for restart files");
    return -1;
  }

  char last_filename[256] = {0};

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
    log_warn("No restart file found with pattern %s*.h5", prefix);
    return 1;
  }

  strncpy(filename, last_filename, 256);
  return 0;
}

static int read_restart_file(const char *base_name) {
  char filename[256];
  int status;
  status = get_last_h5_filename(base_name, filename);
  if (status != 0) {
    return status;
  }
  log_info("Reading restart file: %s", filename);
  status = read_data_hdf5(filename);
  if (status < 0) {
    return -1;
  }
  // parse filename to get hdf_seq_num
  sscanf(filename, "%*[^_]_%*[^_]_%d.h5", &hdf_seq_num);
  log_info("Successfully read restart file %s, hdf_seq_num = %d", filename,
           hdf_seq_num);
  return 0;
}


int initialize_solver(const char *base_name) {
  log_info("Allocating solver");
  if (allocate_solver() != 0) {
    return -1;
  }
  log_info("Trying to read restart file");
  int status = read_restart_file(base_name);
  if (status < 0) {
    return status;
  } else if (status > 0) {
    status = set_domain_decomposition();
    if (status < 0) {
      return -1;
    }
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
