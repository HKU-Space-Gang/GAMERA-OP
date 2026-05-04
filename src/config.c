#include "config.h"

#include <dirent.h>
#include <string.h>
#include <math.h>

#include "log.h"
#include "utils.h"
#include "yaml.h"
#include "curvilinear.h"

FILE *log_file;
int Nt;
double time_stop;
double output_interval;
int restart;
config_t config;
int *proc_coords = config.proc_coords;

int pdm_index_gas = gas_rho_p - gas_rho_h;  // PDM index for gas variables
int pdm_index_gem = mag_bi_p - mag_bi_h;  // PDM index for gem variables

Problem_config_t problem_config;

Norm_t norm_config;

int is, ie, isg, ieg;
int js, je, jsg, jeg;
int ks, ke, ksg, keg;
// time_sim: acumulated simulation time
double dt, dt0, dxx, dyy, dzz, dtout = 0.03, time_sim = 0., time_sim_start = 0.;
double ***x1, ***x2, ***x3;
double ***x1c, ***x2c, ***x3c;
double ***x1ctr, ***x2ctr, ***x3ctr;
double ***dx1, ***dx2, ***dx3;
double ****geo;
double ****rec;
double ****gem;
double *****gas;

int hdf_seq_num = -1;
int log_seq_num = 0;
bool read_restart = true;
// bool read_grid_file = false;
// int read_restart;
// int read_grid_file;
int NF_gem_prim = mags_b3 + 1;
int NF_gas_prim = gas_p_S + 1;
int NF_gem_store = mags_b3_p + 1;
int NF_gas_store = gas_p_S_p + 1;
// Domain boundaries
// The default values are set to 0 and 1, which can be changed in the set_problem_config function in problem files.
double x1min_global = 0.0, x1max_global = 1.0;
double x2min_global = 0.0, x2max_global = 1.0;
double x3min_global = 0.0, x3max_global = 1.0;

double x1min, x1max;
double x2min, x2max;
double x3min, x3max;

int geo_onface_i[KfaceAedgeJ + 1] = {0};
int geo_onface_j[KfaceAedgeJ + 1] = {0};
int geo_onface_k[KfaceAedgeJ + 1] = {0};
int gem_onface_i[magtot_b3 + 1] = {0};
int gem_onface_j[magtot_b3 + 1] = {0};
int gem_onface_k[magtot_b3 + 1] = {0};
int gas_onface_i[bdotv + 1] = {0};
int gas_onface_j[bdotv + 1] = {0};
int gas_onface_k[bdotv + 1] = {0};

// Ringavarage option
// 0 - no ring average; 1 - do ring average
int doRingAverage;

// 0 - not in ring average region; 1 - in ring average region;
int doRingAverage_thisrank_positivepole = 0;
int doRingAverage_thisrank_negativepole = 0;
int Nchunk[32];
int NAVR;

// Physical constants
const double gamma_val = 5.0 / 3.0;
const double CFL = 0.3;
const double PDMB = 1.0;
double CA = 10.0;
const double nlf_cs = 0;
const double nlf_va = 0;
// Set floor value for density and pressure
double rho_floor = 0.00316227; //1e-2.5
double p_floor = 0.000316227; //1e-3.5

// Gravitational field
// Default value is 0. If you want to set gravitational field, you should change
// in file "problem_init.c"
double G01 = 0.0;
double G02 = 0.0;
double G03 = 0.0;

// Divergence of magnetic field
double divB_max = 0.0;

// Factor arrays (Direct initialization)
const double fac8th[8] = {-3.0 / 840, 29.0 / 840, -139.0 / 840, 533.0 / 840,
                          533.0 / 840, -139.0 / 840, 29.0 / 840, -3.0 / 840};
const double facPDMU7[7] = {-1.0 / 140, 5.0 / 84, -101.0 / 420, 319.0 / 420,
                            107.0 / 210, -19.0 / 210, 1.0 / 105};
static int verify_config() {
  if (config.ni_global <= 0 || config.nj_global <= 0 || config.nk_global <= 0) {
    printError("ni_global, nj_global, and nk_global must be positive.\n");
    return -1;
  }
  if (config.proc_dims[0] <= 0 || config.proc_dims[1] <= 0 ||
      config.proc_dims[2] <= 0) {
    printError("proc_dims must be positive.\n");
    return -1;
  }

  if (config.ni_global % config.proc_dims[0] != 0) {
    printError("ni_global must be divisible by ni.\n");
    return -1;
  }
  if (config.nj_global % config.proc_dims[1] != 0) {
    printError("nj_global must be divisible by nj.\n");
    return -1;
  }
  if (config.nk_global % config.proc_dims[2] != 0) {
    printError("nk_global must be divisible by nk.\n");
    return -1;
  }

  return 0;
}

// New function to read the config file
int read_config(const char *filename) {
  // initialize input config to 0
  config.ni_global = 0;
  config.nj_global = 0;
  config.nk_global = 0;
  config.proc_dims[0] = 0;
  config.proc_dims[1] = 0;
  config.proc_dims[2] = 0;

  YamlDocument doc;
  yaml_document_init(&doc);
  if (yaml_parse_file(&doc, filename) != 0) {
    printError("Failed to parse YAML file\n");
    return 1;
  }
  yaml_get_int(&doc, "ni_global", &config.ni_global);
  yaml_get_int(&doc, "nj_global", &config.nj_global);
  yaml_get_int(&doc, "nk_global", &config.nk_global);
  yaml_get_int(&doc, "Nt", &Nt);
  yaml_get_double(&doc, "time_stop", &time_stop);
  yaml_get_int(&doc, "proc_dims_i", &config.proc_dims[0]);
  yaml_get_int(&doc, "proc_dims_j", &config.proc_dims[1]);
  yaml_get_int(&doc, "proc_dims_k", &config.proc_dims[2]);
  yaml_get_double(&doc, "output_interval", &output_interval);
  // Unit normalization
  // Basic Units
  yaml_get_double(&doc, "x_Norm", &norm_config.x_Norm);
  yaml_get_double(&doc, "u_Norm", &norm_config.u_Norm);
  yaml_get_double(&doc, "rho_Norm", &norm_config.rho_Norm);
  //Derived Units
  double mu0 = 4.0 * PI * 1e-7;
  norm_config.Time_Norm = norm_config.x_Norm / norm_config.u_Norm;
  norm_config.p_Norm = norm_config.rho_Norm * norm_config.u_Norm * norm_config.u_Norm;
  norm_config.B_Norm = sqrt(mu0 * norm_config.rho_Norm) * norm_config.u_Norm;
  norm_config.Magnetic_moment_Norm = norm_config.B_Norm * pow(norm_config.x_Norm, 3) / (mu0);
  norm_config.Omega_Norm = 1.0 / norm_config.Time_Norm;
  norm_config.Sigma_Norm = 1.0 / (mu0 * norm_config.u_Norm);

  yaml_document_free(&doc);

  if (verify_config() != 0) {
    printError("Invalid configuration.\n");
    return -1;
  }

  // calculate the derived config variables
  config.ni = config.ni_global / config.proc_dims[0];
  config.nj = config.nj_global / config.proc_dims[1];
  config.nk = config.nk_global / config.proc_dims[2];
  config.NI = config.ni + 2 * NG + 1;
  config.NJ = config.nj + 2 * NG + 1;
  config.NK = config.nk + 2 * NG + 1;

  // initialize ring average config
  if (doRingAverage == 1) {
    init_RingAverage_variables();
  }

  return 0;
}

static void set_boundary_onfaces_geo() {
  // geo_onface_i
  geo_onface_i[edge_idir] = 0;
  geo_onface_i[edge_jdir] = 1;
  geo_onface_i[edge_kdir] = 1;
  geo_onface_i[face_idir] = 1;
  geo_onface_i[face_jdir] = 0;
  geo_onface_i[face_kdir] = 0;
  geo_onface_i[vol_idir] = 1;
  geo_onface_i[vol_jdir] = 0;
  geo_onface_i[vol_kdir] = 0;
  geo_onface_i[vol_center] = 0;
  geo_onface_i[IfaceAedgeJ] = 1;
  geo_onface_i[IfaceAedgeK] = 1;
  geo_onface_i[JfaceAedgeK] = 0;
  geo_onface_i[JfaceAedgeI] = 1;
  geo_onface_i[KfaceAedgeI] = 1;
  geo_onface_i[KfaceAedgeJ] = 0;
  // geo_onface_j
  geo_onface_j[edge_idir] = 1;
  geo_onface_j[edge_jdir] = 0;
  geo_onface_j[edge_kdir] = 1;
  geo_onface_j[face_idir] = 0;
  geo_onface_j[face_jdir] = 1;
  geo_onface_j[face_kdir] = 0;
  geo_onface_j[vol_idir] = 0;
  geo_onface_j[vol_jdir] = 1;
  geo_onface_j[vol_kdir] = 0;
  geo_onface_j[vol_center] = 0;
  geo_onface_j[IfaceAedgeJ] = 1;
  geo_onface_j[IfaceAedgeK] = 0;
  geo_onface_j[JfaceAedgeK] = 1;
  geo_onface_j[JfaceAedgeI] = 1;
  geo_onface_j[KfaceAedgeI] = 0;
  geo_onface_j[KfaceAedgeJ] = 1;
  // geo_onface_k
  geo_onface_k[edge_idir] = 1;
  geo_onface_k[edge_jdir] = 1;
  geo_onface_k[edge_kdir] = 0;
  geo_onface_k[face_idir] = 0;
  geo_onface_k[face_jdir] = 0;
  geo_onface_k[face_kdir] = 1;
  geo_onface_k[vol_idir] = 0;
  geo_onface_k[vol_jdir] = 0;
  geo_onface_k[vol_kdir] = 1;
  geo_onface_k[vol_center] = 0;
  geo_onface_k[IfaceAedgeJ] = 0;
  geo_onface_k[IfaceAedgeK] = 1;
  geo_onface_k[JfaceAedgeK] = 1;
  geo_onface_k[JfaceAedgeI] = 0;
  geo_onface_k[KfaceAedgeI] = 1;
  geo_onface_k[KfaceAedgeJ] = 1;
}

static void set_boundary_onfaces_gem() {
  // gem_onface_i = 1 if var is on the i = i_boundary surface
  gem_onface_i[mag_bi] = 1;
  gem_onface_i[mag0_bi] = 1;
  // gem_onface_j
  gem_onface_j[mag_bj] = 1;
  gem_onface_j[mag0_bj] = 1;
  // gem_onface_k
  gem_onface_k[mag_bk] = 1;
  gem_onface_k[mag0_bk] = 1;
}

static void set_boundary_onfaces_gas() {
  // nothing to do
}

// Function to initialize grid indices
void initialize_grid_indices() {
  isg = 0;
  is = isg + NO2;
  ie = is + config.ni - 1;
  ieg = ie + NO2;
  jsg = 0;
  js = jsg + NO2;
  je = js + config.nj - 1;
  jeg = je + NO2;
  ksg = 0;
  ks = isg + NO2;
  ke = ks + config.nk - 1;
  keg = ke + NO2;
}

int initialize_config(const int rank, const int size, const char *filename) {
  log_info("Reading config file %s for rank %d of %d processes", filename, rank,
           size);
  int ret = read_config(filename);
  if (ret != 0) {
    return -1;
  }

  if (size != config.proc_dims[0] * config.proc_dims[1] * config.proc_dims[2]) {
    printf(
        "The number of processes %d must be equal to the product of the "
        "dimensions %d*%d*%d.\n",
        size, config.proc_dims[0], config.proc_dims[1], config.proc_dims[2]);
    return -1;
  }

  set_boundary_onfaces_geo();
  set_boundary_onfaces_gem();
  set_boundary_onfaces_gas();
  // indices of grid
  // log_info("Initializing grid indices");
  initialize_grid_indices();
  return 0;
}

static int get_last_log_filename(const char *base_name, char *filename) {
  DIR *dir;
  struct dirent *entry;
  char prefix[256];

  snprintf(prefix, sizeof(prefix), out_prefix_pattern, base_name,
           proc_coords[0], proc_coords[1], proc_coords[2]);

  dir = opendir(".");
  if (dir == NULL) {
    log_error("Error opening current directory for log files");
    return -1;
  }

  char last_filename[256] = {0};

  while ((entry = readdir(dir)) != NULL) {
    if (strstr(entry->d_name, prefix) == entry->d_name &&
        strcmp(entry->d_name + strlen(entry->d_name) - 4, ".log") == 0) {
      // Compare to find the lexicographically last filename
      if (strcmp(entry->d_name, last_filename) > 0) {
        strncpy(last_filename, entry->d_name, sizeof(last_filename) - 1);
      }
    }
  }

  closedir(dir);

  if (last_filename[0] == '\0') {
    // log_warn("No log file found with pattern %s*.log", prefix);
    return 1;
  }

  strncpy(filename, last_filename, 256);
  return 0;
}

// Function to initialize logging
int initialize_logging(const char *base_filename, int rank) {
  // Open a log file
  char log_filename[256];
  int log_seq_num = -1;
  if (get_last_log_filename(base_filename, log_filename) == 0) {
    sscanf(log_filename, "%*[^_]_%*[^_]_%d.log", &log_seq_num);
  }
  log_seq_num++;
  get_out_filename(base_filename, log_filename, sizeof(log_filename),
                   log_seq_num, "log");

  log_file = fopen(log_filename, "a");
  if (log_file == NULL) {
    printError("Failed to open log file.\n");
    return -1;
  }
  // Add the log file to the logger, set the log level
  int ret = log_add_fp(log_file, LOG_TRACE);
  if (ret != 0) {
    printError("Failed to add log file.\n");
    return -1;
  }
  // Disable logging to stdout
  log_set_quiet(true);
  log_info("Log initialized for rank %d", rank);
  return 0;
}
