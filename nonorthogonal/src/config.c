#include "config.h"

#include <dirent.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <strings.h>

#include "log.h"
#include "utils.h"
#include "yaml.h"
#include "curvilinear.h"
#include "setup_mpi.h"

FILE *log_file;
int Nt;
double time_stop;
double output_interval;
int analysis_output_enabled;
double restart_interval;
int restart;
config_t config;
int *proc_coords = config.proc_coords;

int pdm_index_gas = gas_rho_p - gas_rho_h;  // PDM index for gas variables
int pdm_index_gem = mag_bi_p - mag_bi_h;  // PDM index for gem variables

Problem_config_t problem_config;

Norm_t norm_config;
Wind_config_t wind_config;
MI_config_t mi_config;
IM_feedback_config_t im_feedback_config;
IM_online_hplus_config_t im_online_hplus_config;
IM_online_hplus_feedback_config_t im_online_hplus_feedback_config;

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
int analysis_seq_num = -1;
int log_seq_num = 0;
bool read_restart = true;
char restart_source_filename[GAMERA_RESTART_PATH_MAX] = {0};
int restart_has_yinyang_hdiv = 0;
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
int ringAverageMode = 1;

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
  if (ringAverageMode != 0 && ringAverageMode != 1) {
    printError("ring_average_mode must be 0 (identity) or 1 (standard).\n");
    return -1;
  }
  if (!(output_interval > 0.0) || !isfinite(output_interval)) {
    printError("output_interval must be finite and positive.\n");
    return -1;
  }
  if (analysis_output_enabled != 0 && analysis_output_enabled != 1) {
    printError("analysis_output_enabled must be 0 or 1.\n");
    return -1;
  }
  if (!isfinite(restart_interval) || restart_interval < 0.0) {
    printError("restart_interval must be finite and nonnegative.\n");
    return -1;
  }
  if (wind_config.enabled &&
      (!(norm_config.x_Norm > 0.0) || !(norm_config.u_Norm > 0.0) ||
       !(norm_config.rho_Norm > 0.0))) {
    printError("x_Norm, u_Norm, and rho_Norm must be positive when wind_file "
               "is used.\n");
    return -1;
  }
  if (wind_config.enabled &&
      (wind_config.enforce_bx_relation < -1 ||
       wind_config.enforce_bx_relation > 1)) {
    printError("wind_enforce_bx_relation must be 0 or 1.\n");
    return -1;
  }
  if (wind_config.enabled &&
      (!(wind_config.velocity_si_scale > 0.0) ||
       !isfinite(wind_config.time_offset) ||
       !isfinite(wind_config.by_coefficient) ||
       !isfinite(wind_config.bz_coefficient) ||
       !isfinite(wind_config.bx_offset) ||
       !isfinite(wind_config.reference[0]) ||
       !isfinite(wind_config.reference[1]) ||
       !isfinite(wind_config.reference[2]) ||
       (wind_config.linear_interpolation != 0 &&
        wind_config.linear_interpolation != 1))) {
    printError("Invalid non-finite solar-wind configuration.\n");
    return -1;
  }
#ifndef GAMERA_TIME_DEPENDENT_WIND
  if (wind_config.enabled) {
    printError("wind_file is unavailable for this problem configuration.\n");
    return -1;
  }
#endif
#ifndef GAMERA_MI_COUPLING
  if (mi_config.enabled) {
    printError("mi_enabled is unavailable for this problem configuration.\n");
    return -1;
  }
#endif
#ifndef GAMERA_NONORTHOGONAL_MATCHED_HPLUS_FEEDBACK
  if (im_feedback_config.enabled) {
    printError("Matched H+ feedback was requested, but this executable was "
               "built without the experimental feedback feature.\n");
    return -1;
  }
#endif
#ifndef GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC
  if (im_online_hplus_config.enabled) {
    printError("Online H+ diagnostics were requested, but this executable "
               "was built without the experimental online provider seam.\n");
    return -1;
  }
#endif
#ifndef GAMERA_NONORTHOGONAL_ONLINE_HPLUS_PRESSURE_FEEDBACK
  if (im_online_hplus_feedback_config.enabled) {
    printError("Live online-H+ pressure feedback was requested, but this "
               "executable was built without that experimental feature.\n");
    return -1;
  }
#endif
  int im_provenance_valid = 1;
  for (size_t byte = 0U; byte < 64U && im_feedback_config.enabled;
       ++byte) {
    const char value = im_feedback_config.provenance_sha256[byte];
    im_provenance_valid &=
        (value >= '0' && value <= '9') || (value >= 'a' && value <= 'f') ||
        (value >= 'A' && value <= 'F');
  }
  if (im_feedback_config.enabled &&
      (im_feedback_config.enabled != 1 ||
       strcmp(im_feedback_config.provider, "prescribed_verification_hdf5") !=
           0 ||
       im_feedback_config.prescribed_state_directory[0] == '\0' ||
       im_feedback_config.prescribed_map_directory[0] == '\0' ||
       strlen(im_feedback_config.provenance_sha256) != 64U ||
       !im_provenance_valid ||
       im_feedback_config.prescribed_end_generation < 0 ||
       !isfinite(im_feedback_config.coupling_interval_s) ||
       im_feedback_config.coupling_interval_s <= 0.0 ||
       !isfinite(im_feedback_config.activation_time_s) ||
       im_feedback_config.activation_time_s < 0.0 ||
       !isfinite(im_feedback_config.relaxation_time_s) ||
       im_feedback_config.relaxation_time_s <= 0.0 ||
       !isfinite(im_feedback_config.matching_time_s) ||
       im_feedback_config.matching_time_s <= 0.0 ||
       !isfinite(im_feedback_config.maximum_fractional_pressure_change) ||
       im_feedback_config.maximum_fractional_pressure_change <= 0.0 ||
       im_feedback_config.maximum_fractional_pressure_change > 1.0 ||
       !isfinite(im_feedback_config.pressure_floor_pa) ||
       im_feedback_config.pressure_floor_pa < 0.0 ||
       !isfinite(im_feedback_config.ion_pressure_fraction) ||
       im_feedback_config.ion_pressure_fraction <= 0.0 ||
       im_feedback_config.ion_pressure_fraction > 1.0 ||
       !isfinite(im_feedback_config.taper_full_l) ||
       im_feedback_config.taper_full_l <= 0.0 ||
       !isfinite(im_feedback_config.taper_zero_l) ||
       im_feedback_config.taper_zero_l <= im_feedback_config.taper_full_l)) {
    printError("Invalid experimental matched-H+ feedback configuration.\n");
    return -1;
  }
  if (im_online_hplus_config.enabled &&
      isfinite(im_online_hplus_config.activation_time_s) &&
      im_online_hplus_config.activation_time_s != 0.0) {
    /* The first production checkpoint schema restores a live provider whose
     * cadence starts at t=0.  Accepting a delayed activation here would
     * create a run that can advance but cannot be restarted. */
    printError("im_hplus_online_activation_time_s must be exactly 0 for the "
               "production checkpoint-v1 path.\n");
    return -1;
  }
  if (im_online_hplus_feedback_config.enabled &&
      (im_online_hplus_feedback_config.enabled != 1 ||
       !im_online_hplus_config.enabled || im_feedback_config.enabled ||
       im_online_hplus_config.activation_time_s != 0.0 ||
       !isfinite(im_online_hplus_feedback_config.relaxation_time_s) ||
       im_online_hplus_feedback_config.relaxation_time_s <= 0.0 ||
       !isfinite(im_online_hplus_feedback_config.matching_time_s) ||
       im_online_hplus_feedback_config.matching_time_s <= 0.0 ||
       !isfinite(im_online_hplus_feedback_config
                     .maximum_fractional_pressure_change) ||
       im_online_hplus_feedback_config.maximum_fractional_pressure_change <=
           0.0 ||
       im_online_hplus_feedback_config.maximum_fractional_pressure_change >
           1.0 ||
       !isfinite(im_online_hplus_feedback_config.pressure_floor_pa) ||
       im_online_hplus_feedback_config.pressure_floor_pa < 0.0 ||
       !isfinite(im_online_hplus_feedback_config.taper_full_l) ||
       im_online_hplus_feedback_config.taper_full_l <= 0.0 ||
       !isfinite(im_online_hplus_feedback_config.taper_zero_l) ||
       im_online_hplus_feedback_config.taper_zero_l <=
           im_online_hplus_feedback_config.taper_full_l)) {
    printError("Invalid live online-H+ pressure-feedback configuration; it "
               "requires the t=0 online provider and cannot coexist with "
               "prescribed feedback.\n");
    return -1;
  }
  const double online_generation =
      im_online_hplus_config.activation_time_s /
      im_online_hplus_config.coupling_interval_s;
  const double online_nearest_generation = nearbyint(online_generation);
  if (im_online_hplus_config.enabled &&
      (im_online_hplus_config.enabled != 1 ||
       im_feedback_config.enabled ||
       strcmp(im_online_hplus_config.provider,
              "online_hplus_callback") != 0 ||
       !isfinite(im_online_hplus_config.coupling_interval_s) ||
       im_online_hplus_config.coupling_interval_s <= 0.0 ||
       !isfinite(im_online_hplus_config.activation_time_s) ||
       !isfinite(online_generation) || online_nearest_generation < 0.0 ||
       fabs(im_online_hplus_config.activation_time_s -
            online_nearest_generation *
                im_online_hplus_config.coupling_interval_s) >
           64.0 * DBL_EPSILON *
               fmax(1.0,
                    fmax(fabs(im_online_hplus_config.activation_time_s),
                         fabs(im_online_hplus_config.coupling_interval_s))) ||
       !isfinite(im_online_hplus_config.cfl_fraction) ||
       im_online_hplus_config.cfl_fraction <= 0.0 ||
       im_online_hplus_config.cfl_fraction > 1.0 ||
       !isfinite(im_online_hplus_config.ion_pressure_fraction) ||
       im_online_hplus_config.ion_pressure_fraction <= 0.0 ||
       im_online_hplus_config.ion_pressure_fraction > 1.0)) {
    printError("Invalid diagnostic-only online H+ configuration.\n");
    return -1;
  }
  if (mi_config.enabled &&
      (mi_config.enabled != 1 ||
       (mi_config.diagnostics_enabled != 0 &&
        mi_config.diagnostics_enabled != 1) ||
       (mi_config.electron_precipitation_enabled != 0 &&
        mi_config.electron_precipitation_enabled != 1) ||
       mi_config.electron_precipitation_hall_model < 0 ||
       mi_config.electron_precipitation_hall_model > 2 ||
       !isfinite(mi_config.electron_precipitation_beta) ||
       !(mi_config.electron_precipitation_beta > 0.0) ||
       !(mi_config.electron_precipitation_ramp_s >= 0.0) ||
       (mi_config.electron_precipitation_conductance_feedback_enabled != 0 &&
        mi_config.electron_precipitation_conductance_feedback_enabled != 1) ||
       !(mi_config.electron_precipitation_conductance_smoothing_s >= 0.0) ||
       !isfinite(mi_config.euv_f107) || !(mi_config.euv_f107 > 0.0) ||
       !isfinite(mi_config.sunward_direction[0]) ||
       !isfinite(mi_config.sunward_direction[1]) ||
       !isfinite(mi_config.sunward_direction[2]) ||
       !(mi_config.sunward_direction[0] * mi_config.sunward_direction[0] +
         mi_config.sunward_direction[1] * mi_config.sunward_direction[1] +
         mi_config.sunward_direction[2] * mi_config.sunward_direction[2] >
         DBL_MIN) ||
       (mi_config.hybrid_dpb_enabled != 0 &&
        mi_config.hybrid_dpb_enabled != 1) ||
       !(mi_config.dpb_absolute_potential_v > 0.0) ||
       !(mi_config.dpb_adaptive_potential_fraction > 0.0) ||
       mi_config.dpb_adaptive_potential_fraction > 1.0 ||
       !(mi_config.dpb_adaptive_potential_floor_v > 0.0) ||
       !(mi_config.dpb_potential_equatorward_offset_deg >= 0.0) ||
       !(mi_config.dpb_fac_equatorward_offset_deg >= 0.0) ||
       !(mi_config.dpb_fac_absolute_floor_a_m2 > 0.0) ||
       !isfinite(mi_config.dpb_minimum_boundary_latitude_deg) ||
       !isfinite(mi_config.dpb_maximum_candidate_latitude_deg) ||
       !isfinite(mi_config.dpb_quiet_initial_nightside_latitude_deg) ||
       !isfinite(mi_config.dpb_oval_center_latitude_deg) ||
       !isfinite(mi_config.dpb_oval_center_mlt_h) ||
       mi_config.dpb_oval_center_mlt_h < 0.0 ||
       mi_config.dpb_oval_center_mlt_h >= 24.0 ||
       !(mi_config.dpb_minimum_boundary_latitude_deg <
         mi_config.dpb_maximum_candidate_latitude_deg) ||
       !(mi_config.dpb_transition_width_deg > 0.0) ||
       !(mi_config.dpb_temporal_timescale_s > 0.0) ||
       !(mi_config.dpb_maximum_slew_deg_per_s > 0.0) ||
       !(mi_config.dpb_minimum_boundary_latitude_deg <
         mi_config.dpb_nightside_poleward_limit_deg) ||
       mi_config.dpb_quiet_initial_nightside_latitude_deg <
           mi_config.dpb_minimum_boundary_latitude_deg ||
       mi_config.dpb_quiet_initial_nightside_latitude_deg >
           mi_config.dpb_nightside_poleward_limit_deg ||
       !(mi_config.dpb_nightside_poleward_limit_deg <
         mi_config.dpb_oval_center_latitude_deg) ||
       !(mi_config.minimum_pedersen_siemens >= 0.0) ||
       !(mi_config.minimum_hall_siemens >= 0.0) ||
       !(mi_config.maximum_hall_to_pedersen > 0.0) ||
       !(mi_config.coupling_interval_s > 0.0) ||
       !(mi_config.diagnostics_interval_s > 0.0) ||
       !(mi_config.ionosphere_radius_re > 0.0) ||
       mi_config.ionosphere_radius_re >= x1min_global ||
       !(mi_config.pedersen_siemens > 0.0) ||
       !(mi_config.hall_siemens >= 0.0) ||
       ((!mi_config.electron_precipitation_enabled ||
         !mi_config.electron_precipitation_conductance_feedback_enabled) &&
        mi_config.hall_siemens != 0.0) || mi_config.current_shell < 0 ||
       mi_config.current_shell >= config.ni_global ||
       !isfinite(mi_config.angular_oversampling) ||
       mi_config.angular_oversampling < 1.0 ||
       mi_config.longitude_count < 4 || mi_config.colatitude_count < 3 ||
       mi_config.maximum_iterations <= 0 ||
       !(mi_config.relative_tolerance > 0.0) ||
       mi_config.absolute_tolerance < 0.0)) {
    printError("Invalid electrostatic M-I precipitation/conductance "
               "configuration.\n");
    return -1;
  }

  return 0;
}

static int read_mi_config(YamlDocument *doc) {
  mi_config = (MI_config_t){
      .enabled = 0,
      .diagnostics_enabled = 1,
      .electron_precipitation_enabled = 0,
      .electron_precipitation_hall_model = 0,
      .electron_precipitation_beta = 0.4362323,
      .electron_precipitation_ramp_s = 600.0,
      .electron_precipitation_conductance_feedback_enabled = 1,
      .electron_precipitation_conductance_smoothing_s = 30.0,
      .euv_f107 = 150.0,
      .sunward_direction = {-1.0, 0.0, 0.0},
      .hybrid_dpb_enabled = 1,
      .dpb_absolute_potential_v = 10000.0,
      .dpb_adaptive_potential_fraction = 0.20,
      .dpb_adaptive_potential_floor_v = 2000.0,
      .dpb_potential_equatorward_offset_deg = 3.0,
      .dpb_fac_equatorward_offset_deg = 2.5,
      .dpb_fac_absolute_floor_a_m2 = 0.03e-6,
      .dpb_minimum_boundary_latitude_deg = 58.0,
      .dpb_maximum_candidate_latitude_deg = 74.0,
      .dpb_transition_width_deg = 3.0,
      .dpb_temporal_timescale_s = 300.0,
      .dpb_maximum_slew_deg_per_s = 1.0 / 120.0,
      .dpb_quiet_initial_nightside_latitude_deg = 64.5,
      .dpb_oval_center_latitude_deg = 85.0,
      .dpb_oval_center_mlt_h = 1.0,
      .dpb_nightside_poleward_limit_deg = 70.0,
      .minimum_pedersen_siemens = 2.0,
      .minimum_hall_siemens = 0.0,
      .maximum_hall_to_pedersen = 6.0,
      .coupling_interval_s = 10.0,
      .diagnostics_interval_s = 10.0,
      .ionosphere_radius_re = 1.0,
      .pedersen_siemens = 5.0,
      .hall_siemens = 0.0,
      .current_shell = 1,
      .angular_oversampling = 1.0,
      .longitude_count = 0,
      .colatitude_count = 0,
      .maximum_iterations = 20000,
      .relative_tolerance = 1.0e-10,
      .absolute_tolerance = 1.0e-8};
  int status = yaml_get_int(doc, "mi_enabled", &mi_config.enabled);
  if (status < 0) {
    return -1;
  }
  const int diagnostics_interval_status = yaml_get_double(
      doc, "mi_diagnostics_interval_s",
      &mi_config.diagnostics_interval_s);
  const int longitude_status =
      yaml_get_int(doc, "mi_longitude_count", &mi_config.longitude_count);
  const int colatitude_status =
      yaml_get_int(doc, "mi_colatitude_count", &mi_config.colatitude_count);
  if (longitude_status < 0 || colatitude_status < 0 ||
      yaml_get_int(doc, "mi_diagnostics_enabled",
                   &mi_config.diagnostics_enabled) < 0 ||
      yaml_get_int(doc, "mi_electron_precipitation_enabled",
                   &mi_config.electron_precipitation_enabled) < 0 ||
      yaml_get_int(doc, "mi_electron_precipitation_hall_model",
                   &mi_config.electron_precipitation_hall_model) < 0 ||
      yaml_get_double(doc, "mi_electron_precipitation_beta",
                      &mi_config.electron_precipitation_beta) < 0 ||
      yaml_get_double(doc, "mi_electron_precipitation_ramp_s",
                      &mi_config.electron_precipitation_ramp_s) < 0 ||
      yaml_get_int(
          doc, "mi_electron_precipitation_conductance_feedback_enabled",
          &mi_config.electron_precipitation_conductance_feedback_enabled) < 0 ||
      yaml_get_double(
          doc, "mi_electron_precipitation_conductance_smoothing_s",
          &mi_config.electron_precipitation_conductance_smoothing_s) < 0 ||
      yaml_get_double(doc, "mi_euv_f107", &mi_config.euv_f107) < 0 ||
      yaml_get_double(doc, "mi_sunward_direction_x",
                      &mi_config.sunward_direction[0]) < 0 ||
      yaml_get_double(doc, "mi_sunward_direction_y",
                      &mi_config.sunward_direction[1]) < 0 ||
      yaml_get_double(doc, "mi_sunward_direction_z",
                      &mi_config.sunward_direction[2]) < 0 ||
      yaml_get_int(doc, "mi_hybrid_dpb_enabled",
                   &mi_config.hybrid_dpb_enabled) < 0 ||
      yaml_get_double(doc, "mi_dpb_absolute_potential_v",
                      &mi_config.dpb_absolute_potential_v) < 0 ||
      yaml_get_double(doc, "mi_dpb_adaptive_potential_fraction",
                      &mi_config.dpb_adaptive_potential_fraction) < 0 ||
      yaml_get_double(doc, "mi_dpb_adaptive_potential_floor_v",
                      &mi_config.dpb_adaptive_potential_floor_v) < 0 ||
      yaml_get_double(doc, "mi_dpb_potential_equatorward_offset_deg",
                      &mi_config.dpb_potential_equatorward_offset_deg) < 0 ||
      yaml_get_double(doc, "mi_dpb_fac_equatorward_offset_deg",
                      &mi_config.dpb_fac_equatorward_offset_deg) < 0 ||
      yaml_get_double(doc, "mi_dpb_fac_absolute_floor_a_m2",
                      &mi_config.dpb_fac_absolute_floor_a_m2) < 0 ||
      yaml_get_double(doc, "mi_dpb_minimum_boundary_latitude_deg",
                      &mi_config.dpb_minimum_boundary_latitude_deg) < 0 ||
      yaml_get_double(doc, "mi_dpb_maximum_candidate_latitude_deg",
                      &mi_config.dpb_maximum_candidate_latitude_deg) < 0 ||
      yaml_get_double(doc, "mi_dpb_transition_width_deg",
                      &mi_config.dpb_transition_width_deg) < 0 ||
      yaml_get_double(doc, "mi_dpb_temporal_timescale_s",
                      &mi_config.dpb_temporal_timescale_s) < 0 ||
      yaml_get_double(doc, "mi_dpb_maximum_slew_deg_per_s",
                      &mi_config.dpb_maximum_slew_deg_per_s) < 0 ||
      yaml_get_double(doc, "mi_dpb_quiet_initial_nightside_latitude_deg",
                      &mi_config.dpb_quiet_initial_nightside_latitude_deg) < 0 ||
      yaml_get_double(doc, "mi_dpb_oval_center_latitude_deg",
                      &mi_config.dpb_oval_center_latitude_deg) < 0 ||
      yaml_get_double(doc, "mi_dpb_oval_center_mlt_h",
                      &mi_config.dpb_oval_center_mlt_h) < 0 ||
      yaml_get_double(doc, "mi_dpb_nightside_poleward_limit_deg",
                      &mi_config.dpb_nightside_poleward_limit_deg) < 0 ||
      yaml_get_double(doc, "mi_minimum_pedersen_siemens",
                      &mi_config.minimum_pedersen_siemens) < 0 ||
      yaml_get_double(doc, "mi_minimum_hall_siemens",
                      &mi_config.minimum_hall_siemens) < 0 ||
      yaml_get_double(doc, "mi_maximum_hall_to_pedersen",
                      &mi_config.maximum_hall_to_pedersen) < 0 ||
      yaml_get_double(doc, "mi_coupling_interval_s",
                      &mi_config.coupling_interval_s) < 0 ||
      diagnostics_interval_status < 0 ||
      yaml_get_double(doc, "mi_ionosphere_radius_re",
                      &mi_config.ionosphere_radius_re) < 0 ||
      yaml_get_double(doc, "mi_pedersen_siemens",
                      &mi_config.pedersen_siemens) < 0 ||
      yaml_get_double(doc, "mi_hall_siemens",
                      &mi_config.hall_siemens) < 0 ||
      yaml_get_int(doc, "mi_current_shell", &mi_config.current_shell) < 0 ||
      yaml_get_double(doc, "mi_angular_oversampling",
                      &mi_config.angular_oversampling) < 0 ||
      yaml_get_int(doc, "mi_maximum_iterations",
                   &mi_config.maximum_iterations) < 0 ||
      yaml_get_double(doc, "mi_relative_tolerance",
                      &mi_config.relative_tolerance) < 0 ||
      yaml_get_double(doc, "mi_absolute_tolerance",
                      &mi_config.absolute_tolerance) < 0) {
    return -1;
  }
  if (mi_config.enabled &&
      (longitude_status > 0 || colatitude_status > 0)) {
    /* Match the polar solver to the nominal Yin--Yang angular cell size.
     * Each patch covers 90 degrees in theta and 270 degrees in phi before
     * the small overset extension.  The ionosphere covers only the polar cap
     * magnetically connected to the inner MHD boundary, so its colatitude
     * count is derived from that cap rather than from a full hemisphere. */
    const double dtheta = 0.5 * PI / (double)config.nj_global;
    const double dphi = 1.5 * PI / (double)config.nk_global;
    const double target_spacing =
        fmin(dtheta, dphi) / mi_config.angular_oversampling;
    if (!(target_spacing > 0.0) || !(x1min_global > 0.0) ||
        !(mi_config.ionosphere_radius_re > 0.0) ||
        mi_config.ionosphere_radius_re >= x1min_global) {
      return -1;
    }
    if (longitude_status > 0) {
      mi_config.longitude_count =
          (int)ceil(2.0 * PI / target_spacing - 1.0e-12);
    }
    if (colatitude_status > 0) {
      const double maximum_colatitude =
          asin(sqrt(mi_config.ionosphere_radius_re / x1min_global));
      mi_config.colatitude_count =
          (int)ceil(maximum_colatitude / target_spacing - 1.0e-12) + 1;
    }
  }
  if (diagnostics_interval_status > 0) {
    /* Keep diagnostic files on the same cadence as the MHD analysis stream
     * unless the user explicitly requests another cadence. */
    mi_config.diagnostics_interval_s = output_interval * norm_config.Time_Norm;
  }
  return 0;
}

static int copy_yaml_text(char *destination, size_t capacity,
                          const char *source) {
  if (destination == NULL || capacity == 0U || source == NULL) {
    return -1;
  }
  size_t begin = 0U;
  size_t length = strlen(source);
  if (length >= 2U &&
      ((source[0] == '"' && source[length - 1U] == '"') ||
       (source[0] == '\'' && source[length - 1U] == '\''))) {
    begin = 1U;
    length -= 2U;
  }
  if (length >= capacity) {
    return -1;
  }
  memcpy(destination, source + begin, length);
  destination[length] = '\0';
  return 0;
}

static int set_path_relative_to_config(char *destination, size_t capacity,
                                       const char *config_filename,
                                       const char *path_value) {
  char path[GAMERA_IM_FEEDBACK_PATH_MAX];
  if (capacity > sizeof(path) ||
      copy_yaml_text(path, sizeof(path), path_value) != 0 ||
      path[0] == '\0') {
    return -1;
  }
  if (path[0] == '/' || config_filename == NULL ||
      strchr(config_filename, '/') == NULL) {
    return copy_yaml_text(destination, capacity, path);
  }
  const char *slash = strrchr(config_filename, '/');
  const size_t directory_length = (size_t)(slash - config_filename);
  const int written = snprintf(destination, capacity, "%.*s/%s",
                               (int)directory_length, config_filename, path);
  return written >= 0 && (size_t)written < capacity ? 0 : -1;
}

static int read_im_feedback_config(YamlDocument *doc,
                                   const char *config_filename) {
  im_feedback_config = (IM_feedback_config_t){
      .enabled = 0,
      .provider = "",
      .prescribed_state_directory = "",
      .prescribed_map_directory = "",
      .provenance_sha256 = "",
      .prescribed_end_generation = -1,
      .coupling_interval_s = 120.0,
      .activation_time_s = 0.0,
      .relaxation_time_s = 3600.0,
      .matching_time_s = 3600.0,
      .maximum_fractional_pressure_change = 0.01,
      .pressure_floor_pa = 0.0,
      .ion_pressure_fraction = 1.0,
      .taper_full_l = 6.0,
      .taper_zero_l = 6.924954264928429};
  if (yaml_get_int(doc, "im_hplus_feedback_enabled",
                   &im_feedback_config.enabled) < 0 ||
      yaml_get_double(doc, "im_hplus_feedback_coupling_interval_s",
                      &im_feedback_config.coupling_interval_s) < 0 ||
      yaml_get_double(doc, "im_hplus_feedback_activation_time_s",
                      &im_feedback_config.activation_time_s) < 0 ||
      yaml_get_double(doc, "im_hplus_feedback_relaxation_time_s",
                      &im_feedback_config.relaxation_time_s) < 0 ||
      yaml_get_double(doc, "im_hplus_feedback_matching_time_s",
                      &im_feedback_config.matching_time_s) < 0 ||
      yaml_get_double(doc,
                      "im_hplus_feedback_max_fractional_pressure_change",
                      &im_feedback_config.maximum_fractional_pressure_change) <
          0 ||
      yaml_get_double(doc, "im_hplus_feedback_pressure_floor_pa",
                      &im_feedback_config.pressure_floor_pa) < 0 ||
      yaml_get_double(doc, "im_hplus_feedback_ion_pressure_fraction",
                      &im_feedback_config.ion_pressure_fraction) < 0 ||
      yaml_get_double(doc, "im_hplus_feedback_taper_full_l",
                      &im_feedback_config.taper_full_l) < 0 ||
      yaml_get_double(doc, "im_hplus_feedback_taper_zero_l",
                      &im_feedback_config.taper_zero_l) < 0) {
    return -1;
  }
  const char *provider =
      yaml_get_string(doc, "im_hplus_feedback_provider");
  if (provider != NULL &&
      copy_yaml_text(im_feedback_config.provider,
                     sizeof(im_feedback_config.provider), provider) != 0) {
    return -1;
  }
  const char *state_directory =
      yaml_get_string(doc, "im_hplus_feedback_state_directory");
  if (state_directory != NULL &&
      set_path_relative_to_config(
          im_feedback_config.prescribed_state_directory,
          sizeof(im_feedback_config.prescribed_state_directory),
          config_filename, state_directory) != 0) {
    return -1;
  }
  const char *provenance =
      yaml_get_string(doc, "im_hplus_feedback_provenance_sha256");
  if (provenance != NULL &&
      copy_yaml_text(im_feedback_config.provenance_sha256,
                     sizeof(im_feedback_config.provenance_sha256),
                     provenance) != 0) {
    return -1;
  }
  const char *map_directory =
      yaml_get_string(doc, "im_hplus_feedback_map_directory");
  if (map_directory != NULL &&
      set_path_relative_to_config(
          im_feedback_config.prescribed_map_directory,
          sizeof(im_feedback_config.prescribed_map_directory),
          config_filename, map_directory) != 0) {
    return -1;
  }
  if (yaml_get_int(doc, "im_hplus_feedback_end_generation",
                   &im_feedback_config.prescribed_end_generation) < 0) {
    return -1;
  }
  return 0;
}

static int read_im_online_hplus_config(YamlDocument *doc) {
  im_online_hplus_config = (IM_online_hplus_config_t){
      .enabled = 0,
      .provider = "",
      .coupling_interval_s = 120.0,
      .activation_time_s = 0.0,
      .cfl_fraction = 0.4,
      .ion_pressure_fraction = 1.0};
  if (yaml_get_int(doc, "im_hplus_online_enabled",
                   &im_online_hplus_config.enabled) < 0 ||
      yaml_get_double(doc, "im_hplus_online_coupling_interval_s",
                      &im_online_hplus_config.coupling_interval_s) < 0 ||
      yaml_get_double(doc, "im_hplus_online_activation_time_s",
                      &im_online_hplus_config.activation_time_s) < 0 ||
      yaml_get_double(doc, "im_hplus_online_cfl_fraction",
                      &im_online_hplus_config.cfl_fraction) < 0 ||
      yaml_get_double(doc, "im_hplus_online_ion_pressure_fraction",
                      &im_online_hplus_config.ion_pressure_fraction) < 0) {
    return -1;
  }
  const char *provider = yaml_get_string(doc, "im_hplus_online_provider");
  return provider == NULL ||
                 copy_yaml_text(im_online_hplus_config.provider,
                                sizeof(im_online_hplus_config.provider),
                                provider) == 0
             ? 0
             : -1;
}

static int read_im_online_hplus_feedback_config(YamlDocument *doc) {
  im_online_hplus_feedback_config = (IM_online_hplus_feedback_config_t){
      .enabled = 0,
      .relaxation_time_s = 3600.0,
      .matching_time_s = 3600.0,
      .maximum_fractional_pressure_change = 0.001,
      .pressure_floor_pa = 0.0,
      .taper_full_l = 6.0,
      .taper_zero_l = 6.924954264928429};
  return yaml_get_int(doc, "im_hplus_online_feedback_enabled",
                      &im_online_hplus_feedback_config.enabled) < 0 ||
                 yaml_get_double(
                     doc, "im_hplus_online_feedback_relaxation_time_s",
                     &im_online_hplus_feedback_config.relaxation_time_s) < 0 ||
                 yaml_get_double(
                     doc, "im_hplus_online_feedback_matching_time_s",
                     &im_online_hplus_feedback_config.matching_time_s) < 0 ||
                 yaml_get_double(
                     doc,
                     "im_hplus_online_feedback_max_fractional_pressure_change",
                     &im_online_hplus_feedback_config
                          .maximum_fractional_pressure_change) < 0 ||
                 yaml_get_double(
                     doc, "im_hplus_online_feedback_pressure_floor_pa",
                     &im_online_hplus_feedback_config.pressure_floor_pa) < 0 ||
                 yaml_get_double(doc,
                                 "im_hplus_online_feedback_taper_full_l",
                                 &im_online_hplus_feedback_config.taper_full_l) <
                     0 ||
                 yaml_get_double(doc,
                                 "im_hplus_online_feedback_taper_zero_l",
                                 &im_online_hplus_feedback_config.taper_zero_l) <
                     0
             ? -1
             : 0;
}

static int set_wind_path(const char *config_filename,
                         const char *wind_filename) {
  char wind_path[GAMERA_WIND_PATH_MAX];
  if (copy_yaml_text(wind_path, sizeof(wind_path), wind_filename) != 0 ||
      wind_path[0] == '\0') {
    return -1;
  }
  if (wind_path[0] == '/' || config_filename == NULL ||
      strchr(config_filename, '/') == NULL) {
    return copy_yaml_text(wind_config.file, sizeof(wind_config.file),
                          wind_path);
  }
  const char *slash = strrchr(config_filename, '/');
  const size_t directory_length = (size_t)(slash - config_filename);
  const int written = snprintf(wind_config.file, sizeof(wind_config.file),
                               "%.*s/%s", (int)directory_length,
                               config_filename, wind_path);
  return written >= 0 && (size_t)written < sizeof(wind_config.file) ? 0 : -1;
}

static int yaml_text_equals(const char *value, const char *expected) {
  char text[64];
  return value != NULL && expected != NULL &&
         copy_yaml_text(text, sizeof(text), value) == 0 &&
         strcasecmp(text, expected) == 0;
}

static int read_wind_config(YamlDocument *doc, const char *config_filename) {
  memset(&wind_config, 0, sizeof(wind_config));
  wind_config.linear_interpolation = 1;
  wind_config.velocity_si_scale = 1.0;
  wind_config.enforce_bx_relation = -1;
  const char *filename = yaml_get_string(doc, "wind_file");
  if (filename == NULL) {
    return 0;
  }
  if (set_wind_path(config_filename, filename) != 0) {
    printError("Invalid or overlong wind_file path.\n");
    return -1;
  }
  wind_config.enabled = 1;

  const char *units = yaml_get_string(doc, "wind_input_units");
  if (units != NULL) {
    if (yaml_text_equals(units, "physical")) {
      wind_config.physical_units = 1;
    } else if (!yaml_text_equals(units, "code")) {
      printError("wind_input_units must be code or physical.\n");
      return -1;
    }
  }
  const char *interpolation =
      yaml_get_string(doc, "wind_interpolation");
  if (interpolation != NULL) {
    if (yaml_text_equals(interpolation, "step")) {
      wind_config.linear_interpolation = 0;
    } else if (!yaml_text_equals(interpolation, "linear")) {
      printError("wind_interpolation must be linear or step.\n");
      return -1;
    }
  }
  const char *velocity_units =
      yaml_get_string(doc, "wind_velocity_units");
  if (velocity_units != NULL) {
    if (yaml_text_equals(velocity_units, "km/s") ||
        yaml_text_equals(velocity_units, "kmps")) {
      wind_config.velocity_si_scale = 1000.0;
    } else if (!yaml_text_equals(velocity_units, "m/s") &&
               !yaml_text_equals(velocity_units, "mps")) {
      printError("wind_velocity_units must be m/s or km/s.\n");
      return -1;
    }
  }
  yaml_get_double(doc, "wind_reference_x", &wind_config.reference[0]);
  yaml_get_double(doc, "wind_reference_y", &wind_config.reference[1]);
  yaml_get_double(doc, "wind_reference_z", &wind_config.reference[2]);
  yaml_get_double(doc, "wind_time_offset", &wind_config.time_offset);
  const int have_by = yaml_get_double(doc, "wind_by_coefficient",
                                      &wind_config.by_coefficient) == 0;
  const int have_bz = yaml_get_double(doc, "wind_bz_coefficient",
                                      &wind_config.bz_coefficient) == 0;
  const int have_bx =
      yaml_get_double(doc, "wind_bx_offset", &wind_config.bx_offset) == 0;
  if (have_by != have_bz) {
    printError("wind_by_coefficient and wind_bz_coefficient must be supplied "
               "together.\n");
    return -1;
  }
  wind_config.coefficient_override = have_by || have_bz || have_bx;
  if (yaml_get_int(doc, "wind_enforce_bx_relation",
                   &wind_config.enforce_bx_relation) < 0) {
    printError("Invalid wind_enforce_bx_relation.\n");
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
  /* Compact, one-file-per-patch analysis is the standard output path.  A
   * configuration may disable routine analysis explicitly, but disabling it
   * must never re-enable the retired per-rank analysis stream. */
  analysis_output_enabled = 1;
  restart_interval = 0.0;
  if (yaml_get_int(&doc, "analysis_output_enabled",
                   &analysis_output_enabled) < 0 ||
      yaml_get_double(&doc, "restart_interval", &restart_interval) < 0) {
    printError("Invalid analysis/restart output configuration.\n");
    yaml_document_free(&doc);
    return -1;
  }
  yaml_get_int(&doc, "ring_average_mode", &ringAverageMode);
#ifdef GAMERA_YINYANG_BACKEND
  /* The historical two-degree extension is less than one angular cell in
   * the production Earth grid.  Permit an explicit overset width so the
   * high-order stencil can be kept several cells away from the ownership
   * switch without changing non-Yin-Yang problems. */
  double yinyang_extension_deg = 0.0;
  const int extension_status = yaml_get_double(
      &doc, "yinyang_angular_extension_deg", &yinyang_extension_deg);
  if (extension_status < 0 ||
      (extension_status == 0 &&
       (!isfinite(yinyang_extension_deg) || yinyang_extension_deg < 0.0 ||
        yinyang_extension_deg >= 45.0))) {
    printError("yinyang_angular_extension_deg must be in [0,45).\n");
    yaml_document_free(&doc);
    return -1;
  }
  if (extension_status == 0) {
    const double extension = yinyang_extension_deg * PI / 180.0;
    x2min_global = 0.25 * PI - extension;
    x2max_global = 0.75 * PI + extension;
    x3min_global = -0.75 * PI - extension;
    x3max_global = 0.75 * PI + extension;
  }
#endif
  // Unit normalization
  // Basic Units
  yaml_get_double(&doc, "x_Norm", &norm_config.x_Norm);
  yaml_get_double(&doc, "u_Norm", &norm_config.u_Norm);
  yaml_get_double(&doc, "rho_Norm", &norm_config.rho_Norm);
  //Derived Units
  double mu0 = 4.0 * PI * 1e-7;
  norm_config.mu0 = mu0;
  norm_config.Time_Norm = norm_config.x_Norm / norm_config.u_Norm;
  norm_config.p_Norm = norm_config.rho_Norm * norm_config.u_Norm * norm_config.u_Norm;
  norm_config.B_Norm = sqrt(mu0 * norm_config.rho_Norm) * norm_config.u_Norm;
  norm_config.Magnetic_moment_Norm = norm_config.B_Norm * pow(norm_config.x_Norm, 3) / (mu0);
  norm_config.Omega_Norm = 1.0 / norm_config.Time_Norm;
  norm_config.Sigma_Norm = 1.0 / (mu0 * norm_config.u_Norm);

  if (read_wind_config(&doc, filename) != 0) {
    yaml_document_free(&doc);
    return -1;
  }
  if (read_mi_config(&doc) != 0) {
    printError("Invalid electrostatic M-I coupling YAML value.\n");
    yaml_document_free(&doc);
    return -1;
  }
  if (read_im_feedback_config(&doc, filename) != 0) {
    printError("Invalid experimental matched-H+ feedback YAML value.\n");
    yaml_document_free(&doc);
    return -1;
  }
  if (read_im_online_hplus_config(&doc) != 0) {
    printError("Invalid diagnostic-only online H+ YAML value.\n");
    yaml_document_free(&doc);
    return -1;
  }
  if (read_im_online_hplus_feedback_config(&doc) != 0) {
    printError("Invalid live online-H+ pressure-feedback YAML value.\n");
    yaml_document_free(&doc);
    return -1;
  }

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

  if (size != patch_count * config.proc_dims[0] * config.proc_dims[1] *
                  config.proc_dims[2]) {
    printf(
        "The number of processes %d must equal patches*proc_dims = "
        "%d*%d*%d*%d.\n",
        size, patch_count, config.proc_dims[0], config.proc_dims[1],
        config.proc_dims[2]);
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

  if (patch_count > 1) {
    snprintf(prefix, sizeof(prefix), "%s_p%d_%02d-%02d-%02d_", base_name,
             patch_id, proc_coords[0], proc_coords[1], proc_coords[2]);
  } else {
    snprintf(prefix, sizeof(prefix), out_prefix_pattern, base_name,
             proc_coords[0], proc_coords[1], proc_coords[2]);
  }

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
    if (patch_count > 1) {
      sscanf(log_filename, "%*[^_]_%*[^_]_%*[^_]_%d.log", &log_seq_num);
    } else {
      sscanf(log_filename, "%*[^_]_%*[^_]_%d.log", &log_seq_num);
    }
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
