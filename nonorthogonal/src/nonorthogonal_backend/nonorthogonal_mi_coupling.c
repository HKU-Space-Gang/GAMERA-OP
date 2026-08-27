#include "nonorthogonal_mi_coupling.h"

#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifndef GAMERA_MI_SNAPSHOT_TESTING
#include "nonorthogonal_mi_collective.h"
#include "nonorthogonal_mi.h"
#include "nonorthogonal_mi_dpb.h"
#include "nonorthogonal_electron_precipitation.h"
#include "nonorthogonal_legacy_adapter.h"
#include "nonorthogonal_mi_restart.h"
#include "nonorthogonal_state.h"
#include "nonorthogonal_yinyang.h"

#include "config.h"
#include "log.h"
#include "problem.h"
#include "setup_mpi.h"

#include <hdf5.h>
#include <mpi.h>
#include <stdio.h>
#endif

typedef struct {
  double *colatitude_rad;
  double *longitude_rad;
  double *north_potential_v;
  size_t theta_points;
  size_t azimuth_points;
  double ionosphere_radius_re;
  double trace_radius_re;
  double epoch_s;
  uint64_t generation;
  int available;
  int ready;
} mi_north_snapshot_state_t;

static mi_north_snapshot_state_t mi_north_snapshot;

static void mi_north_snapshot_clear(void) {
  free(mi_north_snapshot.colatitude_rad);
  free(mi_north_snapshot.longitude_rad);
  free(mi_north_snapshot.north_potential_v);
  memset(&mi_north_snapshot, 0, sizeof(mi_north_snapshot));
}

static int mi_north_snapshot_allocate(size_t theta_points,
                                      size_t azimuth_points) {
  size_t potential_count;
  mi_north_snapshot_clear();
  if (theta_points < 2U || azimuth_points < 2U ||
      theta_points > SIZE_MAX / azimuth_points) {
    return -1;
  }
  potential_count = theta_points * azimuth_points;
  if (theta_points > SIZE_MAX / sizeof(double) ||
      azimuth_points > SIZE_MAX / sizeof(double) ||
      potential_count > SIZE_MAX / sizeof(double)) {
    return -1;
  }
  mi_north_snapshot.colatitude_rad =
      (double *)calloc(theta_points, sizeof(double));
  mi_north_snapshot.longitude_rad =
      (double *)calloc(azimuth_points, sizeof(double));
  mi_north_snapshot.north_potential_v =
      (double *)calloc(potential_count, sizeof(double));
  if (mi_north_snapshot.colatitude_rad == NULL ||
      mi_north_snapshot.longitude_rad == NULL ||
      mi_north_snapshot.north_potential_v == NULL) {
    mi_north_snapshot_clear();
    return -1;
  }
  mi_north_snapshot.theta_points = theta_points;
  mi_north_snapshot.azimuth_points = azimuth_points;
  return 0;
}

static int mi_north_snapshot_configure(double maximum_colatitude_rad,
                                       double ionosphere_radius_re,
                                       double trace_radius_re) {
  if (mi_north_snapshot.colatitude_rad == NULL ||
      mi_north_snapshot.longitude_rad == NULL ||
      mi_north_snapshot.north_potential_v == NULL ||
      !isfinite(maximum_colatitude_rad) ||
      !(maximum_colatitude_rad > 0.0) ||
      !isfinite(ionosphere_radius_re) || !(ionosphere_radius_re > 0.0) ||
      !isfinite(trace_radius_re) || !(trace_radius_re > 0.0)) {
    return -1;
  }
  const double pi = acos(-1.0);
  for (size_t theta = 0U; theta < mi_north_snapshot.theta_points; ++theta) {
    mi_north_snapshot.colatitude_rad[theta] =
        maximum_colatitude_rad * (double)theta /
        (double)(mi_north_snapshot.theta_points - 1U);
  }
  for (size_t azimuth = 0U; azimuth < mi_north_snapshot.azimuth_points;
       ++azimuth) {
    mi_north_snapshot.longitude_rad[azimuth] =
        2.0 * pi * (double)azimuth /
        (double)mi_north_snapshot.azimuth_points;
  }
  mi_north_snapshot.ionosphere_radius_re = ionosphere_radius_re;
  mi_north_snapshot.trace_radius_re = trace_radius_re;
  mi_north_snapshot.available = 1;
  mi_north_snapshot.ready = 0;
  return 0;
}

/* Publication is transactional from a reader's perspective: validate the
 * complete candidate first, copy it into the persistent borrowed-view
 * buffer, then publish epoch/generation and ready last.  The copy happens
 * once per successful M-I update; borrowing the view is zero-copy. */
static int mi_north_snapshot_publish(double epoch_s,
                                     const double *north_potential_v) {
  if (!mi_north_snapshot.available || north_potential_v == NULL ||
      !isfinite(epoch_s) || mi_north_snapshot.generation == UINT64_MAX ||
      mi_north_snapshot.theta_points >
          SIZE_MAX / mi_north_snapshot.azimuth_points) {
    return -1;
  }
  const size_t count = mi_north_snapshot.theta_points *
                       mi_north_snapshot.azimuth_points;
  for (size_t index = 0U; index < count; ++index) {
    if (!isfinite(north_potential_v[index])) {
      return -1;
    }
  }
  memcpy(mi_north_snapshot.north_potential_v, north_potential_v,
         count * sizeof(double));
  mi_north_snapshot.epoch_s = epoch_s;
  ++mi_north_snapshot.generation;
  mi_north_snapshot.ready = 1;
  return 0;
}

gamera_mi_snapshot_status gamera_mi_coupling_borrow_north_snapshot(
    uint64_t expected_generation, gamera_mi_north_snapshot_view *view) {
  if (view == NULL) {
    return GAMERA_MI_SNAPSHOT_INVALID_ARGUMENT;
  }
  memset(view, 0, sizeof(*view));
  if (!mi_north_snapshot.available) {
    return GAMERA_MI_SNAPSHOT_UNAVAILABLE;
  }
  if (!mi_north_snapshot.ready) {
    return GAMERA_MI_SNAPSHOT_UNREADY;
  }
  if (expected_generation != 0U &&
      expected_generation != mi_north_snapshot.generation) {
    return GAMERA_MI_SNAPSHOT_STALE;
  }
  view->ready = 1;
  view->epoch_s = mi_north_snapshot.epoch_s;
  view->generation = mi_north_snapshot.generation;
  view->theta_points = mi_north_snapshot.theta_points;
  view->azimuth_points = mi_north_snapshot.azimuth_points;
  view->ionosphere_radius_re = mi_north_snapshot.ionosphere_radius_re;
  view->trace_radius_re = mi_north_snapshot.trace_radius_re;
  view->colatitude_rad = mi_north_snapshot.colatitude_rad;
  view->longitude_rad = mi_north_snapshot.longitude_rad;
  view->potential_v = mi_north_snapshot.north_potential_v;
  return GAMERA_MI_SNAPSHOT_OK;
}

const char *gamera_mi_snapshot_status_string(
    gamera_mi_snapshot_status status) {
  switch (status) {
    case GAMERA_MI_SNAPSHOT_OK:
      return "ok";
    case GAMERA_MI_SNAPSHOT_UNAVAILABLE:
      return "unavailable";
    case GAMERA_MI_SNAPSHOT_UNREADY:
      return "unready";
    case GAMERA_MI_SNAPSHOT_STALE:
      return "stale";
    case GAMERA_MI_SNAPSHOT_INVALID_ARGUMENT:
      return "invalid argument";
  }
  return "unknown";
}

#ifdef GAMERA_MI_SNAPSHOT_TESTING
int gamera_mi_snapshot_test_prepare(size_t theta_points,
                                    size_t azimuth_points,
                                    double maximum_colatitude_rad,
                                    double ionosphere_radius_re,
                                    double trace_radius_re) {
  if (mi_north_snapshot_allocate(theta_points, azimuth_points) != 0 ||
      mi_north_snapshot_configure(maximum_colatitude_rad,
                                  ionosphere_radius_re,
                                  trace_radius_re) != 0) {
    mi_north_snapshot_clear();
    return -1;
  }
  return 0;
}

int gamera_mi_snapshot_test_publish(double epoch_s,
                                    const double *north_potential_v) {
  return mi_north_snapshot_publish(epoch_s, north_potential_v);
}

void gamera_mi_snapshot_test_finalize(void) {
  mi_north_snapshot_clear();
}
#endif


#ifdef GAMERA_MI_COUPLING
enum {
  MI_NORTH_INDEX = 0,
  MI_SOUTH_INDEX = 1,
  MI_HEMISPHERE_COUNT = 2,
  MI_PATCH_COUNT = 2,
  MI_SOURCE_FAC = 0,
  MI_SOURCE_DENSITY = 1,
  MI_SOURCE_SOUND_SPEED = 2,
  MI_SOURCE_COUNT = 3
};

/*
 * FAC is a derivative of the magnetic field, so even a divergence-preserving
 * overset reconciliation leaves different truncation errors on the rotated
 * Yin and Yang stencils.  A hard owner switch turns that benign candidate
 * difference into a line source for the ionospheric Poisson solve.  Compose
 * the two valid candidates with the same four-cell margin scale used by the
 * overset interface, while suppressing candidates whose bilinear/curl
 * stencil lies in the first fringe cell.
 */
static const double FAC_BLEND_WIDTH_CELLS = 4.0;
static const double FAC_TAPER_ZERO_MARGIN_CELLS = 0.5;
static const double FAC_TAPER_FULL_MARGIN_CELLS = 1.5;

typedef struct {
  const gamera_no_grid *grid;
  const gamera_no_background_data *background;
  double *local_source;
  double *gathered_source;
  double *fac[MI_HEMISPHERE_COUNT];
  double *potential[MI_HEMISPHERE_COUNT];
  double *electron_number_flux[MI_HEMISPHERE_COUNT];
  double *electron_energy_flux[MI_HEMISPHERE_COUNT];
  double *diffuse_precipitation_selector[MI_HEMISPHERE_COUNT];
  double *euv_pedersen_conductance[MI_HEMISPHERE_COUNT];
  double *euv_hall_conductance[MI_HEMISPHERE_COUNT];
  double *auroral_pedersen_conductance[MI_HEMISPHERE_COUNT];
  double *auroral_hall_conductance[MI_HEMISPHERE_COUNT];
  double *pedersen_conductance[MI_HEMISPHERE_COUNT];
  double *hall_conductance[MI_HEMISPHERE_COUNT];
  double *dpb_mask[MI_HEMISPHERE_COUNT];
  double *dpb_boundary_latitude_deg[MI_HEMISPHERE_COUNT];
  gamera_mi_hybrid_dpb_state dpb_state[MI_HEMISPHERE_COUNT];
  gamera_mi_hybrid_dpb_stats dpb_stats[MI_HEMISPHERE_COUNT];
  int conductance_ready[MI_HEMISPHERE_COUNT];
  /* Rank-zero-only audit arrays for the independently interpolated patch
   * candidates and the production composition weight. */
  double *patch_fac[MI_HEMISPHERE_COUNT][MI_PATCH_COUNT];
  double *patch_margin[MI_HEMISPHERE_COUNT][MI_PATCH_COUNT];
  double *patch_fac_delta[MI_HEMISPHERE_COUNT];
  double *patch_yang_weight[MI_HEMISPHERE_COUNT];
  int *patch_valid[MI_HEMISPHERE_COUNT][MI_PATCH_COUNT];
  int *patch_owner[MI_HEMISPHERE_COUNT];
  double *saved_emf[GAMERA_NO_DIM];
  double sample_radius_re;
  double maximum_colatitude;
  double sample_maximum_colatitude;
  double next_update_s;
  double next_diagnostic_s;
  double maximum_cached_emf;
  double maximum_prehold_difference;
  double maximum_posthold_difference;
  long long held_apply_count;
  int update_count;
  int prepared;
  int solution_ready;
} mi_state_t;

static mi_state_t mi_state;

static int write_scalar_double_attribute(hid_t object, const char *name,
                                         double value) {
  hid_t space = H5Screate(H5S_SCALAR);
  if (space < 0) {
    return -1;
  }
  hid_t attribute = H5Acreate2(object, name, H5T_NATIVE_DOUBLE, space,
                               H5P_DEFAULT, H5P_DEFAULT);
  const int status =
      attribute < 0 ||
              H5Awrite(attribute, H5T_NATIVE_DOUBLE, &value) < 0
          ? -1
          : 0;
  if (attribute >= 0) {
    H5Aclose(attribute);
  }
  H5Sclose(space);
  return status;
}

static int write_scalar_int_attribute(hid_t object, const char *name,
                                      int value) {
  hid_t space = H5Screate(H5S_SCALAR);
  if (space < 0) {
    return -1;
  }
  hid_t attribute = H5Acreate2(object, name, H5T_NATIVE_INT, space,
                               H5P_DEFAULT, H5P_DEFAULT);
  const int status =
      attribute < 0 || H5Awrite(attribute, H5T_NATIVE_INT, &value) < 0
          ? -1
          : 0;
  if (attribute >= 0) {
    H5Aclose(attribute);
  }
  H5Sclose(space);
  return status;
}

static int write_double_dataset(hid_t parent, const char *name, int rank_count,
                                const hsize_t dimensions[],
                                const double *values) {
  hid_t space = H5Screate_simple(rank_count, dimensions, NULL);
  if (space < 0) {
    return -1;
  }
  hid_t dataset = H5Dcreate2(parent, name, H5T_NATIVE_DOUBLE, space,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  const int status =
      dataset < 0 ||
              H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                       H5P_DEFAULT, values) < 0
          ? -1
          : 0;
  if (dataset >= 0) {
    H5Dclose(dataset);
  }
  H5Sclose(space);
  return status;
}

static int write_int_dataset(hid_t parent, const char *name, int rank_count,
                             const hsize_t dimensions[], const int *values) {
  hid_t space = H5Screate_simple(rank_count, dimensions, NULL);
  if (space < 0) {
    return -1;
  }
  hid_t dataset = H5Dcreate2(parent, name, H5T_NATIVE_INT, space,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  const int status =
      dataset < 0 ||
              H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                       H5P_DEFAULT, values) < 0
          ? -1
          : 0;
  if (dataset >= 0) {
    H5Dclose(dataset);
  }
  H5Sclose(space);
  return status;
}

static int write_mi_snapshot(double time_code, double time_s) {
  if (!mi_config.diagnostics_enabled || rank != 0) {
    return 0;
  }
  const int write_patch_audit = !analysis_output_enabled;
  char path[96];
  const int path_length = snprintf(path, sizeof(path),
                                   "mi_ionosphere_%06d.h5",
                                   mi_state.update_count);
  if (path_length < 0 || (size_t)path_length >= sizeof(path)) {
    return -1;
  }
  hid_t file = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file < 0) {
    return -1;
  }
  const size_t np = (size_t)mi_config.longitude_count;
  const size_t nt = (size_t)mi_config.colatitude_count;
  int status =
      !mi_north_snapshot.available ||
              mi_north_snapshot.theta_points != nt ||
              mi_north_snapshot.azimuth_points != np
          ? -1
          : 0;
  if (status == 0) {
    const hsize_t longitude_dimension[1] = {(hsize_t)np};
    const hsize_t colatitude_dimension[1] = {(hsize_t)nt};
    status |= write_double_dataset(file, "longitude_rad", 1,
                                   longitude_dimension,
                                   mi_north_snapshot.longitude_rad);
    status |= write_double_dataset(file, "colatitude_rad", 1,
                                   colatitude_dimension,
                                   mi_north_snapshot.colatitude_rad);
  }
  status |= write_scalar_int_attribute(file, "schema_version", 3);
  status |= write_scalar_int_attribute(file, "problem_id", GAMERA_PROBLEM_ID);
  status |= write_scalar_int_attribute(file, "update", mi_state.update_count);
  status |= write_scalar_int_attribute(file, "polar_values_are_nodes", 1);
  /* Both hemispheres are indexed by the same native Cartesian azimuth
   * phi=atan2(y,x).  South must not be reflected by diagnostics. */
  status |= write_scalar_int_attribute(
      file, "hemisphere_longitude_axes_identical", 1);
  status |= write_scalar_int_attribute(file, "fac_uses_residual_field_only", 1);
  /* Compact analysis mode omits the optional patch-composition audit
   * matrices but retains all science fields, including precipitation and
   * conductance. */
  status |= write_scalar_int_attribute(file,
                                       "fac_patch_audit_schema_version",
                                       write_patch_audit ? 2 : 0);
  status |= write_scalar_int_attribute(file, "compact_science_output",
                                       write_patch_audit ? 0 : 1);
  status |= write_scalar_double_attribute(file, "fac_blend_width_cells",
                                          FAC_BLEND_WIDTH_CELLS);
  status |= write_scalar_double_attribute(
      file, "fac_taper_zero_margin_cells", FAC_TAPER_ZERO_MARGIN_CELLS);
  status |= write_scalar_double_attribute(
      file, "fac_taper_full_margin_cells", FAC_TAPER_FULL_MARGIN_CELLS);
#ifdef GAMERA_YINYANG_HDIV_RECONCILE
  status |= write_scalar_int_attribute(file, "fac_composition_mode", 1);
  status |= write_scalar_int_attribute(file,
                                       "owner_patch_is_reference_only", 1);
#else
  status |= write_scalar_int_attribute(file, "fac_composition_mode", 0);
  status |= write_scalar_int_attribute(file,
                                       "owner_patch_is_reference_only", 0);
#endif
  status |= write_scalar_int_attribute(file,
                                       "fac_composition_hard_owner_value", 0);
  status |= write_scalar_int_attribute(file,
                                       "fac_composition_safe_pou_value", 1);
  status |= write_scalar_int_attribute(file, "current_shell_zero_based",
                                       mi_config.current_shell);
  status |= write_scalar_double_attribute(file, "time_code", time_code);
  status |= write_scalar_double_attribute(file, "time_seconds", time_s);
  status |= write_scalar_double_attribute(file, "sample_radius_RE",
                                          mi_state.sample_radius_re);
  status |= write_scalar_double_attribute(file, "ionosphere_radius_RE",
                                          mi_config.ionosphere_radius_re);
  status |= write_scalar_double_attribute(
      file, "target_low_latitude_degrees",
      90.0 - mi_state.maximum_colatitude * 180.0 / PI);
  status |= write_scalar_double_attribute(
      file, "sampled_shell_low_latitude_degrees",
      90.0 - mi_state.sample_maximum_colatitude * 180.0 / PI);
  status |= write_scalar_double_attribute(file, "pedersen_siemens",
                                          mi_config.pedersen_siemens);
  status |= write_scalar_double_attribute(file, "hall_siemens",
                                          mi_config.hall_siemens);
  status |= write_scalar_int_attribute(
      file, "electron_precipitation_enabled",
      mi_config.electron_precipitation_enabled);
  status |= write_scalar_int_attribute(
      file, "electron_precipitation_hall_model",
      mi_config.electron_precipitation_hall_model);
  status |= write_scalar_double_attribute(
      file, "electron_precipitation_beta",
      mi_config.electron_precipitation_beta);
  status |= write_scalar_double_attribute(
      file, "electron_precipitation_ramp_seconds",
      mi_config.electron_precipitation_ramp_s);
  status |= write_scalar_int_attribute(
      file, "conductance_feedback_enabled",
      mi_config.electron_precipitation_conductance_feedback_enabled);
  status |= write_scalar_double_attribute(
      file, "conductance_smoothing_seconds",
      mi_config.electron_precipitation_conductance_smoothing_s);
  status |= write_scalar_double_attribute(file, "euv_f107_sfu",
                                          mi_config.euv_f107);
  status |= write_scalar_double_attribute(file, "sunward_direction_x",
                                          mi_config.sunward_direction[0]);
  status |= write_scalar_double_attribute(file, "sunward_direction_y",
                                          mi_config.sunward_direction[1]);
  status |= write_scalar_double_attribute(file, "sunward_direction_z",
                                          mi_config.sunward_direction[2]);
  status |= write_scalar_int_attribute(file, "hybrid_dpb_enabled",
                                       mi_config.hybrid_dpb_enabled);
  status |= write_scalar_double_attribute(
      file, "dpb_nightside_poleward_limit_degrees",
      mi_config.dpb_nightside_poleward_limit_deg);
  status |= write_scalar_double_attribute(
      file, "dpb_oval_center_latitude_degrees",
      mi_config.dpb_oval_center_latitude_deg);
  status |= write_scalar_double_attribute(file, "dpb_oval_center_mlt_hours",
                                          mi_config.dpb_oval_center_mlt_h);
  status |= write_scalar_double_attribute(
      file, "minimum_pedersen_siemens",
      mi_config.minimum_pedersen_siemens);
  status |= write_scalar_double_attribute(file, "minimum_hall_siemens",
                                          mi_config.minimum_hall_siemens);
  status |= write_scalar_double_attribute(
      file, "maximum_hall_to_pedersen",
      mi_config.maximum_hall_to_pedersen);
  status |= write_scalar_double_attribute(file, "coupling_interval_seconds",
                                          mi_config.coupling_interval_s);
  status |= write_scalar_double_attribute(
      file, "diagnostics_interval_seconds",
      mi_config.diagnostics_interval_s);
  status |= write_scalar_double_attribute(file,
                                          "low_latitude_potential_V", 0.0);
  const hsize_t field_dimensions[2] = {(hsize_t)nt, (hsize_t)np};
  const char *group_name[2] = {"north", "south"};
  for (int hemisphere = 0; hemisphere < MI_HEMISPHERE_COUNT; ++hemisphere) {
    hid_t group = H5Gcreate2(file, group_name[hemisphere], H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
    if (group < 0) {
      status = -1;
      continue;
    }
    status |= write_double_dataset(group, "fac_parallel_A_m2", 2,
                                   field_dimensions,
                                   mi_state.fac[hemisphere]);
    status |= write_double_dataset(group, "potential_V", 2,
                                   field_dimensions,
                                   mi_state.potential[hemisphere]);
    status |= write_double_dataset(
        group, "F_N_cm2_s", 2, field_dimensions,
        mi_state.electron_number_flux[hemisphere]);
    status |= write_double_dataset(
        group, "F_E_erg_cm2_s", 2, field_dimensions,
        mi_state.electron_energy_flux[hemisphere]);
    status |= write_double_dataset(
        group, "diffuse_precipitation_selector", 2, field_dimensions,
        mi_state.diffuse_precipitation_selector[hemisphere]);
    status |= write_double_dataset(
        group, "Sigma_P_S", 2, field_dimensions,
        mi_state.pedersen_conductance[hemisphere]);
    status |= write_double_dataset(
        group, "Sigma_H_S", 2, field_dimensions,
        mi_state.hall_conductance[hemisphere]);
    status |= write_double_dataset(
        group, "Sigma_P_EUV_S", 2, field_dimensions,
        mi_state.euv_pedersen_conductance[hemisphere]);
    status |= write_double_dataset(
        group, "Sigma_H_EUV_S", 2, field_dimensions,
        mi_state.euv_hall_conductance[hemisphere]);
    status |= write_double_dataset(
        group, "DPB_mask", 2, field_dimensions,
        mi_state.dpb_mask[hemisphere]);
    const hsize_t dpb_longitude_dimension[1] = {(hsize_t)np};
    status |= write_double_dataset(
        group, "DPB_boundary_MLAT_deg", 1, dpb_longitude_dimension,
        mi_state.dpb_boundary_latitude_deg[hemisphere]);
    status |= write_scalar_double_attribute(
        group, "dpb_offset_oval_radius_degrees",
        mi_state.dpb_stats[hemisphere].filtered_radius_deg);
    status |= write_scalar_double_attribute(
        group, "dpb_nightside_boundary_degrees",
        mi_state.dpb_stats[hemisphere].nightside_boundary_deg);
    status |= write_scalar_double_attribute(
        group, "dpb_dayside_boundary_degrees",
        mi_state.dpb_stats[hemisphere].dayside_boundary_deg);
    if (write_patch_audit) {
      status |= write_double_dataset(
          group, "fac_yin_A_m2", 2, field_dimensions,
          mi_state.patch_fac[hemisphere][GAMERA_NO_YIN_PATCH]);
      status |= write_double_dataset(
          group, "fac_yang_A_m2", 2, field_dimensions,
          mi_state.patch_fac[hemisphere][GAMERA_NO_YANG_PATCH]);
      status |= write_double_dataset(
          group, "fac_yang_minus_yin_A_m2", 2, field_dimensions,
          mi_state.patch_fac_delta[hemisphere]);
      status |= write_double_dataset(
          group, "fac_weight_yang", 2, field_dimensions,
          mi_state.patch_yang_weight[hemisphere]);
      status |= write_double_dataset(
          group, "margin_yin_cells", 2, field_dimensions,
          mi_state.patch_margin[hemisphere][GAMERA_NO_YIN_PATCH]);
      status |= write_double_dataset(
          group, "margin_yang_cells", 2, field_dimensions,
          mi_state.patch_margin[hemisphere][GAMERA_NO_YANG_PATCH]);
      status |= write_int_dataset(
          group, "valid_yin", 2, field_dimensions,
          mi_state.patch_valid[hemisphere][GAMERA_NO_YIN_PATCH]);
      status |= write_int_dataset(
          group, "valid_yang", 2, field_dimensions,
          mi_state.patch_valid[hemisphere][GAMERA_NO_YANG_PATCH]);
      status |= write_int_dataset(group, "owner_patch", 2,
                                  field_dimensions,
                                  mi_state.patch_owner[hemisphere]);
      status |= write_scalar_int_attribute(group, "owner_unmapped_value",
                                           -1);
      status |= write_scalar_int_attribute(group, "owner_yin_value",
                                           GAMERA_NO_YIN_PATCH);
      status |= write_scalar_int_attribute(group, "owner_yang_value",
                                           GAMERA_NO_YANG_PATCH);
      status |= write_scalar_int_attribute(
          group, "invalid_candidate_values_are_nan", 1);
    }
    status |= write_scalar_int_attribute(
        group, "upward_fac_multiplier",
        hemisphere == MI_NORTH_INDEX ? -1 : 1);
    H5Gclose(group);
  }
  if (H5Fclose(file) < 0) {
    status = -1;
  }
  return status == 0 ? 0 : -1;
}

static double vector_dot(gamera_no_vec3 left, gamera_no_vec3 right) {
  double result = 0.0;
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    result += left.value[component] * right.value[component];
  }
  return result;
}

static gamera_no_vec3 vector_cross(gamera_no_vec3 left,
                                   gamera_no_vec3 right) {
  return (gamera_no_vec3){{
      left.value[1] * right.value[2] - left.value[2] * right.value[1],
      left.value[2] * right.value[0] - left.value[0] * right.value[2],
      left.value[0] * right.value[1] - left.value[1] * right.value[0]}};
}

static double vector_norm(gamera_no_vec3 value) {
  return sqrt(fmax(0.0, vector_dot(value, value)));
}

static int allocate_doubles(size_t count, double **result) {
  if (result == NULL || count == 0U || count > SIZE_MAX / sizeof(double)) {
    return -1;
  }
  *result = (double *)calloc(count, sizeof(double));
  return *result == NULL ? -1 : 0;
}

static int allocate_ints(size_t count, int **result) {
  if (result == NULL || count == 0U || count > SIZE_MAX / sizeof(int)) {
    return -1;
  }
  *result = (int *)calloc(count, sizeof(int));
  return *result == NULL ? -1 : 0;
}

static size_t local_angular_index(int j, int k) {
  return (size_t)j * (size_t)config.nk + (size_t)k;
}

static int global_patch_value(int source, int patch, int global_j,
                              int global_k, double *value) {
  if (value == NULL || source < 0 || source >= MI_SOURCE_COUNT || patch < 0 ||
      patch >= patch_count || global_j < 0 ||
      global_j >= config.nj_global || global_k < 0 ||
      global_k >= config.nk_global) {
    return -1;
  }
  const int radial_coord = mi_config.current_shell / config.ni;
  int coords[3] = {radial_coord, global_j / config.nj,
                   global_k / config.nk};
  int patch_cart_rank;
  if (MPI_Cart_rank(comm_cart, coords, &patch_cart_rank) != MPI_SUCCESS) {
    return -1;
  }
  const int world_rank = patch * patch_size + patch_cart_rank;
  const size_t local_count = (size_t)config.nj * (size_t)config.nk;
  const size_t offset =
      ((size_t)world_rank * MI_SOURCE_COUNT + (size_t)source) * local_count +
      local_angular_index(global_j % config.nj, global_k % config.nk);
  *value = mi_state.gathered_source[offset];
  return isfinite(*value) ? 0 : -1;
}

static void interpolation_coordinate(double logical, int count, int *base,
                                     double *fraction) {
  if (logical <= 0.0) {
    *base = 0;
    *fraction = 0.0;
  } else if (logical >= (double)(count - 1)) {
    *base = count - 2;
    *fraction = 1.0;
  } else {
    *base = (int)floor(logical);
    *fraction = logical - (double)*base;
  }
}

static int angular_margin(int patch, gamera_no_vec3 point, int *valid,
                          double *margin) {
  double radius, theta, phi;
  if (valid == NULL || margin == NULL ||
      gamera_no_yinyang_global_to_logical(patch, point, &radius, &theta,
                                          &phi) != 0) {
    return -1;
  }
  const double theta_spacing =
      (x2max_global - x2min_global) / (double)config.nj_global;
  const double phi_spacing =
      (x3max_global - x3min_global) / (double)config.nk_global;
  const double q_theta = (theta - x2min_global) / theta_spacing - 0.5;
  const double q_phi = (phi - x3min_global) / phi_spacing - 0.5;
  const double tolerance = 2.0e-2;
  *valid = q_theta >= -0.5 - tolerance &&
           q_theta <= (double)config.nj_global - 0.5 + tolerance &&
           q_phi >= -0.5 - tolerance &&
           q_phi <= (double)config.nk_global - 0.5 + tolerance;
  *margin = fmin(fmin(q_theta + 0.5,
                      (double)config.nj_global - 0.5 - q_theta),
                 fmin(q_phi + 0.5,
                      (double)config.nk_global - 0.5 - q_phi));
  return 0;
}

static int owner_patch(gamera_no_vec3 point, int *owner) {
  int valid[2];
  double margin[2];
  if (owner == NULL || angular_margin(0, point, &valid[0], &margin[0]) != 0 ||
      angular_margin(1, point, &valid[1], &margin[1]) != 0 ||
      (!valid[0] && !valid[1])) {
    return -1;
  }
  if (!valid[1]) {
    *owner = 0;
  } else if (!valid[0]) {
    *owner = 1;
  } else {
#ifdef GAMERA_YINYANG_MFE_INTERFACE
    /*
     * The reference Yin-Yang MFE composite keeps a fixed patch priority and
     * changes patches only at a donor-conditioned fringe.  In the
     * ghost-only interface mode, doing the same here avoids taking a hard
     * equal-margin switch between two independently CT-evolved active
     * histories.  Yin's angular boundary cells already use Yang-filled
     * primitive and magnetic ghosts in the current stencil.  Keep writing
     * both candidates and margins to the M-I audit file so this policy is
     * measured rather than hidden.
     */
    *owner = GAMERA_NO_YIN_PATCH;
#else
    *owner = margin[1] > margin[0] + 1.0e-12 ? 1 : 0;
#endif
  }
  return 0;
}

static int interpolate_patch_source(int source, int patch,
                                    gamera_no_vec3 point, double *value) {
  double radius, theta, phi;
  if (value == NULL ||
      gamera_no_yinyang_global_to_logical(patch, point, &radius, &theta,
                                          &phi) != 0) {
    return -1;
  }
  const double q_theta =
      (theta - x2min_global) * (double)config.nj_global /
          (x2max_global - x2min_global) -
      0.5;
  const double q_phi =
      (phi - x3min_global) * (double)config.nk_global /
          (x3max_global - x3min_global) -
      0.5;
  int base_theta, base_phi;
  double fraction_theta, fraction_phi;
  interpolation_coordinate(q_theta, config.nj_global, &base_theta,
                           &fraction_theta);
  interpolation_coordinate(q_phi, config.nk_global, &base_phi,
                           &fraction_phi);
  *value = 0.0;
  for (int dj = 0; dj <= 1; ++dj) {
    for (int dk = 0; dk <= 1; ++dk) {
      double corner;
      if (global_patch_value(source, patch, base_theta + dj, base_phi + dk,
                             &corner) != 0) {
        return -1;
      }
      const double wj = dj == 0 ? 1.0 - fraction_theta : fraction_theta;
      const double wk = dk == 0 ? 1.0 - fraction_phi : fraction_phi;
      *value += wj * wk * corner;
    }
  }
  return isfinite(*value) ? 0 : -1;
}

#ifdef GAMERA_YINYANG_HDIV_RECONCILE
static double fac_stencil_taper(double margin_cells) {
  if (margin_cells <= FAC_TAPER_ZERO_MARGIN_CELLS) {
    return 0.0;
  }
  if (margin_cells >= FAC_TAPER_FULL_MARGIN_CELLS) {
    return 1.0;
  }
  const double coordinate =
      (margin_cells - FAC_TAPER_ZERO_MARGIN_CELLS) /
      (FAC_TAPER_FULL_MARGIN_CELLS - FAC_TAPER_ZERO_MARGIN_CELLS);
  return coordinate * coordinate * (3.0 - 2.0 * coordinate);
}
#endif

static int compose_patch_fac(gamera_no_vec3 point, double *selected,
                             int *owner, double candidate[MI_PATCH_COUNT],
                             int valid[MI_PATCH_COUNT],
                             double margin[MI_PATCH_COUNT],
                             double *weight_yang) {
  if (selected == NULL || owner == NULL || candidate == NULL ||
      valid == NULL || margin == NULL || weight_yang == NULL ||
      owner_patch(point, owner) != 0) {
    return -1;
  }
  for (int patch = 0; patch < MI_PATCH_COUNT; ++patch) {
    candidate[patch] = NAN;
    valid[patch] = 0;
    margin[patch] = NAN;
    int geometrically_valid = 0;
    if (angular_margin(patch, point, &geometrically_valid, &margin[patch]) !=
        0) {
      return -1;
    }
    if (geometrically_valid) {
      /* A geometrically present donor with a non-finite/missing FAC is data
       * corruption, not an absent patch.  Fail collectively instead of
       * silently hiding a one-patch instability behind the other candidate. */
      if (interpolate_patch_source(MI_SOURCE_FAC, patch, point,
                                   &candidate[patch]) != 0) {
        return -1;
      }
      valid[patch] = 1;
    }
  }
  if (!valid[GAMERA_NO_YIN_PATCH] && !valid[GAMERA_NO_YANG_PATCH]) {
    return -1;
  }
  if (!valid[GAMERA_NO_YANG_PATCH]) {
    *weight_yang = 0.0;
    *selected = candidate[GAMERA_NO_YIN_PATCH];
  } else if (!valid[GAMERA_NO_YIN_PATCH]) {
    *weight_yang = 1.0;
    *selected = candidate[GAMERA_NO_YANG_PATCH];
  } else {
#ifdef GAMERA_YINYANG_HDIV_RECONCILE
    const double geometric_weight =
        0.5 * (1.0 + tanh((margin[GAMERA_NO_YANG_PATCH] -
                           margin[GAMERA_NO_YIN_PATCH]) /
                          FAC_BLEND_WIDTH_CELLS));
    const double yin_weight =
        (1.0 - geometric_weight) *
        fac_stencil_taper(margin[GAMERA_NO_YIN_PATCH]);
    const double yang_weight =
        geometric_weight *
        fac_stencil_taper(margin[GAMERA_NO_YANG_PATCH]);
    const double weight_sum = yin_weight + yang_weight;
    if (weight_sum > 0.0 && isfinite(weight_sum)) {
      *weight_yang = yang_weight / weight_sum;
    } else {
      /* Reintroducing a hard switch here would silently restore the FAC seam
       * on an overlap too narrow for either curl/bilinear stencil. */
      log_error("No stencil-safe Yin/Yang FAC candidate at xyz=(%.9g,%.9g,"
                "%.9g), margins=(%.9g,%.9g) cells",
                point.value[0], point.value[1], point.value[2],
                margin[GAMERA_NO_YIN_PATCH],
                margin[GAMERA_NO_YANG_PATCH]);
      return -1;
    }
#else
    *weight_yang = *owner == GAMERA_NO_YANG_PATCH ? 1.0 : 0.0;
#endif
    *selected = (1.0 - *weight_yang) * candidate[GAMERA_NO_YIN_PATCH] +
                *weight_yang * candidate[GAMERA_NO_YANG_PATCH];
  }
  if (!isfinite(*weight_yang) || *weight_yang < 0.0 ||
      *weight_yang > 1.0) {
    return -1;
  }
  return isfinite(*selected) ? 0 : -1;
}

/* Reuse the exact valid-patch set and production composition weight selected
 * for FAC so density and sound speed cannot acquire a different Yin/Yang seam
 * before entering the precipitation model. */
static int compose_patch_source(gamera_no_vec3 point, int source,
                                const int valid[MI_PATCH_COUNT],
                                double weight_yang, double *selected) {
  if (source <= MI_SOURCE_FAC || source >= MI_SOURCE_COUNT || valid == NULL ||
      selected == NULL || !isfinite(weight_yang) || weight_yang < 0.0 ||
      weight_yang > 1.0) {
    return -1;
  }
  double candidate[MI_PATCH_COUNT] = {NAN, NAN};
  for (int patch = 0; patch < MI_PATCH_COUNT; ++patch) {
    if (valid[patch] &&
        interpolate_patch_source(source, patch, point, &candidate[patch]) !=
            0) {
      return -1;
    }
  }
  if (!valid[GAMERA_NO_YANG_PATCH]) {
    *selected = candidate[GAMERA_NO_YIN_PATCH];
  } else if (!valid[GAMERA_NO_YIN_PATCH]) {
    *selected = candidate[GAMERA_NO_YANG_PATCH];
  } else {
    *selected =
        (1.0 - weight_yang) * candidate[GAMERA_NO_YIN_PATCH] +
        weight_yang * candidate[GAMERA_NO_YANG_PATCH];
  }
  return isfinite(*selected) ? 0 : -1;
}

static int reset_fac_patch_audit(size_t polar_count) {
  if (!mi_config.diagnostics_enabled || analysis_output_enabled || rank != 0) {
    return 0;
  }
  for (int hemisphere = 0; hemisphere < MI_HEMISPHERE_COUNT; ++hemisphere) {
    if (mi_state.patch_owner[hemisphere] == NULL ||
        mi_state.patch_fac_delta[hemisphere] == NULL ||
        mi_state.patch_yang_weight[hemisphere] == NULL) {
      return -1;
    }
    for (size_t index = 0; index < polar_count; ++index) {
      mi_state.patch_owner[hemisphere][index] = -1;
      mi_state.patch_fac_delta[hemisphere][index] = NAN;
      mi_state.patch_yang_weight[hemisphere][index] = NAN;
    }
    for (int patch = 0; patch < MI_PATCH_COUNT; ++patch) {
      if (mi_state.patch_fac[hemisphere][patch] == NULL ||
          mi_state.patch_margin[hemisphere][patch] == NULL ||
          mi_state.patch_valid[hemisphere][patch] == NULL) {
        return -1;
      }
      for (size_t index = 0; index < polar_count; ++index) {
        mi_state.patch_fac[hemisphere][patch][index] = NAN;
        mi_state.patch_margin[hemisphere][patch][index] = NAN;
        mi_state.patch_valid[hemisphere][patch][index] = 0;
      }
    }
  }
  return 0;
}

/* Record the exact candidates and weight used by the production composition.
 * A missing candidate is represented by valid_patch=0 and NaN. */
static void capture_fac_patch_audit(int hemisphere, size_t index,
                                    int selected_owner,
                                    const double candidate[MI_PATCH_COUNT],
                                    const int valid[MI_PATCH_COUNT],
                                    const double margin[MI_PATCH_COUNT],
                                    double weight_yang) {
  if (!mi_config.diagnostics_enabled || analysis_output_enabled || rank != 0 ||
      hemisphere < 0 || hemisphere >= MI_HEMISPHERE_COUNT) {
    return;
  }
  mi_state.patch_owner[hemisphere][index] = selected_owner;
  mi_state.patch_yang_weight[hemisphere][index] = weight_yang;
  for (int patch = 0; patch < MI_PATCH_COUNT; ++patch) {
    mi_state.patch_margin[hemisphere][patch][index] = margin[patch];
    if (valid[patch]) {
      mi_state.patch_fac[hemisphere][patch][index] = candidate[patch];
      mi_state.patch_valid[hemisphere][patch][index] = valid[patch];
    }
  }
  if (mi_state.patch_valid[hemisphere][GAMERA_NO_YIN_PATCH][index] != 0 &&
      mi_state.patch_valid[hemisphere][GAMERA_NO_YANG_PATCH][index] != 0) {
    mi_state.patch_fac_delta[hemisphere][index] =
        mi_state.patch_fac[hemisphere][GAMERA_NO_YANG_PATCH][index] -
        mi_state.patch_fac[hemisphere][GAMERA_NO_YIN_PATCH][index];
  }
}

static int fill_local_source(void) {
  const size_t local_count = (size_t)config.nj * (size_t)config.nk;
  memset(mi_state.local_source, 0,
         MI_SOURCE_COUNT * local_count * sizeof(*mi_state.local_source));
  const int radial_coord = mi_config.current_shell / config.ni;
  if (proc_coords[0] != radial_coord) {
    return 0;
  }
  const size_t i = (size_t)(is + mi_config.current_shell % config.ni);
  const gamera_no_storage *storage = gamera_no_legacy_storage();
  if (storage == NULL) {
    return -1;
  }
  const double current_norm =
      norm_config.B_Norm /
      (norm_config.mu0 > 0.0 ? norm_config.mu0 : 4.0e-7 * PI) /
      norm_config.x_Norm;
  int failed = 0;
#pragma omp parallel for collapse(2) reduction(| : failed) schedule(static)
  for (int local_j = 0; local_j < config.nj; ++local_j) {
    for (int local_k = 0; local_k < config.nk; ++local_k) {
      const size_t j = (size_t)(js + local_j);
      const size_t k = (size_t)(ks + local_k);
      const size_t cell =
          gamera_no_index3(mi_state.grid->cell_extent, i, j, k);
      gamera_no_primitive primitive;
      gamera_no_vec3 current;
      if (gamera_no_conserved_to_primitive(
              &storage->conserved[cell * GAMERA_NO_FLUX_COUNT], gamma_val,
              rho_floor, p_floor, &primitive) != 0 ||
          !(primitive.density > 0.0) || !(primitive.pressure > 0.0) ||
          gamera_no_cell_current_from_residual(
              mi_state.grid, storage->cell_magnetic, i, j, k, &current) != 0) {
        failed = 1;
        continue;
      }
      const gamera_no_vec3 b0 = mi_state.background->cell_magnetic[cell];
      const double b0_magnitude = vector_norm(b0);
      const gamera_no_vec3 point = mi_state.grid->cell[cell].centroid;
      const double radius = vector_norm(point);
      if (!(b0_magnitude > 0.0) || !(radius > 0.0)) {
        failed = 1;
        continue;
      }
      const double polar_colatitude =
          acos(fmax(-1.0, fmin(1.0, fabs(point.value[2]) / radius)));
      const double field_ratio = gamera_mi_dipole_field_ratio(
          polar_colatitude, radius, mi_config.ionosphere_radius_re);
      const double parallel = vector_dot(current, b0) / b0_magnitude;
      const size_t angular = local_angular_index(local_j, local_k);
      mi_state.local_source[MI_SOURCE_FAC * local_count + angular] =
          parallel * field_ratio * current_norm;
      mi_state.local_source[MI_SOURCE_DENSITY * local_count + angular] =
          primitive.density * norm_config.rho_Norm;
      mi_state.local_source[MI_SOURCE_SOUND_SPEED * local_count + angular] =
          sqrt(gamma_val * primitive.pressure / primitive.density) *
          norm_config.u_Norm;
      if (!isfinite(
              mi_state.local_source[MI_SOURCE_FAC * local_count + angular]) ||
          !(mi_state.local_source[MI_SOURCE_DENSITY * local_count + angular] >
            0.0) ||
          !(mi_state.local_source[MI_SOURCE_SOUND_SPEED * local_count +
                                  angular] >= 0.0)) {
        failed = 1;
      }
    }
  }
  return failed ? -1 : 0;
}

static int gather_source(void) {
  const size_t local_count = (size_t)config.nj * (size_t)config.nk;
  if (local_count > (size_t)INT32_MAX / MI_SOURCE_COUNT ||
      fill_local_source() != 0) {
    return -1;
  }
  const int packed_count = (int)(MI_SOURCE_COUNT * local_count);
  return MPI_Allgather(mi_state.local_source, packed_count, MPI_DOUBLE,
                       mi_state.gathered_source, packed_count, MPI_DOUBLE,
                       MPI_COMM_WORLD) == MPI_SUCCESS
             ? 0
             : -1;
}

static int map_mhd_source_to_ionosphere(double time_s) {
  const size_t np = (size_t)mi_config.longitude_count;
  const size_t nt = (size_t)mi_config.colatitude_count;
  const size_t polar_count = np * nt;
  if (reset_fac_patch_audit(polar_count) != 0) {
    return -1;
  }
  gamera_mi_fedder95_config precipitation_config =
      gamera_mi_fedder95_default_config();
  precipitation_config.beta = mi_config.electron_precipitation_beta;
  precipitation_config.hall_model =
      (gamera_mi_hall_model)mi_config.electron_precipitation_hall_model;
  const double ramp_factor =
      mi_config.electron_precipitation_ramp_s > 0.0
          ? fmin(1.0, fmax(0.0, time_s /
                                    mi_config.electron_precipitation_ramp_s))
          : 1.0;
  gamera_mi_hybrid_dpb_config dpb_config =
      gamera_mi_hybrid_dpb_default_config();
  dpb_config.absolute_potential_v = mi_config.dpb_absolute_potential_v;
  dpb_config.adaptive_potential_fraction =
      mi_config.dpb_adaptive_potential_fraction;
  dpb_config.adaptive_potential_floor_v =
      mi_config.dpb_adaptive_potential_floor_v;
  dpb_config.potential_equatorward_offset_deg =
      mi_config.dpb_potential_equatorward_offset_deg;
  dpb_config.fac_equatorward_offset_deg =
      mi_config.dpb_fac_equatorward_offset_deg;
  dpb_config.fac_absolute_floor_a_m2 =
      mi_config.dpb_fac_absolute_floor_a_m2;
  dpb_config.minimum_boundary_latitude_deg =
      mi_config.dpb_minimum_boundary_latitude_deg;
  dpb_config.maximum_candidate_latitude_deg =
      mi_config.dpb_maximum_candidate_latitude_deg;
  dpb_config.transition_width_deg = mi_config.dpb_transition_width_deg;
  dpb_config.temporal_timescale_s = mi_config.dpb_temporal_timescale_s;
  dpb_config.maximum_slew_deg_per_s =
      mi_config.dpb_maximum_slew_deg_per_s;
  dpb_config.quiet_initial_nightside_latitude_deg =
      mi_config.dpb_quiet_initial_nightside_latitude_deg;
  dpb_config.oval_center_latitude_deg =
      mi_config.dpb_oval_center_latitude_deg;
  dpb_config.oval_center_mlt_h = mi_config.dpb_oval_center_mlt_h;
  dpb_config.nightside_poleward_limit_deg =
      mi_config.dpb_nightside_poleward_limit_deg;
  const double dtheta = mi_state.maximum_colatitude / (double)(nt - 1U);
  const double dphi = 2.0 * PI / (double)np;
  const double sunward_norm =
      sqrt(mi_config.sunward_direction[0] * mi_config.sunward_direction[0] +
           mi_config.sunward_direction[1] * mi_config.sunward_direction[1] +
           mi_config.sunward_direction[2] * mi_config.sunward_direction[2]);
  const double subsolar_colatitude =
      acos(fmax(-1.0, fmin(1.0,
                           mi_config.sunward_direction[2] / sunward_norm)));
  const double subsolar_longitude =
      atan2(mi_config.sunward_direction[1],
            mi_config.sunward_direction[0]);
  for (int hemisphere_index = 0; hemisphere_index < 2;
       ++hemisphere_index) {
    for (size_t i = 0U; i < nt; ++i) {
      const double colatitude = (double)i * dtheta;
      for (size_t j = 0U; j < np; ++j) {
        const size_t polar = i * np + j;
        mi_state.electron_number_flux[hemisphere_index][polar] = 0.0;
        mi_state.electron_energy_flux[hemisphere_index][polar] = 0.0;
        mi_state.diffuse_precipitation_selector[hemisphere_index][polar] = 0.0;
        mi_state.auroral_pedersen_conductance[hemisphere_index][polar] = 0.0;
        mi_state.auroral_hall_conductance[hemisphere_index][polar] = 0.0;
        mi_state.dpb_mask[hemisphere_index][polar] = 0.0;
        if (mi_config.electron_precipitation_enabled) {
          const double longitude = (double)j * dphi;
          /* Compute illumination from the Yin--Yang model's native Cartesian
           * sunward vector, not from another model's longitude convention. */
          const double solar_zenith = gamera_mi_solar_zenith_angle(
              colatitude, longitude, subsolar_colatitude,
              subsolar_longitude);
          double euv_pedersen;
          double euv_hall;
          if (gamera_mi_lompe_euv_conductance(
                  solar_zenith, mi_config.euv_f107,
                  &euv_pedersen, &euv_hall) != 0) {
            return -1;
          }
          mi_state.euv_pedersen_conductance[hemisphere_index][polar] =
              euv_pedersen;
          mi_state.euv_hall_conductance[hemisphere_index][polar] =
              euv_hall;
        } else {
          mi_state.euv_pedersen_conductance[hemisphere_index][polar] =
              mi_config.pedersen_siemens;
          mi_state.euv_hall_conductance[hemisphere_index][polar] =
              mi_config.hall_siemens;
          mi_state.pedersen_conductance[hemisphere_index][polar] =
              fmax(mi_config.minimum_pedersen_siemens,
                   mi_config.pedersen_siemens);
          mi_state.hall_conductance[hemisphere_index][polar] =
              fmax(mi_config.minimum_hall_siemens, mi_config.hall_siemens);
        }
      }
    }
  }
  for (int hemisphere_index = 0; hemisphere_index < 2;
       ++hemisphere_index) {
    const double z_sign = hemisphere_index == MI_NORTH_INDEX ? 1.0 : -1.0;
    for (size_t i = 0; i < nt; ++i) {
      const double ionosphere_colatitude = (double)i * dtheta;
      if (ionosphere_colatitude >
          mi_state.sample_maximum_colatitude + 4.0e-14) {
        memset(&mi_state.fac[hemisphere_index][i * np], 0,
               np * sizeof(double));
        continue;
      }
      double argument = sin(ionosphere_colatitude) *
                        sqrt(mi_state.sample_radius_re /
                             mi_config.ionosphere_radius_re);
      argument = fmax(0.0, fmin(1.0, argument));
      const double mhd_colatitude = asin(argument);
      for (size_t j = 0; j < np; ++j) {
        const double phi = (double)j * dphi;
        const gamera_no_vec3 point = {{
            mi_state.sample_radius_re * sin(mhd_colatitude) * cos(phi),
            mi_state.sample_radius_re * sin(mhd_colatitude) * sin(phi),
            z_sign * mi_state.sample_radius_re * cos(mhd_colatitude)}};
        int patch;
        double candidate[MI_PATCH_COUNT];
        int valid[MI_PATCH_COUNT];
        double margin[MI_PATCH_COUNT];
        double weight_yang;
        const size_t polar = i * np + j;
        if (compose_patch_fac(point, &mi_state.fac[hemisphere_index][polar],
                              &patch, candidate, valid, margin,
                              &weight_yang) != 0) {
          return -1;
        }
        if (mi_config.electron_precipitation_enabled) {
          double mass_density_kg_m3;
          double sound_speed_m_s;
          gamera_mi_electron_precipitation precipitation;
          const int hemisphere = hemisphere_index == MI_NORTH_INDEX ? -1 : 1;
          if (compose_patch_source(point, MI_SOURCE_DENSITY, valid,
                                   weight_yang, &mass_density_kg_m3) != 0 ||
              compose_patch_source(point, MI_SOURCE_SOUND_SPEED, valid,
                                   weight_yang, &sound_speed_m_s) != 0 ||
              gamera_mi_fedder95_electron_precipitation(
                  &precipitation_config, mass_density_kg_m3, sound_speed_m_s,
                  mi_state.fac[hemisphere_index][polar],
                  fmax(mi_config.minimum_pedersen_siemens,
                       mi_state.euv_pedersen_conductance[hemisphere_index]
                                                                 [polar]),
                  ramp_factor, hemisphere,
                  &precipitation) != 0) {
            return -1;
          }
          mi_state.electron_number_flux[hemisphere_index][polar] =
              precipitation.number_flux_cm2_s;
          mi_state.electron_energy_flux[hemisphere_index][polar] =
              precipitation.energy_flux_erg_cm2_s;
          mi_state.auroral_pedersen_conductance[hemisphere_index][polar] =
              precipitation.pedersen_siemens;
          mi_state.auroral_hall_conductance[hemisphere_index][polar] =
              precipitation.hall_siemens;
          mi_state.diffuse_precipitation_selector[hemisphere_index][polar] =
              precipitation.is_diffuse ? 1.0 : 0.0;
        }
        capture_fac_patch_audit(hemisphere_index, polar, patch, candidate,
                                valid, margin, weight_yang);
      }
    }
  }
  if (mi_config.electron_precipitation_enabled) {
    const double elapsed_s = mi_config.coupling_interval_s;
    const double smoothing_alpha =
        mi_config.electron_precipitation_conductance_smoothing_s > 0.0
            ? 1.0 - exp(-elapsed_s /
                        mi_config.electron_precipitation_conductance_smoothing_s)
            : 1.0;
    for (int hemisphere_index = 0; hemisphere_index < 2;
         ++hemisphere_index) {
      if (mi_config.hybrid_dpb_enabled) {
        const double upward_multiplier =
            hemisphere_index == MI_NORTH_INDEX ? -1.0 : 1.0;
        if (gamera_mi_hybrid_dpb_update(
                &dpb_config, nt, np, mi_state.maximum_colatitude,
                mi_state.sample_maximum_colatitude,
                mi_state.fac[hemisphere_index], upward_multiplier,
                mi_state.potential[hemisphere_index], elapsed_s,
                &mi_state.dpb_state[hemisphere_index],
                mi_state.dpb_boundary_latitude_deg[hemisphere_index],
                mi_state.dpb_mask[hemisphere_index],
                &mi_state.dpb_stats[hemisphere_index]) != 0) {
          return -1;
        }
      } else {
        for (size_t polar = 0U; polar < polar_count; ++polar) {
          mi_state.dpb_mask[hemisphere_index][polar] = 1.0;
        }
        for (size_t j = 0U; j < np; ++j) {
          mi_state.dpb_boundary_latitude_deg[hemisphere_index][j] =
              mi_config.dpb_minimum_boundary_latitude_deg;
        }
        mi_state.dpb_stats[hemisphere_index] =
            (gamera_mi_hybrid_dpb_stats){0};
      }
      for (size_t polar = 0U; polar < polar_count; ++polar) {
        const double diffuse_factor =
            mi_state.diffuse_precipitation_selector[hemisphere_index][polar] >
                    0.5
                ? fmax(0.0,
                       fmin(1.0,
                            mi_state.dpb_mask[hemisphere_index][polar]))
                : 1.0;
        const double conductance_factor = sqrt(diffuse_factor);
        mi_state.electron_number_flux[hemisphere_index][polar] *=
            diffuse_factor;
        mi_state.electron_energy_flux[hemisphere_index][polar] *=
            diffuse_factor;
        double target_pedersen;
        double target_hall;
        if (gamera_mi_combine_conductance(
                mi_state.euv_pedersen_conductance[hemisphere_index][polar],
                mi_state.euv_hall_conductance[hemisphere_index][polar],
                conductance_factor *
                    mi_state.auroral_pedersen_conductance[hemisphere_index]
                                                        [polar],
                conductance_factor *
                    mi_state.auroral_hall_conductance[hemisphere_index]
                                                   [polar],
                mi_config.minimum_pedersen_siemens,
                mi_config.minimum_hall_siemens,
                mi_config.maximum_hall_to_pedersen,
                &target_pedersen, &target_hall) != 0) {
          return -1;
        }
        if (!mi_state.conductance_ready[hemisphere_index]) {
          mi_state.pedersen_conductance[hemisphere_index][polar] =
              target_pedersen;
          mi_state.hall_conductance[hemisphere_index][polar] = target_hall;
        } else {
          mi_state.pedersen_conductance[hemisphere_index][polar] +=
              smoothing_alpha *
              (target_pedersen -
               mi_state.pedersen_conductance[hemisphere_index][polar]);
          mi_state.hall_conductance[hemisphere_index][polar] +=
              smoothing_alpha *
              (target_hall -
               mi_state.hall_conductance[hemisphere_index][polar]);
        }
      }
      mi_state.conductance_ready[hemisphere_index] = 1;
    }
  }
  return 0;
}

static int interpolate_potential(gamera_no_vec3 point, double *potential_v) {
  if (potential_v == NULL || !mi_state.solution_ready) {
    return -1;
  }
  const double radius = vector_norm(point);
  if (!(radius >= mi_config.ionosphere_radius_re)) {
    return -1;
  }
  const int hemisphere = point.value[2] >= 0.0 ? MI_NORTH_INDEX
                                               : MI_SOUTH_INDEX;
  const double polar_colatitude =
      acos(fmax(-1.0, fmin(1.0, fabs(point.value[2]) / radius)));
  const double ionosphere_colatitude = gamera_mi_mapped_colatitude(
      polar_colatitude, radius, mi_config.ionosphere_radius_re);
  if (!isfinite(ionosphere_colatitude)) {
    return -1;
  }
  if (ionosphere_colatitude >= mi_state.maximum_colatitude) {
    *potential_v = 0.0;
    return 0;
  }
  const size_t np = (size_t)mi_config.longitude_count;
  const size_t nt = (size_t)mi_config.colatitude_count;
  const double theta_logical = ionosphere_colatitude /
      mi_state.maximum_colatitude * (double)(nt - 1U);
  int theta_base;
  double theta_fraction;
  interpolation_coordinate(theta_logical, (int)nt, &theta_base,
                           &theta_fraction);
  double phi = atan2(point.value[1], point.value[0]);
  if (phi < 0.0) {
    phi += 2.0 * PI;
  }
  const double phi_logical = phi * (double)np / (2.0 * PI);
  const int phi_base = (int)floor(phi_logical) % (int)np;
  const int phi_next = (phi_base + 1) % (int)np;
  const double phi_fraction = phi_logical - floor(phi_logical);
  const double *potential = mi_state.potential[hemisphere];
  const double lower =
      (1.0 - phi_fraction) *
          potential[(size_t)theta_base * np + (size_t)phi_base] +
      phi_fraction *
          potential[(size_t)theta_base * np + (size_t)phi_next];
  const double upper =
      (1.0 - phi_fraction) *
          potential[(size_t)(theta_base + 1) * np + (size_t)phi_base] +
      phi_fraction *
          potential[(size_t)(theta_base + 1) * np + (size_t)phi_next];
  *potential_v =
      (1.0 - theta_fraction) * lower + theta_fraction * upper;
  return isfinite(*potential_v) ? 0 : -1;
}

static int edge_endpoints(const gamera_no_grid *grid, int direction,
                          size_t i, size_t j, size_t k,
                          gamera_no_vec3 endpoint[2]) {
  if (grid == NULL || endpoint == NULL || direction < 0 ||
      direction >= GAMERA_NO_DIM || i >= grid->vertex_extent[0] ||
      j >= grid->vertex_extent[1] || k >= grid->vertex_extent[2]) {
    return -1;
  }
  size_t end[3] = {i, j, k};
  ++end[direction];
  if (end[direction] >= grid->vertex_extent[direction]) {
    return -1;
  }
  endpoint[0] = grid->vertex[
      gamera_no_index3(grid->vertex_extent, i, j, k)];
  endpoint[1] = grid->vertex[gamera_no_index3(
      grid->vertex_extent, end[0], end[1], end[2])];
  return 0;
}

static int cache_boundary_emf(void) {
  const size_t wall_i = (size_t)is;
  const double potential_norm = norm_config.u_Norm * norm_config.B_Norm *
                                norm_config.x_Norm;
  int failed = !(potential_norm > 0.0);
  double local_maximum = 0.0;
  for (int direction = GAMERA_NO_J;
       !failed && direction <= GAMERA_NO_K; ++direction) {
    double direction_maximum = 0.0;
#pragma omp parallel for collapse(2) reduction(| : failed) \
    reduction(max : direction_maximum) schedule(static)
    for (size_t j = 0; j < mi_state.grid->edge[direction].extent[1]; ++j) {
      for (size_t k = 0; k < mi_state.grid->edge[direction].extent[2]; ++k) {
        const size_t edge = gamera_no_index3(
            mi_state.grid->edge[direction].extent, wall_i, j, k);
        if (mi_state.grid->edge[direction].valid[edge] == 0U) {
          continue;
        }
        gamera_no_vec3 endpoint[2];
        double potential[2];
        if (edge_endpoints(mi_state.grid, direction, wall_i, j, k,
                           endpoint) != 0 ||
            interpolate_potential(endpoint[0], &potential[0]) != 0 ||
            interpolate_potential(endpoint[1], &potential[1]) != 0) {
          failed = 1;
          continue;
        }
        mi_state.saved_emf[direction][edge] =
            -(potential[1] - potential[0]) / potential_norm;
        direction_maximum =
            fmax(direction_maximum,
                 fabs(mi_state.saved_emf[direction][edge]));
      }
    }
    local_maximum = fmax(local_maximum, direction_maximum);
  }
  return gamera_mi_collective_maximum_if_all_valid(
      MPI_COMM_WORLD, failed, local_maximum,
      &mi_state.maximum_cached_emf);
}

static int mi_nearly_equal(double first, double second) {
  const double scale = fmax(1.0, fmax(fabs(first), fabs(second)));
  return isfinite(first) && isfinite(second) &&
         fabs(first - second) <= 256.0 * DBL_EPSILON * scale;
}

int gamera_mi_coupling_export_restart(
    struct gamera_mi_restart_state **state) {
  if (state == NULL) {
    return -1;
  }
  *state = NULL;
  int local_failed =
      !mi_config.enabled || !mi_state.prepared || !mi_state.solution_ready ||
      !mi_north_snapshot.available || !mi_north_snapshot.ready ||
      mi_state.update_count <= 0 || mi_state.held_apply_count < 0 ||
      mi_north_snapshot.generation != (uint64_t)mi_state.update_count ||
      mi_north_snapshot.theta_points !=
          (size_t)mi_config.colatitude_count ||
      mi_north_snapshot.azimuth_points !=
          (size_t)mi_config.longitude_count ||
      (rank == 0 && mi_config.electron_precipitation_enabled &&
       (!mi_state.conductance_ready[MI_NORTH_INDEX] ||
        !mi_state.conductance_ready[MI_SOUTH_INDEX]));
  int global_failed = 0;
  if (MPI_Allreduce(&local_failed, &global_failed, 1, MPI_INT, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS ||
      global_failed) {
    return -1;
  }
  gamera_mi_restart_state *candidate = NULL;
  local_failed = gamera_mi_restart_create(
                     (size_t)mi_config.colatitude_count,
                     (size_t)mi_config.longitude_count, &candidate) != 0;
  if (MPI_Allreduce(&local_failed, &global_failed, 1, MPI_INT, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS ||
      global_failed) {
    gamera_mi_restart_destroy(candidate);
    return -1;
  }
  /* next_diagnostic_s advances on rank zero only because only rank zero
   * writes diagnostics.  Normalize it together with the rank-local held-EMF
   * audit window so every rank hashes the same detached restart. */
  double local_maximum[3] = {mi_state.maximum_prehold_difference,
                             mi_state.maximum_posthold_difference,
                             mi_state.next_diagnostic_s};
  double global_maximum[3] = {0.0, 0.0, 0.0};
  long long local_apply_count = mi_state.held_apply_count;
  long long global_apply_count = 0;
  if (MPI_Allreduce(local_maximum, global_maximum, 3, MPI_DOUBLE, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Allreduce(&local_apply_count, &global_apply_count, 1,
                    MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS ||
      global_apply_count < 0) {
    gamera_mi_restart_destroy(candidate);
    return -1;
  }
  const size_t polar_count = candidate->theta_points *
                             candidate->azimuth_points;
  memcpy(candidate->potential_v, mi_state.potential[MI_NORTH_INDEX],
         polar_count * sizeof(double));
  memcpy(candidate->potential_v + polar_count,
         mi_state.potential[MI_SOUTH_INDEX],
         polar_count * sizeof(double));
  if (rank == 0 && mi_config.electron_precipitation_enabled) {
    candidate->precipitation_state_ready = 1U;
    for (int hemisphere = 0; hemisphere < MI_HEMISPHERE_COUNT;
         ++hemisphere) {
      candidate->dpb_initialized[hemisphere] =
          (uint32_t)(mi_state.dpb_state[hemisphere].initialized != 0);
      candidate->dpb_radius_deg[hemisphere] =
          mi_state.dpb_state[hemisphere].radius_deg;
      memcpy(candidate->pedersen_siemens +
                 (size_t)hemisphere * polar_count,
             mi_state.pedersen_conductance[hemisphere],
             polar_count * sizeof(double));
      memcpy(candidate->hall_siemens + (size_t)hemisphere * polar_count,
             mi_state.hall_conductance[hemisphere],
             polar_count * sizeof(double));
    }
  }
  uint32_t precipitation_metadata[3] = {
      candidate->precipitation_state_ready, candidate->dpb_initialized[0],
      candidate->dpb_initialized[1]};
  if (polar_count > (size_t)INT_MAX / MI_HEMISPHERE_COUNT ||
      MPI_Bcast(precipitation_metadata, 3, MPI_UINT32_T, 0,
                MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(candidate->dpb_radius_deg, MI_HEMISPHERE_COUNT, MPI_DOUBLE, 0,
                MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(candidate->pedersen_siemens,
                (int)(MI_HEMISPHERE_COUNT * polar_count), MPI_DOUBLE, 0,
                MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(candidate->hall_siemens,
                (int)(MI_HEMISPHERE_COUNT * polar_count), MPI_DOUBLE, 0,
                MPI_COMM_WORLD) != MPI_SUCCESS) {
    gamera_mi_restart_destroy(candidate);
    return -1;
  }
  candidate->precipitation_state_ready = precipitation_metadata[0];
  candidate->dpb_initialized[0] = precipitation_metadata[1];
  candidate->dpb_initialized[1] = precipitation_metadata[2];
  candidate->solution_ready = 1U;
  candidate->snapshot_ready = 1U;
  candidate->next_update_s = mi_state.next_update_s;
  candidate->next_diagnostic_s = global_maximum[2];
  candidate->update_count = (uint64_t)mi_state.update_count;
  candidate->snapshot_epoch_s = mi_north_snapshot.epoch_s;
  candidate->snapshot_generation = mi_north_snapshot.generation;
  candidate->maximum_cached_emf = mi_state.maximum_cached_emf;
  candidate->maximum_prehold_difference = global_maximum[0];
  candidate->maximum_posthold_difference = global_maximum[1];
  candidate->held_apply_count = (uint64_t)global_apply_count;
  if (memcmp(mi_north_snapshot.north_potential_v,
             mi_state.potential[MI_NORTH_INDEX],
             polar_count * sizeof(double)) != 0 ||
      gamera_mi_restart_validate(candidate) != 0) {
    gamera_mi_restart_destroy(candidate);
    return -1;
  }
  *state = candidate;
  return 0;
}

int gamera_mi_coupling_restore_restart(
    const struct gamera_mi_restart_state *state, double restart_time_code) {
  int local_failed =
      !mi_config.enabled || !mi_state.prepared || mi_state.solution_ready ||
      gamera_mi_restart_validate(state) != 0 ||
      !isfinite(restart_time_code) || restart_time_code < 0.0 ||
      !(norm_config.Time_Norm > 0.0);
  if (!local_failed) {
    local_failed =
        state->theta_points != (size_t)mi_config.colatitude_count ||
        state->azimuth_points != (size_t)mi_config.longitude_count ||
        state->update_count > (uint64_t)INT_MAX ||
        state->held_apply_count > (uint64_t)LLONG_MAX ||
        state->precipitation_state_ready !=
            (uint32_t)(mi_config.electron_precipitation_enabled != 0) ||
        (mi_config.electron_precipitation_enabled &&
         mi_config.hybrid_dpb_enabled &&
         (state->dpb_initialized[0] == 0U ||
          state->dpb_initialized[1] == 0U));
  }
  int global_failed = 0;
  if (MPI_Allreduce(&local_failed, &global_failed, 1, MPI_INT, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS ||
      global_failed) {
    return -1;
  }
  const double restart_time_s = restart_time_code * norm_config.Time_Norm;
  const double tolerance_s =
      512.0 * DBL_EPSILON * fmax(1.0, fabs(restart_time_s));
  local_failed =
      !isfinite(restart_time_s) ||
      state->snapshot_epoch_s > restart_time_s + tolerance_s ||
      !(state->next_update_s > restart_time_s + tolerance_s) ||
      !mi_nearly_equal(state->next_update_s,
                       state->snapshot_epoch_s +
                           mi_config.coupling_interval_s) ||
      (mi_config.diagnostics_enabled &&
       !(state->next_diagnostic_s > restart_time_s + tolerance_s)) ||
      (!mi_config.diagnostics_enabled && state->next_diagnostic_s != 0.0);
  if (MPI_Allreduce(&local_failed, &global_failed, 1, MPI_INT, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS ||
      global_failed) {
    return -1;
  }
  const size_t polar_count = state->theta_points * state->azimuth_points;
  memcpy(mi_state.potential[MI_NORTH_INDEX], state->potential_v,
         polar_count * sizeof(double));
  memcpy(mi_state.potential[MI_SOUTH_INDEX],
         state->potential_v + polar_count, polar_count * sizeof(double));
  for (int hemisphere = 0; hemisphere < MI_HEMISPHERE_COUNT;
       ++hemisphere) {
    memcpy(mi_state.pedersen_conductance[hemisphere],
           state->pedersen_siemens + (size_t)hemisphere * polar_count,
           polar_count * sizeof(double));
    memcpy(mi_state.hall_conductance[hemisphere],
           state->hall_siemens + (size_t)hemisphere * polar_count,
           polar_count * sizeof(double));
    mi_state.conductance_ready[hemisphere] =
        state->precipitation_state_ready != 0U;
    mi_state.dpb_state[hemisphere].initialized =
        state->dpb_initialized[hemisphere] != 0U;
    mi_state.dpb_state[hemisphere].radius_deg =
        state->dpb_radius_deg[hemisphere];
  }
  memcpy(mi_north_snapshot.north_potential_v, state->potential_v,
         polar_count * sizeof(double));
  mi_state.next_update_s = state->next_update_s;
  mi_state.next_diagnostic_s = state->next_diagnostic_s;
  mi_state.update_count = (int)state->update_count;
  mi_state.maximum_prehold_difference =
      state->maximum_prehold_difference;
  mi_state.maximum_posthold_difference =
      state->maximum_posthold_difference;
  mi_state.held_apply_count = (long long)state->held_apply_count;
  mi_north_snapshot.epoch_s = state->snapshot_epoch_s;
  mi_north_snapshot.generation = state->snapshot_generation;
  mi_north_snapshot.ready = 1;
  mi_state.solution_ready = 1;
  if (cache_boundary_emf() != 0 ||
      !mi_nearly_equal(mi_state.maximum_cached_emf,
                       state->maximum_cached_emf)) {
    mi_state.solution_ready = 0;
    mi_north_snapshot.ready = 0;
    return -1;
  }
  return 0;
}

int gamera_mi_coupling_prepare(
    const gamera_no_grid *grid,
    const gamera_no_background_data *background) {
  if (!mi_config.enabled) {
    return 0;
  }
  if (grid == NULL || background == NULL || mi_state.prepared ||
      !(norm_config.x_Norm > 0.0) || !(norm_config.B_Norm > 0.0) ||
      !(norm_config.u_Norm > 0.0) || !(norm_config.rho_Norm > 0.0) ||
      !(norm_config.p_Norm > 0.0)) {
    return -1;
  }
  memset(&mi_state, 0, sizeof(mi_state));
  mi_state.grid = grid;
  mi_state.background = background;
  const size_t local_count = (size_t)config.nj * (size_t)config.nk;
  const size_t polar_count = (size_t)mi_config.longitude_count *
                             (size_t)mi_config.colatitude_count;
  int local_allocation_failed =
      local_count > SIZE_MAX / MI_SOURCE_COUNT ||
      MI_SOURCE_COUNT * local_count > SIZE_MAX / (size_t)size ||
      allocate_doubles(MI_SOURCE_COUNT * local_count,
                       &mi_state.local_source) != 0 ||
      allocate_doubles(MI_SOURCE_COUNT * local_count * (size_t)size,
                       &mi_state.gathered_source) != 0;
  if (!local_allocation_failed) {
    for (int hemisphere = 0; hemisphere < MI_HEMISPHERE_COUNT; ++hemisphere) {
      if (allocate_doubles(polar_count, &mi_state.fac[hemisphere]) != 0 ||
          allocate_doubles(polar_count, &mi_state.potential[hemisphere]) !=
              0 ||
          allocate_doubles(
              polar_count,
              &mi_state.electron_number_flux[hemisphere]) != 0 ||
          allocate_doubles(
              polar_count,
              &mi_state.electron_energy_flux[hemisphere]) != 0 ||
          allocate_doubles(
              polar_count,
              &mi_state.diffuse_precipitation_selector[hemisphere]) != 0 ||
          allocate_doubles(
              polar_count,
              &mi_state.euv_pedersen_conductance[hemisphere]) != 0 ||
          allocate_doubles(
              polar_count,
              &mi_state.euv_hall_conductance[hemisphere]) != 0 ||
          allocate_doubles(
              polar_count,
              &mi_state.auroral_pedersen_conductance[hemisphere]) != 0 ||
          allocate_doubles(
              polar_count,
              &mi_state.auroral_hall_conductance[hemisphere]) != 0 ||
          allocate_doubles(
              polar_count,
              &mi_state.pedersen_conductance[hemisphere]) != 0 ||
          allocate_doubles(
              polar_count,
              &mi_state.hall_conductance[hemisphere]) != 0 ||
          allocate_doubles(
              polar_count,
              &mi_state.dpb_mask[hemisphere]) != 0 ||
          allocate_doubles(
              (size_t)mi_config.longitude_count,
              &mi_state.dpb_boundary_latitude_deg[hemisphere]) != 0) {
        local_allocation_failed = 1;
        break;
      }
    }
  }
  if (!local_allocation_failed &&
      mi_north_snapshot_allocate(
          (size_t)mi_config.colatitude_count,
          (size_t)mi_config.longitude_count) != 0) {
    local_allocation_failed = 1;
  }
  /* Only rank zero maps and writes FAC diagnostics.  Keeping the audit
   * storage rank-local avoids replicating six floating-point and three
   * integer polar fields on every MPI rank. */
  if (!local_allocation_failed && mi_config.diagnostics_enabled &&
      !analysis_output_enabled && rank == 0) {
    for (int hemisphere = 0; hemisphere < MI_HEMISPHERE_COUNT; ++hemisphere) {
      if (allocate_doubles(polar_count,
                           &mi_state.patch_fac_delta[hemisphere]) != 0 ||
          allocate_doubles(polar_count,
                           &mi_state.patch_yang_weight[hemisphere]) != 0 ||
          allocate_ints(polar_count, &mi_state.patch_owner[hemisphere]) != 0) {
        local_allocation_failed = 1;
        break;
      }
      for (int patch = 0; patch < MI_PATCH_COUNT; ++patch) {
        if (allocate_doubles(
                polar_count, &mi_state.patch_fac[hemisphere][patch]) != 0 ||
            allocate_doubles(
                polar_count, &mi_state.patch_margin[hemisphere][patch]) != 0 ||
            allocate_ints(
                polar_count, &mi_state.patch_valid[hemisphere][patch]) != 0) {
          local_allocation_failed = 1;
          break;
        }
      }
      if (local_allocation_failed) {
        break;
      }
    }
  }
  if (!local_allocation_failed) {
    for (int direction = GAMERA_NO_J; direction <= GAMERA_NO_K; ++direction) {
      if (allocate_doubles(
              gamera_no_element_count3(grid->edge[direction].extent),
              &mi_state.saved_emf[direction]) != 0) {
        local_allocation_failed = 1;
        break;
      }
    }
  }
  int global_allocation_failed = 0;
  if (MPI_Allreduce(&local_allocation_failed, &global_allocation_failed, 1,
                    MPI_INT, MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS ||
      global_allocation_failed) {
    gamera_mi_coupling_finalize();
    return -1;
  }

  const int radial_coord = mi_config.current_shell / config.ni;
  double radius_sum = 0.0;
  long long radius_count = 0;
  if (proc_coords[0] == radial_coord) {
    const size_t i = (size_t)(is + mi_config.current_shell % config.ni);
    for (int local_j = 0; local_j < config.nj; ++local_j) {
      for (int local_k = 0; local_k < config.nk; ++local_k) {
        const size_t cell = gamera_no_index3(
            grid->cell_extent, i, (size_t)(js + local_j),
            (size_t)(ks + local_k));
        radius_sum += vector_norm(grid->cell[cell].centroid);
        ++radius_count;
      }
    }
  }
  double global_radius_sum = 0.0;
  long long global_radius_count = 0;
  if (MPI_Allreduce(&radius_sum, &global_radius_sum, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Allreduce(&radius_count, &global_radius_count, 1, MPI_LONG_LONG,
                    MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS ||
      global_radius_count <= 0) {
    gamera_mi_coupling_finalize();
    return -1;
  }
  mi_state.sample_radius_re =
      global_radius_sum / (double)global_radius_count;
  mi_state.maximum_colatitude =
      asin(sqrt(mi_config.ionosphere_radius_re / x1min_global));
  mi_state.sample_maximum_colatitude =
      asin(sqrt(mi_config.ionosphere_radius_re /
                mi_state.sample_radius_re));
  if (mi_north_snapshot_configure(mi_state.maximum_colatitude,
                                  mi_config.ionosphere_radius_re,
                                  mi_state.sample_radius_re) != 0) {
    gamera_mi_coupling_finalize();
    return -1;
  }
  /* A restart reconstructs the electrostatic solution from the checkpointed
   * MHD state.  Preserve the diagnostic sequence implied by physical time so
   * that the recomputed restart-time snapshot and all following snapshots do
   * not restart at 000001 or overwrite unrelated early-time output. */
  if (read_restart && time_sim > 0.0) {
    const double completed_cadences =
        floor((time_sim * norm_config.Time_Norm) /
                  mi_config.coupling_interval_s +
              1.0e-9);
    if (completed_cadences > (double)(INT32_MAX - 1)) {
      gamera_mi_coupling_finalize();
      return -1;
    }
    mi_state.update_count = (int)completed_cadences;
  }
  /* Seed the publication lineage from the same cadence-derived update index
   * used by restart diagnostics.  The recomputed restart-time solution then
   * receives the same generation it would have had without a process
   * restart (the t=0 solve is generation 1). */
  mi_north_snapshot.generation = (uint64_t)mi_state.update_count;
  mi_state.prepared = 1;
  mi_state.next_update_s = 0.0;
  mi_state.next_diagnostic_s = 0.0;
  log_info("Prepared electrostatic M-I coupling: FAC shell=%d r=%.9g RE, "
           "ionosphere=%dx%d per hemisphere to %.6g deg MLAT "
           "(automatic angular oversampling %.6g), "
           "SigmaP=%.6g S SigmaH=%.6g S, electron precipitation=%s "
           "(Hall model %d, beta=%.7g, EUV F10.7=%.6g sfu, "
           "sunward=(%.6g,%.6g,%.6g), "
           "hybrid DPB=%s, "
           "conductance feedback=%s), "
           "coupling cadence=%.6g s, "
           "diagnostic cadence=%.6g s",
           mi_config.current_shell + 1, mi_state.sample_radius_re,
           mi_config.longitude_count, mi_config.colatitude_count,
           90.0 - mi_state.maximum_colatitude * 180.0 / PI,
           mi_config.angular_oversampling,
           mi_config.pedersen_siemens,
           mi_config.hall_siemens,
           mi_config.electron_precipitation_enabled ? "on" : "off",
           mi_config.electron_precipitation_hall_model,
           mi_config.electron_precipitation_beta,
           mi_config.euv_f107,
           mi_config.sunward_direction[0],
           mi_config.sunward_direction[1],
           mi_config.sunward_direction[2],
           mi_config.hybrid_dpb_enabled ? "on" : "off",
           mi_config.electron_precipitation_conductance_feedback_enabled
               ? "on" : "off",
           mi_config.coupling_interval_s,
           mi_config.diagnostics_interval_s);
  return 0;
}

void gamera_mi_coupling_finalize(void) {
  mi_north_snapshot_clear();
  free(mi_state.local_source);
  free(mi_state.gathered_source);
  for (int hemisphere = 0; hemisphere < MI_HEMISPHERE_COUNT; ++hemisphere) {
    free(mi_state.fac[hemisphere]);
    free(mi_state.potential[hemisphere]);
    free(mi_state.electron_number_flux[hemisphere]);
    free(mi_state.electron_energy_flux[hemisphere]);
    free(mi_state.diffuse_precipitation_selector[hemisphere]);
    free(mi_state.euv_pedersen_conductance[hemisphere]);
    free(mi_state.euv_hall_conductance[hemisphere]);
    free(mi_state.auroral_pedersen_conductance[hemisphere]);
    free(mi_state.auroral_hall_conductance[hemisphere]);
    free(mi_state.pedersen_conductance[hemisphere]);
    free(mi_state.hall_conductance[hemisphere]);
    free(mi_state.dpb_mask[hemisphere]);
    free(mi_state.dpb_boundary_latitude_deg[hemisphere]);
    free(mi_state.patch_fac_delta[hemisphere]);
    free(mi_state.patch_yang_weight[hemisphere]);
    free(mi_state.patch_owner[hemisphere]);
    for (int patch = 0; patch < MI_PATCH_COUNT; ++patch) {
      free(mi_state.patch_fac[hemisphere][patch]);
      free(mi_state.patch_margin[hemisphere][patch]);
      free(mi_state.patch_valid[hemisphere][patch]);
    }
  }
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    free(mi_state.saved_emf[direction]);
  }
  memset(&mi_state, 0, sizeof(mi_state));
}

int gamera_mi_coupling_maybe_update(double time_code) {
  if (!mi_config.enabled) {
    return 0;
  }
  if (!mi_state.prepared || !isfinite(time_code)) {
    return -1;
  }
  const double time_s = time_code * norm_config.Time_Norm;
  if (mi_state.solution_ready &&
      time_s + 2.0e-12 * fmax(1.0, fabs(time_s)) < mi_state.next_update_s) {
    return 0;
  }
  const double update_wall_start = MPI_Wtime();
  const double reduce_wall_start = MPI_Wtime();
  double held_local[2] = {mi_state.maximum_prehold_difference,
                          mi_state.maximum_posthold_difference};
  double held_global[2] = {0.0, 0.0};
  long long held_apply_count = 0;
  if (MPI_Reduce(held_local, held_global, 2, MPI_DOUBLE, MPI_MAX, 0,
                 MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Reduce(&mi_state.held_apply_count, &held_apply_count, 1,
                 MPI_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD) != MPI_SUCCESS) {
    return -1;
  }
  const double reduce_wall_s = MPI_Wtime() - reduce_wall_start;
  mi_state.maximum_prehold_difference = 0.0;
  mi_state.maximum_posthold_difference = 0.0;
  mi_state.held_apply_count = 0;
  const double gather_wall_start = MPI_Wtime();
  if (gather_source() != 0) {
    return -1;
  }
  const double gather_wall_s = MPI_Wtime() - gather_wall_start;
  const double map_wall_start = MPI_Wtime();
  int update_status =
      rank == 0 ? map_mhd_source_to_ionosphere(time_s) : 0;
  const double map_wall_s = MPI_Wtime() - map_wall_start;
  const double status_bcast_wall_start = MPI_Wtime();
  if (MPI_Bcast(&update_status, 1, MPI_INT, 0, MPI_COMM_WORLD) !=
          MPI_SUCCESS ||
      update_status != 0) {
    return -1;
  }
  const double status_bcast_wall_s =
      MPI_Wtime() - status_bcast_wall_start;
  const size_t polar_count = (size_t)mi_config.longitude_count *
                             (size_t)mi_config.colatitude_count;
  gamera_mi_solver_stats stats[2] = {{0}};
  double solve_wall_s[2] = {0.0, 0.0};
  double solve_total_wall_s = 0.0;
  int diagnostic_status = 0;
  if (rank == 0) {
    const double solve_total_start = MPI_Wtime();
    for (int hemisphere_index = 0; hemisphere_index < 2;
         ++hemisphere_index) {
      const gamera_mi_solver_config solver = {
          .longitude_count = (size_t)mi_config.longitude_count,
          .colatitude_count = (size_t)mi_config.colatitude_count,
          .maximum_colatitude = mi_state.maximum_colatitude,
          .ionosphere_radius_m =
              mi_config.ionosphere_radius_re * norm_config.x_Norm,
          .pedersen_siemens = mi_config.pedersen_siemens,
          .hall_siemens = mi_config.hall_siemens,
          .low_latitude_potential_v = 0.0,
          .hemisphere = hemisphere_index == MI_NORTH_INDEX
                            ? GAMERA_MI_NORTH
                            : GAMERA_MI_SOUTH,
          .maximum_iterations = mi_config.maximum_iterations,
          .relative_tolerance = mi_config.relative_tolerance,
          .absolute_tolerance = mi_config.absolute_tolerance};
      const double solve_start = MPI_Wtime();
      const int solve_status =
          mi_config.electron_precipitation_enabled &&
                  mi_config.electron_precipitation_conductance_feedback_enabled
              ? gamera_mi_solve_conductance_tensor(
                    &solver, mi_state.fac[hemisphere_index],
                    mi_state.pedersen_conductance[hemisphere_index],
                    mi_state.hall_conductance[hemisphere_index],
                    mi_state.potential[hemisphere_index],
                    &stats[hemisphere_index])
              : gamera_mi_solve_constant_pedersen(
                    &solver, mi_state.fac[hemisphere_index],
                    mi_state.potential[hemisphere_index],
                    &stats[hemisphere_index]);
      solve_wall_s[hemisphere_index] = MPI_Wtime() - solve_start;
      if (solve_status != 0) {
        log_error("M-I potential solve failed in %s hemisphere: iter=%d "
                  "residual %.6e -> %.6e",
                  hemisphere_index == MI_NORTH_INDEX ? "North" : "South",
                  stats[hemisphere_index].iterations,
                  stats[hemisphere_index].initial_residual,
                  stats[hemisphere_index].final_residual);
        update_status = -1;
        break;
      }
    }
    solve_total_wall_s = MPI_Wtime() - solve_total_start;
  }
  if (MPI_Bcast(&update_status, 1, MPI_INT, 0, MPI_COMM_WORLD) !=
          MPI_SUCCESS ||
      update_status != 0) {
    return -1;
  }
  const double potential_bcast_wall_start = MPI_Wtime();
  for (int hemisphere_index = 0; hemisphere_index < 2;
       ++hemisphere_index) {
    if (MPI_Bcast(mi_state.potential[hemisphere_index], (int)polar_count,
                  MPI_DOUBLE, 0, MPI_COMM_WORLD) != MPI_SUCCESS) {
      return -1;
    }
  }
  const double potential_bcast_wall_s =
      MPI_Wtime() - potential_bcast_wall_start;
  mi_state.solution_ready = 1;
  const double cache_wall_start = MPI_Wtime();
  if (cache_boundary_emf() != 0) {
    return -1;
  }
  const double cache_wall_s = MPI_Wtime() - cache_wall_start;
  const double update_core_wall_s = MPI_Wtime() - update_wall_start;
  ++mi_state.update_count;
  const double cadence_tolerance_s =
      2.0e-12 * fmax(1.0, fabs(time_s));
  const double completed_cadences = floor(
      (time_s + cadence_tolerance_s) / mi_config.coupling_interval_s);
  mi_state.next_update_s =
      (completed_cadences + 1.0) * mi_config.coupling_interval_s;
  if (rank == 0) {
    double cpcp[2] = {0.0, 0.0};
    double fac_minimum[2] = {0.0, 0.0};
    double fac_maximum[2] = {0.0, 0.0};
    double fac_rms[2] = {0.0, 0.0};
    for (int hemisphere_index = 0; hemisphere_index < 2;
         ++hemisphere_index) {
      double minimum = mi_state.potential[hemisphere_index][0];
      double maximum = minimum;
      fac_minimum[hemisphere_index] = mi_state.fac[hemisphere_index][0];
      fac_maximum[hemisphere_index] = mi_state.fac[hemisphere_index][0];
      double fac_square_sum = 0.0;
      for (size_t index = 1; index < polar_count; ++index) {
        minimum = fmin(minimum, mi_state.potential[hemisphere_index][index]);
        maximum = fmax(maximum, mi_state.potential[hemisphere_index][index]);
      }
      for (size_t index = 0; index < polar_count; ++index) {
        const double fac_value = mi_state.fac[hemisphere_index][index];
        fac_minimum[hemisphere_index] =
            fmin(fac_minimum[hemisphere_index], fac_value);
        fac_maximum[hemisphere_index] =
            fmax(fac_maximum[hemisphere_index], fac_value);
        fac_square_sum += fac_value * fac_value;
      }
      fac_rms[hemisphere_index] =
          sqrt(fac_square_sum / (double)polar_count);
      cpcp[hemisphere_index] = (maximum - minimum) * 1.0e-3;
    }
    if (mi_config.diagnostics_enabled &&
        time_s + cadence_tolerance_s >= mi_state.next_diagnostic_s) {
      if (write_mi_snapshot(time_code, time_s) != 0) {
        log_error("Failed to write M-I diagnostic snapshot %d",
                  mi_state.update_count);
        diagnostic_status = -1;
      } else {
        const double completed_diagnostic_cadences = floor(
            (time_s + cadence_tolerance_s) /
            mi_config.diagnostics_interval_s);
        mi_state.next_diagnostic_s =
            (completed_diagnostic_cadences + 1.0) *
            mi_config.diagnostics_interval_s;
      }
    }
    log_info("M-I update %d at %.6f s: CPCP N/S=%.6g/%.6g kV, "
             "solver iterations N/S=%d/%d, solver wall N/S/total="
             "%.6f/%.6f/%.6f s, update wall "
             "reduce/gather/map/status_bcast/potential_bcast/cache/core="
             "%.6f/%.6f/%.6f/%.6f/%.6f/%.6f/%.6f s, next=%.6f s",
             mi_state.update_count, time_s, cpcp[0], cpcp[1],
             stats[0].iterations, stats[1].iterations, solve_wall_s[0],
             solve_wall_s[1], solve_total_wall_s, reduce_wall_s,
             gather_wall_s, map_wall_s, status_bcast_wall_s,
             potential_bcast_wall_s, cache_wall_s, update_core_wall_s,
             mi_state.next_update_s);
    if (mi_config.electron_precipitation_enabled &&
        mi_config.hybrid_dpb_enabled) {
      log_info("M-I hybrid DPB N/S nightside=%.4g/%.4g MLAT, "
               "dayside=%.4g/%.4g MLAT, radius=%.4g/%.4g deg, "
               "evidence=%.3g/%.3g",
               mi_state.dpb_stats[0].nightside_boundary_deg,
               mi_state.dpb_stats[1].nightside_boundary_deg,
               mi_state.dpb_stats[0].dayside_boundary_deg,
               mi_state.dpb_stats[1].dayside_boundary_deg,
               mi_state.dpb_stats[0].filtered_radius_deg,
               mi_state.dpb_stats[1].filtered_radius_deg,
               mi_state.dpb_stats[0].evidence_fraction,
               mi_state.dpb_stats[1].evidence_fraction);
    }
    log_info("M-I FAC parallel N/S [min,max,rms]="
             "[%.6g,%.6g,%.6g]/[%.6g,%.6g,%.6g] uA/m^2; "
             "held EMF max=%.6e code, previous interval applies=%lld "
             "pre/post overwrite max differences=%.6e/%.6e",
             1.0e6 * fac_minimum[0], 1.0e6 * fac_maximum[0],
             1.0e6 * fac_rms[0], 1.0e6 * fac_minimum[1],
             1.0e6 * fac_maximum[1], 1.0e6 * fac_rms[1],
             mi_state.maximum_cached_emf, held_apply_count,
             held_global[0], held_global[1]);
  }
  if (MPI_Bcast(&diagnostic_status, 1, MPI_INT, 0, MPI_COMM_WORLD) !=
          MPI_SUCCESS ||
      diagnostic_status != 0) {
    return -1;
  }
  if (mi_north_snapshot_publish(
          time_s, mi_state.potential[MI_NORTH_INDEX]) != 0) {
    return -1;
  }
  return 1;
}

double gamera_mi_coupling_limit_timestep(double time_code,
                                         double proposed_dt_code) {
  if (!mi_config.enabled || !mi_state.prepared || !mi_state.solution_ready ||
      !isfinite(time_code) || !isfinite(proposed_dt_code) ||
      !(proposed_dt_code > 0.0) || !(norm_config.Time_Norm > 0.0)) {
    return proposed_dt_code;
  }
  const double time_s = time_code * norm_config.Time_Norm;
  const double remaining_s = mi_state.next_update_s - time_s;
  const double tolerance_s =
      2.0e-12 * fmax(1.0, fabs(mi_state.next_update_s));
  if (remaining_s > tolerance_s &&
      proposed_dt_code * norm_config.Time_Norm > remaining_s) {
    return remaining_s / norm_config.Time_Norm;
  }
  return proposed_dt_code;
}

int gamera_mi_coupling_apply_held_emf(
    gamera_no_storage *storage, const gamera_no_grid *grid,
    const size_t active_lower[3], const size_t active_upper[3]) {
  (void)active_upper;
  if (!mi_config.enabled || !mi_state.solution_ready) {
    return 0;
  }
  if (storage == NULL || grid != mi_state.grid || active_lower == NULL) {
    return -1;
  }
  if (proc_coords[0] != 0) {
    return 0;
  }
  const size_t wall_i = active_lower[GAMERA_NO_I];
  double maximum_prehold_difference = 0.0;
  double maximum_posthold_difference = 0.0;
  for (int direction = GAMERA_NO_J; direction <= GAMERA_NO_K; ++direction) {
#pragma omp parallel for collapse(2) \
    reduction(max : maximum_prehold_difference, maximum_posthold_difference) \
    schedule(static)
    for (size_t j = 0; j < grid->edge[direction].extent[1]; ++j) {
      for (size_t k = 0; k < grid->edge[direction].extent[2]; ++k) {
        const size_t edge = gamera_no_index3(
            grid->edge[direction].extent, wall_i, j, k);
        if (grid->edge[direction].valid[edge] != 0U) {
          maximum_prehold_difference = fmax(
              maximum_prehold_difference,
              fabs(storage->edge_emf[direction][edge] -
                   mi_state.saved_emf[direction][edge]));
          storage->edge_emf[direction][edge] =
              mi_state.saved_emf[direction][edge];
          maximum_posthold_difference = fmax(
              maximum_posthold_difference,
              fabs(storage->edge_emf[direction][edge] -
                   mi_state.saved_emf[direction][edge]));
        }
      }
    }
  }
  mi_state.maximum_prehold_difference =
      fmax(mi_state.maximum_prehold_difference, maximum_prehold_difference);
  mi_state.maximum_posthold_difference =
      fmax(mi_state.maximum_posthold_difference, maximum_posthold_difference);
  ++mi_state.held_apply_count;
  return 0;
}

int gamera_mi_coupling_ghost_velocity(gamera_no_vec3 point,
                                      gamera_no_vec3 source_velocity,
                                      gamera_no_vec3 *ghost_velocity) {
  if (ghost_velocity == NULL || !mi_state.solution_ready) {
    return -1;
  }
  const double radius = vector_norm(point);
  if (!(radius > 0.0)) {
    return -1;
  }
  const double epsilon = 1.0e-4 * fmax(1.0, radius);
  const double electric_norm = norm_config.u_Norm * norm_config.B_Norm;
  gamera_no_vec3 electric = {{0.0, 0.0, 0.0}};
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    gamera_no_vec3 plus = point;
    gamera_no_vec3 minus = point;
    plus.value[component] += epsilon;
    minus.value[component] -= epsilon;
    double plus_potential, minus_potential;
    if (interpolate_potential(plus, &plus_potential) != 0 ||
        interpolate_potential(minus, &minus_potential) != 0) {
      return -1;
    }
    electric.value[component] =
        -(plus_potential - minus_potential) /
        (2.0 * epsilon * norm_config.x_Norm * electric_norm);
  }
  gamera_no_vec3 background;
  if (problem_nonorthogonal_background_field(point, NULL, &background) != 0) {
    return -1;
  }
  const double background_square = vector_dot(background, background);
  if (!(background_square > 0.0)) {
    return -1;
  }
  gamera_no_vec3 drift = vector_cross(electric, background);
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    drift.value[component] /= background_square;
  }
  const gamera_no_vec3 radial = {{point.value[0] / radius,
                                  point.value[1] / radius,
                                  point.value[2] / radius}};
  const double drift_radial = vector_dot(drift, radial);
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    drift.value[component] -= drift_radial * radial.value[component];
  }
  /* Kaiju stores V_EB directly because its boundary object consumes a wall
   * state.  GAMERA-OP's mirrored ghost reconstruction instead needs
   * 2*V_EB-v_active to place V_EB at the physical wall. */
  return gamera_mi_compose_ghost_velocity(source_velocity, background, drift,
                                          ghost_velocity);
}

int gamera_mi_coupling_ready(void) {
  return mi_state.solution_ready;
}
#else
int gamera_mi_coupling_prepare(
    const gamera_no_grid *grid,
    const gamera_no_background_data *background) {
  (void)grid;
  (void)background;
  return 0;
}

void gamera_mi_coupling_finalize(void) {}

int gamera_mi_coupling_maybe_update(double time_code) {
  (void)time_code;
  return 0;
}

double gamera_mi_coupling_limit_timestep(double time_code,
                                         double proposed_dt_code) {
  (void)time_code;
  return proposed_dt_code;
}

int gamera_mi_coupling_apply_held_emf(
    gamera_no_storage *storage, const gamera_no_grid *grid,
    const size_t active_lower[3], const size_t active_upper[3]) {
  (void)storage;
  (void)grid;
  (void)active_lower;
  (void)active_upper;
  return 0;
}

int gamera_mi_coupling_ghost_velocity(gamera_no_vec3 point,
                                      gamera_no_vec3 source_velocity,
                                      gamera_no_vec3 *ghost_velocity) {
  (void)point;
  (void)source_velocity;
  (void)ghost_velocity;
  return -1;
}

int gamera_mi_coupling_ready(void) { return 0; }

int gamera_mi_coupling_export_restart(
    struct gamera_mi_restart_state **state) {
  if (state != NULL) {
    *state = NULL;
  }
  return -1;
}

int gamera_mi_coupling_restore_restart(
    const struct gamera_mi_restart_state *state, double restart_time_code) {
  (void)state;
  (void)restart_time_code;
  return -1;
}
#endif
