#ifndef GAMERA_NONORTHOGONAL_MI_DPB_H
#define GAMERA_NONORTHOGONAL_MI_DPB_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Dynamic diffuse-precipitation boundary (DPB).  The boundary is represented
 * by the Zhang-style offset small circle used by the accepted offline
 * diagnostic.  Its radius is inferred from the convection-potential edge and
 * an R2/R1 opposite-sign FAC pair, then filtered in time.  Only the nightside
 * 20--04 MLT sector is prevented from contracting poleward of 70 MLAT; the
 * dayside remains free to follow the offset-oval geometry.
 */
typedef struct {
  double absolute_potential_v;
  double adaptive_potential_fraction;
  double adaptive_potential_floor_v;
  double potential_equatorward_offset_deg;
  double fac_equatorward_offset_deg;
  double fac_absolute_floor_a_m2;
  double minimum_boundary_latitude_deg;
  double maximum_candidate_latitude_deg;
  double transition_width_deg;
  double temporal_timescale_s;
  double maximum_slew_deg_per_s;
  double quiet_initial_nightside_latitude_deg;
  double oval_center_latitude_deg;
  double oval_center_mlt_h;
  double nightside_poleward_limit_deg;
} gamera_mi_hybrid_dpb_config;

typedef struct {
  int initialized;
  double radius_deg;
} gamera_mi_hybrid_dpb_state;

typedef struct {
  double target_radius_deg;
  double filtered_radius_deg;
  double nightside_boundary_deg;
  double dayside_boundary_deg;
  double cpcp_v;
  double adaptive_potential_threshold_v;
  double fac_scale_a_m2;
  double evidence_fraction;
  int nightside_limit_active;
} gamera_mi_hybrid_dpb_stats;

gamera_mi_hybrid_dpb_config gamera_mi_hybrid_dpb_default_config(void);

/* Longitudes are uniformly spaced on [0,2*pi).  FAC uses the MHD Cartesian
 * sign; upward_fac_multiplier is -1 North and +1 South.  potential_v is the
 * most recently available ionospheric solution (one coupling cadence old),
 * avoiding an implicit DPB/conductance/potential iteration. */
int gamera_mi_hybrid_dpb_update(
    const gamera_mi_hybrid_dpb_config *config, size_t colatitude_count,
    size_t longitude_count, double maximum_colatitude_rad,
    double mapped_maximum_colatitude_rad, const double *fac_a_m2,
    double upward_fac_multiplier, const double *potential_v,
    double elapsed_s, gamera_mi_hybrid_dpb_state *state,
    double *boundary_latitude_deg, double *mask,
    gamera_mi_hybrid_dpb_stats *stats);

/* Public geometry helper used by deterministic tests and diagnostics. */
int gamera_mi_hybrid_dpb_boundary_from_radius(
    size_t longitude_count, double radius_deg,
    double center_latitude_deg, double center_mlt_h,
    double *boundary_latitude_deg);

#ifdef __cplusplus
}
#endif

#endif
