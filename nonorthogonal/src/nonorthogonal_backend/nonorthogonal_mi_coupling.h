#ifndef GAMERA_NONORTHOGONAL_MI_COUPLING_H
#define GAMERA_NONORTHOGONAL_MI_COUPLING_H

#include "nonorthogonal_background.h"
#include "nonorthogonal_grid.h"
#include "nonorthogonal_storage.h"

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

struct gamera_mi_restart_state;

int gamera_mi_coupling_prepare(
    const gamera_no_grid *grid,
    const gamera_no_background_data *background);
void gamera_mi_coupling_finalize(void);

/* Return 1 after an update, 0 when the held solution remains current. */
int gamera_mi_coupling_maybe_update(double time_code);

/* Shorten an MHD step so that it lands exactly on the next coupling time. */
double gamera_mi_coupling_limit_timestep(double time_code,
                                         double proposed_dt_code);

/* Reapply the cached first-active-shell tangential EMFs every MHD step. */
int gamera_mi_coupling_apply_held_emf(
    gamera_no_storage *storage, const gamera_no_grid *grid,
    const size_t active_lower[3], const size_t active_upper[3]);

/* Inner ghost velocity: reflected parallel flow and imposed ExB drift. */
int gamera_mi_coupling_ghost_velocity(gamera_no_vec3 point,
                                      gamera_no_vec3 source_velocity,
                                      gamera_no_vec3 *ghost_velocity);

int gamera_mi_coupling_ready(void);

/* Export the electrostatic scheduler/potential state and restore it before
 * the first post-restart maybe_update.  Held edge EMFs are reconstructed
 * from the restored N/S potentials on each decomposition-local grid. */
int gamera_mi_coupling_export_restart(
    struct gamera_mi_restart_state **state);
int gamera_mi_coupling_restore_restart(
    const struct gamera_mi_restart_state *state, double restart_time_code);

/*
 * Immutable, borrowed view of the last completely published North M-I
 * solution.  The polar arrays use the solver/writer's native layout:
 *
 *   potential_v[theta_index * azimuth_points + azimuth_index]
 *
 * colatitude_rad contains theta nodes from the magnetic pole toward the
 * low-latitude boundary.  longitude_rad contains periodic azimuth nodes;
 * phi=0 is 00 MLT (midnight) and increases toward 06 MLT.  The pointed-to
 * storage remains owned by the M-I coupling and must not be modified or
 * freed.  It is valid only until the next successful M-I publication or
 * gamera_mi_coupling_finalize().
 */
typedef struct gamera_mi_north_snapshot_view {
  int ready;
  double epoch_s;
  uint64_t generation;
  size_t theta_points;
  size_t azimuth_points;
  double ionosphere_radius_re;
  double trace_radius_re;
  const double *colatitude_rad;
  const double *longitude_rad;
  const double *potential_v;
} gamera_mi_north_snapshot_view;

typedef enum gamera_mi_snapshot_status {
  GAMERA_MI_SNAPSHOT_OK = 0,
  GAMERA_MI_SNAPSHOT_UNAVAILABLE = 1,
  GAMERA_MI_SNAPSHOT_UNREADY = 2,
  GAMERA_MI_SNAPSHOT_STALE = 3,
  GAMERA_MI_SNAPSHOT_INVALID_ARGUMENT = 4
} gamera_mi_snapshot_status;

/*
 * Borrow the latest view.  expected_generation=0 selects the current
 * generation; a nonzero value must match exactly or STALE is returned.
 * Every non-OK return clears view, so stale/unready data cannot be consumed
 * accidentally.  This routine performs no allocation, copy, or collective.
 */
gamera_mi_snapshot_status gamera_mi_coupling_borrow_north_snapshot(
    uint64_t expected_generation, gamera_mi_north_snapshot_view *view);

const char *gamera_mi_snapshot_status_string(
    gamera_mi_snapshot_status status);

#ifdef GAMERA_MI_SNAPSHOT_TESTING
/* Narrow lifecycle seam used only by the focused snapshot unit test. */
int gamera_mi_snapshot_test_prepare(size_t theta_points,
                                    size_t azimuth_points,
                                    double maximum_colatitude_rad,
                                    double ionosphere_radius_re,
                                    double trace_radius_re);
int gamera_mi_snapshot_test_publish(double epoch_s,
                                    const double *north_potential_v);
void gamera_mi_snapshot_test_finalize(void);
#endif

#ifdef __cplusplus
}
#endif

#endif
