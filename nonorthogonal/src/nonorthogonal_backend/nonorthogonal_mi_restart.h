#ifndef GAMERA_NONORTHOGONAL_MI_RESTART_H
#define GAMERA_NONORTHOGONAL_MI_RESTART_H

#include <mpi.h>

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define GAMERA_MI_RESTART_SCHEMA 2U
#define GAMERA_MI_RESTART_DIGEST_BYTES 32U

/*
 * Detached electrostatic M-I restart state.  The held inner-wall EMF is a
 * deterministic derivative of the two potential arrays and is deliberately
 * rebuilt on every rank after restore rather than checkpointed as a
 * decomposition-dependent edge array.
 */
typedef struct gamera_mi_restart_state {
  uint32_t schema_version;
  uint32_t solution_ready;
  uint32_t snapshot_ready;
  size_t theta_points;
  size_t azimuth_points;
  double next_update_s;
  double next_diagnostic_s;
  uint64_t update_count;
  double snapshot_epoch_s;
  uint64_t snapshot_generation;
  double maximum_cached_emf;
  double maximum_prehold_difference;
  double maximum_posthold_difference;
  uint64_t held_apply_count;
  uint32_t precipitation_state_ready;
  uint32_t dpb_initialized[2];
  double dpb_radius_deg[2];
  /* Native solver layout: [hemisphere][theta][azimuth], North then South. */
  double *potential_v;
  /* Smoothed total conductance used by the next variable-tensor solve. */
  double *pedersen_siemens;
  double *hall_siemens;
} gamera_mi_restart_state;

int gamera_mi_restart_create(size_t theta_points, size_t azimuth_points,
                             gamera_mi_restart_state **state);
void gamera_mi_restart_destroy(gamera_mi_restart_state *state);
int gamera_mi_restart_validate(const gamera_mi_restart_state *state);
int gamera_mi_restart_digest(
    const gamera_mi_restart_state *state,
    unsigned char digest[GAMERA_MI_RESTART_DIGEST_BYTES]);

/* Append/read the strict /mi payload in the composite online-H+ sidecar. */
int gamera_mi_restart_append_hdf5(const char *path,
                                  const gamera_mi_restart_state *state);
/* Return 1 when /mi exists, 0 for a readable legacy file without it, and -1
 * when the HDF5 file itself cannot be inspected. */
int gamera_mi_restart_hdf5_present(const char *path);
int gamera_mi_restart_read_hdf5(const char *path,
                                gamera_mi_restart_state **state);

/* Root owns *state on entry; every rank owns an independent copy on return. */
int gamera_mi_restart_broadcast(MPI_Comm communicator, int root,
                                gamera_mi_restart_state **state);

#ifdef __cplusplus
}
#endif

#endif
