#ifndef GAMERA_NONORTHOGONAL_ADVANCE_H
#define GAMERA_NONORTHOGONAL_ADVANCE_H

#include "nonorthogonal_step.h"
#include "nonorthogonal_sweep.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  gamera_no_sweep_options stress;
  gamera_no_emf_options emf;
  gamera_no_update_options update;
  /* Optional conservative source evaluated from the half-time predictor. */
  int (*cell_source)(gamera_no_vec3 point,
                     const double predicted[GAMERA_NO_FLUX_COUNT],
                     double time, void *context,
                     double rate[GAMERA_NO_FLUX_COUNT]);
  void *source_context;
  /*
   * Optional additive source with O(1) access to the full local storage-cell
   * index.  This is intended for sparse/cached couplings such as ring-current
   * pressure feedback; it is evaluated independently of cell_source and the
   * two returned conservative rates are summed.
   */
  int (*indexed_cell_source)(
      size_t cell, gamera_no_vec3 point,
      const double predicted[GAMERA_NO_FLUX_COUNT], double time,
      void *context, double rate[GAMERA_NO_FLUX_COUNT]);
  void *indexed_source_context;
  double source_time;
  /*
   * Optional physical-boundary correction after all three fluid/Maxwell
   * sweeps and before their divergence is accumulated into cell rates.
   */
  int (*fluid_flux_sync)(gamera_no_storage *storage,
                         const gamera_no_grid *grid,
                         const size_t active_lower[3],
                         const size_t active_upper[3], void *context);
  void *fluid_flux_context;
  /*
   * Optional overset-grid synchronization after all three edge-EMF sweeps
   * and before the CT curl updates any face flux.  The callback may replace
   * receiver-edge line integrals, but must not advance the magnetic state.
   */
  int (*edge_emf_sync)(gamera_no_storage *storage,
                       const gamera_no_grid *grid,
                       const size_t active_lower[3],
                       const size_t active_upper[3], void *context);
  void *edge_emf_context;
} gamera_no_advance_options;

/*
 * Complete single-fluid local-rank GAMERA step.  predictor_ratio is
 * (0.5*dt)/(time_current-time_old).  MPI exchange, physical boundaries,
 * Ring Average, and output remain responsibilities of the outer backend.
 */
int gamera_no_advance(
    gamera_no_storage *storage, const gamera_no_grid *grid,
    const size_t active_lower[3], const size_t active_upper[3],
    double predictor_ratio, double dt,
    const gamera_no_advance_options *options,
    const gamera_no_background_field *background);

#ifdef __cplusplus
}
#endif

#endif
