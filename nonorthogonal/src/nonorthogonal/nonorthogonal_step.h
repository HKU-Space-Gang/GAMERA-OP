#ifndef GAMERA_NONORTHOGONAL_STEP_H
#define GAMERA_NONORTHOGONAL_STEP_H

#include "nonorthogonal_grid.h"
#include "nonorthogonal_storage.h"

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

double gamera_no_cell_net_flux(const gamera_no_storage *storage, size_t i,
                               size_t j, size_t k);

/*
 * Kaiju earthcmi.F90:IonFlux hard-wall guard.  At the inner radial face,
 * prevent the numerical reconstruction from injecting mass and positive
 * energy into the active domain.  Momentum and Maxwell fluxes are untouched.
 */
int gamera_no_trap_inner_outward_mass_energy_flux(
    gamera_no_storage *storage, const size_t active_lower[3],
    const size_t active_upper[3], size_t *clamped_face_count);

/* Kaiju msphutils.F90:ChillOut pressure control for Boris magnetospheres. */
int gamera_no_apply_kaiju_chillout(
    gamera_no_storage *storage, const gamera_no_grid *grid,
    const size_t active_lower[3], const size_t active_upper[3], double dt,
    double gamma, double density_floor, double pressure_floor,
    double light_speed, double low_density_threshold,
    double low_density_sound_speed, size_t *low_density_cell_count,
    size_t *high_sound_speed_cell_count, double *maximum_sound_speed,
    double *maximum_pressure);

/* Apply the grid-level E2Flux curl to the current face-flux arrays. */
int gamera_no_advance_ct(gamera_no_storage *storage, double dt);

/* CT update restricted to faces bounding an active cell box. */
int gamera_no_advance_ct_active(gamera_no_storage *storage,
                                const size_t active_lower[3],
                                const size_t active_upper[3], double dt);

/* Recover Cartesian cell B for an arbitrary matching set of face-flux arrays. */
int gamera_no_recover_magnetic_field(
    const gamera_no_grid *grid, const double *const face_flux[3],
    gamera_no_vec3 *cell_magnetic);

/* Save current gas/B/face-flux into the old time level. */
int gamera_no_save_current_as_old(gamera_no_storage *storage);

/* Single-fluid primitive AB predictor plus face-flux predictor and B recovery. */
int gamera_no_predict_storage(gamera_no_storage *storage,
                              const gamera_no_grid *grid, double ratio,
                              double gamma, double density_floor,
                              double pressure_floor);

typedef struct {
  double gamma;
  double density_floor;
  double pressure_floor;
  double light_speed;
  bool use_mhd;
  bool use_boris;
  bool use_background;
} gamera_no_update_options;

typedef struct {
  double gamma;
  double density_floor;
  double pressure_floor;
  double cfl;
  double light_speed;
  bool use_mhd;
  bool use_boris;
  bool use_background;
} gamera_no_timestep_options;

/* Local-rank CFL minimum; the caller performs the MPI minimum reduction. */
int gamera_no_local_timestep(
    const gamera_no_grid *grid, const double *conserved,
    const gamera_no_vec3 *cell_magnetic, const size_t active_lower[3],
    const size_t active_upper[3],
    const gamera_no_timestep_options *options,
    const gamera_no_vec3 *background_cell_magnetic, double *local_dt);

/* Apply precomputed hydro/Maxwell rates and edge EMF to the active box. */
int gamera_no_apply_active_update(
    gamera_no_storage *storage, const gamera_no_grid *grid,
    const size_t active_lower[3], const size_t active_upper[3], double dt,
    const gamera_no_update_options *options,
    const gamera_no_vec3 *background_cell_magnetic);

#ifdef __cplusplus
}
#endif

#endif
