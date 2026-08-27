#ifndef GAMERA_NONORTHOGONAL_SWEEP_H
#define GAMERA_NONORTHOGONAL_SWEEP_H

#include "nonorthogonal_grid.h"
#include "nonorthogonal_storage.h"

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  double gamma;
  double density_floor;
  double pressure_floor;
  double pdm_coefficient;
  double hydro_hogs_coefficient;
  double magnetic_hogs_coefficient;
  double light_speed;
  bool use_mhd;
  bool use_hogs;
  bool use_boris;
  bool use_background;
} gamera_no_sweep_options;

/*
 * Optional split-field data.  The face arrays use the same extents as the
 * corresponding grid faces.  cell_force is the precomputed dpB0 term added
 * after the Maxwell-stress divergence.
 */
typedef struct {
  const gamera_no_vec3 *cell_magnetic;
  const gamera_no_vec3 *face_magnetic[3];
  const double *face_flux[3];
  const gamera_no_vec3 *edge_magnetic[3];
  const gamera_no_vec3 *cell_force;
} gamera_no_background_field;

typedef struct {
  double density_floor;
  double pdm_coefficient;
  double diffusion_coefficient;
  double light_speed;
  double cfl;
  double dt;
  bool use_boris;
  bool use_background;
} gamera_no_emf_options;

/*
 * Reconstruct one logical-direction set of active faces and write area-scaled
 * Reynolds and Maxwell fluxes into storage.  active_upper is exclusive for
 * cells; the sweep includes both bounding faces in direction.
 *
 * The direction being reconstructed requires four ghost cells on each side:
 * active_lower[direction] >= 4 and
 * active_upper[direction] + 3 < grid->cell_extent[direction].
 */
int gamera_no_sweep_face_fluxes(
    const gamera_no_grid *grid, gamera_no_storage *storage,
    const double *conserved, const gamera_no_vec3 *cell_magnetic,
    const double *const face_magnetic_flux[3], int direction,
    const size_t active_lower[3], const size_t active_upper[3],
    const gamera_no_sweep_options *options,
    const gamera_no_background_field *background);

/* Convert the three sets of area-scaled face fluxes into cell rates. */
int gamera_no_calculate_stress_rates(
    const gamera_no_grid *grid, gamera_no_storage *storage,
    const size_t active_lower[3], const size_t active_upper[3], bool use_mhd,
    const gamera_no_vec3 *background_cell_force);

/* Two-stage center-to-face-to-edge velocity interpolation and GetCornerB. */
int gamera_no_sweep_edge_emf(
    const gamera_no_grid *grid, gamera_no_storage *storage,
    const double *conserved, const double *const face_magnetic_flux[3],
    int edge_direction, const size_t active_lower[3],
    const size_t active_upper[3], const gamera_no_emf_options *options,
    const gamera_no_background_field *background);

#ifdef __cplusplus
}
#endif

#endif
