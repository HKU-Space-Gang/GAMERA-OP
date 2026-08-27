#ifndef GAMERA_NONORTHOGONAL_BACKGROUND_H
#define GAMERA_NONORTHOGONAL_BACKGROUND_H

#include "nonorthogonal_sweep.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Cartesian vector field evaluated at a physical point. */
typedef int (*gamera_no_background_function)(gamera_no_vec3 point,
                                             void *context,
                                             gamera_no_vec3 *field);

/*
 * Owned form of gamera_no_background_field.  The public field member can be
 * passed directly to the sweep, CT, timestep, and update kernels.  The data
 * are static after construction and are therefore safe to share between
 * OpenMP threads.
 */
typedef struct {
  gamera_no_background_field field;
  gamera_no_vec3 *cell_magnetic;
  gamera_no_vec3 *face_magnetic[GAMERA_NO_DIM];
  double *face_flux[GAMERA_NO_DIM];
  gamera_no_vec3 *edge_magnetic[GAMERA_NO_DIM];
  gamera_no_vec3 *cell_force;
} gamera_no_background_data;

/*
 * Reproduce GAMERA background.F90:AddB0 using 12-point Gauss-Legendre
 * quadrature: volume-averaged cell B0, area-averaged Cartesian face B0,
 * oriented face flux, parametric edge-averaged B0, and the volume-normalized
 * pure-B0 Maxwell-stress divergence.
 */
int gamera_no_background_create(const gamera_no_grid *grid,
                                gamera_no_background_function function,
                                void *context,
                                gamera_no_background_data *background);
void gamera_no_background_destroy(gamera_no_background_data *background);

typedef struct {
  /* Physical/code-unit dipole moment vector. */
  gamera_no_vec3 moment;
} gamera_no_dipole;

/* B = 3 (m dot r) r / r^5 - m / r^3. */
int gamera_no_dipole_field(gamera_no_vec3 point, void *context,
                           gamera_no_vec3 *field);

#ifdef __cplusplus
}
#endif

#endif
