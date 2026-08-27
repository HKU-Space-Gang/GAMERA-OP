#ifndef GAMERA_NONORTHOGONAL_INITIALIZATION_H
#define GAMERA_NONORTHOGONAL_INITIALIZATION_H

#include "nonorthogonal_grid.h"
#include "nonorthogonal_storage.h"

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef int (*gamera_no_primitive_function)(
    gamera_no_vec3 point, void *context, gamera_no_primitive *primitive);
typedef int (*gamera_no_vector_potential_function)(
    gamera_no_vec3 point, void *context, gamera_no_vec3 *potential);

/* Initialize cell averages either at the metric centroid or by 12^3 Gauss. */
int gamera_no_initialize_primitives(
    const gamera_no_grid *grid, gamera_no_storage *storage,
    gamera_no_primitive_function function, void *context,
    bool volume_average, double gamma, double density_floor,
    double pressure_floor);

/*
 * Initialize unique face magnetic fluxes from line-integrated vector
 * potential and a discrete Stokes curl.  The result is divergence-free to
 * roundoff independently of grid orthogonality.
 */
int gamera_no_initialize_flux_from_vector_potential(
    const gamera_no_grid *grid, gamera_no_storage *storage,
    gamera_no_vector_potential_function function, void *context);

typedef struct {
  double density;
  double pressure;
  double magnetic_amplitude;
} gamera_no_orszag_tang_options;

void gamera_no_orszag_tang_fortran_defaults(
    gamera_no_orszag_tang_options *options);
/* Zhang et al. (2019), Section 4.4 / Equation (108) normalization. */
void gamera_no_orszag_tang_paper_defaults(
    gamera_no_orszag_tang_options *options);
int gamera_no_initialize_orszag_tang(
    const gamera_no_grid *grid, gamera_no_storage *storage,
    const gamera_no_orszag_tang_options *options, double gamma,
    double density_floor, double pressure_floor);

typedef struct {
  gamera_no_vec3 center;
  double radius;
  double ambient_density;
  double ambient_pressure;
  double density_ratio;
  double pressure_ratio;
  gamera_no_vec3 velocity;
  gamera_no_vec3 magnetic;
  bool volume_average;
} gamera_no_blast_options;

void gamera_no_spherical_blast_defaults(gamera_no_blast_options *options);
int gamera_no_initialize_blast(
    const gamera_no_grid *grid, gamera_no_storage *storage,
    const gamera_no_blast_options *options, double gamma,
    double density_floor, double pressure_floor);

#ifdef __cplusplus
}
#endif

#endif
