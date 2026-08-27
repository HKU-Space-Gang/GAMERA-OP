#ifndef GAMERA_NONORTHOGONAL_FLUX_H
#define GAMERA_NONORTHOGONAL_FLUX_H

#include "nonorthogonal_geometry.h"

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

enum {
  GAMERA_NO_FLUX_DENSITY = 0,
  GAMERA_NO_FLUX_MOMENTUM_X,
  GAMERA_NO_FLUX_MOMENTUM_Y,
  GAMERA_NO_FLUX_MOMENTUM_Z,
  GAMERA_NO_FLUX_ENERGY,
  GAMERA_NO_FLUX_COUNT
};

typedef struct {
  double density;
  gamera_no_vec3 velocity;
  double pressure;
} gamera_no_primitive;

typedef struct {
  double conserved[GAMERA_NO_FLUX_COUNT];
  double conserved_jump[GAMERA_NO_FLUX_COUNT];
  gamera_no_vec3 velocity_jump;
  double normal_velocity[2];
} gamera_no_fluid_flux;

typedef struct {
  gamera_no_vec3 momentum;
  double alfven_diffusion_speed;
  double magnetic_pressure_sum;
} gamera_no_maxwell_flux;

/* Fortran stress.F90:GasKinFlux, summed over LEFT and RIGHT half spaces. */
int gamera_no_kinetic_fluid_flux(const gamera_no_primitive state[2],
                                 double gamma,
                                 const gamera_no_face_geometry *face,
                                 gamera_no_fluid_flux *flux);

/*
 * Fortran stress.F90:MagKinFlux, summed over LEFT and RIGHT half spaces.
 * magnetic[] is the reconstructed residual Cartesian field. face_normal_field
 * is the continuous residual face-normal B from face flux / face area.
 */
int gamera_no_kinetic_maxwell_flux(
    const gamera_no_primitive state[2], const gamera_no_vec3 magnetic[2],
    double face_normal_field, const gamera_no_face_geometry *face,
    bool use_background, gamera_no_vec3 background,
    double background_face_normal_field, bool use_boris, double light_speed,
    const double normal_velocity[2], gamera_no_maxwell_flux *flux);

/* Apply the default high-order gradient stabilization terms. */
int gamera_no_apply_hogs(gamera_no_fluid_flux *fluid,
                         gamera_no_maxwell_flux *maxwell,
                         double hydro_coefficient,
                         double magnetic_coefficient, bool use_boris,
                         double light_speed);

/*
 * Emergency interface diffusion: if the RMS reconstructed interface speed
 * reaches 1.5 times the Boris light speed, add c*(U_R-U_L) cell-centered
 * diffusion.
 */
int gamera_no_apply_emergency_interface_diffusion(
    gamera_no_fluid_flux *fluid, const gamera_no_primitive interface[2],
    const double lower_cell[GAMERA_NO_FLUX_COUNT],
    const double upper_cell[GAMERA_NO_FLUX_COUNT], bool use_boris,
    double light_speed, bool *applied);

#ifdef __cplusplus
}
#endif

#endif
