#ifndef GAMERA_NONORTHOGONAL_STORAGE_H
#define GAMERA_NONORTHOGONAL_STORAGE_H

#include "nonorthogonal_flux.h"

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  size_t cell_extent[3];
  size_t face_extent[3][3];
  size_t edge_extent[3][3];
  double *conserved;
  double *old_conserved;
  double *predicted_conserved;
  gamera_no_vec3 *cell_magnetic;
  gamera_no_vec3 *old_cell_magnetic;
  gamera_no_vec3 *predicted_cell_magnetic;
  double *face_flux[3];
  double *old_face_flux[3];
  double *predicted_face_flux[3];
  double *edge_emf[3];
  double *fluid_face_flux[3];
  gamera_no_vec3 *maxwell_face_flux[3];
  double *hydro_rate;
  gamera_no_vec3 *maxwell_rate;
  size_t emergency_diffusion_face_count;
  double emergency_diffusion_max_interface_speed;
  size_t inner_wall_clamped_face_count[2];
  double inner_wall_positive_mass_max[2];
  double inner_wall_positive_energy_max[2];
  double inner_wall_fluid_momentum_max[2];
  double inner_wall_maxwell_momentum_max[2];
  double inner_wall_density_max[2];
  double inner_wall_pressure_max[2];
  double inner_wall_pressure_gradient_max[2];
  double inner_wall_speed_max[2];
  double inner_wall_residual_magnetic_max[2];
  gamera_no_vec3 inner_wall_residual_magnetic_location[2];
  double inner_wall_lowlat_pressure_max[2];
  double inner_wall_lowlat_pressure_gradient_max[2];
  double inner_wall_lowlat_speed_max[2];
  double inner_wall_lowlat_residual_magnetic_max[2];
} gamera_no_storage;

int gamera_no_storage_create(const size_t cell_extent[3],
                             gamera_no_storage *storage);
void gamera_no_storage_destroy(gamera_no_storage *storage);

#ifdef __cplusplus
}
#endif

#endif
