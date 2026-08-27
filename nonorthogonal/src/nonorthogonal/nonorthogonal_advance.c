#include "nonorthogonal_advance.h"

#include <math.h>
#include <stddef.h>

static int same_value(double left, double right) {
  const double scale = fmax(1.0, fmax(fabs(left), fabs(right)));
  return fabs(left - right) <= 8.0e-15 * scale;
}

static int consistent_options(const gamera_no_advance_options *options) {
  if (options == NULL ||
      !same_value(options->stress.gamma, options->update.gamma) ||
      !same_value(options->stress.density_floor,
                  options->update.density_floor) ||
      !same_value(options->stress.pressure_floor,
                  options->update.pressure_floor) ||
      options->stress.use_mhd != options->update.use_mhd ||
      options->stress.use_boris != options->update.use_boris ||
      options->stress.use_background != options->update.use_background ||
      options->emf.use_boris != options->update.use_boris ||
      options->emf.use_background != options->update.use_background ||
      !same_value(options->emf.density_floor,
                  options->update.density_floor) ||
      !isfinite(options->source_time)) {
    return 0;
  }
  return !options->update.use_boris ||
         (same_value(options->stress.light_speed,
                     options->update.light_speed) &&
          same_value(options->emf.light_speed,
                     options->update.light_speed));
}

static int add_cell_source_rates(
    gamera_no_storage *storage, const gamera_no_grid *grid,
    const size_t active_lower[3], const size_t active_upper[3],
    const gamera_no_advance_options *options) {
  if (options->cell_source == NULL && options->indexed_cell_source == NULL) {
    return 0;
  }
  int failed = 0;
#pragma omp parallel for collapse(3) reduction(| : failed) schedule(static)
  for (size_t i = active_lower[0]; i < active_upper[0]; ++i) {
    for (size_t j = active_lower[1]; j < active_upper[1]; ++j) {
      for (size_t k = active_lower[2]; k < active_upper[2]; ++k) {
        const size_t cell =
            gamera_no_index3(storage->cell_extent, i, j, k);
        double rate[GAMERA_NO_FLUX_COUNT] = {0.0};
        if (options->cell_source != NULL &&
            options->cell_source(
                grid->cell[cell].centroid,
                &storage->predicted_conserved[cell * GAMERA_NO_FLUX_COUNT],
                options->source_time, options->source_context, rate) != 0) {
          failed = 1;
          continue;
        }
        if (options->indexed_cell_source != NULL) {
          double indexed_rate[GAMERA_NO_FLUX_COUNT] = {0.0};
          if (options->indexed_cell_source(
                  cell, grid->cell[cell].centroid,
                  &storage->predicted_conserved
                       [cell * GAMERA_NO_FLUX_COUNT],
                  options->source_time, options->indexed_source_context,
                  indexed_rate) != 0) {
            failed = 1;
            continue;
          }
          for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT;
               ++variable) {
            rate[variable] += indexed_rate[variable];
          }
        }
        int finite_rate = 1;
        for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
          finite_rate = finite_rate && isfinite(rate[variable]);
        }
        if (!finite_rate) {
          failed = 1;
          continue;
        }
        for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
          storage->hydro_rate[cell * GAMERA_NO_FLUX_COUNT +
                              (size_t)variable] += rate[variable];
        }
      }
    }
  }
  return failed ? -1 : 0;
}

int gamera_no_advance(
    gamera_no_storage *storage, const gamera_no_grid *grid,
    const size_t active_lower[3], const size_t active_upper[3],
    double predictor_ratio, double dt,
    const gamera_no_advance_options *options,
    const gamera_no_background_field *background) {
  if (storage == NULL || grid == NULL || active_lower == NULL ||
      active_upper == NULL || !isfinite(predictor_ratio) || !isfinite(dt) ||
      dt <= 0.0 || !consistent_options(options) ||
      (options->update.use_background && background == NULL)) {
    return -1;
  }
  storage->emergency_diffusion_face_count = 0U;
  storage->emergency_diffusion_max_interface_speed = 0.0;
  for (int hemisphere = 0; hemisphere < 2; ++hemisphere) {
    storage->inner_wall_clamped_face_count[hemisphere] = 0U;
    storage->inner_wall_positive_mass_max[hemisphere] = 0.0;
    storage->inner_wall_positive_energy_max[hemisphere] = 0.0;
    storage->inner_wall_fluid_momentum_max[hemisphere] = 0.0;
    storage->inner_wall_maxwell_momentum_max[hemisphere] = 0.0;
    storage->inner_wall_density_max[hemisphere] = 0.0;
    storage->inner_wall_pressure_max[hemisphere] = 0.0;
    storage->inner_wall_pressure_gradient_max[hemisphere] = 0.0;
    storage->inner_wall_speed_max[hemisphere] = 0.0;
    storage->inner_wall_residual_magnetic_max[hemisphere] = 0.0;
    storage->inner_wall_residual_magnetic_location[hemisphere] =
        (gamera_no_vec3){{0.0, 0.0, 0.0}};
    storage->inner_wall_lowlat_pressure_max[hemisphere] = 0.0;
    storage->inner_wall_lowlat_pressure_gradient_max[hemisphere] = 0.0;
    storage->inner_wall_lowlat_speed_max[hemisphere] = 0.0;
    storage->inner_wall_lowlat_residual_magnetic_max[hemisphere] = 0.0;
  }
  if (gamera_no_predict_storage(
          storage, grid, predictor_ratio, options->update.gamma,
          options->update.density_floor,
          options->update.pressure_floor) != 0) {
    return -1;
  }

  const double *predicted_face_flux[3] = {
      storage->predicted_face_flux[GAMERA_NO_I],
      storage->predicted_face_flux[GAMERA_NO_J],
      storage->predicted_face_flux[GAMERA_NO_K]};
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    if (gamera_no_sweep_face_fluxes(
            grid, storage, storage->predicted_conserved,
            storage->predicted_cell_magnetic, predicted_face_flux, direction,
            active_lower, active_upper, &options->stress, background) != 0) {
      return -1;
    }
  }
  if (options->fluid_flux_sync != NULL &&
      options->fluid_flux_sync(storage, grid, active_lower, active_upper,
                               options->fluid_flux_context) != 0) {
    return -1;
  }

  if (options->update.use_mhd) {
    gamera_no_emf_options emf_options = options->emf;
    emf_options.dt = dt;
    for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
      if (gamera_no_sweep_edge_emf(
              grid, storage, storage->predicted_conserved,
              predicted_face_flux, direction, active_lower, active_upper,
              &emf_options, background) != 0) {
        return -1;
      }
    }
    if (options->edge_emf_sync != NULL &&
        options->edge_emf_sync(storage, grid, active_lower, active_upper,
                               options->edge_emf_context) != 0) {
      return -1;
    }
  }

  const gamera_no_vec3 *background_force =
      background == NULL ? NULL : background->cell_force;
  if (gamera_no_calculate_stress_rates(
          grid, storage, active_lower, active_upper,
          options->update.use_mhd, background_force) != 0) {
    return -1;
  }
  if (add_cell_source_rates(storage, grid, active_lower, active_upper,
                            options) != 0) {
    return -1;
  }
  const gamera_no_vec3 *background_magnetic =
      background == NULL ? NULL : background->cell_magnetic;
  return gamera_no_apply_active_update(
      storage, grid, active_lower, active_upper, dt, &options->update,
      background_magnetic);
}
