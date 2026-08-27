#include "nonorthogonal_step.h"

#include "nonorthogonal_operators.h"
#include "nonorthogonal_state.h"

#include <math.h>
#include <float.h>
#include <stddef.h>
#include <string.h>

int gamera_no_trap_inner_outward_mass_energy_flux(
    gamera_no_storage *storage, const size_t active_lower[3],
    const size_t active_upper[3], size_t *clamped_face_count) {
  if (storage == NULL || active_lower == NULL || active_upper == NULL ||
      storage->fluid_face_flux[GAMERA_NO_I] == NULL ||
      active_lower[GAMERA_NO_I] >= storage->face_extent[GAMERA_NO_I][0] ||
      active_lower[GAMERA_NO_J] >= active_upper[GAMERA_NO_J] ||
      active_lower[GAMERA_NO_K] >= active_upper[GAMERA_NO_K] ||
      active_upper[GAMERA_NO_J] > storage->face_extent[GAMERA_NO_I][1] ||
      active_upper[GAMERA_NO_K] > storage->face_extent[GAMERA_NO_I][2]) {
    return -1;
  }

  size_t count = 0U;
  const size_t i = active_lower[GAMERA_NO_I];
  for (size_t j = active_lower[GAMERA_NO_J];
       j < active_upper[GAMERA_NO_J]; ++j) {
    for (size_t k = active_lower[GAMERA_NO_K];
         k < active_upper[GAMERA_NO_K]; ++k) {
      const size_t face = gamera_no_index3(
          storage->face_extent[GAMERA_NO_I], i, j, k);
      double *flux =
          &storage->fluid_face_flux[GAMERA_NO_I]
                                   [face * GAMERA_NO_FLUX_COUNT];
      if (flux[GAMERA_NO_FLUX_DENSITY] > 0.0) {
        flux[GAMERA_NO_FLUX_DENSITY] = 0.0;
        if (flux[GAMERA_NO_FLUX_ENERGY] > 0.0) {
          flux[GAMERA_NO_FLUX_ENERGY] = 0.0;
        }
        ++count;
      }
    }
  }
  if (clamped_face_count != NULL) {
    *clamped_face_count = count;
  }
  return 0;
}

int gamera_no_apply_pressure_control(
    gamera_no_storage *storage, const gamera_no_grid *grid,
    const size_t active_lower[3], const size_t active_upper[3], double dt,
    double gamma, double density_floor, double pressure_floor,
    double light_speed, double low_density_threshold,
    double low_density_sound_speed, size_t *low_density_cell_count,
    size_t *high_sound_speed_cell_count, double *maximum_sound_speed,
    double *maximum_pressure) {
  if (storage == NULL || grid == NULL || active_lower == NULL ||
      active_upper == NULL || storage->conserved == NULL || !(dt > 0.0) ||
      !(gamma > 1.0) || !(density_floor > 0.0) ||
      !(pressure_floor > 0.0) || !(light_speed > 0.0) ||
      !(low_density_threshold > 0.0) ||
      !(low_density_sound_speed >= 0.0)) {
    return -1;
  }
  size_t low_count = 0U;
  size_t hot_count = 0U;
  double max_sound = 0.0;
  double max_pressure_local = 0.0;
  int failed = 0;
#pragma omp parallel for collapse(3) reduction(| : failed) \
    reduction(+ : low_count, hot_count) reduction(max : max_sound, max_pressure_local) \
    schedule(static)
  for (size_t i = active_lower[0]; i < active_upper[0]; ++i) {
    for (size_t j = active_lower[1]; j < active_upper[1]; ++j) {
      for (size_t k = active_lower[2]; k < active_upper[2]; ++k) {
        const size_t cell =
            gamera_no_index3(storage->cell_extent, i, j, k);
        double *conserved =
            &storage->conserved[cell * GAMERA_NO_FLUX_COUNT];
        gamera_no_primitive primitive;
        if (gamera_no_conserved_to_primitive(
                conserved, gamma, density_floor, pressure_floor,
                &primitive) != 0) {
          failed = 1;
          continue;
        }
        int changed = 0;
        if (primitive.density <= low_density_threshold) {
          primitive.pressure = fmax(
              pressure_floor,
              primitive.density * low_density_sound_speed *
                  low_density_sound_speed / gamma);
          ++low_count;
          changed = 1;
        }
        const double sound_speed =
            sqrt(gamma * primitive.pressure / primitive.density);
        max_sound = fmax(max_sound, sound_speed);
        max_pressure_local = fmax(max_pressure_local, primitive.pressure);
        if (sound_speed > 1.5 * light_speed) {
          const gamera_no_vec3 point = grid->cell[cell].centroid;
          const double radius_square =
              point.value[0] * point.value[0] +
              point.value[1] * point.value[1] +
              point.value[2] * point.value[2];
          const double cylindrical_square =
              point.value[0] * point.value[0] +
              point.value[1] * point.value[1];
          if (!(radius_square > 0.0) || !(cylindrical_square > DBL_MIN)) {
            failed = 1;
            continue;
          }
          const double radius = sqrt(radius_square);
          const double dipole_l = radius * radius_square / cylindrical_square;
          const double target_pressure =
              primitive.density * (1.5 * light_speed) *
              (1.5 * light_speed) / gamma;
          const double bounce_time = dipole_l / sound_speed;
          primitive.pressure -=
              (dt / bounce_time) * (primitive.pressure - target_pressure);
          ++hot_count;
          changed = 1;
        }
        if (changed &&
            gamera_no_primitive_to_conserved(
                &primitive, gamma, density_floor, pressure_floor,
                conserved) != 0) {
          failed = 1;
        }
      }
    }
  }
  if (low_density_cell_count != NULL) {
    *low_density_cell_count = low_count;
  }
  if (high_sound_speed_cell_count != NULL) {
    *high_sound_speed_cell_count = hot_count;
  }
  if (maximum_sound_speed != NULL) {
    *maximum_sound_speed = max_sound;
  }
  if (maximum_pressure != NULL) {
    *maximum_pressure = max_pressure_local;
  }
  return failed ? -1 : 0;
}

static int matching_extents(const gamera_no_storage *storage,
                            const gamera_no_grid *grid) {
  if (storage == NULL || grid == NULL) {
    return 0;
  }
  for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
    if (storage->cell_extent[axis] != grid->cell_extent[axis]) {
      return 0;
    }
  }
  return 1;
}

double gamera_no_cell_net_flux(const gamera_no_storage *storage, size_t i,
                               size_t j, size_t k) {
  if (storage == NULL || i >= storage->cell_extent[0] ||
      j >= storage->cell_extent[1] || k >= storage->cell_extent[2]) {
    return NAN;
  }
  const double lower_i =
      storage->face_flux[GAMERA_NO_I][gamera_no_index3(
          storage->face_extent[GAMERA_NO_I], i, j, k)];
  const double upper_i =
      storage->face_flux[GAMERA_NO_I][gamera_no_index3(
          storage->face_extent[GAMERA_NO_I], i + 1, j, k)];
  const double lower_j =
      storage->face_flux[GAMERA_NO_J][gamera_no_index3(
          storage->face_extent[GAMERA_NO_J], i, j, k)];
  const double upper_j =
      storage->face_flux[GAMERA_NO_J][gamera_no_index3(
          storage->face_extent[GAMERA_NO_J], i, j + 1, k)];
  const double lower_k =
      storage->face_flux[GAMERA_NO_K][gamera_no_index3(
          storage->face_extent[GAMERA_NO_K], i, j, k)];
  const double upper_k =
      storage->face_flux[GAMERA_NO_K][gamera_no_index3(
          storage->face_extent[GAMERA_NO_K], i, j, k + 1)];
  return upper_i - lower_i + upper_j - lower_j + upper_k - lower_k;
}

static int valid_extent_range(const size_t extent[3], const size_t lower[3],
                              const size_t upper[3]) {
  if (extent == NULL || lower == NULL || upper == NULL) {
    return 0;
  }
  for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
    if (lower[axis] >= upper[axis] || upper[axis] > extent[axis]) {
      return 0;
    }
  }
  return 1;
}

static int valid_active_range(const gamera_no_storage *storage,
                              const size_t lower[3],
                              const size_t upper[3]) {
  return storage != NULL &&
         valid_extent_range(storage->cell_extent, lower, upper);
}

int gamera_no_advance_ct_active(gamera_no_storage *storage,
                                const size_t active_lower[3],
                                const size_t active_upper[3], double dt) {
  if (!valid_active_range(storage, active_lower, active_upper) ||
      !isfinite(dt)) {
    return -1;
  }

#pragma omp parallel for collapse(3) schedule(static)
  for (size_t i = active_lower[0]; i <= active_upper[0]; ++i) {
    for (size_t j = active_lower[1]; j < active_upper[1]; ++j) {
      for (size_t k = active_lower[2]; k < active_upper[2]; ++k) {
        const double curl =
            storage->edge_emf[GAMERA_NO_K][gamera_no_index3(
                storage->edge_extent[GAMERA_NO_K], i, j + 1, k)] -
            storage->edge_emf[GAMERA_NO_K][gamera_no_index3(
                storage->edge_extent[GAMERA_NO_K], i, j, k)] -
            storage->edge_emf[GAMERA_NO_J][gamera_no_index3(
                storage->edge_extent[GAMERA_NO_J], i, j, k + 1)] +
            storage->edge_emf[GAMERA_NO_J][gamera_no_index3(
                storage->edge_extent[GAMERA_NO_J], i, j, k)];
        storage->face_flux[GAMERA_NO_I][gamera_no_index3(
            storage->face_extent[GAMERA_NO_I], i, j, k)] -= dt * curl;
      }
    }
  }

#pragma omp parallel for collapse(3) schedule(static)
  for (size_t i = active_lower[0]; i < active_upper[0]; ++i) {
    for (size_t j = active_lower[1]; j <= active_upper[1]; ++j) {
      for (size_t k = active_lower[2]; k < active_upper[2]; ++k) {
        const double curl =
            storage->edge_emf[GAMERA_NO_I][gamera_no_index3(
                storage->edge_extent[GAMERA_NO_I], i, j, k + 1)] -
            storage->edge_emf[GAMERA_NO_I][gamera_no_index3(
                storage->edge_extent[GAMERA_NO_I], i, j, k)] -
            storage->edge_emf[GAMERA_NO_K][gamera_no_index3(
                storage->edge_extent[GAMERA_NO_K], i + 1, j, k)] +
            storage->edge_emf[GAMERA_NO_K][gamera_no_index3(
                storage->edge_extent[GAMERA_NO_K], i, j, k)];
        storage->face_flux[GAMERA_NO_J][gamera_no_index3(
            storage->face_extent[GAMERA_NO_J], i, j, k)] -= dt * curl;
      }
    }
  }

#pragma omp parallel for collapse(3) schedule(static)
  for (size_t i = active_lower[0]; i < active_upper[0]; ++i) {
    for (size_t j = active_lower[1]; j < active_upper[1]; ++j) {
      for (size_t k = active_lower[2]; k <= active_upper[2]; ++k) {
        const double curl =
            storage->edge_emf[GAMERA_NO_J][gamera_no_index3(
                storage->edge_extent[GAMERA_NO_J], i + 1, j, k)] -
            storage->edge_emf[GAMERA_NO_J][gamera_no_index3(
                storage->edge_extent[GAMERA_NO_J], i, j, k)] -
            storage->edge_emf[GAMERA_NO_I][gamera_no_index3(
                storage->edge_extent[GAMERA_NO_I], i, j + 1, k)] +
            storage->edge_emf[GAMERA_NO_I][gamera_no_index3(
                storage->edge_extent[GAMERA_NO_I], i, j, k)];
        storage->face_flux[GAMERA_NO_K][gamera_no_index3(
            storage->face_extent[GAMERA_NO_K], i, j, k)] -= dt * curl;
      }
    }
  }
  return 0;
}

int gamera_no_advance_ct(gamera_no_storage *storage, double dt) {
  if (storage == NULL) {
    return -1;
  }
  const size_t lower[3] = {0, 0, 0};
  const size_t upper[3] = {storage->cell_extent[0], storage->cell_extent[1],
                           storage->cell_extent[2]};
  return gamera_no_advance_ct_active(storage, lower, upper, dt);
}

int gamera_no_recover_magnetic_field(
    const gamera_no_grid *grid, const double *const face_flux[3],
    gamera_no_vec3 *cell_magnetic) {
  if (grid == NULL || face_flux == NULL || cell_magnetic == NULL) {
    return -1;
  }
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    if (face_flux[direction] == NULL) {
      return -1;
    }
  }

  int failed = 0;
#pragma omp parallel for collapse(3) reduction(| : failed) schedule(static)
  for (size_t i = 0; i < grid->cell_extent[0]; ++i) {
    for (size_t j = 0; j < grid->cell_extent[1]; ++j) {
      for (size_t k = 0; k < grid->cell_extent[2]; ++k) {
        double local_flux[3][2];
        local_flux[GAMERA_NO_I][GAMERA_NO_LOWER] =
            face_flux[GAMERA_NO_I][gamera_no_index3(
                grid->face[GAMERA_NO_I].extent, i, j, k)];
        local_flux[GAMERA_NO_I][GAMERA_NO_UPPER] =
            face_flux[GAMERA_NO_I][gamera_no_index3(
                grid->face[GAMERA_NO_I].extent, i + 1, j, k)];
        local_flux[GAMERA_NO_J][GAMERA_NO_LOWER] =
            face_flux[GAMERA_NO_J][gamera_no_index3(
                grid->face[GAMERA_NO_J].extent, i, j, k)];
        local_flux[GAMERA_NO_J][GAMERA_NO_UPPER] =
            face_flux[GAMERA_NO_J][gamera_no_index3(
                grid->face[GAMERA_NO_J].extent, i, j + 1, k)];
        local_flux[GAMERA_NO_K][GAMERA_NO_LOWER] =
            face_flux[GAMERA_NO_K][gamera_no_index3(
                grid->face[GAMERA_NO_K].extent, i, j, k)];
        local_flux[GAMERA_NO_K][GAMERA_NO_UPPER] =
            face_flux[GAMERA_NO_K][gamera_no_index3(
                grid->face[GAMERA_NO_K].extent, i, j, k + 1)];
        const size_t cell_index =
            gamera_no_index3(grid->cell_extent, i, j, k);
        if (gamera_no_flux_to_cell_field(&grid->cell[cell_index], local_flux,
                                         &cell_magnetic[cell_index], NULL) !=
            0) {
          failed = 1;
        }
      }
    }
  }
  return failed ? -1 : 0;
}

int gamera_no_save_current_as_old(gamera_no_storage *storage) {
  if (storage == NULL) {
    return -1;
  }
  const size_t cell_count = gamera_no_element_count3(storage->cell_extent);
  memcpy(storage->old_conserved, storage->conserved,
         cell_count * GAMERA_NO_FLUX_COUNT * sizeof(*storage->conserved));
  memcpy(storage->old_cell_magnetic, storage->cell_magnetic,
         cell_count * sizeof(*storage->cell_magnetic));
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t face_count =
        gamera_no_element_count3(storage->face_extent[direction]);
    memcpy(storage->old_face_flux[direction], storage->face_flux[direction],
           face_count * sizeof(*storage->face_flux[direction]));
  }
  return 0;
}

int gamera_no_predict_storage(gamera_no_storage *storage,
                              const gamera_no_grid *grid, double ratio,
                              double gamma, double density_floor,
                              double pressure_floor) {
  if (!matching_extents(storage, grid) || !isfinite(ratio)) {
    return -1;
  }
  const size_t cell_count = gamera_no_element_count3(storage->cell_extent);
  int failed = 0;
#pragma omp parallel for reduction(| : failed) schedule(static)
  for (size_t cell = 0; cell < cell_count; ++cell) {
    const size_t offset = cell * GAMERA_NO_FLUX_COUNT;
    if (gamera_no_predict_cell(&storage->old_conserved[offset],
                               &storage->conserved[offset], ratio, gamma,
                               density_floor, pressure_floor,
                               &storage->predicted_conserved[offset]) != 0) {
      failed = 1;
    }
  }
  if (failed) {
    return -1;
  }
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t face_count =
        gamera_no_element_count3(storage->face_extent[direction]);
#pragma omp parallel for schedule(static)
    for (size_t face = 0; face < face_count; ++face) {
      storage->predicted_face_flux[direction][face] =
          storage->face_flux[direction][face] +
          ratio * (storage->face_flux[direction][face] -
                   storage->old_face_flux[direction][face]);
    }
  }
  const double *predicted_flux[3] = {
      storage->predicted_face_flux[GAMERA_NO_I],
      storage->predicted_face_flux[GAMERA_NO_J],
      storage->predicted_face_flux[GAMERA_NO_K]};
  return gamera_no_recover_magnetic_field(
      grid, predicted_flux, storage->predicted_cell_magnetic);
}

static int valid_update_options(const gamera_no_update_options *options) {
  if (options == NULL || !isfinite(options->gamma) || options->gamma <= 1.0 ||
      !isfinite(options->density_floor) || options->density_floor <= 0.0 ||
      !isfinite(options->pressure_floor) || options->pressure_floor <= 0.0) {
    return 0;
  }
  return !options->use_boris ||
         (isfinite(options->light_speed) && options->light_speed > DBL_MIN);
}

static gamera_no_vec3 add_vectors(gamera_no_vec3 left,
                                  gamera_no_vec3 right) {
  gamera_no_vec3 result;
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    result.value[component] =
        left.value[component] + right.value[component];
  }
  return result;
}

static double vector_norm(gamera_no_vec3 vector) {
  double square = 0.0;
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    square += vector.value[component] * vector.value[component];
  }
  return sqrt(square);
}

int gamera_no_local_timestep(
    const gamera_no_grid *grid, const double *conserved,
    const gamera_no_vec3 *cell_magnetic, const size_t active_lower[3],
    const size_t active_upper[3],
    const gamera_no_timestep_options *options,
    const gamera_no_vec3 *background_cell_magnetic, double *local_dt) {
  if (grid == NULL || conserved == NULL || options == NULL ||
      local_dt == NULL || !isfinite(options->gamma) ||
      options->gamma <= 1.0 || !isfinite(options->density_floor) ||
      options->density_floor <= 0.0 || !isfinite(options->pressure_floor) ||
      options->pressure_floor <= 0.0 || !isfinite(options->cfl) ||
      options->cfl <= 0.0 ||
      (options->use_mhd && cell_magnetic == NULL) ||
      (options->use_background && background_cell_magnetic == NULL) ||
      (options->use_boris &&
       (!isfinite(options->light_speed) ||
        options->light_speed <= DBL_MIN))) {
    return -1;
  }
  if (!valid_extent_range(grid->cell_extent, active_lower, active_upper)) {
    return -1;
  }

  double minimum = DBL_MAX;
  int failed = 0;
#pragma omp parallel for collapse(3) reduction(min : minimum) \
    reduction(| : failed) schedule(static)
  for (size_t i = active_lower[0]; i < active_upper[0]; ++i) {
    for (size_t j = active_lower[1]; j < active_upper[1]; ++j) {
      for (size_t k = active_lower[2]; k < active_upper[2]; ++k) {
        const size_t cell = gamera_no_index3(grid->cell_extent, i, j, k);
        gamera_no_primitive primitive;
        if (gamera_no_conserved_to_primitive(
                &conserved[cell * GAMERA_NO_FLUX_COUNT], options->gamma,
                options->density_floor, options->pressure_floor,
                &primitive) != 0) {
          failed = 1;
          continue;
        }
        const double sound_speed =
            sqrt(options->gamma * primitive.pressure / primitive.density);
        double alfven_speed = 0.0;
        if (options->use_mhd) {
          gamera_no_vec3 magnetic = cell_magnetic[cell];
          if (options->use_background) {
            magnetic =
                add_vectors(magnetic, background_cell_magnetic[cell]);
          }
          alfven_speed = vector_norm(magnetic) / sqrt(primitive.density);
          if (options->use_boris) {
            alfven_speed =
                options->light_speed * alfven_speed /
                hypot(options->light_speed, alfven_speed);
          }
        }
        const double characteristic_speed =
            vector_norm(primitive.velocity) +
            hypot(sound_speed, alfven_speed);
        double length = grid->cell[cell].cfl_length[0];
        for (int direction = 1; direction < GAMERA_NO_DIM; ++direction) {
          length = fmin(length, grid->cell[cell].cfl_length[direction]);
        }
        const double candidate = options->cfl * length / characteristic_speed;
        if (!isfinite(candidate) || candidate <= 0.0) {
          failed = 1;
          continue;
        }
        minimum = fmin(minimum, candidate);
      }
    }
  }
  if (failed) {
    return -1;
  }
  *local_dt = minimum;
  return 0;
}

int gamera_no_apply_active_update(
    gamera_no_storage *storage, const gamera_no_grid *grid,
    const size_t active_lower[3], const size_t active_upper[3], double dt,
    const gamera_no_update_options *options,
    const gamera_no_vec3 *background_cell_magnetic) {
  if (!matching_extents(storage, grid) ||
      !valid_active_range(storage, active_lower, active_upper) ||
      !isfinite(dt) || !valid_update_options(options) ||
      (options->use_background && background_cell_magnetic == NULL)) {
    return -1;
  }
  if (gamera_no_save_current_as_old(storage) != 0) {
    return -1;
  }

#pragma omp parallel for collapse(3) schedule(static)
  for (size_t i = active_lower[0]; i < active_upper[0]; ++i) {
    for (size_t j = active_lower[1]; j < active_upper[1]; ++j) {
      for (size_t k = active_lower[2]; k < active_upper[2]; ++k) {
        const size_t cell =
            gamera_no_index3(storage->cell_extent, i, j, k);
        gamera_no_apply_reynolds(
            &storage->conserved[cell * GAMERA_NO_FLUX_COUNT],
            &storage->hydro_rate[cell * GAMERA_NO_FLUX_COUNT], dt);
      }
    }
  }
  if (!options->use_mhd) {
    return 0;
  }

  if (gamera_no_advance_ct_active(storage, active_lower, active_upper, dt) !=
      0) {
    return -1;
  }
  const double *current_flux[3] = {storage->face_flux[GAMERA_NO_I],
                                   storage->face_flux[GAMERA_NO_J],
                                   storage->face_flux[GAMERA_NO_K]};
  if (gamera_no_recover_magnetic_field(grid, current_flux,
                                       storage->cell_magnetic) != 0) {
    return -1;
  }

  int failed = 0;
#pragma omp parallel for collapse(3) reduction(| : failed) schedule(static)
  for (size_t i = active_lower[0]; i < active_upper[0]; ++i) {
    for (size_t j = active_lower[1]; j < active_upper[1]; ++j) {
      for (size_t k = active_lower[2]; k < active_upper[2]; ++k) {
        const size_t cell =
            gamera_no_index3(storage->cell_extent, i, j, k);
        double *updated =
            &storage->conserved[cell * GAMERA_NO_FLUX_COUNT];
        if (options->use_boris) {
          gamera_no_vec3 new_magnetic = storage->cell_magnetic[cell];
          gamera_no_vec3 old_magnetic = storage->old_cell_magnetic[cell];
          if (options->use_background) {
            new_magnetic =
                add_vectors(new_magnetic, background_cell_magnetic[cell]);
            old_magnetic =
                add_vectors(old_magnetic, background_cell_magnetic[cell]);
          }
          if (gamera_no_apply_boris(
                  updated,
                  &storage->old_conserved[cell * GAMERA_NO_FLUX_COUNT],
                  new_magnetic, old_magnetic,
                  &storage->hydro_rate[cell * GAMERA_NO_FLUX_COUNT],
                  storage->maxwell_rate[cell], options->light_speed, dt,
                  options->gamma, options->density_floor,
                  options->pressure_floor) != 0) {
            failed = 1;
          }
        } else if (gamera_no_apply_maxwell(
                       updated, storage->maxwell_rate[cell], dt,
                       options->gamma, options->density_floor,
                       options->pressure_floor) != 0) {
          failed = 1;
        }
      }
    }
  }
  return failed ? -1 : 0;
}
