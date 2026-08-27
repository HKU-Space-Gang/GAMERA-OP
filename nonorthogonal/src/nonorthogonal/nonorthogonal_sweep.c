#include "nonorthogonal_sweep.h"

#include "nonorthogonal_flux.h"
#include "nonorthogonal_operators.h"
#include "nonorthogonal_reconstruction.h"
#include "nonorthogonal_state.h"

#include <float.h>
#include <math.h>
#include <stddef.h>

static int matching_extents(const gamera_no_grid *grid,
                            const gamera_no_storage *storage) {
  if (grid == NULL || storage == NULL) {
    return 0;
  }
  for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
    if (grid->cell_extent[axis] != storage->cell_extent[axis]) {
      return 0;
    }
  }
  return 1;
}

static int valid_active_range(const gamera_no_grid *grid,
                              const size_t lower[3],
                              const size_t upper[3]) {
  if (grid == NULL || lower == NULL || upper == NULL) {
    return 0;
  }
  for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
    if (lower[axis] >= upper[axis] || upper[axis] > grid->cell_extent[axis]) {
      return 0;
    }
  }
  return 1;
}

static int valid_options(const gamera_no_sweep_options *options) {
  if (options == NULL || !isfinite(options->gamma) || options->gamma <= 1.0 ||
      !isfinite(options->density_floor) || options->density_floor <= 0.0 ||
      !isfinite(options->pressure_floor) || options->pressure_floor <= 0.0 ||
      !isfinite(options->pdm_coefficient) ||
      options->pdm_coefficient < 0.0 ||
      !isfinite(options->hydro_hogs_coefficient) ||
      options->hydro_hogs_coefficient < 0.0 ||
      !isfinite(options->magnetic_hogs_coefficient) ||
      options->magnetic_hogs_coefficient < 0.0) {
    return 0;
  }
  return !options->use_boris ||
         (isfinite(options->light_speed) && options->light_speed > DBL_MIN);
}

static double primitive_component(gamera_no_primitive primitive,
                                  int variable) {
  if (variable == GAMERA_NO_FLUX_DENSITY) {
    return primitive.density;
  }
  if (variable == GAMERA_NO_FLUX_ENERGY) {
    return primitive.pressure;
  }
  return primitive.velocity.value[variable - GAMERA_NO_FLUX_MOMENTUM_X];
}

static void set_primitive_component(gamera_no_primitive *primitive,
                                    int variable, double value) {
  if (variable == GAMERA_NO_FLUX_DENSITY) {
    primitive->density = value;
  } else if (variable == GAMERA_NO_FLUX_ENERGY) {
    primitive->pressure = value;
  } else {
    primitive->velocity.value[variable - GAMERA_NO_FLUX_MOMENTUM_X] = value;
  }
}

static int reconstruct_face_states(
    const gamera_no_grid *grid, const double *conserved,
    const gamera_no_vec3 *cell_magnetic, int direction,
    const size_t face_coordinate[3], const gamera_no_sweep_options *options,
    gamera_no_primitive primitive[2], gamera_no_vec3 magnetic[2]) {
  double volume[GAMERA_NO_RECON_STENCIL];
  gamera_no_primitive primitive_stencil[GAMERA_NO_RECON_STENCIL];
  gamera_no_vec3 magnetic_stencil[GAMERA_NO_RECON_STENCIL];

  for (size_t n = 0; n < GAMERA_NO_RECON_STENCIL; ++n) {
    size_t cell_coordinate[3] = {face_coordinate[0], face_coordinate[1],
                                 face_coordinate[2]};
    cell_coordinate[direction] = face_coordinate[direction] + n - 4U;
    const size_t cell = gamera_no_index3(
        grid->cell_extent, cell_coordinate[0], cell_coordinate[1],
        cell_coordinate[2]);
    volume[n] = grid->cell[cell].volume;
    if (gamera_no_conserved_to_primitive(
            &conserved[cell * GAMERA_NO_FLUX_COUNT], options->gamma,
            options->density_floor, options->pressure_floor,
            &primitive_stencil[n]) != 0) {
      return -1;
    }
    if (options->use_mhd) {
      magnetic_stencil[n] = cell_magnetic[cell];
    }
  }

  for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
    double stencil[GAMERA_NO_RECON_STENCIL];
    for (int n = 0; n < GAMERA_NO_RECON_STENCIL; ++n) {
      stencil[n] = primitive_component(primitive_stencil[n], variable);
    }
    double left = 0.0;
    double right = 0.0;
    if (gamera_no_reconstruct_up7_pdm(volume, stencil,
                                      options->pdm_coefficient, &left,
                                      &right) != 0) {
      return -1;
    }
    set_primitive_component(&primitive[GAMERA_NO_LOWER], variable, left);
    set_primitive_component(&primitive[GAMERA_NO_UPPER], variable, right);
  }

  if (options->use_mhd) {
    for (int component = 0; component < GAMERA_NO_DIM; ++component) {
      double stencil[GAMERA_NO_RECON_STENCIL];
      for (int n = 0; n < GAMERA_NO_RECON_STENCIL; ++n) {
        stencil[n] = magnetic_stencil[n].value[component];
      }
      if (gamera_no_reconstruct_up7_pdm(
              volume, stencil, options->pdm_coefficient,
              &magnetic[GAMERA_NO_LOWER].value[component],
              &magnetic[GAMERA_NO_UPPER].value[component]) != 0) {
        return -1;
      }
    }
  } else {
    magnetic[GAMERA_NO_LOWER] = (gamera_no_vec3){{0.0, 0.0, 0.0}};
    magnetic[GAMERA_NO_UPPER] = (gamera_no_vec3){{0.0, 0.0, 0.0}};
  }
  return 0;
}

static int valid_background(const gamera_no_sweep_options *options,
                            const gamera_no_background_field *background,
                            int direction) {
  if (!options->use_background) {
    return 1;
  }
  return background != NULL && background->face_magnetic[direction] != NULL &&
         background->face_flux[direction] != NULL;
}

static double vector_dot(gamera_no_vec3 left, gamera_no_vec3 right) {
  double result = 0.0;
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    result += left.value[component] * right.value[component];
  }
  return result;
}

static int valid_emf_options(const gamera_no_emf_options *options) {
  if (options == NULL || !isfinite(options->density_floor) ||
      options->density_floor <= 0.0 ||
      !isfinite(options->pdm_coefficient) ||
      options->pdm_coefficient < 0.0 ||
      !isfinite(options->diffusion_coefficient) || !isfinite(options->cfl) ||
      options->cfl < 0.0 || !isfinite(options->dt) ||
      options->dt <= DBL_MIN) {
    return 0;
  }
  return !options->use_boris ||
         (isfinite(options->light_speed) && options->light_speed > DBL_MIN);
}

static int face_velocity(const gamera_no_grid *grid, const double *conserved,
                         int face_direction,
                         const size_t face_coordinate[3],
                         double density_floor, gamera_no_vec3 *velocity) {
  double volume[GAMERA_NO_RECON_STENCIL];
  double weighted_velocity[GAMERA_NO_DIM][GAMERA_NO_RECON_STENCIL];
  for (size_t n = 0; n < GAMERA_NO_RECON_STENCIL; ++n) {
    size_t cell_coordinate[3] = {face_coordinate[0], face_coordinate[1],
                                 face_coordinate[2]};
    cell_coordinate[face_direction] =
        face_coordinate[face_direction] + n - 4U;
    const size_t cell = gamera_no_index3(
        grid->cell_extent, cell_coordinate[0], cell_coordinate[1],
        cell_coordinate[2]);
    volume[n] = grid->cell[cell].volume;
    const double density =
        fmax(conserved[cell * GAMERA_NO_FLUX_COUNT +
                       GAMERA_NO_FLUX_DENSITY],
             density_floor);
    for (int component = 0; component < GAMERA_NO_DIM; ++component) {
      weighted_velocity[component][n] =
          volume[n] *
          conserved[cell * GAMERA_NO_FLUX_COUNT +
                    (size_t)(GAMERA_NO_FLUX_MOMENTUM_X + component)] /
          density;
    }
  }
  const double interface_volume = gamera_no_central8(volume);
  if (!isfinite(interface_volume) || fabs(interface_volume) <= DBL_MIN) {
    return -1;
  }
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    velocity->value[component] =
        gamera_no_central8(weighted_velocity[component]) / interface_volume;
  }
  return 0;
}

static int edge_velocity(const gamera_no_grid *grid, const double *conserved,
                         int first_transverse, int second_transverse,
                         const size_t edge_coordinate[3],
                         double density_floor, gamera_no_vec3 *velocity) {
  double area[GAMERA_NO_RECON_STENCIL];
  gamera_no_vec3 face_value[GAMERA_NO_RECON_STENCIL];
  for (size_t n = 0; n < GAMERA_NO_RECON_STENCIL; ++n) {
    size_t face_coordinate[3] = {edge_coordinate[0], edge_coordinate[1],
                                 edge_coordinate[2]};
    face_coordinate[second_transverse] =
        edge_coordinate[second_transverse] + n - 4U;
    const size_t face = gamera_no_index3(
        grid->face[first_transverse].extent, face_coordinate[0],
        face_coordinate[1], face_coordinate[2]);
    area[n] = grid->face[first_transverse].value[face].area;
    if (face_velocity(grid, conserved, first_transverse, face_coordinate,
                      density_floor, &face_value[n]) != 0) {
      return -1;
    }
  }
  const double edge_area = gamera_no_central8(area);
  if (!isfinite(edge_area) || fabs(edge_area) <= DBL_MIN) {
    return -1;
  }
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    double weighted[GAMERA_NO_RECON_STENCIL];
    for (int n = 0; n < GAMERA_NO_RECON_STENCIL; ++n) {
      weighted[n] = area[n] * face_value[n].value[component];
    }
    velocity->value[component] = gamera_no_central8(weighted) / edge_area;
  }
  return 0;
}

static int edge_density(const gamera_no_grid *grid, const double *conserved,
                        int first_transverse, int second_transverse,
                        const size_t edge_coordinate[3], double *density) {
  *density = 0.0;
  for (size_t first_side = 0; first_side < 2U; ++first_side) {
    for (size_t second_side = 0; second_side < 2U; ++second_side) {
      size_t cell_coordinate[3] = {edge_coordinate[0], edge_coordinate[1],
                                   edge_coordinate[2]};
      cell_coordinate[first_transverse] -= first_side;
      cell_coordinate[second_transverse] -= second_side;
      const size_t cell = gamera_no_index3(
          grid->cell_extent, cell_coordinate[0], cell_coordinate[1],
          cell_coordinate[2]);
      *density +=
          0.25 * conserved[cell * GAMERA_NO_FLUX_COUNT +
                           GAMERA_NO_FLUX_DENSITY];
    }
  }
  return isfinite(*density) && *density > DBL_MIN ? 0 : -1;
}

static int transverse_face_field(
    const gamera_no_grid *grid, const double *face_magnetic_flux,
    int face_direction, int sweep_direction,
    const size_t edge_coordinate[3], double pdm_coefficient, double *left,
    double *right, double transverse_normal[2]) {
  double area[GAMERA_NO_RECON_STENCIL];
  double field[GAMERA_NO_RECON_STENCIL];
  gamera_no_vec3 normal[GAMERA_NO_RECON_STENCIL];
  for (size_t n = 0; n < GAMERA_NO_RECON_STENCIL; ++n) {
    size_t face_coordinate[3] = {edge_coordinate[0], edge_coordinate[1],
                                 edge_coordinate[2]};
    face_coordinate[sweep_direction] =
        edge_coordinate[sweep_direction] + n - 4U;
    const size_t face = gamera_no_index3(
        grid->face[face_direction].extent, face_coordinate[0],
        face_coordinate[1], face_coordinate[2]);
    const gamera_no_face_geometry *geometry =
        &grid->face[face_direction].value[face];
    area[n] = geometry->area;
    field[n] = face_magnetic_flux[face] / geometry->area;
    normal[n] = geometry->normal;
  }
  const size_t edge = gamera_no_index3(
      grid->edge[3 - face_direction - sweep_direction].extent,
      edge_coordinate[0], edge_coordinate[1], edge_coordinate[2]);
  const gamera_no_edge_geometry *edge_geometry =
      &grid->edge[3 - face_direction - sweep_direction].value[edge];
  if (gamera_no_reconstruct_up7_pdm(area, field, pdm_coefficient, left,
                                    right) != 0 ||
      gamera_no_interpolate_face_normal_to_edge(
          area, normal, edge_geometry, transverse_normal) != 0) {
    return -1;
  }
  return 0;
}

int gamera_no_sweep_face_fluxes(
    const gamera_no_grid *grid, gamera_no_storage *storage,
    const double *conserved, const gamera_no_vec3 *cell_magnetic,
    const double *const face_magnetic_flux[3], int direction,
    const size_t active_lower[3], const size_t active_upper[3],
    const gamera_no_sweep_options *options,
    const gamera_no_background_field *background) {
  if (!matching_extents(grid, storage) || conserved == NULL ||
      !valid_active_range(grid, active_lower, active_upper) ||
      !valid_options(options) || direction < 0 || direction >= GAMERA_NO_DIM ||
      grid->cell_extent[direction] < 4U ||
      active_lower[direction] < 4U ||
      active_upper[direction] > grid->cell_extent[direction] - 4U ||
      (options->use_mhd &&
       (cell_magnetic == NULL || face_magnetic_flux == NULL ||
        face_magnetic_flux[direction] == NULL)) ||
      !valid_background(options, background, direction)) {
    return -1;
  }

  size_t face_upper[3] = {active_upper[0], active_upper[1], active_upper[2]};
  ++face_upper[direction];
  int failed = 0;
  size_t nuclear_count = 0U;
  double nuclear_maximum_speed = 0.0;
#pragma omp parallel for collapse(3) reduction(| : failed) \
    reduction(+ : nuclear_count) reduction(max : nuclear_maximum_speed) \
    schedule(static)
  for (size_t i = active_lower[0]; i < face_upper[0]; ++i) {
    for (size_t j = active_lower[1]; j < face_upper[1]; ++j) {
      for (size_t k = active_lower[2]; k < face_upper[2]; ++k) {
        const size_t coordinate[3] = {i, j, k};
        const size_t face = gamera_no_index3(grid->face[direction].extent, i,
                                             j, k);
        const gamera_no_face_geometry *geometry =
            &grid->face[direction].value[face];
        gamera_no_primitive primitive[2];
        gamera_no_vec3 magnetic[2];
        if (reconstruct_face_states(grid, conserved, cell_magnetic, direction,
                                    coordinate, options, primitive,
                                    magnetic) != 0) {
          failed = 1;
          continue;
        }

        gamera_no_fluid_flux fluid;
        if (gamera_no_kinetic_fluid_flux(primitive, options->gamma, geometry,
                                         &fluid) != 0) {
          failed = 1;
          continue;
        }
        gamera_no_maxwell_flux maxwell = {
            {{0.0, 0.0, 0.0}}, 0.0, 0.0};
        if (options->use_mhd) {
          const double normal_field =
              face_magnetic_flux[direction][face] / geometry->area;
          gamera_no_vec3 background_field = {{0.0, 0.0, 0.0}};
          double background_normal_field = 0.0;
          if (options->use_background) {
            background_field = background->face_magnetic[direction][face];
            background_normal_field =
                background->face_flux[direction][face] / geometry->area;
          }
          if (gamera_no_kinetic_maxwell_flux(
                  primitive, magnetic, normal_field, geometry,
                  options->use_background, background_field,
                  background_normal_field, options->use_boris,
                  options->light_speed, fluid.normal_velocity, &maxwell) != 0) {
            failed = 1;
            continue;
          }
          if (options->use_hogs &&
              gamera_no_apply_hogs(
                  &fluid, &maxwell, options->hydro_hogs_coefficient,
                  options->magnetic_hogs_coefficient, options->use_boris,
                  options->light_speed) != 0) {
            failed = 1;
            continue;
          }
          if (options->use_hogs && options->use_boris) {
            size_t lower_coordinate[3] = {i, j, k};
            size_t upper_coordinate[3] = {i, j, k};
            --lower_coordinate[direction];
            const size_t lower_cell = gamera_no_index3(
                grid->cell_extent, lower_coordinate[0], lower_coordinate[1],
                lower_coordinate[2]);
            const size_t upper_cell = gamera_no_index3(
                grid->cell_extent, upper_coordinate[0], upper_coordinate[1],
                upper_coordinate[2]);
            bool applied = false;
            if (gamera_no_apply_nuclear_hogs(
                    &fluid, primitive,
                    &conserved[lower_cell * GAMERA_NO_FLUX_COUNT],
                    &conserved[upper_cell * GAMERA_NO_FLUX_COUNT], true,
                    options->light_speed, &applied) != 0) {
              failed = 1;
              continue;
            }
            if (applied) {
              double speed_square_sum = 0.0;
              for (int side = 0; side < 2; ++side) {
                for (int component = 0; component < GAMERA_NO_DIM;
                     ++component) {
                  const double value =
                      primitive[side].velocity.value[component];
                  speed_square_sum += value * value;
                }
              }
              nuclear_maximum_speed =
                  fmax(nuclear_maximum_speed, sqrt(0.5 * speed_square_sum));
              ++nuclear_count;
            }
          }
        }

        for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
          storage->fluid_face_flux[direction]
                                  [face * GAMERA_NO_FLUX_COUNT +
                                   (size_t)variable] =
              geometry->area * fluid.conserved[variable];
        }
        for (int component = 0; component < GAMERA_NO_DIM; ++component) {
          storage->maxwell_face_flux[direction][face].value[component] =
              geometry->area * maxwell.momentum.value[component];
        }
      }
    }
  }
  storage->nuclear_hogs_face_count += nuclear_count;
  storage->nuclear_hogs_max_interface_speed =
      fmax(storage->nuclear_hogs_max_interface_speed,
           nuclear_maximum_speed);
  return failed ? -1 : 0;
}

int gamera_no_calculate_stress_rates(
    const gamera_no_grid *grid, gamera_no_storage *storage,
    const size_t active_lower[3], const size_t active_upper[3], bool use_mhd,
    const gamera_no_vec3 *background_cell_force) {
  if (!matching_extents(grid, storage) ||
      !valid_active_range(grid, active_lower, active_upper)) {
    return -1;
  }
#pragma omp parallel for collapse(3) schedule(static)
  for (size_t i = active_lower[0]; i < active_upper[0]; ++i) {
    for (size_t j = active_lower[1]; j < active_upper[1]; ++j) {
      for (size_t k = active_lower[2]; k < active_upper[2]; ++k) {
        const size_t coordinate[3] = {i, j, k};
        const size_t cell = gamera_no_index3(grid->cell_extent, i, j, k);
        const double inverse_volume = 1.0 / grid->cell[cell].volume;
        for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
          double difference = 0.0;
          for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
            size_t upper_coordinate[3] = {coordinate[0], coordinate[1],
                                          coordinate[2]};
            ++upper_coordinate[direction];
            const size_t lower_face = gamera_no_index3(
                grid->face[direction].extent, coordinate[0], coordinate[1],
                coordinate[2]);
            const size_t upper_face = gamera_no_index3(
                grid->face[direction].extent, upper_coordinate[0],
                upper_coordinate[1], upper_coordinate[2]);
            difference +=
                storage->fluid_face_flux[direction]
                                        [lower_face * GAMERA_NO_FLUX_COUNT +
                                         (size_t)variable] -
                storage->fluid_face_flux[direction]
                                        [upper_face * GAMERA_NO_FLUX_COUNT +
                                         (size_t)variable];
          }
          storage->hydro_rate[cell * GAMERA_NO_FLUX_COUNT +
                              (size_t)variable] =
              difference * inverse_volume;
        }

        gamera_no_vec3 maxwell_rate = {{0.0, 0.0, 0.0}};
        if (use_mhd) {
          for (int component = 0; component < GAMERA_NO_DIM; ++component) {
            double difference = 0.0;
            for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
              size_t upper_coordinate[3] = {coordinate[0], coordinate[1],
                                            coordinate[2]};
              ++upper_coordinate[direction];
              const size_t lower_face = gamera_no_index3(
                  grid->face[direction].extent, coordinate[0], coordinate[1],
                  coordinate[2]);
              const size_t upper_face = gamera_no_index3(
                  grid->face[direction].extent, upper_coordinate[0],
                  upper_coordinate[1], upper_coordinate[2]);
              difference +=
                  storage->maxwell_face_flux[direction][lower_face]
                      .value[component] -
                  storage->maxwell_face_flux[direction][upper_face]
                      .value[component];
            }
            maxwell_rate.value[component] = difference * inverse_volume;
            if (background_cell_force != NULL) {
              maxwell_rate.value[component] +=
                  background_cell_force[cell].value[component];
            }
          }
        }
        storage->maxwell_rate[cell] = maxwell_rate;
      }
    }
  }
  return 0;
}

int gamera_no_sweep_edge_emf(
    const gamera_no_grid *grid, gamera_no_storage *storage,
    const double *conserved, const double *const face_magnetic_flux[3],
    int edge_direction, const size_t active_lower[3],
    const size_t active_upper[3], const gamera_no_emf_options *options,
    const gamera_no_background_field *background) {
  if (!matching_extents(grid, storage) || conserved == NULL ||
      face_magnetic_flux == NULL ||
      !valid_active_range(grid, active_lower, active_upper) ||
      !valid_emf_options(options) || edge_direction < 0 ||
      edge_direction >= GAMERA_NO_DIM) {
    return -1;
  }
  const int first_transverse = (edge_direction + 1) % GAMERA_NO_DIM;
  const int second_transverse = (edge_direction + 2) % GAMERA_NO_DIM;
  if (face_magnetic_flux[first_transverse] == NULL ||
      face_magnetic_flux[second_transverse] == NULL ||
      grid->cell_extent[first_transverse] < 4U ||
      grid->cell_extent[second_transverse] < 4U ||
      active_lower[first_transverse] < 4U ||
      active_lower[second_transverse] < 4U ||
      active_upper[first_transverse] >
          grid->cell_extent[first_transverse] - 4U ||
      active_upper[second_transverse] >
          grid->cell_extent[second_transverse] - 4U ||
      (options->use_background &&
       (background == NULL ||
        background->edge_magnetic[edge_direction] == NULL))) {
    return -1;
  }

  size_t edge_upper[3] = {active_upper[0], active_upper[1], active_upper[2]};
  ++edge_upper[first_transverse];
  ++edge_upper[second_transverse];
  int failed = 0;
#pragma omp parallel for collapse(3) reduction(| : failed) schedule(static)
  for (size_t i = active_lower[0]; i < edge_upper[0]; ++i) {
    for (size_t j = active_lower[1]; j < edge_upper[1]; ++j) {
      for (size_t k = active_lower[2]; k < edge_upper[2]; ++k) {
        const size_t coordinate[3] = {i, j, k};
        const size_t edge = gamera_no_index3(
            grid->edge[edge_direction].extent, i, j, k);
        if (grid->edge[edge_direction].valid[edge] == 0U) {
          failed = 1;
          continue;
        }
        const gamera_no_edge_geometry *edge_geometry =
            &grid->edge[edge_direction].value[edge];

        gamera_no_vec3 velocity;
        if (edge_velocity(grid, conserved, first_transverse,
                          second_transverse, coordinate,
                          options->density_floor, &velocity) != 0) {
          failed = 1;
          continue;
        }
        const double velocity_tangent1 =
            vector_dot(velocity, edge_geometry->tangent1);
        const double velocity_tangent2 =
            vector_dot(velocity, edge_geometry->tangent2);

        double first_left = 0.0;
        double first_right = 0.0;
        double second_left = 0.0;
        double second_right = 0.0;
        double first_normal[2];
        double second_normal[2];
        if (transverse_face_field(
                grid, face_magnetic_flux[first_transverse], first_transverse,
                second_transverse, coordinate, options->pdm_coefficient,
                &first_left, &first_right, first_normal) != 0 ||
            transverse_face_field(
                grid, face_magnetic_flux[second_transverse],
                second_transverse, first_transverse, coordinate,
                options->pdm_coefficient, &second_left, &second_right,
                second_normal) != 0) {
          failed = 1;
          continue;
        }
        double edge_field[2];
        if (gamera_no_solve_edge_field(
                first_normal, second_normal,
                0.5 * (first_left + first_right),
                0.5 * (second_left + second_right), edge_field) != 0) {
          failed = 1;
          continue;
        }
        if (options->use_background) {
          const gamera_no_vec3 background_field =
              background->edge_magnetic[edge_direction][edge];
          edge_field[0] +=
              vector_dot(background_field, edge_geometry->tangent1);
          edge_field[1] +=
              vector_dot(background_field, edge_geometry->tangent2);
        }
        const double current_jump =
            (second_right - second_left) - (first_right - first_left);
        double density = 0.0;
        if (edge_density(grid, conserved, first_transverse,
                         second_transverse, coordinate, &density) != 0) {
          failed = 1;
          continue;
        }
        double diffusion_speed = 0.0;
        if (gamera_no_compute_edge_emf(
                velocity_tangent1, velocity_tangent2, edge_field[0],
                edge_field[1], current_jump, density, edge_geometry,
                options->diffusion_coefficient, options->use_boris,
                options->light_speed, options->cfl, options->dt,
                &storage->edge_emf[edge_direction][edge],
                &diffusion_speed) != 0) {
          failed = 1;
        }
      }
    }
  }
  return failed ? -1 : 0;
}
