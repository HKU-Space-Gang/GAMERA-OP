#include "nonorthogonal_initialization.h"

#include "nonorthogonal_state.h"
#include "nonorthogonal_step.h"

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

enum { GAMERA_NO_INIT_GAUSS_POINTS = 12 };

static const double gauss_node[GAMERA_NO_INIT_GAUSS_POINTS] = {
    0.1252334085114689154724,  0.3678314989981801937527,
    0.5873179542866174472967,  0.7699026741943046870369,
    0.9041172563704748566785,  0.9815606342467192506906,
    -0.1252334085114689154724, -0.3678314989981801937527,
    -0.5873179542866174472967, -0.7699026741943046870369,
    -0.9041172563704748566785, -0.9815606342467192506906};

static const double gauss_weight[GAMERA_NO_INIT_GAUSS_POINTS] = {
    0.2491470458134027850006, 0.2334925365383548087610,
    0.2031674267230659217490, 0.1600783285433462263350,
    0.1069393259953184309603, 0.0471753363865118271946,
    0.2491470458134027850006, 0.2334925365383548087610,
    0.2031674267230659217490, 0.1600783285433462263350,
    0.1069393259953184309603, 0.0471753363865118271946};

static const int logical_sign[8][GAMERA_NO_DIM] = {
    {-1, -1, -1}, {+1, -1, -1}, {+1, +1, -1}, {-1, +1, -1},
    {-1, -1, +1}, {+1, -1, +1}, {+1, +1, +1}, {-1, +1, +1}};

static gamera_no_vec3 add(gamera_no_vec3 a, gamera_no_vec3 b) {
  gamera_no_vec3 result;
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    result.value[d] = a.value[d] + b.value[d];
  }
  return result;
}

static gamera_no_vec3 subtract(gamera_no_vec3 a, gamera_no_vec3 b) {
  gamera_no_vec3 result;
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    result.value[d] = a.value[d] - b.value[d];
  }
  return result;
}

static gamera_no_vec3 scale(gamera_no_vec3 value, double factor) {
  gamera_no_vec3 result;
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    result.value[d] = factor * value.value[d];
  }
  return result;
}

static double dot(gamera_no_vec3 a, gamera_no_vec3 b) {
  double result = 0.0;
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    result += a.value[d] * b.value[d];
  }
  return result;
}

static gamera_no_vec3 cross(gamera_no_vec3 a, gamera_no_vec3 b) {
  return (gamera_no_vec3){{
      a.value[1] * b.value[2] - a.value[2] * b.value[1],
      a.value[2] * b.value[0] - a.value[0] * b.value[2],
      a.value[0] * b.value[1] - a.value[1] * b.value[0]}};
}

static void cell_corners(const gamera_no_grid *grid, size_t i, size_t j,
                         size_t k, gamera_no_vec3 corners[8]) {
  const size_t offset[8][3] = {
      {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
      {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
  for (int corner = 0; corner < 8; ++corner) {
    corners[corner] =
        grid->vertex[gamera_no_index3(
            grid->vertex_extent, i + offset[corner][0],
            j + offset[corner][1], k + offset[corner][2])];
  }
}

static int average_primitive(const gamera_no_vec3 corners[8], double volume,
                             gamera_no_primitive_function function,
                             void *context,
                             gamera_no_primitive *average) {
  *average = (gamera_no_primitive){0};
  double integrated_volume = 0.0;
  for (int i = 0; i < GAMERA_NO_INIT_GAUSS_POINTS; ++i) {
    const double xi = gauss_node[i];
    for (int j = 0; j < GAMERA_NO_INIT_GAUSS_POINTS; ++j) {
      const double eta = gauss_node[j];
      for (int k = 0; k < GAMERA_NO_INIT_GAUSS_POINTS; ++k) {
        const double zeta = gauss_node[k];
        gamera_no_vec3 point = {{0.0, 0.0, 0.0}};
        gamera_no_vec3 derivative[3] = {{{0.0, 0.0, 0.0}},
                                         {{0.0, 0.0, 0.0}},
                                         {{0.0, 0.0, 0.0}}};
        for (int corner = 0; corner < 8; ++corner) {
          const double sx = (double)logical_sign[corner][0];
          const double sy = (double)logical_sign[corner][1];
          const double sz = (double)logical_sign[corner][2];
          const double fx = 1.0 + sx * xi;
          const double fy = 1.0 + sy * eta;
          const double fz = 1.0 + sz * zeta;
          const double shape = 0.125 * fx * fy * fz;
          const double dshape[3] = {0.125 * sx * fy * fz,
                                    0.125 * fx * sy * fz,
                                    0.125 * fx * fy * sz};
          point = add(point, scale(corners[corner], shape));
          for (int d = 0; d < GAMERA_NO_DIM; ++d) {
            derivative[d] =
                add(derivative[d], scale(corners[corner], dshape[d]));
          }
        }
        const double jacobian = fabs(dot(derivative[0],
                                         cross(derivative[1], derivative[2])));
        const double weight = jacobian * gauss_weight[i] * gauss_weight[j] *
                              gauss_weight[k];
        gamera_no_primitive primitive;
        if (!isfinite(weight) || function(point, context, &primitive) != 0) {
          return -1;
        }
        integrated_volume += weight;
        average->density += weight * primitive.density;
        average->pressure += weight * primitive.pressure;
        for (int d = 0; d < GAMERA_NO_DIM; ++d) {
          average->velocity.value[d] +=
              weight * primitive.velocity.value[d];
        }
      }
    }
  }
  if (!isfinite(integrated_volume) || integrated_volume <= 0.0 ||
      fabs(integrated_volume - volume) >
          5.0e-12 * fmax(integrated_volume, volume)) {
    return -1;
  }
  average->density /= integrated_volume;
  average->pressure /= integrated_volume;
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    average->velocity.value[d] /= integrated_volume;
  }
  return 0;
}

int gamera_no_initialize_primitives(
    const gamera_no_grid *grid, gamera_no_storage *storage,
    gamera_no_primitive_function function, void *context,
    bool volume_average, double gamma, double density_floor,
    double pressure_floor) {
  if (grid == NULL || storage == NULL || function == NULL ||
      storage->conserved == NULL) {
    return -1;
  }
  int failed = 0;
#pragma omp parallel for collapse(3) reduction(| : failed) schedule(static)
  for (size_t i = 0; i < grid->cell_extent[0]; ++i) {
    for (size_t j = 0; j < grid->cell_extent[1]; ++j) {
      for (size_t k = 0; k < grid->cell_extent[2]; ++k) {
        const size_t cell = gamera_no_index3(grid->cell_extent, i, j, k);
        gamera_no_primitive primitive;
        int status;
        if (volume_average) {
          gamera_no_vec3 corners[8];
          cell_corners(grid, i, j, k, corners);
          status = average_primitive(corners, grid->cell[cell].volume,
                                     function, context, &primitive);
        } else {
          status = function(grid->cell[cell].centroid, context, &primitive);
        }
        if (status != 0 ||
            gamera_no_primitive_to_conserved(
                &primitive, gamma, density_floor, pressure_floor,
                &storage->conserved[cell * GAMERA_NO_FLUX_COUNT]) != 0) {
          failed = 1;
        }
      }
    }
  }
  return failed ? -1 : 0;
}

static int allocate_doubles(size_t count, double **array) {
  if (array == NULL || count == 0 || count > SIZE_MAX / sizeof(**array)) {
    return -1;
  }
  *array = (double *)calloc(count, sizeof(**array));
  return *array == NULL ? -1 : 0;
}

static int integrate_edge(gamera_no_vec3 start, gamera_no_vec3 end,
                          gamera_no_vector_potential_function function,
                          void *context, double *integral) {
  const gamera_no_vec3 displacement = subtract(end, start);
  const gamera_no_vec3 midpoint = scale(add(start, end), 0.5);
  *integral = 0.0;
  for (int n = 0; n < GAMERA_NO_INIT_GAUSS_POINTS; ++n) {
    const gamera_no_vec3 point =
        add(midpoint, scale(displacement, 0.5 * gauss_node[n]));
    gamera_no_vec3 potential;
    if (function(point, context, &potential) != 0) {
      return -1;
    }
    *integral += 0.5 * gauss_weight[n] * dot(potential, displacement);
  }
  return isfinite(*integral) ? 0 : -1;
}

static gamera_no_vec3 vertex_at(const gamera_no_grid *grid, size_t i,
                                size_t j, size_t k) {
  return grid->vertex[gamera_no_index3(grid->vertex_extent, i, j, k)];
}

int gamera_no_initialize_flux_from_vector_potential(
    const gamera_no_grid *grid, gamera_no_storage *storage,
    gamera_no_vector_potential_function function, void *context) {
  if (grid == NULL || storage == NULL || function == NULL) {
    return -1;
  }
  double *line[3] = {NULL, NULL, NULL};
  int status = -1;
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t count =
        gamera_no_element_count3(grid->edge[direction].extent);
    if (storage->face_flux[direction] == NULL ||
        allocate_doubles(count, &line[direction]) != 0) {
      goto cleanup;
    }
    int failed = 0;
#pragma omp parallel for collapse(3) reduction(| : failed) schedule(static)
    for (size_t i = 0; i < grid->edge[direction].extent[0]; ++i) {
      for (size_t j = 0; j < grid->edge[direction].extent[1]; ++j) {
        for (size_t k = 0; k < grid->edge[direction].extent[2]; ++k) {
          size_t end[3] = {i, j, k};
          ++end[direction];
          const size_t edge = gamera_no_index3(
              grid->edge[direction].extent, i, j, k);
          if (integrate_edge(vertex_at(grid, i, j, k),
                             vertex_at(grid, end[0], end[1], end[2]),
                             function, context, &line[direction][edge]) != 0) {
            failed = 1;
          }
        }
      }
    }
    if (failed) {
      goto cleanup;
    }
  }

#pragma omp parallel for collapse(3) schedule(static)
  for (size_t i = 0; i < grid->face[GAMERA_NO_I].extent[0]; ++i) {
    for (size_t j = 0; j < grid->face[GAMERA_NO_I].extent[1]; ++j) {
      for (size_t k = 0; k < grid->face[GAMERA_NO_I].extent[2]; ++k) {
        const size_t face = gamera_no_index3(
            grid->face[GAMERA_NO_I].extent, i, j, k);
        const size_t j0 = gamera_no_index3(
            grid->edge[GAMERA_NO_J].extent, i, j, k);
        const size_t j1 = gamera_no_index3(
            grid->edge[GAMERA_NO_J].extent, i, j, k + 1U);
        const size_t k0 = gamera_no_index3(
            grid->edge[GAMERA_NO_K].extent, i, j, k);
        const size_t k1 = gamera_no_index3(
            grid->edge[GAMERA_NO_K].extent, i, j + 1U, k);
        storage->face_flux[GAMERA_NO_I][face] =
            line[GAMERA_NO_J][j0] - line[GAMERA_NO_J][j1] +
            line[GAMERA_NO_K][k1] - line[GAMERA_NO_K][k0];
      }
    }
  }
#pragma omp parallel for collapse(3) schedule(static)
  for (size_t i = 0; i < grid->face[GAMERA_NO_J].extent[0]; ++i) {
    for (size_t j = 0; j < grid->face[GAMERA_NO_J].extent[1]; ++j) {
      for (size_t k = 0; k < grid->face[GAMERA_NO_J].extent[2]; ++k) {
        const size_t face = gamera_no_index3(
            grid->face[GAMERA_NO_J].extent, i, j, k);
        const size_t i0 = gamera_no_index3(
            grid->edge[GAMERA_NO_I].extent, i, j, k);
        const size_t i1 = gamera_no_index3(
            grid->edge[GAMERA_NO_I].extent, i, j, k + 1U);
        const size_t k0 = gamera_no_index3(
            grid->edge[GAMERA_NO_K].extent, i, j, k);
        const size_t k1 = gamera_no_index3(
            grid->edge[GAMERA_NO_K].extent, i + 1U, j, k);
        storage->face_flux[GAMERA_NO_J][face] =
            -line[GAMERA_NO_I][i0] + line[GAMERA_NO_I][i1] -
            line[GAMERA_NO_K][k1] + line[GAMERA_NO_K][k0];
      }
    }
  }
#pragma omp parallel for collapse(3) schedule(static)
  for (size_t i = 0; i < grid->face[GAMERA_NO_K].extent[0]; ++i) {
    for (size_t j = 0; j < grid->face[GAMERA_NO_K].extent[1]; ++j) {
      for (size_t k = 0; k < grid->face[GAMERA_NO_K].extent[2]; ++k) {
        const size_t face = gamera_no_index3(
            grid->face[GAMERA_NO_K].extent, i, j, k);
        const size_t i0 = gamera_no_index3(
            grid->edge[GAMERA_NO_I].extent, i, j, k);
        const size_t i1 = gamera_no_index3(
            grid->edge[GAMERA_NO_I].extent, i, j + 1U, k);
        const size_t j0 = gamera_no_index3(
            grid->edge[GAMERA_NO_J].extent, i, j, k);
        const size_t j1 = gamera_no_index3(
            grid->edge[GAMERA_NO_J].extent, i + 1U, j, k);
        storage->face_flux[GAMERA_NO_K][face] =
            line[GAMERA_NO_I][i0] - line[GAMERA_NO_I][i1] +
            line[GAMERA_NO_J][j1] - line[GAMERA_NO_J][j0];
      }
    }
  }

  {
    const double *face_flux[3] = {storage->face_flux[0],
                                  storage->face_flux[1],
                                  storage->face_flux[2]};
    if (gamera_no_recover_magnetic_field(grid, face_flux,
                                         storage->cell_magnetic) != 0) {
      goto cleanup;
    }
  }
  status = 0;

cleanup:
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    free(line[direction]);
  }
  return status;
}

static int orszag_tang_primitive(gamera_no_vec3 point, void *context,
                                 gamera_no_primitive *primitive) {
  const gamera_no_orszag_tang_options *options =
      (const gamera_no_orszag_tang_options *)context;
  const double pi = acos(-1.0);
  primitive->density = options->density;
  primitive->pressure = options->pressure;
  primitive->velocity =
      (gamera_no_vec3){{-sin(2.0 * pi * point.value[1]),
                        sin(2.0 * pi * point.value[0]), 0.0}};
  return 0;
}

static int orszag_tang_potential(gamera_no_vec3 point, void *context,
                                 gamera_no_vec3 *potential) {
  const gamera_no_orszag_tang_options *options =
      (const gamera_no_orszag_tang_options *)context;
  const double pi = acos(-1.0);
  *potential = (gamera_no_vec3){{
      0.0, 0.0,
      options->magnetic_amplitude *
          (cos(4.0 * pi * point.value[0]) / (4.0 * pi) +
           cos(2.0 * pi * point.value[1]) / (2.0 * pi))}};
  return 0;
}

void gamera_no_orszag_tang_fortran_defaults(
    gamera_no_orszag_tang_options *options) {
  if (options == NULL) {
    return;
  }
  const double pi = acos(-1.0);
  options->density = 25.0 / (36.0 * pi);
  options->pressure = 5.0 / (12.0 * pi);
  options->magnetic_amplitude = 1.0 / sqrt(4.0 * pi);
}

void gamera_no_orszag_tang_paper_defaults(
    gamera_no_orszag_tang_options *options) {
  if (options == NULL) {
    return;
  }
  const double pi = acos(-1.0);
  options->density = 25.0 * pi / 36.0;
  options->pressure = 5.0 * pi / 12.0;
  options->magnetic_amplitude = 1.0;
}

int gamera_no_initialize_orszag_tang(
    const gamera_no_grid *grid, gamera_no_storage *storage,
    const gamera_no_orszag_tang_options *options, double gamma,
    double density_floor, double pressure_floor) {
  if (options == NULL || !isfinite(options->density) ||
      !isfinite(options->pressure) ||
      !isfinite(options->magnetic_amplitude) || options->density <= 0.0 ||
      options->pressure <= 0.0) {
    return -1;
  }
  if (gamera_no_initialize_primitives(
          grid, storage, orszag_tang_primitive, (void *)options, false, gamma,
          density_floor, pressure_floor) != 0 ||
      gamera_no_initialize_flux_from_vector_potential(
          grid, storage, orszag_tang_potential, (void *)options) != 0 ||
      gamera_no_save_current_as_old(storage) != 0) {
    return -1;
  }
  return 0;
}

static int blast_primitive(gamera_no_vec3 point, void *context,
                           gamera_no_primitive *primitive) {
  const gamera_no_blast_options *options =
      (const gamera_no_blast_options *)context;
  const gamera_no_vec3 offset = subtract(point, options->center);
  const double radius = sqrt(dot(offset, offset));
  const int inside = radius <= options->radius;
  primitive->density =
      options->ambient_density * (inside ? options->density_ratio : 1.0);
  primitive->pressure =
      options->ambient_pressure * (inside ? options->pressure_ratio : 1.0);
  primitive->velocity = options->velocity;
  return 0;
}

static int uniform_field_potential(gamera_no_vec3 point, void *context,
                                   gamera_no_vec3 *potential) {
  const gamera_no_blast_options *options =
      (const gamera_no_blast_options *)context;
  *potential = scale(cross(options->magnetic, point), 0.5);
  return 0;
}

void gamera_no_spherical_blast_defaults(gamera_no_blast_options *options) {
  if (options == NULL) {
    return;
  }
  *options = (gamera_no_blast_options){
      .center = {{1.0, 0.0, 0.0}},
      .radius = 0.1,
      .ambient_density = 1.0,
      .ambient_pressure = 0.1,
      .density_ratio = 1.0,
      .pressure_ratio = 100.0,
      .velocity = {{0.0, 0.0, 0.0}},
      .magnetic = {{0.70710678118654752440,
                    0.70710678118654752440, 0.0}},
      .volume_average = true};
}

int gamera_no_initialize_blast(
    const gamera_no_grid *grid, gamera_no_storage *storage,
    const gamera_no_blast_options *options, double gamma,
    double density_floor, double pressure_floor) {
  if (options == NULL || !isfinite(options->radius) ||
      !isfinite(options->ambient_density) ||
      !isfinite(options->ambient_pressure) ||
      !isfinite(options->density_ratio) ||
      !isfinite(options->pressure_ratio) || options->radius <= 0.0 ||
      options->ambient_density <= 0.0 || options->ambient_pressure <= 0.0 ||
      options->density_ratio <= 0.0 || options->pressure_ratio <= 0.0) {
    return -1;
  }
  if (gamera_no_initialize_primitives(
          grid, storage, blast_primitive, (void *)options,
          options->volume_average, gamma, density_floor, pressure_floor) !=
          0 ||
      gamera_no_initialize_flux_from_vector_potential(
          grid, storage, uniform_field_potential, (void *)options) != 0 ||
      gamera_no_save_current_as_old(storage) != 0) {
    return -1;
  }
  return 0;
}
