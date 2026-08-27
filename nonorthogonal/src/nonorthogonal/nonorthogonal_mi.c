#include "nonorthogonal_mi.h"

#include "nonorthogonal_operators.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>


static double vector_dot(gamera_no_vec3 left, gamera_no_vec3 right) {
  double result = 0.0;
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    result += left.value[component] * right.value[component];
  }
  return result;
}

int gamera_mi_compose_ghost_velocity(gamera_no_vec3 source_velocity,
                                     gamera_no_vec3 background_magnetic,
                                     gamera_no_vec3 drift_velocity,
                                     gamera_no_vec3 *ghost_velocity) {
  if (ghost_velocity == NULL) {
    return -1;
  }
  const double background_square =
      vector_dot(background_magnetic, background_magnetic);
  if (!(background_square > 0.0) || !isfinite(background_square)) {
    return -1;
  }
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    ghost_velocity->value[component] =
        2.0 * drift_velocity.value[component] -
        source_velocity.value[component];
    if (!isfinite(ghost_velocity->value[component])) {
      return -1;
    }
  }
  return 0;
}

static int edge_integral(const gamera_no_grid *grid,
                         const gamera_no_vec3 *magnetic, int direction,
                         size_t i, size_t j, size_t k, double *result) {
  if (grid == NULL || magnetic == NULL || result == NULL || direction < 0 ||
      direction >= GAMERA_NO_DIM ||
      i >= grid->edge[direction].extent[0] ||
      j >= grid->edge[direction].extent[1] ||
      k >= grid->edge[direction].extent[2]) {
    return -1;
  }
  const size_t edge =
      gamera_no_index3(grid->edge[direction].extent, i, j, k);
  if (grid->edge[direction].valid[edge] == 0U) {
    return -1;
  }

  const size_t cell_coordinates[4][3] = {
      {i, j, k}, {i, j, k}, {i, j, k}, {i, j, k}};
  size_t stencil[4][3];
  memcpy(stencil, cell_coordinates, sizeof(stencil));
  if (direction == GAMERA_NO_I) {
    --stencil[1][1];
    --stencil[2][2];
    --stencil[3][1];
    --stencil[3][2];
  } else if (direction == GAMERA_NO_J) {
    --stencil[1][0];
    --stencil[2][2];
    --stencil[3][0];
    --stencil[3][2];
  } else {
    --stencil[1][0];
    --stencil[2][1];
    --stencil[3][0];
    --stencil[3][1];
  }

  gamera_no_vec3 average = {{0.0, 0.0, 0.0}};
  for (int point = 0; point < 4; ++point) {
    for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
      if (stencil[point][axis] >= grid->cell_extent[axis]) {
        return -1;
      }
    }
    const gamera_no_vec3 value = magnetic[gamera_no_index3(
        grid->cell_extent, stencil[point][0], stencil[point][1],
        stencil[point][2])];
    for (int component = 0; component < GAMERA_NO_DIM; ++component) {
      average.value[component] += 0.25 * value.value[component];
    }
  }
  const gamera_no_edge_geometry *geometry =
      &grid->edge[direction].value[edge];
  *result = vector_dot(average, geometry->normal) * geometry->length;
  return isfinite(*result) ? 0 : -1;
}

static int face_current_flux(const gamera_no_grid *grid,
                             const gamera_no_vec3 *magnetic, int direction,
                             size_t i, size_t j, size_t k, double *result) {
  double first, second, third, fourth;
  if (direction == GAMERA_NO_I) {
    if (edge_integral(grid, magnetic, GAMERA_NO_J, i, j, k, &first) != 0 ||
        edge_integral(grid, magnetic, GAMERA_NO_K, i, j + 1U, k,
                      &second) != 0 ||
        edge_integral(grid, magnetic, GAMERA_NO_J, i, j, k + 1U,
                      &third) != 0 ||
        edge_integral(grid, magnetic, GAMERA_NO_K, i, j, k, &fourth) != 0) {
      return -1;
    }
  } else if (direction == GAMERA_NO_J) {
    if (edge_integral(grid, magnetic, GAMERA_NO_K, i, j, k, &first) != 0 ||
        edge_integral(grid, magnetic, GAMERA_NO_I, i, j, k + 1U,
                      &second) != 0 ||
        edge_integral(grid, magnetic, GAMERA_NO_K, i + 1U, j, k,
                      &third) != 0 ||
        edge_integral(grid, magnetic, GAMERA_NO_I, i, j, k, &fourth) != 0) {
      return -1;
    }
  } else if (direction == GAMERA_NO_K) {
    if (edge_integral(grid, magnetic, GAMERA_NO_I, i, j, k, &first) != 0 ||
        edge_integral(grid, magnetic, GAMERA_NO_J, i + 1U, j, k,
                      &second) != 0 ||
        edge_integral(grid, magnetic, GAMERA_NO_I, i, j + 1U, k,
                      &third) != 0 ||
        edge_integral(grid, magnetic, GAMERA_NO_J, i, j, k, &fourth) != 0) {
      return -1;
    }
  } else {
    return -1;
  }
  *result = first + second - third - fourth;
  return isfinite(*result) ? 0 : -1;
}

int gamera_no_cell_current_from_residual(
    const gamera_no_grid *grid, const gamera_no_vec3 *residual_magnetic,
    size_t i, size_t j, size_t k, gamera_no_vec3 *current) {
  if (grid == NULL || residual_magnetic == NULL || current == NULL ||
      i == 0U || j == 0U || k == 0U ||
      i + 1U >= grid->cell_extent[0] ||
      j + 1U >= grid->cell_extent[1] ||
      k + 1U >= grid->cell_extent[2]) {
    return -1;
  }
  double flux[GAMERA_NO_DIM][2];
  if (face_current_flux(grid, residual_magnetic, GAMERA_NO_I, i, j, k,
                        &flux[GAMERA_NO_I][GAMERA_NO_LOWER]) != 0 ||
      face_current_flux(grid, residual_magnetic, GAMERA_NO_I, i + 1U, j, k,
                        &flux[GAMERA_NO_I][GAMERA_NO_UPPER]) != 0 ||
      face_current_flux(grid, residual_magnetic, GAMERA_NO_J, i, j, k,
                        &flux[GAMERA_NO_J][GAMERA_NO_LOWER]) != 0 ||
      face_current_flux(grid, residual_magnetic, GAMERA_NO_J, i, j + 1U, k,
                        &flux[GAMERA_NO_J][GAMERA_NO_UPPER]) != 0 ||
      face_current_flux(grid, residual_magnetic, GAMERA_NO_K, i, j, k,
                        &flux[GAMERA_NO_K][GAMERA_NO_LOWER]) != 0 ||
      face_current_flux(grid, residual_magnetic, GAMERA_NO_K, i, j, k + 1U,
                        &flux[GAMERA_NO_K][GAMERA_NO_UPPER]) != 0) {
    return -1;
  }
  const gamera_no_cell_geometry *cell =
      &grid->cell[gamera_no_index3(grid->cell_extent, i, j, k)];
  return gamera_no_flux_to_cell_field(cell, flux, current, NULL);
}

double gamera_mi_mapped_colatitude(double mhd_colatitude,
                                   double mhd_radius,
                                   double ionosphere_radius) {
  if (!isfinite(mhd_colatitude) || !isfinite(mhd_radius) ||
      !isfinite(ionosphere_radius) || mhd_radius <= 0.0 ||
      ionosphere_radius <= 0.0 || ionosphere_radius > mhd_radius) {
    return NAN;
  }
  double argument =
      sin(mhd_colatitude) * sqrt(ionosphere_radius / mhd_radius);
  argument = fmax(-1.0, fmin(1.0, argument));
  return asin(argument);
}

double gamera_mi_dipole_field_ratio(double mhd_colatitude,
                                    double mhd_radius,
                                    double ionosphere_radius) {
  const double ionosphere_colatitude = gamera_mi_mapped_colatitude(
      mhd_colatitude, mhd_radius, ionosphere_radius);
  if (!isfinite(ionosphere_colatitude)) {
    return NAN;
  }
  const double radial_ratio = mhd_radius / ionosphere_radius;
  const double numerator =
      sqrt(1.0 + 3.0 * pow(cos(ionosphere_colatitude), 2.0));
  const double denominator =
      sqrt(1.0 + 3.0 * pow(cos(mhd_colatitude), 2.0));
  return radial_ratio * radial_ratio * radial_ratio * numerator /
         denominator;
}

double gamera_mi_dipole_cos_inclination(double ionosphere_colatitude,
                                        int hemisphere) {
  if (!isfinite(ionosphere_colatitude) || ionosphere_colatitude < 0.0 ||
      ionosphere_colatitude > 0.5 * acos(-1.0) ||
      (hemisphere != GAMERA_MI_NORTH && hemisphere != GAMERA_MI_SOUTH)) {
    return NAN;
  }
  const double cosine = cos(ionosphere_colatitude);
  return (double)hemisphere * 2.0 * cosine /
         sqrt(1.0 + 3.0 * cosine * cosine);
}

static int valid_solver_config(const gamera_mi_solver_config *config) {
  return config != NULL && config->longitude_count >= 4U &&
         config->colatitude_count >= 3U &&
         isfinite(config->maximum_colatitude) &&
         config->maximum_colatitude > 0.0 &&
         config->maximum_colatitude <= 0.5 * acos(-1.0) &&
         isfinite(config->ionosphere_radius_m) &&
         config->ionosphere_radius_m > 0.0 &&
         isfinite(config->pedersen_siemens) &&
         config->pedersen_siemens > 0.0 &&
         config->hall_siemens == 0.0 &&
         isfinite(config->low_latitude_potential_v) &&
         (config->hemisphere == GAMERA_MI_NORTH ||
          config->hemisphere == GAMERA_MI_SOUTH) &&
         config->maximum_iterations > 0 &&
         isfinite(config->relative_tolerance) &&
         config->relative_tolerance > 0.0 &&
         isfinite(config->absolute_tolerance) &&
         config->absolute_tolerance >= 0.0;
}

static size_t polar_index(const gamera_mi_solver_config *config, size_t i,
                          size_t j) {
  return i * config->longitude_count + j;
}

static double theta_coefficient(const gamera_mi_solver_config *config,
                                double theta_half, double theta_center,
                                double delta_theta) {
  const double cosine = cos(theta_half);
  const double cosine_inclination_square =
      4.0 * cosine * cosine / (1.0 + 3.0 * cosine * cosine);
  return sin(theta_half) * config->pedersen_siemens /
         (cosine_inclination_square * sin(theta_center) *
          delta_theta * delta_theta);
}

int gamera_mi_apply_constant_pedersen(
    const gamera_mi_solver_config *config, const double *potential_v,
    double *result) {
  if (!valid_solver_config(config) || potential_v == NULL || result == NULL) {
    return -1;
  }
  const size_t np = config->longitude_count;
  const size_t nt = config->colatitude_count;
  const double dtheta = config->maximum_colatitude / (double)(nt - 1U);
  const double dphi = 2.0 * acos(-1.0) / (double)np;
  double first_ring_mean = 0.0;
  for (size_t j = 0; j < np; ++j) {
    first_ring_mean += potential_v[polar_index(config, 1U, j)] / (double)np;
  }
  for (size_t j = 0; j < np; ++j) {
    result[polar_index(config, 0U, j)] =
        potential_v[polar_index(config, 0U, j)] - first_ring_mean;
  }
  for (size_t i = 1U; i + 1U < nt; ++i) {
    const double theta = (double)i * dtheta;
    const double lower =
        theta_coefficient(config, theta - 0.5 * dtheta, theta, dtheta);
    const double upper =
        theta_coefficient(config, theta + 0.5 * dtheta, theta, dtheta);
    const double longitude = config->pedersen_siemens /
        (sin(theta) * sin(theta) * dphi * dphi);
    for (size_t j = 0; j < np; ++j) {
      const size_t jm = j == 0U ? np - 1U : j - 1U;
      const size_t jp = j + 1U == np ? 0U : j + 1U;
      const size_t center = polar_index(config, i, j);
      result[center] =
          lower * (potential_v[center] -
                   potential_v[polar_index(config, i - 1U, j)]) +
          upper * (potential_v[center] -
                   potential_v[polar_index(config, i + 1U, j)]) +
          longitude *
              (2.0 * potential_v[center] -
               potential_v[polar_index(config, i, jm)] -
               potential_v[polar_index(config, i, jp)]);
    }
  }
  for (size_t j = 0; j < np; ++j) {
    result[polar_index(config, nt - 1U, j)] =
        potential_v[polar_index(config, nt - 1U, j)];
  }
  return 0;
}

static double array_dot(const double *left, const double *right,
                        size_t count) {
  double result = 0.0;
  for (size_t i = 0; i < count; ++i) {
    result += left[i] * right[i];
  }
  return result;
}

/*
 * There is only one physical node at the pole.  Keeping one copy per
 * longitude and constraining every copy to the first-ring mean produces a
 * nonsymmetric algebraic system even though the constant-Pedersen operator
 * is elliptic.  The old BiCGStab solve occasionally stagnated on that
 * representation at production resolution.
 *
 * Use one pole unknown and multiply each interior equation by sin(theta).
 * The resulting reduced operator is symmetric positive definite: adjacent
 * theta rows share the same face coefficient and the pole row is scaled so
 * that its coupling to every first-ring node has that coefficient.  This is
 * algebraically equivalent to the full operator once the single pole value
 * is copied to every longitude.
 */
static size_t reduced_index(const gamera_mi_solver_config *config, size_t i,
                            size_t j) {
  return 1U + (i - 1U) * config->longitude_count + j;
}

static double pole_edge_coefficient(
    const gamera_mi_solver_config *config) {
  const double dtheta = config->maximum_colatitude /
                        (double)(config->colatitude_count - 1U);
  const double theta_half = 0.5 * dtheta;
  const double cosine = cos(theta_half);
  const double cosine_inclination_square =
      4.0 * cosine * cosine / (1.0 + 3.0 * cosine * cosine);
  return sin(theta_half) * config->pedersen_siemens /
         (cosine_inclination_square * dtheta * dtheta);
}

static int build_reduced_diagonal(const gamera_mi_solver_config *config,
                                  double *diagonal) {
  const size_t np = config->longitude_count;
  const size_t nt = config->colatitude_count;
  const double dtheta = config->maximum_colatitude / (double)(nt - 1U);
  const double dphi = 2.0 * acos(-1.0) / (double)np;
  diagonal[0] = (double)np * pole_edge_coefficient(config);
  for (size_t i = 1U; i + 1U < nt; ++i) {
    const double theta = (double)i * dtheta;
    const double weight = sin(theta);
    const double lower =
        theta_coefficient(config, theta - 0.5 * dtheta, theta, dtheta);
    const double upper =
        theta_coefficient(config, theta + 0.5 * dtheta, theta, dtheta);
    const double longitude = config->pedersen_siemens /
        (sin(theta) * sin(theta) * dphi * dphi);
    for (size_t j = 0; j < np; ++j) {
      diagonal[reduced_index(config, i, j)] =
          weight * (lower + upper + 2.0 * longitude);
    }
  }
  return 0;
}

static void apply_preconditioner(const double *diagonal,
                                 const double *input, double *output,
                                 size_t count) {
  for (size_t i = 0; i < count; ++i) {
    output[i] = input[i] / diagonal[i];
  }
}

static int apply_reduced_constant_pedersen(
    const gamera_mi_solver_config *config, const double *potential_v,
    double *result) {
  const size_t np = config->longitude_count;
  const size_t nt = config->colatitude_count;
  const double dtheta = config->maximum_colatitude / (double)(nt - 1U);
  const double dphi = 2.0 * acos(-1.0) / (double)np;
  const double pole_edge = pole_edge_coefficient(config);
  double first_ring_sum = 0.0;
  for (size_t j = 0; j < np; ++j) {
    first_ring_sum += potential_v[reduced_index(config, 1U, j)];
  }
  result[0] = pole_edge *
              ((double)np * potential_v[0] - first_ring_sum);
  for (size_t i = 1U; i + 1U < nt; ++i) {
    const double theta = (double)i * dtheta;
    const double weight = sin(theta);
    const double lower =
        theta_coefficient(config, theta - 0.5 * dtheta, theta, dtheta);
    const double upper =
        theta_coefficient(config, theta + 0.5 * dtheta, theta, dtheta);
    const double longitude = config->pedersen_siemens /
        (sin(theta) * sin(theta) * dphi * dphi);
    for (size_t j = 0; j < np; ++j) {
      const size_t jm = j == 0U ? np - 1U : j - 1U;
      const size_t jp = j + 1U == np ? 0U : j + 1U;
      const size_t center = reduced_index(config, i, j);
      const double lower_value =
          i == 1U ? potential_v[0]
                  : potential_v[reduced_index(config, i - 1U, j)];
      const double upper_value =
          i + 2U < nt
              ? potential_v[reduced_index(config, i + 1U, j)]
              : 0.0;
      result[center] =
          weight *
          (lower * (potential_v[center] - lower_value) +
           upper * (potential_v[center] - upper_value) +
           longitude *
               (2.0 * potential_v[center] -
                potential_v[reduced_index(config, i, jm)] -
                potential_v[reduced_index(config, i, jp)]));
    }
  }
  return 0;
}

static double reduced_original_residual_norm(
    const gamera_mi_solver_config *config, const double *residual) {
  const size_t np = config->longitude_count;
  const size_t nt = config->colatitude_count;
  const double dtheta = config->maximum_colatitude / (double)(nt - 1U);
  const double pole_scale = (double)np * pole_edge_coefficient(config);
  const double pole_residual = residual[0] / pole_scale;
  double norm_square = (double)np * pole_residual * pole_residual;
  for (size_t i = 1U; i + 1U < nt; ++i) {
    const double weight = sin((double)i * dtheta);
    for (size_t j = 0; j < np; ++j) {
      const double value = residual[reduced_index(config, i, j)] / weight;
      norm_square += value * value;
    }
  }
  return sqrt(fmax(0.0, norm_square));
}

int gamera_mi_solve_constant_pedersen(
    const gamera_mi_solver_config *config, const double *fac_a_m2,
    double *potential_v, gamera_mi_solver_stats *stats) {
  if (!valid_solver_config(config) || fac_a_m2 == NULL ||
      potential_v == NULL) {
    return -1;
  }
  const size_t np = config->longitude_count;
  const size_t nt = config->colatitude_count;
  if (np > (SIZE_MAX - 1U) / (nt - 2U)) {
    return -1;
  }
  const size_t count = 1U + np * (nt - 2U);
  double *work = (double *)calloc(7U * count, sizeof(double));
  if (work == NULL) {
    return -1;
  }
  double *right = work;
  double *diagonal = right + count;
  double *solution = diagonal + count;
  double *residual = solution + count;
  double *preconditioned = residual + count;
  double *search = preconditioned + count;
  double *operator_search = search + count;

  const double dtheta = config->maximum_colatitude / (double)(nt - 1U);
  const double radius_square =
      config->ionosphere_radius_m * config->ionosphere_radius_m;
  double right_norm_square =
      (double)np * config->low_latitude_potential_v *
      config->low_latitude_potential_v;
  for (size_t i = 1U; i + 1U < nt; ++i) {
    const double theta = (double)i * dtheta;
    const double weight = sin(theta);
    const double inclination =
        gamera_mi_dipole_cos_inclination(theta, config->hemisphere);
    const double upper =
        theta_coefficient(config, theta + 0.5 * dtheta, theta, dtheta);
    for (size_t j = 0; j < np; ++j) {
      const size_t full_index = polar_index(config, i, j);
      const size_t index = reduced_index(config, i, j);
      if (!isfinite(fac_a_m2[full_index])) {
        free(work);
        return -1;
      }
      const double unweighted_right =
          -radius_square * fac_a_m2[full_index] * inclination;
      right[index] = weight * unweighted_right;
      right_norm_square += unweighted_right * unweighted_right;
      if (i + 2U == nt) {
        right[index] += weight * upper *
                        config->low_latitude_potential_v;
      }
    }
  }
  double pole_mean = 0.0;
  for (size_t j = 0; j < np; ++j) {
    const double value = potential_v[polar_index(config, 0U, j)];
    if (isfinite(value)) {
      pole_mean += value / (double)np;
    }
  }
  solution[0] = pole_mean;
  for (size_t i = 1U; i + 1U < nt; ++i) {
    for (size_t j = 0; j < np; ++j) {
      const double value = potential_v[polar_index(config, i, j)];
      solution[reduced_index(config, i, j)] =
          isfinite(value) ? value : 0.0;
    }
  }
  build_reduced_diagonal(config, diagonal);
  if (apply_reduced_constant_pedersen(config, solution,
                                      operator_search) != 0) {
    free(work);
    return -1;
  }
  for (size_t i = 0; i < count; ++i) {
    residual[i] = right[i] - operator_search[i];
  }
  const double initial = reduced_original_residual_norm(config, residual);
  const double target = fmax(
      config->absolute_tolerance,
      config->relative_tolerance *
          fmax(sqrt(fmax(0.0, right_norm_square)), 1.0));
  if (stats != NULL) {
    *stats = (gamera_mi_solver_stats){0, initial, initial, initial <= target};
  }
  if (initial <= target) {
    for (size_t j = 0; j < np; ++j) {
      potential_v[polar_index(config, 0U, j)] = solution[0];
      potential_v[polar_index(config, nt - 1U, j)] =
          config->low_latitude_potential_v;
    }
    for (size_t i = 1U; i + 1U < nt; ++i) {
      for (size_t j = 0; j < np; ++j) {
        potential_v[polar_index(config, i, j)] =
            solution[reduced_index(config, i, j)];
      }
    }
    free(work);
    return 0;
  }

  apply_preconditioner(diagonal, residual, preconditioned, count);
  memcpy(search, preconditioned, count * sizeof(*search));
  double rho = array_dot(residual, preconditioned, count);
  double final = initial;
  int iterations = 0;
  int converged = 0;
  for (int iteration = 1; iteration <= config->maximum_iterations;
       ++iteration) {
    if (!isfinite(rho) || rho <= DBL_MIN ||
        apply_reduced_constant_pedersen(config, search,
                                        operator_search) != 0) {
      break;
    }
    const double denominator = array_dot(search, operator_search, count);
    if (!isfinite(denominator) || denominator <= DBL_MIN) {
      break;
    }
    const double alpha = rho / denominator;
    for (size_t i = 0; i < count; ++i) {
      solution[i] += alpha * search[i];
      residual[i] -= alpha * operator_search[i];
    }
    final = reduced_original_residual_norm(config, residual);
    iterations = iteration;
    if (final <= target || iteration % 1024 == 0) {
      /* Verify with a true residual before accepting the tight 1e-10 solve.
       * If recurrence drift is visible, this also gives PCG a clean restart.
       */
      if (apply_reduced_constant_pedersen(config, solution,
                                          operator_search) != 0) {
        break;
      }
      for (size_t i = 0; i < count; ++i) {
        residual[i] = right[i] - operator_search[i];
      }
      final = reduced_original_residual_norm(config, residual);
      if (final <= target) {
        converged = 1;
        break;
      }
      apply_preconditioner(diagonal, residual, preconditioned, count);
      memcpy(search, preconditioned, count * sizeof(*search));
      rho = array_dot(residual, preconditioned, count);
      continue;
    }
    apply_preconditioner(diagonal, residual, preconditioned, count);
    const double rho_new = array_dot(residual, preconditioned, count);
    if (!isfinite(rho_new) || rho_new <= DBL_MIN) {
      break;
    }
    const double beta = rho_new / rho;
    for (size_t i = 0; i < count; ++i) {
      search[i] = preconditioned[i] + beta * search[i];
    }
    rho = rho_new;
  }
  for (size_t j = 0; j < np; ++j) {
    potential_v[polar_index(config, 0U, j)] = solution[0];
    potential_v[polar_index(config, nt - 1U, j)] =
        config->low_latitude_potential_v;
  }
  for (size_t i = 1U; i + 1U < nt; ++i) {
    for (size_t j = 0; j < np; ++j) {
      potential_v[polar_index(config, i, j)] =
          solution[reduced_index(config, i, j)];
    }
  }
  if (stats != NULL) {
    *stats =
        (gamera_mi_solver_stats){iterations, initial, final, converged};
  }
  free(work);
  return converged ? 0 : -1;
}

static int valid_tensor_config(const gamera_mi_solver_config *config,
                               const double *pedersen,
                               const double *hall) {
  if (config == NULL || pedersen == NULL || hall == NULL ||
      config->longitude_count < 4U || config->colatitude_count < 3U ||
      !isfinite(config->maximum_colatitude) ||
      !(config->maximum_colatitude > 0.0) ||
      config->maximum_colatitude > 0.5 * acos(-1.0) ||
      !isfinite(config->ionosphere_radius_m) ||
      !(config->ionosphere_radius_m > 0.0) ||
      !isfinite(config->low_latitude_potential_v) ||
      fabs(config->low_latitude_potential_v) > DBL_MIN ||
      (config->hemisphere != GAMERA_MI_NORTH &&
       config->hemisphere != GAMERA_MI_SOUTH) ||
      config->maximum_iterations <= 0 ||
      !isfinite(config->relative_tolerance) ||
      !(config->relative_tolerance > 0.0) ||
      !isfinite(config->absolute_tolerance) ||
      config->absolute_tolerance < 0.0 ||
      config->longitude_count > SIZE_MAX / config->colatitude_count) {
    return 0;
  }
  const size_t count = config->longitude_count * config->colatitude_count;
  for (size_t index = 0U; index < count; ++index) {
    if (!isfinite(pedersen[index]) || !(pedersen[index] > 0.0) ||
        !isfinite(hall[index]) || hall[index] < 0.0) {
      return 0;
    }
  }
  return 1;
}

static double inclination_square(double theta) {
  const double cosine = cos(theta);
  return 4.0 * cosine * cosine / (1.0 + 3.0 * cosine * cosine);
}

static double tensor_f12(const gamera_mi_solver_config *config,
                         const double *hall, size_t i, size_t j) {
  const double dtheta = config->maximum_colatitude /
                        (double)(config->colatitude_count - 1U);
  const double inclination = gamera_mi_dipole_cos_inclination(
      (double)i * dtheta, config->hemisphere);
  return -hall[polar_index(config, i, j)] / inclination;
}

int gamera_mi_apply_conductance_tensor(
    const gamera_mi_solver_config *config, const double *pedersen,
    const double *hall, const double *potential_v, double *result) {
  if (!valid_tensor_config(config, pedersen, hall) || potential_v == NULL ||
      result == NULL) {
    return -1;
  }
  const size_t np = config->longitude_count;
  const size_t nt = config->colatitude_count;
  const double dtheta = config->maximum_colatitude / (double)(nt - 1U);
  const double dphi = 2.0 * acos(-1.0) / (double)np;
  double first_ring_mean = 0.0;
  for (size_t j = 0U; j < np; ++j) {
    first_ring_mean += potential_v[polar_index(config, 1U, j)] / (double)np;
  }
  for (size_t j = 0U; j < np; ++j) {
    result[polar_index(config, 0U, j)] =
        potential_v[polar_index(config, 0U, j)] - first_ring_mean;
  }
  for (size_t i = 1U; i + 1U < nt; ++i) {
    const double theta = (double)i * dtheta;
    const double sine = sin(theta);
    const double lower_theta = theta - 0.5 * dtheta;
    const double upper_theta = theta + 0.5 * dtheta;
    for (size_t j = 0U; j < np; ++j) {
      const size_t jm = j == 0U ? np - 1U : j - 1U;
      const size_t jp = j + 1U == np ? 0U : j + 1U;
      const size_t center = polar_index(config, i, j);
      const double lower =
          sin(lower_theta) *
          0.5 * (pedersen[center] +
                 pedersen[polar_index(config, i - 1U, j)]) /
          (inclination_square(lower_theta) * sine * dtheta * dtheta);
      const double upper =
          sin(upper_theta) *
          0.5 * (pedersen[center] +
                 pedersen[polar_index(config, i + 1U, j)]) /
          (inclination_square(upper_theta) * sine * dtheta * dtheta);
      const double longitude_lower =
          0.5 * (pedersen[center] +
                 pedersen[polar_index(config, i, jm)]) /
          (sine * sine * dphi * dphi);
      const double longitude_upper =
          0.5 * (pedersen[center] +
                 pedersen[polar_index(config, i, jp)]) /
          (sine * sine * dphi * dphi);
      const double pedersen_term =
          lower * (potential_v[center] -
                   potential_v[polar_index(config, i - 1U, j)]) +
          upper * (potential_v[center] -
                   potential_v[polar_index(config, i + 1U, j)]) +
          longitude_lower *
              (potential_v[center] -
               potential_v[polar_index(config, i, jm)]) +
          longitude_upper *
              (potential_v[center] -
               potential_v[polar_index(config, i, jp)]);
      const double df12_dtheta =
          (tensor_f12(config, hall, i + 1U, j) -
           tensor_f12(config, hall, i - 1U, j)) /
          (2.0 * dtheta);
      const double df12_dphi =
          (tensor_f12(config, hall, i, jp) -
           tensor_f12(config, hall, i, jm)) /
          (2.0 * dphi);
      const double dpotential_dtheta =
          (potential_v[polar_index(config, i + 1U, j)] -
           potential_v[polar_index(config, i - 1U, j)]) /
          (2.0 * dtheta);
      const double dpotential_dphi =
          (potential_v[polar_index(config, i, jp)] -
           potential_v[polar_index(config, i, jm)]) /
          (2.0 * dphi);
      const double hall_term =
          (-df12_dtheta * dpotential_dphi +
           df12_dphi * dpotential_dtheta) /
          sine;
      result[center] = pedersen_term + hall_term;
    }
  }
  for (size_t j = 0U; j < np; ++j) {
    result[polar_index(config, nt - 1U, j)] =
        potential_v[polar_index(config, nt - 1U, j)];
  }
  return 0;
}

static double variable_pole_edge_coefficient(
    const gamera_mi_solver_config *config, const double *pedersen,
    size_t j, double pole_pedersen_mean) {
  const double dtheta = config->maximum_colatitude /
                        (double)(config->colatitude_count - 1U);
  const double theta_half = 0.5 * dtheta;
  const double face_pedersen =
      0.5 * (pole_pedersen_mean +
             pedersen[polar_index(config, 1U, j)]);
  return sin(theta_half) * face_pedersen /
         (inclination_square(theta_half) * dtheta * dtheta);
}

typedef struct {
  const gamera_mi_solver_config *config;
  double *diagonal;
  double *lower;
  double *upper;
  double *longitude_lower;
  double *longitude_upper;
  double *pole_edge;
  double *preconditioner_work;
  double *storage;
} gamera_mi_tensor_stencil;

static void destroy_tensor_stencil(gamera_mi_tensor_stencil *stencil) {
  if (stencil != NULL) {
    free(stencil->storage);
    *stencil = (gamera_mi_tensor_stencil){0};
  }
}

static int build_tensor_stencil(
    const gamera_mi_solver_config *config, const double *pedersen,
    const double *hall, double *diagonal,
    gamera_mi_tensor_stencil *stencil) {
  if (stencil == NULL || diagonal == NULL) {
    return -1;
  }
  const size_t np = config->longitude_count;
  const size_t nt = config->colatitude_count;
  const size_t count = 1U + np * (nt - 2U);
  if (count > (SIZE_MAX - np) / 5U) {
    return -1;
  }
  double *storage =
      (double *)calloc(5U * count + np, sizeof(*storage));
  if (storage == NULL) {
    return -1;
  }
  *stencil = (gamera_mi_tensor_stencil){
      .config = config,
      .diagonal = diagonal,
      .lower = storage,
      .upper = storage + count,
      .longitude_lower = storage + 2U * count,
      .longitude_upper = storage + 3U * count,
      .pole_edge = storage + 4U * count,
      .preconditioner_work = storage + 4U * count + np,
      .storage = storage};
  const double dtheta = config->maximum_colatitude / (double)(nt - 1U);
  const double dphi = 2.0 * acos(-1.0) / (double)np;
  double pole_pedersen_mean = 0.0;
  for (size_t j = 0U; j < np; ++j) {
    pole_pedersen_mean += pedersen[polar_index(config, 0U, j)] / (double)np;
  }
  diagonal[0] = 0.0;
  for (size_t j = 0U; j < np; ++j) {
    stencil->pole_edge[j] = variable_pole_edge_coefficient(
        config, pedersen, j, pole_pedersen_mean);
    diagonal[0] += stencil->pole_edge[j];
  }
  for (size_t i = 1U; i + 1U < nt; ++i) {
    const double theta = (double)i * dtheta;
    const double sine = sin(theta);
    const double lower_theta = theta - 0.5 * dtheta;
    const double upper_theta = theta + 0.5 * dtheta;
    for (size_t j = 0U; j < np; ++j) {
      const size_t jm = j == 0U ? np - 1U : j - 1U;
      const size_t jp = j + 1U == np ? 0U : j + 1U;
      const size_t full = polar_index(config, i, j);
      const double lower =
          sin(lower_theta) *
          0.5 * (pedersen[full] +
                 pedersen[polar_index(config, i - 1U, j)]) /
          (inclination_square(lower_theta) * dtheta * dtheta);
      const double upper =
          sin(upper_theta) *
          0.5 * (pedersen[full] +
                 pedersen[polar_index(config, i + 1U, j)]) /
          (inclination_square(upper_theta) * dtheta * dtheta);
      const double longitude_lower =
          0.5 * (pedersen[full] +
                 pedersen[polar_index(config, i, jm)]) /
          (sine * dphi * dphi);
      const double longitude_upper =
          0.5 * (pedersen[full] +
                 pedersen[polar_index(config, i, jp)]) /
          (sine * dphi * dphi);
      const size_t center = reduced_index(config, i, j);
      diagonal[center] =
          lower + upper + longitude_lower + longitude_upper;
      const double df12_dtheta =
          (tensor_f12(config, hall, i + 1U, j) -
           tensor_f12(config, hall, i - 1U, j)) /
          (2.0 * dtheta);
      const double df12_dphi =
          (tensor_f12(config, hall, i, jp) -
           tensor_f12(config, hall, i, jm)) /
          (2.0 * dphi);
      stencil->lower[center] =
          -lower - df12_dphi / (2.0 * dtheta);
      stencil->upper[center] =
          -upper + df12_dphi / (2.0 * dtheta);
      stencil->longitude_lower[center] =
          -longitude_lower + df12_dtheta / (2.0 * dphi);
      stencil->longitude_upper[center] =
          -longitude_upper - df12_dtheta / (2.0 * dphi);
    }
  }
  return 0;
}

static void apply_tensor_stencil(const gamera_mi_tensor_stencil *stencil,
                                 const double *potential_v,
                                 double *result) {
  const gamera_mi_solver_config *config = stencil->config;
  const size_t np = config->longitude_count;
  const size_t nt = config->colatitude_count;
  result[0] = stencil->diagonal[0] * potential_v[0];
  for (size_t j = 0U; j < np; ++j) {
    result[0] -= stencil->pole_edge[j] *
                 potential_v[reduced_index(config, 1U, j)];
  }
  for (size_t i = 1U; i + 1U < nt; ++i) {
    for (size_t j = 0U; j < np; ++j) {
      const size_t jm = j == 0U ? np - 1U : j - 1U;
      const size_t jp = j + 1U == np ? 0U : j + 1U;
      const size_t center = reduced_index(config, i, j);
      const double lower_value =
          i == 1U ? potential_v[0]
                  : potential_v[reduced_index(config, i - 1U, j)];
      const double upper_value =
          i + 2U < nt
              ? potential_v[reduced_index(config, i + 1U, j)]
              : 0.0;
      result[center] =
          stencil->diagonal[center] * potential_v[center] +
          stencil->lower[center] * lower_value +
          stencil->upper[center] * upper_value +
          stencil->longitude_lower[center] *
              potential_v[reduced_index(config, i, jm)] +
          stencil->longitude_upper[center] *
              potential_v[reduced_index(config, i, jp)];
    }
  }
}

/*
 * One symmetric Gauss-Seidel/SSOR application to the cached tensor matrix.
 * The former diagonal Jacobi preconditioner leaves the near-pole longitude
 * stiffness almost untouched, which makes the iteration count grow rapidly
 * with angular resolution.  This triangular preconditioner retains the
 * complete five-point Pedersen/Hall stencil, including the periodic seam and
 * the single reduced pole node, without changing the solved matrix.
 */
static void apply_tensor_preconditioner(gamera_mi_tensor_stencil *stencil,
                                        const double *input,
                                        double *output) {
  const gamera_mi_solver_config *config = stencil->config;
  const size_t np = config->longitude_count;
  const size_t nt = config->colatitude_count;
  double *forward = stencil->preconditioner_work;

  forward[0] = input[0] / stencil->diagonal[0];
  for (size_t i = 1U; i + 1U < nt; ++i) {
    for (size_t j = 0U; j < np; ++j) {
      const size_t center = reduced_index(config, i, j);
      double value = input[center];
      const size_t lower =
          i == 1U ? 0U : reduced_index(config, i - 1U, j);
      value -= stencil->lower[center] * forward[lower];
      if (j > 0U) {
        value -= stencil->longitude_lower[center] *
                 forward[reduced_index(config, i, j - 1U)];
      }
      if (j + 1U == np) {
        value -= stencil->longitude_upper[center] *
                 forward[reduced_index(config, i, 0U)];
      }
      forward[center] = value / stencil->diagonal[center];
    }
  }

  for (size_t i = nt - 2U; i > 0U; --i) {
    for (size_t j_reverse = np; j_reverse > 0U; --j_reverse) {
      const size_t j = j_reverse - 1U;
      const size_t center = reduced_index(config, i, j);
      double value = stencil->diagonal[center] * forward[center];
      if (i + 2U < nt) {
        value -= stencil->upper[center] *
                 output[reduced_index(config, i + 1U, j)];
      }
      if (j + 1U < np) {
        value -= stencil->longitude_upper[center] *
                 output[reduced_index(config, i, j + 1U)];
      }
      if (j == 0U) {
        value -= stencil->longitude_lower[center] *
                 output[reduced_index(config, i, np - 1U)];
      }
      output[center] = value / stencil->diagonal[center];
    }
  }
  double pole_value = stencil->diagonal[0] * forward[0];
  for (size_t j = 0U; j < np; ++j) {
    pole_value += stencil->pole_edge[j] *
                  output[reduced_index(config, 1U, j)];
  }
  output[0] = pole_value / stencil->diagonal[0];
}

int gamera_mi_solve_conductance_tensor(
    const gamera_mi_solver_config *config, const double *fac_a_m2,
    const double *pedersen, const double *hall, double *potential_v,
    gamera_mi_solver_stats *stats) {
  if (!valid_tensor_config(config, pedersen, hall) || fac_a_m2 == NULL ||
      potential_v == NULL) {
    return -1;
  }
  const size_t np = config->longitude_count;
  const size_t nt = config->colatitude_count;
  if (np > (SIZE_MAX - 1U) / (nt - 2U)) {
    return -1;
  }
  const size_t count = 1U + np * (nt - 2U);
  double *work = (double *)calloc(11U * count, sizeof(double));
  if (work == NULL) {
    return -1;
  }
  double *right = work;
  double *diagonal = right + count;
  double *solution = diagonal + count;
  double *residual = solution + count;
  double *shadow = residual + count;
  double *search = shadow + count;
  double *operator_search = search + count;
  double *preconditioned_search = operator_search + count;
  double *intermediate = preconditioned_search + count;
  double *preconditioned_intermediate = intermediate + count;
  double *operator_intermediate = preconditioned_intermediate + count;

  const double dtheta = config->maximum_colatitude / (double)(nt - 1U);
  const double radius_square =
      config->ionosphere_radius_m * config->ionosphere_radius_m;
  double right_norm_square = 0.0;
  for (size_t i = 1U; i + 1U < nt; ++i) {
    const double theta = (double)i * dtheta;
    const double weight = sin(theta);
    const double inclination =
        gamera_mi_dipole_cos_inclination(theta, config->hemisphere);
    for (size_t j = 0U; j < np; ++j) {
      const size_t full = polar_index(config, i, j);
      const size_t index = reduced_index(config, i, j);
      if (!isfinite(fac_a_m2[full])) {
        free(work);
        return -1;
      }
      const double unweighted = -radius_square * fac_a_m2[full] * inclination;
      right[index] = weight * unweighted;
      right_norm_square += unweighted * unweighted;
    }
  }
  double pole_mean = 0.0;
  for (size_t j = 0U; j < np; ++j) {
    const double value = potential_v[polar_index(config, 0U, j)];
    if (isfinite(value)) {
      pole_mean += value / (double)np;
    }
  }
  solution[0] = pole_mean;
  for (size_t i = 1U; i + 1U < nt; ++i) {
    for (size_t j = 0U; j < np; ++j) {
      const double value = potential_v[polar_index(config, i, j)];
      solution[reduced_index(config, i, j)] =
          isfinite(value) ? value : 0.0;
    }
  }
  gamera_mi_tensor_stencil stencil = {0};
  if (build_tensor_stencil(config, pedersen, hall, diagonal, &stencil) != 0) {
    free(work);
    return -1;
  }
  apply_tensor_stencil(&stencil, solution, operator_search);
  for (size_t index = 0U; index < count; ++index) {
    residual[index] = right[index] - operator_search[index];
    shadow[index] = residual[index];
  }
  const double initial = reduced_original_residual_norm(config, residual);
  const double target = fmax(
      config->absolute_tolerance,
      config->relative_tolerance *
          fmax(sqrt(fmax(0.0, right_norm_square)), 1.0));
  if (stats != NULL) {
    *stats = (gamera_mi_solver_stats){0, initial, initial, initial <= target};
  }
  double final = initial;
  int iterations = 0;
  int converged = initial <= target;
  double rho_previous = 1.0;
  double alpha = 1.0;
  double omega = 1.0;
  memset(search, 0, count * sizeof(*search));
  memset(operator_search, 0, count * sizeof(*operator_search));
  for (int iteration = 1; !converged && iteration <= config->maximum_iterations;
       ++iteration) {
    const double rho = array_dot(shadow, residual, count);
    if (!isfinite(rho) || fabs(rho) <= DBL_MIN ||
        !isfinite(omega) || fabs(omega) <= DBL_MIN) {
      break;
    }
    const double beta = (rho / rho_previous) * (alpha / omega);
    for (size_t index = 0U; index < count; ++index) {
      search[index] = residual[index] +
                      beta * (search[index] -
                              omega * operator_search[index]);
    }
    apply_tensor_preconditioner(&stencil, search, preconditioned_search);
    apply_tensor_stencil(&stencil, preconditioned_search, operator_search);
    const double alpha_denominator =
        array_dot(shadow, operator_search, count);
    if (!isfinite(alpha_denominator) ||
        fabs(alpha_denominator) <= DBL_MIN) {
      break;
    }
    alpha = rho / alpha_denominator;
    for (size_t index = 0U; index < count; ++index) {
      intermediate[index] = residual[index] -
                            alpha * operator_search[index];
    }
    final = reduced_original_residual_norm(config, intermediate);
    if (final <= target) {
      for (size_t index = 0U; index < count; ++index) {
        solution[index] += alpha * preconditioned_search[index];
      }
      iterations = iteration;
      converged = 1;
      break;
    }
    apply_tensor_preconditioner(&stencil, intermediate,
                                preconditioned_intermediate);
    apply_tensor_stencil(&stencil, preconditioned_intermediate,
                         operator_intermediate);
    const double omega_denominator =
        array_dot(operator_intermediate, operator_intermediate, count);
    if (!isfinite(omega_denominator) || omega_denominator <= DBL_MIN) {
      break;
    }
    omega = array_dot(operator_intermediate, intermediate, count) /
            omega_denominator;
    if (!isfinite(omega) || fabs(omega) <= DBL_MIN) {
      break;
    }
    for (size_t index = 0U; index < count; ++index) {
      solution[index] += alpha * preconditioned_search[index] +
                         omega * preconditioned_intermediate[index];
      residual[index] = intermediate[index] -
                        omega * operator_intermediate[index];
    }
    final = reduced_original_residual_norm(config, residual);
    iterations = iteration;
    if (final <= target || iteration % 256 == 0) {
      apply_tensor_stencil(&stencil, solution, operator_intermediate);
      for (size_t index = 0U; index < count; ++index) {
        residual[index] = right[index] - operator_intermediate[index];
      }
      final = reduced_original_residual_norm(config, residual);
      if (final <= target) {
        converged = 1;
        break;
      }
      memcpy(shadow, residual, count * sizeof(*shadow));
      memset(search, 0, count * sizeof(*search));
      memset(operator_search, 0, count * sizeof(*operator_search));
      rho_previous = 1.0;
      alpha = 1.0;
      omega = 1.0;
      continue;
    }
    rho_previous = rho;
  }
  for (size_t j = 0U; j < np; ++j) {
    potential_v[polar_index(config, 0U, j)] = solution[0];
    potential_v[polar_index(config, nt - 1U, j)] = 0.0;
  }
  for (size_t i = 1U; i + 1U < nt; ++i) {
    for (size_t j = 0U; j < np; ++j) {
      potential_v[polar_index(config, i, j)] =
          solution[reduced_index(config, i, j)];
    }
  }
  if (stats != NULL) {
    *stats = (gamera_mi_solver_stats){iterations, initial, final, converged};
  }
  destroy_tensor_stencil(&stencil);
  free(work);
  return converged ? 0 : -1;
}
