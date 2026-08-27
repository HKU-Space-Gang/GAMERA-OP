#include "nonorthogonal_grid.h"
#include "nonorthogonal_mesh.h"
#include "nonorthogonal_mi.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


static int failures;

static void expect_near(const char *name, double actual, double expected,
                        double tolerance) {
  if (!isfinite(actual) || fabs(actual - expected) > tolerance) {
    fprintf(stderr,
            "FAIL %s: actual=%.17g expected=%.17g tolerance=%.3e\n",
            name, actual, expected, tolerance);
    ++failures;
  }
}

static void test_finite_volume_current(void) {
  const size_t cells[3] = {7U, 7U, 7U};
  const size_t vertices_extent[3] = {8U, 8U, 8U};
  const double lower[3] = {-3.5, -3.5, -3.5};
  const double upper[3] = {3.5, 3.5, 3.5};
  gamera_no_vec3 *vertices = (gamera_no_vec3 *)calloc(
      gamera_no_element_count3(vertices_extent), sizeof(*vertices));
  gamera_no_vec3 *magnetic = (gamera_no_vec3 *)calloc(
      gamera_no_element_count3(cells), sizeof(*magnetic));
  gamera_no_grid grid = {0};
  if (vertices == NULL || magnetic == NULL ||
      gamera_no_generate_cartesian_vertices(cells, lower, upper, vertices) !=
          0 ||
      gamera_no_grid_create(cells, vertices, &grid) != 0) {
    fprintf(stderr, "FAIL finite-volume current setup\n");
    ++failures;
    free(vertices);
    free(magnetic);
    gamera_no_grid_destroy(&grid);
    return;
  }
  for (size_t i = 0; i < cells[0]; ++i) {
    for (size_t j = 0; j < cells[1]; ++j) {
      for (size_t k = 0; k < cells[2]; ++k) {
        const size_t cell = gamera_no_index3(cells, i, j, k);
        const gamera_no_vec3 point = grid.cell[cell].centroid;
        magnetic[cell] = (gamera_no_vec3){{
            -0.5 * point.value[1], 0.5 * point.value[0], 0.0}};
      }
    }
  }
  gamera_no_vec3 current;
  if (gamera_no_cell_current_from_residual(
          &grid, magnetic, 3U, 3U, 3U, &current) != 0) {
    fprintf(stderr, "FAIL finite-volume current returned an error\n");
    ++failures;
  } else {
    expect_near("solid-rotation curl x", current.value[0], 0.0, 2.0e-13);
    expect_near("solid-rotation curl y", current.value[1], 0.0, 2.0e-13);
    expect_near("solid-rotation curl z", current.value[2], 1.0, 2.0e-13);
  }
  free(vertices);
  free(magnetic);
  gamera_no_grid_destroy(&grid);
}

static void test_dipole_mapping(void) {
  const double pi = acos(-1.0);
  const double mapped =
      gamera_mi_mapped_colatitude(0.5 * pi, 4.0, 1.0);
  expect_near("4 RE equator maps to 60 MLAT", mapped, pi / 6.0,
              2.0e-15);
  expect_near("north inclination sign",
              gamera_mi_dipole_cos_inclination(mapped, GAMERA_MI_NORTH),
              -sqrt(3.0) / sqrt(3.25), 2.0e-15);
  expect_near("south inclination sign",
              gamera_mi_dipole_cos_inclination(mapped, GAMERA_MI_SOUTH),
              sqrt(3.0) / sqrt(3.25), 2.0e-15);
  expect_near("dipole field ratio",
              gamera_mi_dipole_field_ratio(0.5 * pi, 4.0, 1.0),
              64.0 * sqrt(3.25), 2.0e-13);
}

static void test_mi_ghost_velocity(void) {
  const gamera_no_vec3 source = {{1.0, 2.0, 3.0}};
  const gamera_no_vec3 background = {{0.0, 0.0, 2.0}};
  const gamera_no_vec3 drift = {{4.0, 5.0, 0.0}};
  gamera_no_vec3 ghost = {{0.0, 0.0, 0.0}};
  if (gamera_mi_compose_ghost_velocity(source, background, drift, &ghost) !=
      0) {
    fprintf(stderr, "FAIL M-I ghost velocity composition returned an error\n");
    ++failures;
    return;
  }
  expect_near("M-I ghost transverse x", ghost.value[0], 7.0, 1.0e-15);
  expect_near("M-I ghost transverse y", ghost.value[1], 8.0, 1.0e-15);
  expect_near("M-I ghost reflected parallel", ghost.value[2], -3.0,
              1.0e-15);
}

static void test_manufactured_potential(void) {
  /* Match the automatically derived M-I grid used by the 256x96x288
   * production run.  This directly exercises the pole reduction that avoids
   * the production-resolution BiCGStab stagnation. */
  const size_t np = 384U;
  const size_t nt = 39U;
  const double pi = acos(-1.0);
  const gamera_mi_solver_config config = {
      .longitude_count = np,
      .colatitude_count = nt,
      .maximum_colatitude = pi / 6.0,
      .ionosphere_radius_m = 6.371e6,
      .pedersen_siemens = 5.0,
      .hall_siemens = 0.0,
      .low_latitude_potential_v = 0.0,
      .hemisphere = GAMERA_MI_NORTH,
      .maximum_iterations = 10000,
      .relative_tolerance = 1.0e-11,
      .absolute_tolerance = 1.0e-10};
  const size_t count = np * nt;
  double *exact = (double *)calloc(count, sizeof(*exact));
  double *equation = (double *)calloc(count, sizeof(*equation));
  double *fac = (double *)calloc(count, sizeof(*fac));
  double *solution = (double *)calloc(count, sizeof(*solution));
  if (exact == NULL || equation == NULL || fac == NULL || solution == NULL) {
    fprintf(stderr, "FAIL manufactured potential allocation\n");
    ++failures;
    free(exact);
    free(equation);
    free(fac);
    free(solution);
    return;
  }
  const double dtheta = config.maximum_colatitude / (double)(nt - 1U);
  const double dphi = 2.0 * pi / (double)np;
  const double cap_sine_square =
      pow(sin(config.maximum_colatitude), 2.0);
  for (size_t i = 0; i < nt; ++i) {
    const double theta = (double)i * dtheta;
    const double radial = pow(sin(theta), 2.0) *
                          (1.0 - pow(sin(theta), 2.0) / cap_sine_square);
    for (size_t j = 0; j < np; ++j) {
      exact[i * np + j] = 30000.0 * radial *
                          (sin((double)j * dphi) +
                           0.25 * cos(2.0 * (double)j * dphi));
    }
  }
  if (gamera_mi_apply_constant_pedersen(&config, exact, equation) != 0) {
    fprintf(stderr, "FAIL manufactured operator application\n");
    ++failures;
  } else {
    const double radius_square = config.ionosphere_radius_m *
                                 config.ionosphere_radius_m;
    for (size_t i = 1U; i + 1U < nt; ++i) {
      const double inclination = gamera_mi_dipole_cos_inclination(
          (double)i * dtheta, config.hemisphere);
      for (size_t j = 0; j < np; ++j) {
        const size_t index = i * np + j;
        fac[index] = -equation[index] / (radius_square * inclination);
      }
    }
    gamera_mi_solver_stats stats;
    if (gamera_mi_solve_constant_pedersen(
            &config, fac, solution, &stats) != 0 || !stats.converged) {
      fprintf(stderr,
              "FAIL manufactured solver did not converge: iter=%d "
              "r0=%.6e r=%.6e\n",
              stats.iterations, stats.initial_residual,
              stats.final_residual);
      ++failures;
    } else {
      double maximum_error = 0.0;
      for (size_t index = 0; index < count; ++index) {
        maximum_error =
            fmax(maximum_error, fabs(solution[index] - exact[index]));
      }
      expect_near("manufactured potential maximum error", maximum_error,
                  0.0, 2.0e-5);
      if (stats.final_residual >= stats.initial_residual) {
        fprintf(stderr, "FAIL manufactured residual did not decrease\n");
        ++failures;
      }
      gamera_mi_solver_config south_config = config;
      south_config.hemisphere = GAMERA_MI_SOUTH;
      for (size_t index = 0; index < count; ++index) {
        fac[index] = -fac[index];
        solution[index] = 0.0;
      }
      if (gamera_mi_solve_constant_pedersen(
              &south_config, fac, solution, &stats) != 0 ||
          !stats.converged) {
        fprintf(stderr,
                "FAIL South manufactured solver did not converge: iter=%d "
                "r0=%.6e r=%.6e\n",
                stats.iterations, stats.initial_residual,
                stats.final_residual);
        ++failures;
      } else {
        double south_maximum_error = 0.0;
        for (size_t index = 0; index < count; ++index) {
          south_maximum_error =
              fmax(south_maximum_error,
                   fabs(solution[index] - exact[index]));
        }
        expect_near("North/South manufactured potential symmetry",
                    south_maximum_error, 0.0, 2.0e-5);
      }
    }
  }
  free(exact);
  free(equation);
  free(fac);
  free(solution);
}

static void test_variable_conductance_tensor(void) {
  const size_t np = 96U;
  const size_t nt = 25U;
  const size_t count = np * nt;
  const double pi = acos(-1.0);
  const gamera_mi_solver_config config = {
      .longitude_count = np,
      .colatitude_count = nt,
      .maximum_colatitude = pi / 6.0,
      .ionosphere_radius_m = 6.371e6,
      .pedersen_siemens = 5.0,
      .hall_siemens = 0.0,
      .low_latitude_potential_v = 0.0,
      .hemisphere = GAMERA_MI_NORTH,
      .maximum_iterations = 20000,
      .relative_tolerance = 1.0e-10,
      .absolute_tolerance = 1.0e-8};
  double *exact = (double *)calloc(count, sizeof(*exact));
  double *constant_result =
      (double *)calloc(count, sizeof(*constant_result));
  double *tensor_result = (double *)calloc(count, sizeof(*tensor_result));
  double *fac = (double *)calloc(count, sizeof(*fac));
  double *solution = (double *)calloc(count, sizeof(*solution));
  double *pedersen = (double *)calloc(count, sizeof(*pedersen));
  double *hall = (double *)calloc(count, sizeof(*hall));
  if (exact == NULL || constant_result == NULL || tensor_result == NULL ||
      fac == NULL || solution == NULL || pedersen == NULL || hall == NULL) {
    fprintf(stderr, "FAIL variable conductance allocation\n");
    ++failures;
    free(exact);
    free(constant_result);
    free(tensor_result);
    free(fac);
    free(solution);
    free(pedersen);
    free(hall);
    return;
  }
  const double dtheta = config.maximum_colatitude / (double)(nt - 1U);
  const double dphi = 2.0 * pi / (double)np;
  const double cap_sine_square = pow(sin(config.maximum_colatitude), 2.0);
  for (size_t i = 0U; i < nt; ++i) {
    const double theta = (double)i * dtheta;
    const double radial = pow(sin(theta), 2.0) *
                          (1.0 - pow(sin(theta), 2.0) /
                                     cap_sine_square);
    for (size_t j = 0U; j < np; ++j) {
      const double phi = (double)j * dphi;
      const size_t index = i * np + j;
      exact[index] = 25000.0 * radial *
                     (sin(phi) + 0.2 * cos(2.0 * phi));
      pedersen[index] = 5.0;
      hall[index] = 0.0;
    }
  }
  if (gamera_mi_apply_constant_pedersen(
          &config, exact, constant_result) != 0 ||
      gamera_mi_apply_conductance_tensor(
          &config, pedersen, hall, exact, tensor_result) != 0) {
    fprintf(stderr, "FAIL constant/tensor operator identity application\n");
    ++failures;
  } else {
    double maximum_difference = 0.0;
    for (size_t index = 0U; index < count; ++index) {
      maximum_difference = fmax(
          maximum_difference,
          fabs(constant_result[index] - tensor_result[index]));
    }
    expect_near("constant/tensor operator identity", maximum_difference,
                0.0, 2.0e-8);
  }
  for (size_t i = 0U; i < nt; ++i) {
    const double theta = (double)i * dtheta;
    for (size_t j = 0U; j < np; ++j) {
      const double phi = (double)j * dphi;
      const size_t index = i * np + j;
      pedersen[index] = 4.5 + 1.2 * pow(sin(theta), 2.0) *
                                  (1.0 + 0.25 * cos(phi));
      hall[index] = 1.0 + 0.6 * pow(sin(theta), 2.0) *
                             (1.0 - 0.3 * sin(phi));
    }
  }
  if (gamera_mi_apply_conductance_tensor(
          &config, pedersen, hall, exact, tensor_result) != 0) {
    fprintf(stderr, "FAIL variable tensor operator application\n");
    ++failures;
  } else {
    const double radius_square = config.ionosphere_radius_m *
                                 config.ionosphere_radius_m;
    for (size_t i = 1U; i + 1U < nt; ++i) {
      const double inclination = gamera_mi_dipole_cos_inclination(
          (double)i * dtheta, config.hemisphere);
      for (size_t j = 0U; j < np; ++j) {
        const size_t index = i * np + j;
        fac[index] = -tensor_result[index] /
                     (radius_square * inclination);
      }
    }
    gamera_mi_solver_stats stats;
    if (gamera_mi_solve_conductance_tensor(
            &config, fac, pedersen, hall, solution, &stats) != 0 ||
        !stats.converged) {
      fprintf(stderr,
              "FAIL variable tensor solver did not converge: iter=%d "
              "r0=%.6e r=%.6e\n",
              stats.iterations, stats.initial_residual,
              stats.final_residual);
      ++failures;
    } else {
      double maximum_error = 0.0;
      for (size_t index = 0U; index < count; ++index) {
        maximum_error =
            fmax(maximum_error, fabs(solution[index] - exact[index]));
      }
      expect_near("variable tensor manufactured solution", maximum_error,
                  0.0, 2.0e-3);
    }
  }
  /* Repeat with a production-like EUV background plus a sharp, three-degree
   * auroral oval transition and a large Hall/Pedersen ratio. */
  memset(solution, 0, count * sizeof(*solution));
  for (size_t i = 0U; i < nt; ++i) {
    const double latitude_deg =
        90.0 - 180.0 * (double)i * dtheta / pi;
    for (size_t j = 0U; j < np; ++j) {
      const double phi = (double)j * dphi;
      const size_t index = i * np + j;
      const double boundary_deg = 65.0 + 4.0 * cos(phi);
      const double coordinate =
          fmax(0.0, fmin(1.0, (latitude_deg - boundary_deg) / 3.0));
      const double oval_mask =
          coordinate * coordinate * (3.0 - 2.0 * coordinate);
      const double euv_pedersen = 2.0 + 12.0 * fmax(0.0, cos(phi));
      const double euv_hall = 0.8 * euv_pedersen;
      const double auroral_pedersen = 20.0 * oval_mask;
      const double auroral_hall = 3.0 * auroral_pedersen;
      pedersen[index] = hypot(euv_pedersen, auroral_pedersen);
      hall[index] = fmin(6.0 * pedersen[index],
                         hypot(euv_hall, auroral_hall));
    }
  }
  if (gamera_mi_apply_conductance_tensor(
          &config, pedersen, hall, exact, tensor_result) != 0) {
    fprintf(stderr, "FAIL sharp-oval tensor operator application\n");
    ++failures;
  } else {
    const double radius_square = config.ionosphere_radius_m *
                                 config.ionosphere_radius_m;
    memset(fac, 0, count * sizeof(*fac));
    for (size_t i = 1U; i + 1U < nt; ++i) {
      const double inclination = gamera_mi_dipole_cos_inclination(
          (double)i * dtheta, config.hemisphere);
      for (size_t j = 0U; j < np; ++j) {
        const size_t index = i * np + j;
        fac[index] = -tensor_result[index] /
                     (radius_square * inclination);
      }
    }
    gamera_mi_solver_stats stats;
    if (gamera_mi_solve_conductance_tensor(
            &config, fac, pedersen, hall, solution, &stats) != 0 ||
        !stats.converged) {
      fprintf(stderr,
              "FAIL sharp-oval tensor solver did not converge: iter=%d "
              "r0=%.6e r=%.6e\n",
              stats.iterations, stats.initial_residual,
              stats.final_residual);
      ++failures;
    } else {
      double maximum_error = 0.0;
      for (size_t index = 0U; index < count; ++index) {
        maximum_error =
            fmax(maximum_error, fabs(solution[index] - exact[index]));
      }
      expect_near("sharp-oval tensor manufactured solution", maximum_error,
                  0.0, 1.0e-1);
    }
  }
  free(exact);
  free(constant_result);
  free(tensor_result);
  free(fac);
  free(solution);
  free(pedersen);
  free(hall);
}

int main(void) {
  test_finite_volume_current();
  test_dipole_mapping();
  test_mi_ghost_velocity();
  test_manufactured_potential();
  test_variable_conductance_tensor();
  if (failures != 0) {
    fprintf(stderr, "%d M-I coupling kernel check(s) failed\n", failures);
    return EXIT_FAILURE;
  }
  printf("all M-I coupling kernel checks passed\n");
  return EXIT_SUCCESS;
}
