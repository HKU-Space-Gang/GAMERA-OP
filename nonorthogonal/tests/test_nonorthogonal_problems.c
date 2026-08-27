#include "nonorthogonal_initialization.h"
#include "nonorthogonal_mesh.h"
#include "nonorthogonal_state.h"
#include "nonorthogonal_step.h"
#include "nonorthogonal_yinyang.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures;
static void expect_near(const char *name, double actual, double expected,
                        double tolerance);

static void test_yinyang_coordinates(void) {
  const double pi = acos(-1.0);
  const double logical[][3] = {
      {1.0, 0.5 * pi, 0.0},
      {1.2, 0.25 * pi, -0.7 * pi},
      {0.7, 0.72 * pi, 0.63 * pi},
  };
  for (int patch = GAMERA_NO_YIN_PATCH;
       patch <= GAMERA_NO_YANG_PATCH; ++patch) {
    for (size_t point = 0; point < sizeof(logical) / sizeof(logical[0]);
         ++point) {
      gamera_no_vec3 global;
      double radius;
      double theta;
      double phi;
      if (gamera_no_yinyang_logical_to_global(
              patch, logical[point][0], logical[point][1], logical[point][2],
              &global) != 0 ||
          gamera_no_yinyang_global_to_logical(
              patch, global, &radius, &theta, &phi) != 0) {
        fprintf(stderr, "FAIL Yin-Yang coordinate round trip setup\n");
        ++failures;
        continue;
      }
      expect_near("Yin-Yang radius round trip", radius, logical[point][0],
                  3.0e-15);
      expect_near("Yin-Yang theta round trip", theta, logical[point][1],
                  3.0e-15);
      expect_near("Yin-Yang phi round trip", phi, logical[point][2],
                  3.0e-15);
    }
  }

  gamera_no_vec3 north;
  if (gamera_no_yinyang_logical_to_global(
          GAMERA_NO_YANG_PATCH, 1.0, 0.5 * pi, 0.5 * pi, &north) != 0) {
    fprintf(stderr, "FAIL Yang north-pole mapping\n");
    ++failures;
  } else {
    expect_near("Yang covers global north x", north.value[0], 0.0, 2.0e-15);
    expect_near("Yang covers global north y", north.value[1], 0.0, 2.0e-15);
    expect_near("Yang covers global north z", north.value[2], 1.0, 2.0e-15);
  }
}

static void expect_near(const char *name, double actual, double expected,
                        double tolerance) {
  if (!isfinite(actual) || fabs(actual - expected) > tolerance) {
    fprintf(stderr,
            "FAIL %s: actual=%.17g expected=%.17g tolerance=%.3e\n",
            name, actual, expected, tolerance);
    ++failures;
  }
}

static gamera_no_vec3 *allocate_vertices(const size_t cells[3]) {
  const size_t extent[3] = {cells[0] + 1U, cells[1] + 1U,
                            cells[2] + 1U};
  return (gamera_no_vec3 *)calloc(gamera_no_element_count3(extent),
                                  sizeof(gamera_no_vec3));
}

static void test_warp_conventions(void) {
  const size_t cells[3] = {4U, 4U, 1U};
  const size_t vertices_extent[3] = {5U, 5U, 2U};
  const double lower[3] = {0.0, 0.0, -0.5};
  const double upper[3] = {1.0, 1.0, 0.5};
  gamera_no_vec3 *vertices = allocate_vertices(cells);
  if (vertices == NULL) {
    ++failures;
    return;
  }
  if (gamera_no_generate_warped_cartesian_vertices(
          cells, lower, upper, 0.1, 1, GAMERA_NO_WARP_PAPER, vertices) != 0) {
    fprintf(stderr, "FAIL paper warped mesh generation\n");
    ++failures;
    free(vertices);
    return;
  }
  size_t center = gamera_no_index3(vertices_extent, 2U, 2U, 0U);
  expect_near("paper warp center x", vertices[center].value[0], 0.6,
              2.0e-15);
  expect_near("paper warp center y", vertices[center].value[1], 0.6,
              2.0e-15);
  expect_near("paper warp boundary x", vertices[0].value[0], 0.0,
              2.0e-15);

  if (gamera_no_generate_warped_cartesian_vertices(
          cells, lower, upper, 0.1, 1, GAMERA_NO_WARP_FORTRAN, vertices) !=
      0) {
    fprintf(stderr, "FAIL Fortran warped mesh generation\n");
    ++failures;
  } else {
    expect_near("Fortran warp center x", vertices[center].value[0], 0.6,
                2.0e-15);
    expect_near("Fortran warp center y", vertices[center].value[1], 0.4,
                2.0e-15);
  }
  free(vertices);
}

static void test_orszag_tang_presets(void) {
  const double pi = acos(-1.0);
  gamera_no_orszag_tang_options paper;
  gamera_no_orszag_tang_options fortran;
  gamera_no_orszag_tang_paper_defaults(&paper);
  gamera_no_orszag_tang_fortran_defaults(&fortran);
  expect_near("paper OT density", paper.density, 25.0 * pi / 36.0,
              2.0e-15);
  expect_near("paper OT pressure", paper.pressure, 5.0 * pi / 12.0,
              2.0e-15);
  expect_near("paper OT magnetic amplitude", paper.magnetic_amplitude, 1.0,
              0.0);
  expect_near("Fortran OT density", fortran.density, 25.0 / (36.0 * pi),
              2.0e-15);
  expect_near("Fortran OT pressure", fortran.pressure, 5.0 / (12.0 * pi),
              2.0e-15);
  expect_near("Fortran OT magnetic amplitude", fortran.magnetic_amplitude,
              1.0 / sqrt(4.0 * pi), 2.0e-15);
}

static void test_warped_orszag_tang(void) {
  const size_t cells[3] = {16U, 16U, 1U};
  const double lower[3] = {0.0, 0.0, -0.5};
  const double upper[3] = {1.0, 1.0, 0.5};
  gamera_no_vec3 *vertices = allocate_vertices(cells);
  gamera_no_grid grid = {0};
  gamera_no_storage storage = {0};
  if (vertices == NULL ||
      gamera_no_generate_warped_cartesian_vertices(
          cells, lower, upper, 0.1, 1, GAMERA_NO_WARP_FORTRAN, vertices) !=
          0 ||
      gamera_no_grid_create(cells, vertices, &grid) != 0 ||
      gamera_no_storage_create(cells, &storage) != 0) {
    fprintf(stderr, "FAIL warped OT setup\n");
    ++failures;
    free(vertices);
    gamera_no_grid_destroy(&grid);
    gamera_no_storage_destroy(&storage);
    return;
  }
  gamera_no_orszag_tang_options options;
  gamera_no_orszag_tang_fortran_defaults(&options);
  if (gamera_no_initialize_orszag_tang(&grid, &storage, &options, 5.0 / 3.0,
                                       1.0e-12, 1.0e-12) != 0) {
    fprintf(stderr, "FAIL warped OT initialization\n");
    ++failures;
  } else {
    double max_net_flux = 0.0;
    for (size_t i = 0; i < cells[0]; ++i) {
      for (size_t j = 0; j < cells[1]; ++j) {
        const double net = fabs(gamera_no_cell_net_flux(&storage, i, j, 0U));
        max_net_flux = fmax(max_net_flux, net);
      }
    }
    expect_near("warped OT discrete divB", max_net_flux, 0.0, 2.0e-15);

    const size_t i = 5U;
    const size_t j = 7U;
    const size_t cell = gamera_no_index3(cells, i, j, 0U);
    const double fortran_gas_gold[GAMERA_NO_FLUX_COUNT] = {
        2.21048532072076864e-1, -1.49444576078564928e-1,
        9.32883522910203239e-2, 0.0, 2.69146357212345144e-1};
    for (size_t variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
      expect_near("warped OT Fortran gas gold",
                  storage.conserved[cell * GAMERA_NO_FLUX_COUNT + variable],
                  fortran_gas_gold[variable], 8.0e-14);
    }
    gamera_no_primitive primitive;
    if (gamera_no_conserved_to_primitive(
            &storage.conserved[cell * GAMERA_NO_FLUX_COUNT], 5.0 / 3.0,
            1.0e-12, 1.0e-12, &primitive) != 0) {
      fprintf(stderr, "FAIL warped OT conserved-to-primitive\n");
      ++failures;
    } else {
      const gamera_no_vec3 point = grid.cell[cell].centroid;
      const double pi = acos(-1.0);
      expect_near("warped OT density", primitive.density, options.density,
                  2.0e-15);
      expect_near("warped OT pressure", primitive.pressure, options.pressure,
                  2.0e-15);
      expect_near("warped OT vx", primitive.velocity.value[0],
                  -sin(2.0 * pi * point.value[1]), 3.0e-15);
      expect_near("warped OT vy", primitive.velocity.value[1],
                  sin(2.0 * pi * point.value[0]), 3.0e-15);
    }
  }
  free(vertices);
  gamera_no_grid_destroy(&grid);
  gamera_no_storage_destroy(&storage);
}

static void test_spherical_blast(void) {
  const size_t cells[3] = {8U, 8U, 8U};
  const double pi = acos(-1.0);
  const double lower[3] = {0.5, 0.3 * pi, -0.2 * pi};
  const double upper[3] = {1.5, 0.7 * pi, 0.2 * pi};
  gamera_no_vec3 *vertices = allocate_vertices(cells);
  gamera_no_grid grid = {0};
  gamera_no_storage storage = {0};
  if (vertices == NULL ||
      gamera_no_generate_spherical_vertices(cells, lower, upper, vertices) !=
          0 ||
      gamera_no_grid_create(cells, vertices, &grid) != 0 ||
      gamera_no_storage_create(cells, &storage) != 0) {
    fprintf(stderr, "FAIL spherical blast setup\n");
    ++failures;
    free(vertices);
    gamera_no_grid_destroy(&grid);
    gamera_no_storage_destroy(&storage);
    return;
  }
  gamera_no_blast_options options;
  gamera_no_spherical_blast_defaults(&options);
  if (gamera_no_initialize_blast(&grid, &storage, &options, 5.0 / 3.0,
                                 1.0e-12, 1.0e-12) != 0) {
    fprintf(stderr, "FAIL spherical blast initialization\n");
    ++failures;
  } else {
    double max_net_flux = 0.0;
    double maximum_pressure = 0.0;
    double minimum_pressure = INFINITY;
    double maximum_field_error = 0.0;
    for (size_t i = 0; i < cells[0]; ++i) {
      for (size_t j = 0; j < cells[1]; ++j) {
        for (size_t k = 0; k < cells[2]; ++k) {
          const size_t cell = gamera_no_index3(cells, i, j, k);
          gamera_no_primitive primitive;
          if (gamera_no_conserved_to_primitive(
                  &storage.conserved[cell * GAMERA_NO_FLUX_COUNT],
                  5.0 / 3.0, 1.0e-12, 1.0e-12, &primitive) != 0) {
            ++failures;
            continue;
          }
          maximum_pressure = fmax(maximum_pressure, primitive.pressure);
          minimum_pressure = fmin(minimum_pressure, primitive.pressure);
          max_net_flux =
              fmax(max_net_flux,
                   fabs(gamera_no_cell_net_flux(&storage, i, j, k)));
          for (int d = 0; d < GAMERA_NO_DIM; ++d) {
            maximum_field_error =
                fmax(maximum_field_error,
                     fabs(storage.cell_magnetic[cell].value[d] -
                          options.magnetic.value[d]));
          }
        }
      }
    }
    expect_near("spherical blast ambient pressure", minimum_pressure, 0.1,
                3.0e-14);
    if (!(maximum_pressure > 0.1 && maximum_pressure <= 10.0)) {
      fprintf(stderr, "FAIL spherical blast pressure range: %.17g\n",
              maximum_pressure);
      ++failures;
    }
    expect_near("spherical blast discrete divB", max_net_flux, 0.0,
                2.0e-15);
    expect_near("spherical blast uniform B recovery", maximum_field_error,
                0.0, 2.0e-12);

    const size_t inner_cell_a = gamera_no_index3(cells, 3U, 3U, 3U);
    const size_t inner_cell_b = gamera_no_index3(cells, 4U, 3U, 3U);
    const double pressure_a =
        (5.0 / 3.0 - 1.0) *
        storage.conserved[inner_cell_a * GAMERA_NO_FLUX_COUNT +
                          GAMERA_NO_FLUX_ENERGY];
    const double pressure_b =
        (5.0 / 3.0 - 1.0) *
        storage.conserved[inner_cell_b * GAMERA_NO_FLUX_COUNT +
                          GAMERA_NO_FLUX_ENERGY];
    expect_near("spherical blast Fortran volume-average gold A", pressure_a,
                1.82086481422018154, 2.0e-13);
    expect_near("spherical blast Fortran volume-average gold B", pressure_b,
                1.78292014323528525, 2.0e-13);
  }
  free(vertices);
  gamera_no_grid_destroy(&grid);
  gamera_no_storage_destroy(&storage);
}

int main(void) {
  test_yinyang_coordinates();
  test_warp_conventions();
  test_orszag_tang_presets();
  test_warped_orszag_tang();
  test_spherical_blast();
  if (failures != 0) {
    fprintf(stderr, "%d non-orthogonal problem test(s) failed\n", failures);
    return EXIT_FAILURE;
  }
  printf("all non-orthogonal problem tests passed\n");
  return EXIT_SUCCESS;
}
