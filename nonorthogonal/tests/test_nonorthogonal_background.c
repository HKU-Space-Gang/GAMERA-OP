#include "nonorthogonal_advance.h"
#include "nonorthogonal_background.h"
#include "nonorthogonal_grid.h"
#include "nonorthogonal_mesh.h"
#include "nonorthogonal_state.h"
#include "nonorthogonal_storage.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures;

static double dot(gamera_no_vec3 left, gamera_no_vec3 right) {
  return left.value[0] * right.value[0] +
         left.value[1] * right.value[1] +
         left.value[2] * right.value[2];
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

static int constant_field(gamera_no_vec3 point, void *context,
                          gamera_no_vec3 *field) {
  (void)point;
  if (context == NULL || field == NULL) {
    return -1;
  }
  *field = *(const gamera_no_vec3 *)context;
  return 0;
}

static gamera_no_vec3 *allocate_vertices(const size_t cells[3]) {
  const size_t extent[3] = {cells[0] + 1U, cells[1] + 1U,
                            cells[2] + 1U};
  return (gamera_no_vec3 *)calloc(gamera_no_element_count3(extent),
                                  sizeof(gamera_no_vec3));
}

static void test_constant_background_quadrature(void) {
  const size_t cells[3] = {3U, 3U, 3U};
  const size_t vertices_extent[3] = {4U, 4U, 4U};
  gamera_no_vec3 *vertices = allocate_vertices(cells);
  gamera_no_grid grid = {0};
  gamera_no_background_data background = {0};
  const gamera_no_vec3 expected = {{0.37, -0.22, 0.41}};
  if (vertices == NULL) {
    ++failures;
    return;
  }
  for (size_t i = 0; i < vertices_extent[0]; ++i) {
    for (size_t j = 0; j < vertices_extent[1]; ++j) {
      for (size_t k = 0; k < vertices_extent[2]; ++k) {
        const double fi = (double)i;
        const double fj = (double)j;
        const double fk = (double)k;
        vertices[gamera_no_index3(vertices_extent, i, j, k)] =
            (gamera_no_vec3){{fi + 0.07 * fj * fk,
                              fj + 0.04 * fi * fk,
                              fk + 0.03 * fi * fj}};
      }
    }
  }
  if (gamera_no_grid_create(cells, vertices, &grid) != 0 ||
      gamera_no_background_create(&grid, constant_field, (void *)&expected,
                                  &background) != 0) {
    fprintf(stderr, "FAIL constant background setup\n");
    ++failures;
    free(vertices);
    gamera_no_background_destroy(&background);
    gamera_no_grid_destroy(&grid);
    return;
  }

  double maximum_force = 0.0;
  for (size_t cell = 0; cell < gamera_no_element_count3(cells); ++cell) {
    for (int component = 0; component < GAMERA_NO_DIM; ++component) {
      expect_near("constant cell average",
                  background.cell_magnetic[cell].value[component],
                  expected.value[component], 2.0e-14);
      maximum_force =
          fmax(maximum_force,
               fabs(background.cell_force[cell].value[component]));
    }
  }
  expect_near("constant pure-B0 force", maximum_force, 0.0, 8.0e-14);

  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t face_count =
        gamera_no_element_count3(grid.face[direction].extent);
    for (size_t face = 0; face < face_count; ++face) {
      for (int component = 0; component < GAMERA_NO_DIM; ++component) {
        expect_near("constant face average",
                    background.face_magnetic[direction][face]
                        .value[component],
                    expected.value[component], 2.0e-14);
      }
      expect_near("constant oriented face flux",
                  background.face_flux[direction][face],
                  dot(expected,
                      grid.face[direction].value[face].area_vector),
                  8.0e-14);
    }
    const size_t edge_count =
        gamera_no_element_count3(grid.edge[direction].extent);
    for (size_t edge = 0; edge < edge_count; ++edge) {
      for (int component = 0; component < GAMERA_NO_DIM; ++component) {
        expect_near("constant edge average",
                    background.edge_magnetic[direction][edge]
                        .value[component],
                    expected.value[component], 2.0e-14);
      }
    }
  }

  free(vertices);
  gamera_no_background_destroy(&background);
  gamera_no_grid_destroy(&grid);
}

static void test_dipole_divergence_and_equilibrium(void) {
  const size_t cells[3] = {9U, 9U, 9U};
  const double lower[3] = {2.0, 0.30 * acos(-1.0), -0.30 * acos(-1.0)};
  const double upper[3] = {5.0, 0.70 * acos(-1.0), 0.30 * acos(-1.0)};
  gamera_no_vec3 *vertices = allocate_vertices(cells);
  gamera_no_grid grid = {0};
  gamera_no_storage storage = {0};
  gamera_no_background_data background = {0};
  const gamera_no_dipole dipole = {{{0.0, 0.0, -1.0}}};
  if (vertices == NULL ||
      gamera_no_generate_spherical_vertices(cells, lower, upper, vertices) !=
          0 ||
      gamera_no_grid_create(cells, vertices, &grid) != 0 ||
      gamera_no_storage_create(cells, &storage) != 0 ||
      gamera_no_background_create(&grid, gamera_no_dipole_field,
                                  (void *)&dipole, &background) != 0) {
    fprintf(stderr, "FAIL dipole background setup\n");
    ++failures;
    free(vertices);
    gamera_no_background_destroy(&background);
    gamera_no_storage_destroy(&storage);
    gamera_no_grid_destroy(&grid);
    return;
  }

  double maximum_flux_divergence = 0.0;
  double maximum_force = 0.0;
  for (size_t i = 0; i < cells[0]; ++i) {
    for (size_t j = 0; j < cells[1]; ++j) {
      for (size_t k = 0; k < cells[2]; ++k) {
        double net_flux = 0.0;
        for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
          size_t upper_coordinate[3] = {i, j, k};
          ++upper_coordinate[direction];
          const size_t lower_face = gamera_no_index3(
              grid.face[direction].extent, i, j, k);
          const size_t upper_face = gamera_no_index3(
              grid.face[direction].extent, upper_coordinate[0],
              upper_coordinate[1], upper_coordinate[2]);
          net_flux += background.face_flux[direction][upper_face] -
                      background.face_flux[direction][lower_face];
        }
        maximum_flux_divergence = fmax(maximum_flux_divergence,
                                       fabs(net_flux));
        const size_t cell = gamera_no_index3(cells, i, j, k);
        maximum_force =
            fmax(maximum_force,
                 sqrt(dot(background.cell_force[cell],
                          background.cell_force[cell])));
      }
    }
  }
  expect_near("dipole integrated divB", maximum_flux_divergence, 0.0,
              3.0e-13);
  expect_near("dipole integrated force-free balance", maximum_force, 0.0,
              3.0e-12);

  const gamera_no_primitive primitive = {1.0, {{0.0, 0.0, 0.0}}, 0.8};
  double conserved[GAMERA_NO_FLUX_COUNT];
  (void)gamera_no_primitive_to_conserved(&primitive, 5.0 / 3.0, 1.0e-12,
                                         1.0e-12, conserved);
  const size_t cell_count = gamera_no_element_count3(cells);
  for (size_t cell = 0; cell < cell_count; ++cell) {
    memcpy(&storage.conserved[cell * GAMERA_NO_FLUX_COUNT], conserved,
           sizeof(conserved));
    memcpy(&storage.old_conserved[cell * GAMERA_NO_FLUX_COUNT], conserved,
           sizeof(conserved));
  }
  const size_t active_lower[3] = {4U, 4U, 4U};
  const size_t active_upper[3] = {5U, 5U, 5U};
  const gamera_no_advance_options options = {
      {5.0 / 3.0, 1.0e-12, 1.0e-12, 1.0, 0.15, 0.20, 4.0,
       true,        true,      true,      true},
      {1.0e-12, 1.0, 0.6, 4.0, 0.8, 0.001, true, true},
      {5.0 / 3.0, 1.0e-12, 1.0e-12, 4.0, true, true, true},
      NULL, NULL, NULL, NULL, 0.0};
  if (gamera_no_advance(&storage, &grid, active_lower, active_upper, 0.5,
                        0.001, &options, &background.field) != 0) {
    fprintf(stderr, "FAIL dipole background split-field advance\n");
    ++failures;
  } else {
    const size_t active = gamera_no_index3(cells, 4U, 4U, 4U);
    for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
      expect_near("dipole background equilibrium conserved",
                  storage.conserved[active * GAMERA_NO_FLUX_COUNT +
                                    (size_t)variable],
                  conserved[variable], 8.0e-12);
    }
    for (int component = 0; component < GAMERA_NO_DIM; ++component) {
      expect_near("dipole residual magnetic field",
                  storage.cell_magnetic[active].value[component], 0.0,
                  3.0e-13);
    }
  }

  free(vertices);
  gamera_no_background_destroy(&background);
  gamera_no_storage_destroy(&storage);
  gamera_no_grid_destroy(&grid);
}

int main(void) {
  test_constant_background_quadrature();
  test_dipole_divergence_and_equilibrium();
  if (failures != 0) {
    fprintf(stderr, "%d background-field check(s) failed\n", failures);
    return EXIT_FAILURE;
  }
  printf("all non-orthogonal background-field checks passed\n");
  return EXIT_SUCCESS;
}
