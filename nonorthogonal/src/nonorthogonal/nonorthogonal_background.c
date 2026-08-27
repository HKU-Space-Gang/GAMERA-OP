#include "nonorthogonal_background.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

enum { GAMERA_NO_BACKGROUND_GAUSS_POINTS = 12 };

static const double gauss_node[GAMERA_NO_BACKGROUND_GAUSS_POINTS] = {
    0.1252334085114689154724,  0.3678314989981801937527,
    0.5873179542866174472967,  0.7699026741943046870369,
    0.9041172563704748566785,  0.9815606342467192506906,
    -0.1252334085114689154724, -0.3678314989981801937527,
    -0.5873179542866174472967, -0.7699026741943046870369,
    -0.9041172563704748566785, -0.9815606342467192506906};

static const double gauss_weight[GAMERA_NO_BACKGROUND_GAUSS_POINTS] = {
    0.2491470458134027850006, 0.2334925365383548087610,
    0.2031674267230659217490, 0.1600783285433462263350,
    0.1069393259953184309603, 0.0471753363865118271946,
    0.2491470458134027850006, 0.2334925365383548087610,
    0.2031674267230659217490, 0.1600783285433462263350,
    0.1069393259953184309603, 0.0471753363865118271946};

static const int logical_sign[8][GAMERA_NO_DIM] = {
    {-1, -1, -1}, {+1, -1, -1}, {+1, +1, -1}, {-1, +1, -1},
    {-1, -1, +1}, {+1, -1, +1}, {+1, +1, +1}, {-1, +1, +1}};

static const size_t corner_offset[8][GAMERA_NO_DIM] = {
    {0U, 0U, 0U}, {1U, 0U, 0U}, {1U, 1U, 0U}, {0U, 1U, 0U},
    {0U, 0U, 1U}, {1U, 0U, 1U}, {1U, 1U, 1U}, {0U, 1U, 1U}};

/* Same positive-logical-normal ordering as nonorthogonal_geometry.c. */
static const int face_corner[GAMERA_NO_DIM][2][4] = {
    {{0, 3, 4, 7}, {1, 2, 5, 6}},
    {{0, 4, 1, 5}, {3, 7, 2, 6}},
    {{0, 1, 3, 2}, {4, 5, 7, 6}}};

static gamera_no_vec3 zero_vector(void) {
  return (gamera_no_vec3){{0.0, 0.0, 0.0}};
}

static gamera_no_vec3 add(gamera_no_vec3 left, gamera_no_vec3 right) {
  gamera_no_vec3 result;
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    result.value[component] =
        left.value[component] + right.value[component];
  }
  return result;
}

static gamera_no_vec3 subtract(gamera_no_vec3 left,
                               gamera_no_vec3 right) {
  gamera_no_vec3 result;
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    result.value[component] =
        left.value[component] - right.value[component];
  }
  return result;
}

static gamera_no_vec3 scale(gamera_no_vec3 vector, double factor) {
  gamera_no_vec3 result;
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    result.value[component] = factor * vector.value[component];
  }
  return result;
}

static double dot(gamera_no_vec3 left, gamera_no_vec3 right) {
  double result = 0.0;
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    result += left.value[component] * right.value[component];
  }
  return result;
}

static gamera_no_vec3 cross(gamera_no_vec3 left, gamera_no_vec3 right) {
  return (gamera_no_vec3){{
      left.value[1] * right.value[2] - left.value[2] * right.value[1],
      left.value[2] * right.value[0] - left.value[0] * right.value[2],
      left.value[0] * right.value[1] - left.value[1] * right.value[0]}};
}

static double norm(gamera_no_vec3 vector) { return sqrt(dot(vector, vector)); }

static int finite_vector(gamera_no_vec3 vector) {
  return isfinite(vector.value[0]) && isfinite(vector.value[1]) &&
         isfinite(vector.value[2]);
}

static gamera_no_vec3 vertex(const gamera_no_grid *grid, size_t i, size_t j,
                             size_t k) {
  return grid->vertex[gamera_no_index3(grid->vertex_extent, i, j, k)];
}

static void get_cell_corners(const gamera_no_grid *grid, size_t i, size_t j,
                             size_t k, gamera_no_vec3 corners[8]) {
  for (int corner = 0; corner < 8; ++corner) {
    corners[corner] =
        vertex(grid, i + corner_offset[corner][0],
               j + corner_offset[corner][1], k + corner_offset[corner][2]);
  }
}

static gamera_no_vec3 bilinear_point(const gamera_no_vec3 face[4],
                                     double eta, double psi) {
  const gamera_no_vec3 d_eta = subtract(face[1], face[0]);
  const gamera_no_vec3 d_psi = subtract(face[2], face[0]);
  const gamera_no_vec3 d_cross =
      add(subtract(face[3], face[2]), subtract(face[0], face[1]));
  return add(face[0],
             add(scale(d_eta, eta),
                 add(scale(d_psi, psi), scale(d_cross, eta * psi))));
}

static int integrate_face(const gamera_no_vec3 face[4], double scalar_area,
                          gamera_no_background_function function,
                          void *context, gamera_no_vec3 *average,
                          double *flux, gamera_no_vec3 *stress) {
  const gamera_no_vec3 d_eta = subtract(face[1], face[0]);
  const gamera_no_vec3 d_psi = subtract(face[2], face[0]);
  const gamera_no_vec3 d_cross =
      add(subtract(face[3], face[2]), subtract(face[0], face[1]));
  gamera_no_vec3 integral = zero_vector();
  gamera_no_vec3 stress_integral = zero_vector();
  double flux_integral = 0.0;
  if (!(scalar_area > DBL_MIN) || !isfinite(scalar_area)) {
    return -1;
  }

  for (int i = 0; i < GAMERA_NO_BACKGROUND_GAUSS_POINTS; ++i) {
    const double eta = 0.5 * (1.0 + gauss_node[i]);
    for (int j = 0; j < GAMERA_NO_BACKGROUND_GAUSS_POINTS; ++j) {
      const double psi = 0.5 * (1.0 + gauss_node[j]);
      const double weight = 0.25 * gauss_weight[i] * gauss_weight[j];
      const gamera_no_vec3 tangent_eta =
          add(d_eta, scale(d_cross, psi));
      const gamera_no_vec3 tangent_psi =
          add(d_psi, scale(d_cross, eta));
      const gamera_no_vec3 oriented_area =
          cross(tangent_eta, tangent_psi);
      const double area = norm(oriented_area);
      const gamera_no_vec3 point = bilinear_point(face, eta, psi);
      gamera_no_vec3 field;
      if (!(area > DBL_MIN) || !isfinite(area) ||
          function(point, context, &field) != 0 || !finite_vector(field)) {
        return -1;
      }
      integral = add(integral, scale(field, weight * area));
      const double normal_field = dot(field, oriented_area);
      flux_integral += weight * normal_field;
      const double magnetic_pressure = 0.5 * dot(field, field);
      stress_integral =
          add(stress_integral,
              scale(subtract(scale(oriented_area, magnetic_pressure),
                             scale(field, normal_field)),
                    weight));
    }
  }
  *average = scale(integral, 1.0 / scalar_area);
  *flux = flux_integral;
  *stress = stress_integral;
  return finite_vector(*average) && isfinite(*flux) && finite_vector(*stress)
             ? 0
             : -1;
}

static int integrate_cell(const gamera_no_vec3 corners[8], double volume,
                          gamera_no_background_function function,
                          void *context, gamera_no_vec3 *average) {
  gamera_no_vec3 integral = zero_vector();
  double integrated_volume = 0.0;
  if (!(volume > DBL_MIN) || !isfinite(volume)) {
    return -1;
  }
  for (int i = 0; i < GAMERA_NO_BACKGROUND_GAUSS_POINTS; ++i) {
    const double xi = gauss_node[i];
    for (int j = 0; j < GAMERA_NO_BACKGROUND_GAUSS_POINTS; ++j) {
      const double eta = gauss_node[j];
      for (int k = 0; k < GAMERA_NO_BACKGROUND_GAUSS_POINTS; ++k) {
        const double zeta = gauss_node[k];
        gamera_no_vec3 point = zero_vector();
        gamera_no_vec3 derivative[3] = {zero_vector(), zero_vector(),
                                        zero_vector()};
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
          for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
            derivative[direction] =
                add(derivative[direction],
                    scale(corners[corner], dshape[direction]));
          }
        }
        const double jacobian =
            fabs(dot(derivative[0], cross(derivative[1], derivative[2])));
        const double weight = jacobian * gauss_weight[i] * gauss_weight[j] *
                              gauss_weight[k];
        gamera_no_vec3 field;
        if (!isfinite(weight) || function(point, context, &field) != 0 ||
            !finite_vector(field)) {
          return -1;
        }
        integrated_volume += weight;
        integral = add(integral, scale(field, weight));
      }
    }
  }
  if (!(integrated_volume > DBL_MIN) ||
      fabs(integrated_volume - volume) > 2.0e-11 * volume) {
    return -1;
  }
  *average = scale(integral, 1.0 / volume);
  return finite_vector(*average) ? 0 : -1;
}

static int integrate_edge(gamera_no_vec3 start, gamera_no_vec3 end,
                          gamera_no_background_function function,
                          void *context, gamera_no_vec3 *average) {
  const gamera_no_vec3 displacement = subtract(end, start);
  gamera_no_vec3 result = zero_vector();
  if (!(norm(displacement) > DBL_MIN)) {
    return -1;
  }
  for (int n = 0; n < GAMERA_NO_BACKGROUND_GAUSS_POINTS; ++n) {
    const double parameter = 0.5 * (1.0 + gauss_node[n]);
    gamera_no_vec3 field;
    if (function(add(start, scale(displacement, parameter)), context,
                 &field) != 0 ||
        !finite_vector(field)) {
      return -1;
    }
    result = add(result, scale(field, 0.5 * gauss_weight[n]));
  }
  *average = result;
  return finite_vector(*average) ? 0 : -1;
}

static int allocate_items(size_t count, size_t size, void **result) {
  if (count == 0U || size == 0U || count > SIZE_MAX / size ||
      result == NULL) {
    return -1;
  }
  *result = calloc(count, size);
  return *result == NULL ? -1 : 0;
}

void gamera_no_background_destroy(gamera_no_background_data *background) {
  if (background == NULL) {
    return;
  }
  free(background->cell_magnetic);
  free(background->cell_force);
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    free(background->face_magnetic[direction]);
    free(background->face_flux[direction]);
    free(background->edge_magnetic[direction]);
  }
  memset(background, 0, sizeof(*background));
}

int gamera_no_background_create(const gamera_no_grid *grid,
                                gamera_no_background_function function,
                                void *context,
                                gamera_no_background_data *background) {
  if (grid == NULL || function == NULL || background == NULL ||
      grid->vertex == NULL || grid->cell == NULL) {
    return -1;
  }
  memset(background, 0, sizeof(*background));
  const size_t cell_count = gamera_no_element_count3(grid->cell_extent);
  gamera_no_vec3 *face_stress[GAMERA_NO_DIM] = {NULL, NULL, NULL};
  if (allocate_items(cell_count, sizeof(*background->cell_magnetic),
                     (void **)&background->cell_magnetic) != 0 ||
      allocate_items(cell_count, sizeof(*background->cell_force),
                     (void **)&background->cell_force) != 0) {
    gamera_no_background_destroy(background);
    return -1;
  }
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t face_count =
        gamera_no_element_count3(grid->face[direction].extent);
    const size_t edge_count =
        gamera_no_element_count3(grid->edge[direction].extent);
    if (allocate_items(face_count, sizeof(*background->face_magnetic[direction]),
                       (void **)&background->face_magnetic[direction]) != 0 ||
        allocate_items(face_count, sizeof(*background->face_flux[direction]),
                       (void **)&background->face_flux[direction]) != 0 ||
        allocate_items(edge_count, sizeof(*background->edge_magnetic[direction]),
                       (void **)&background->edge_magnetic[direction]) != 0 ||
        allocate_items(face_count, sizeof(*face_stress[direction]),
                       (void **)&face_stress[direction]) != 0) {
      for (int d = 0; d < GAMERA_NO_DIM; ++d) {
        free(face_stress[d]);
      }
      gamera_no_background_destroy(background);
      return -1;
    }
  }

  int failed = 0;
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
#pragma omp parallel for collapse(3) reduction(| : failed) schedule(static)
    for (size_t i = 0; i < grid->face[direction].extent[0]; ++i) {
      for (size_t j = 0; j < grid->face[direction].extent[1]; ++j) {
        for (size_t k = 0; k < grid->face[direction].extent[2]; ++k) {
          size_t cell_coordinate[3] = {i, j, k};
          int side = GAMERA_NO_LOWER;
          if (cell_coordinate[direction] == grid->cell_extent[direction]) {
            --cell_coordinate[direction];
            side = GAMERA_NO_UPPER;
          }
          gamera_no_vec3 corners[8];
          gamera_no_vec3 face[4];
          get_cell_corners(grid, cell_coordinate[0], cell_coordinate[1],
                           cell_coordinate[2], corners);
          for (int corner = 0; corner < 4; ++corner) {
            face[corner] = corners[face_corner[direction][side][corner]];
          }
          const size_t index = gamera_no_index3(
              grid->face[direction].extent, i, j, k);
          if (integrate_face(face, grid->face[direction].value[index].area,
                             function, context,
                             &background->face_magnetic[direction][index],
                             &background->face_flux[direction][index],
                             &face_stress[direction][index]) != 0) {
            failed = 1;
          }
        }
      }
    }
  }

#pragma omp parallel for collapse(3) reduction(| : failed) schedule(static)
  for (size_t i = 0; i < grid->cell_extent[0]; ++i) {
    for (size_t j = 0; j < grid->cell_extent[1]; ++j) {
      for (size_t k = 0; k < grid->cell_extent[2]; ++k) {
        gamera_no_vec3 corners[8];
        get_cell_corners(grid, i, j, k, corners);
        const size_t cell =
            gamera_no_index3(grid->cell_extent, i, j, k);
        if (integrate_cell(corners, grid->cell[cell].volume, function,
                           context, &background->cell_magnetic[cell]) != 0) {
          failed = 1;
          continue;
        }
        gamera_no_vec3 upper_minus_lower = zero_vector();
        for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
          size_t lower_coordinate[3] = {i, j, k};
          size_t upper_coordinate[3] = {i, j, k};
          ++upper_coordinate[direction];
          const size_t lower_face = gamera_no_index3(
              grid->face[direction].extent, lower_coordinate[0],
              lower_coordinate[1], lower_coordinate[2]);
          const size_t upper_face = gamera_no_index3(
              grid->face[direction].extent, upper_coordinate[0],
              upper_coordinate[1], upper_coordinate[2]);
          upper_minus_lower =
              add(upper_minus_lower,
                  subtract(face_stress[direction][upper_face],
                           face_stress[direction][lower_face]));
        }
        background->cell_force[cell] =
            scale(upper_minus_lower, -1.0 / grid->cell[cell].volume);
        if (!finite_vector(background->cell_force[cell])) {
          failed = 1;
        }
      }
    }
  }

  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
#pragma omp parallel for collapse(3) reduction(| : failed) schedule(static)
    for (size_t i = 0; i < grid->edge[direction].extent[0]; ++i) {
      for (size_t j = 0; j < grid->edge[direction].extent[1]; ++j) {
        for (size_t k = 0; k < grid->edge[direction].extent[2]; ++k) {
          size_t end_coordinate[3] = {i, j, k};
          ++end_coordinate[direction];
          const size_t edge = gamera_no_index3(
              grid->edge[direction].extent, i, j, k);
          if (integrate_edge(
                  vertex(grid, i, j, k),
                  vertex(grid, end_coordinate[0], end_coordinate[1],
                         end_coordinate[2]),
                  function, context,
                  &background->edge_magnetic[direction][edge]) != 0) {
            failed = 1;
          }
        }
      }
    }
  }

  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    free(face_stress[direction]);
    background->field.face_magnetic[direction] =
        background->face_magnetic[direction];
    background->field.face_flux[direction] = background->face_flux[direction];
    background->field.edge_magnetic[direction] =
        background->edge_magnetic[direction];
  }
  background->field.cell_magnetic = background->cell_magnetic;
  background->field.cell_force = background->cell_force;
  if (failed) {
    gamera_no_background_destroy(background);
    return -1;
  }
  return 0;
}

int gamera_no_dipole_field(gamera_no_vec3 point, void *context,
                           gamera_no_vec3 *field) {
  const gamera_no_dipole *dipole = (const gamera_no_dipole *)context;
  if (dipole == NULL || field == NULL || !finite_vector(point) ||
      !finite_vector(dipole->moment)) {
    return -1;
  }
  const double radius_squared = dot(point, point);
  if (!(radius_squared > DBL_MIN) || !isfinite(radius_squared)) {
    return -1;
  }
  const double radius = sqrt(radius_squared);
  const double inverse_radius_cubed =
      1.0 / (radius_squared * radius);
  const double inverse_radius_fifth =
      inverse_radius_cubed / radius_squared;
  const double projection = dot(dipole->moment, point);
  *field = subtract(scale(point, 3.0 * projection * inverse_radius_fifth),
                    scale(dipole->moment, inverse_radius_cubed));
  return finite_vector(*field) ? 0 : -1;
}
