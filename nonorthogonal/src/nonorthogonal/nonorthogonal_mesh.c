#include "nonorthogonal_mesh.h"

#include "nonorthogonal_grid.h"

#include <math.h>
#include <stddef.h>

static int valid_mesh_arguments(const size_t cell_extent[3],
                                const double lower[3],
                                const double upper[3],
                                const gamera_no_vec3 *vertices) {
  if (cell_extent == NULL || lower == NULL || upper == NULL ||
      vertices == NULL) {
    return 0;
  }
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    if (cell_extent[direction] == 0 || !isfinite(lower[direction]) ||
        !isfinite(upper[direction]) ||
        upper[direction] <= lower[direction]) {
      return 0;
    }
  }
  return 1;
}

static void vertex_extents(const size_t cell_extent[3], size_t extent[3]) {
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    extent[direction] = cell_extent[direction] + 1U;
  }
}

int gamera_no_generate_cartesian_vertices(
    const size_t cell_extent[3], const double lower[3], const double upper[3],
    gamera_no_vec3 *vertices) {
  if (!valid_mesh_arguments(cell_extent, lower, upper, vertices)) {
    return -1;
  }
  size_t extent[3];
  vertex_extents(cell_extent, extent);
  for (size_t i = 0; i < extent[0]; ++i) {
    const double u = (double)i / (double)cell_extent[0];
    for (size_t j = 0; j < extent[1]; ++j) {
      const double v = (double)j / (double)cell_extent[1];
      for (size_t k = 0; k < extent[2]; ++k) {
        const double w = (double)k / (double)cell_extent[2];
        const size_t index = gamera_no_index3(extent, i, j, k);
        vertices[index].value[0] =
            lower[0] + (upper[0] - lower[0]) * u;
        vertices[index].value[1] =
            lower[1] + (upper[1] - lower[1]) * v;
        vertices[index].value[2] =
            lower[2] + (upper[2] - lower[2]) * w;
      }
    }
  }
  return 0;
}

int gamera_no_generate_warped_cartesian_vertices(
    const size_t cell_extent[3], const double lower[3], const double upper[3],
    double w0, int wave_number, gamera_no_warp_convention convention,
    gamera_no_vec3 *vertices) {
  if (!valid_mesh_arguments(cell_extent, lower, upper, vertices) ||
      !isfinite(w0) || w0 < 0.0 || wave_number <= 0 ||
      (convention != GAMERA_NO_WARP_PAPER &&
       convention != GAMERA_NO_WARP_FORTRAN)) {
    return -1;
  }
  const double pi = acos(-1.0);
  size_t extent[3];
  vertex_extents(cell_extent, extent);
  for (size_t i = 0; i < extent[0]; ++i) {
    const double u = (double)i / (double)cell_extent[0];
    for (size_t j = 0; j < extent[1]; ++j) {
      const double v = (double)j / (double)cell_extent[1];
      const double displacement =
          w0 * sin((double)wave_number * pi * u) *
          sin((double)wave_number * pi * v);
      for (size_t k = 0; k < extent[2]; ++k) {
        const double w = (double)k / (double)cell_extent[2];
        const size_t index = gamera_no_index3(extent, i, j, k);
        vertices[index].value[0] =
            lower[0] + (upper[0] - lower[0]) * (u + displacement);
        vertices[index].value[1] =
            lower[1] + (upper[1] - lower[1]) *
                           (v + (double)convention * displacement);
        vertices[index].value[2] =
            lower[2] + (upper[2] - lower[2]) * w;
      }
    }
  }
  return 0;
}

int gamera_no_generate_spherical_vertices(
    const size_t cell_extent[3], const double lower[3], const double upper[3],
    gamera_no_vec3 *vertices) {
  if (!valid_mesh_arguments(cell_extent, lower, upper, vertices) ||
      lower[0] <= 0.0 || lower[1] < 0.0 || upper[1] > acos(-1.0)) {
    return -1;
  }
  size_t extent[3];
  vertex_extents(cell_extent, extent);
  for (size_t i = 0; i < extent[0]; ++i) {
    const double u = (double)i / (double)cell_extent[0];
    const double radius = lower[0] + (upper[0] - lower[0]) * u;
    for (size_t j = 0; j < extent[1]; ++j) {
      const double v = (double)j / (double)cell_extent[1];
      const double theta = lower[1] + (upper[1] - lower[1]) * v;
      for (size_t k = 0; k < extent[2]; ++k) {
        const double w = (double)k / (double)cell_extent[2];
        const double phi = lower[2] + (upper[2] - lower[2]) * w;
        const size_t index = gamera_no_index3(extent, i, j, k);
        vertices[index].value[0] = radius * sin(theta) * cos(phi);
        vertices[index].value[1] = radius * sin(theta) * sin(phi);
        vertices[index].value[2] = radius * cos(theta);
      }
    }
  }
  return 0;
}
