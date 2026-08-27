#include "nonorthogonal_grid.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

size_t gamera_no_index3(const size_t extent[3], size_t i, size_t j,
                        size_t k) {
  return (i * extent[1] + j) * extent[2] + k;
}

size_t gamera_no_element_count3(const size_t extent[3]) {
  if (extent == NULL || extent[0] == 0 || extent[1] == 0 || extent[2] == 0 ||
      extent[0] > SIZE_MAX / extent[1] ||
      extent[0] * extent[1] > SIZE_MAX / extent[2]) {
    return 0;
  }
  return extent[0] * extent[1] * extent[2];
}

static int allocate_array(size_t count, size_t item_size, void **array) {
  if (array == NULL || count == 0 || item_size == 0 ||
      count > SIZE_MAX / item_size) {
    return -1;
  }
  *array = calloc(count, item_size);
  return *array == NULL ? -1 : 0;
}

static gamera_no_vec3 vertex_at(const gamera_no_grid *grid, size_t i,
                                size_t j, size_t k) {
  return grid->vertex[gamera_no_index3(grid->vertex_extent, i, j, k)];
}

static gamera_no_vec3 vector_add(gamera_no_vec3 a, gamera_no_vec3 b) {
  gamera_no_vec3 result;
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    result.value[d] = a.value[d] + b.value[d];
  }
  return result;
}

static gamera_no_vec3 vector_subtract(gamera_no_vec3 a, gamera_no_vec3 b) {
  gamera_no_vec3 result;
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    result.value[d] = a.value[d] - b.value[d];
  }
  return result;
}

static gamera_no_vec3 vector_scale(gamera_no_vec3 a, double scale) {
  gamera_no_vec3 result;
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    result.value[d] = scale * a.value[d];
  }
  return result;
}

static int build_cells(gamera_no_grid *grid) {
  const size_t corner_offset[8][3] = {
      {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
      {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
  int failed = 0;
#pragma omp parallel for collapse(3) reduction(| : failed) schedule(static)
  for (size_t i = 0; i < grid->cell_extent[0]; ++i) {
    for (size_t j = 0; j < grid->cell_extent[1]; ++j) {
      for (size_t k = 0; k < grid->cell_extent[2]; ++k) {
        gamera_no_vec3 corners[8];
        for (int corner = 0; corner < 8; ++corner) {
          corners[corner] =
              vertex_at(grid, i + corner_offset[corner][0],
                        j + corner_offset[corner][1],
                        k + corner_offset[corner][2]);
        }
        gamera_no_cell_geometry *cell =
            &grid->cell[gamera_no_index3(grid->cell_extent, i, j, k)];
        if (gamera_no_compute_cell_geometry(corners, cell) != 0) {
          failed = 1;
        }
      }
    }
  }
  return failed ? -1 : 0;
}

static void set_staggered_extents(gamera_no_grid *grid) {
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
      grid->face[direction].extent[axis] =
          grid->cell_extent[axis] + (axis == direction ? 1U : 0U);
      grid->edge[direction].extent[axis] =
          grid->cell_extent[axis] + (axis == direction ? 0U : 1U);
    }
  }
}

static int build_faces(gamera_no_grid *grid) {
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t count =
        gamera_no_element_count3(grid->face[direction].extent);
    if (allocate_array(count, sizeof(*grid->face[direction].value),
                       (void **)&grid->face[direction].value) != 0) {
      return -1;
    }

#pragma omp parallel for collapse(3) schedule(static)
    for (size_t i = 0; i < grid->face[direction].extent[0]; ++i) {
      for (size_t j = 0; j < grid->face[direction].extent[1]; ++j) {
        for (size_t k = 0; k < grid->face[direction].extent[2]; ++k) {
          size_t cell_index[3] = {i, j, k};
          int side = GAMERA_NO_LOWER;
          if (cell_index[direction] == grid->cell_extent[direction]) {
            --cell_index[direction];
            side = GAMERA_NO_UPPER;
          }
          const gamera_no_cell_geometry *cell =
              &grid->cell[gamera_no_index3(grid->cell_extent, cell_index[0],
                                           cell_index[1], cell_index[2])];
          grid->face[direction]
              .value[gamera_no_index3(grid->face[direction].extent, i, j, k)] =
              cell->face[direction][side];
        }
      }
    }
  }
  return 0;
}

static int edge_is_interior(int direction, const gamera_no_grid *grid,
                            size_t i, size_t j, size_t k) {
  if (direction == GAMERA_NO_I) {
    return j > 0 && j < grid->cell_extent[GAMERA_NO_J] && k > 0 &&
           k < grid->cell_extent[GAMERA_NO_K];
  }
  if (direction == GAMERA_NO_J) {
    return i > 0 && i < grid->cell_extent[GAMERA_NO_I] && k > 0 &&
           k < grid->cell_extent[GAMERA_NO_K];
  }
  return i > 0 && i < grid->cell_extent[GAMERA_NO_I] && j > 0 &&
         j < grid->cell_extent[GAMERA_NO_J];
}

static gamera_no_vec3 edge_transverse_average(int direction,
                                               const gamera_no_grid *grid,
                                               size_t i, size_t j, size_t k) {
  gamera_no_vec3 first;
  gamera_no_vec3 second;
  if (direction == GAMERA_NO_I) {
    first = vector_subtract(vertex_at(grid, i, j, k + 1),
                            vertex_at(grid, i, j, k - 1));
    second = vector_subtract(vertex_at(grid, i + 1, j, k + 1),
                             vertex_at(grid, i + 1, j, k - 1));
  } else if (direction == GAMERA_NO_J) {
    first = vector_subtract(vertex_at(grid, i + 1, j, k),
                            vertex_at(grid, i - 1, j, k));
    second = vector_subtract(vertex_at(grid, i + 1, j + 1, k),
                             vertex_at(grid, i - 1, j + 1, k));
  } else {
    first = vector_subtract(vertex_at(grid, i, j + 1, k),
                            vertex_at(grid, i, j - 1, k));
    second = vector_subtract(vertex_at(grid, i, j + 1, k + 1),
                             vertex_at(grid, i, j - 1, k + 1));
  }
  return vector_scale(vector_add(first, second), 0.5);
}

static int build_edges(gamera_no_grid *grid) {
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t count =
        gamera_no_element_count3(grid->edge[direction].extent);
    if (allocate_array(count, sizeof(*grid->edge[direction].value),
                       (void **)&grid->edge[direction].value) != 0 ||
        allocate_array(count, sizeof(*grid->edge[direction].valid),
                       (void **)&grid->edge[direction].valid) != 0) {
      return -1;
    }

    int failed = 0;
#pragma omp parallel for collapse(3) reduction(| : failed) schedule(static)
    for (size_t i = 0; i < grid->edge[direction].extent[0]; ++i) {
      for (size_t j = 0; j < grid->edge[direction].extent[1]; ++j) {
        for (size_t k = 0; k < grid->edge[direction].extent[2]; ++k) {
          const size_t index =
              gamera_no_index3(grid->edge[direction].extent, i, j, k);
          if (!edge_is_interior(direction, grid, i, j, k)) {
            continue;
          }
          size_t end_index[3] = {i, j, k};
          ++end_index[direction];
          const gamera_no_vec3 start = vertex_at(grid, i, j, k);
          const gamera_no_vec3 end =
              vertex_at(grid, end_index[0], end_index[1], end_index[2]);
          const gamera_no_vec3 average =
              edge_transverse_average(direction, grid, i, j, k);
          if (gamera_no_compute_edge_geometry(
                  start, end, average, &grid->edge[direction].value[index]) !=
              0) {
            failed = 1;
            continue;
          }
          grid->edge[direction].valid[index] = 1U;
        }
      }
    }
    if (failed) {
      return -1;
    }
  }
  return 0;
}

int gamera_no_grid_create(const size_t cell_extent[3],
                          const gamera_no_vec3 *vertices,
                          gamera_no_grid *grid) {
  if (cell_extent == NULL || vertices == NULL || grid == NULL) {
    return -1;
  }
  memset(grid, 0, sizeof(*grid));
  for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
    if (cell_extent[axis] == 0 || cell_extent[axis] == SIZE_MAX) {
      return -1;
    }
    grid->cell_extent[axis] = cell_extent[axis];
    grid->vertex_extent[axis] = cell_extent[axis] + 1U;
  }
  set_staggered_extents(grid);

  const size_t vertex_count = gamera_no_element_count3(grid->vertex_extent);
  const size_t cell_count = gamera_no_element_count3(grid->cell_extent);
  if (allocate_array(vertex_count, sizeof(*grid->vertex),
                     (void **)&grid->vertex) != 0 ||
      allocate_array(cell_count, sizeof(*grid->cell),
                     (void **)&grid->cell) != 0) {
    gamera_no_grid_destroy(grid);
    return -1;
  }
  memcpy(grid->vertex, vertices, vertex_count * sizeof(*grid->vertex));

  if (build_cells(grid) != 0 || build_faces(grid) != 0 ||
      build_edges(grid) != 0) {
    gamera_no_grid_destroy(grid);
    return -1;
  }
  return 0;
}

void gamera_no_grid_destroy(gamera_no_grid *grid) {
  if (grid == NULL) {
    return;
  }
  free(grid->vertex);
  free(grid->cell);
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    free(grid->face[direction].value);
    free(grid->edge[direction].value);
    free(grid->edge[direction].valid);
  }
  memset(grid, 0, sizeof(*grid));
}
