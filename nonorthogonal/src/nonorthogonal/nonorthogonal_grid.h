#ifndef GAMERA_NONORTHOGONAL_GRID_H
#define GAMERA_NONORTHOGONAL_GRID_H

#include "nonorthogonal_geometry.h"

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  size_t extent[3];
  gamera_no_face_geometry *value;
} gamera_no_face_grid;

typedef struct {
  size_t extent[3];
  gamera_no_edge_geometry *value;
  unsigned char *valid;
} gamera_no_edge_grid;

typedef struct {
  size_t cell_extent[3];
  size_t vertex_extent[3];
  gamera_no_vec3 *vertex;
  gamera_no_cell_geometry *cell;
  gamera_no_face_grid face[3];
  gamera_no_edge_grid edge[3];
} gamera_no_grid;

/* k is the contiguous dimension, matching GAMERA-OP's existing arrays. */
size_t gamera_no_index3(const size_t extent[3], size_t i, size_t j,
                        size_t k);
size_t gamera_no_element_count3(const size_t extent[3]);

/* Build owned geometry from a (ni+1)*(nj+1)*(nk+1) Cartesian vertex array. */
int gamera_no_grid_create(const size_t cell_extent[3],
                          const gamera_no_vec3 *vertices,
                          gamera_no_grid *grid);
void gamera_no_grid_destroy(gamera_no_grid *grid);

#ifdef __cplusplus
}
#endif

#endif
