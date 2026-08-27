#ifndef GAMERA_NONORTHOGONAL_GEOMETRY_H
#define GAMERA_NONORTHOGONAL_GEOMETRY_H

#ifdef __cplusplus
extern "C" {
#endif

enum {
  GAMERA_NO_I = 0,
  GAMERA_NO_J = 1,
  GAMERA_NO_K = 2,
  GAMERA_NO_DIM = 3,
  GAMERA_NO_LOWER = 0,
  GAMERA_NO_UPPER = 1
};

typedef struct {
  double value[GAMERA_NO_DIM];
} gamera_no_vec3;

typedef struct {
  double area;
  gamera_no_vec3 area_vector;
  gamera_no_vec3 centroid;
  gamera_no_vec3 normal;
  gamera_no_vec3 tangent1;
  gamera_no_vec3 tangent2;
} gamera_no_face_geometry;

typedef struct {
  double volume;
  gamera_no_vec3 centroid;
  gamera_no_face_geometry face[GAMERA_NO_DIM][2];
  double cfl_length[GAMERA_NO_DIM];
} gamera_no_cell_geometry;

typedef struct {
  double length;
  gamera_no_vec3 normal;
  gamera_no_vec3 tangent1;
  gamera_no_vec3 tangent2;
} gamera_no_edge_geometry;

/*
 * Compute the geometry of a trilinearly mapped hexahedral cell using the
 * 12-point-per-dimension Gauss-Legendre rule used by GAMERA.
 *
 * Corner order:
 *   0: (0,0,0), 1: (1,0,0), 2: (1,1,0), 3: (0,1,0),
 *   4: (0,0,1), 5: (1,0,1), 6: (1,1,1), 7: (0,1,1).
 *
 * Face normals and area vectors point in the positive logical direction for
 * both the lower and upper face. Returns 0 on success and -1 for a degenerate
 * or non-finite cell.
 */
int gamera_no_compute_cell_geometry(
    const gamera_no_vec3 corners[8], gamera_no_cell_geometry *geometry);

/* Fortran init.F90 edge frame: N=edge, T1=eAvg x N, T2=N x T1. */
int gamera_no_compute_edge_geometry(gamera_no_vec3 start, gamera_no_vec3 end,
                                    gamera_no_vec3 transverse_average,
                                    gamera_no_edge_geometry *geometry);

#ifdef __cplusplus
}
#endif

#endif
