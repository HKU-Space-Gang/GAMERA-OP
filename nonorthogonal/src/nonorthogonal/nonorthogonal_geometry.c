#include "nonorthogonal_geometry.h"

#include <float.h>
#include <math.h>
#include <stddef.h>

enum { GAMERA_NO_GAUSS_POINTS = 12 };

static const double gauss_node[GAMERA_NO_GAUSS_POINTS] = {
    0.1252334085114689154724,  0.3678314989981801937527,
    0.5873179542866174472967,  0.7699026741943046870369,
    0.9041172563704748566785,  0.9815606342467192506906,
    -0.1252334085114689154724, -0.3678314989981801937527,
    -0.5873179542866174472967, -0.7699026741943046870369,
    -0.9041172563704748566785, -0.9815606342467192506906};

static const double gauss_weight[GAMERA_NO_GAUSS_POINTS] = {
    0.2491470458134027850006, 0.2334925365383548087610,
    0.2031674267230659217490, 0.1600783285433462263350,
    0.1069393259953184309603, 0.0471753363865118271946,
    0.2491470458134027850006, 0.2334925365383548087610,
    0.2031674267230659217490, 0.1600783285433462263350,
    0.1069393259953184309603, 0.0471753363865118271946};

static const int logical_sign[8][GAMERA_NO_DIM] = {
    {-1, -1, -1}, {+1, -1, -1}, {+1, +1, -1}, {-1, +1, -1},
    {-1, -1, +1}, {+1, -1, +1}, {+1, +1, +1}, {-1, +1, +1}};

/* f0=SW, f1=SE, f2=NW, f3=NE, matching metric.F90:faceCoords. */
static const int face_corner[GAMERA_NO_DIM][2][4] = {
    {{0, 3, 4, 7}, {1, 2, 5, 6}},
    {{0, 4, 1, 5}, {3, 7, 2, 6}},
    {{0, 1, 3, 2}, {4, 5, 7, 6}}};

static gamera_no_vec3 vec_zero(void) {
  gamera_no_vec3 result = {{0.0, 0.0, 0.0}};
  return result;
}

static gamera_no_vec3 vec_add(gamera_no_vec3 a, gamera_no_vec3 b) {
  gamera_no_vec3 result;
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    result.value[d] = a.value[d] + b.value[d];
  }
  return result;
}

static gamera_no_vec3 vec_sub(gamera_no_vec3 a, gamera_no_vec3 b) {
  gamera_no_vec3 result;
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    result.value[d] = a.value[d] - b.value[d];
  }
  return result;
}

static gamera_no_vec3 vec_scale(gamera_no_vec3 vector, double scale) {
  gamera_no_vec3 result;
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    result.value[d] = vector.value[d] * scale;
  }
  return result;
}

static double vec_dot(gamera_no_vec3 a, gamera_no_vec3 b) {
  double result = 0.0;
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    result += a.value[d] * b.value[d];
  }
  return result;
}

static gamera_no_vec3 vec_cross(gamera_no_vec3 a, gamera_no_vec3 b) {
  gamera_no_vec3 result = {{a.value[1] * b.value[2] - a.value[2] * b.value[1],
                            a.value[2] * b.value[0] - a.value[0] * b.value[2],
                            a.value[0] * b.value[1] - a.value[1] * b.value[0]}};
  return result;
}

static double vec_norm(gamera_no_vec3 vector) {
  return sqrt(vec_dot(vector, vector));
}

static int vec_normalize(gamera_no_vec3 vector, gamera_no_vec3 *unit) {
  const double norm = vec_norm(vector);
  if (!isfinite(norm) || norm <= DBL_MIN) {
    return -1;
  }
  *unit = vec_scale(vector, 1.0 / norm);
  return 0;
}

static gamera_no_vec3 bilinear_point(const gamera_no_vec3 face[4], double eta,
                                     double psi) {
  const gamera_no_vec3 d_eta = vec_sub(face[1], face[0]);
  const gamera_no_vec3 d_psi = vec_sub(face[2], face[0]);
  const gamera_no_vec3 d_cross =
      vec_add(vec_sub(face[3], face[2]), vec_sub(face[0], face[1]));
  return vec_add(
      face[0],
      vec_add(vec_scale(d_eta, eta),
              vec_add(vec_scale(d_psi, psi),
                      vec_scale(d_cross, eta * psi))));
}

static int compute_face_geometry(const gamera_no_vec3 face[4],
                                 gamera_no_face_geometry *geometry) {
  const gamera_no_vec3 d_eta = vec_sub(face[1], face[0]);
  const gamera_no_vec3 d_psi = vec_sub(face[2], face[0]);
  const gamera_no_vec3 d_cross =
      vec_add(vec_sub(face[3], face[2]), vec_sub(face[0], face[1]));

  geometry->area = 0.0;
  geometry->area_vector = vec_zero();
  geometry->centroid = vec_zero();

  for (int i = 0; i < GAMERA_NO_GAUSS_POINTS; ++i) {
    const double eta = 0.5 * (1.0 + gauss_node[i]);
    for (int j = 0; j < GAMERA_NO_GAUSS_POINTS; ++j) {
      const double psi = 0.5 * (1.0 + gauss_node[j]);
      const double weight = 0.25 * gauss_weight[i] * gauss_weight[j];
      const gamera_no_vec3 tangent_eta =
          vec_add(d_eta, vec_scale(d_cross, psi));
      const gamera_no_vec3 tangent_psi =
          vec_add(d_psi, vec_scale(d_cross, eta));
      const gamera_no_vec3 oriented_area =
          vec_cross(tangent_eta, tangent_psi);
      const double area_weight = weight * vec_norm(oriented_area);
      const gamera_no_vec3 point = bilinear_point(face, eta, psi);

      geometry->area += area_weight;
      geometry->area_vector =
          vec_add(geometry->area_vector, vec_scale(oriented_area, weight));
      geometry->centroid =
          vec_add(geometry->centroid, vec_scale(point, area_weight));
    }
  }

  if (!isfinite(geometry->area) || geometry->area <= DBL_MIN) {
    return -1;
  }
  geometry->centroid = vec_scale(geometry->centroid, 1.0 / geometry->area);

  /* The production Fortran path constructs Tf from these centerline vectors. */
  const gamera_no_vec3 center_eta =
      vec_scale(vec_add(vec_sub(face[3], face[2]),
                        vec_sub(face[1], face[0])),
                0.5);
  const gamera_no_vec3 center_psi =
      vec_scale(vec_add(vec_sub(face[3], face[1]),
                        vec_sub(face[2], face[0])),
                0.5);
  if (vec_normalize(vec_cross(center_eta, center_psi), &geometry->normal) != 0 ||
      vec_normalize(center_psi, &geometry->tangent2) != 0 ||
      vec_normalize(vec_cross(geometry->tangent2, geometry->normal),
                    &geometry->tangent1) != 0) {
    return -1;
  }

  return 0;
}

static int compute_volume_geometry(const gamera_no_vec3 corners[8],
                                   double *volume,
                                   gamera_no_vec3 *centroid) {
  *volume = 0.0;
  *centroid = vec_zero();

  for (int i = 0; i < GAMERA_NO_GAUSS_POINTS; ++i) {
    const double xi = gauss_node[i];
    for (int j = 0; j < GAMERA_NO_GAUSS_POINTS; ++j) {
      const double eta = gauss_node[j];
      for (int k = 0; k < GAMERA_NO_GAUSS_POINTS; ++k) {
        const double zeta = gauss_node[k];
        gamera_no_vec3 point = vec_zero();
        gamera_no_vec3 derivative[3] = {vec_zero(), vec_zero(), vec_zero()};

        for (int corner = 0; corner < 8; ++corner) {
          const double sx = (double)logical_sign[corner][GAMERA_NO_I];
          const double sy = (double)logical_sign[corner][GAMERA_NO_J];
          const double sz = (double)logical_sign[corner][GAMERA_NO_K];
          const double fx = 1.0 + sx * xi;
          const double fy = 1.0 + sy * eta;
          const double fz = 1.0 + sz * zeta;
          const double shape = 0.125 * fx * fy * fz;
          const double dshape[3] = {0.125 * sx * fy * fz,
                                    0.125 * fx * sy * fz,
                                    0.125 * fx * fy * sz};

          point = vec_add(point, vec_scale(corners[corner], shape));
          for (int d = 0; d < GAMERA_NO_DIM; ++d) {
            derivative[d] =
                vec_add(derivative[d],
                        vec_scale(corners[corner], dshape[d]));
          }
        }

        const double determinant =
            vec_dot(derivative[0], vec_cross(derivative[1], derivative[2]));
        const double weighted_volume =
            fabs(determinant) * gauss_weight[i] * gauss_weight[j] *
            gauss_weight[k];
        *volume += weighted_volume;
        *centroid = vec_add(*centroid, vec_scale(point, weighted_volume));
      }
    }
  }

  if (!isfinite(*volume) || *volume <= DBL_MIN) {
    return -1;
  }
  *centroid = vec_scale(*centroid, 1.0 / *volume);
  return 0;
}

int gamera_no_compute_cell_geometry(
    const gamera_no_vec3 corners[8], gamera_no_cell_geometry *geometry) {
  if (corners == NULL || geometry == NULL) {
    return -1;
  }

  if (compute_volume_geometry(corners, &geometry->volume,
                              &geometry->centroid) != 0) {
    return -1;
  }

  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    for (int side = 0; side < 2; ++side) {
      gamera_no_vec3 face[4];
      for (int corner = 0; corner < 4; ++corner) {
        face[corner] = corners[face_corner[direction][side][corner]];
      }
      if (compute_face_geometry(face, &geometry->face[direction][side]) != 0) {
        return -1;
      }
    }
    const double maximum_area =
        fmax(geometry->face[direction][GAMERA_NO_LOWER].area,
             geometry->face[direction][GAMERA_NO_UPPER].area);
    geometry->cfl_length[direction] = geometry->volume / maximum_area;
  }

  return 0;
}

int gamera_no_compute_edge_geometry(gamera_no_vec3 start, gamera_no_vec3 end,
                                    gamera_no_vec3 transverse_average,
                                    gamera_no_edge_geometry *geometry) {
  if (geometry == NULL) {
    return -1;
  }
  const gamera_no_vec3 edge = vec_sub(end, start);
  geometry->length = vec_norm(edge);
  if (!isfinite(geometry->length) || geometry->length <= DBL_MIN ||
      vec_normalize(edge, &geometry->normal) != 0 ||
      vec_normalize(vec_cross(transverse_average, geometry->normal),
                    &geometry->tangent1) != 0 ||
      vec_normalize(vec_cross(geometry->normal, geometry->tangent1),
                    &geometry->tangent2) != 0) {
    return -1;
  }
  return 0;
}
