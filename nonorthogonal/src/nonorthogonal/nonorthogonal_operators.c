#include "nonorthogonal_operators.h"

#include "nonorthogonal_reconstruction.h"

#include <float.h>
#include <math.h>
#include <stddef.h>

static double vector_dot(gamera_no_vec3 a, gamera_no_vec3 b) {
  double result = 0.0;
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    result += a.value[d] * b.value[d];
  }
  return result;
}

int gamera_no_flux_to_cell_field(const gamera_no_cell_geometry *geometry,
                                 const double face_flux[GAMERA_NO_DIM][2],
                                 gamera_no_vec3 *field,
                                 double *net_flux) {
  if (geometry == NULL || face_flux == NULL || field == NULL ||
      !isfinite(geometry->volume) || geometry->volume <= DBL_MIN) {
    return -1;
  }

  double divergence = 0.0;
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    field->value[component] = 0.0;
  }

  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const double lower_flux = face_flux[direction][GAMERA_NO_LOWER];
    const double upper_flux = face_flux[direction][GAMERA_NO_UPPER];
    divergence += upper_flux - lower_flux;
    for (int component = 0; component < GAMERA_NO_DIM; ++component) {
      field->value[component] +=
          upper_flux * geometry->face[direction][GAMERA_NO_UPPER]
                           .centroid.value[component] -
          lower_flux * geometry->face[direction][GAMERA_NO_LOWER]
                           .centroid.value[component];
    }
  }

  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    field->value[component] =
        (field->value[component] -
         divergence * geometry->centroid.value[component]) /
        geometry->volume;
  }
  if (net_flux != NULL) {
    *net_flux = divergence;
  }
  return 0;
}

int gamera_no_solve_edge_field(const double normal1[2],
                               const double normal2[2],
                               double face_field1, double face_field2,
                               double edge_field[2]) {
  if (normal1 == NULL || normal2 == NULL || edge_field == NULL) {
    return -1;
  }
  const double determinant =
      normal1[0] * normal2[1] - normal2[0] * normal1[1];
  const double matrix_scale =
      fmax(fmax(fabs(normal1[0]), fabs(normal1[1])),
           fmax(fabs(normal2[0]), fabs(normal2[1])));
  if (!isfinite(determinant) || matrix_scale <= DBL_MIN ||
      fabs(determinant) <= 64.0 * DBL_EPSILON * matrix_scale * matrix_scale) {
    return -1;
  }

  const double inverse_determinant = 1.0 / determinant;
  edge_field[0] = inverse_determinant *
                  (normal2[1] * face_field1 - normal1[1] * face_field2);
  edge_field[1] = inverse_determinant *
                  (-normal2[0] * face_field1 + normal1[0] * face_field2);
  return 0;
}

int gamera_no_interpolate_face_normal_to_edge(
    const double face_area[8], const gamera_no_vec3 face_normal[8],
    const gamera_no_edge_geometry *edge, double transverse_normal[2]) {
  if (face_area == NULL || face_normal == NULL || edge == NULL ||
      transverse_normal == NULL) {
    return -1;
  }

  const double interpolated_area = gamera_no_central6(face_area);
  if (!isfinite(interpolated_area) || fabs(interpolated_area) <= DBL_MIN) {
    return -1;
  }

  gamera_no_vec3 normal_at_edge = {{0.0, 0.0, 0.0}};
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    double weighted_normal[8];
    for (int n = 0; n < 8; ++n) {
      weighted_normal[n] = face_area[n] * face_normal[n].value[component];
    }
    normal_at_edge.value[component] =
        gamera_no_central6(weighted_normal) / interpolated_area;
  }

  const double projected1 = vector_dot(normal_at_edge, edge->tangent1);
  const double projected2 = vector_dot(normal_at_edge, edge->tangent2);
  const double projected_norm = hypot(projected1, projected2);
  if (!isfinite(projected_norm) || projected_norm <= DBL_MIN) {
    return -1;
  }
  transverse_normal[0] = projected1 / projected_norm;
  transverse_normal[1] = projected2 / projected_norm;
  return 0;
}

void gamera_no_ct_face_increments(const double emf_i[2][2],
                                  const double emf_j[2][2],
                                  const double emf_k[2][2], double dt,
                                  double face_increment[GAMERA_NO_DIM][2]) {
  for (int i_side = 0; i_side < 2; ++i_side) {
    const double curl = emf_k[i_side][GAMERA_NO_UPPER] -
                        emf_k[i_side][GAMERA_NO_LOWER] -
                        emf_j[i_side][GAMERA_NO_UPPER] +
                        emf_j[i_side][GAMERA_NO_LOWER];
    face_increment[GAMERA_NO_I][i_side] = -dt * curl;
  }

  for (int j_side = 0; j_side < 2; ++j_side) {
    const double curl = emf_i[j_side][GAMERA_NO_UPPER] -
                        emf_i[j_side][GAMERA_NO_LOWER] -
                        emf_k[GAMERA_NO_UPPER][j_side] +
                        emf_k[GAMERA_NO_LOWER][j_side];
    face_increment[GAMERA_NO_J][j_side] = -dt * curl;
  }

  for (int k_side = 0; k_side < 2; ++k_side) {
    const double curl = emf_j[GAMERA_NO_UPPER][k_side] -
                        emf_j[GAMERA_NO_LOWER][k_side] -
                        emf_i[GAMERA_NO_UPPER][k_side] +
                        emf_i[GAMERA_NO_LOWER][k_side];
    face_increment[GAMERA_NO_K][k_side] = -dt * curl;
  }
}

int gamera_no_compute_edge_emf(
    double velocity_tangent1, double velocity_tangent2,
    double magnetic_tangent1, double magnetic_tangent2,
    double magnetic_current_jump, double edge_density,
    const gamera_no_edge_geometry *edge, double diffusion_coefficient,
    bool use_boris, double light_speed, double cfl, double dt,
    double *line_integrated_emf, double *diffusion_speed) {
  if (edge == NULL || line_integrated_emf == NULL ||
      diffusion_speed == NULL || !isfinite(velocity_tangent1) ||
      !isfinite(velocity_tangent2) || !isfinite(magnetic_tangent1) ||
      !isfinite(magnetic_tangent2) ||
      !isfinite(magnetic_current_jump) || !isfinite(edge_density) ||
      edge_density <= DBL_MIN || !isfinite(edge->length) ||
      edge->length <= DBL_MIN || !isfinite(diffusion_coefficient) ||
      !isfinite(cfl) || cfl < 0.0 || !isfinite(dt) || dt <= DBL_MIN ||
      (use_boris && (!isfinite(light_speed) || light_speed <= DBL_MIN))) {
    return -1;
  }

  const double transverse_flow_speed =
      hypot(velocity_tangent1, velocity_tangent2);
  double alfven_speed =
      hypot(magnetic_tangent1, magnetic_tangent2) / sqrt(edge_density);
  if (use_boris) {
    alfven_speed = light_speed * alfven_speed /
                   hypot(light_speed, alfven_speed);
  }
  *diffusion_speed = transverse_flow_speed + alfven_speed;
  *diffusion_speed = fmin(*diffusion_speed, cfl * edge->length / dt);

  const double ideal_electric =
      -(velocity_tangent1 * magnetic_tangent2 -
        velocity_tangent2 * magnetic_tangent1);
  *line_integrated_emf =
      (ideal_electric + diffusion_coefficient * (*diffusion_speed) *
                            magnetic_current_jump) *
      edge->length;
  return 0;
}
