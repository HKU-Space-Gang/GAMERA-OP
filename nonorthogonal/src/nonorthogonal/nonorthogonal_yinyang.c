#include "nonorthogonal_yinyang.h"

#include <math.h>
#include <stddef.h>

static int valid_patch(int patch) {
  return patch == GAMERA_NO_YIN_PATCH || patch == GAMERA_NO_YANG_PATCH;
}

static int finite_vector(gamera_no_vec3 vector) {
  return isfinite(vector.value[0]) && isfinite(vector.value[1]) &&
         isfinite(vector.value[2]);
}

int gamera_no_yinyang_local_to_global(int patch, gamera_no_vec3 local,
                                      gamera_no_vec3 *global) {
  if (!valid_patch(patch) || !finite_vector(local) || global == NULL) {
    return -1;
  }
  if (patch == GAMERA_NO_YIN_PATCH) {
    *global = local;
  } else {
    *global = (gamera_no_vec3){{-local.value[0], local.value[2],
                                local.value[1]}};
  }
  return 0;
}

int gamera_no_yinyang_global_to_local(int patch, gamera_no_vec3 global,
                                      gamera_no_vec3 *local) {
  /* The Yang rotation matrix is symmetric and therefore self-inverse. */
  return gamera_no_yinyang_local_to_global(patch, global, local);
}

int gamera_no_yinyang_logical_to_global(int patch, double radius,
                                        double theta, double phi,
                                        gamera_no_vec3 *global) {
  if (!valid_patch(patch) || global == NULL || !isfinite(radius) ||
      !isfinite(theta) || !isfinite(phi) || !(radius > 0.0)) {
    return -1;
  }
  const gamera_no_vec3 local = {
      {radius * sin(theta) * cos(phi), radius * sin(theta) * sin(phi),
       radius * cos(theta)}};
  return gamera_no_yinyang_local_to_global(patch, local, global);
}

int gamera_no_yinyang_global_to_logical(int patch, gamera_no_vec3 global,
                                        double *radius, double *theta,
                                        double *phi) {
  if (!valid_patch(patch) || !finite_vector(global) || radius == NULL ||
      theta == NULL || phi == NULL) {
    return -1;
  }
  gamera_no_vec3 local;
  if (gamera_no_yinyang_global_to_local(patch, global, &local) != 0) {
    return -1;
  }
  const double r = sqrt(local.value[0] * local.value[0] +
                        local.value[1] * local.value[1] +
                        local.value[2] * local.value[2]);
  if (!(r > 0.0) || !isfinite(r)) {
    return -1;
  }
  double cosine = local.value[2] / r;
  if (cosine > 1.0) {
    cosine = 1.0;
  } else if (cosine < -1.0) {
    cosine = -1.0;
  }
  *radius = r;
  *theta = acos(cosine);
  *phi = atan2(local.value[1], local.value[0]);
  return 0;
}

static int valid_angular_domain(
    const gamera_no_yinyang_angular_domain *domain) {
  return domain != NULL && isfinite(domain->theta_min) &&
         isfinite(domain->theta_max) && isfinite(domain->phi_min) &&
         isfinite(domain->phi_max) &&
         domain->theta_max > domain->theta_min &&
         domain->phi_max > domain->phi_min && domain->theta_cells > 0U &&
         domain->phi_cells > 0U &&
         isfinite(domain->boundary_tolerance_cells) &&
         domain->boundary_tolerance_cells >= 0.0;
}

int gamera_no_yinyang_angular_margin(
    int patch, gamera_no_vec3 point,
    const gamera_no_yinyang_angular_domain *domain, int *valid,
    double *margin_cells) {
  double radius;
  double theta;
  double phi;
  if (!valid_patch(patch) || !valid_angular_domain(domain) || valid == NULL ||
      margin_cells == NULL ||
      gamera_no_yinyang_global_to_logical(patch, point, &radius, &theta,
                                          &phi) != 0) {
    return -1;
  }
  const double theta_spacing =
      (domain->theta_max - domain->theta_min) /
      (double)domain->theta_cells;
  const double phi_spacing =
      (domain->phi_max - domain->phi_min) / (double)domain->phi_cells;
  const double q_theta = (theta - domain->theta_min) / theta_spacing - 0.5;
  const double q_phi = (phi - domain->phi_min) / phi_spacing - 0.5;
  const double tolerance = domain->boundary_tolerance_cells;
  *valid = q_theta >= -0.5 - tolerance &&
           q_theta <= (double)domain->theta_cells - 0.5 + tolerance &&
           q_phi >= -0.5 - tolerance &&
           q_phi <= (double)domain->phi_cells - 0.5 + tolerance;
  *margin_cells =
      fmin(fmin(q_theta + 0.5,
                (double)domain->theta_cells - 0.5 - q_theta),
           fmin(q_phi + 0.5,
                (double)domain->phi_cells - 0.5 - q_phi));
  return isfinite(*margin_cells) ? 0 : -1;
}

int gamera_no_yinyang_composite_owner(
    gamera_no_vec3 point, const gamera_no_yinyang_angular_domain *domain,
    int *owner_patch, int *in_overlap, double margin_cells[2]) {
  int valid[2];
  double margin[2];
  if (owner_patch == NULL || in_overlap == NULL ||
      gamera_no_yinyang_angular_margin(GAMERA_NO_YIN_PATCH, point, domain,
                                       &valid[0], &margin[0]) != 0 ||
      gamera_no_yinyang_angular_margin(GAMERA_NO_YANG_PATCH, point, domain,
                                       &valid[1], &margin[1]) != 0 ||
      (!valid[0] && !valid[1])) {
    return -1;
  }
  *in_overlap = valid[0] && valid[1];
  if (!valid[1]) {
    *owner_patch = GAMERA_NO_YIN_PATCH;
  } else if (!valid[0]) {
    *owner_patch = GAMERA_NO_YANG_PATCH;
  } else {
    const double tie_tolerance = 1.0e-12;
    *owner_patch = margin[1] > margin[0] + tie_tolerance
                       ? GAMERA_NO_YANG_PATCH
                       : GAMERA_NO_YIN_PATCH;
  }
  if (margin_cells != NULL) {
    margin_cells[0] = margin[0];
    margin_cells[1] = margin[1];
  }
  return 0;
}
