#include "nonorthogonal_reconstruction.h"

#include <float.h>
#include <math.h>
#include <stddef.h>

static const double central8_coefficient[GAMERA_NO_RECON_STENCIL] = {
    -3.0 / 840.0, 29.0 / 840.0,  -139.0 / 840.0, 533.0 / 840.0,
    533.0 / 840.0, -139.0 / 840.0, 29.0 / 840.0,  -3.0 / 840.0};

static const double central6_coefficient[GAMERA_NO_RECON_STENCIL] = {
    0.0, 1.0 / 60.0, -8.0 / 60.0, 37.0 / 60.0,
    37.0 / 60.0, -8.0 / 60.0, 1.0 / 60.0, 0.0};

static const double upwind7_coefficient[7] = {
    -3.0 / 420.0, 25.0 / 420.0, -101.0 / 420.0, 319.0 / 420.0,
    214.0 / 420.0, -38.0 / 420.0, 4.0 / 420.0};

static double dot8(const double values[GAMERA_NO_RECON_STENCIL],
                   const double coefficients[GAMERA_NO_RECON_STENCIL]) {
  double result = 0.0;
  for (int n = 0; n < GAMERA_NO_RECON_STENCIL; ++n) {
    result += coefficients[n] * values[n];
  }
  return result;
}

static double upwind7(const double values[7]) {
  double result = 0.0;
  for (int n = 0; n < 7; ++n) {
    result += upwind7_coefficient[n] * values[n];
  }
  return result;
}

/* Fortran SIGN(1.0, x), including its nonzero result at x=0. */
static double fortran_unit_sign(double value) { return copysign(1.0, value); }

static double pdm_left(double q0, double q1, double q2,
                       double interpolated, double pdmb_code) {
  const double bounded = fmax(fmin(q1, q2), fmin(interpolated, fmax(q1, q2)));
  const double delta0 = pdmb_code * (q1 - q0);
  const double delta1 = pdmb_code * (q2 - q1);
  const double sign0 = fortran_unit_sign(delta0);
  const double sign1 = fortran_unit_sign(delta1);
  const double same_sign_factor = fabs(sign0 + sign1);
  const double local_slope = bounded - q1;
  return bounded -
         sign1 * fmax(0.0, fabs(local_slope) -
                               same_sign_factor * fabs(delta0));
}

double gamera_no_central8(const double stencil[GAMERA_NO_RECON_STENCIL]) {
  return dot8(stencil, central8_coefficient);
}

double gamera_no_central6(const double stencil[GAMERA_NO_RECON_STENCIL]) {
  return dot8(stencil, central6_coefficient);
}

int gamera_no_reconstruct_up7_pdm(
    const double volume[GAMERA_NO_RECON_STENCIL],
    const double primitive[GAMERA_NO_RECON_STENCIL], double pdmb_code,
    double *left, double *right) {
  if (volume == NULL || primitive == NULL || left == NULL || right == NULL ||
      !isfinite(pdmb_code) || pdmb_code < 0.0) {
    return -1;
  }

  double weighted[GAMERA_NO_RECON_STENCIL];
  for (int n = 0; n < GAMERA_NO_RECON_STENCIL; ++n) {
    weighted[n] = volume[n] * primitive[n];
  }
  const double interface_volume = gamera_no_central8(volume);
  if (!isfinite(interface_volume) || fabs(interface_volume) <= DBL_MIN) {
    return -1;
  }

  double left_stencil[7];
  double right_stencil[7];
  for (int n = 0; n < 7; ++n) {
    left_stencil[n] = weighted[n];
    right_stencil[n] = weighted[7 - n];
  }
  const double high_left = upwind7(left_stencil) / interface_volume;
  const double high_right = upwind7(right_stencil) / interface_volume;

  *left = pdm_left(primitive[2], primitive[3], primitive[4], high_left,
                   pdmb_code);
  *right = pdm_left(primitive[5], primitive[4], primitive[3], high_right,
                    pdmb_code);
  return 0;
}
