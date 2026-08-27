#include "nonorthogonal_mi_dpb.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define REQUIRE(condition)                                                     \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "requirement failed at %s:%d: %s\n", __FILE__,         \
              __LINE__, #condition);                                           \
      return EXIT_FAILURE;                                                     \
    }                                                                          \
  } while (0)

static double sector_mean(const double *boundary, size_t np, int night) {
  double sum = 0.0;
  size_t count = 0U;
  for (size_t j = 0U; j < np; ++j) {
    const double mlt = 24.0 * (double)j / (double)np;
    const int selected = night ? (mlt >= 20.0 || mlt <= 4.0)
                               : (mlt >= 10.0 && mlt <= 14.0);
    if (selected) {
      sum += boundary[j];
      ++count;
    }
  }
  return sum / (double)count;
}

static int test_offset_oval_and_nightside_limit(void) {
  const size_t np = 192U;
  double *boundary = (double *)calloc(np, sizeof(*boundary));
  REQUIRE(boundary != NULL);
  REQUIRE(gamera_mi_hybrid_dpb_boundary_from_radius(
              np, 20.5, 85.0, 1.0, boundary) == 0);
  REQUIRE(sector_mean(boundary, np, 0) > sector_mean(boundary, np, 1));
  free(boundary);
  return EXIT_SUCCESS;
}

static int test_quiet_initialization_and_mask(void) {
  const size_t np = 192U;
  const size_t nt = 20U;
  const size_t count = np * nt;
  const double pi = acos(-1.0);
  double *fac = (double *)calloc(count, sizeof(*fac));
  double *potential = (double *)calloc(count, sizeof(*potential));
  double *boundary = (double *)calloc(np, sizeof(*boundary));
  double *mask = (double *)calloc(count, sizeof(*mask));
  REQUIRE(fac != NULL && potential != NULL && boundary != NULL &&
          mask != NULL);
  const gamera_mi_hybrid_dpb_config config =
      gamera_mi_hybrid_dpb_default_config();
  gamera_mi_hybrid_dpb_state state = {0};
  gamera_mi_hybrid_dpb_stats stats;
  REQUIRE(gamera_mi_hybrid_dpb_update(
              &config, nt, np, pi / 6.0, pi / 6.0, fac, -1.0,
              potential, 120.0, &state, boundary, mask, &stats) == 0);
  REQUIRE(state.initialized);
  REQUIRE(isfinite(stats.nightside_boundary_deg));
  REQUIRE(stats.nightside_boundary_deg <= 70.0 + 1.0e-10);
  REQUIRE(stats.dayside_boundary_deg > stats.nightside_boundary_deg);
  for (size_t j = 0U; j < np; ++j) {
    REQUIRE(mask[j] > 0.999999);
    REQUIRE(mask[(nt - 1U) * np + j] < 1.0e-12);
  }
  free(fac);
  free(potential);
  free(boundary);
  free(mask);
  return EXIT_SUCCESS;
}

static int test_hybrid_evidence_and_slew(void) {
  const size_t np = 192U;
  const size_t nt = 25U;
  const size_t count = np * nt;
  const double pi = acos(-1.0);
  double *fac = (double *)calloc(count, sizeof(*fac));
  double *potential = (double *)calloc(count, sizeof(*potential));
  double *boundary = (double *)calloc(np, sizeof(*boundary));
  double *prior = (double *)calloc(np, sizeof(*prior));
  double *mask = (double *)calloc(count, sizeof(*mask));
  REQUIRE(fac != NULL && potential != NULL && boundary != NULL &&
          prior != NULL && mask != NULL);
  const double maximum_colatitude = 35.0 * pi / 180.0;
  for (size_t i = 0U; i < nt; ++i) {
    const double latitude = 90.0 -
        180.0 * maximum_colatitude * (double)i /
            ((double)(nt - 1U) * pi);
    for (size_t j = 0U; j < np; ++j) {
      const double phi = 2.0 * pi * (double)j / (double)np;
      const double cell = sin(phi) > 0.0 ? 1.0 : -1.0;
      potential[i * np + j] =
          25000.0 * cell /
          (1.0 + exp(-(latitude - 66.0) / 0.8));
      fac[i * np + j] =
          -1.0e-6 * cos(phi) * exp(-0.5 * pow((latitude - 67.0) / 1.0, 2.0)) +
           0.8e-6 * cos(phi) * exp(-0.5 * pow((latitude - 74.0) / 1.0, 2.0));
    }
  }
  const gamera_mi_hybrid_dpb_config config =
      gamera_mi_hybrid_dpb_default_config();
  gamera_mi_hybrid_dpb_state state = {0};
  gamera_mi_hybrid_dpb_stats stats;
  REQUIRE(gamera_mi_hybrid_dpb_update(
              &config, nt, np, maximum_colatitude, maximum_colatitude,
              fac, 1.0, potential, 120.0, &state, boundary, mask,
              &stats) == 0);
  REQUIRE(stats.evidence_fraction > 0.25);
  for (size_t j = 0U; j < np; ++j) {
    prior[j] = boundary[j];
  }
  for (size_t index = 0U; index < count; ++index) {
    potential[index] *= 0.1;
    fac[index] *= 0.1;
  }
  REQUIRE(gamera_mi_hybrid_dpb_update(
              &config, nt, np, maximum_colatitude, maximum_colatitude,
              fac, 1.0, potential, 120.0, &state, boundary, mask,
              &stats) == 0);
  double maximum_step = 0.0;
  for (size_t j = 0U; j < np; ++j) {
    maximum_step = fmax(maximum_step, fabs(boundary[j] - prior[j]));
  }
  REQUIRE(maximum_step <= 1.0 + 1.0e-10);
  REQUIRE(stats.nightside_boundary_deg <= 70.0 + 1.0e-10);
  free(fac);
  free(potential);
  free(boundary);
  free(prior);
  free(mask);
  return EXIT_SUCCESS;
}

int main(void) {
  if (test_offset_oval_and_nightside_limit() != EXIT_SUCCESS ||
      test_quiet_initialization_and_mask() != EXIT_SUCCESS ||
      test_hybrid_evidence_and_slew() != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }
  puts("nonorthogonal hybrid DPB tests passed");
  return EXIT_SUCCESS;
}
