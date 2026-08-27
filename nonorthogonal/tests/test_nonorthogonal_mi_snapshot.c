#include "nonorthogonal_mi_coupling.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures;

static void expect_true(const char *name, int condition) {
  if (!condition) {
    fprintf(stderr, "FAIL %s\n", name);
    ++failures;
  }
}

static void expect_near(const char *name, double actual, double expected,
                        double tolerance) {
  if (!isfinite(actual) || fabs(actual - expected) > tolerance) {
    fprintf(stderr,
            "FAIL %s: actual=%.17g expected=%.17g tolerance=%.3e\n",
            name, actual, expected, tolerance);
    ++failures;
  }
}

static int cleared_view(const gamera_mi_north_snapshot_view *view) {
  return view->ready == 0 && view->epoch_s == 0.0 &&
         view->generation == 0U && view->theta_points == 0U &&
         view->azimuth_points == 0U && view->ionosphere_radius_re == 0.0 &&
         view->trace_radius_re == 0.0 && view->colatitude_rad == NULL &&
         view->longitude_rad == NULL && view->potential_v == NULL;
}

static void fill_potential(double *values, size_t nt, size_t np,
                           double offset) {
  for (size_t theta = 0U; theta < nt; ++theta) {
    for (size_t azimuth = 0U; azimuth < np; ++azimuth) {
      values[theta * np + azimuth] =
          offset + 100.0 * (double)theta + (double)azimuth;
    }
  }
}

static void test_snapshot_lifecycle(void) {
  const size_t nt = 4U;
  const size_t np = 8U;
  const size_t count = nt * np;
  const double pi = acos(-1.0);
  double first[32];
  double second[32];
  double invalid[32];
  gamera_mi_north_snapshot_view view;

  memset(&view, 0xa5, sizeof(view));
  expect_true("disabled/unprepared snapshot is unavailable",
              gamera_mi_coupling_borrow_north_snapshot(0U, &view) ==
                  GAMERA_MI_SNAPSHOT_UNAVAILABLE);
  expect_true("unavailable clears view", cleared_view(&view));
  expect_true("null view is invalid",
              gamera_mi_coupling_borrow_north_snapshot(0U, NULL) ==
                  GAMERA_MI_SNAPSHOT_INVALID_ARGUMENT);

  expect_true("snapshot test prepare",
              gamera_mi_snapshot_test_prepare(nt, np, pi / 6.0, 1.0,
                                              4.24795) == 0);
  memset(&view, 0xa5, sizeof(view));
  expect_true("prepared snapshot is explicitly unready",
              gamera_mi_coupling_borrow_north_snapshot(0U, &view) ==
                  GAMERA_MI_SNAPSHOT_UNREADY);
  expect_true("unready clears view", cleared_view(&view));

  fill_potential(first, nt, np, 1000.0);
  expect_true("first successful publication",
              gamera_mi_snapshot_test_publish(123.25, first) == 0);
  expect_true("borrow first generation",
              gamera_mi_coupling_borrow_north_snapshot(0U, &view) ==
                  GAMERA_MI_SNAPSHOT_OK);
  expect_true("ready flag", view.ready == 1);
  expect_true("first generation", view.generation == UINT64_C(1));
  expect_near("exact first epoch", view.epoch_s, 123.25, 0.0);
  expect_true("native dimensions",
              view.theta_points == nt && view.azimuth_points == np);
  expect_near("ionosphere radius", view.ionosphere_radius_re, 1.0, 0.0);
  expect_near("trace radius", view.trace_radius_re, 4.24795, 0.0);
  expect_near("theta pole node", view.colatitude_rad[0], 0.0, 0.0);
  expect_near("theta low-latitude node", view.colatitude_rad[nt - 1U],
              pi / 6.0, 2.0e-16);
  expect_near("phi=0 is first/midnight node", view.longitude_rad[0], 0.0,
              0.0);
  expect_near("phi quarter-turn node", view.longitude_rad[np / 4U],
              pi / 2.0, 2.0e-16);
  expect_near("theta-major potential layout",
              view.potential_v[2U * np + 3U], 1203.0, 0.0);
  const double *const persistent_pointer = view.potential_v;

  memcpy(invalid, first, count * sizeof(double));
  invalid[7] = NAN;
  expect_true("invalid/failed update rejected",
              gamera_mi_snapshot_test_publish(124.0, invalid) != 0);
  expect_true("failed update leaves publication available",
              gamera_mi_coupling_borrow_north_snapshot(1U, &view) ==
                  GAMERA_MI_SNAPSHOT_OK);
  expect_true("failed update does not advance generation",
              view.generation == UINT64_C(1));
  expect_near("failed update does not advance epoch", view.epoch_s, 123.25,
              0.0);
  expect_true("failed update does not alter data",
              memcmp(view.potential_v, first,
                     count * sizeof(double)) == 0);

  expect_true("nonfinite epoch rejected",
              gamera_mi_snapshot_test_publish(NAN, first) != 0);
  fill_potential(second, nt, np, 2000.0);
  expect_true("second successful publication",
              gamera_mi_snapshot_test_publish(180.0, second) == 0);
  memset(&view, 0xa5, sizeof(view));
  expect_true("old generation is explicitly stale",
              gamera_mi_coupling_borrow_north_snapshot(1U, &view) ==
                  GAMERA_MI_SNAPSHOT_STALE);
  expect_true("stale clears view", cleared_view(&view));
  expect_true("borrow exact second generation",
              gamera_mi_coupling_borrow_north_snapshot(2U, &view) ==
                  GAMERA_MI_SNAPSHOT_OK);
  expect_true("generation advances monotonically",
              view.generation == UINT64_C(2));
  expect_near("exact second epoch", view.epoch_s, 180.0, 0.0);
  expect_true("borrow is zero-copy persistent storage",
              view.potential_v == persistent_pointer);
  expect_true("second potential published exactly",
              memcmp(view.potential_v, second,
                     count * sizeof(double)) == 0);

  gamera_mi_snapshot_test_finalize();
  memset(&view, 0xa5, sizeof(view));
  expect_true("finalize makes snapshot unavailable",
              gamera_mi_coupling_borrow_north_snapshot(0U, &view) ==
                  GAMERA_MI_SNAPSHOT_UNAVAILABLE);
  expect_true("finalize clears returned view", cleared_view(&view));
}

static void test_status_strings(void) {
  expect_true("status ok string",
              strcmp(gamera_mi_snapshot_status_string(
                         GAMERA_MI_SNAPSHOT_OK),
                     "ok") == 0);
  expect_true("status stale string",
              strcmp(gamera_mi_snapshot_status_string(
                         GAMERA_MI_SNAPSHOT_STALE),
                     "stale") == 0);
}

int main(void) {
  test_snapshot_lifecycle();
  test_status_strings();
  if (failures != 0) {
    fprintf(stderr, "%d M-I snapshot check(s) failed\n", failures);
    return EXIT_FAILURE;
  }
  printf("all M-I snapshot checks passed\n");
  return EXIT_SUCCESS;
}
