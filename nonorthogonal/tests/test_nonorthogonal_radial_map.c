#include "nonorthogonal_radial_map.h"

#include <math.h>
#include <stdio.h>

static int failures;

static void expect_near(const char *name, double actual, double expected,
                        double tolerance) {
  if (!isfinite(actual) || fabs(actual - expected) > tolerance) {
    fprintf(stderr,
            "FAIL %s: actual=%.17g expected=%.17g tolerance=%.3e\n",
            name, actual, expected, tolerance);
    ++failures;
  }
}

static void expect_true(const char *name, int condition) {
  if (!condition) {
    fprintf(stderr, "FAIL %s\n", name);
    ++failures;
  }
}

static double evaluate(const gamera_no_radial_map *map, double logical) {
  double radius = NAN;
  if (gamera_no_radial_map_forward(map, logical, &radius) != 0) {
    fprintf(stderr, "FAIL radial map evaluation at %.17g\n", logical);
    ++failures;
  }
  return radius;
}

static void test_legacy_exponential(void) {
  gamera_no_radial_map map;
  expect_true("legacy init",
              gamera_no_radial_map_init(
                  &map, GAMERA_NO_RADIAL_MAP_LEGACY_EXPONENTIAL, 4.0,
                  100.0, 3.0) == 0);
  for (int index = -4; index <= 68; ++index) {
    const double logical = (double)index / 64.0;
    const double expected =
        4.0 + 96.0 * expm1(3.0 * logical) / expm1(3.0);
    const double radius = evaluate(&map, logical);
    expect_near("legacy forward parity", radius, expected, 2.0e-13);
    double inverse = NAN;
    expect_true("legacy inverse status",
                gamera_no_radial_map_inverse(&map, radius, &inverse) == 0);
    expect_near("legacy inverse round trip", inverse, logical, 3.0e-14);
  }
}

static void test_earth_production_v2_profile(void) {
  gamera_no_radial_map map;
  expect_true("production init",
              gamera_no_radial_map_init(
                  &map, GAMERA_NO_RADIAL_MAP_EARTH_PRODUCTION, 4.0, 100.0,
                  3.0) == 0);
  expect_near("production inner endpoint", evaluate(&map, 0.0), 4.0,
              1.0e-14);
  expect_near("production 15 RE knot", evaluate(&map, 11.0 / 32.0), 15.0,
              2.0e-14);
  expect_near("production 50 RE knot", evaluate(&map, 25.0 / 32.0), 50.0,
              5.0e-14);
  expect_near("production outer endpoint", evaluate(&map, 1.0), 100.0,
              1.0e-13);
  expect_true("production inner halo stays positive",
              evaluate(&map, -4.0 / 64.0) > 0.0);
  expect_true("production coarse-grid inner halo stays positive",
              evaluate(&map, -4.0 / 24.0) > 0.0);

  double maximum_adjacent_ratio = 1.0;
  double previous_width = NAN;
  for (int index = 0; index < 64; ++index) {
    const double left = evaluate(&map, (double)index / 64.0);
    const double right = evaluate(&map, (double)(index + 1) / 64.0);
    const double width = right - left;
    expect_true("production width positive", width > 0.0);
    if (index < 22) {
      expect_near("production inner width", width, 0.5, 4.0e-14);
    }
    if (index > 0) {
      const double ratio = fmax(width / previous_width,
                                previous_width / width);
      maximum_adjacent_ratio = fmax(maximum_adjacent_ratio, ratio);
    }
    previous_width = width;
  }
  expect_true("production adjacent ratio bounded",
              maximum_adjacent_ratio < 1.061);
  expect_true("production stronger far-tail stretch",
              previous_width > 4.7 && previous_width < 4.9);

  for (int index = -4; index <= 68; ++index) {
    const double logical = (double)index / 64.0;
    const double radius = evaluate(&map, logical);
    double inverse = NAN;
    expect_true("production inverse status",
                gamera_no_radial_map_inverse(&map, radius, &inverse) == 0);
    expect_near("production inverse round trip", inverse, logical, 8.0e-14);
  }

  const double h = 1.0e-7;
  const double knots[] = {11.0 / 32.0, 25.0 / 32.0};
  for (size_t knot = 0; knot < sizeof(knots) / sizeof(knots[0]); ++knot) {
    const double center = evaluate(&map, knots[knot]);
    const double left_derivative =
        (center - evaluate(&map, knots[knot] - h)) / h;
    const double right_derivative =
        (evaluate(&map, knots[knot] + h) - center) / h;
    expect_near("production C1 join", right_derivative, left_derivative,
                8.0e-5);
  }

  double parameters[GAMERA_NO_RADIAL_MAP_PARAMETER_COUNT];
  expect_true("production identity status",
              gamera_no_radial_map_parameters(&map, parameters) == 0);
  expect_near("production identity inner logical knot", parameters[0],
              11.0 / 32.0, 0.0);
  expect_near("production identity outer logical knot", parameters[1],
              25.0 / 32.0, 0.0);
  expect_near("production identity inner physical knot", parameters[2],
              15.0, 0.0);
  expect_near("production identity outer physical knot", parameters[3],
              50.0, 0.0);
  expect_near("production middle stretch", parameters[4],
              1.618788125264683, 3.0e-14);
  expect_near("production outer stretch", parameters[5],
              0.6586298120651575, 3.0e-14);
}

static void test_earth_production_v3_profile(void) {
  gamera_no_radial_map map;
  expect_true("production v3 init",
              gamera_no_radial_map_init(
                  &map, GAMERA_NO_RADIAL_MAP_EARTH_PRODUCTION_V3, 3.0,
                  200.0, 3.0) == 0);
  expect_near("production v3 inner endpoint", evaluate(&map, 0.0), 3.0,
              1.0e-14);
  expect_near("production v3 15 RE knot", evaluate(&map, 3.0 / 8.0),
              15.0, 2.0e-14);
  expect_near("production v3 50 RE knot", evaluate(&map, 3.0 / 4.0),
              50.0, 8.0e-14);
  expect_near("production v3 outer endpoint", evaluate(&map, 1.0), 200.0,
              3.0e-13);
  expect_true("production v3 inner halo stays positive",
              evaluate(&map, -4.0 / 64.0) > 0.0);

  double maximum_adjacent_ratio = 1.0;
  double previous_width = NAN;
  for (int index = 0; index < 64; ++index) {
    const double left = evaluate(&map, (double)index / 64.0);
    const double right = evaluate(&map, (double)(index + 1) / 64.0);
    const double width = right - left;
    expect_true("production v3 width positive", width > 0.0);
    if (index < 24) {
      expect_near("production v3 inner width", width, 0.5, 5.0e-14);
    }
    if (index > 0) {
      const double ratio = fmax(width / previous_width,
                                previous_width / width);
      maximum_adjacent_ratio = fmax(maximum_adjacent_ratio, ratio);
    }
    previous_width = width;
  }
  expect_true("production v3 adjacent ratio bounded",
              maximum_adjacent_ratio < 1.124);
  expect_true("production v3 far-domain stretch",
              previous_width > 19.4 && previous_width < 19.6);

  for (int index = -4; index <= 68; ++index) {
    const double logical = (double)index / 64.0;
    const double radius = evaluate(&map, logical);
    double inverse = NAN;
    expect_true("production v3 inverse status",
                gamera_no_radial_map_inverse(&map, radius, &inverse) == 0);
    expect_near("production v3 inverse round trip", inverse, logical,
                1.2e-13);
  }

  const double h = 1.0e-7;
  const double knots[] = {3.0 / 8.0, 3.0 / 4.0};
  for (size_t knot = 0; knot < sizeof(knots) / sizeof(knots[0]); ++knot) {
    const double center = evaluate(&map, knots[knot]);
    const double left_derivative =
        (center - evaluate(&map, knots[knot] - h)) / h;
    const double right_derivative =
        (evaluate(&map, knots[knot] + h) - center) / h;
    expect_near("production v3 C1 join", right_derivative,
                left_derivative, 2.0e-4);
  }

  double parameters[GAMERA_NO_RADIAL_MAP_PARAMETER_COUNT];
  expect_true("production v3 identity status",
              gamera_no_radial_map_parameters(&map, parameters) == 0);
  expect_near("production v3 identity inner logical knot", parameters[0],
              3.0 / 8.0, 0.0);
  expect_near("production v3 identity outer logical knot", parameters[1],
              3.0 / 4.0, 0.0);
  expect_near("production v3 identity inner physical knot", parameters[2],
              15.0, 0.0);
  expect_near("production v3 identity outer physical knot", parameters[3],
              50.0, 0.0);
  expect_near("production v3 middle stretch", parameters[4],
              1.8603600760888819, 4.0e-14);
  expect_near("production v3 outer stretch", parameters[5],
              1.8609665578258845, 4.0e-14);
}

static void test_earth_production_v4_profile(void) {
  gamera_no_radial_map map;
  expect_true("production v4 init",
              gamera_no_radial_map_init(
                  &map, GAMERA_NO_RADIAL_MAP_EARTH_PRODUCTION_V4, 2.5,
                  200.0, 3.0) == 0);
  expect_near("production v4 inner endpoint", evaluate(&map, 0.0), 2.5,
              1.0e-14);
  expect_near("production v4 14 RE knot", evaluate(&map, 23.0 / 72.0),
              14.0, 3.0e-14);
  expect_near("production v4 outer endpoint", evaluate(&map, 1.0), 200.0,
              3.0e-13);
  expect_true("production v4 inner halo stays positive",
              evaluate(&map, -4.0 / 72.0) > 0.0);

  double maximum_adjacent_ratio = 1.0;
  double previous_width = NAN;
  for (int index = 0; index < 72; ++index) {
    const double left = evaluate(&map, (double)index / 72.0);
    const double right = evaluate(&map, (double)(index + 1) / 72.0);
    const double width = right - left;
    expect_true("production v4 width positive", width > 0.0);
    if (index < 23) {
      expect_near("production v4 inner width", width, 0.5, 7.0e-14);
    }
    if (index > 0) {
      const double ratio =
          fmax(width / previous_width, previous_width / width);
      maximum_adjacent_ratio = fmax(maximum_adjacent_ratio, ratio);
    }
    previous_width = width;
  }
  expect_true("production v4 adjacent ratio bounded",
              maximum_adjacent_ratio < 1.069);
  expect_true("production v4 far-domain stretch",
              previous_width > 12.3 && previous_width < 12.5);

  for (int index = -4; index <= 76; ++index) {
    const double logical = (double)index / 72.0;
    const double radius = evaluate(&map, logical);
    double inverse = NAN;
    expect_true("production v4 inverse status",
                gamera_no_radial_map_inverse(&map, radius, &inverse) == 0);
    expect_near("production v4 inverse round trip", inverse, logical,
                1.5e-13);
  }

  const double h = 1.0e-7;
  const double knot = 23.0 / 72.0;
  const double center = evaluate(&map, knot);
  const double left_derivative =
      (center - evaluate(&map, knot - h)) / h;
  const double right_derivative =
      (evaluate(&map, knot + h) - center) / h;
  expect_near("production v4 C1 join", right_derivative, left_derivative,
              2.0e-4);

  double parameters[GAMERA_NO_RADIAL_MAP_PARAMETER_COUNT];
  expect_true("production v4 identity status",
              gamera_no_radial_map_parameters(&map, parameters) == 0);
  expect_near("production v4 identity inner logical knot", parameters[0],
              23.0 / 72.0, 0.0);
  expect_near("production v4 identity outer logical knot", parameters[1],
              1.0, 0.0);
  expect_near("production v4 identity inner physical knot", parameters[2],
              14.0, 0.0);
  expect_near("production v4 identity outer physical knot", parameters[3],
              200.0, 0.0);
  expect_near("production v4 stretch", parameters[4],
              3.243548585696993, 5.0e-14);
  expect_near("production v4 unused outer stretch", parameters[5], 0.0,
              0.0);
}

static void test_invalid_profiles(void) {
  gamera_no_radial_map map;
  expect_true("unknown version rejected",
              gamera_no_radial_map_init(&map, 17, 4.0, 100.0, 3.0) != 0);
  expect_true("production outer below knot rejected",
              gamera_no_radial_map_init(
                  &map, GAMERA_NO_RADIAL_MAP_EARTH_PRODUCTION, 4.0, 45.0,
                  3.0) != 0);
  expect_true("legacy zero stretch rejected",
              gamera_no_radial_map_init(
                  &map, GAMERA_NO_RADIAL_MAP_LEGACY_EXPONENTIAL, 4.0,
                  100.0, 0.0) != 0);
  expect_true("production v3 outer below knot rejected",
              gamera_no_radial_map_init(
                  &map, GAMERA_NO_RADIAL_MAP_EARTH_PRODUCTION_V3, 3.0,
                  45.0, 3.0) != 0);
  expect_true("production v4 outer below knot rejected",
              gamera_no_radial_map_init(
                  &map, GAMERA_NO_RADIAL_MAP_EARTH_PRODUCTION_V4, 2.5,
                  13.5, 3.0) != 0);
}

int main(void) {
  test_legacy_exponential();
  test_earth_production_v2_profile();
  test_earth_production_v3_profile();
  test_earth_production_v4_profile();
  test_invalid_profiles();
  if (failures != 0) {
    fprintf(stderr, "%d radial-map assertions failed\n", failures);
    return 1;
  }
  printf("all radial-map tests passed\n");
  return 0;
}
