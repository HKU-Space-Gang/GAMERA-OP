#include "solar_wind.h"

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

static gamera_solar_wind_series make_series(void) {
  const double time[3] = {0.0, 10.0, 20.0};
  gamera_solar_wind_state state[3] = {0};
  for (int sample = 0; sample < 3; ++sample) {
    state[sample].density = 1.0 + 2.0 * sample;
    state[sample].pressure = 0.5 + sample;
    state[sample].velocity[0] = 2.0;
    state[sample].velocity[1] = 0.25;
    state[sample].velocity[2] = -0.5;
    state[sample].magnetic[0] = -99.0;
    state[sample].magnetic[1] = 1.0 + sample;
    state[sample].magnetic[2] = 0.5 * sample;
  }
  gamera_solar_wind_series result;
  gamera_solar_wind_init(&result);
  result.by_coefficient = 0.5;
  result.bz_coefficient = -0.25;
  result.bx_offset = 0.1;
  result.enforce_bx_relation = 1;
  if (gamera_solar_wind_set(&result, 3U, time, state) != 0) {
    fprintf(stderr, "FAIL valid time-series setup\n");
    ++failures;
  }
  return result;
}

static void test_time_sampling(void) {
  gamera_solar_wind_series series = make_series();
  gamera_solar_wind_state sample;
  if (gamera_solar_wind_sample_time(&series, 5.0, &sample) != 0) {
    fprintf(stderr, "FAIL linear time sample\n");
    ++failures;
  } else {
    expect_near("linear density", sample.density, 2.0, 1.0e-15);
    expect_near("linear pressure", sample.pressure, 1.0, 1.0e-15);
    expect_near("Lyon Bx relation", sample.magnetic[0], 0.7875, 1.0e-15);
  }
  series.linear_interpolation = 0;
  gamera_solar_wind_sample_time(&series, 9.99, &sample);
  expect_near("step density", sample.density, 1.0, 0.0);
  gamera_solar_wind_sample_time(&series, -100.0, &sample);
  expect_near("lower clamp", sample.density, 1.0, 0.0);
  gamera_solar_wind_sample_time(&series, 100.0, &sample);
  expect_near("upper clamp", sample.density, 5.0, 0.0);
  gamera_solar_wind_destroy(&series);
}

static void test_ballistic_mapping(void) {
  gamera_solar_wind_series series = make_series();
  series.linear_interpolation = 1;
  series.by_coefficient = 1.0;
  series.bz_coefficient = 0.0;
  series.reference[0] = 0.0;
  series.reference[1] = 0.0;
  series.reference[2] = 0.0;
  const double point[3] = {4.0, 2.0, 0.0};
  gamera_solar_wind_state sample;
  double delay = 0.0;
  if (gamera_solar_wind_sample_at(&series, point, 6.0, &sample, &delay) != 0) {
    fprintf(stderr, "FAIL ballistic sample\n");
    ++failures;
  } else {
    /* delay=(x-ay*y-az*z)/Vx=(4-2)/2=1 */
    expect_near("tilted-front delay", delay, 1.0, 2.0e-15);
    expect_near("tilted-front phase", sample.density, 2.0, 2.0e-15);
  }
  series.time_offset = 2.0;
  if (gamera_solar_wind_sample_at(&series, point, 4.0, &sample, &delay) != 0) {
    ++failures;
  } else {
    expect_near("monitor time offset", sample.density, 2.0, 2.0e-15);
  }
  gamera_solar_wind_destroy(&series);
}

static void test_boundary_weight(void) {
  gamera_solar_wind_series series = make_series();
  double weight = -1.0;
  const double upstream[3] = {-1.0, 0.0, 0.0};
  const double terminator[3] = {0.0, 1.0, 0.0};
  const double downstream[3] = {1.0, 0.0, 0.0};
  const double flank[3] = {1.0, 1.0, 0.0};
  /* Remove transverse velocity so these angles have exact reference values. */
  for (size_t sample = 0; sample < series.count; ++sample) {
    series.state[sample].velocity[1] = 0.0;
    series.state[sample].velocity[2] = 0.0;
  }
  gamera_solar_wind_weight(&series, upstream, 0.0, &weight);
  expect_near("upstream boundary weight", weight, 1.0, 0.0);
  gamera_solar_wind_weight(&series, terminator, 0.0, &weight);
  expect_near("terminator boundary weight", weight, 1.0, 0.0);
  gamera_solar_wind_weight(&series, downstream, 0.0, &weight);
  expect_near("downstream boundary weight", weight, 0.0, 0.0);
  gamera_solar_wind_weight(&series, flank, 0.0, &weight);
  expect_near("smooth tail ramp", weight, exp(-1.0 / 3.0), 2.0e-15);
  gamera_solar_wind_destroy(&series);
}

int main(void) {
  test_time_sampling();
  test_ballistic_mapping();
  test_boundary_weight();
  if (failures != 0) {
    fprintf(stderr, "%d solar-wind tests failed\n", failures);
    return 1;
  }
  printf("all solar-wind tests passed\n");
  return 0;
}
