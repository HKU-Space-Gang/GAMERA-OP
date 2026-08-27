#include "nonorthogonal_absolute_cadence.h"

#include <float.h>
#include <math.h>
#include <stdio.h>

#define REQUIRE(condition)                                                     \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #condition);    \
      return 1;                                                                \
    }                                                                          \
  } while (0)

static int nearly_equal(double first, double second) {
  return fabs(first - second) <=
         64.0 * DBL_EPSILON *
             fmax(1.0, fmax(fabs(first), fabs(second)));
}

int main(void) {
  gamera_no_absolute_cadence cadence = {0};
  REQUIRE(gamera_no_absolute_cadence_init(10.0, 10.0, 0.0, 0,
                                          &cadence) == 0);
  REQUIRE(cadence.next_due_s == 10.0);

  double limited = 0.0;
  REQUIRE(gamera_no_absolute_cadence_limit_timestep(
              &cadence, 0.75, 1.0, &limited) == 0);
  REQUIRE(nearly_equal(limited, 0.25));

  int due = -1;
  double due_time_s = -1.0;
  REQUIRE(gamera_no_absolute_cadence_poll(
              &cadence, 0.9, &due, &due_time_s) == 0);
  REQUIRE(due == 0);
  REQUIRE(due_time_s == 0.0);
  REQUIRE(gamera_no_absolute_cadence_poll(
              &cadence, 1.0, &due, &due_time_s) == 0);
  REQUIRE(due == 1);
  REQUIRE(due_time_s == 10.0);

  /* A failed physics refresh simply omits commit and leaves the due event. */
  REQUIRE(cadence.next_due_s == 10.0);
  REQUIRE(gamera_no_absolute_cadence_poll(
              &cadence, 1.0, &due, &due_time_s) == 0);
  REQUIRE(due == 1 && due_time_s == 10.0);
  REQUIRE(gamera_no_absolute_cadence_commit(&cadence, due_time_s) == 0);
  REQUIRE(cadence.next_due_s == 20.0);

  /* Reproduce the 10800 -> 10920 s gate landing: floating addition is one
   * ULP short, but both a continuous cadence and a split stop canonicalize to
   * the same exact code time before poll/commit. */
  const double gate_norm_s = 63.71;
  const double gate_start_code = 10800.0 / gate_norm_s;
  const double gate_stop_code = 10920.0 / gate_norm_s;
  gamera_no_absolute_cadence gate_cadence = {0};
  REQUIRE(gamera_no_absolute_cadence_init(
              120.0, gate_norm_s, gate_start_code, 1,
              &gate_cadence) == 0);
  REQUIRE(gamera_no_absolute_cadence_limit_timestep(
              &gate_cadence, gate_start_code, 10.0, &limited) == 0);
  const double raw_gate_end = gate_start_code + limited;
  REQUIRE(raw_gate_end < gate_stop_code);
  double continuous_end = -1.0;
  double split_end = -1.0;
  REQUIRE(gamera_no_absolute_cadence_canonicalize_step_end(
              &gate_cadence, raw_gate_end, gate_stop_code + 10.0,
              &continuous_end) == 0);
  REQUIRE(gamera_no_absolute_cadence_canonicalize_step_end(
              &gate_cadence, raw_gate_end, gate_stop_code,
              &split_end) == 0);
  REQUIRE(continuous_end == gate_stop_code && split_end == gate_stop_code);
  double exact_end = -1.0;
  double successor_end = -1.0;
  REQUIRE(gamera_no_absolute_cadence_canonicalize_step_end(
              &gate_cadence, gate_stop_code, gate_stop_code + 10.0,
              &exact_end) == 0);
  REQUIRE(gamera_no_absolute_cadence_canonicalize_step_end(
              &gate_cadence, nextafter(gate_stop_code, INFINITY),
              gate_stop_code + 10.0, &successor_end) == 0);
  REQUIRE(exact_end == gate_stop_code && successor_end == gate_stop_code);
  double successor_split_end = -1.0;
  REQUIRE(gamera_no_absolute_cadence_canonicalize_step_end(
              &gate_cadence, nextafter(gate_stop_code, INFINITY),
              gate_stop_code, &successor_split_end) == 0);
  REQUIRE(successor_split_end == gate_stop_code);
  double lower_stop_end = -1.0;
  double upper_stop_end = -1.0;
  REQUIRE(gamera_no_absolute_cadence_canonicalize_step_end(
              &gate_cadence, raw_gate_end,
              nextafter(gate_stop_code, -INFINITY),
              &lower_stop_end) == 0);
  REQUIRE(gamera_no_absolute_cadence_canonicalize_step_end(
              &gate_cadence, raw_gate_end,
              nextafter(gate_stop_code, INFINITY),
              &upper_stop_end) == 0);
  REQUIRE(lower_stop_end == gate_stop_code &&
          upper_stop_end == gate_stop_code);
  REQUIRE(gamera_no_absolute_cadence_poll(
              &gate_cadence, continuous_end, &due, &due_time_s) == 0);
  REQUIRE(due == 1 && due_time_s == 10920.0);
  REQUIRE(gamera_no_absolute_cadence_commit(&gate_cadence,
                                             due_time_s) == 0);
  REQUIRE(gamera_no_absolute_cadence_limit_timestep(
              &gate_cadence, continuous_end, 1.0, &limited) == 0);
  REQUIRE(limited == 1.0);

  double unchanged_end = -1.0;
  REQUIRE(gamera_no_absolute_cadence_init(
              120.0, gate_norm_s, gate_start_code, 1,
              &gate_cadence) == 0);
  REQUIRE(gamera_no_absolute_cadence_canonicalize_step_end(
              &gate_cadence, gate_stop_code - 1.0e-9,
              gate_stop_code + 10.0, &unchanged_end) == 0);
  REQUIRE(unchanged_end == gate_stop_code - 1.0e-9);
  REQUIRE(gamera_no_absolute_cadence_canonicalize_step_end(
              &gate_cadence, gate_stop_code + 1.0e-9,
              gate_stop_code + 10.0, &unchanged_end) != 0);

  /* Crossing a cadence boundary is rejected instead of silently skipping it. */
  REQUIRE(gamera_no_absolute_cadence_poll(
              &cadence, 2.01, &due, &due_time_s) != 0);
  REQUIRE(gamera_no_absolute_cadence_commit(&cadence, 19.0) != 0);
  REQUIRE(cadence.next_due_s == 20.0);

  int aligned = -1;
  REQUIRE(gamera_no_absolute_cadence_time_is_aligned(
              2.0, 10.0, 10.0, &aligned) == 0);
  REQUIRE(aligned == 1);
  REQUIRE(gamera_no_absolute_cadence_time_is_aligned(
              1.5, 10.0, 10.0, &aligned) == 0);
  REQUIRE(aligned == 0);
  REQUIRE(gamera_no_absolute_cadence_init(10.0, 10.0, 1.5, 1,
                                          &cadence) != 0);
  REQUIRE(gamera_no_absolute_cadence_init(10.0, 10.0, 2.0, 1,
                                          &cadence) == 0);
  REQUIRE(cadence.next_due_s == 30.0);

  /* Values within floating-point roundoff of a boundary are aligned. */
  REQUIRE(gamera_no_absolute_cadence_init(
              120.0, 63.7, 240.0 / 63.7, 1, &cadence) == 0);
  REQUIRE(cadence.next_due_s == 360.0);

  double saved_interval = 0.0;
  double saved_norm = 0.0;
  double saved_next = 0.0;
  REQUIRE(gamera_no_absolute_cadence_export(
              &cadence, &saved_interval, &saved_norm, &saved_next) == 0);
  REQUIRE(saved_interval == 120.0 && saved_norm == 63.7 &&
          saved_next == 360.0);
  gamera_no_absolute_cadence restored = {1.0, 1.0, 1.0};
  REQUIRE(gamera_no_absolute_cadence_restore(
              saved_interval, saved_norm, saved_next, 300.0 / 63.7,
              &restored) == 0);
  REQUIRE(restored.interval_s == cadence.interval_s &&
          restored.time_norm_s == cadence.time_norm_s &&
          restored.next_due_s == cadence.next_due_s);
  REQUIRE(gamera_no_absolute_cadence_restore(
              saved_interval, saved_norm, 350.0, 300.0 / 63.7,
              &restored) != 0);
  REQUIRE(restored.next_due_s == 360.0);
  REQUIRE(gamera_no_absolute_cadence_restore(
              saved_interval, saved_norm, saved_next, 360.0 / 63.7,
              &restored) != 0);
  REQUIRE(restored.next_due_s == 360.0);

  REQUIRE(gamera_no_absolute_cadence_init(0.0, 1.0, 0.0, 0,
                                          &cadence) != 0);
  REQUIRE(gamera_no_absolute_cadence_init(1.0, INFINITY, 0.0, 0,
                                          &cadence) != 0);
  REQUIRE(gamera_no_absolute_cadence_limit_timestep(
              &cadence, 0.0, NAN, &limited) != 0);

  puts("PASS target-independent absolute coupling cadence");
  return 0;
}
