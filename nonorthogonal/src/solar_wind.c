#include "solar_wind.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static int finite_state(const gamera_solar_wind_state *state) {
  if (state == NULL || !isfinite(state->density) ||
      !isfinite(state->pressure) || state->density <= 0.0 ||
      state->pressure <= 0.0) {
    return 0;
  }
  for (int component = 0; component < GAMERA_WIND_DIM; ++component) {
    if (!isfinite(state->velocity[component]) ||
        !isfinite(state->magnetic[component])) {
      return 0;
    }
  }
  /* The planar-front propagation map is singular when monitor Vx is zero. */
  return fabs(state->velocity[GAMERA_WIND_X]) > DBL_MIN;
}

void gamera_solar_wind_init(gamera_solar_wind_series *series) {
  if (series == NULL) {
    return;
  }
  memset(series, 0, sizeof(*series));
  series->linear_interpolation = 1;
}

void gamera_solar_wind_destroy(gamera_solar_wind_series *series) {
  if (series == NULL) {
    return;
  }
  free(series->time);
  free(series->state);
  gamera_solar_wind_init(series);
}

void gamera_solar_wind_apply_bx_relation(gamera_solar_wind_series *series) {
  if (series == NULL || !series->enforce_bx_relation) {
    return;
  }
  for (size_t sample = 0; sample < series->count; ++sample) {
    series->state[sample].magnetic[GAMERA_WIND_X] =
        series->bx_offset +
        series->by_coefficient *
            series->state[sample].magnetic[GAMERA_WIND_Y] +
        series->bz_coefficient *
            series->state[sample].magnetic[GAMERA_WIND_Z];
  }
}

int gamera_solar_wind_set(gamera_solar_wind_series *series, size_t count,
                          const double *time,
                          const gamera_solar_wind_state *state) {
  if (series == NULL || count == 0U || time == NULL || state == NULL) {
    return -1;
  }
  for (size_t sample = 0; sample < count; ++sample) {
    if (!isfinite(time[sample]) || !finite_state(&state[sample]) ||
        (sample > 0U && !(time[sample] > time[sample - 1U]))) {
      return -1;
    }
  }
  double *new_time = (double *)malloc(count * sizeof(*new_time));
  gamera_solar_wind_state *new_state =
      (gamera_solar_wind_state *)malloc(count * sizeof(*new_state));
  if (new_time == NULL || new_state == NULL) {
    free(new_time);
    free(new_state);
    return -1;
  }
  memcpy(new_time, time, count * sizeof(*new_time));
  memcpy(new_state, state, count * sizeof(*new_state));
  free(series->time);
  free(series->state);
  series->time = new_time;
  series->state = new_state;
  series->count = count;
  gamera_solar_wind_apply_bx_relation(series);
  return 0;
}

static size_t lower_sample(const gamera_solar_wind_series *series,
                           double time) {
  size_t lower = 0U;
  size_t upper = series->count - 1U;
  while (upper - lower > 1U) {
    const size_t middle = lower + (upper - lower) / 2U;
    if (series->time[middle] <= time) {
      lower = middle;
    } else {
      upper = middle;
    }
  }
  return lower;
}

static gamera_solar_wind_state interpolate(
    const gamera_solar_wind_state *left,
    const gamera_solar_wind_state *right, double fraction) {
  gamera_solar_wind_state result;
  result.density = left->density + fraction * (right->density - left->density);
  result.pressure =
      left->pressure + fraction * (right->pressure - left->pressure);
  for (int component = 0; component < GAMERA_WIND_DIM; ++component) {
    result.velocity[component] =
        left->velocity[component] +
        fraction * (right->velocity[component] - left->velocity[component]);
    result.magnetic[component] =
        left->magnetic[component] +
        fraction * (right->magnetic[component] - left->magnetic[component]);
  }
  return result;
}

int gamera_solar_wind_sample_time(const gamera_solar_wind_series *series,
                                  double time,
                                  gamera_solar_wind_state *state) {
  if (series == NULL || series->count == 0U || series->time == NULL ||
      series->state == NULL || state == NULL || !isfinite(time)) {
    return -1;
  }
  if (series->count == 1U || time <= series->time[0]) {
    *state = series->state[0];
    return 0;
  }
  if (time >= series->time[series->count - 1U]) {
    *state = series->state[series->count - 1U];
    return 0;
  }
  const size_t lower = lower_sample(series, time);
  if (!series->linear_interpolation) {
    *state = series->state[lower];
    return 0;
  }
  const double fraction =
      (time - series->time[lower]) /
      (series->time[lower + 1U] - series->time[lower]);
  *state = interpolate(&series->state[lower], &series->state[lower + 1U],
                       fraction);
  return finite_state(state) ? 0 : -1;
}

int gamera_solar_wind_sample_at(const gamera_solar_wind_series *series,
                                const double point[GAMERA_WIND_DIM],
                                double simulation_time,
                                gamera_solar_wind_state *state,
                                double *delay) {
  if (series == NULL || point == NULL || state == NULL ||
      !isfinite(simulation_time)) {
    return -1;
  }
  const double monitor_time = simulation_time + series->time_offset;
  gamera_solar_wind_state current;
  if (gamera_solar_wind_sample_time(series, monitor_time, &current) != 0) {
    return -1;
  }
  const double coefficient_square =
      series->by_coefficient * series->by_coefficient +
      series->bz_coefficient * series->bz_coefficient;
  const double factor = current.velocity[GAMERA_WIND_X] /
                        (1.0 + coefficient_square);
  const double front_velocity[GAMERA_WIND_DIM] = {
      factor, -factor * series->by_coefficient,
      -factor * series->bz_coefficient};
  double speed_square = 0.0;
  double projection = 0.0;
  for (int component = 0; component < GAMERA_WIND_DIM; ++component) {
    if (!isfinite(point[component])) {
      return -1;
    }
    speed_square += front_velocity[component] * front_velocity[component];
    projection += (point[component] - series->reference[component]) *
                  front_velocity[component];
  }
  if (!(speed_square > DBL_MIN) || !isfinite(speed_square)) {
    return -1;
  }
  const double local_delay = projection / speed_square;
  if (delay != NULL) {
    *delay = local_delay;
  }
  return gamera_solar_wind_sample_time(series, monitor_time - local_delay,
                                       state);
}

int gamera_solar_wind_weight(const gamera_solar_wind_series *series,
                             const double point[GAMERA_WIND_DIM],
                             double simulation_time, double *weight) {
  if (series == NULL || point == NULL || weight == NULL) {
    return -1;
  }
  gamera_solar_wind_state state;
  if (gamera_solar_wind_sample_at(series, point, simulation_time, &state,
                                  NULL) != 0) {
    return -1;
  }
  double radius_square = 0.0;
  double speed_square = 0.0;
  double cosine = 0.0;
  for (int component = 0; component < GAMERA_WIND_DIM; ++component) {
    radius_square += point[component] * point[component];
    speed_square += state.velocity[component] * state.velocity[component];
    cosine -= point[component] * state.velocity[component];
  }
  if (!(radius_square > DBL_MIN) || !(speed_square > DBL_MIN)) {
    return -1;
  }
  cosine /= sqrt(radius_square * speed_square);
  cosine = fmax(-1.0, fmin(1.0, cosine));
  const double angle = acos(cosine);
  const double half_pi = 0.5 * acos(-1.0);
  if (angle <= half_pi) {
    *weight = 1.0;
  } else if (angle >= 2.0 * half_pi) {
    *weight = 0.0;
  } else {
    const double scaled = (angle - half_pi) / half_pi;
    *weight = exp(1.0 - 1.0 / (1.0 - scaled * scaled));
  }
  return 0;
}
