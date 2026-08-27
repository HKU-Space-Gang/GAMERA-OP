#include "nonorthogonal_absolute_cadence.h"

#include <float.h>
#include <math.h>
#include <stddef.h>

static int finite_positive(double value) {
  return isfinite(value) && value > 0.0;
}

static double cadence_tolerance(double time_s, double interval_s) {
  return 64.0 * DBL_EPSILON *
         fmax(1.0, fmax(fabs(time_s), fabs(interval_s)));
}

int gamera_no_absolute_cadence_time_is_aligned(
    double time_code, double time_norm_s, double interval_s, int *aligned) {
  if (aligned == NULL || !isfinite(time_code) || time_code < 0.0 ||
      !finite_positive(time_norm_s) || !finite_positive(interval_s)) {
    return -1;
  }
  const double time_s = time_code * time_norm_s;
  if (!isfinite(time_s)) {
    return -1;
  }
  const double nearest = nearbyint(time_s / interval_s) * interval_s;
  if (!isfinite(nearest)) {
    return -1;
  }
  *aligned =
      fabs(time_s - nearest) <= cadence_tolerance(time_s, interval_s);
  return 0;
}

int gamera_no_absolute_cadence_init(
    double interval_s, double time_norm_s, double initial_time_code,
    int require_alignment, gamera_no_absolute_cadence *cadence) {
  if (cadence == NULL || !finite_positive(interval_s) ||
      !finite_positive(time_norm_s) || !isfinite(initial_time_code) ||
      initial_time_code < 0.0 ||
      (require_alignment != 0 && require_alignment != 1)) {
    return -1;
  }
  int aligned = 0;
  if (gamera_no_absolute_cadence_time_is_aligned(
          initial_time_code, time_norm_s, interval_s, &aligned) != 0 ||
      (require_alignment && !aligned)) {
    return -1;
  }
  const double time_s = initial_time_code * time_norm_s;
  const double tolerance = cadence_tolerance(time_s, interval_s);
  const double next_due_s =
      interval_s * (floor((time_s + tolerance) / interval_s) + 1.0);
  if (!finite_positive(next_due_s) || next_due_s <= time_s) {
    return -1;
  }
  *cadence = (gamera_no_absolute_cadence){
      interval_s, time_norm_s, next_due_s};
  return 0;
}

int gamera_no_absolute_cadence_limit_timestep(
    const gamera_no_absolute_cadence *cadence, double time_code,
    double proposed_dt_code, double *limited_dt_code) {
  if (cadence == NULL || limited_dt_code == NULL ||
      !finite_positive(cadence->interval_s) ||
      !finite_positive(cadence->time_norm_s) ||
      !finite_positive(cadence->next_due_s) || !isfinite(time_code) ||
      time_code < 0.0 || !finite_positive(proposed_dt_code)) {
    return -1;
  }
  const double time_s = time_code * cadence->time_norm_s;
  if (!isfinite(time_s)) {
    return -1;
  }
  const double remaining_s = cadence->next_due_s - time_s;
  const double tolerance = cadence_tolerance(time_s, cadence->interval_s);
  if (remaining_s <= tolerance) {
    return -1;
  }
  const double landing_dt_code = remaining_s / cadence->time_norm_s;
  if (!finite_positive(landing_dt_code)) {
    return -1;
  }
  *limited_dt_code = fmin(proposed_dt_code, landing_dt_code);
  return 0;
}

int gamera_no_absolute_cadence_canonicalize_step_end(
    const gamera_no_absolute_cadence *cadence, double time_code,
    double stop_time_code, double *canonical_time_code) {
  if (cadence == NULL || canonical_time_code == NULL ||
      !finite_positive(cadence->interval_s) ||
      !finite_positive(cadence->time_norm_s) ||
      !finite_positive(cadence->next_due_s) || !isfinite(time_code) ||
      !isfinite(stop_time_code) || time_code < 0.0 ||
      stop_time_code < 0.0) {
    return -1;
  }
  const double initial_stop_tolerance =
      64.0 * DBL_EPSILON *
      fmax(1.0, fmax(fabs(time_code), fabs(stop_time_code)));
  if (time_code > stop_time_code + initial_stop_tolerance) {
    return -1;
  }
  double canonical = time_code;
  const double time_s = time_code * cadence->time_norm_s;
  const double due_code = cadence->next_due_s / cadence->time_norm_s;
  const double tolerance_s =
      cadence_tolerance(time_s, cadence->interval_s);
  if (!isfinite(time_s) || !isfinite(due_code) ||
      time_s > cadence->next_due_s + tolerance_s) {
    return -1;
  }
  const int cadence_snapped =
      fabs(time_s - cadence->next_due_s) <= tolerance_s;
  if (cadence_snapped) {
    canonical = due_code;
  }
  const double stop_tolerance =
      64.0 * DBL_EPSILON *
      fmax(1.0, fmax(fabs(canonical), fabs(stop_time_code)));
  if (cadence_snapped) {
    if (canonical > stop_time_code + stop_tolerance) {
      return -1;
    }
  } else if (fabs(stop_time_code - canonical) <= stop_tolerance) {
    canonical = stop_time_code;
  } else if (canonical > stop_time_code) {
    return -1;
  }
  *canonical_time_code = canonical;
  return 0;
}

int gamera_no_absolute_cadence_poll(
    const gamera_no_absolute_cadence *cadence, double time_code, int *due,
    double *due_time_s) {
  if (cadence == NULL || due == NULL || due_time_s == NULL ||
      !finite_positive(cadence->interval_s) ||
      !finite_positive(cadence->time_norm_s) ||
      !finite_positive(cadence->next_due_s) || !isfinite(time_code) ||
      time_code < 0.0) {
    return -1;
  }
  const double time_s = time_code * cadence->time_norm_s;
  if (!isfinite(time_s)) {
    return -1;
  }
  const double tolerance = cadence_tolerance(time_s, cadence->interval_s);
  if (time_s > cadence->next_due_s + tolerance) {
    return -1;
  }
  *due = time_s + tolerance >= cadence->next_due_s;
  *due_time_s = *due ? cadence->next_due_s : 0.0;
  return 0;
}

int gamera_no_absolute_cadence_commit(
    gamera_no_absolute_cadence *cadence, double due_time_s) {
  if (cadence == NULL || !finite_positive(cadence->interval_s) ||
      !finite_positive(cadence->time_norm_s) ||
      !finite_positive(cadence->next_due_s) || !isfinite(due_time_s)) {
    return -1;
  }
  const double tolerance =
      cadence_tolerance(due_time_s, cadence->interval_s);
  if (fabs(due_time_s - cadence->next_due_s) > tolerance) {
    return -1;
  }
  const double next_due_s = cadence->next_due_s + cadence->interval_s;
  if (!isfinite(next_due_s) || next_due_s <= cadence->next_due_s) {
    return -1;
  }
  cadence->next_due_s = next_due_s;
  return 0;
}

int gamera_no_absolute_cadence_export(
    const gamera_no_absolute_cadence *cadence, double *interval_s,
    double *time_norm_s, double *next_due_s) {
  if (cadence == NULL || interval_s == NULL || time_norm_s == NULL ||
      next_due_s == NULL || interval_s == time_norm_s ||
      interval_s == next_due_s || time_norm_s == next_due_s ||
      !finite_positive(cadence->interval_s) ||
      !finite_positive(cadence->time_norm_s) ||
      !finite_positive(cadence->next_due_s)) {
    return -1;
  }
  *interval_s = cadence->interval_s;
  *time_norm_s = cadence->time_norm_s;
  *next_due_s = cadence->next_due_s;
  return 0;
}

int gamera_no_absolute_cadence_restore(
    double interval_s, double time_norm_s, double next_due_s,
    double restart_time_code, gamera_no_absolute_cadence *cadence) {
  if (cadence == NULL || !finite_positive(interval_s) ||
      !finite_positive(time_norm_s) || !finite_positive(next_due_s) ||
      !isfinite(restart_time_code) || restart_time_code < 0.0) {
    return -1;
  }
  const double restart_time_s = restart_time_code * time_norm_s;
  if (!isfinite(restart_time_s) ||
      next_due_s - restart_time_s <=
          cadence_tolerance(restart_time_s, interval_s)) {
    return -1;
  }
  const double quotient = next_due_s / interval_s;
  const double nearest = nearbyint(quotient);
  if (!isfinite(quotient) || nearest < 1.0 ||
      fabs(next_due_s - nearest * interval_s) >
          cadence_tolerance(next_due_s, interval_s)) {
    return -1;
  }
  *cadence = (gamera_no_absolute_cadence){interval_s, time_norm_s,
                                          next_due_s};
  return 0;
}
