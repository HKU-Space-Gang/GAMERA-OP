#include "nonorthogonal_mi_dpb.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static const double DPB_PI = 3.141592653589793238462643383279502884;

static double clamp_value(double value, double lower, double upper) {
  return fmax(lower, fmin(upper, value));
}

static int compare_double(const void *left, const void *right) {
  const double a = *(const double *)left;
  const double b = *(const double *)right;
  return (a > b) - (a < b);
}

static double percentile_copy(const double *values, size_t count,
                              double fraction, double *scratch) {
  if (values == NULL || scratch == NULL || count == 0U) {
    return NAN;
  }
  memcpy(scratch, values, count * sizeof(*scratch));
  qsort(scratch, count, sizeof(*scratch), compare_double);
  const double logical = clamp_value(fraction, 0.0, 1.0) * (double)(count - 1U);
  const size_t lower = (size_t)floor(logical);
  const size_t upper = lower + 1U < count ? lower + 1U : lower;
  const double weight = logical - (double)lower;
  return (1.0 - weight) * scratch[lower] + weight * scratch[upper];
}

static double latitude_deg(size_t theta, size_t theta_count,
                           double maximum_colatitude_rad) {
  const double colatitude = maximum_colatitude_rad * (double)theta /
                            (double)(theta_count - 1U);
  return 90.0 - 180.0 * colatitude / DPB_PI;
}

static double mlt_h(size_t longitude, size_t longitude_count) {
  return 24.0 * (double)longitude / (double)longitude_count;
}

static size_t index2(size_t theta, size_t longitude,
                     size_t longitude_count) {
  return theta * longitude_count + longitude;
}

static int valid_config(const gamera_mi_hybrid_dpb_config *config) {
  return config != NULL && isfinite(config->absolute_potential_v) &&
         config->absolute_potential_v > 0.0 &&
         isfinite(config->adaptive_potential_fraction) &&
         config->adaptive_potential_fraction > 0.0 &&
         config->adaptive_potential_fraction <= 1.0 &&
         isfinite(config->adaptive_potential_floor_v) &&
         config->adaptive_potential_floor_v > 0.0 &&
         isfinite(config->potential_equatorward_offset_deg) &&
         config->potential_equatorward_offset_deg >= 0.0 &&
         isfinite(config->fac_equatorward_offset_deg) &&
         config->fac_equatorward_offset_deg >= 0.0 &&
         isfinite(config->fac_absolute_floor_a_m2) &&
         config->fac_absolute_floor_a_m2 > 0.0 &&
         isfinite(config->minimum_boundary_latitude_deg) &&
         isfinite(config->maximum_candidate_latitude_deg) &&
         config->minimum_boundary_latitude_deg <
             config->maximum_candidate_latitude_deg &&
         isfinite(config->transition_width_deg) &&
         config->transition_width_deg > 0.0 &&
         isfinite(config->temporal_timescale_s) &&
         config->temporal_timescale_s > 0.0 &&
         isfinite(config->maximum_slew_deg_per_s) &&
         config->maximum_slew_deg_per_s > 0.0 &&
         isfinite(config->quiet_initial_nightside_latitude_deg) &&
         isfinite(config->oval_center_latitude_deg) &&
         isfinite(config->oval_center_mlt_h) &&
         isfinite(config->nightside_poleward_limit_deg) &&
         config->minimum_boundary_latitude_deg <
             config->nightside_poleward_limit_deg &&
         config->nightside_poleward_limit_deg <
             config->oval_center_latitude_deg;
}

gamera_mi_hybrid_dpb_config gamera_mi_hybrid_dpb_default_config(void) {
  const gamera_mi_hybrid_dpb_config config = {
      .absolute_potential_v = 10000.0,
      .adaptive_potential_fraction = 0.20,
      .adaptive_potential_floor_v = 2000.0,
      .potential_equatorward_offset_deg = 3.0,
      .fac_equatorward_offset_deg = 2.5,
      .fac_absolute_floor_a_m2 = 0.03e-6,
      .minimum_boundary_latitude_deg = 58.0,
      .maximum_candidate_latitude_deg = 74.0,
      .transition_width_deg = 3.0,
      .temporal_timescale_s = 300.0,
      .maximum_slew_deg_per_s = 1.0 / 120.0,
      .quiet_initial_nightside_latitude_deg = 64.5,
      .oval_center_latitude_deg = 85.0,
      .oval_center_mlt_h = 1.0,
      .nightside_poleward_limit_deg = 70.0};
  return config;
}

int gamera_mi_hybrid_dpb_boundary_from_radius(
    size_t longitude_count, double radius_deg,
    double center_latitude_deg, double center_mlt_h,
    double *boundary_latitude_deg) {
  if (longitude_count < 4U || !isfinite(radius_deg) || radius_deg <= 0.0 ||
      !isfinite(center_latitude_deg) || !isfinite(center_mlt_h) ||
      boundary_latitude_deg == NULL) {
    return -1;
  }
  const double center_colatitude =
      (90.0 - center_latitude_deg) * DPB_PI / 180.0;
  const double center_longitude = center_mlt_h * DPB_PI / 12.0;
  const double radius = radius_deg * DPB_PI / 180.0;
  for (size_t j = 0U; j < longitude_count; ++j) {
    const double longitude = 2.0 * DPB_PI * (double)j /
                             (double)longitude_count;
    const double delta = longitude - center_longitude;
    const double coefficient_cos = cos(center_colatitude);
    const double coefficient_sin =
        sin(center_colatitude) * cos(delta);
    const double amplitude = hypot(coefficient_cos, coefficient_sin);
    if (!(amplitude > DBL_MIN)) {
      return -1;
    }
    const double phase = atan2(coefficient_sin, coefficient_cos);
    const double argument = clamp_value(cos(radius) / amplitude, -1.0, 1.0);
    const double colatitude = phase + acos(argument);
    boundary_latitude_deg[j] = 90.0 - 180.0 * colatitude / DPB_PI;
    if (!isfinite(boundary_latitude_deg[j])) {
      return -1;
    }
  }
  return 0;
}

static double sector_median(const double *boundary, size_t longitude_count,
                            int nightside, double *scratch) {
  size_t count = 0U;
  for (size_t j = 0U; j < longitude_count; ++j) {
    const double mlt = mlt_h(j, longitude_count);
    const int selected = nightside ? (mlt >= 20.0 || mlt <= 4.0)
                                   : (mlt >= 10.0 && mlt <= 14.0);
    if (selected && isfinite(boundary[j])) {
      scratch[count++] = boundary[j];
    }
  }
  if (count == 0U) {
    return NAN;
  }
  qsort(scratch, count, sizeof(*scratch), compare_double);
  return count % 2U == 0U
             ? 0.5 * (scratch[count / 2U - 1U] + scratch[count / 2U])
             : scratch[count / 2U];
}

static double radius_for_nightside_limit(
    const gamera_mi_hybrid_dpb_config *config, size_t longitude_count,
    double *boundary, double *scratch) {
  const double center_colatitude = 90.0 - config->oval_center_latitude_deg;
  double lower = center_colatitude + 0.25;
  double upper = 89.0;
  for (int iteration = 0; iteration < 64; ++iteration) {
    const double middle = 0.5 * (lower + upper);
    gamera_mi_hybrid_dpb_boundary_from_radius(
        longitude_count, middle, config->oval_center_latitude_deg,
        config->oval_center_mlt_h, boundary);
    const double nightside =
        sector_median(boundary, longitude_count, 1, scratch);
    if (nightside > config->nightside_poleward_limit_deg) {
      lower = middle;
    } else {
      upper = middle;
    }
  }
  return upper;
}

static void smooth_periodic_3x3(const double *input, size_t theta_count,
                                size_t longitude_count, double *output) {
  static const double weight[3] = {1.0, 2.0, 1.0};
  for (size_t i = 0U; i < theta_count; ++i) {
    for (size_t j = 0U; j < longitude_count; ++j) {
      double sum = 0.0;
      double total = 0.0;
      for (int di = -1; di <= 1; ++di) {
        size_t ii;
        if (di < 0 && i == 0U) {
          ii = 0U;
        } else if (di > 0 && i + 1U == theta_count) {
          ii = theta_count - 1U;
        } else {
          ii = (size_t)((long long)i + di);
        }
        for (int dj = -1; dj <= 1; ++dj) {
          const size_t jj = (size_t)((long long)j + dj +
                                     (long long)longitude_count) %
                            longitude_count;
          const double w = weight[di + 1] * weight[dj + 1];
          sum += w * input[index2(ii, jj, longitude_count)];
          total += w;
        }
      }
      output[index2(i, j, longitude_count)] = sum / total;
    }
  }
}

static double circle_radius(double longitude_rad, double boundary_deg,
                            double center_latitude_deg,
                            double center_mlt_h) {
  const double latitude = boundary_deg * DPB_PI / 180.0;
  const double center_latitude = center_latitude_deg * DPB_PI / 180.0;
  const double center_longitude = center_mlt_h * DPB_PI / 12.0;
  const double cosine =
      sin(center_latitude) * sin(latitude) +
      cos(center_latitude) * cos(latitude) *
          cos(longitude_rad - center_longitude);
  return 180.0 * acos(clamp_value(cosine, -1.0, 1.0)) / DPB_PI;
}

static double robust_weighted_location(const double *value,
                                       const double *weight, size_t count,
                                       double *scratch) {
  size_t valid = 0U;
  for (size_t i = 0U; i < count; ++i) {
    if (isfinite(value[i]) && isfinite(weight[i]) && weight[i] > 0.0) {
      scratch[valid++] = value[i];
    }
  }
  if (valid < 5U) {
    return NAN;
  }
  qsort(scratch, valid, sizeof(*scratch), compare_double);
  double location = scratch[valid / 2U];
  for (int iteration = 0; iteration < 8; ++iteration) {
    size_t residual_count = 0U;
    for (size_t i = 0U; i < count; ++i) {
      if (isfinite(value[i]) && isfinite(weight[i]) && weight[i] > 0.0) {
        scratch[residual_count++] = fabs(value[i] - location);
      }
    }
    qsort(scratch, residual_count, sizeof(*scratch), compare_double);
    const double scale = 1.4826 * scratch[residual_count / 2U];
    if (!isfinite(scale) || scale < 0.05) {
      break;
    }
    double numerator = 0.0;
    double denominator = 0.0;
    for (size_t i = 0U; i < count; ++i) {
      if (isfinite(value[i]) && isfinite(weight[i]) && weight[i] > 0.0) {
        const double normalized = fabs(value[i] - location) / (1.5 * scale);
        const double robust = normalized <= 1.0 ? 1.0 : 1.0 / normalized;
        numerator += weight[i] * robust * value[i];
        denominator += weight[i] * robust;
      }
    }
    if (!(denominator > 0.0)) {
      break;
    }
    location = numerator / denominator;
  }
  return location;
}

int gamera_mi_hybrid_dpb_update(
    const gamera_mi_hybrid_dpb_config *config, size_t theta_count,
    size_t longitude_count, double maximum_colatitude_rad,
    double mapped_maximum_colatitude_rad, const double *fac_a_m2,
    double upward_fac_multiplier, const double *potential_v,
    double elapsed_s, gamera_mi_hybrid_dpb_state *state,
    double *boundary_latitude_deg, double *mask,
    gamera_mi_hybrid_dpb_stats *stats) {
  if (!valid_config(config) || theta_count < 5U || longitude_count < 8U ||
      !isfinite(maximum_colatitude_rad) || maximum_colatitude_rad <= 0.0 ||
      !isfinite(mapped_maximum_colatitude_rad) ||
      mapped_maximum_colatitude_rad <= 0.0 ||
      !isfinite(upward_fac_multiplier) ||
      fabs(fabs(upward_fac_multiplier) - 1.0) > 16.0 * DBL_EPSILON ||
      fac_a_m2 == NULL || potential_v == NULL || !isfinite(elapsed_s) ||
      elapsed_s <= 0.0 || state == NULL || boundary_latitude_deg == NULL ||
      mask == NULL || theta_count > SIZE_MAX / longitude_count) {
    return -1;
  }
  const size_t field_count = theta_count * longitude_count;
  /* Full-field storage: centered potential, smoothed amplitude, smoothed
   * upward FAC, gradient scratch, and FAC-statistics scratch.  The latter
   * must be field-sized because it receives one sample per selected
   * latitude and longitude, not merely one value per longitude. */
  if (field_count > SIZE_MAX / (5U * sizeof(double)) ||
      longitude_count > SIZE_MAX / (4U * sizeof(double))) {
    return -1;
  }
  double *field = (double *)calloc(5U * field_count, sizeof(double));
  double *column = (double *)calloc(4U * longitude_count, sizeof(double));
  if (field == NULL || column == NULL) {
    free(field);
    free(column);
    return -1;
  }
  double *centered = field;
  double *smooth_potential = centered + field_count;
  double *smooth_fac = smooth_potential + field_count;
  double *sample = smooth_fac + field_count;
  double *statistics_scratch = sample + field_count;
  double *candidate_radius = column;
  double *candidate_weight = candidate_radius + longitude_count;
  double *candidate_boundary = candidate_weight + longitude_count;
  double *geometry_scratch = candidate_boundary + longitude_count;

  size_t mapped_count = 0U;
  for (size_t i = 0U; i < theta_count; ++i) {
    const double theta = maximum_colatitude_rad * (double)i /
                         (double)(theta_count - 1U);
    if (theta <= mapped_maximum_colatitude_rad + 4.0e-14) {
      for (size_t j = 0U; j < longitude_count; ++j) {
        const size_t index = index2(i, j, longitude_count);
        if (!isfinite(potential_v[index]) || !isfinite(fac_a_m2[index])) {
          free(field);
          free(column);
          return -1;
        }
        sample[mapped_count++] = potential_v[index];
      }
    }
  }
  if (mapped_count < longitude_count) {
    free(field);
    free(column);
    return -1;
  }
  const double low = percentile_copy(sample, mapped_count, 0.01,
                                     centered);
  const double high = percentile_copy(sample, mapped_count, 0.99,
                                      centered);
  const double gauge = 0.5 * (low + high);
  const double cpcp = high - low;
  const double adaptive_threshold = fmin(
      config->absolute_potential_v,
      fmax(config->adaptive_potential_floor_v,
           config->adaptive_potential_fraction * cpcp));
  for (size_t i = 0U; i < theta_count; ++i) {
    const double theta = maximum_colatitude_rad * (double)i /
                         (double)(theta_count - 1U);
    for (size_t j = 0U; j < longitude_count; ++j) {
      const size_t index = index2(i, j, longitude_count);
      const int mapped = theta <= mapped_maximum_colatitude_rad + 4.0e-14;
      centered[index] = mapped ? fabs(potential_v[index] - gauge) : 0.0;
      sample[index] = mapped ? upward_fac_multiplier * fac_a_m2[index] : 0.0;
    }
  }
  smooth_periodic_3x3(centered, theta_count, longitude_count,
                      smooth_potential);
  smooth_periodic_3x3(sample, theta_count, longitude_count, smooth_fac);

  size_t gradient_count = 0U;
  size_t fac_count = 0U;
  for (size_t i = 1U; i + 1U < theta_count; ++i) {
    const double lat = latitude_deg(i, theta_count, maximum_colatitude_rad);
    const double theta = maximum_colatitude_rad * (double)i /
                         (double)(theta_count - 1U);
    if (lat < 59.0 || lat > 82.0 ||
        theta > mapped_maximum_colatitude_rad + 4.0e-14) {
      continue;
    }
    const double latitude_span =
        latitude_deg(i - 1U, theta_count, maximum_colatitude_rad) -
        latitude_deg(i + 1U, theta_count, maximum_colatitude_rad);
    for (size_t j = 0U; j < longitude_count; ++j) {
      sample[gradient_count++] = fmax(
          0.0,
          (smooth_potential[index2(i - 1U, j, longitude_count)] -
           smooth_potential[index2(i + 1U, j, longitude_count)]) /
              latitude_span);
      statistics_scratch[fac_count++] =
          fabs(smooth_fac[index2(i, j, longitude_count)]);
    }
  }
  const double gradient_scale = gradient_count > 0U
                                    ? percentile_copy(sample, gradient_count,
                                                      0.95, centered)
                                    : 0.0;
  const double fac_scale = fac_count > 0U
                               ? percentile_copy(statistics_scratch, fac_count,
                                                 0.95, centered)
                               : 0.0;
  memset(candidate_radius, 0, longitude_count * sizeof(double));
  memset(candidate_weight, 0, longitude_count * sizeof(double));
  size_t evidence_count = 0U;
  for (size_t j = 0U; j < longitude_count; ++j) {
    double potential_boundary = NAN;
    double potential_weight = 0.0;
    double best_gradient = 0.0;
    double local_gradient_max = 0.0;
    double local_potential_max = 0.0;
    for (size_t i = 1U; i + 1U < theta_count; ++i) {
      const double lat = latitude_deg(i, theta_count, maximum_colatitude_rad);
      const double theta = maximum_colatitude_rad * (double)i /
                           (double)(theta_count - 1U);
      if (lat < 59.0 || lat > 82.0 ||
          theta > mapped_maximum_colatitude_rad + 4.0e-14) {
        continue;
      }
      const double latitude_span =
          latitude_deg(i - 1U, theta_count, maximum_colatitude_rad) -
          latitude_deg(i + 1U, theta_count, maximum_colatitude_rad);
      const double gradient = fmax(
          0.0,
          (smooth_potential[index2(i - 1U, j, longitude_count)] -
           smooth_potential[index2(i + 1U, j, longitude_count)]) /
              latitude_span);
      local_gradient_max = fmax(local_gradient_max, gradient);
      local_potential_max = fmax(
          local_potential_max,
          smooth_potential[index2(i, j, longitude_count)]);
    }
    /* Scan equatorward to poleward (large theta to small theta) and keep the
     * first strong local potential-gradient peak. */
    for (size_t reverse = theta_count - 2U; reverse > 0U; --reverse) {
      const size_t i = reverse;
      const double lat = latitude_deg(i, theta_count, maximum_colatitude_rad);
      const double theta = maximum_colatitude_rad * (double)i /
                           (double)(theta_count - 1U);
      if (lat < 59.0 || lat > 82.0 ||
          theta > mapped_maximum_colatitude_rad + 4.0e-14) {
        continue;
      }
      const double latitude_span =
          latitude_deg(i - 1U, theta_count, maximum_colatitude_rad) -
          latitude_deg(i + 1U, theta_count, maximum_colatitude_rad);
      const double gradient = fmax(
          0.0,
          (smooth_potential[index2(i - 1U, j, longitude_count)] -
           smooth_potential[index2(i + 1U, j, longitude_count)]) /
              latitude_span);
      if (gradient >= 0.55 * local_gradient_max &&
          gradient >= 0.15 * gradient_scale &&
          local_potential_max >= config->adaptive_potential_floor_v) {
        potential_boundary = lat - config->potential_equatorward_offset_deg;
        best_gradient = gradient;
        potential_weight = fmin(1.0, gradient /
                                     fmax(gradient_scale, DBL_MIN));
        break;
      }
    }
    if (!isfinite(potential_boundary)) {
      for (size_t reverse = theta_count - 1U; reverse > 0U; --reverse) {
        const size_t i = reverse;
        const double lat =
            latitude_deg(i, theta_count, maximum_colatitude_rad);
        const double theta = maximum_colatitude_rad * (double)i /
                             (double)(theta_count - 1U);
        if (lat >= 59.0 && lat <= 82.0 &&
            theta <= mapped_maximum_colatitude_rad + 4.0e-14 &&
            smooth_potential[index2(i, j, longitude_count)] >=
                adaptive_threshold) {
          potential_boundary = lat;
          potential_weight = 0.35 * fmax(
              0.15, fmin(1.0, local_potential_max /
                                  adaptive_threshold - 1.0));
          break;
        }
      }
    }

    double fac_boundary = NAN;
    double fac_weight = 0.0;
    double local_fac_max = 0.0;
    for (size_t i = 1U; i + 1U < theta_count; ++i) {
      const double lat = latitude_deg(i, theta_count, maximum_colatitude_rad);
      const double theta = maximum_colatitude_rad * (double)i /
                           (double)(theta_count - 1U);
      if (lat >= 59.0 && lat <= 82.0 &&
          theta <= mapped_maximum_colatitude_rad + 4.0e-14) {
        local_fac_max = fmax(
            local_fac_max,
            fabs(smooth_fac[index2(i, j, longitude_count)]));
      }
    }
    const double fac_height = fmax(
        config->fac_absolute_floor_a_m2,
        fmax(0.18 * fac_scale, 0.22 * local_fac_max));
    size_t peak_index[32];
    size_t peak_count = 0U;
    for (size_t reverse = theta_count - 2U; reverse > 0U; --reverse) {
      const size_t i = reverse;
      const double lat = latitude_deg(i, theta_count, maximum_colatitude_rad);
      const double theta = maximum_colatitude_rad * (double)i /
                           (double)(theta_count - 1U);
      const double amplitude =
          fabs(smooth_fac[index2(i, j, longitude_count)]);
      if (lat >= 59.0 && lat <= 82.0 &&
          theta <= mapped_maximum_colatitude_rad + 4.0e-14 &&
          amplitude >= fac_height &&
          amplitude >= fabs(smooth_fac[index2(i - 1U, j, longitude_count)]) &&
          amplitude >= fabs(smooth_fac[index2(i + 1U, j, longitude_count)]) &&
          peak_count < sizeof(peak_index) / sizeof(peak_index[0])) {
        peak_index[peak_count++] = i;
      }
    }
    double best_pair_score = 0.0;
    size_t best_equatorward = 0U;
    for (size_t first = 0U; first < peak_count; ++first) {
      for (size_t second = first + 1U; second < peak_count; ++second) {
        const size_t lower = peak_index[first];
        const size_t upper = peak_index[second];
        const double lower_lat =
            latitude_deg(lower, theta_count, maximum_colatitude_rad);
        const double upper_lat =
            latitude_deg(upper, theta_count, maximum_colatitude_rad);
        const double separation = upper_lat - lower_lat;
        const double lower_value =
            smooth_fac[index2(lower, j, longitude_count)];
        const double upper_value =
            smooth_fac[index2(upper, j, longitude_count)];
        if (separation >= 2.0 && separation <= 15.0 &&
            lower_value * upper_value < 0.0) {
          const double strength = fmin(fabs(lower_value), fabs(upper_value));
          const double score = strength *
              exp(-0.5 * pow((separation - 7.0) / 4.0, 2.0));
          if (score > best_pair_score) {
            best_pair_score = score;
            best_equatorward = lower;
          }
        }
      }
    }
    if (best_pair_score > 0.0) {
      fac_boundary =
          latitude_deg(best_equatorward, theta_count,
                       maximum_colatitude_rad) -
          config->fac_equatorward_offset_deg;
      fac_weight = fmax(0.35, fmin(1.0, best_pair_score /
                                           fmax(fac_scale,
                                                config->fac_absolute_floor_a_m2)));
    }

    double weighted_boundary = 0.0;
    double total_weight = 0.0;
    if (isfinite(potential_boundary)) {
      weighted_boundary += potential_weight * potential_boundary;
      total_weight += potential_weight;
    }
    if (isfinite(fac_boundary) && fac_weight >= 0.35) {
      weighted_boundary += 2.0 * fac_weight * fac_boundary;
      total_weight += 2.0 * fac_weight;
    }
    if (total_weight > 0.0) {
      candidate_boundary[j] = clamp_value(
          weighted_boundary / total_weight,
          config->minimum_boundary_latitude_deg,
          config->maximum_candidate_latitude_deg);
      candidate_radius[j] = circle_radius(
          2.0 * DPB_PI * (double)j / (double)longitude_count,
          candidate_boundary[j], config->oval_center_latitude_deg,
          config->oval_center_mlt_h);
      candidate_weight[j] = total_weight;
      ++evidence_count;
    } else {
      candidate_boundary[j] = NAN;
      candidate_radius[j] = NAN;
      candidate_weight[j] = 0.0;
    }
    (void)best_gradient;
  }

  double target_radius = robust_weighted_location(
      candidate_radius, candidate_weight, longitude_count,
      statistics_scratch);
  const double center_colatitude = 90.0 - config->oval_center_latitude_deg;
  const double minimum_radius = radius_for_nightside_limit(
      config, longitude_count, boundary_latitude_deg, geometry_scratch);
  const double maximum_radius = fmax(
      minimum_radius,
      90.0 - config->minimum_boundary_latitude_deg - center_colatitude);
  if (!isfinite(target_radius)) {
    target_radius = state->initialized
                        ? state->radius_deg
                        : fabs(config->oval_center_latitude_deg -
                               config->quiet_initial_nightside_latitude_deg);
  }
  const double unconstrained_target = target_radius;
  target_radius = clamp_value(target_radius, minimum_radius, maximum_radius);
  double filtered_radius = target_radius;
  if (state->initialized) {
    const double alpha = 1.0 - exp(-elapsed_s / config->temporal_timescale_s);
    const double maximum_increment =
        config->maximum_slew_deg_per_s * elapsed_s;
    filtered_radius = state->radius_deg + clamp_value(
        alpha * (target_radius - state->radius_deg),
        -maximum_increment, maximum_increment);
  }
  if (gamera_mi_hybrid_dpb_boundary_from_radius(
          longitude_count, filtered_radius,
          config->oval_center_latitude_deg, config->oval_center_mlt_h,
          boundary_latitude_deg) != 0) {
    free(field);
    free(column);
    return -1;
  }
  if (state->initialized) {
    gamera_mi_hybrid_dpb_boundary_from_radius(
        longitude_count, state->radius_deg,
        config->oval_center_latitude_deg, config->oval_center_mlt_h,
        candidate_boundary);
    double maximum_step = 0.0;
    for (size_t j = 0U; j < longitude_count; ++j) {
      maximum_step = fmax(maximum_step,
                          fabs(boundary_latitude_deg[j] -
                               candidate_boundary[j]));
    }
    const double maximum_increment =
        config->maximum_slew_deg_per_s * elapsed_s;
    if (maximum_step > maximum_increment && maximum_step > 0.0) {
      const double fraction = maximum_increment / maximum_step;
      filtered_radius = state->radius_deg +
                        fraction * (filtered_radius - state->radius_deg);
      gamera_mi_hybrid_dpb_boundary_from_radius(
          longitude_count, filtered_radius,
          config->oval_center_latitude_deg, config->oval_center_mlt_h,
          boundary_latitude_deg);
    }
  }
  state->initialized = 1;
  state->radius_deg = filtered_radius;
  for (size_t i = 0U; i < theta_count; ++i) {
    const double latitude =
        latitude_deg(i, theta_count, maximum_colatitude_rad);
    const double theta = maximum_colatitude_rad * (double)i /
                         (double)(theta_count - 1U);
    for (size_t j = 0U; j < longitude_count; ++j) {
      const double coordinate = clamp_value(
          (latitude - boundary_latitude_deg[j]) /
              config->transition_width_deg,
          0.0, 1.0);
      mask[index2(i, j, longitude_count)] =
          theta <= mapped_maximum_colatitude_rad + 4.0e-14
              ? coordinate * coordinate * (3.0 - 2.0 * coordinate)
              : 0.0;
    }
  }
  if (stats != NULL) {
    stats->target_radius_deg = target_radius;
    stats->filtered_radius_deg = filtered_radius;
    stats->nightside_boundary_deg = sector_median(
        boundary_latitude_deg, longitude_count, 1, geometry_scratch);
    stats->dayside_boundary_deg = sector_median(
        boundary_latitude_deg, longitude_count, 0, geometry_scratch);
    stats->cpcp_v = cpcp;
    stats->adaptive_potential_threshold_v = adaptive_threshold;
    stats->fac_scale_a_m2 = fac_scale;
    stats->evidence_fraction =
        (double)evidence_count / (double)longitude_count;
    stats->nightside_limit_active =
        unconstrained_target < minimum_radius ? 1 : 0;
  }
  free(field);
  free(column);
  return 0;
}
