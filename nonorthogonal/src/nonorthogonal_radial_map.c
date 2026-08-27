#include "nonorthogonal_radial_map.h"

#include <math.h>

static int solve_positive_exponential_stretch(double normalized_start_slope,
                                              double *stretch) {
  if (stretch == NULL || !isfinite(normalized_start_slope) ||
      !(normalized_start_slope > 0.0) ||
      !(normalized_start_slope < 1.0)) {
    return -1;
  }

  double lower = 0.0;
  double upper = 1.0;
  while (upper / expm1(upper) > normalized_start_slope) {
    upper *= 2.0;
    if (!(upper < 128.0)) {
      return -1;
    }
  }
  for (int iteration = 0; iteration < 100; ++iteration) {
    const double midpoint = 0.5 * (lower + upper);
    if (midpoint / expm1(midpoint) > normalized_start_slope) {
      lower = midpoint;
    } else {
      upper = midpoint;
    }
  }
  *stretch = 0.5 * (lower + upper);
  return isfinite(*stretch) && *stretch > 0.0 ? 0 : -1;
}

static int is_piecewise_earth_map(int version) {
  return version == GAMERA_NO_RADIAL_MAP_EARTH_PRODUCTION_V2 ||
         version == GAMERA_NO_RADIAL_MAP_EARTH_PRODUCTION_V3 ||
         version == GAMERA_NO_RADIAL_MAP_EARTH_PRODUCTION_V4;
}

int gamera_no_radial_map_init(gamera_no_radial_map *map, int version,
                              double inner_radius, double outer_radius,
                              double legacy_stretch) {
  if (map == NULL || !isfinite(inner_radius) || !isfinite(outer_radius) ||
      !(inner_radius > 0.0) || !(outer_radius > inner_radius)) {
    return -1;
  }

  *map = (gamera_no_radial_map){0};
  map->version = version;
  map->inner_radius = inner_radius;
  map->outer_radius = outer_radius;
  map->legacy_stretch = legacy_stretch;

  if (version == GAMERA_NO_RADIAL_MAP_LEGACY_EXPONENTIAL) {
    return isfinite(legacy_stretch) && legacy_stretch > 0.0 ? 0 : -1;
  }
  if (!is_piecewise_earth_map(version)) {
    return -1;
  }

  if (version == GAMERA_NO_RADIAL_MAP_EARTH_PRODUCTION_V4) {
    map->inner_logical_knot = 23.0 / 72.0;
    map->outer_logical_knot = 1.0;
    map->inner_physical_knot = 14.0;
    map->outer_physical_knot = outer_radius;
    if (!(inner_radius < map->inner_physical_knot) ||
        !(map->inner_physical_knot < outer_radius)) {
      return -1;
    }
    const double inner_derivative =
        (map->inner_physical_knot - inner_radius) /
        map->inner_logical_knot;
    const double outer_span = 1.0 - map->inner_logical_knot;
    const double outer_range = outer_radius - map->inner_physical_knot;
    const double normalized_start =
        inner_derivative * outer_span / outer_range;
    if (solve_positive_exponential_stretch(normalized_start,
                                           &map->middle_stretch) != 0) {
      return -1;
    }
    map->outer_stretch = 0.0;
    return 0;
  }

  if (version == GAMERA_NO_RADIAL_MAP_EARTH_PRODUCTION_V2) {
    map->inner_logical_knot = 11.0 / 32.0;
    map->outer_logical_knot = 25.0 / 32.0;
  } else {
    map->inner_logical_knot = 3.0 / 8.0;
    map->outer_logical_knot = 3.0 / 4.0;
  }
  map->inner_physical_knot = 15.0;
  map->outer_physical_knot = 50.0;
  if (!(inner_radius < map->inner_physical_knot) ||
      !(map->inner_physical_knot < map->outer_physical_knot) ||
      !(map->outer_physical_knot < outer_radius)) {
    return -1;
  }

  const double inner_span = map->inner_logical_knot;
  const double middle_span =
      map->outer_logical_knot - map->inner_logical_knot;
  const double outer_span = 1.0 - map->outer_logical_knot;
  const double inner_derivative =
      (map->inner_physical_knot - inner_radius) / inner_span;
  const double middle_range =
      map->outer_physical_knot - map->inner_physical_knot;
  const double middle_normalized_start =
      inner_derivative * middle_span / middle_range;
  if (solve_positive_exponential_stretch(middle_normalized_start,
                                         &map->middle_stretch) != 0) {
    return -1;
  }

  const double middle_end_derivative =
      middle_range / middle_span * map->middle_stretch *
      exp(map->middle_stretch) / expm1(map->middle_stretch);
  const double outer_range = outer_radius - map->outer_physical_knot;
  const double outer_normalized_start =
      middle_end_derivative * outer_span / outer_range;
  if (solve_positive_exponential_stretch(outer_normalized_start,
                                         &map->outer_stretch) != 0) {
    return -1;
  }
  return 0;
}

int gamera_no_radial_map_forward(const gamera_no_radial_map *map,
                                 double logical_fraction, double *radius) {
  if (map == NULL || radius == NULL || !isfinite(logical_fraction)) {
    return -1;
  }
  if (map->version == GAMERA_NO_RADIAL_MAP_LEGACY_EXPONENTIAL) {
    *radius = map->inner_radius +
              (map->outer_radius - map->inner_radius) *
                  expm1(map->legacy_stretch * logical_fraction) /
                  expm1(map->legacy_stretch);
  } else if (map->version == GAMERA_NO_RADIAL_MAP_EARTH_PRODUCTION_V4) {
    if (logical_fraction < 0.0) {
      const double inner_derivative =
          (map->inner_physical_knot - map->inner_radius) /
          map->inner_logical_knot;
      *radius = map->inner_radius *
                exp(inner_derivative / map->inner_radius *
                    logical_fraction);
    } else if (logical_fraction <= map->inner_logical_knot) {
      *radius = map->inner_radius +
                (map->inner_physical_knot - map->inner_radius) *
                    logical_fraction / map->inner_logical_knot;
    } else {
      const double fraction =
          (logical_fraction - map->inner_logical_knot) /
          (1.0 - map->inner_logical_knot);
      *radius = map->inner_physical_knot +
                (map->outer_radius - map->inner_physical_knot) *
                    expm1(map->middle_stretch * fraction) /
                    expm1(map->middle_stretch);
    }
  } else if (is_piecewise_earth_map(map->version)) {
    if (logical_fraction < 0.0) {
      const double inner_derivative =
          (map->inner_physical_knot - map->inner_radius) /
          map->inner_logical_knot;
      *radius = map->inner_radius *
                exp(inner_derivative / map->inner_radius *
                    logical_fraction);
    } else if (logical_fraction <= map->inner_logical_knot) {
      *radius = map->inner_radius +
                (map->inner_physical_knot - map->inner_radius) *
                    logical_fraction / map->inner_logical_knot;
    } else if (logical_fraction <= map->outer_logical_knot) {
      const double fraction =
          (logical_fraction - map->inner_logical_knot) /
          (map->outer_logical_knot - map->inner_logical_knot);
      *radius = map->inner_physical_knot +
                (map->outer_physical_knot - map->inner_physical_knot) *
                    expm1(map->middle_stretch * fraction) /
                    expm1(map->middle_stretch);
    } else {
      const double fraction =
          (logical_fraction - map->outer_logical_knot) /
          (1.0 - map->outer_logical_knot);
      *radius = map->outer_physical_knot +
                (map->outer_radius - map->outer_physical_knot) *
                    expm1(map->outer_stretch * fraction) /
                    expm1(map->outer_stretch);
    }
  } else {
    return -1;
  }
  return isfinite(*radius) ? 0 : -1;
}

int gamera_no_radial_map_inverse(const gamera_no_radial_map *map,
                                 double radius, double *logical_fraction) {
  if (map == NULL || logical_fraction == NULL || !isfinite(radius)) {
    return -1;
  }
  if (map->version == GAMERA_NO_RADIAL_MAP_LEGACY_EXPONENTIAL) {
    const double scaled =
        (radius - map->inner_radius) /
        (map->outer_radius - map->inner_radius) *
        expm1(map->legacy_stretch);
    if (!(scaled > -1.0)) {
      return -1;
    }
    *logical_fraction = log1p(scaled) / map->legacy_stretch;
  } else if (map->version == GAMERA_NO_RADIAL_MAP_EARTH_PRODUCTION_V4) {
    if (radius < map->inner_radius) {
      if (!(radius > 0.0)) {
        return -1;
      }
      const double inner_derivative =
          (map->inner_physical_knot - map->inner_radius) /
          map->inner_logical_knot;
      *logical_fraction =
          log(radius / map->inner_radius) * map->inner_radius /
          inner_derivative;
    } else if (radius <= map->inner_physical_knot) {
      *logical_fraction =
          (radius - map->inner_radius) /
          (map->inner_physical_knot - map->inner_radius) *
          map->inner_logical_knot;
    } else {
      const double scaled =
          (radius - map->inner_physical_knot) /
          (map->outer_radius - map->inner_physical_knot) *
          expm1(map->middle_stretch);
      if (!(scaled > -1.0)) {
        return -1;
      }
      *logical_fraction =
          map->inner_logical_knot +
          (1.0 - map->inner_logical_knot) *
              log1p(scaled) / map->middle_stretch;
    }
  } else if (is_piecewise_earth_map(map->version)) {
    if (radius < map->inner_radius) {
      if (!(radius > 0.0)) {
        return -1;
      }
      const double inner_derivative =
          (map->inner_physical_knot - map->inner_radius) /
          map->inner_logical_knot;
      *logical_fraction =
          log(radius / map->inner_radius) * map->inner_radius /
          inner_derivative;
    } else if (radius <= map->inner_physical_knot) {
      *logical_fraction =
          (radius - map->inner_radius) /
          (map->inner_physical_knot - map->inner_radius) *
          map->inner_logical_knot;
    } else if (radius <= map->outer_physical_knot) {
      const double scaled =
          (radius - map->inner_physical_knot) /
          (map->outer_physical_knot - map->inner_physical_knot) *
          expm1(map->middle_stretch);
      if (!(scaled > -1.0)) {
        return -1;
      }
      *logical_fraction =
          map->inner_logical_knot +
          (map->outer_logical_knot - map->inner_logical_knot) *
              log1p(scaled) / map->middle_stretch;
    } else {
      const double scaled =
          (radius - map->outer_physical_knot) /
          (map->outer_radius - map->outer_physical_knot) *
          expm1(map->outer_stretch);
      if (!(scaled > -1.0)) {
        return -1;
      }
      *logical_fraction =
          map->outer_logical_knot +
          (1.0 - map->outer_logical_knot) * log1p(scaled) /
              map->outer_stretch;
    }
  } else {
    return -1;
  }
  return isfinite(*logical_fraction) ? 0 : -1;
}

int gamera_no_radial_map_parameters(
    const gamera_no_radial_map *map,
    double parameters[GAMERA_NO_RADIAL_MAP_PARAMETER_COUNT]) {
  if (map == NULL || parameters == NULL) {
    return -1;
  }
  if (map->version == GAMERA_NO_RADIAL_MAP_LEGACY_EXPONENTIAL) {
    parameters[0] = map->legacy_stretch;
    for (size_t index = 1; index < GAMERA_NO_RADIAL_MAP_PARAMETER_COUNT;
         ++index) {
      parameters[index] = 0.0;
    }
  } else if (is_piecewise_earth_map(map->version)) {
    parameters[0] = map->inner_logical_knot;
    parameters[1] = map->outer_logical_knot;
    parameters[2] = map->inner_physical_knot;
    parameters[3] = map->outer_physical_knot;
    parameters[4] = map->middle_stretch;
    parameters[5] = map->outer_stretch;
  } else {
    return -1;
  }
  for (size_t index = 0; index < GAMERA_NO_RADIAL_MAP_PARAMETER_COUNT;
       ++index) {
    if (!isfinite(parameters[index])) {
      return -1;
    }
  }
  return 0;
}
