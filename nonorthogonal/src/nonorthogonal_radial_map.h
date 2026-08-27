#ifndef NONORTHOGONAL_RADIAL_MAP_H
#define NONORTHOGONAL_RADIAL_MAP_H

#include <stddef.h>

#define GAMERA_NO_RADIAL_MAP_LEGACY_EXPONENTIAL 1
#define GAMERA_NO_RADIAL_MAP_EARTH_PRODUCTION_V2 2
#define GAMERA_NO_RADIAL_MAP_EARTH_PRODUCTION_V3 3
#define GAMERA_NO_RADIAL_MAP_EARTH_PRODUCTION_V4 4
/* Source compatibility for archived v2 callers and restart tests. */
#define GAMERA_NO_RADIAL_MAP_EARTH_PRODUCTION \
  GAMERA_NO_RADIAL_MAP_EARTH_PRODUCTION_V2
#define GAMERA_NO_RADIAL_MAP_PARAMETER_COUNT 6U

typedef struct {
  int version;
  double inner_radius;
  double outer_radius;
  double legacy_stretch;
  double inner_logical_knot;
  double outer_logical_knot;
  double inner_physical_knot;
  double outer_physical_knot;
  double middle_stretch;
  double outer_stretch;
} gamera_no_radial_map;

/*
 * Version 1 is the historical single exponential.  Version 2 preserves the
 * accepted 4--100 RE Earth profile for restart compatibility.  Version 3 is
 * the 3--200 RE production profile: 3--15 RE is uniform, 15--50 RE is mildly
 * stretched, and the stronger far-domain stretch begins beyond 50 RE.  Both
 * joins are C1 continuous.  The v3 logical knots are 3/8 and 3/4 so that the
 * 64/128/256-cell production families place faces exactly at 15 and 50 RE;
 * the 64-cell floor has 24 uniform 0.5 RE cells inside 15 RE.  Version 4 is
 * the 2.5--200 RE candidate production profile: 2.5--14 RE is uniform, and a
 * single C1 exponential segment spans 14--200 RE.  Its logical join is 23/72,
 * giving exactly 23 uniform 0.5 RE cells on the 72-cell quick-test grid.
 */
int gamera_no_radial_map_init(gamera_no_radial_map *map, int version,
                              double inner_radius, double outer_radius,
                              double legacy_stretch);

/* Continue into halos; v2 uses a positive C1 exponential inside the wall. */
int gamera_no_radial_map_forward(const gamera_no_radial_map *map,
                                 double logical_fraction, double *radius);
int gamera_no_radial_map_inverse(const gamera_no_radial_map *map,
                                 double radius, double *logical_fraction);

/* Version-specific checkpoint identity, excluding separately stored bounds. */
int gamera_no_radial_map_parameters(
    const gamera_no_radial_map *map,
    double parameters[GAMERA_NO_RADIAL_MAP_PARAMETER_COUNT]);

#endif
