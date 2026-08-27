#ifndef GAMERA_SOLAR_WIND_H
#define GAMERA_SOLAR_WIND_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

enum {
  GAMERA_WIND_X = 0,
  GAMERA_WIND_Y = 1,
  GAMERA_WIND_Z = 2,
  GAMERA_WIND_DIM = 3
};

typedef struct {
  double density;
  double pressure;
  double velocity[GAMERA_WIND_DIM];
  double magnetic[GAMERA_WIND_DIM];
} gamera_solar_wind_state;

typedef struct {
  size_t count;
  double *time;
  gamera_solar_wind_state *state;
  double reference[GAMERA_WIND_DIM];
  double time_offset;
  double by_coefficient;
  double bz_coefficient;
  double bx_offset;
  int enforce_bx_relation;
  int linear_interpolation;
} gamera_solar_wind_series;

typedef struct {
  int physical_units;
  double time_norm;
  double density_norm;
  double velocity_norm;
  double pressure_norm;
  double magnetic_norm;
  /* Multiplier taking the input velocity to m/s (1 for Kaiju files). */
  double velocity_si_scale;
} gamera_solar_wind_hdf5_units;

void gamera_solar_wind_init(gamera_solar_wind_series *series);
void gamera_solar_wind_destroy(gamera_solar_wind_series *series);

/* Copy a validated, monotonically increasing time series. */
int gamera_solar_wind_set(gamera_solar_wind_series *series, size_t count,
                          const double *time,
                          const gamera_solar_wind_state *state);
void gamera_solar_wind_apply_bx_relation(gamera_solar_wind_series *series);

/* Read /T,/D,/Vx,/Vy,/Vz,/Bx,/By,/Bz and /P or /Temp. */
int gamera_solar_wind_load_hdf5(
    const char *filename, const gamera_solar_wind_hdf5_units *units,
    gamera_solar_wind_series *series);

/* Sample only the monitor time series, with endpoint clamping. */
int gamera_solar_wind_sample_time(const gamera_solar_wind_series *series,
                                  double time,
                                  gamera_solar_wind_state *state);

/*
 * Lyon et al. (2004), section 5.1 / Kaiju wind.F90 ballistic mapping.
 * The current monitor Vx defines a single-valued front velocity, then the
 * monitor is sampled at time - dot(x-x_ref,v_front)/|v_front|^2.
 */
int gamera_solar_wind_sample_at(const gamera_solar_wind_series *series,
                                const double point[GAMERA_WIND_DIM],
                                double simulation_time,
                                gamera_solar_wind_state *state,
                                double *delay);

/* Kaiju's smooth upstream-hemisphere-to-tail boundary weight. */
int gamera_solar_wind_weight(const gamera_solar_wind_series *series,
                             const double point[GAMERA_WIND_DIM],
                             double simulation_time, double *weight);

#ifdef __cplusplus
}
#endif

#endif
