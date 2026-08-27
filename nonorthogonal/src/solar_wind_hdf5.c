#include "solar_wind.h"

#include "hdf5.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int dataset_length(hid_t file, const char *name, size_t *length) {
  if (H5Lexists(file, name, H5P_DEFAULT) <= 0) {
    return 1;
  }
  hid_t dataset = H5Dopen2(file, name, H5P_DEFAULT);
  if (dataset < 0) {
    return -1;
  }
  hid_t space = H5Dget_space(dataset);
  if (space < 0) {
    H5Dclose(dataset);
    return -1;
  }
  const int rank = H5Sget_simple_extent_ndims(space);
  hsize_t extent = 0;
  const int status = rank == 1 ? H5Sget_simple_extent_dims(space, &extent, NULL)
                               : -1;
  H5Sclose(space);
  H5Dclose(dataset);
  if (status < 0 || extent == 0 || extent > (hsize_t)SIZE_MAX) {
    return -1;
  }
  *length = (size_t)extent;
  return 0;
}

static int read_array(hid_t file, const char *name, size_t expected,
                      double *values) {
  size_t length = 0U;
  if (dataset_length(file, name, &length) != 0 || length != expected) {
    return -1;
  }
  hid_t dataset = H5Dopen2(file, name, H5P_DEFAULT);
  if (dataset < 0) {
    return -1;
  }
  const herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                                H5P_DEFAULT, values);
  H5Dclose(dataset);
  return status >= 0 ? 0 : -1;
}

static int read_optional_scalar(hid_t file, const char *name, double *value) {
  if (H5Lexists(file, name, H5P_DEFAULT) <= 0) {
    return 1;
  }
  hid_t dataset = H5Dopen2(file, name, H5P_DEFAULT);
  if (dataset < 0) {
    return -1;
  }
  hid_t space = H5Dget_space(dataset);
  if (space < 0) {
    H5Dclose(dataset);
    return -1;
  }
  const int rank = H5Sget_simple_extent_ndims(space);
  hsize_t extent = 1U;
  if (rank == 1 && H5Sget_simple_extent_dims(space, &extent, NULL) < 0) {
    extent = 0U;
  }
  const int valid = rank == 0 || (rank == 1 && extent == 1U);
  const herr_t status = valid ? H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
                                        H5S_ALL, H5P_DEFAULT, value)
                              : -1;
  H5Sclose(space);
  H5Dclose(dataset);
  return status >= 0 && isfinite(*value) ? 0 : -1;
}

static double *allocate_array(size_t count) {
  return (double *)malloc(count * sizeof(double));
}

int gamera_solar_wind_load_hdf5(
    const char *filename, const gamera_solar_wind_hdf5_units *units,
    gamera_solar_wind_series *series) {
  if (filename == NULL || units == NULL || series == NULL ||
      !(units->time_norm > 0.0) || !(units->density_norm > 0.0) ||
      !(units->velocity_norm > 0.0) || !(units->pressure_norm > 0.0) ||
      !(units->magnetic_norm > 0.0) ||
      !(units->velocity_si_scale > 0.0)) {
    return -1;
  }
  hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) {
    return -1;
  }
  size_t count = 0U;
  if (dataset_length(file, "/T", &count) != 0) {
    H5Fclose(file);
    return -1;
  }
  double *time = allocate_array(count);
  double *density = allocate_array(count);
  double *pressure = allocate_array(count);
  double *temperature = allocate_array(count);
  double *velocity[3] = {allocate_array(count), allocate_array(count),
                         allocate_array(count)};
  double *magnetic[3] = {allocate_array(count), allocate_array(count),
                         allocate_array(count)};
  int failed = time == NULL || density == NULL || pressure == NULL ||
               temperature == NULL || velocity[0] == NULL ||
               velocity[1] == NULL || velocity[2] == NULL ||
               magnetic[0] == NULL || magnetic[1] == NULL ||
               magnetic[2] == NULL;
  const char *velocity_name[3] = {"/Vx", "/Vy", "/Vz"};
  const char *magnetic_name[3] = {"/Bx", "/By", "/Bz"};
  if (!failed &&
      (read_array(file, "/T", count, time) != 0 ||
       read_array(file, "/D", count, density) != 0)) {
    failed = 1;
  }
  for (int component = 0; component < 3 && !failed; ++component) {
    if (read_array(file, velocity_name[component], count,
                   velocity[component]) != 0 ||
        read_array(file, magnetic_name[component], count,
                   magnetic[component]) != 0) {
      failed = 1;
    }
  }
  int have_pressure = 0;
  if (!failed && H5Lexists(file, "/P", H5P_DEFAULT) > 0) {
    have_pressure = read_array(file, "/P", count, pressure) == 0;
    failed = !have_pressure;
  } else if (!failed && units->physical_units &&
             H5Lexists(file, "/Temp", H5P_DEFAULT) > 0) {
    failed = read_array(file, "/Temp", count, temperature) != 0;
  } else if (!failed) {
    failed = 1;
  }

  gamera_solar_wind_state *state = NULL;
  if (!failed) {
    state = (gamera_solar_wind_state *)malloc(count * sizeof(*state));
    failed = state == NULL;
  }
  const double proton_mass = 1.67262192369e-27;
  const double boltzmann = 1.380649e-23;
  for (size_t sample = 0; sample < count && !failed; ++sample) {
    if (units->physical_units) {
      const double number_density_si = density[sample] * 1.0e6;
      time[sample] /= units->time_norm;
      state[sample].density =
          number_density_si * proton_mass / units->density_norm;
      state[sample].pressure =
          (have_pressure ? pressure[sample] * 1.0e-9
                         : number_density_si * boltzmann * temperature[sample]) /
          units->pressure_norm;
      for (int component = 0; component < 3; ++component) {
        state[sample].velocity[component] =
            velocity[component][sample] * units->velocity_si_scale /
            units->velocity_norm;
        state[sample].magnetic[component] =
            magnetic[component][sample] * 1.0e-9 / units->magnetic_norm;
      }
    } else {
      state[sample].density = density[sample];
      state[sample].pressure = pressure[sample];
      for (int component = 0; component < 3; ++component) {
        state[sample].velocity[component] = velocity[component][sample];
        state[sample].magnetic[component] = magnetic[component][sample];
      }
    }
  }

  double by_coefficient = 0.0;
  double bz_coefficient = 0.0;
  double bx_offset = 0.0;
  const int have_by =
      !failed ? read_optional_scalar(file, "/ByC", &by_coefficient) : -1;
  const int have_bz =
      !failed ? read_optional_scalar(file, "/BzC", &bz_coefficient) : -1;
  const int have_bx0 =
      !failed ? read_optional_scalar(file, "/Bx0", &bx_offset) : -1;
  if (have_by < 0 || have_bz < 0 || have_bx0 < 0 ||
      ((have_by == 0) != (have_bz == 0))) {
    failed = 1;
  }
  if (!failed && units->physical_units && have_bx0 == 0) {
    bx_offset = bx_offset * 1.0e-9 / units->magnetic_norm;
  }
  H5Fclose(file);

  int result = -1;
  if (!failed) {
    series->by_coefficient = by_coefficient;
    series->bz_coefficient = bz_coefficient;
    series->bx_offset = bx_offset;
    /* Preserve the measured Bx until the caller chooses the Lyon relation. */
    series->enforce_bx_relation = 0;
    result = gamera_solar_wind_set(series, count, time, state);
    series->enforce_bx_relation = have_by == 0 && have_bz == 0;
  }
  free(time);
  free(density);
  free(pressure);
  free(temperature);
  free(state);
  for (int component = 0; component < 3; ++component) {
    free(velocity[component]);
    free(magnetic[component]);
  }
  return result;
}
