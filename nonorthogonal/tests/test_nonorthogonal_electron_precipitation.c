#include "nonorthogonal_electron_precipitation.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define REQUIRE(condition)                                                     \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "requirement failed at %s:%d: %s\n", __FILE__,         \
              __LINE__, #condition);                                           \
      return EXIT_FAILURE;                                                     \
    }                                                                          \
  } while (0)

static int close_relative(double value, double expected, double tolerance) {
  return isfinite(value) && isfinite(expected) &&
         fabs(value - expected) <=
             tolerance * fmax(1.0, fmax(fabs(value), fabs(expected)));
}

static int test_robinson_reference(void) {
  const double pedersen = gamera_mi_robinson_pedersen(5.0, 1.0);
  const double hall = gamera_mi_robinson_hall(
      5.0, pedersen, GAMERA_MI_HALL_ROBINSON_1987);
  REQUIRE(close_relative(pedersen, 4.878048780487805, 2.0e-15));
  REQUIRE(close_relative(hall, 8.621506429373024, 2.0e-15));
  REQUIRE(isnan(gamera_mi_robinson_pedersen(-1.0, 1.0)));
  REQUIRE(isnan(gamera_mi_robinson_hall(
      1.0, -1.0, GAMERA_MI_HALL_ROBINSON_1987)));
  return EXIT_SUCCESS;
}

static int test_fedder_reference(void) {
  gamera_mi_fedder95_config config = gamera_mi_fedder95_default_config();
  gamera_mi_electron_precipitation north;
  gamera_mi_electron_precipitation south;
  const double density = 1.16 * 1.67262192369e-21;
  REQUIRE(gamera_mi_fedder95_electron_precipitation(
              &config, density, 4.0e5, -1.0e-7, 2.0, 1.0, -1,
              &north) == 0);
  REQUIRE(gamera_mi_fedder95_electron_precipitation(
              &config, density, 4.0e5, 1.0e-7, 2.0, 1.0, 1,
              &south) == 0);
  REQUIRE(close_relative(north.source_average_energy_kev,
                         2.0020245096450697, 2.0e-14));
  REQUIRE(north.is_diffuse == 0);
  REQUIRE(close_relative(north.field_aligned_energy_change_kev,
                         3.0818962787574886, 2.0e-14));
  REQUIRE(close_relative(north.average_energy_kev, 5.083920788402558,
                         2.0e-14));
  REQUIRE(close_relative(north.number_flux_cm2_s, 42247350.33567401,
                         2.0e-14));
  REQUIRE(close_relative(north.energy_flux_erg_cm2_s,
                         0.3441189944036332, 2.0e-14));
  REQUIRE(close_relative(north.pedersen_siemens, 2.850732964986577,
                         2.0e-14));
  REQUIRE(close_relative(north.hall_siemens, 5.110201328884505,
                         2.0e-14));
  REQUIRE(close_relative(south.number_flux_cm2_s,
                         north.number_flux_cm2_s, 2.0e-14));
  REQUIRE(close_relative(south.energy_flux_erg_cm2_s,
                         north.energy_flux_erg_cm2_s, 2.0e-14));
  REQUIRE(close_relative(south.pedersen_siemens,
                         north.pedersen_siemens, 2.0e-14));
  REQUIRE(close_relative(south.hall_siemens,
                         north.hall_siemens, 2.0e-14));
  return EXIT_SUCCESS;
}

static int test_beta_scaling(void) {
  gamera_mi_fedder95_config baseline_config =
      gamera_mi_fedder95_default_config();
  gamera_mi_fedder95_config doubled_config = baseline_config;
  doubled_config.beta *= 2.0;
  gamera_mi_electron_precipitation baseline;
  gamera_mi_electron_precipitation doubled;
  const double density = 1.16 * 1.67262192369e-21;
  REQUIRE(gamera_mi_fedder95_electron_precipitation(
              &baseline_config, density, 4.0e5, -1.0e-7, 2.0, 1.0, -1,
              &baseline) == 0);
  REQUIRE(gamera_mi_fedder95_electron_precipitation(
              &doubled_config, density, 4.0e5, -1.0e-7, 2.0, 1.0, -1,
              &doubled) == 0);
  REQUIRE(close_relative(doubled.source_average_energy_kev,
                         baseline.source_average_energy_kev, 2.0e-14));
  REQUIRE(close_relative(doubled.field_aligned_energy_change_kev,
                         baseline.field_aligned_energy_change_kev,
                         2.0e-14));
  REQUIRE(close_relative(doubled.average_energy_kev,
                         baseline.average_energy_kev, 2.0e-14));
  REQUIRE(close_relative(doubled.number_flux_cm2_s,
                         2.0 * baseline.number_flux_cm2_s, 2.0e-14));
  REQUIRE(close_relative(doubled.energy_flux_erg_cm2_s,
                         2.0 * baseline.energy_flux_erg_cm2_s, 2.0e-14));
  REQUIRE(close_relative(doubled.pedersen_siemens,
                         sqrt(2.0) * baseline.pedersen_siemens, 2.0e-14));
  REQUIRE(close_relative(doubled.hall_siemens,
                         sqrt(2.0) * baseline.hall_siemens, 2.0e-14));
  return EXIT_SUCCESS;
}

static int test_drop_cap_and_hall_options(void) {
  gamera_mi_fedder95_config config = gamera_mi_fedder95_default_config();
  gamera_mi_electron_precipitation result;
  const double density = 1.16 * 1.67262192369e-22;
  REQUIRE(gamera_mi_fedder95_electron_precipitation(
              &config, density, 4.0e5, -1.0e-5, 2.0, 1.0, -1,
              &result) == 0);
  REQUIRE(result.field_aligned_energy_change_kev > 20.0);
  REQUIRE(result.is_diffuse == 0);
  REQUIRE(close_relative(result.average_energy_kev,
                         result.source_average_energy_kev + 20.0,
                         2.0e-14));

  const double legacy = gamera_mi_robinson_hall(
      30.0, 5.0, GAMERA_MI_HALL_ROBINSON_1987);
  const double capped = gamera_mi_robinson_hall(
      30.0, 5.0, GAMERA_MI_HALL_ROBINSON_ENERGY_CAP);
  const double kaeppler = gamera_mi_robinson_hall(
      30.0, 5.0, GAMERA_MI_HALL_KAEPPLER_2015);
  REQUIRE(legacy > capped);
  REQUIRE(capped > 0.0);
  REQUIRE(kaeppler > 0.0);
  REQUIRE(gamera_mi_fedder95_electron_precipitation(
              &config, density, 4.0e5, 0.0, 2.0, 1.0, -1,
              &result) == 0);
  REQUIRE(result.is_diffuse == 1);
  return EXIT_SUCCESS;
}

static int test_conductance_combination(void) {
  double pedersen;
  double hall;
  REQUIRE(gamera_mi_combine_conductance(2.0, 1.0, 3.0, 20.0, 2.0,
                                        1.0, 4.0, &pedersen, &hall) == 0);
  REQUIRE(close_relative(pedersen, sqrt(13.0), 2.0e-15));
  REQUIRE(close_relative(hall, 4.0 * sqrt(13.0), 2.0e-15));
  REQUIRE(gamera_mi_combine_conductance(-1.0, 1.0, 3.0, 4.0, 2.0,
                                        1.0, 4.0, &pedersen, &hall) != 0);
  return EXIT_SUCCESS;
}

static int test_lompe_euv(void) {
  double pedersen;
  double hall;
  double pedersen_89;
  double hall_89;
  double pedersen_90;
  double hall_90;
  double pedersen_100;
  double hall_100;
  double pedersen_119;
  double hall_119;
  const double f107 = 150.0;
  REQUIRE(gamera_mi_lompe_euv_conductance(
              0.0, f107, &pedersen, &hall) == 0);
  REQUIRE(close_relative(pedersen, pow(f107, 0.49) * 1.27, 2.0e-15));
  REQUIRE(close_relative(hall, pow(f107, 0.53) * 1.35, 2.0e-15));
  REQUIRE(gamera_mi_lompe_euv_conductance(
              89.0 * acos(-1.0) / 180.0, f107,
              &pedersen_89, &hall_89) == 0);
  REQUIRE(gamera_mi_lompe_euv_conductance(
              90.0 * acos(-1.0) / 180.0, f107,
              &pedersen_90, &hall_90) == 0);
  REQUIRE(gamera_mi_lompe_euv_conductance(
              100.0 * acos(-1.0) / 180.0, f107,
              &pedersen_100, &hall_100) == 0);
  REQUIRE(gamera_mi_lompe_euv_conductance(
              119.0 * acos(-1.0) / 180.0, f107,
              &pedersen_119, &hall_119) == 0);
  REQUIRE(pedersen_89 > pedersen_90);
  REQUIRE(pedersen_90 > pedersen_100);
  REQUIRE(pedersen_100 > pedersen_119);
  REQUIRE(pedersen_119 > 0.0);
  REQUIRE(hall_89 > hall_90);
  REQUIRE(hall_90 > hall_100);
  REQUIRE(hall_100 > hall_119);
  REQUIRE(hall_119 > 0.0);
  REQUIRE(gamera_mi_lompe_euv_conductance(
              120.0 * acos(-1.0) / 180.0, f107,
              &pedersen, &hall) == 0);
  REQUIRE(close_relative(pedersen, 0.0, 2.0e-15));
  REQUIRE(close_relative(hall, 0.0, 2.0e-15));
  REQUIRE(gamera_mi_lompe_euv_conductance(
              2.0 * acos(-1.0) / 3.0, f107,
              &pedersen, &hall) == 0);
  REQUIRE(close_relative(pedersen, 0.0, 2.0e-15));
  REQUIRE(close_relative(hall, 0.0, 2.0e-15));
  REQUIRE(gamera_mi_lompe_euv_conductance(
              0.0, 0.0, &pedersen, &hall) != 0);
  return EXIT_SUCCESS;
}

static int test_mi_solar_zenith_mapping(void) {
  const double pi = acos(-1.0);
  const double colatitude = pi / 6.0;
  const double subsolar_colatitude = pi / 2.0;
  const double subsolar_longitude = pi;
  /* Stored phi=pi is 12 MLT; its zenith angle is the magnetic latitude. */
  REQUIRE(close_relative(
      gamera_mi_solar_zenith_angle(colatitude, pi, subsolar_colatitude,
                                   subsolar_longitude),
      pi / 3.0, 2.0e-15));
  /* Stored phi=0 is 00 MLT and must be on the nightside. */
  REQUIRE(close_relative(
      gamera_mi_solar_zenith_angle(colatitude, 0.0, subsolar_colatitude,
                                   subsolar_longitude),
      2.0 * pi / 3.0, 2.0e-15));
  REQUIRE(close_relative(
      gamera_mi_solar_zenith_angle(colatitude, pi / 2.0,
                                   subsolar_colatitude,
                                   subsolar_longitude),
      pi / 2.0, 2.0e-15));
  REQUIRE(isnan(gamera_mi_solar_zenith_angle(
      -1.0, 0.0, subsolar_colatitude, subsolar_longitude)));
  return EXIT_SUCCESS;
}

int main(void) {
  if (test_robinson_reference() != EXIT_SUCCESS ||
      test_fedder_reference() != EXIT_SUCCESS ||
      test_beta_scaling() != EXIT_SUCCESS ||
      test_drop_cap_and_hall_options() != EXIT_SUCCESS ||
      test_conductance_combination() != EXIT_SUCCESS ||
      test_lompe_euv() != EXIT_SUCCESS ||
      test_mi_solar_zenith_mapping() != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }
  puts("nonorthogonal electron precipitation/conductance tests passed");
  return EXIT_SUCCESS;
}
