#ifndef GAMERA_NONORTHOGONAL_ELECTRON_PRECIPITATION_H
#define GAMERA_NONORTHOGONAL_ELECTRON_PRECIPITATION_H

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  GAMERA_MI_HALL_ROBINSON_1987 = 0,
  GAMERA_MI_HALL_KAIJU_FEDDER_CAP = 1,
  GAMERA_MI_HALL_KAEPPLER_2015 = 2
} gamera_mi_hall_model;

typedef struct {
  double alpha;
  double beta;
  double current_voltage_scale;
  double outward_resistance;
  double inward_resistance;
  double maximum_acceleration_kev;
  double minimum_average_energy_kev;
  double helium_mass_factor;
  /* Kaiju limits only the density used by the current-voltage relation to
   * coefficient * background Pedersen conductance.  A non-positive value
   * disables that optional compatibility limit. */
  double fac_density_cap_kg_m3_per_siemens;
  gamera_mi_hall_model hall_model;
} gamera_mi_fedder95_config;

typedef struct {
  int is_diffuse;
  double source_average_energy_kev;
  double field_aligned_energy_change_kev;
  double average_energy_kev;
  double number_flux_cm2_s;
  double energy_flux_erg_cm2_s;
  double pedersen_siemens;
  double hall_siemens;
} gamera_mi_electron_precipitation;

/* Kaiju/REMIX Fedder defaults, with the published Robinson-1987 Hall fit. */
gamera_mi_fedder95_config gamera_mi_fedder95_default_config(void);

/*
 * Evaluate the Fedder-1995 causal electron precipitation model.
 *
 * mass_density_kg_m3 and sound_speed_m_s are sampled from the MHD mapping
 * shell. fac_a_m2 uses the Cartesian sign stored by the MHD code; hemisphere
 * is -1 for North and +1 for South, matching the existing M-I mapping.
 * background_pedersen_siemens is used only by the optional Kaiju density cap.
 * ramp_factor is dimensionless and normally one after spin-up.
 */
int gamera_mi_fedder95_electron_precipitation(
    const gamera_mi_fedder95_config *config, double mass_density_kg_m3,
    double sound_speed_m_s, double fac_a_m2,
    double background_pedersen_siemens, double ramp_factor, int hemisphere,
    gamera_mi_electron_precipitation *result);

double gamera_mi_robinson_pedersen(double average_energy_kev,
                                   double energy_flux_erg_cm2_s);
double gamera_mi_robinson_hall(double average_energy_kev,
                               double pedersen_siemens,
                               gamera_mi_hall_model model);

/* Spherical dot-product solar zenith angle.  All angles are radians in the
 * model's native M-I coordinate convention.  The caller obtains the
 * subsolar angles from the native Cartesian sunward vector. */
double gamera_mi_solar_zenith_angle(double colatitude_rad,
                                    double longitude_rad,
                                    double subsolar_colatitude_rad,
                                    double subsolar_longitude_rad);

/* Quadrature addition used by REMIX for background plus auroral conductance. */
int gamera_mi_combine_conductance(double background_pedersen_siemens,
                                  double background_hall_siemens,
                                  double auroral_pedersen_siemens,
                                  double auroral_hall_siemens,
                                  double minimum_pedersen_siemens,
                                  double minimum_hall_siemens,
                                  double maximum_hall_to_pedersen,
                                  double *pedersen_siemens,
                                  double *hall_siemens);

/* Kaiju/REMIX LOMPE solar-EUV conductance.  The tabulated q' fit and its
 * twilight tail are used through 120 degrees solar zenith angle, exactly as
 * in Kaiju's InterpQP routine; the EUV term is zero at and beyond 120 degrees.
 * F10.7 is in sfu. */
int gamera_mi_lompe_euv_conductance(double solar_zenith_rad, double f107,
                                    double *pedersen_siemens,
                                    double *hall_siemens);

#ifdef __cplusplus
}
#endif

#endif
