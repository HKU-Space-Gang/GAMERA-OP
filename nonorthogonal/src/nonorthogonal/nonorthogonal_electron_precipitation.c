#include "nonorthogonal_electron_precipitation.h"

#include <float.h>
#include <math.h>
#include <stddef.h>

static const double PROTON_MASS_G = 1.67262192369e-24;
static const double ELEMENTARY_CHARGE_C = 1.60217663e-19;
static const double KEV_TO_ERG = 1.602176634e-9;

/* Published LOMPE q'(solar zenith angle) fit at one-degree spacing. */
static const double LOMPE_QP[121] = {
    1.000000e+00, 9.998500e-01, 9.994000e-01, 9.986500e-01,
    9.976000e-01, 9.962510e-01, 9.946030e-01, 9.926570e-01,
    9.904130e-01, 9.878720e-01, 9.850340e-01, 9.819010e-01,
    9.784740e-01, 9.747540e-01, 9.707420e-01, 9.664380e-01,
    9.618460e-01, 9.569650e-01, 9.517980e-01, 9.463470e-01,
    9.406120e-01, 9.345960e-01, 9.283010e-01, 9.217290e-01,
    9.148810e-01, 9.077610e-01, 9.003700e-01, 8.927110e-01,
    8.847860e-01, 8.765980e-01, 8.681500e-01, 8.594440e-01,
    8.504830e-01, 8.412700e-01, 8.318080e-01, 8.221000e-01,
    8.121500e-01, 8.019600e-01, 7.915350e-01, 7.808770e-01,
    7.699900e-01, 7.588780e-01, 7.475450e-01, 7.359940e-01,
    7.242290e-01, 7.122560e-01, 7.000760e-01, 6.876960e-01,
    6.751190e-01, 6.623500e-01, 6.493930e-01, 6.362540e-01,
    6.229360e-01, 6.094460e-01, 5.957870e-01, 5.819660e-01,
    5.679880e-01, 5.538590e-01, 5.395830e-01, 5.251680e-01,
    5.106190e-01, 4.959420e-01, 4.811460e-01, 4.662350e-01,
    4.512190e-01, 4.361030e-01, 4.208970e-01, 4.056090e-01,
    3.902480e-01, 3.748240e-01, 3.593460e-01, 3.438260e-01,
    3.282750e-01, 3.127080e-01, 2.971370e-01, 2.815790e-01,
    2.660500e-01, 2.505700e-01, 2.351610e-01, 2.198450e-01,
    2.046500e-01, 1.896070e-01, 1.747500e-01, 1.601180e-01,
    1.457550e-01, 1.317130e-01, 1.180480e-01, 1.048230e-01,
    9.211030e-02, 7.998570e-02, 6.853160e-02, 5.783310e-02,
    4.797420e-02, 3.903300e-02, 3.107490e-02, 2.414490e-02,
    1.826040e-02, 1.340510e-02, 9.526090e-03, 6.535640e-03,
    4.318240e-03, 2.741420e-03, 1.668790e-03, 9.722760e-04,
    5.412880e-04, 2.875250e-04, 1.455210e-04, 7.008060e-05,
    3.207010e-05, 1.392600e-05, 5.729740e-06, 2.230220e-06,
    8.198480e-07, 2.841240e-07, 9.264680e-08, 2.836570e-08,
    8.136250e-09, 2.181100e-09, 5.450290e-10, 1.266040e-10,
    2.725590e-11};

static int finite_nonnegative(double value) {
  return isfinite(value) && value >= 0.0;
}

double gamera_mi_solar_zenith_angle(double colatitude_rad,
                                    double longitude_rad,
                                    double subsolar_colatitude_rad,
                                    double subsolar_longitude_rad) {
  const double pi = acos(-1.0);
  if (!isfinite(colatitude_rad) || colatitude_rad < 0.0 ||
      colatitude_rad > pi || !isfinite(longitude_rad) ||
      !isfinite(subsolar_colatitude_rad) ||
      subsolar_colatitude_rad < 0.0 || subsolar_colatitude_rad > pi ||
      !isfinite(subsolar_longitude_rad)) {
    return NAN;
  }
  const double cosine_zenith =
      cos(colatitude_rad) * cos(subsolar_colatitude_rad) +
      sin(colatitude_rad) * sin(subsolar_colatitude_rad) *
          cos(longitude_rad - subsolar_longitude_rad);
  return acos(fmax(-1.0, fmin(1.0, cosine_zenith)));
}

int gamera_mi_lompe_euv_conductance(double solar_zenith_rad, double f107,
                                    double *pedersen_siemens,
                                    double *hall_siemens) {
  if (!isfinite(solar_zenith_rad) || solar_zenith_rad < 0.0 ||
      !isfinite(f107) || !(f107 > 0.0) || pedersen_siemens == NULL ||
      hall_siemens == NULL) {
    return -1;
  }
  const double zenith_deg = 180.0 * solar_zenith_rad / acos(-1.0);
  double qprime;
  if (zenith_deg <= 0.0) {
    qprime = 1.0;
  } else if (zenith_deg >=
             120.0 - 64.0 * DBL_EPSILON * fmax(1.0, zenith_deg)) {
    /* The tabulated twilight fit terminates at 120 degrees. */
    qprime = 0.0;
  } else {
    const size_t lower = (size_t)floor(zenith_deg);
    const double fraction = zenith_deg - (double)lower;
    qprime = (1.0 - fraction) * LOMPE_QP[lower] +
             fraction * LOMPE_QP[lower + 1U];
  }
  qprime = fmax(0.0, qprime);
  *pedersen_siemens =
      pow(f107, 0.49) * (0.34 * qprime + 0.93 * sqrt(qprime));
  *hall_siemens =
      pow(f107, 0.53) * (0.81 * qprime + 0.54 * sqrt(qprime));
  return isfinite(*pedersen_siemens) && isfinite(*hall_siemens) ? 0 : -1;
}

gamera_mi_fedder95_config gamera_mi_fedder95_default_config(void) {
  const gamera_mi_fedder95_config config = {
      .alpha = 1.0332467,
      .beta = 0.4362323,
      .current_voltage_scale = 0.083567956,
      .outward_resistance = 6.0,
      .inward_resistance = 1.2,
      .maximum_acceleration_kev = 20.0,
      .minimum_average_energy_kev = 1.0e-12,
      .helium_mass_factor = 1.16,
      .fac_density_cap_kg_m3_per_siemens = 1.65e-21,
      .hall_model = GAMERA_MI_HALL_ROBINSON_1987};
  return config;
}

static int valid_config(const gamera_mi_fedder95_config *config) {
  return config != NULL && isfinite(config->alpha) && config->alpha > 0.0 &&
         isfinite(config->beta) && config->beta >= 0.0 &&
         isfinite(config->current_voltage_scale) &&
         config->current_voltage_scale >= 0.0 &&
         isfinite(config->outward_resistance) &&
         config->outward_resistance > 0.0 &&
         isfinite(config->inward_resistance) &&
         config->inward_resistance > 0.0 &&
         isfinite(config->maximum_acceleration_kev) &&
         config->maximum_acceleration_kev >= 0.0 &&
         isfinite(config->minimum_average_energy_kev) &&
         config->minimum_average_energy_kev > 0.0 &&
         isfinite(config->helium_mass_factor) &&
         config->helium_mass_factor > 0.0 &&
         isfinite(config->fac_density_cap_kg_m3_per_siemens) &&
         config->fac_density_cap_kg_m3_per_siemens >= 0.0 &&
         config->hall_model >= GAMERA_MI_HALL_ROBINSON_1987 &&
         config->hall_model <= GAMERA_MI_HALL_KAEPPLER_2015;
}

double gamera_mi_robinson_pedersen(double average_energy_kev,
                                   double energy_flux_erg_cm2_s) {
  if (!finite_nonnegative(average_energy_kev) ||
      !finite_nonnegative(energy_flux_erg_cm2_s)) {
    return NAN;
  }
  const double denominator = 16.0 + average_energy_kev * average_energy_kev;
  return 40.0 * average_energy_kev * sqrt(energy_flux_erg_cm2_s) /
         denominator;
}

double gamera_mi_robinson_hall(double average_energy_kev,
                               double pedersen_siemens,
                               gamera_mi_hall_model model) {
  if (!finite_nonnegative(average_energy_kev) ||
      !finite_nonnegative(pedersen_siemens)) {
    return NAN;
  }
  if (model == GAMERA_MI_HALL_ROBINSON_1987) {
    return 0.45 * pedersen_siemens * pow(average_energy_kev, 0.85);
  }
  if (model == GAMERA_MI_HALL_ROBINSON_ENERGY_CAP) {
    return 0.45 * pedersen_siemens * pow(average_energy_kev, 0.85) /
           (1.0 + 0.0025 * average_energy_kev * average_energy_kev);
  }
  if (model == GAMERA_MI_HALL_KAEPPLER_2015) {
    return 0.57 * pedersen_siemens * pow(average_energy_kev, 0.53);
  }
  return NAN;
}

int gamera_mi_fedder95_electron_precipitation(
    const gamera_mi_fedder95_config *config, double mass_density_kg_m3,
    double sound_speed_m_s, double fac_a_m2,
    double background_pedersen_siemens, double ramp_factor, int hemisphere,
    gamera_mi_electron_precipitation *result) {
  if (!valid_config(config) || result == NULL ||
      !isfinite(mass_density_kg_m3) || mass_density_kg_m3 <= 0.0 ||
      !isfinite(sound_speed_m_s) || sound_speed_m_s < 0.0 ||
      !isfinite(fac_a_m2) || !finite_nonnegative(background_pedersen_siemens) ||
      !isfinite(ramp_factor) || ramp_factor < 0.0 ||
      (hemisphere != -1 && hemisphere != 1)) {
    return -1;
  }

  const double density_g_cm3 = 1.0e-3 * mass_density_kg_m3;
  const double sound_speed_cm_s = 100.0 * sound_speed_m_s;
  const double effective_mass_g =
      config->helium_mass_factor * PROTON_MASS_G;
  const double erg_to_kev = 1.0 / KEV_TO_ERG;
  const double source_average_energy_kev =
      config->alpha * effective_mass_g * erg_to_kev *
      sound_speed_cm_s * sound_speed_cm_s * ramp_factor;
  if (!finite_nonnegative(source_average_energy_kev)) {
    return -1;
  }

  const double base_number_flux =
      sqrt(KEV_TO_ERG) / pow(effective_mass_g, 1.5) * config->beta *
      density_g_cm3 * sqrt(source_average_energy_kev) * ramp_factor;
  if (!finite_nonnegative(base_number_flux)) {
    return -1;
  }

  double current_voltage_density_g_cm3 = density_g_cm3;
  if (config->fac_density_cap_kg_m3_per_siemens > 0.0 &&
      background_pedersen_siemens > 0.0) {
    const double cap_g_cm3 =
        1.0e-3 * config->fac_density_cap_kg_m3_per_siemens *
        background_pedersen_siemens;
    current_voltage_density_g_cm3 =
        fmin(current_voltage_density_g_cm3, cap_g_cm3);
  }
  if (!(current_voltage_density_g_cm3 > 0.0) ||
      !isfinite(current_voltage_density_g_cm3)) {
    return -1;
  }

  const double signed_upward_fac = (double)hemisphere * fac_a_m2;
  const double resistance =
      signed_upward_fac >= 0.0 ? 2.0 * config->outward_resistance
                               : 2.0 * config->inward_resistance;
  const double field_aligned_energy_change_kev =
      pow(effective_mass_g, 1.5) / ELEMENTARY_CHARGE_C * 1.0e-4 *
      sqrt(erg_to_kev) * config->current_voltage_scale * resistance *
      signed_upward_fac * sqrt(source_average_energy_kev) /
      current_voltage_density_g_cm3;
  if (!isfinite(field_aligned_energy_change_kev)) {
    return -1;
  }
  const double clipped_change_kev =
      fmin(field_aligned_energy_change_kev,
           config->maximum_acceleration_kev);
  const double average_energy_kev =
      fmax(source_average_energy_kev + clipped_change_kev,
           config->minimum_average_energy_kev);

  double number_flux_cm2_s;
  if (clipped_change_kev > 0.0 && source_average_energy_kev > DBL_MIN) {
    number_flux_cm2_s =
        base_number_flux *
        (8.0 - 7.0 * exp(-clipped_change_kev /
                         (7.0 * source_average_energy_kev)));
  } else if (source_average_energy_kev > DBL_MIN) {
    number_flux_cm2_s =
        base_number_flux * exp(clipped_change_kev /
                               source_average_energy_kev);
  } else {
    number_flux_cm2_s = 0.0;
  }
  const double energy_flux_erg_cm2_s =
      average_energy_kev * number_flux_cm2_s * KEV_TO_ERG;
  const double pedersen_siemens = gamera_mi_robinson_pedersen(
      average_energy_kev, energy_flux_erg_cm2_s);
  const double hall_siemens = gamera_mi_robinson_hall(
      average_energy_kev, pedersen_siemens, config->hall_model);
  if (!finite_nonnegative(number_flux_cm2_s) ||
      !finite_nonnegative(energy_flux_erg_cm2_s) ||
      !finite_nonnegative(pedersen_siemens) ||
      !finite_nonnegative(hall_siemens)) {
    return -1;
  }

  result->is_diffuse =
      field_aligned_energy_change_kev < source_average_energy_kev;
  result->source_average_energy_kev = source_average_energy_kev;
  result->field_aligned_energy_change_kev =
      field_aligned_energy_change_kev;
  result->average_energy_kev = average_energy_kev;
  result->number_flux_cm2_s = number_flux_cm2_s;
  result->energy_flux_erg_cm2_s = energy_flux_erg_cm2_s;
  result->pedersen_siemens = pedersen_siemens;
  result->hall_siemens = hall_siemens;
  return 0;
}

int gamera_mi_combine_conductance(double background_pedersen_siemens,
                                  double background_hall_siemens,
                                  double auroral_pedersen_siemens,
                                  double auroral_hall_siemens,
                                  double minimum_pedersen_siemens,
                                  double minimum_hall_siemens,
                                  double maximum_hall_to_pedersen,
                                  double *pedersen_siemens,
                                  double *hall_siemens) {
  if (!finite_nonnegative(background_pedersen_siemens) ||
      !finite_nonnegative(background_hall_siemens) ||
      !finite_nonnegative(auroral_pedersen_siemens) ||
      !finite_nonnegative(auroral_hall_siemens) ||
      !finite_nonnegative(minimum_pedersen_siemens) ||
      !finite_nonnegative(minimum_hall_siemens) ||
      !isfinite(maximum_hall_to_pedersen) ||
      maximum_hall_to_pedersen <= 0.0 || pedersen_siemens == NULL ||
      hall_siemens == NULL) {
    return -1;
  }
  *pedersen_siemens =
      fmax(minimum_pedersen_siemens,
           hypot(background_pedersen_siemens, auroral_pedersen_siemens));
  *hall_siemens =
      fmax(minimum_hall_siemens,
           hypot(background_hall_siemens, auroral_hall_siemens));
  *hall_siemens =
      fmin(*hall_siemens, maximum_hall_to_pedersen * *pedersen_siemens);
  return isfinite(*pedersen_siemens) && isfinite(*hall_siemens) ? 0 : -1;
}
