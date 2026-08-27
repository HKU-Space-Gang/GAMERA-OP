#if defined(PROBLEM_NONORTHOGONAL_YINYANG_MHD_OBSTACLE) || \
    defined(PROBLEM_NONORTHOGONAL_YINYANG_IONOSPHERE) || \
    defined(PROBLEM_NONORTHOGONAL_YINYANG_CRUSTAL) || \
    defined(PROBLEM_NONORTHOGONAL_EARTH_MAGNETOSPHERE) || \
    defined(PROBLEM_NONORTHOGONAL_YINYANG_BOW_SHOCK) || \
    defined(PROBLEM_NONORTHOGONAL_EARTH_MAGNETOSPHERE_STRETCHED)

#include "problem.h"

#include "config.h"
#include "log.h"
#include "nonorthogonal_initialization.h"
#include "nonorthogonal_background.h"
#include "nonorthogonal_driver.h"
#include "nonorthogonal_legacy_adapter.h"
#include "nonorthogonal_radial_map.h"
#include "nonorthogonal_state.h"
#include "nonorthogonal_step.h"
#include "nonorthogonal_yinyang.h"
#ifdef GAMERA_MI_COUPLING
#include "nonorthogonal_mi_coupling.h"
#endif
#include "setup_mpi.h"
#include "solar_wind.h"

#include <math.h>

#ifdef GAMERA_BOW_SHOCK
#ifndef GAMERA_BOW_SHOCK_INNER_RADIUS_RE
#ifdef PROBLEM_NONORTHOGONAL_YINYANG_BOW_SHOCK
#define GAMERA_BOW_SHOCK_INNER_RADIUS_RE 4.0
#else
#define GAMERA_BOW_SHOCK_INNER_RADIUS_RE 3.0
#endif
#endif
#ifndef GAMERA_BOW_SHOCK_OUTER_RADIUS_RE
#ifdef PROBLEM_NONORTHOGONAL_YINYANG_BOW_SHOCK
#define GAMERA_BOW_SHOCK_OUTER_RADIUS_RE 100.0
#else
#define GAMERA_BOW_SHOCK_OUTER_RADIUS_RE 200.0
#endif
#endif
#ifndef GAMERA_BOW_SHOCK_RADIAL_STRETCH
#define GAMERA_BOW_SHOCK_RADIAL_STRETCH 3.0
#endif
#ifndef GAMERA_EARTH_RADIAL_MAP_VERSION
#define GAMERA_EARTH_RADIAL_MAP_VERSION 4
#endif
static const double obstacle_radius = GAMERA_BOW_SHOCK_INNER_RADIUS_RE;
static const double outer_radius = GAMERA_BOW_SHOCK_OUTER_RADIUS_RE;
static const double inflow_density = 5.0;
static const double inflow_speed = 4.0;
static const double inflow_mach = 10.0;
static const gamera_no_vec3 upstream_magnetic = {{0.0, 0.0, 0.0}};
static const double radial_stretch = GAMERA_BOW_SHOCK_RADIAL_STRETCH;
static gamera_no_radial_map radial_map;
static int radial_map_ready;
#elif defined(GAMERA_EARTH_DIPOLE_BACKGROUND)
static const double obstacle_radius = 2.5;
static const double outer_radius = 20.0;
static const double inflow_density = 1.0;
static const double inflow_speed = 1.0;
static const double inflow_mach = 2.0;
static const gamera_no_vec3 upstream_magnetic = {{0.0, 0.10, 0.05}};
#else
static const double obstacle_radius = 1.0;
static const double outer_radius = 4.0;
static const double inflow_density = 1.0;
static const double inflow_speed = 1.0;
static const double inflow_mach = 2.0;
static const gamera_no_vec3 upstream_magnetic = {{0.0, 0.10, 0.05}};
#endif

#ifdef GAMERA_BOW_SHOCK
static int requested_radial_map_version(void) {
#ifdef GAMERA_EARTH_MAGNETOSPHERE_STRETCHED
  return GAMERA_EARTH_RADIAL_MAP_VERSION;
#else
  return GAMERA_NO_RADIAL_MAP_LEGACY_EXPONENTIAL;
#endif
}

static int ensure_radial_map(void) {
  if (!radial_map_ready) {
    if (gamera_no_radial_map_init(&radial_map,
                                  requested_radial_map_version(),
                                  obstacle_radius, outer_radius,
                                  radial_stretch) != 0) {
      return -1;
    }
    radial_map_ready = 1;
  }
  return 0;
}
#endif
#ifndef GAMERA_YINYANG_ANGULAR_EXTENSION_DEG
#define GAMERA_YINYANG_ANGULAR_EXTENSION_DEG 2.0
#endif
/* Keep the historical two-degree overlap by default.  The compile-time
 * value is intentionally exposed for overset-fringe convergence tests: a
 * seventh-order reconstruction reaches several cells beyond the ownership
 * switch, so the Earth inner-wall audit also exercises wider overlaps. */
static const double angular_extension =
    GAMERA_YINYANG_ANGULAR_EXTENSION_DEG * PI / 180.0;

static gamera_solar_wind_series time_dependent_wind;
static int time_dependent_wind_ready;
static double initial_density = 1.0;
static double initial_pressure;
static gamera_no_vec3 initial_velocity = {{1.0, 0.0, 0.0}};
static gamera_no_vec3 initial_magnetic = {{0.0, 0.10, 0.05}};

#ifdef GAMERA_EARTH_DIPOLE_BACKGROUND
/* Use a southward code dipole moment corresponding to 0.31 G. */
static const double earth_equatorial_surface_field_tesla = 3.10e-5;
static gamera_no_dipole earth_dipole = {{{0.0, 0.0, -1.0}}};
/* Tenuous-magnetosphere floors in the accepted production code units. The
 * legacy C defaults are intentionally much larger and are inappropriate for
 * this Earth M-I problem. */
static const double earth_density_floor = 1.0e-6;
static const double earth_pressure_floor = 1.0e-8;
#endif

#ifdef GAMERA_EARTH_UPSTREAM_STARTUP
/* Accepted Earth startup parameters with an explicit upstream launch plane. */
static const double earth_cut_radius = 16.0;
static const double earth_cut_length = 8.0;
static const double startup_front_x = -30.0;
static const double startup_front_width = 1.5;
static const double startup_ambient_density = 1.0;
static const double startup_ambient_pressure = 0.001;

static double cubic_ramp_down(double coordinate, double start,
                              double length) {
  const double scaled = (coordinate - start) / length;
  if (scaled <= 0.0) {
    return 1.0;
  }
  if (scaled >= 1.0) {
    return 0.0;
  }
  return 1.0 - 3.0 * scaled * scaled + 2.0 * scaled * scaled * scaled;
}

static double startup_wind_fraction(double x) {
  return 0.5 *
         (1.0 - tanh((x - startup_front_x) / startup_front_width));
}

/*
 * Integral of startup_wind_fraction with a zero downstream gauge.  Its
 * derivative is the wind fraction, so Ay=Bz*F and Az=-By*F make a planar,
 * divergence-free tangential IMF front without filling the inner domain.
 */
static double startup_front_antiderivative(double x) {
  const double q = (x - startup_front_x) / startup_front_width;
  const double absolute_q = fabs(q);
  const double log_cosh =
      absolute_q + log1p(exp(-2.0 * absolute_q)) - log(2.0);
  return 0.5 *
         (x - startup_front_width * log_cosh - startup_front_x -
          startup_front_width * log(2.0));
}
#endif

#ifdef GAMERA_IONOSPHERE_SOURCE_LOSS
static const double ionosphere_scale_height = 0.12;
static const double photoionization_rate = 0.08;
static const double recombination_coefficient = 0.02;
static const double neutral_drag_frequency = 0.15;
static const double injected_specific_internal_energy = 0.075;
#endif

#ifdef GAMERA_CRUSTAL_FIELD
static const gamera_no_vec3 crustal_center = {{0.0, -0.35, -0.25}};
static const gamera_no_vec3 crustal_moment = {{0.0, 0.0, 1.0}};
static const double crustal_amplitude = 0.008;
static const double crustal_softening = 0.15;

static gamera_no_vec3 crustal_vector_potential(gamera_no_vec3 point) {
  const gamera_no_vec3 displacement =
      {{point.value[0] - crustal_center.value[0],
        point.value[1] - crustal_center.value[1],
        point.value[2] - crustal_center.value[2]}};
  const double square = displacement.value[0] * displacement.value[0] +
                        displacement.value[1] * displacement.value[1] +
                        displacement.value[2] * displacement.value[2] +
                        crustal_softening * crustal_softening;
  const double factor = crustal_amplitude / pow(square, 1.5);
  return (gamera_no_vec3){{
      factor * (crustal_moment.value[1] * displacement.value[2] -
                crustal_moment.value[2] * displacement.value[1]),
      factor * (crustal_moment.value[2] * displacement.value[0] -
                crustal_moment.value[0] * displacement.value[2]),
      factor * (crustal_moment.value[0] * displacement.value[1] -
                crustal_moment.value[1] * displacement.value[0])}};
}

static gamera_no_vec3 crustal_magnetic(gamera_no_vec3 point) {
  const gamera_no_vec3 displacement =
      {{point.value[0] - crustal_center.value[0],
        point.value[1] - crustal_center.value[1],
        point.value[2] - crustal_center.value[2]}};
  const double radius_squared =
      displacement.value[0] * displacement.value[0] +
      displacement.value[1] * displacement.value[1] +
      displacement.value[2] * displacement.value[2];
  const double softened_squared =
      radius_squared + crustal_softening * crustal_softening;
  const double inverse_power =
      crustal_amplitude / pow(softened_squared, 2.5);
  double moment_dot_displacement = 0.0;
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    moment_dot_displacement +=
        crustal_moment.value[component] * displacement.value[component];
  }
  gamera_no_vec3 result;
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    result.value[component] =
        inverse_power *
        (3.0 * displacement.value[component] * moment_dot_displacement +
         crustal_moment.value[component] *
             (2.0 * crustal_softening * crustal_softening -
              radius_squared));
  }
  return result;
}
#endif

static double inflow_pressure(void) {
  if (initial_pressure > 0.0) {
    return initial_pressure;
  }
  return inflow_density * inflow_speed * inflow_speed /
         (gamma_val * inflow_mach * inflow_mach);
}

static int inflow_primitive(gamera_no_vec3 point, void *context,
                            gamera_no_primitive *primitive) {
  (void)context;
  if (primitive == NULL) {
    return -1;
  }
#ifdef GAMERA_EARTH_UPSTREAM_STARTUP
  const double wind_fraction = startup_wind_fraction(point.value[0]);
  primitive->density =
      startup_ambient_density +
      wind_fraction * (initial_density - startup_ambient_density);
  primitive->pressure =
      startup_ambient_pressure +
      wind_fraction * (inflow_pressure() - startup_ambient_pressure);
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    primitive->velocity.value[component] =
        wind_fraction * initial_velocity.value[component];
  }
#else
  primitive->density = initial_density;
  primitive->velocity = initial_velocity;
  primitive->pressure = inflow_pressure();
#endif
#if defined(GAMERA_EARTH_DIPOLE_BACKGROUND) || defined(GAMERA_BOW_SHOCK)
  /*
   * Start the hard-wall neighborhood at rest and introduce the upstream flow
   * smoothly.  A uniform 400 km/s flow right against a reflecting sphere
   * creates an avoidable startup rarefaction before the bow shock can form.
   */
  const double radius = sqrt(point.value[0] * point.value[0] +
                             point.value[1] * point.value[1] +
                             point.value[2] * point.value[2]);
#ifdef GAMERA_BOW_SHOCK
  const double startup_width = 2.0;
#else
  const double startup_width = 2.5;
#endif
  double fraction = (radius - obstacle_radius) / startup_width;
  fraction = fmax(0.0, fmin(1.0, fraction));
  const double ramp = fraction * fraction * (3.0 - 2.0 * fraction);
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    primitive->velocity.value[component] *= ramp;
  }
#else
  (void)point;
#endif
  return 0;
}

/*
 * A = 0.5*(1-R^3/r^3) B_inf x r.  This is a uniform IMF plus the potential
 * field of an induced dipole m=-0.5 R^3 B_inf.  It is divergence-free,
 * approaches B_inf, and gives exactly zero line potential (hence zero CT
 * normal flux) on the conducting r=R surface.
 */
static int conducting_sphere_potential(gamera_no_vec3 point, void *context,
                                       gamera_no_vec3 *potential) {
  (void)context;
  if (potential == NULL) {
    return -1;
  }
  const double radius = sqrt(point.value[0] * point.value[0] +
                             point.value[1] * point.value[1] +
                             point.value[2] * point.value[2]);
  if (!(radius > 0.0)) {
    return -1;
  }
  *potential = (gamera_no_vec3){{0.0, 0.0, 0.0}};
#ifdef GAMERA_EARTH_UPSTREAM_STARTUP
  /*
   * Initialize total B before subtracting B0 face fluxes: a radially cut
   * dipole around Earth plus an IMF front at x=-30 RE.
   */
  const double dipole_weight =
      cubic_ramp_down(radius, earth_cut_radius, earth_cut_length);
  const gamera_no_vec3 moment_cross_position = {{
      earth_dipole.moment.value[1] * point.value[2] -
          earth_dipole.moment.value[2] * point.value[1],
      earth_dipole.moment.value[2] * point.value[0] -
          earth_dipole.moment.value[0] * point.value[2],
      earth_dipole.moment.value[0] * point.value[1] -
          earth_dipole.moment.value[1] * point.value[0]}};
  const double inverse_radius_cubed = 1.0 / (radius * radius * radius);
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    potential->value[component] =
        dipole_weight * inverse_radius_cubed *
        moment_cross_position.value[component];
  }
  /* A constant normal Bx must cross a planar front continuously. */
  potential->value[1] -= 0.5 * initial_magnetic.value[0] * point.value[2];
  potential->value[2] += 0.5 * initial_magnetic.value[0] * point.value[1];
  const double front_integral =
      startup_front_antiderivative(point.value[0]);
  potential->value[1] += initial_magnetic.value[2] * front_integral;
  potential->value[2] -= initial_magnetic.value[1] * front_integral;
  return 0;
#endif
  if (radius > obstacle_radius * (1.0 + 8.0e-15)) {
    const double factor =
        0.5 * (1.0 - pow(obstacle_radius / radius, 3.0));
    potential->value[0] +=
        factor * (initial_magnetic.value[1] * point.value[2] -
                  initial_magnetic.value[2] * point.value[1]);
    potential->value[1] +=
        factor * (initial_magnetic.value[2] * point.value[0] -
                  initial_magnetic.value[0] * point.value[2]);
    potential->value[2] +=
        factor * (initial_magnetic.value[0] * point.value[1] -
                  initial_magnetic.value[1] * point.value[0]);
  }
#ifdef GAMERA_CRUSTAL_FIELD
  const gamera_no_vec3 crustal = crustal_vector_potential(point);
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    potential->value[component] += crustal.value[component];
  }
#endif
  return 0;
}

static gamera_no_vec3 analytic_magnetic(gamera_no_vec3 point) {
  const double radius = sqrt(point.value[0] * point.value[0] +
                             point.value[1] * point.value[1] +
                             point.value[2] * point.value[2]);
  const double inverse_radius = 1.0 / radius;
  const gamera_no_vec3 normal = {{point.value[0] * inverse_radius,
                                  point.value[1] * inverse_radius,
                                  point.value[2] * inverse_radius}};
  const double dipole_scale =
      -0.5 * pow(obstacle_radius * inverse_radius, 3.0);
  double b_dot_n = 0.0;
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    b_dot_n += initial_magnetic.value[component] * normal.value[component];
  }
  gamera_no_vec3 result;
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    result.value[component] =
        initial_magnetic.value[component] +
        dipole_scale *
            (3.0 * normal.value[component] * b_dot_n -
             upstream_magnetic.value[component]);
  }
#ifdef GAMERA_CRUSTAL_FIELD
  const gamera_no_vec3 crustal = crustal_magnetic(point);
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    result.value[component] += crustal.value[component];
  }
#endif
  return result;
}

void set_problem_config(void) {
  patch_count = 2;
  x1min_global = obstacle_radius;
  x1max_global = outer_radius;
  x2min_global = 0.25 * PI - angular_extension;
  x2max_global = 0.75 * PI + angular_extension;
  x3min_global = -0.75 * PI - angular_extension;
  x3max_global = 0.75 * PI + angular_extension;
  CA = 10.0;
#ifdef GAMERA_EARTH_DIPOLE_BACKGROUND
  rho_floor = earth_density_floor;
  p_floor = earth_pressure_floor;
#endif
  doRingAverage = 0;
  problem_config.Boundary_i = BC_NON_PERIODIC;
  problem_config.Boundary_j = BC_NON_PERIODIC;
  problem_config.Boundary_k = BC_NON_PERIODIC;
}

int problem_runtime_init(void) {
  gamera_solar_wind_init(&time_dependent_wind);
  time_dependent_wind_ready = 0;
  initial_density = inflow_density;
  initial_pressure = inflow_density * inflow_speed * inflow_speed /
                     (gamma_val * inflow_mach * inflow_mach);
  initial_velocity = (gamera_no_vec3){{inflow_speed, 0.0, 0.0}};
  initial_magnetic = upstream_magnetic;
  log_info("Yin-Yang angular extension %.6g deg; cell widths "
           "dtheta=%.6g deg dphi=%.6g deg on patch %d",
           (0.25 * PI - x2min_global) * 180.0 / PI,
           (x2max_global - x2min_global) * 180.0 /
               (PI * (double)config.nj_global),
           (x3max_global - x3min_global) * 180.0 /
               (PI * (double)config.nk_global),
           patch_id);
#ifdef GAMERA_EARTH_DIPOLE_BACKGROUND
  log_info("Earth magnetosphere floors rho=%.9g p=%.9g code "
           "(accepted production values)",
           rho_floor, p_floor);
  if (!(norm_config.B_Norm > 0.0) || !isfinite(norm_config.B_Norm)) {
    log_error("Earth background dipole requires a finite positive B_Norm");
    return -1;
  }
  earth_dipole.moment = (gamera_no_vec3){{
      0.0, 0.0,
      -earth_equatorial_surface_field_tesla / norm_config.B_Norm}};
  log_info("Earth background dipole moment Mz=%.9g code units "
           "(surface equatorial field 0.31 G)",
           earth_dipole.moment.value[2]);
#endif
  if (!wind_config.enabled) {
    log_info("Using constant analytic solar-wind boundary on patch %d",
             patch_id);
    return 0;
  }

  const gamera_solar_wind_hdf5_units units = {
      .physical_units = wind_config.physical_units,
      .time_norm = norm_config.Time_Norm,
      .density_norm = norm_config.rho_Norm,
      .velocity_norm = norm_config.u_Norm,
      .pressure_norm = norm_config.p_Norm,
      .magnetic_norm = norm_config.B_Norm,
      .velocity_si_scale = wind_config.velocity_si_scale};
  if (gamera_solar_wind_load_hdf5(wind_config.file, &units,
                                  &time_dependent_wind) != 0) {
    log_error("Unable to load solar-wind time series %s", wind_config.file);
    gamera_solar_wind_destroy(&time_dependent_wind);
    return -1;
  }
  const int file_relation = time_dependent_wind.enforce_bx_relation;
  for (int component = 0; component < GAMERA_WIND_DIM; ++component) {
    time_dependent_wind.reference[component] =
        wind_config.reference[component];
  }
  time_dependent_wind.time_offset = wind_config.time_offset;
  time_dependent_wind.linear_interpolation =
      wind_config.linear_interpolation;
  if (wind_config.coefficient_override) {
    time_dependent_wind.by_coefficient = wind_config.by_coefficient;
    time_dependent_wind.bz_coefficient = wind_config.bz_coefficient;
    time_dependent_wind.bx_offset = wind_config.bx_offset;
  }
  time_dependent_wind.enforce_bx_relation =
      wind_config.enforce_bx_relation >= 0
          ? wind_config.enforce_bx_relation
          : (wind_config.coefficient_override ? 1 : file_relation);
  gamera_solar_wind_apply_bx_relation(&time_dependent_wind);

  gamera_solar_wind_state initial;
  if (gamera_solar_wind_sample_time(
          &time_dependent_wind,
          time_sim + time_dependent_wind.time_offset, &initial) != 0) {
    log_error("Unable to sample initial solar-wind state");
    gamera_solar_wind_destroy(&time_dependent_wind);
    return -1;
  }
  initial_density = initial.density;
  initial_pressure = initial.pressure;
  for (int component = 0; component < GAMERA_WIND_DIM; ++component) {
    initial_velocity.value[component] = initial.velocity[component];
    initial_magnetic.value[component] = initial.magnetic[component];
  }
  time_dependent_wind_ready = 1;
  log_info("Loaded %zu solar-wind samples from %s: t=[%.9g, %.9g], "
           "interpolation=%s, front normal=(1,%.6g,%.6g), Bx relation=%s",
           time_dependent_wind.count, wind_config.file,
           time_dependent_wind.time[0],
           time_dependent_wind.time[time_dependent_wind.count - 1U],
           time_dependent_wind.linear_interpolation ? "linear" : "step",
           -time_dependent_wind.by_coefficient,
           -time_dependent_wind.bz_coefficient,
           time_dependent_wind.enforce_bx_relation ? "on" : "off");
  return 0;
}

void problem_runtime_finalize(void) {
  gamera_solar_wind_destroy(&time_dependent_wind);
  time_dependent_wind_ready = 0;
}

int problem_grid_init(void) {
#ifdef GAMERA_BOW_SHOCK
  if (!isfinite(obstacle_radius) || !isfinite(outer_radius) ||
      !isfinite(radial_stretch) || !(obstacle_radius > 0.0) ||
      !(outer_radius > obstacle_radius) || !(radial_stretch > 0.0) ||
      ensure_radial_map() != 0) {
    log_error("Invalid stretched radial map: inner=%.17g outer=%.17g "
              "version=%d legacy_stretch=%.17g",
              obstacle_radius, outer_radius, requested_radial_map_version(),
              radial_stretch);
    return -1;
  }
#endif
  for (int i = 0; i < config.NI; ++i) {
    const int global_i = proc_coords[0] * config.ni + i - NG;
#ifdef GAMERA_BOW_SHOCK
    /*
     * Exponential vertex map.  Extending the same analytic map into the
     * geometry halos keeps the metric smooth and leaves all inner ghost
     * radii positive for the supported coarse grids.
     */
    const double radial_fraction =
        (double)global_i / (double)config.ni_global;
    double radius = 0.0;
    if (gamera_no_radial_map_forward(&radial_map, radial_fraction,
                                     &radius) != 0) {
      log_error("Failed to evaluate radial map at logical fraction %.17g",
                radial_fraction);
      return -1;
    }
#else
    const double radius =
        x1min_global + (x1max_global - x1min_global) *
                           (double)global_i / (double)config.ni_global;
#endif
    if (!(radius > 0.0)) {
      log_error("Yin-Yang IMF-obstacle radial ghost vertex is non-positive");
      return -1;
    }
    for (int j = 0; j < config.NJ; ++j) {
      const int global_j = proc_coords[1] * config.nj + j - NG;
      const double theta =
          x2min_global + (x2max_global - x2min_global) *
                             (double)global_j / (double)config.nj_global;
      for (int k = 0; k < config.NK; ++k) {
        const int global_k = proc_coords[2] * config.nk + k - NG;
        const double phi =
            x3min_global + (x3max_global - x3min_global) *
                               (double)global_k / (double)config.nk_global;
        gamera_no_vec3 point;
        if (gamera_no_yinyang_logical_to_global(
                patch_id, radius, theta, phi, &point) != 0) {
          return -1;
        }
        x1[i][j][k] = point.value[0];
        x2[i][j][k] = point.value[1];
        x3[i][j][k] = point.value[2];
      }
    }
  }
  return 0;
}

#ifdef GAMERA_BOW_SHOCK
int problem_nonorthogonal_radial_logical(double radius, double *logical) {
  if (logical == NULL || !isfinite(radius) || ensure_radial_map() != 0) {
    return -1;
  }
  double logical_fraction = 0.0;
  if (gamera_no_radial_map_inverse(&radial_map, radius,
                                   &logical_fraction) != 0) {
    return -1;
  }
  *logical = (double)config.ni_global * logical_fraction - 0.5;
  return isfinite(*logical) ? 0 : -1;
}

int problem_nonorthogonal_radial_map_version(void) {
  return requested_radial_map_version();
}

double problem_nonorthogonal_radial_stretch(void) {
  return radial_stretch;
}

int problem_nonorthogonal_radial_map_parameters(
    double parameters[GAMERA_NO_RADIAL_MAP_PARAMETER_COUNT]) {
  return ensure_radial_map() == 0
             ? gamera_no_radial_map_parameters(&radial_map, parameters)
             : -1;
}
#endif

static void set_inflow_cell(int i, int j, int k, int old_state) {
  const int offset = old_state ? gas_rho_p - gas_rho : 0;
  gas[0][gas_rho + offset][i][j][k] = inflow_density;
  gas[0][gas_v1 + offset][i][j][k] = inflow_speed;
  gas[0][gas_v2 + offset][i][j][k] = 0.0;
  gas[0][gas_v3 + offset][i][j][k] = 0.0;
  gas[0][gas_p + offset][i][j][k] = inflow_pressure();
  gas[0][gas_p_S + offset][i][j][k] = inflow_pressure();
}

static void copy_outer_outflow(int destination, int j, int k,
                               int old_state) {
  const int offset = old_state ? gas_rho_p - gas_rho : 0;
  for (int field = gas_rho; field <= gas_p_S; ++field) {
    gas[0][field + offset][destination][j][k] =
        gas[0][field + offset][ie][j][k];
  }
}

static double old_boundary_time(void) {
  return isfinite(dt0) && dt0 > 0.0 ? time_sim - dt0 : time_sim;
}

static int wind_state_and_weight(gamera_no_vec3 point, double time,
                                 gamera_solar_wind_state *state,
                                 double *weight) {
  const double position[3] = {point.value[0], point.value[1], point.value[2]};
  return gamera_solar_wind_sample_at(&time_dependent_wind, position, time,
                                     state, NULL) == 0 &&
                 gamera_solar_wind_weight(&time_dependent_wind, position,
                                          time, weight) == 0
             ? 0
             : -1;
}

static int set_time_dependent_outer_cell(int destination, int j, int k,
                                         int old_state, double time) {
  const int offset = old_state ? gas_rho_p - gas_rho : 0;
  const gamera_no_vec3 point = {{x1ctr[destination][j][k],
                                 x2ctr[destination][j][k],
                                 x3ctr[destination][j][k]}};
  gamera_solar_wind_state wind;
  double wind_weight;
  if (wind_state_and_weight(point, time, &wind, &wind_weight) != 0) {
    return -1;
  }

  double out_density = gas[0][gas_rho + offset][ie][j][k];
  double out_pressure = gas[0][gas_p + offset][ie][j][k];
  double out_velocity[3] = {
      gas[0][gas_v1 + offset][ie][j][k],
      gas[0][gas_v2 + offset][ie][j][k],
      gas[0][gas_v3 + offset][ie][j][k]};
  const double radius = sqrt(point.value[0] * point.value[0] +
                             point.value[1] * point.value[1] +
                             point.value[2] * point.value[2]);
  if (!(radius > 0.0) || !(out_density > 0.0) || !(out_pressure > 0.0)) {
    return -1;
  }
  const double normal[3] = {point.value[0] / radius,
                            point.value[1] / radius,
                            point.value[2] / radius};
  const double normal_velocity = out_velocity[0] * normal[0] +
                                 out_velocity[1] * normal[1] +
                                 out_velocity[2] * normal[2];
  const double sound_speed = sqrt(gamma_val * out_pressure / out_density);
  if (normal_velocity < sound_speed) {
    for (int component = 0; component < 3; ++component) {
      out_velocity[component] +=
          (sound_speed - normal_velocity) * normal[component];
    }
  }

  const double out_weight = 1.0 - wind_weight;
  gas[0][gas_rho + offset][destination][j][k] =
      wind_weight * wind.density + out_weight * out_density;
  gas[0][gas_p + offset][destination][j][k] =
      wind_weight * wind.pressure + out_weight * out_pressure;
  gas[0][gas_p_S + offset][destination][j][k] =
      gas[0][gas_p + offset][destination][j][k];
  for (int component = 0; component < 3; ++component) {
    gas[0][gas_v1 + component + offset][destination][j][k] =
        wind_weight * wind.velocity[component] +
        out_weight * out_velocity[component];
  }
  return 0;
}

static void reflect_inner_cell(int destination, int source, int j, int k,
                               int old_state) {
  const int offset = old_state ? gas_rho_p - gas_rho : 0;
  const double x = x1ctr[destination][j][k];
  const double y = x2ctr[destination][j][k];
  const double z = x3ctr[destination][j][k];
  const double inverse_radius = 1.0 / sqrt(x * x + y * y + z * z);
  const double normal[3] = {x * inverse_radius, y * inverse_radius,
                            z * inverse_radius};
  const double velocity[3] = {
      gas[0][gas_v1 + offset][source][j][k],
      gas[0][gas_v2 + offset][source][j][k],
      gas[0][gas_v3 + offset][source][j][k]};
  const double normal_velocity = velocity[0] * normal[0] +
                                 velocity[1] * normal[1] +
                                 velocity[2] * normal[2];
  gas[0][gas_rho + offset][destination][j][k] =
      gas[0][gas_rho + offset][source][j][k];
  gas[0][gas_p + offset][destination][j][k] =
      gas[0][gas_p + offset][source][j][k];
  gas[0][gas_p_S + offset][destination][j][k] =
      gas[0][gas_p_S + offset][source][j][k];
#ifdef GAMERA_MI_COUPLING
  if (gamera_mi_coupling_ready()) {
    const gamera_no_vec3 point = {{x, y, z}};
    const gamera_no_vec3 source_velocity = {
        {velocity[0], velocity[1], velocity[2]}};
    gamera_no_vec3 ghost_velocity;
    if (gamera_mi_coupling_ghost_velocity(
            point, source_velocity, &ghost_velocity) == 0) {
#ifdef GAMERA_STRICT_INNER_GAS_WALL
      /* The helper returns 2*V_EB-v_active for the mirrored-ghost convention.
       * Recover the directly stored wall state V_EB for strict-wall mode. */
      for (int component = 0; component < GAMERA_NO_DIM; ++component) {
        ghost_velocity.value[component] =
            0.5 * (ghost_velocity.value[component] + velocity[component]);
      }
#endif
      gas[0][gas_v1 + offset][destination][j][k] =
          ghost_velocity.value[0];
      gas[0][gas_v2 + offset][destination][j][k] =
          ghost_velocity.value[1];
      gas[0][gas_v3 + offset][destination][j][k] =
          ghost_velocity.value[2];
      return;
    }
  }
#endif
#ifdef GAMERA_STRICT_INNER_GAS_WALL
  gas[0][gas_v1 + offset][destination][j][k] = 0.0;
  gas[0][gas_v2 + offset][destination][j][k] = 0.0;
  gas[0][gas_v3 + offset][destination][j][k] = 0.0;
  return;
#endif
  gas[0][gas_v1 + offset][destination][j][k] =
      velocity[0] - 2.0 * normal_velocity * normal[0];
  gas[0][gas_v2 + offset][destination][j][k] =
      velocity[1] - 2.0 * normal_velocity * normal[1];
  gas[0][gas_v3 + offset][destination][j][k] =
      velocity[2] - 2.0 * normal_velocity * normal[2];
}

void boundary_conditions(void) {
  for (int i = 0; i < config.NI - 1; ++i) {
    if (i >= is && i <= ie) {
      continue;
    }
    const int inner = i < is;
    if ((inner && proc_coords[0] != 0) ||
        (!inner && proc_coords[0] != config.proc_dims[0] - 1)) {
      continue;
    }
    const int inner_source = 2 * is - 1 - i;
    for (int j = 0; j < config.NJ - 1; ++j) {
      for (int k = 0; k < config.NK - 1; ++k) {
        if (inner) {
          reflect_inner_cell(i, inner_source, j, k, 0);
          reflect_inner_cell(i, inner_source, j, k, 1);
          continue;
        }
        if (time_dependent_wind_ready) {
          if (set_time_dependent_outer_cell(i, j, k, 0, time_sim) != 0 ||
              set_time_dependent_outer_cell(i, j, k, 1,
                                            old_boundary_time()) != 0) {
            log_error("Failed to sample time-dependent wind at outer ghost "
                      "(%d,%d,%d)",
                      i, j, k);
          }
          continue;
        }
        const double x = x1ctr[i][j][k];
        const double y = x2ctr[i][j][k];
        const double z = x3ctr[i][j][k];
        const double outward_x = x / sqrt(x * x + y * y + z * z);
        if (inflow_speed * outward_x < 0.0) {
          set_inflow_cell(i, j, k, 0);
          set_inflow_cell(i, j, k, 1);
        } else {
          copy_outer_outflow(i, j, k, 0);
          copy_outer_outflow(i, j, k, 1);
        }
      }
    }
  }
}

static gamera_no_vec3 wind_edge_midpoint(const gamera_no_grid *grid,
                                         int direction, size_t i, size_t j,
                                         size_t k) {
  size_t upper[3] = {i, j, k};
  ++upper[direction];
  const gamera_no_vec3 start =
      grid->vertex[gamera_no_index3(grid->vertex_extent, i, j, k)];
  const gamera_no_vec3 end = grid->vertex[gamera_no_index3(
      grid->vertex_extent, upper[0], upper[1], upper[2])];
  gamera_no_vec3 midpoint;
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    midpoint.value[component] =
        0.5 * (start.value[component] + end.value[component]);
  }
  return midpoint;
}

static int face_flux_density(const gamera_no_grid *grid, const double *flux,
                             int direction, size_t i, size_t j, size_t k,
                             double *result) {
  if (grid == NULL || flux == NULL || result == NULL ||
      i >= grid->face[direction].extent[0] ||
      j >= grid->face[direction].extent[1] ||
      k >= grid->face[direction].extent[2]) {
    return -1;
  }
  const size_t face =
      gamera_no_index3(grid->face[direction].extent, i, j, k);
  const double area = grid->face[direction].value[face].area;
  if (!(area > 0.0)) {
    return -1;
  }
  *result = flux[face] / area;
  return isfinite(*result) ? 0 : -1;
}

/* Prevent transverse magnetic field from hanging on the outer shell. */
static int add_wind_diffusive_emf(
    gamera_no_storage *storage, const gamera_no_grid *grid,
    const size_t active_lower[3], const size_t active_upper[3], double time) {
  const size_t i = active_upper[GAMERA_NO_I];
  if (i == 0U || !(dt > 0.0)) {
    return -1;
  }
  int failed = 0;
  for (size_t j = active_lower[GAMERA_NO_J];
       j <= active_upper[GAMERA_NO_J]; ++j) {
    for (size_t k = active_lower[GAMERA_NO_K];
         k <= active_upper[GAMERA_NO_K]; ++k) {
      if (j == 0U || k == 0U ||
          i >= grid->edge[GAMERA_NO_J].extent[0] ||
          j >= grid->edge[GAMERA_NO_J].extent[1] ||
          k >= grid->edge[GAMERA_NO_J].extent[2] ||
          i >= grid->edge[GAMERA_NO_K].extent[0] ||
          j >= grid->edge[GAMERA_NO_K].extent[1] ||
          k >= grid->edge[GAMERA_NO_K].extent[2]) {
        continue;
      }
      const size_t edge_j = gamera_no_index3(
          grid->edge[GAMERA_NO_J].extent, i, j, k);
      const size_t edge_k = gamera_no_index3(
          grid->edge[GAMERA_NO_K].extent, i, j, k);
      if (grid->edge[GAMERA_NO_J].valid[edge_j] == 0U ||
          grid->edge[GAMERA_NO_K].valid[edge_k] == 0U) {
        continue;
      }
      double bi_k_plus, bi_k_minus, bk_i_plus, bk_i_minus;
      double bj_i_plus, bj_i_minus, bi_j_plus, bi_j_minus;
      if (face_flux_density(
              grid, storage->predicted_face_flux[GAMERA_NO_I], GAMERA_NO_I,
              i, j, k, &bi_k_plus) != 0 ||
          face_flux_density(
              grid, storage->predicted_face_flux[GAMERA_NO_I], GAMERA_NO_I,
              i, j, k - 1U, &bi_k_minus) != 0 ||
          face_flux_density(
              grid, storage->predicted_face_flux[GAMERA_NO_K], GAMERA_NO_K,
              i, j, k, &bk_i_plus) != 0 ||
          face_flux_density(
              grid, storage->predicted_face_flux[GAMERA_NO_K], GAMERA_NO_K,
              i - 1U, j, k, &bk_i_minus) != 0 ||
          face_flux_density(
              grid, storage->predicted_face_flux[GAMERA_NO_J], GAMERA_NO_J,
              i, j, k, &bj_i_plus) != 0 ||
          face_flux_density(
              grid, storage->predicted_face_flux[GAMERA_NO_J], GAMERA_NO_J,
              i - 1U, j, k, &bj_i_minus) != 0 ||
          face_flux_density(
              grid, storage->predicted_face_flux[GAMERA_NO_I], GAMERA_NO_I,
              i, j, k, &bi_j_plus) != 0 ||
          face_flux_density(
              grid, storage->predicted_face_flux[GAMERA_NO_I], GAMERA_NO_I,
              i, j - 1U, k, &bi_j_minus) != 0) {
        failed = 1;
        continue;
      }
      const double current_j =
          (bi_k_plus - bi_k_minus) - (bk_i_plus - bk_i_minus);
      const double current_k =
          (bj_i_plus - bj_i_minus) - (bi_j_plus - bi_j_minus);
      const double dl =
          sqrt(grid->edge[GAMERA_NO_J].value[edge_j].length *
               grid->edge[GAMERA_NO_K].value[edge_k].length);
      const gamera_no_vec3 point_j =
          wind_edge_midpoint(grid, GAMERA_NO_J, i, j, k);
      const gamera_no_vec3 point_k =
          wind_edge_midpoint(grid, GAMERA_NO_K, i, j, k);
      gamera_solar_wind_state wind_j, wind_k;
      double weight_j, weight_k;
      if (!(dl > 0.0) ||
          wind_state_and_weight(point_j, time, &wind_j, &weight_j) != 0 ||
          wind_state_and_weight(point_k, time, &wind_k, &weight_k) != 0) {
        failed = 1;
        continue;
      }
      const double coefficient_square =
          time_dependent_wind.by_coefficient *
              time_dependent_wind.by_coefficient +
          time_dependent_wind.bz_coefficient *
              time_dependent_wind.bz_coefficient;
      const double front_speed =
          fmax(fabs(wind_j.velocity[GAMERA_WIND_X]),
               fabs(wind_k.velocity[GAMERA_WIND_X])) /
          sqrt(1.0 + coefficient_square);
      const double diffusion_speed = fmin(front_speed, CFL * dl / dt);
      storage->edge_emf[GAMERA_NO_J][edge_j] +=
          weight_j * diffusion_speed * current_j * dl;
      storage->edge_emf[GAMERA_NO_K][edge_k] +=
          weight_k * diffusion_speed * current_k * dl;
    }
  }
  return failed ? -1 : 0;
}

static int impose_time_dependent_wind_emf(
    gamera_no_storage *storage, const gamera_no_grid *grid,
    const size_t active_lower[3], const size_t active_upper[3], double time) {
  if (!time_dependent_wind_ready ||
      proc_coords[0] != config.proc_dims[0] - 1) {
    return 0;
  }
  const size_t outer_i = active_upper[GAMERA_NO_I];
  const size_t first_i = outer_i > 3U ? outer_i - 3U : 0U;
  int failed = 0;
  for (int direction = GAMERA_NO_I; direction <= GAMERA_NO_K; ++direction) {
    const size_t last_i =
        outer_i < grid->edge[direction].extent[0]
            ? outer_i + 1U
            : grid->edge[direction].extent[0];
#pragma omp parallel for collapse(3) reduction(| : failed) schedule(static)
    for (size_t i = first_i; i < last_i; ++i) {
      for (size_t j = 0; j < grid->edge[direction].extent[1]; ++j) {
        for (size_t k = 0; k < grid->edge[direction].extent[2]; ++k) {
          const size_t edge = gamera_no_index3(
              grid->edge[direction].extent, i, j, k);
          if (grid->edge[direction].valid[edge] == 0U) {
            continue;
          }
          const gamera_no_vec3 point =
              wind_edge_midpoint(grid, direction, i, j, k);
          gamera_solar_wind_state wind;
          double wind_weight;
          if (wind_state_and_weight(point, time, &wind, &wind_weight) != 0) {
            failed = 1;
            continue;
          }
          if (i < outer_i) {
            wind_weight /= 1.0 + (double)(outer_i - i);
          }
          const double electric[3] = {
              -(wind.velocity[1] * wind.magnetic[2] -
                wind.velocity[2] * wind.magnetic[1]),
              -(wind.velocity[2] * wind.magnetic[0] -
                wind.velocity[0] * wind.magnetic[2]),
              -(wind.velocity[0] * wind.magnetic[1] -
                wind.velocity[1] * wind.magnetic[0])};
          const gamera_no_edge_geometry *geometry =
              &grid->edge[direction].value[edge];
          double line_integral = 0.0;
          for (int component = 0; component < GAMERA_NO_DIM; ++component) {
            line_integral += electric[component] *
                             geometry->normal.value[component] *
                             geometry->length;
          }
          storage->edge_emf[direction][edge] =
              (1.0 - wind_weight) * storage->edge_emf[direction][edge] +
              wind_weight * line_integral;
        }
      }
    }
  }
  if (!failed && add_wind_diffusive_emf(storage, grid, active_lower,
                                        active_upper, time) != 0) {
    failed = 1;
  }
  return failed ? -1 : 0;
}

int problem_nonorthogonal_edge_emf_boundary(
    gamera_no_storage *storage, const gamera_no_grid *grid,
    const size_t active_lower[3], const size_t active_upper[3]) {
  if (storage == NULL || grid == NULL) {
    return -1;
  }
  if (impose_time_dependent_wind_emf(
          storage, grid, active_lower, active_upper,
          time_sim + 0.5 * dt) != 0) {
    return -1;
  }
  if (proc_coords[0] != 0) {
    return 0;
  }
  const size_t wall_i = active_lower[GAMERA_NO_I];
  for (int direction = GAMERA_NO_J; direction <= GAMERA_NO_K; ++direction) {
    for (size_t j = 0; j < grid->edge[direction].extent[1]; ++j) {
      for (size_t k = 0; k < grid->edge[direction].extent[2]; ++k) {
        const size_t edge = gamera_no_index3(
            grid->edge[direction].extent, wall_i, j, k);
        if (grid->edge[direction].valid[edge] != 0U) {
          storage->edge_emf[direction][edge] = 0.0;
        }
      }
    }
  }
  return 0;
}

#ifdef GAMERA_NONORTHOGONAL_HAS_FLUID_FLUX_BOUNDARY
int problem_nonorthogonal_fluid_flux_boundary(
    gamera_no_storage *storage, const gamera_no_grid *grid,
    const size_t active_lower[3], const size_t active_upper[3]) {
  if (storage == NULL || grid == NULL) {
    return -1;
  }
  if (proc_coords[0] != 0) {
    return 0;
  }
  const size_t i = active_lower[GAMERA_NO_I];
  for (size_t j = active_lower[GAMERA_NO_J];
       j < active_upper[GAMERA_NO_J]; ++j) {
    for (size_t k = active_lower[GAMERA_NO_K];
         k < active_upper[GAMERA_NO_K]; ++k) {
      const size_t face = gamera_no_index3(
          storage->face_extent[GAMERA_NO_I], i, j, k);
      const gamera_no_face_geometry *geometry =
          &grid->face[GAMERA_NO_I].value[face];
      if (!(geometry->area > 0.0)) {
        return -1;
      }
      const int hemisphere = geometry->centroid.value[2] >= 0.0 ? 0 : 1;
      const double *fluid =
          &storage->fluid_face_flux[GAMERA_NO_I]
                                   [face * GAMERA_NO_FLUX_COUNT];
      const double inverse_area = 1.0 / geometry->area;
      if (fluid[GAMERA_NO_FLUX_DENSITY] > 0.0) {
        ++storage->inner_wall_clamped_face_count[hemisphere];
        storage->inner_wall_positive_mass_max[hemisphere] =
            fmax(storage->inner_wall_positive_mass_max[hemisphere],
                 fluid[GAMERA_NO_FLUX_DENSITY] * inverse_area);
        storage->inner_wall_positive_energy_max[hemisphere] =
            fmax(storage->inner_wall_positive_energy_max[hemisphere],
                 fmax(0.0, fluid[GAMERA_NO_FLUX_ENERGY]) * inverse_area);
      }
      double fluid_momentum_square = 0.0;
      double maxwell_momentum_square = 0.0;
      for (int component = 0; component < GAMERA_NO_DIM; ++component) {
        const double fluid_value =
            fluid[GAMERA_NO_FLUX_MOMENTUM_X + component] * inverse_area;
        const double maxwell_value =
            storage->maxwell_face_flux[GAMERA_NO_I][face].value[component] *
            inverse_area;
        fluid_momentum_square += fluid_value * fluid_value;
        maxwell_momentum_square += maxwell_value * maxwell_value;
      }
      storage->inner_wall_fluid_momentum_max[hemisphere] =
          fmax(storage->inner_wall_fluid_momentum_max[hemisphere],
               sqrt(fluid_momentum_square));
      storage->inner_wall_maxwell_momentum_max[hemisphere] =
          fmax(storage->inner_wall_maxwell_momentum_max[hemisphere],
               sqrt(maxwell_momentum_square));

      const size_t cell =
          gamera_no_index3(storage->cell_extent, i, j, k);
      const size_t outer_cell =
          gamera_no_index3(storage->cell_extent, i + 1U, j, k);
      gamera_no_primitive primitive;
      gamera_no_primitive outer_primitive;
      if (gamera_no_conserved_to_primitive(
              &storage->conserved[cell * GAMERA_NO_FLUX_COUNT], gamma_val,
              rho_floor, p_floor, &primitive) != 0 ||
          gamera_no_conserved_to_primitive(
              &storage->conserved[outer_cell * GAMERA_NO_FLUX_COUNT],
              gamma_val, rho_floor, p_floor, &outer_primitive) != 0) {
        return -1;
      }
      double speed_square = 0.0;
      double residual_magnetic_square = 0.0;
      for (int component = 0; component < GAMERA_NO_DIM; ++component) {
        speed_square += primitive.velocity.value[component] *
                        primitive.velocity.value[component];
        residual_magnetic_square +=
            storage->cell_magnetic[cell].value[component] *
            storage->cell_magnetic[cell].value[component];
      }
      const gamera_no_vec3 center = grid->cell[cell].centroid;
      const gamera_no_vec3 outer_center = grid->cell[outer_cell].centroid;
      const double radius = sqrt(center.value[0] * center.value[0] +
                                 center.value[1] * center.value[1] +
                                 center.value[2] * center.value[2]);
      const double outer_radius_local =
          sqrt(outer_center.value[0] * outer_center.value[0] +
               outer_center.value[1] * outer_center.value[1] +
               outer_center.value[2] * outer_center.value[2]);
      const double radial_spacing = fabs(outer_radius_local - radius);
      storage->inner_wall_density_max[hemisphere] =
          fmax(storage->inner_wall_density_max[hemisphere],
               primitive.density);
      storage->inner_wall_pressure_max[hemisphere] =
          fmax(storage->inner_wall_pressure_max[hemisphere],
               primitive.pressure);
      if (radial_spacing > 0.0) {
        storage->inner_wall_pressure_gradient_max[hemisphere] = fmax(
            storage->inner_wall_pressure_gradient_max[hemisphere],
            fabs(outer_primitive.pressure - primitive.pressure) /
                radial_spacing);
      }
      storage->inner_wall_speed_max[hemisphere] =
          fmax(storage->inner_wall_speed_max[hemisphere],
               sqrt(speed_square));
      const double residual_magnetic_norm =
          sqrt(residual_magnetic_square);
      if (residual_magnetic_norm >
          storage->inner_wall_residual_magnetic_max[hemisphere]) {
        storage->inner_wall_residual_magnetic_max[hemisphere] =
            residual_magnetic_norm;
        storage->inner_wall_residual_magnetic_location[hemisphere] = center;
      }
      /* The ionospheric low-latitude boundary maps to the magnetic-equatorial
       * part of the inner MHD shell.  Keep a separate +/-15 degree diagnostic
       * so polar/plasma-sheet maxima cannot hide a boundary-localized pulse. */
      if (fabs(center.value[2]) <= radius * sin(15.0 * PI / 180.0)) {
        storage->inner_wall_lowlat_pressure_max[hemisphere] =
            fmax(storage->inner_wall_lowlat_pressure_max[hemisphere],
                 primitive.pressure);
        if (radial_spacing > 0.0) {
          storage->inner_wall_lowlat_pressure_gradient_max[hemisphere] = fmax(
              storage->inner_wall_lowlat_pressure_gradient_max[hemisphere],
              fabs(outer_primitive.pressure - primitive.pressure) /
                  radial_spacing);
        }
        storage->inner_wall_lowlat_speed_max[hemisphere] =
            fmax(storage->inner_wall_lowlat_speed_max[hemisphere],
                 sqrt(speed_square));
        storage->inner_wall_lowlat_residual_magnetic_max[hemisphere] = fmax(
            storage->inner_wall_lowlat_residual_magnetic_max[hemisphere],
            residual_magnetic_norm);
      }
    }
  }
  return gamera_no_trap_inner_outward_mass_energy_flux(
      storage, active_lower, active_upper, NULL);
}
#endif

static gamera_no_vec3 reflected_magnetic(gamera_no_vec3 magnetic,
                                         gamera_no_vec3 point) {
  const double radius = sqrt(point.value[0] * point.value[0] +
                             point.value[1] * point.value[1] +
                             point.value[2] * point.value[2]);
  double normal_component = 0.0;
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    normal_component +=
        magnetic.value[component] * point.value[component] / radius;
  }
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    magnetic.value[component] -=
        2.0 * normal_component * point.value[component] / radius;
  }
  return magnetic;
}

static int time_dependent_wind_face_flux(
    const gamera_no_storage *storage, const gamera_no_grid *grid,
    int direction, size_t i, size_t j, size_t k, int old_state, double time,
    double *result) {
  if (storage == NULL || grid == NULL || result == NULL) {
    return -1;
  }
  const size_t face = gamera_no_index3(grid->face[direction].extent, i, j, k);
  const gamera_no_face_geometry *geometry =
      &grid->face[direction].value[face];
  gamera_solar_wind_state wind;
  double wind_weight;
  if (wind_state_and_weight(geometry->centroid, time, &wind, &wind_weight) !=
      0) {
    return -1;
  }
  double wind_flux = 0.0;
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    wind_flux += wind.magnetic[component] *
                 geometry->area_vector.value[component];
  }
#ifdef GAMERA_EARTH_DIPOLE_BACKGROUND
  const gamera_no_background_field *background =
      gamera_no_driver_background_field();
  if (background == NULL || background->face_flux[direction] == NULL) {
    return -1;
  }
  /* The evolved CT field is residual: dB = B_total - B0. */
  wind_flux -= background->face_flux[direction][face];
#endif
  size_t source_i = direction == GAMERA_NO_I ? (size_t)(ie + 1)
                                              : (size_t)ie;
  if (source_i >= grid->face[direction].extent[0]) {
    source_i = grid->face[direction].extent[0] - 1U;
  }
  const size_t source = gamera_no_index3(grid->face[direction].extent,
                                         source_i, j, k);
  const double outflow_flux = old_state
                                  ? storage->old_face_flux[direction][source]
                                  : storage->face_flux[direction][source];
  *result = wind_weight * wind_flux + (1.0 - wind_weight) * outflow_flux;
  return isfinite(*result) ? 0 : -1;
}

int problem_nonorthogonal_magnetic_boundary(gamera_no_storage *storage,
                                             const gamera_no_grid *grid) {
  if (storage == NULL || grid == NULL) {
    return -1;
  }
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    for (size_t i = 0; i < grid->face[direction].extent[0]; ++i) {
      const int inner = direction == GAMERA_NO_I ? i < (size_t)is
                                                 : i < (size_t)is;
      const int outer = direction == GAMERA_NO_I ? i > (size_t)(ie + 1)
                                                 : i > (size_t)ie;
      if ((!inner || proc_coords[0] != 0) &&
          (!outer || proc_coords[0] != config.proc_dims[0] - 1)) {
        continue;
      }
      for (size_t j = 0; j < grid->face[direction].extent[1]; ++j) {
        for (size_t k = 0; k < grid->face[direction].extent[2]; ++k) {
          const size_t face = gamera_no_index3(
              grid->face[direction].extent, i, j, k);
          const gamera_no_face_geometry *geometry =
              &grid->face[direction].value[face];
          if (!inner && time_dependent_wind_ready) {
            if (time_dependent_wind_face_flux(
                    storage, grid, direction, i, j, k, 0, time_sim,
                    &storage->face_flux[direction][face]) != 0 ||
                time_dependent_wind_face_flux(
                    storage, grid, direction, i, j, k, 1,
                    old_boundary_time(),
                    &storage->old_face_flux[direction][face]) != 0) {
              return -1;
            }
            continue;
          }
          gamera_no_vec3 field[2];
          if (inner) {
#ifdef GAMERA_STRICT_INNER_MAGNETIC_WALL
            /* Copy the first active-shell face-flux density into every radial
             * ghost with the face-area ratio. Preserve the CT face value
             * directly instead of reconstructing it from Cartesian B. */
            size_t source_coordinate[3] = {(size_t)is, j, k};
            for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
              if (source_coordinate[axis] >=
                  grid->face[direction].extent[axis]) {
                source_coordinate[axis] =
                    grid->face[direction].extent[axis] - 1U;
              }
            }
            const size_t source_face = gamera_no_index3(
                grid->face[direction].extent, source_coordinate[0],
                source_coordinate[1], source_coordinate[2]);
            const double source_area =
                grid->face[direction].value[source_face].area;
            if (!(source_area > 0.0)) {
              return -1;
            }
            const double area_ratio = geometry->area / source_area;
            storage->face_flux[direction][face] =
                area_ratio * storage->face_flux[direction][source_face];
            storage->old_face_flux[direction][face] =
                area_ratio * storage->old_face_flux[direction][source_face];
            continue;
#else
            size_t source[3] = {(size_t)(2 * is - 1 - (int)i), j, k};
            for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
              if (source[axis] >= grid->cell_extent[axis]) {
                source[axis] = grid->cell_extent[axis] - 1U;
              }
            }
            const size_t cell = gamera_no_index3(
                grid->cell_extent, source[0], source[1], source[2]);
#ifdef GAMERA_CRUSTAL_FIELD
            const gamera_no_vec3 source_crustal =
                crustal_magnetic(grid->cell[cell].centroid);
            const gamera_no_vec3 face_crustal =
                crustal_magnetic(geometry->centroid);
            field[0] = storage->cell_magnetic[cell];
            field[1] = storage->old_cell_magnetic[cell];
            for (int component = 0; component < GAMERA_NO_DIM; ++component) {
              field[0].value[component] -= source_crustal.value[component];
              field[1].value[component] -= source_crustal.value[component];
            }
            field[0] = reflected_magnetic(field[0], geometry->centroid);
            field[1] = reflected_magnetic(field[1], geometry->centroid);
            for (int component = 0; component < GAMERA_NO_DIM; ++component) {
              field[0].value[component] += face_crustal.value[component];
              field[1].value[component] += face_crustal.value[component];
            }
#else
#ifdef GAMERA_BOW_SHOCK
            /*
             * The strict magnetosphere wall copies Cartesian B into the
             * conjugate ghost and reconstructs every face flux from that
             * field. This is dB/dn=0, not a polar-vector reflection.
             */
            field[0] = storage->cell_magnetic[cell];
            field[1] = storage->old_cell_magnetic[cell];
#else
            field[0] = reflected_magnetic(storage->cell_magnetic[cell],
                                          geometry->centroid);
            field[1] = reflected_magnetic(storage->old_cell_magnetic[cell],
                                          geometry->centroid);
#endif
#endif
#endif
          } else {
            field[0] = analytic_magnetic(geometry->centroid);
            field[1] = field[0];
          }
          storage->face_flux[direction][face] =
              field[0].value[0] * geometry->area_vector.value[0] +
              field[0].value[1] * geometry->area_vector.value[1] +
              field[0].value[2] * geometry->area_vector.value[2];
          storage->old_face_flux[direction][face] =
              field[1].value[0] * geometry->area_vector.value[0] +
              field[1].value[1] * geometry->area_vector.value[1] +
              field[1].value[2] * geometry->area_vector.value[2];
        }
      }
    }
  }
  const double *current[3] = {storage->face_flux[0], storage->face_flux[1],
                              storage->face_flux[2]};
  const double *old[3] = {storage->old_face_flux[0],
                          storage->old_face_flux[1],
                          storage->old_face_flux[2]};
  if (gamera_no_recover_magnetic_field(grid, current,
                                       storage->cell_magnetic) != 0 ||
      gamera_no_recover_magnetic_field(grid, old,
                                       storage->old_cell_magnetic) != 0) {
    return -1;
  }
#ifdef GAMERA_BOW_SHOCK
  if (proc_coords[0] == 0) {
    /* Explicitly impose zero Cartesian-field gradient in predictor/ghost
     * storage as well as on the reconstructed face fluxes. */
    for (int i = 0; i < is; ++i) {
#ifdef GAMERA_STRICT_INNER_MAGNETIC_WALL
      const size_t source_i = (size_t)is;
#else
      const size_t source_i = (size_t)(2 * is - 1 - i);
#endif
      for (size_t j = 0; j < grid->cell_extent[1]; ++j) {
        for (size_t k = 0; k < grid->cell_extent[2]; ++k) {
          const size_t ghost = gamera_no_index3(
              grid->cell_extent, (size_t)i, j, k);
          const size_t source = gamera_no_index3(
              grid->cell_extent, source_i, j, k);
          storage->cell_magnetic[ghost] = storage->cell_magnetic[source];
          storage->old_cell_magnetic[ghost] =
              storage->old_cell_magnetic[source];
        }
      }
    }
  }
#endif
  return 0;
}

#ifdef GAMERA_IONOSPHERE_SOURCE_LOSS
int problem_nonorthogonal_cell_source(
    gamera_no_vec3 point,
    const double predicted[GAMERA_NO_FLUX_COUNT], double time,
    void *context, double rate[GAMERA_NO_FLUX_COUNT]) {
  (void)time;
  (void)context;
  if (predicted == NULL || rate == NULL ||
      !(predicted[GAMERA_NO_FLUX_DENSITY] > 0.0)) {
    return -1;
  }
  const double radius = sqrt(point.value[0] * point.value[0] +
                             point.value[1] * point.value[1] +
                             point.value[2] * point.value[2]);
  if (radius < obstacle_radius ||
      radius > obstacle_radius + 5.0 * ionosphere_scale_height) {
    return 0;
  }
  const double profile =
      exp(-(radius - obstacle_radius) / ionosphere_scale_height);
  /* The solar wind enters from -x, so the illuminated normal points -x. */
  const double solar_zenith_cosine = fmax(0.0, -point.value[0] / radius);
  const double illumination = 0.15 + 0.85 * solar_zenith_cosine;
  const double density = predicted[GAMERA_NO_FLUX_DENSITY];
  const double photo = photoionization_rate * profile * illumination;
  const double recombination =
      recombination_coefficient * profile * density * density;
  const double drag = neutral_drag_frequency * profile;
  rate[GAMERA_NO_FLUX_DENSITY] = photo - recombination;

  double momentum_squared = 0.0;
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    const int variable = GAMERA_NO_FLUX_MOMENTUM_X + component;
    const double momentum = predicted[variable];
    rate[variable] = -(recombination / density + drag) * momentum;
    momentum_squared += momentum * momentum;
  }
  const double energy = predicted[GAMERA_NO_FLUX_ENERGY];
  rate[GAMERA_NO_FLUX_ENERGY] =
      photo * injected_specific_internal_energy -
      recombination * energy / density - drag * momentum_squared / density;
  return isfinite(rate[GAMERA_NO_FLUX_ENERGY]) ? 0 : -1;
}
#endif

void problem_init(void) {
#ifdef GAMERA_EARTH_UPSTREAM_STARTUP
  log_info("Initializing upstream-launched Earth magnetosphere: stationary "
           "ambient (rho=%.3g,p=%.3g), wind/IMF front x=%.1f RE width %.1f "
           "RE, cut dipole rCut=%.1f lCut=%.1f, r=[%.1f,%.1f] RE, "
           "radial map v%d (legacy stretch %.8g), patch %d",
           startup_ambient_density, startup_ambient_pressure,
           startup_front_x, startup_front_width, earth_cut_radius,
           earth_cut_length, obstacle_radius, outer_radius,
           requested_radial_map_version(), radial_stretch, patch_id);
#elif defined(GAMERA_EARTH_MAGNETOSPHERE_STRETCHED)
  log_info("Initializing stretched-grid Earth magnetosphere: Mach %.1f "
           "wind IMF=(%.6g,%.6g,%.6g) nT plus 0.31 G split background dipole, "
           "r=[%.1f,%.1f] RE, radial map v%d (legacy stretch %.2f), "
           "Boris CA=%.1f on patch %d",
           inflow_mach,
           initial_magnetic.value[0] * norm_config.B_Norm * 1.0e9,
           initial_magnetic.value[1] * norm_config.B_Norm * 1.0e9,
           initial_magnetic.value[2] * norm_config.B_Norm * 1.0e9,
           obstacle_radius, outer_radius, requested_radial_map_version(),
           radial_stretch, CA, patch_id);
#elif defined(GAMERA_BOW_SHOCK)
  log_info("Initializing Mach %.1f unmagnetized Yin-Yang bow shock: "
           "r=[%.1f,%.1f], radial stretch %.2f, Boris CA=%.1f on patch %d",
           inflow_mach, obstacle_radius, outer_radius, radial_stretch, CA,
           patch_id);
#elif defined(GAMERA_EARTH_DIPOLE_BACKGROUND)
  log_info("Initializing Yin-Yang Earth magnetosphere: solar wind plus "
           "12-point background dipole on patch %d", patch_id);
#elif defined(GAMERA_IONOSPHERE_SOURCE_LOSS)
#ifdef GAMERA_CRUSTAL_FIELD
  log_info("Initializing Mach %.1f Yin-Yang ionosphere with localized "
           "crustal field on patch %d", inflow_mach, patch_id);
#else
  log_info("Initializing Mach %.1f Yin-Yang IMF obstacle with ionosphere "
           "source/loss on patch %d", inflow_mach, patch_id);
#endif
#else
  log_info("Initializing Mach %.1f Yin-Yang IMF sphere obstacle on patch %d",
           inflow_mach, patch_id);
#endif
  if (gamera_no_legacy_adapter_create_empty() != 0) {
    log_error("Failed to create Yin-Yang IMF-obstacle adapter");
    return;
  }
  gamera_no_grid *grid = gamera_no_legacy_grid();
  gamera_no_storage *storage = gamera_no_legacy_storage();
  if (gamera_no_initialize_primitives(
          grid, storage, inflow_primitive, NULL, false, gamma_val, rho_floor,
          p_floor) != 0 ||
      gamera_no_initialize_flux_from_vector_potential(
          grid, storage, conducting_sphere_potential, NULL) != 0 ||
      gamera_no_save_current_as_old(storage) != 0 ||
      gamera_no_legacy_export() != 0) {
    log_error("Failed to initialize Yin-Yang IMF sphere obstacle");
    gamera_no_legacy_adapter_destroy();
  }
}

#ifdef GAMERA_EARTH_DIPOLE_BACKGROUND
int problem_nonorthogonal_background_field(gamera_no_vec3 point,
                                            void *context,
                                            gamera_no_vec3 *field) {
  (void)context;
  if (gamera_no_dipole_field(point, &earth_dipole, field) != 0) {
    return -1;
  }
#ifdef GAMERA_EARTH_UPSTREAM_STARTUP
  const double radius = sqrt(point.value[0] * point.value[0] +
                             point.value[1] * point.value[1] +
                             point.value[2] * point.value[2]);
  const double weight =
      cubic_ramp_down(radius, earth_cut_radius, earth_cut_length);
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    field->value[component] *= weight;
  }
#endif
  return 0;
}

int problem_nonorthogonal_adjust_background(
    const gamera_no_grid *grid, gamera_no_background_data *background) {
  if (grid == NULL || background == NULL || background->cell_force == NULL) {
    return -1;
  }
  /* Zero dpB0 in the force-free inner part of the cut dipole. */
  const size_t count = gamera_no_element_count3(grid->cell_extent);
  for (size_t cell = 0; cell < count; ++cell) {
#ifdef GAMERA_EARTH_UPSTREAM_STARTUP
    const gamera_no_vec3 point = grid->cell[cell].centroid;
    const double radius = sqrt(point.value[0] * point.value[0] +
                               point.value[1] * point.value[1] +
                               point.value[2] * point.value[2]);
    if (radius > 0.75 * earth_cut_radius) {
      continue;
    }
#endif
    background->cell_force[cell] =
        (gamera_no_vec3){{0.0, 0.0, 0.0}};
  }
  return 0;
}

#ifdef GAMERA_NONORTHOGONAL_HAS_BACKGROUND_RESIDUAL_INIT
int problem_nonorthogonal_initialize_background_residual(
    const gamera_no_grid *grid, gamera_no_storage *storage,
    const gamera_no_background_data *background) {
  if (grid == NULL || storage == NULL || background == NULL) {
    return -1;
  }
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    if (storage->face_flux[direction] == NULL ||
        storage->old_face_flux[direction] == NULL ||
        background->face_flux[direction] == NULL) {
      return -1;
    }
    const size_t face_count =
        gamera_no_element_count3(grid->face[direction].extent);
#pragma omp parallel for schedule(static)
    for (size_t face = 0; face < face_count; ++face) {
      storage->face_flux[direction][face] -=
          background->face_flux[direction][face];
      storage->old_face_flux[direction][face] -=
          background->face_flux[direction][face];
    }
  }
  const double *current_flux[GAMERA_NO_DIM] = {
      storage->face_flux[GAMERA_NO_I], storage->face_flux[GAMERA_NO_J],
      storage->face_flux[GAMERA_NO_K]};
  const double *old_flux[GAMERA_NO_DIM] = {
      storage->old_face_flux[GAMERA_NO_I],
      storage->old_face_flux[GAMERA_NO_J],
      storage->old_face_flux[GAMERA_NO_K]};
  if (gamera_no_recover_magnetic_field(grid, current_flux,
                                       storage->cell_magnetic) != 0 ||
      gamera_no_recover_magnetic_field(grid, old_flux,
                                       storage->old_cell_magnetic) != 0) {
    return -1;
  }
  log_info("Converted initial total CT flux to residual B1=Btotal-B0 on "
           "patch %d",
           patch_id);
  return 0;
}
#endif
#endif

#endif
