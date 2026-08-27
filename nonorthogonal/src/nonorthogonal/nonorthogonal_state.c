#include "nonorthogonal_state.h"

#include <float.h>
#include <math.h>
#include <stddef.h>

static double vector_dot(gamera_no_vec3 a, gamera_no_vec3 b) {
  double result = 0.0;
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    result += a.value[d] * b.value[d];
  }
  return result;
}

static gamera_no_vec3 vector_add(gamera_no_vec3 a, gamera_no_vec3 b) {
  gamera_no_vec3 result;
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    result.value[d] = a.value[d] + b.value[d];
  }
  return result;
}

static gamera_no_vec3 vector_scale(gamera_no_vec3 value, double scale) {
  gamera_no_vec3 result;
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    result.value[d] = scale * value.value[d];
  }
  return result;
}

static gamera_no_vec3 vector_parallel(gamera_no_vec3 value,
                                      gamera_no_vec3 unit) {
  return vector_scale(unit, vector_dot(value, unit));
}

static gamera_no_vec3 vector_perpendicular(gamera_no_vec3 value,
                                           gamera_no_vec3 unit) {
  gamera_no_vec3 parallel = vector_parallel(value, unit);
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    value.value[d] -= parallel.value[d];
  }
  return value;
}

static gamera_no_vec3 normalized_or_zero(gamera_no_vec3 value) {
  const double norm = sqrt(vector_dot(value, value));
  if (norm <= 1.0e-12) {
    return (gamera_no_vec3){{0.0, 0.0, 0.0}};
  }
  return vector_scale(value, 1.0 / norm);
}

static int valid_conversion_parameters(double gamma, double density_floor,
                                       double pressure_floor) {
  return isfinite(gamma) && gamma > 1.0 && isfinite(density_floor) &&
         density_floor > 0.0 && isfinite(pressure_floor) &&
         pressure_floor > 0.0;
}

int gamera_no_conserved_to_primitive(
    const double conserved[GAMERA_NO_FLUX_COUNT], double gamma,
    double density_floor, double pressure_floor,
    gamera_no_primitive *primitive) {
  if (conserved == NULL || primitive == NULL ||
      !valid_conversion_parameters(gamma, density_floor, pressure_floor) ||
      !isfinite(conserved[GAMERA_NO_FLUX_DENSITY]) ||
      conserved[GAMERA_NO_FLUX_DENSITY] <= DBL_MIN) {
    return -1;
  }

  const double raw_density = conserved[GAMERA_NO_FLUX_DENSITY];
  double velocity_squared = 0.0;
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    const double momentum =
        conserved[GAMERA_NO_FLUX_MOMENTUM_X + d];
    if (!isfinite(momentum)) {
      return -1;
    }
    primitive->velocity.value[d] = momentum / raw_density;
    velocity_squared += primitive->velocity.value[d] *
                        primitive->velocity.value[d];
  }
  if (!isfinite(conserved[GAMERA_NO_FLUX_ENERGY])) {
    return -1;
  }
  const double kinetic = 0.5 * raw_density * velocity_squared;
  primitive->density = fmax(raw_density, density_floor);
  primitive->pressure =
      fmax((gamma - 1.0) *
               (conserved[GAMERA_NO_FLUX_ENERGY] - kinetic),
           pressure_floor);
  return 0;
}

int gamera_no_primitive_to_conserved(
    const gamera_no_primitive *primitive, double gamma, double density_floor,
    double pressure_floor, double conserved[GAMERA_NO_FLUX_COUNT]) {
  if (primitive == NULL || conserved == NULL ||
      !valid_conversion_parameters(gamma, density_floor, pressure_floor) ||
      !isfinite(primitive->density) || !isfinite(primitive->pressure)) {
    return -1;
  }

  const double density = fmax(primitive->density, density_floor);
  const double pressure = fmax(primitive->pressure, pressure_floor);
  double velocity_squared = 0.0;
  conserved[GAMERA_NO_FLUX_DENSITY] = density;
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    if (!isfinite(primitive->velocity.value[d])) {
      return -1;
    }
    conserved[GAMERA_NO_FLUX_MOMENTUM_X + d] =
        density * primitive->velocity.value[d];
    velocity_squared += primitive->velocity.value[d] *
                        primitive->velocity.value[d];
  }
  conserved[GAMERA_NO_FLUX_ENERGY] =
      pressure / (gamma - 1.0) + 0.5 * density * velocity_squared;
  return 0;
}

int gamera_no_predict_cell(
    const double old_conserved[GAMERA_NO_FLUX_COUNT],
    const double current_conserved[GAMERA_NO_FLUX_COUNT], double ratio,
    double gamma, double density_floor, double pressure_floor,
    double predicted_conserved[GAMERA_NO_FLUX_COUNT]) {
  if (old_conserved == NULL || current_conserved == NULL ||
      predicted_conserved == NULL || !isfinite(ratio)) {
    return -1;
  }
  gamera_no_primitive old_primitive;
  gamera_no_primitive current_primitive;
  if (gamera_no_conserved_to_primitive(old_conserved, gamma, density_floor,
                                       pressure_floor, &old_primitive) != 0 ||
      gamera_no_conserved_to_primitive(current_conserved, gamma, density_floor,
                                       pressure_floor,
                                       &current_primitive) != 0) {
    return -1;
  }

  gamera_no_primitive predicted = current_primitive;
  predicted.density = fmax(current_primitive.density +
                               ratio * (current_primitive.density -
                                        old_primitive.density),
                           density_floor);
  predicted.pressure = fmax(current_primitive.pressure +
                                ratio * (current_primitive.pressure -
                                         old_primitive.pressure),
                            pressure_floor);
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    predicted.velocity.value[d] =
        current_primitive.velocity.value[d] +
        ratio * (current_primitive.velocity.value[d] -
                 old_primitive.velocity.value[d]);
  }
  return gamera_no_primitive_to_conserved(&predicted, gamma, density_floor,
                                          pressure_floor,
                                          predicted_conserved);
}

void gamera_no_apply_reynolds(
    double conserved[GAMERA_NO_FLUX_COUNT],
    const double hydro_rate[GAMERA_NO_FLUX_COUNT], double dt) {
  for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
    conserved[variable] += dt * hydro_rate[variable];
  }
}

int gamera_no_apply_maxwell(double conserved[GAMERA_NO_FLUX_COUNT],
                            gamera_no_vec3 magnetic_momentum_rate, double dt,
                            double gamma, double density_floor,
                            double pressure_floor) {
  if (conserved == NULL || !isfinite(dt) ||
      !valid_conversion_parameters(gamma, density_floor, pressure_floor)) {
    return -1;
  }
  const double density =
      fmax(conserved[GAMERA_NO_FLUX_DENSITY], density_floor);
  double old_momentum_squared = 0.0;
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    const double momentum =
        conserved[GAMERA_NO_FLUX_MOMENTUM_X + d];
    old_momentum_squared += momentum * momentum;
  }
  const double pressure =
      fmax((gamma - 1.0) *
               (conserved[GAMERA_NO_FLUX_ENERGY] -
                0.5 * old_momentum_squared / density),
           pressure_floor);

  double new_momentum_squared = 0.0;
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    conserved[GAMERA_NO_FLUX_MOMENTUM_X + d] +=
        dt * magnetic_momentum_rate.value[d];
    const double momentum =
        conserved[GAMERA_NO_FLUX_MOMENTUM_X + d];
    new_momentum_squared += momentum * momentum;
  }
  conserved[GAMERA_NO_FLUX_DENSITY] = density;
  conserved[GAMERA_NO_FLUX_ENERGY] =
      0.5 * new_momentum_squared / density + pressure / (gamma - 1.0);
  return 0;
}

int gamera_no_apply_boris(
    double hydro_updated[GAMERA_NO_FLUX_COUNT],
    const double old_conserved[GAMERA_NO_FLUX_COUNT],
    gamera_no_vec3 new_total_magnetic, gamera_no_vec3 old_total_magnetic,
    const double hydro_rate[GAMERA_NO_FLUX_COUNT],
    gamera_no_vec3 magnetic_momentum_rate, double light_speed, double dt,
    double gamma, double density_floor, double pressure_floor) {
  if (hydro_updated == NULL || old_conserved == NULL || hydro_rate == NULL ||
      !isfinite(light_speed) || light_speed <= DBL_MIN || !isfinite(dt)) {
    return -1;
  }

  gamera_no_primitive updated_primitive;
  if (gamera_no_conserved_to_primitive(
          hydro_updated, gamma, density_floor, pressure_floor,
          &updated_primitive) != 0 ||
      old_conserved[GAMERA_NO_FLUX_DENSITY] <= DBL_MIN) {
    return -1;
  }
  const double updated_density = fmax(updated_primitive.density, density_floor);
  const double updated_pressure =
      fmax(updated_primitive.pressure, pressure_floor);
  const double old_density = old_conserved[GAMERA_NO_FLUX_DENSITY];
  gamera_no_vec3 old_velocity;
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    old_velocity.value[d] =
        old_conserved[GAMERA_NO_FLUX_MOMENTUM_X + d] /
        fmax(old_density, density_floor);
  }

  const gamera_no_vec3 half_magnetic =
      vector_scale(vector_add(new_total_magnetic, old_total_magnetic), 0.5);
  const double alpha = vector_dot(half_magnetic, half_magnetic) /
                       (updated_density * light_speed * light_speed);
  const double density_increment =
      dt * hydro_rate[GAMERA_NO_FLUX_DENSITY];
  gamera_no_vec3 gas_momentum_increment;
  gamera_no_vec3 magnetic_momentum_increment;
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    gas_momentum_increment.value[d] =
        dt * hydro_rate[GAMERA_NO_FLUX_MOMENTUM_X + d];
    magnetic_momentum_increment.value[d] =
        dt * magnetic_momentum_rate.value[d];
  }

  gamera_no_vec3 momentum_increment;
  if (alpha > 1.0e-12) {
    const double plasma_fraction = 1.0 / (1.0 + alpha);
    const gamera_no_vec3 field_unit = normalized_or_zero(half_magnetic);
    const gamera_no_vec3 total_stress_increment =
        vector_add(gas_momentum_increment, magnetic_momentum_increment);
    momentum_increment = vector_add(
        vector_parallel(gas_momentum_increment, field_unit),
        vector_add(
            vector_scale(vector_perpendicular(total_stress_increment,
                                              field_unit),
                         plasma_fraction),
            vector_scale(vector_perpendicular(old_velocity, field_unit),
                         alpha * plasma_fraction * density_increment)));
  } else {
    momentum_increment =
        vector_add(gas_momentum_increment, magnetic_momentum_increment);
  }

  gamera_no_primitive result = updated_primitive;
  result.density = updated_density;
  result.pressure = updated_pressure;
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    result.velocity.value[d] =
        (old_density * old_velocity.value[d] +
         momentum_increment.value[d]) /
        updated_density;
  }
  return gamera_no_primitive_to_conserved(&result, gamma, density_floor,
                                          pressure_floor, hydro_updated);
}
