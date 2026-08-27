#include "nonorthogonal_flux.h"

#include <float.h>
#include <math.h>
#include <stddef.h>

static const double gamera_no_pi = 3.14159265358979323846264338327950288;

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

static int valid_primitive(gamera_no_primitive state) {
  if (!isfinite(state.density) || !isfinite(state.pressure) ||
      state.density <= DBL_MIN || state.pressure <= DBL_MIN) {
    return 0;
  }
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    if (!isfinite(state.velocity.value[d])) {
      return 0;
    }
  }
  return 1;
}

int gamera_no_kinetic_fluid_flux(const gamera_no_primitive state[2],
                                 double gamma,
                                 const gamera_no_face_geometry *face,
                                 gamera_no_fluid_flux *flux) {
  if (state == NULL || face == NULL || flux == NULL || !isfinite(gamma) ||
      gamma <= 1.0 || !valid_primitive(state[0]) ||
      !valid_primitive(state[1])) {
    return -1;
  }

  for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
    flux->conserved[variable] = 0.0;
    flux->conserved_jump[variable] = 0.0;
  }
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    flux->velocity_jump.value[d] = 0.0;
  }

  const double side_sign[2] = {-1.0, 1.0};
  for (int side = 0; side < 2; ++side) {
    const double density = state[side].density;
    const double pressure = state[side].pressure;
    const gamera_no_vec3 velocity = state[side].velocity;
    const double speed_squared = vector_dot(velocity, velocity);
    const double energy =
        0.5 * density * speed_squared + pressure / (gamma - 1.0);
    const double lambda = density / (2.0 * pressure);
    const double normal_velocity = vector_dot(face->normal, velocity);
    const double tangent1_velocity = vector_dot(face->tangent1, velocity);
    const double tangent2_velocity = vector_dot(face->tangent2, velocity);
    const double sign = side_sign[side];
    const double zeroth_moment =
        0.5 * erfc(sign * sqrt(lambda) * normal_velocity);
    const double first_moment =
        normal_velocity * zeroth_moment -
        0.5 * sign * exp(-lambda * normal_velocity * normal_velocity) /
            sqrt(gamera_no_pi * lambda);

    flux->normal_velocity[side] = normal_velocity;
    flux->conserved[GAMERA_NO_FLUX_DENSITY] += density * first_moment;
    flux->conserved[GAMERA_NO_FLUX_ENERGY] +=
        (energy + 0.5 * pressure) * first_moment +
        0.5 * pressure * zeroth_moment * normal_velocity;

    const double normal_momentum =
        density * normal_velocity * first_moment + pressure * zeroth_moment;
    const double tangent1_momentum =
        density * tangent1_velocity * first_moment;
    const double tangent2_momentum =
        density * tangent2_velocity * first_moment;
    for (int d = 0; d < GAMERA_NO_DIM; ++d) {
      flux->conserved[GAMERA_NO_FLUX_MOMENTUM_X + d] +=
          face->normal.value[d] * normal_momentum +
          face->tangent1.value[d] * tangent1_momentum +
          face->tangent2.value[d] * tangent2_momentum;
      flux->velocity_jump.value[d] += sign * velocity.value[d];
    }

    flux->conserved_jump[GAMERA_NO_FLUX_DENSITY] += sign * density;
    for (int d = 0; d < GAMERA_NO_DIM; ++d) {
      flux->conserved_jump[GAMERA_NO_FLUX_MOMENTUM_X + d] +=
          sign * density * velocity.value[d];
    }
    flux->conserved_jump[GAMERA_NO_FLUX_ENERGY] += sign * energy;
  }

  return 0;
}

int gamera_no_kinetic_maxwell_flux(
    const gamera_no_primitive state[2], const gamera_no_vec3 magnetic[2],
    double face_normal_field, const gamera_no_face_geometry *face,
    bool use_background, gamera_no_vec3 background,
    double background_face_normal_field, bool use_boris, double light_speed,
    const double normal_velocity[2], gamera_no_maxwell_flux *flux) {
  if (state == NULL || magnetic == NULL || face == NULL ||
      normal_velocity == NULL || flux == NULL || !isfinite(face_normal_field) ||
      (use_boris && (!isfinite(light_speed) || light_speed <= DBL_MIN))) {
    return -1;
  }
  for (int side = 0; side < 2; ++side) {
    if (!valid_primitive(state[side])) {
      return -1;
    }
    for (int d = 0; d < GAMERA_NO_DIM; ++d) {
      if (!isfinite(magnetic[side].value[d]) ||
          (use_background && !isfinite(background.value[d]))) {
        return -1;
      }
    }
  }

  flux->momentum = (gamera_no_vec3){{0.0, 0.0, 0.0}};
  flux->alfven_diffusion_speed = 0.0;
  flux->magnetic_pressure_sum = 0.0;

  const double side_sign[2] = {-1.0, 1.0};
  for (int side = 0; side < 2; ++side) {
    const double residual_squared =
        vector_dot(magnetic[side], magnetic[side]);
    const double residual_pressure = 0.5 * residual_squared;
    const gamera_no_vec3 total =
        use_background ? vector_add(magnetic[side], background) : magnetic[side];
    const double total_squared = vector_dot(total, total);
    const double alfven_squared = total_squared / state[side].density;
    const double corrected_alfven_squared =
        use_boris
            ? alfven_squared /
                  (1.0 + alfven_squared / (light_speed * light_speed))
            : alfven_squared;
    const double lambda =
        1.0 / (2.0 * state[side].pressure / state[side].density +
               corrected_alfven_squared);
    const double zeroth_moment =
        0.5 * erfc(side_sign[side] * sqrt(lambda) * normal_velocity[side]);

    flux->alfven_diffusion_speed +=
        0.5 * sqrt(corrected_alfven_squared);
    flux->magnetic_pressure_sum += 0.5 * total_squared;

    const double background_dot_residual =
        use_background ? vector_dot(background, magnetic[side]) : 0.0;
    for (int d = 0; d < GAMERA_NO_DIM; ++d) {
      double stress = residual_pressure * face->normal.value[d] -
                      magnetic[side].value[d] * face_normal_field;
      if (use_background) {
        stress += -background.value[d] * face_normal_field -
                  magnetic[side].value[d] *
                      background_face_normal_field +
                  background_dot_residual * face->normal.value[d];
      }
      flux->momentum.value[d] += zeroth_moment * stress;
    }
  }
  return 0;
}

int gamera_no_apply_hogs(gamera_no_fluid_flux *fluid,
                         gamera_no_maxwell_flux *maxwell,
                         double hydro_coefficient,
                         double magnetic_coefficient, bool use_boris,
                         double light_speed) {
  if (fluid == NULL || maxwell == NULL || !isfinite(hydro_coefficient) ||
      !isfinite(magnetic_coefficient) || hydro_coefficient < 0.0 ||
      magnetic_coefficient < 0.0 ||
      (use_boris && (!isfinite(light_speed) || light_speed <= DBL_MIN))) {
    return -1;
  }

  for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
    fluid->conserved[variable] -=
        hydro_coefficient * maxwell->alfven_diffusion_speed *
        fluid->conserved_jump[variable];
  }
  if (use_boris) {
    const double inverse_light_speed_squared = 1.0 / (light_speed * light_speed);
    for (int d = 0; d < GAMERA_NO_DIM; ++d) {
      const double magnetic_momentum_jump =
          maxwell->magnetic_pressure_sum * fluid->velocity_jump.value[d] *
          inverse_light_speed_squared;
      maxwell->momentum.value[d] -=
          magnetic_coefficient * maxwell->alfven_diffusion_speed *
          magnetic_momentum_jump;
    }
  }
  return 0;
}

int gamera_no_apply_emergency_interface_diffusion(
    gamera_no_fluid_flux *fluid, const gamera_no_primitive interface[2],
    const double lower_cell[GAMERA_NO_FLUX_COUNT],
    const double upper_cell[GAMERA_NO_FLUX_COUNT], bool use_boris,
    double light_speed, bool *applied) {
  if (fluid == NULL || interface == NULL || lower_cell == NULL ||
      upper_cell == NULL || applied == NULL ||
      (use_boris && (!isfinite(light_speed) || light_speed <= DBL_MIN))) {
    return -1;
  }
  *applied = false;
  if (!use_boris) {
    return 0;
  }
  double speed_square_sum = 0.0;
  for (int side = 0; side < 2; ++side) {
    if (!valid_primitive(interface[side])) {
      return -1;
    }
    speed_square_sum +=
        vector_dot(interface[side].velocity, interface[side].velocity);
  }
  const double interface_speed = sqrt(0.5 * speed_square_sum);
  if (interface_speed < 1.5 * light_speed) {
    return 0;
  }
  for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
    if (!isfinite(lower_cell[variable]) || !isfinite(upper_cell[variable])) {
      return -1;
    }
    fluid->conserved[variable] -=
        light_speed * (upper_cell[variable] - lower_cell[variable]);
  }
  *applied = true;
  return 0;
}
