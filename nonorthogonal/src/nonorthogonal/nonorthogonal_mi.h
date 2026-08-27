#ifndef GAMERA_NONORTHOGONAL_MI_H
#define GAMERA_NONORTHOGONAL_MI_H

#include "nonorthogonal_grid.h"

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

enum {
  GAMERA_MI_NORTH = -1,
  GAMERA_MI_SOUTH = 1
};

/*
 * Kaiju GetShellJ finite-volume/Stokes current from a Cartesian residual B.
 * The target cell needs one magnetic-cell halo in every direction.
 */
int gamera_no_cell_current_from_residual(
    const gamera_no_grid *grid, const gamera_no_vec3 *residual_magnetic,
    size_t i, size_t j, size_t k, gamera_no_vec3 *current);

/* Dipole J/B mapping helpers. Colatitude is measured from either pole. */
double gamera_mi_mapped_colatitude(double mhd_colatitude,
                                   double mhd_radius,
                                   double ionosphere_radius);
double gamera_mi_dipole_field_ratio(double mhd_colatitude,
                                    double mhd_radius,
                                    double ionosphere_radius);
double gamera_mi_dipole_cos_inclination(double ionosphere_colatitude,
                                        int hemisphere);

/*
 * Compose the MHD inner-ghost velocity from the ionospheric E x B drift.
 * GAMERA-OP reconstructs the prescribed wall drift from a conjugate active
 * state, so the ghost value is 2*V_EB-v_active.
 */
int gamera_mi_compose_ghost_velocity(gamera_no_vec3 source_velocity,
                                     gamera_no_vec3 background_magnetic,
                                     gamera_no_vec3 drift_velocity,
                                     gamera_no_vec3 *ghost_velocity);

typedef struct {
  size_t longitude_count;
  size_t colatitude_count; /* includes pole and low-latitude boundary */
  double maximum_colatitude;
  double ionosphere_radius_m;
  double pedersen_siemens;
  double hall_siemens; /* first implementation requires exactly zero */
  double low_latitude_potential_v;
  int hemisphere;
  int maximum_iterations;
  double relative_tolerance;
  double absolute_tolerance;
} gamera_mi_solver_config;

typedef struct {
  int iterations;
  double initial_residual;
  double final_residual;
  int converged;
} gamera_mi_solver_stats;

/*
 * Apply the positive-sign discrete constant-Pedersen operator. The pole row
 * imposes Psi_pole(phi)=mean(Psi(first ring)); the outer row is Dirichlet.
 */
int gamera_mi_apply_constant_pedersen(
    const gamera_mi_solver_config *config, const double *potential_v,
    double *result);

/*
 * FAC is in A/m^2, potential is in volts, and arrays are theta-major with
 * longitude contiguous. potential_v is also the iterative initial guess;
 * non-finite entries are replaced with zero.  The constant-Pedersen solve
 * uses one physical pole degree of freedom and area-weighted PCG, then
 * expands the pole value back to every longitude slot.
 */
int gamera_mi_solve_constant_pedersen(
    const gamera_mi_solver_config *config, const double *fac_a_m2,
    double *potential_v, gamera_mi_solver_stats *stats);

/*
 * Apply/solve the spatially varying height-integrated conductance tensor.
 * SigmaP and SigmaH use the same theta-major node layout as potential.  The
 * Hall tensor is antisymmetric, so its mixed second derivatives cancel and
 * only conductance-gradient first-derivative terms remain.  The resulting
 * operator is nonsymmetric and is solved with diagonally preconditioned
 * BiCGStab.  Constant SigmaP with SigmaH=0 is exactly the legacy
 * constant-Pedersen operator; nonzero Hall conductance retains the dipole
 * inclination metric used by REMIX.
 */
int gamera_mi_apply_conductance_tensor(
    const gamera_mi_solver_config *config, const double *pedersen_siemens,
    const double *hall_siemens, const double *potential_v, double *result);

int gamera_mi_solve_conductance_tensor(
    const gamera_mi_solver_config *config, const double *fac_a_m2,
    const double *pedersen_siemens, const double *hall_siemens,
    double *potential_v, gamera_mi_solver_stats *stats);

#ifdef __cplusplus
}
#endif

#endif
