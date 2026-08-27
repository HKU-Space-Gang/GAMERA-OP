#ifndef GAMERA_NONORTHOGONAL_STATE_H
#define GAMERA_NONORTHOGONAL_STATE_H

#include "nonorthogonal_flux.h"

#ifdef __cplusplus
extern "C" {
#endif

int gamera_no_conserved_to_primitive(
    const double conserved[GAMERA_NO_FLUX_COUNT], double gamma,
    double density_floor, double pressure_floor,
    gamera_no_primitive *primitive);

int gamera_no_primitive_to_conserved(
    const gamera_no_primitive *primitive, double gamma, double density_floor,
    double pressure_floor, double conserved[GAMERA_NO_FLUX_COUNT]);

/* Fortran mhdgroup.F90:CellPredictor for a single-fluid cell. */
int gamera_no_predict_cell(
    const double old_conserved[GAMERA_NO_FLUX_COUNT],
    const double current_conserved[GAMERA_NO_FLUX_COUNT], double ratio,
    double gamma, double density_floor, double pressure_floor,
    double predicted_conserved[GAMERA_NO_FLUX_COUNT]);

void gamera_no_apply_reynolds(
    double conserved[GAMERA_NO_FLUX_COUNT],
    const double hydro_rate[GAMERA_NO_FLUX_COUNT], double dt);

/* Fortran mhdgroup.F90:CellMaxwell for a single-fluid cell. */
int gamera_no_apply_maxwell(double conserved[GAMERA_NO_FLUX_COUNT],
                            gamera_no_vec3 magnetic_momentum_rate, double dt,
                            double gamma, double density_floor,
                            double pressure_floor);

/* Fortran mhdgroup.F90:CellBoris, single-fluid branch. */
int gamera_no_apply_boris(
    double hydro_updated[GAMERA_NO_FLUX_COUNT],
    const double old_conserved[GAMERA_NO_FLUX_COUNT],
    gamera_no_vec3 new_total_magnetic, gamera_no_vec3 old_total_magnetic,
    const double hydro_rate[GAMERA_NO_FLUX_COUNT],
    gamera_no_vec3 magnetic_momentum_rate, double light_speed, double dt,
    double gamma, double density_floor, double pressure_floor);

#ifdef __cplusplus
}
#endif

#endif
