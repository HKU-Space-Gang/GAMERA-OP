#ifndef GAMERA_NONORTHOGONAL_OPERATORS_H
#define GAMERA_NONORTHOGONAL_OPERATORS_H

#include "nonorthogonal_geometry.h"

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* GAMERA face-flux-to-cell-Cartesian-field recovery, including Div correction. */
int gamera_no_flux_to_cell_field(const gamera_no_cell_geometry *geometry,
                                 const double face_flux[GAMERA_NO_DIM][2],
                                 gamera_no_vec3 *field,
                                 double *net_flux);

/* Solve the non-orthogonal 2x2 edge-normal system used by GetCornerB. */
int gamera_no_solve_edge_field(const double normal1[2],
                               const double normal2[2],
                               double face_field1, double face_field2,
                               double edge_field[2]);

/*
 * Fortran init.F90:ebGeom: area-weight a face-normal stencil with the current
 * 6th-order centered edge interpolation, project it into (T1,T2), normalize,
 * and return one row of the edge 2x2 magnetic transform.
 */
int gamera_no_interpolate_face_normal_to_edge(
    const double face_area[8], const gamera_no_vec3 face_normal[8],
    const gamera_no_edge_geometry *edge, double transverse_normal[2]);

/*
 * Return CT increments for the six faces of one logical cell.
 * emf_i[j_side][k_side], emf_j[i_side][k_side], and
 * emf_k[i_side][j_side] are already line-integrated edge EMFs.
 */
void gamera_no_ct_face_increments(const double emf_i[2][2],
                                  const double emf_j[2][2],
                                  const double emf_k[2][2], double dt,
                                  double face_increment[GAMERA_NO_DIM][2]);

/* Fortran fields.F90 local edge EMF after interpolation/reconstruction. */
int gamera_no_compute_edge_emf(
    double velocity_tangent1, double velocity_tangent2,
    double magnetic_tangent1, double magnetic_tangent2,
    double magnetic_current_jump, double edge_density,
    const gamera_no_edge_geometry *edge, double diffusion_coefficient,
    bool use_boris, double light_speed, double cfl, double dt,
    double *line_integrated_emf, double *diffusion_speed);

#ifdef __cplusplus
}
#endif

#endif
