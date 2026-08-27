#ifndef GAMERA_NONORTHOGONAL_YINYANG_H
#define GAMERA_NONORTHOGONAL_YINYANG_H

#include "nonorthogonal_geometry.h"

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

enum { GAMERA_NO_YIN_PATCH = 0, GAMERA_NO_YANG_PATCH = 1 };

/*
 * Both patches store vectors in the same global Cartesian xyz basis.  These
 * rotations are used only to map physical points to/from a patch's logical
 * spherical coordinates.
 *
 * Yang local -> global: (x,y,z) = (-xl,zl,yl).  The rotation is self-inverse.
 */
int gamera_no_yinyang_local_to_global(int patch, gamera_no_vec3 local,
                                      gamera_no_vec3 *global);
int gamera_no_yinyang_global_to_local(int patch, gamera_no_vec3 global,
                                      gamera_no_vec3 *local);

int gamera_no_yinyang_logical_to_global(int patch, double radius,
                                        double theta, double phi,
                                        gamera_no_vec3 *global);
int gamera_no_yinyang_global_to_logical(int patch, gamera_no_vec3 global,
                                        double *radius, double *theta,
                                        double *phi);

typedef struct gamera_no_yinyang_angular_domain {
  double theta_min;
  double theta_max;
  double phi_min;
  double phi_max;
  size_t theta_cells;
  size_t phi_cells;
  /* Straight Cartesian cell faces may bow slightly past analytic bounds. */
  double boundary_tolerance_cells;
} gamera_no_yinyang_angular_domain;

/*
 * Return the signed distance, in logical cells, from point to the nearest
 * angular edge of patch.  valid accepts the configured geometric tolerance;
 * margin itself is not clipped and remains useful for deterministic ownership.
 */
int gamera_no_yinyang_angular_margin(
    int patch, gamera_no_vec3 point,
    const gamera_no_yinyang_angular_domain *domain, int *valid,
    double *margin_cells);

/*
 * Hard composite owner used for physical integrals and ledgers: choose the
 * only valid patch, or the patch farther from its angular fringe.  Exact ties
 * go to Yin.  MFE may still evolve/apply sources to both active copies; this
 * owner exists to count the overlapping physical volume exactly once.
 */
int gamera_no_yinyang_composite_owner(
    gamera_no_vec3 point, const gamera_no_yinyang_angular_domain *domain,
    int *owner_patch, int *in_overlap, double margin_cells[2]);

#ifdef __cplusplus
}
#endif

#endif
