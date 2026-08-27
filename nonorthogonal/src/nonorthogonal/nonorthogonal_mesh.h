#ifndef GAMERA_NONORTHOGONAL_MESH_H
#define GAMERA_NONORTHOGONAL_MESH_H

#include "nonorthogonal_geometry.h"

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Fill a caller-owned (ni+1)*(nj+1)*(nk+1) Cartesian vertex array.  The
 * array uses gamera_no_index3 ordering (k contiguous).
 */
int gamera_no_generate_cartesian_vertices(
    const size_t cell_extent[3], const double lower[3], const double upper[3],
    gamera_no_vec3 *vertices);

typedef enum {
  /* Zhang et al. (2019), Equations (102)-(103): x and y warp together. */
  GAMERA_NO_WARP_PAPER = 1,
  /* kaiju/src/gamera/init.F90: x += dsp and y -= dsp. */
  GAMERA_NO_WARP_FORTRAN = -1
} gamera_no_warp_convention;

/*
 * Generate the distorted Cartesian mapping used by the GAMERA tests:
 *
 *   x = xmin + Lx [u + w0 sin(n*pi*u) sin(n*pi*v)]
 *   y = ymin + Ly [v + s*w0 sin(n*pi*u) sin(n*pi*v)]
 *   z = zmin + Lz*w
 *
 * where s=+1 is the paper mapping and s=-1 is the current Fortran mapping.
 * The amplitude w0 is a fraction of the corresponding domain length.
 */
int gamera_no_generate_warped_cartesian_vertices(
    const size_t cell_extent[3], const double lower[3], const double upper[3],
    double w0, int wave_number, gamera_no_warp_convention convention,
    gamera_no_vec3 *vertices);

/* Map logical (r,theta,phi) vertices into global Cartesian coordinates. */
int gamera_no_generate_spherical_vertices(
    const size_t cell_extent[3], const double lower[3], const double upper[3],
    gamera_no_vec3 *vertices);

#ifdef __cplusplus
}
#endif

#endif
