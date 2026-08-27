#ifndef GAMERA_NONORTHOGONAL_RECONSTRUCTION_H
#define GAMERA_NONORTHOGONAL_RECONSTRUCTION_H

#ifdef __cplusplus
extern "C" {
#endif

enum { GAMERA_NO_RECON_STENCIL = 8 };

double gamera_no_central8(const double stencil[GAMERA_NO_RECON_STENCIL]);
double gamera_no_central6(const double stencil[GAMERA_NO_RECON_STENCIL]);

/*
 * Fortran recon.F90:Up7LRs for one primitive variable and one interface.
 * volume and primitive run from the four cells below the interface through
 * the four cells above it. pdmb_code=1 corresponds to A=4 in the paper.
 */
int gamera_no_reconstruct_up7_pdm(
    const double volume[GAMERA_NO_RECON_STENCIL],
    const double primitive[GAMERA_NO_RECON_STENCIL], double pdmb_code,
    double *left, double *right);

#ifdef __cplusplus
}
#endif

#endif
