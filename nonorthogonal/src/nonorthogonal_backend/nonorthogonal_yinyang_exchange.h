#ifndef GAMERA_NONORTHOGONAL_YINYANG_EXCHANGE_H
#define GAMERA_NONORTHOGONAL_YINYANG_EXCHANGE_H

#include "nonorthogonal_grid.h"
#include "nonorthogonal_storage.h"

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Exchange strict-HD Cartesian conserved state across Yin/Yang boundaries. */
int gamera_no_yinyang_exchange_hd(void);
/* Refresh only angular HD ghosts, preserving checkpointed active overlap. */
int gamera_no_yinyang_exchange_hd_ghosts_only(void);
/* Synchronize receptor-edge Cartesian electric-field line integrals for CT. */
int gamera_no_yinyang_sync_edge_emf(
    gamera_no_storage *storage, const gamera_no_grid *grid,
    const size_t active_lower[3], const size_t active_upper[3],
    void *context);
/* Fill only angular ghost face fluxes from donor Cartesian magnetic fields. */
int gamera_no_yinyang_exchange_magnetic_ghosts(void);
/* Composite owner for physical integrals; independent of MFE state updates. */
int gamera_no_yinyang_physical_owner(gamera_no_vec3 point,
                                     int *owner_patch, int *in_overlap,
                                     double margin_cells[2]);
void gamera_no_yinyang_exchange_destroy(void);

extern int gamera_no_yinyang_receptor_count;
extern int gamera_no_yinyang_active_receptor_count;
extern double gamera_no_yinyang_max_donor_extrapolation;
extern int gamera_no_yinyang_edge_receptor_count;
extern double gamera_no_yinyang_max_emf_donor_extrapolation;
extern int gamera_no_yinyang_magnetic_face_receptor_count;
extern double gamera_no_yinyang_max_magnetic_donor_extrapolation;
extern int gamera_no_yinyang_active_magnetic_face_receptor_count;
extern double gamera_no_yinyang_hdiv_max_before;
extern double gamera_no_yinyang_hdiv_max_after;
extern double gamera_no_yinyang_hdiv_max_correction;
extern int gamera_no_yinyang_hdiv_max_iterations;
extern double gamera_no_yinyang_hdiv_min_weight;
extern double gamera_no_yinyang_hdiv_max_weight;
extern double gamera_no_yinyang_hdiv_min_metric_cosine;

#ifdef __cplusplus
}
#endif

#endif
