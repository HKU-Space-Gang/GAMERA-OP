#ifndef GAMERA_NONORTHOGONAL_LEGACY_ADAPTER_H
#define GAMERA_NONORTHOGONAL_LEGACY_ADAPTER_H

#include "nonorthogonal_grid.h"
#include "nonorthogonal_storage.h"

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Adapter for GAMERA-OP's global gas/gem/x1/x2/x3 arrays.  Under the future
 * non-orthogonal coordinate type, x1/x2/x3 are Cartesian vertex coordinates,
 * gas velocities/momenta and mags_b1..b3 are Cartesian components, while
 * mag_bi/bj/bk remain face-normal values at the legacy boundary/I/O surface.
 */
int gamera_no_legacy_adapter_create(bool import_restart_history);
/* Build geometry/storage without importing legacy state (fresh IC path). */
int gamera_no_legacy_adapter_create_empty(void);
void gamera_no_legacy_adapter_destroy(void);

/* Import current primitive/face-normal state after MPI exchange and BC. */
int gamera_no_legacy_import_current(void);

/* Export current and old state, EMF, geometry diagnostics and divB. */
int gamera_no_legacy_export(void);

gamera_no_grid *gamera_no_legacy_grid(void);
gamera_no_storage *gamera_no_legacy_storage(void);

#ifdef __cplusplus
}
#endif

#endif
