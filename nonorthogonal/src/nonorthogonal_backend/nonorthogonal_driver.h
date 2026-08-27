#ifndef GAMERA_NONORTHOGONAL_DRIVER_H
#define GAMERA_NONORTHOGONAL_DRIVER_H

#ifdef GAMERA_NONORTHOGONAL_HAS_BACKGROUND_FIELD
#include "nonorthogonal_sweep.h"
#endif

/* Prepare/free problem-specific numerical data after the adapter exists. */
int gamera_no_driver_prepare(void);
void gamera_no_driver_finalize(void);

#ifdef GAMERA_NONORTHOGONAL_HAS_BACKGROUND_FIELD
const gamera_no_background_field *gamera_no_driver_background_field(void);
#endif

int solve_nonorthogonal(int maximum_steps);

#endif
