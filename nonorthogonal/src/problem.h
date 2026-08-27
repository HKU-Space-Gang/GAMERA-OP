#ifndef PROBLEM_H
#define PROBLEM_H

#ifdef GAMERA_NONORTHOGONAL_HAS_CELL_SOURCE
#include "nonorthogonal_flux.h"
#endif
#if defined(GAMERA_NONORTHOGONAL_HAS_FLUID_FLUX_BOUNDARY) || \
    defined(GAMERA_NONORTHOGONAL_HAS_EDGE_EMF_BOUNDARY) || \
    defined(GAMERA_NONORTHOGONAL_HAS_MAGNETIC_BOUNDARY)
#include "nonorthogonal_grid.h"
#include "nonorthogonal_storage.h"
#endif
#ifdef GAMERA_NONORTHOGONAL_HAS_BACKGROUND_FIELD
#include "nonorthogonal_background.h"
#endif

// Sets the configuration parameters for the problem.
void set_problem_config(void);

// Initializes the problem, including setting up initial conditions.
void problem_init(void);

#ifdef GAMERA_TIME_DEPENDENT_WIND
/* Load/free optional time-dependent outer-boundary data around restarts too. */
int problem_runtime_init(void);
void problem_runtime_finalize(void);
#endif

// Applies boundary conditions to the problem domain.
void boundary_conditions(void);
#ifdef GAMERA_NONORTHOGONAL_BACKEND
/* Fill local Cartesian vertex arrays, including four geometry halos. */
int problem_grid_init(void);
#endif
#ifdef GAMERA_NONORTHOGONAL_HAS_RADIAL_MAP
#include <stddef.h>

#include "nonorthogonal_radial_map.h"

/* Map a physical radius to the global fractional cell-center index. */
int problem_nonorthogonal_radial_logical(double radius, double *logical);
/* Analytic-map identity persisted in checkpoints to prevent mixed grids. */
int problem_nonorthogonal_radial_map_version(void);
double problem_nonorthogonal_radial_stretch(void);
int problem_nonorthogonal_radial_map_parameters(
    double parameters[GAMERA_NO_RADIAL_MAP_PARAMETER_COUNT]);
#endif
#ifdef GAMERA_NONORTHOGONAL_HAS_CELL_SOURCE
/* Conservative rate evaluated from the non-orthogonal half-time state. */
int problem_nonorthogonal_cell_source(
    gamera_no_vec3 point,
    const double predicted[GAMERA_NO_FLUX_COUNT], double time,
    void *context, double rate[GAMERA_NO_FLUX_COUNT]);
#endif
#ifdef GAMERA_NONORTHOGONAL_HAS_FLUID_FLUX_BOUNDARY
int problem_nonorthogonal_fluid_flux_boundary(
    gamera_no_storage *storage, const gamera_no_grid *grid,
    const size_t active_lower[3], const size_t active_upper[3]);
#endif
#ifdef GAMERA_NONORTHOGONAL_HAS_EDGE_EMF_BOUNDARY
int problem_nonorthogonal_edge_emf_boundary(
    gamera_no_storage *storage, const gamera_no_grid *grid,
    const size_t active_lower[3], const size_t active_upper[3]);
#endif
#ifdef GAMERA_NONORTHOGONAL_HAS_MAGNETIC_BOUNDARY
int problem_nonorthogonal_magnetic_boundary(gamera_no_storage *storage,
                                             const gamera_no_grid *grid);
#endif
#ifdef GAMERA_NONORTHOGONAL_HAS_BACKGROUND_FIELD
int problem_nonorthogonal_background_field(gamera_no_vec3 point,
                                            void *context,
                                            gamera_no_vec3 *field);
#ifdef GAMERA_NONORTHOGONAL_HAS_BACKGROUND_ADJUST
int problem_nonorthogonal_adjust_background(
    const gamera_no_grid *grid, gamera_no_background_data *background);
#endif
#ifdef GAMERA_NONORTHOGONAL_HAS_BACKGROUND_RESIDUAL_INIT
/* Convert an initialized total CT flux to B1=Btotal-B0 after B0 is built. */
int problem_nonorthogonal_initialize_background_residual(
    const gamera_no_grid *grid, gamera_no_storage *storage,
    const gamera_no_background_data *background);
#endif
#endif
#endif
