#ifndef SOLVER_H
#define SOLVER_H

#include <mpi.h>
#include <stdbool.h>

#include "config.h"
#include "common.h"

void AinitB(double ***magi, double ***magj, double ***magk, Field A1, Field A2,
            Field A3, Field B10, Field B20, Field B30);
double GaussianLineIntegral(double (*fx)(double, double, double), double xa,
                            double ya, double za, double xb, double yb,
                            double zb);

// Initialize the solver with base configuration
int initialize_solver(const char *base_name);

// Allocate necessary memory and resources for the solver
int allocate_solver();

// Set up domain decomposition for parallel processing
int set_domain_decomposition();

// Free allocated resources and memory
int finalize_solver();

// Calculate time step based on CFL condition
double get_dt();

// Get derived variables
void get_derived_variables();

// Get variables needed for AB2 time-stepping scheme
void get_AB2_variables();

// Reconstruct variables for magnetic field variables
void reconstruct_3dv_gem(dir_t dir);

// Reconstruct variables for gas variables
void reconstruct_3dv_gas(int dir);

// Calculate electric fields based on reconstructed variables
void get_e_fields();

// Calculate electric fields in i-direction
void get_e_fields_i();

// Calculate electric fields in j-direction
void get_e_fields_j();

// Calculate electric fields in k-direction
void get_e_fields_k();

// Calculate fluid flux using Rusanov's scheme
void get_fluid_flux_rusanov(dir_t dir);

// Calculate magnetic stress using Rusanov's scheme
void get_magnetic_stress_rusanov(dir_t dir);

// Calculate fluid flux using Gk scheme
void get_fluid_flux_gk(int dir);

// main function
int solve(int nt);

// Set background magnetic field
void set_background_field();

// Check the positivity of density and pressure, and reset them if necessary
void reset_rho();
void reset_p();
void check_positivity();

// Math utility functions
double factorial(int k);
double nchoosek(int j, int k);
double minmod(double x, double y);
double DipoleL(double r, double theta, double phi);

// Ring average functions
void PPM(double fm2, double fm1, double f, double fp1, double fp2, double* vm, double* vp);
void FilterBaseMode(double* Q, int N, double* Q0);
void RingBlockAvg(double* Q, int N, int Nchunk);

#endif
