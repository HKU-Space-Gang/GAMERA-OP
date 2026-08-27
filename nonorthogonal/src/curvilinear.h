#ifndef CURVILINEAR_H
#define CURVILINEAR_H

#include <stdbool.h>
#include <stdio.h>
#include "config.h"

// Ring related arrays and variables
extern int Ring_Ni, Ring_Nj, Ring_Nk;

extern double ****geo_ring;
enum {
  face_idir_ring = 0,
  face_jdir_ring,
  face_kdir_ring,
  NF_geo_ring
};

extern double ****gem_ring;
enum {
  mag_bi_ring = 0,
  mag_bj_ring,
  mag_bk_ring,

  mags_b1_ring,
  mags_b2_ring,
  mags_b3_ring,

  mag_bi_flux_ring,
  mag_bj_flux_ring,
  mag_bk_flux_ring,

  // Here, efield_tackle_ei_ring is ring average of electric field used in tackle_Efield function
  efield_tackle_ei_ring,
  efield_tackle_ej_ring,
  efield_tackle_ek_ring,

  // Here, efield_ring is ring average of electric field
  efield_ei_ring,
  efield_ej_ring,
  efield_ek_ring,
  NF_gem_ring
};

extern double *****gas_ring;
enum {
  gas_rho_ring = 0,
  gas_v1_ring,
  gas_v2_ring,
  gas_v3_ring,
  gas_p_ring,
  gas_p_S_ring,
  NF_gas_ring
};

// Lame coefficients
double H1(double x1, double x2, double x3);
double H2(double x1, double x2, double x3);
double H3(double x1, double x2, double x3);

// Calculate geometry arrays for curvilinear coordinates
int set_geometry_arrays();

// Calculate reconstruction weights for curvilinear coordinates
void getweights();

// Calculate source term for curvilinear coordinates
void Source_Term(dir_t dir);

// Ring average related functions in curvilinear coordinates
void gather_Etackle_RingAverage_arrays();
void broadcast_Etackle_RingAverage_arrays();
int RingAverage_gather_data_3d_array(double ***field_array_scr, double ***field_array_dst);
int RingAverage_broadcast_data_3d_array(double ***field_array_scr, double ***field_array_dst);

// Special operation in pole regions for curvilinear coordinates
void tackle_Efield_pole();
void tackle_Magfield_pole();

// Set configuration for ring average in curvilinear coordinates
void init_RingAverage_variables();
void set_RingAverage_geo_data();
// Spherical: ring average in k direction. The chunk number depends on j direction
void set_RingAverage_config();

// Get time step constraint from ring average in curvilinear coordinates
double get_dt_ring();

// Ring average function
void gather_value_RingAverage_arrays();
void HydroRingAverage();
void MagneticRingAverage();
void broadcast_value_RingAverage_arrays();

// specific boundary conditions for curvilinear coordinates
void gas_bc_pole_j_low();
void gas_bc_pole_j_high();
void gem_bc_pole_j_low();
void gem_bc_pole_j_high();
void gas_bc_solar_wind();
void gem_bc_solar_wind();

void InnerEFix();
void WindEFix();

// coordinate transformation functions
void cartesian_to_curvilinear_coord(double x, double y, double z,
                                    double* x1, double* x2, double* x3);

void curvilinear_to_cartesian_coord(double x1, double x2, double x3,
                                    double* x, double* y, double* z);

void cartesian_to_curvilinear_vector(double x, double y, double z,
                                     double vx, double vy, double vz,
                                     double* v1, double* v2, double* v3);

void curvilinear_to_cartesian_vector(double x1, double x2, double x3,
                                     double v1, double v2, double v3,
                                     double* vx, double* vy, double* vz);

#endif // CURVILINEAR_H
