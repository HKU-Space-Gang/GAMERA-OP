#ifndef COMMON_H
#define COMMON_H

// For common functions of problems
// Field functions
typedef double (*Field)(double x, double y, double z);
extern Field A1, A2, A3, B10, B20, B30, Source1_B0, Source2_B0, Source3_B0;

// Function to define initialized magnetic field
double Field_Default(double x, double y, double z);

// Function to define background magnetic field
double B10_uniform(double x, double y, double z);
double B20_uniform(double x, double y, double z);
double B30_uniform(double x, double y, double z);

//============================================================
// Common boundary conditions
void gas_bc_symmetric_i_low();
void gas_bc_symmetric_i_high();
void gas_bc_symmetric_j_low();
void gas_bc_symmetric_j_high();
void gas_bc_symmetric_k_low();
void gas_bc_symmetric_k_high();

void gas_bc_reflective_i_low();
void gas_bc_reflective_i_high();
void gas_bc_reflective_j_low();
void gas_bc_reflective_j_high();
void gas_bc_reflective_k_low();
void gas_bc_reflective_k_high();

void gas_bc_fixed_i_low();
void gas_bc_fixed_i_high();
void gas_bc_fixed_j_low();
void gas_bc_fixed_j_high();
void gas_bc_fixed_k_low();
void gas_bc_fixed_k_high();

void gas_bc_extrapolated_i_low();
void gas_bc_extrapolated_i_high();
void gas_bc_extrapolated_j_low();
void gas_bc_extrapolated_j_high();
void gas_bc_extrapolated_k_low();
void gas_bc_extrapolated_k_high();

void gas_bc_none_i_low();
void gas_bc_none_i_high();
void gas_bc_none_j_low();
void gas_bc_none_j_high();
void gas_bc_none_k_low();
void gas_bc_none_k_high();


void gem_bc_symmetric_i_low();
void gem_bc_symmetric_i_high();
void gem_bc_symmetric_j_low();
void gem_bc_symmetric_j_high();
void gem_bc_symmetric_k_low();
void gem_bc_symmetric_k_high();

void gem_bc_reflective_i_low();
void gem_bc_reflective_i_high();
void gem_bc_reflective_j_low();
void gem_bc_reflective_j_high();
void gem_bc_reflective_k_low();
void gem_bc_reflective_k_high();

void gem_bc_fixed_i_low();
void gem_bc_fixed_i_high();
void gem_bc_fixed_j_low();
void gem_bc_fixed_j_high();
void gem_bc_fixed_k_low();
void gem_bc_fixed_k_high();

void gem_bc_extrapolated_i_low();
void gem_bc_extrapolated_i_high();
void gem_bc_extrapolated_j_low();
void gem_bc_extrapolated_j_high();
void gem_bc_extrapolated_k_low();
void gem_bc_extrapolated_k_high();

void gem_bc_none_i_low();
void gem_bc_none_i_high();
void gem_bc_none_j_low();
void gem_bc_none_j_high();
void gem_bc_none_k_low();
void gem_bc_none_k_high();

#endif
