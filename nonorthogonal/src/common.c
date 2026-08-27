#include "common.h"

#include <math.h>

#include "config.h"

// Default field fuction: 0
double Field_Default(double x, double y, double z) { return 0.0; }
// Field configuration including initialized magnetic field and background
// magnetic field Default value is 0. If you want to set other magnetic field,
// you should change in file "problem_init.c"
Field A1 = Field_Default;
Field A2 = Field_Default;
Field A3 = Field_Default;
Field B10 = Field_Default;
Field B20 = Field_Default;
Field B30 = Field_Default;
Field Source1_B0 = Field_Default;
Field Source2_B0 = Field_Default;
Field Source3_B0 = Field_Default;

//============================================================
// Function to define background magnetic field
// Uniform B field
double B10_uniform(double x, double y, double z) { return 1 / sqrt(2); }
double B20_uniform(double x, double y, double z) { return 1 / sqrt(2); }
double B30_uniform(double x, double y, double z) { return 1 / sqrt(2); }
// Dipole field
double B10_dipole(double x, double y, double z) {
  double r_squared = x * x + y * y + z * z;
  return (3 * x * z) / pow(r_squared, 2.5);
}
double B20_dipole(double x, double y, double z) {
  double r_squared = x * x + y * y + z * z;
  return (3 * y * z) / pow(r_squared, 2.5);
}
double B30_dipole(double x, double y, double z) {
  double r_squared = x * x + y * y + z * z;
  return (2 * z * z - x * x - y * y) / pow(r_squared, 2.5);
}
