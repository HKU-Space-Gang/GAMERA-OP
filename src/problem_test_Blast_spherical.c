#ifdef PROBLEM_Blast_spherical
#include <math.h>

#include "log.h"

#include "common.h"
#include "config.h"
#include "problem.h"
#include "curvilinear.h"
#include "solver.h"

void set_problem_config() {
  // basic settings
  // Spherical Domain: x1 is the radial coordinate, x2 is the polar angle (theta), x3 is the azimuthal angle (phi)
  x1min_global = 0.5;
  x1max_global = 1.5;
  x2min_global = 3.0/10.0*PI;
  x2max_global = 7.0/10.0*PI;
  x3min_global = -1.0/5.0*PI;
  x3max_global = 1.0/5.0*PI;

  CA = 10.0; // Light speed in the normalized unit
  // Ringavarage option
  // 0 - no ring average; 1 - do ring average
  // Note that if you want to use pole_boundary in spherical coordinate, you must set doRingAverage = 1
  doRingAverage = 0;

  problem_config.Boundary_i = BC_NON_PERIODIC;
  problem_config.Boundary_j = BC_NON_PERIODIC;
  problem_config.Boundary_k = BC_NON_PERIODIC;
}

void boundary_conditions() {
  // if problem_config.Boundary_i == BC_PERIODIC, then we don't need to set boundary conditions for i direction. The same fot j and k direction.
  gas_bc_symmetric_i_low();
  gas_bc_symmetric_i_high();
  gas_bc_symmetric_j_low();
  gas_bc_symmetric_j_high();
  gas_bc_symmetric_k_low();
  gas_bc_symmetric_k_high();

  gem_bc_symmetric_i_low();
  gem_bc_symmetric_i_high();
  gem_bc_symmetric_j_low();
  gem_bc_symmetric_j_high();
  gem_bc_symmetric_k_low();
  gem_bc_symmetric_k_high();
}


double Field_A1(double x1, double x2, double x3) {
  double A1 = 1.0/sqrt(2.0)*x1*sin(x2)*cos(x2)*(sin(x3)-cos(x3));
  return A1;
}

double Field_A2(double x1, double x2, double x3) {
  double A2 = -1.0/sqrt(2.0)*x1*sin(x2)*sin(x2)*(sin(x3)-cos(x3));
  return A2;
}

void problem_init() {
  log_info("Initializing Blast problem");
  A1 = Field_A1;
  A2 = Field_A2;
  Field B10_init = Field_Default;
  Field B20_init = Field_Default;
  Field B30_init = Field_Default;
  for (int s = 0; s < NS1; s++) {
    for (int i = is; i <= ie; i++) {
      for (int j = js; j <= je; j++) {
        for (int k = ks; k <= ke; k++) {
          gas[s][gas_rho][i][j][k] = 1.0;

          gas[s][gas_v1][i][j][k] = 0.0;
          gas[s][gas_v2][i][j][k] = 0.0;
          gas[s][gas_v3][i][j][k] = 0.0;

          // Calculating the radius
          double x_temp = x1ctr[i][j][k]*sin(x2ctr[i][j][k])*cos(x3ctr[i][j][k]);
          double y_temp = x1ctr[i][j][k]*sin(x2ctr[i][j][k])*sin(x3ctr[i][j][k]);
          double z_temp = x1ctr[i][j][k]*cos(x2ctr[i][j][k]);
          double radius =
              sqrt((x_temp-1.0) * (x_temp-1.0) +
                   (y_temp) * (y_temp) +
                   (z_temp) * (z_temp));

          if (radius < 0.1) {
            // High pressure inside the blast radius
            gas[s][gas_p][i][j][k] = 10.0;
          } else {
            // Low pressure outside the blast radius
            gas[s][gas_p][i][j][k] = 0.1;
          }
        }
      }
    }
  }
  AinitB(gem[mag_bi], gem[mag_bj], gem[mag_bk], A1, A2, A3, B10_init, B20_init, B30_init);
}
#endif