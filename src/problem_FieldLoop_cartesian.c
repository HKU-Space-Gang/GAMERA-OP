#ifdef PROBLEM_FieldLoop_cartesian
#include <math.h>

#include "log.h"

#include "common.h"
#include "config.h"
#include "problem.h"
#include "solver.h"

void set_problem_config() {
  // basic settings
  // Cartesian Domain: x1, x2, x3 are the three Cartesian coordinates
  x1min_global = -2.0;
  x1max_global = 2.0;
  x2min_global = -2.0;
  x2max_global = 2.0;
  x3min_global = -2.0;
  x3max_global = 2.0;

  CA = 10.0; // Light speed in the normalized unit
  // Ringavarage option
  // 0 - no ring average; 1 - do ring average
  // Note that if you want to use pole_boundary in spherical coordinate, you must set doRingAverage = 1
  doRingAverage = 0;

  problem_config.Boundary_i = BC_PERIODIC;
  problem_config.Boundary_j = BC_PERIODIC;
  problem_config.Boundary_k = BC_PERIODIC;
}


void boundary_conditions() {
  // If problem_config.Boundary_i == BC_PERIODIC, then we don't need to set boundary conditions for i direction. The same fot j and k direction.
//  gas_bc_symmetric_i_low();
//  gas_bc_symmetric_i_high();
//  gas_bc_symmetric_j_low();
//  gas_bc_symmetric_j_high();
//  gas_bc_symmetric_k_low();
//  gas_bc_symmetric_k_high();
//
//  gem_bc_symmetric_i_low();
//  gem_bc_symmetric_i_high();
//  gem_bc_symmetric_j_low();
//  gem_bc_symmetric_j_high();
//  gem_bc_symmetric_k_low();
//  gem_bc_symmetric_k_high();
}


double Field_Az(double x, double y, double z) {
  double B0 = 1e-3;
  double sigma = 0.2, R = 0.5;
  double x0 = -sqrt(2)/2.0, y0 = 0, z0 = sqrt(2.0)/2.0;
  double Az = B0*exp(-(z-z0)*(z-z0)/sigma/sigma)*fmax(0, R - sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)));
  return Az;
}

void problem_init() {
  log_info("Initializing Field Loop problem");
  A3 = Field_Az;
  Field B10_init = Field_Default;
  Field B20_init = Field_Default;
  Field B30_init = Field_Default;
  for (int s = 0; s < NS1; s++) {
    for (int i = is; i <= ie; i++) {
      for (int j = js; j <= je; j++) {
        for (int k = ks; k <= ke; k++) {
          gas[s][gas_rho][i][j][k] = 1.0;
          // Assuming vx refers to the x-direction velocity
          gas[s][gas_v1][i][j][k] = 1.0;
          // Assuming vy refers to the y-direction velocity
          gas[s][gas_v2][i][j][k] = 0.0;
          // Assuming vz refers to the z-direction velocity
          gas[s][gas_v3][i][j][k] = 0.0;

          gas[s][gas_p][i][j][k] = 0.5;
        }
      }
    }
  }

  AinitB(gem[mag_bi], gem[mag_bj], gem[mag_bk], A1, A2, A3, B10_init, B20_init, B30_init);
}

#endif