#ifdef PROBLEM_Blast_cartesian
#include <math.h>

#include "log.h"

#include "common.h"
#include "config.h"
#include "problem.h"

void set_problem_config() {
  // basic settings
  // Cartesian Domain: x1, x2, x3 are the three Cartesian coordinates
  x1min_global = -0.5;
  x1max_global = 0.5;
  x2min_global = -0.5;
  x2max_global = 0.5;
  x3min_global = -0.5;
  x3max_global = 0.5;

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
  // if problem_config.Boundary_i == BC_PERIODIC, then we don't need to set boundary conditions for i direction. The same fot j and k direction.
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

void problem_init() {
  log_info("Initializing Blast problem");
  for (int s = 0; s < NS1; s++) {
    for (int i = is; i <= ie; i++) {
      for (int j = js; j <= je; j++) {
        for (int k = ks; k <= ke; k++) {
          gas[s][gas_rho][i][j][k] = 1.0;
          // Assuming vx refers to the x-direction velocity
          gas[s][gas_v1][i][j][k] = 0.0;
          // Assuming vy refers to the y-direction velocity
          gas[s][gas_v2][i][j][k] = 0.0;
          // Assuming vz refers to the z-direction velocity
          gas[s][gas_v3][i][j][k] = 0.0;

          // Calculating the radius in the x-y plane
          double radius =
              sqrt((x1c[i][j][k] - 0.5) * (x1c[i][j][k] - 0.5) +
                   (x2c[i][j][k] - 0.5) * (x2c[i][j][k] - 0.5) +
                   (x3c[i][j][k] - 0.5) * (x3c[i][j][k] - 0.5));

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

  for (int i = is; i <= ie + 1; i++) {
    for (int j = js; j <= je; j++) {
      for (int k = ks; k <= ke; k++) {
        gem[mag_bi][i][j][k] = 1.0 / sqrt(2);
      }
    }
  }

  for (int i = is; i <= ie; i++) {
    for (int j = js; j <= je + 1; j++) {
      for (int k = ks; k <= ke; k++) {
        gem[mag_bj][i][j][k] = 1.0 / sqrt(2);
      }
    }
  }

  for (int i = is; i <= ie; i++) {
    for (int j = js; j <= je; j++) {
      for (int k = ks; k <= ke + 1; k++) {
        gem[mag_bk][i][j][k] = 0.0;
      }
    }
  }
}

#endif