#ifdef PROBLEM_OT_cartesian
#include <math.h>

#include "log.h"

#include "common.h"
#include "config.h"
#include "problem.h"
#include "solver.h"

void set_problem_config() {
  // basic settings
  // Cartesian Domain: x1, x2, x3 are the three Cartesian coordinates
  x1min_global = 0.0;
  x1max_global = 1.0;
  x2min_global = 0.0;
  x2max_global = 1.0;
  x3min_global = 0.0;
  x3max_global = 1.0;

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

// OT2 Az field function
double Field_OT2_Az(double x, double y, double z) {
  return (cos(4 * PI * x) / PI / 4 + cos(2 * PI * y) / PI / 2) / sqrt(PI * 4);
}

void problem_init() {
  // the Orszag-Tang problem, x-y plane
  log_info("Initializing Orszag-Tang problem");
  A3 = Field_OT2_Az;
  Field B10_init = Field_Default;
  Field B20_init = Field_Default;
  Field B30_init = Field_Default;
  for (int s = 0; s < NS1; s++) {
    for (int i = is; i <= ie; i++) {
      for (int j = js; j <= je; j++) {
        for (int k = ks; k <= ke; k++) {
          gas[s][gas_rho][i][j][k] = 25.0 / (PI * 36.0);
          gas[s][gas_p][i][j][k] = 5.0 / (PI * 12.0);
          // /* for i-j
          gas[s][gas_v1][i][j][k] = -sin(2 * PI * x2c[i][j][k]);
          gas[s][gas_v2][i][j][k] = sin(2 * PI * x1c[i][j][k]);
          gas[s][gas_v3][i][j][k] = 0.0;
          //*/
          /*for j-k
          gas[s][gas_v2][i][j][k] = -sin(2*PI*(*x3c)[i][j][k]);
          gas[s][gas_v3][i][j][k] =  sin(2*PI*x2c[i][j][k]);
          gas[s][gas_v1][i][j][k] = 0.0;
          */
          /* for k-i
          gas[s][gas_v3][i][j][k] = -sin(2*PI*x1c[i][j][k]);
          gas[s][gas_v1][i][j][k] =  sin(2*PI*(*x3c)[i][j][k]);
          gas[s][gas_v2][i][j][k] = 0.0;
          */
        }
      }
    }
  }

  for (int i = is; i <= ie + 1; i++) {
    for (int j = js; j <= je; j++) {
      for (int k = ks; k <= ke; k++) {
        // for i-j
        gem[mag_bi][i][j][k] =
            -1.0 * sin(2.0 * PI * x2c[i][j][k]) / sqrt(4 * PI);
        // for j-k
        // gem[mag_bi][i][j][k] = 0.0;
        //  for k-i
        // gem[mag_bi][i][j][k]
        //= 1.0*sin(4.0*PI*(*x3c)[i][j][k])/sqrt(4*PI);
      }
    }
  }
  for (int i = is; i <= ie; i++) {
    for (int j = js; j <= je + 1; j++) {
      for (int k = ks; k <= ke; k++) {
        // for i-j
        gem[mag_bj][i][j][k] =
            1.0 * sin(4.0 * PI * x1c[i][j][k]) / sqrt(4 * PI);
        // for j-k
        // gem[mag_bj][i][j][k] =
        //-1.0*sin(2.0*PI*(*x3c)[i][j][k])/sqrt(4*PI);
        //  for k-i
        // gem[mag_bj][i][j][k] = 0.0;
      }
    }
  }
  for (int i = is; i <= ie; i++) {
    for (int j = js; j <= je; j++) {
      for (int k = ks; k <= ke + 1; k++) {
        // for i-j
        gem[mag_bk][i][j][k] = 0.0;
        // for j-k
        // gem[mag_bk][i][j][k]
        //= 1.0*sin(4.0*PI*x2c[i][j][k])/sqrt(4*PI);
        //  for k-i
        // gem[mag_bk][i][j][k] =
        //-1.0*sin(2.0*PI*x1c[i][j][k])/sqrt(4*PI);
      }
    }
  }
  AinitB(gem[mag_bi], gem[mag_bj], gem[mag_bk], A1, A2, A3, B10_init, B20_init, B30_init);
}

#endif