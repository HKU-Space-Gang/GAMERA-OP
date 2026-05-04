#ifdef PROBLEM_FieldLoop_spherical
#include <math.h>

#include "log.h"

#include "common.h"
#include "config.h"
#include "problem.h"
#include "solver.h"
#include "curvilinear.h"

void set_problem_config() {
  // basic settings
  // Spherical Domain: x1 is the radial coordinate, x2 is the polar angle (theta), x3 is the azimuthal angle (phi)
  x1min_global = 0.1;
  x1max_global = 2.0;
  x2min_global = 0.0/10.0*PI;
  x2max_global = 5.0/10.0*PI;
  x3min_global = 0.0/5.0*PI;
  x3max_global = 10.0/5.0*PI;

  CA = 10.0; // Light speed in the normalized unit
  // Ringavarage option
  // 0 - no ring average; 1 - do ring average
  // Note that if you want to use pole_boundary in spherical coordinate, you must set doRingAverage = 1
  doRingAverage = 1;

  problem_config.Boundary_i = BC_NON_PERIODIC;
  problem_config.Boundary_j = BC_NON_PERIODIC;
  problem_config.Boundary_k = BC_PERIODIC;
}


void boundary_conditions() {
  // If problem_config.Boundary_i == BC_PERIODIC, then we don't need to set boundary conditions for i direction. The same fot j and k direction.
  // If you want to fix the boudary condition as initialization, we don't need to set boundary conditions

//  gas_bc_symmetric_i_low();
//  gas_bc_symmetric_i_high();
//  gem_bc_symmetric_i_low();
//  gem_bc_symmetric_i_high();

  gas_bc_pole_j_low();
  gem_bc_pole_j_low();
//  gas_bc_pole_j_high();
//  gem_bc_pole_j_high();
  // gas_bc_symmetric_j_high();
  // gem_bc_symmetric_j_high();
}


double Field_Az(double x, double y, double z) {
  double B0 = 1e-3;
  double sigma = 0.2, R = 0.5;
  double x0 = -sqrt(2)/2.0, y0 = 0, z0 = sqrt(2.0)/2.0;
  double Az = B0*exp(-(z-z0)*(z-z0)/sigma/sigma)*fmax(0, R - sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)));
  return Az;
}

double Field_A1(double x1, double x2, double x3) {
  double x, y, z, A1, A2, A3, Ax, Ay, Az;
  curvilinear_to_cartesian_coord(x1, x2, x3, &x, &y, &z);
  Ax = 0;
  Ay = 0;
  Az = Field_Az(x, y, z);
  cartesian_to_curvilinear_vector(x, y, z, Ax, Ay, Az, &A1, &A2, &A3);
  return A1;
}


double Field_A2(double x1, double x2, double x3) {
  double x, y, z, A1, A2, A3, Ax, Ay, Az;
  curvilinear_to_cartesian_coord(x1, x2, x3, &x, &y, &z);
  Ax = 0;
  Ay = 0;
  Az = Field_Az(x, y, z);
  cartesian_to_curvilinear_vector(x, y, z, Ax, Ay, Az, &A1, &A2, &A3);
  return A2;
}


double Field_A3(double x1, double x2, double x3) {
  double x, y, z, A1, A2, A3, Ax, Ay, Az;
  curvilinear_to_cartesian_coord(x1, x2, x3, &x, &y, &z);
  Ax = 0;
  Ay = 0;
  Az = Field_Az(x, y, z);
  cartesian_to_curvilinear_vector(x, y, z, Ax, Ay, Az, &A1, &A2, &A3);
  return A3;
}

void problem_init() {
  log_info("Initializing Field Loop problem");
  A1 = Field_A1;
  A2 = Field_A2;
  A3 = Field_A3;
  Field B10_init = Field_Default;
  Field B20_init = Field_Default;
  Field B30_init = Field_Default;
  for (int s = 0; s < NS1; s++) {
    for (int i = isg; i <= ieg; i++) {
      for (int j = jsg; j <= jeg; j++) {
        for (int k = ksg; k <= keg; k++) {
          gas[s][gas_rho][i][j][k] = 1.0;

          double x, y, z, vx, vy, vz;
          vx = 1.0;
          vy = 0.0;
          vz = 0.0;
          // vx = 1.0 in Cartesian coordinates, we need to convert it to curvilinear coordinates
          curvilinear_to_cartesian_coord(x1ctr[i][j][k], x2ctr[i][j][k], x3ctr[i][j][k], &x, &y, &z);
          cartesian_to_curvilinear_vector(x, y, z, vx, vy, vz, &gas[s][gas_v1][i][j][k], &gas[s][gas_v2][i][j][k], &gas[s][gas_v3][i][j][k]);

          gas[s][gas_p][i][j][k] = 0.5;
        }
      }
    }
  }

  AinitB(gem[mag_bi], gem[mag_bj], gem[mag_bk], A1, A2, A3, B10_init, B20_init, B30_init);
}

#endif