#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "common.h"
#include "log.h"
#include "solver.h"
#include "config.h"

#include "curvilinear.h"
#include "spsolver.h"
#include "lu_solver.h"

void AinitB(double ***magi, double ***magj, double ***magk, Field A1, Field A2,
            Field A3, Field B10, Field B20, Field B30)  {
  for (int i = isg; i <= ieg; i++) {
    for (int j = jsg; j <= jeg + 1; j++) {
      for (int k = ksg; k <= keg + 1; k++) {
        gem[LAi][i][j][k] =
            geo[edge_idir][i][j][k] *
            GaussianLineIntegral(*A1, x1[i + 1][j][k], x2[i + 1][j][k],
                                 x3[i + 1][j][k], x1[i][j][k], x2[i][j][k],
                                 x3[i][j][k]);
      }
    }
  }
  for (int i = isg; i <= ieg + 1; i++) {
    for (int j = jsg; j <= jeg; j++) {
      for (int k = ksg; k <= keg + 1; k++) {
        gem[LAj][i][j][k] =
            geo[edge_jdir][i][j][k] *
            GaussianLineIntegral(*A2, x1[i][j + 1][k], x2[i][j + 1][k],
                                 x3[i][j + 1][k], x1[i][j][k], x2[i][j][k],
                                 x3[i][j][k]);
      }
    }
  }
  for (int i = isg; i <= ieg + 1; i++) {
    for (int j = jsg; j <= jeg + 1; j++) {
      for (int k = ksg; k <= keg; k++) {
        gem[LAk][i][j][k] =
            geo[edge_kdir][i][j][k] *
            GaussianLineIntegral(*A3, x1[i][j][k + 1], x2[i][j][k + 1],
                                 x3[i][j][k + 1], x1[i][j][k], x2[i][j][k],
                                 x3[i][j][k]);
      }
    }
  }

  for (int i = isg; i <= ieg + 1; i++) {
    for (int j = jsg; j <= jeg; j++) {
      for (int k = ksg; k <= keg; k++) {
        if (geo[face_idir][i][j][k] < 1e-10) {
//          magi[i][j][k] = 0.0;
          magi[i][j][k] = GaussianLineIntegral(B10, x1[i][j][k], x2[i][j][k], x3[i][j][k], x1[i][j][k+1], x2[i][j][k+1], x3[i][j][k+1]);
        } else{
          magi[i][j][k] = (gem[LAj][i][j][k] - gem[LAj][i][j][k + 1] +
              gem[LAk][i][j + 1][k] - gem[LAk][i][j][k]) /
              geo[face_idir][i][j][k];
        }
      }
    }
  }

  for (int i = isg; i <= ieg; i++) {
    for (int j = jsg; j <= jeg + 1; j++) {
      for (int k = ksg; k <= keg; k++) {
        if (geo[face_jdir][i][j][k] < 1e-10) {
//          magj[i][j][k] = 0.0;
          magj[i][j][k] = GaussianLineIntegral(B20, x1[i][j][k], x2[i][j][k], x3[i][j][k], x1[i+1][j][k], x2[i+1][j][k], x3[i+1][j][k]);
        } else {
          magj[i][j][k] = -(gem[LAi][i][j][k] - gem[LAi][i][j][k + 1] +
              gem[LAk][i + 1][j][k] - gem[LAk][i][j][k]) /
              geo[face_jdir][i][j][k];
        }
      }
    }
  }

  for (int i = isg; i <= ieg; i++) {
    for (int j = jsg; j <= jeg; j++) {
      for (int k = ksg; k <= keg + 1; k++) {
        if (geo[face_kdir][i][j][k] < 1e-10) {
//          magk[i][j][k] = 0.0;
          magk[i][j][k] = GaussianLineIntegral(B30, x1[i][j][k], x2[i][j][k], x3[i][j][k], x1[i][j+1][k], x2[i][j+1][k], x3[i][j+1][k]);
        } else {
          magk[i][j][k] = (gem[LAi][i][j][k] - gem[LAi][i][j + 1][k] +
              gem[LAj][i + 1][j][k] - gem[LAj][i][j][k]) /
              geo[face_kdir][i][j][k];
        }
      }
    }
  }
}

double GaussianLineIntegral(double (*fx)(double, double, double), double xa,
                            double ya, double za, double xb, double yb,
                            double zb) {
  // Positive zeros of 12th order Legendre polynomial
  double A[] = {0.1252334085, 0.3678314989, 0.5873179542,
                0.7699026741, 0.9041172563, 0.9815606342};
  // Gaussian Integration coefficients for a 12th order polynomial
  double WT[] = {0.2491470458, 0.2334925365, 0.2031674267,
                 0.1600783285, 0.1069393259, 0.0471753363};

  double dx = (xb - xa) / 2.0;
  double dy = (yb - ya) / 2.0;
  double dz = (zb - za) / 2.0;
  double xbar = (xa + xb) / 2.0;
  double ybar = (ya + yb) / 2.0;
  double zbar = (za + zb) / 2.0;

  double sum = 0;
  int k;

  for (k = 0; k < sizeof(A) / sizeof(A[0]); k++) {
    sum += WT[k] * (fx(xbar + A[k] * dx, ybar + A[k] * dy, zbar + A[k] * dz) +
                    fx(xbar - A[k] * dx, ybar - A[k] * dy, zbar - A[k] * dz));
  }

  sum /= 2.0;

  return sum;
}

// one dimensional integration, no sequence requirement for a and b point
double g1int_ortho(dir_t dir, double (*f)(double, double, double), double xa,
                   double ya, double za, double xb, double yb,
                   double zb, double edge_len) {
  // Positive zeros of 12th order Legendre polynomial
  double A[] = {0.125233408500000, 0.367831498900000, 0.587317954200000,
                0.769902674100000, 0.904117256300000, 0.981560634200000};
  // Gaussian Integration coefficients for a 12th order polynomial
  double WT[] = {0.249147045800000, 0.233492536500000, 0.203167426700000,
                 0.160078328500000, 0.106939325900000, 0.047175336300000};

  double dx = (xb - xa) / 2.0;
  double dy = (yb - ya) / 2.0;
  double dz = (zb - za) / 2.0;
  double xbar = (xa + xb) / 2.0;
  double ybar = (ya + yb) / 2.0;
  double zbar = (za + zb) / 2.0;

  double sum = 0.0;
  int k;

  switch (dir) {
    case i_dir:
      for (k = 0; k < sizeof(A) / sizeof(A[0]); k++) {
        sum += WT[k] * (f(xbar + A[k] * dx, ybar + A[k] * dy, zbar + A[k] * dz)* H1(xbar + A[k] * dx, ybar + A[k] * dy, zbar + A[k] * dz) +
            f(xbar - A[k] * dx, ybar - A[k] * dy, zbar - A[k] * dz)* H1(xbar - A[k] * dx, ybar - A[k] * dy, zbar - A[k] * dz));
      }
      break;
    case j_dir:
      for (k = 0; k < sizeof(A) / sizeof(A[0]); k++) {
        sum += WT[k] * (f(xbar + A[k] * dx, ybar + A[k] * dy, zbar + A[k] * dz)*H2(xbar + A[k] * dx, ybar + A[k] * dy, zbar + A[k] * dz) +
            f(xbar - A[k] * dx, ybar - A[k] * dy, zbar - A[k] * dz)*H2(xbar - A[k] * dx, ybar - A[k] * dy, zbar - A[k] * dz));
      }
      break;
    case k_dir:
      for (k = 0; k < sizeof(A) / sizeof(A[0]); k++) {
        sum += WT[k] * (f(xbar + A[k] * dx, ybar + A[k] * dy, zbar + A[k] * dz)*H3(xbar + A[k] * dx, ybar + A[k] * dy, zbar + A[k] * dz) +
            f(xbar - A[k] * dx, ybar - A[k] * dy, zbar - A[k] * dz)*H3(xbar - A[k] * dx, ybar - A[k] * dy, zbar - A[k] * dz));
      }
      break;
  }

  sum = sum* sqrt(dx*dx + dy*dy + dz*dz);
  if (edge_len > 1e-10){
    sum = sum / edge_len;
  } else {
    sum = 0.0;
  }

  return sum;
}


// Factorial calculation function
double factorial(int k) {
  if (k < 0) return 0.0;
  double result = 1.0;
  for(int i = 2; i <= k; ++i) {
    result *= i;
  }
  return result;
}

// Combination number calculation function
double nchoosek(int j, int k) {
  if(k < 0 || k > j) return 0.0;
  if(k == 0 || k == j) return 1.0;

  // Optimization calculation: C(j,k) = C(j, j-k)
  k = (k < j - k) ? k : j - k;

  double result = 1.0;
  for(int i = 1; i <= k; ++i) {
    result *= (j - k + i) / (double)i;
  }
  return result;
}

static double fsign(double x) {
  if (x > 0.0) {
    return 1.0;
  } else if (x < 0.0) {
    return -1.0;
  } else {
    return 0.0;
  }
}

static inline double PDMU7(double fh0[], double fp0[], double wgtPDMU[], double faceim1, double facei, double vol, double dx) {
  //double facPDMU7[7] = {-1.0/140, 5.0/84, -101.0/420, 319.0/420, 107.0/210, -19.0/210, 1.0/105};
  int size = 7, im1, i0, ip1;
  im1 = 2; i0 = im1+1; ip1 = i0+1;
  double f_itp = 0.0;

  for (int i = 0; i<size; i++){
    f_itp += fh0[i]*wgtPDMU[i];
  }

  double maxf = fmax(fp0[i0],fp0[ip1]);
  double minf = fmin(fp0[i0],fp0[ip1]);

  double f = fmax(minf,fmin(f_itp,maxf));

  double df0 = PDMB*(fp0[i0] - fp0[im1]);
  double df1 = PDMB*(fp0[ip1] - fp0[i0]);

  double s0 = fsign(df0);
  double s1 = fsign(df1);

  df0 = fabs(df0);

  double q0 = fabs(s0 + s1);

  double df_left = f - fp0[i0];

  double q = vol/dx;
  double e = 1.0/(2*PDMB+1.0);
  double df0c = fabs( (q-e*facei)/(e*facei)*fp0[i0] + (e*faceim1-q)/(e*facei)*fp0[im1] )/2.0;

  double a1 = PDMB/(PDMB+1.0);
  double a2 = 1.0-a1;
  df0c = a1*df0c + a2*df0;

  double f_pdm = f - s1*fmax(0 ,fabs(df_left)-q0*df0c);
  return f_pdm;
//   return f_itp;
}

double get_dt_normal() {
  double Vfluid, Vsound, B_tot_sq, Va, Va_eff, VCFL, dtCFL, minDtCFL = INFINITY;
  int i, j, k;

  // Assuming gas has dimensions [4][keg+1][jeg+1][ieg+1] and indexed as
  // gas[variable][k][j][i]
  #pragma omp parallel for collapse(4) schedule(static) reduction(min:minDtCFL)
  for (int s = 0; s < NS1; s++) {
    for (i = is; i <= ie; i++) {
      for (j = js; j <= je; j++) {
        for (k = ks; k <= ke; k++) {
          B_tot_sq = gem[magtot_b1][i][j][k] * gem[magtot_b1][i][j][k] +
                     gem[magtot_b2][i][j][k] * gem[magtot_b2][i][j][k] +
                     gem[magtot_b3][i][j][k] * gem[magtot_b3][i][j][k];

          Va = sqrt(B_tot_sq / gas[s][gas_rho][i][j][k]);
          Va_eff = Va * CA / (sqrt(Va * Va + CA * CA));

          Vfluid = sqrt(pow(gas[s][gas_v1][i][j][k], 2) +
                        pow(gas[s][gas_v2][i][j][k], 2) +
                        pow(gas[s][gas_v3][i][j][k], 2));
          Vsound = sqrt(gamma_val * gas[s][gas_p][i][j][k] /
                        gas[s][gas_rho][i][j][k]);

          VCFL = Vfluid + sqrt(Vsound * Vsound + Va_eff * Va_eff);
//          dtCFL = CFL / (fmax(VCFL / dx1[i][j][k],
//                              fmax(VCFL / dx2[i][j][k], VCFL / dx3[i][j][k])));
           dtCFL = CFL / (VCFL / geo[dxi][i][j][k] + VCFL / geo[dxj][i][j][k] +
                          VCFL / geo[dxk][i][j][k]);
          //          if (problem_config.initial_condition == IC_JUPITER) {
          //            double r = sqrt(x1c[i][j][k] * x1c[i][j][k] +
          //            x2c[i][j][k] * x2c[i][j][k] + x3c[i][j][k] *
          //            x3c[i][j][k]); if (r < planet_config.r_inner) {
          //              dtCFL = 100;
          //            }
          //          }
          minDtCFL = fmin(minDtCFL, dtCFL);
        }
      }
    }
  }

  return minDtCFL;
}

double get_dt() {
  if (doRingAverage == 1) {
    return get_dt_ring();
  } else {
    return get_dt_normal();
  }
}

void get_derived_variables() {
  // update derived magnetic field
  #pragma omp parallel for collapse(3) schedule(static)
  for (int i = isg; i <= ieg; i++) {
    for (int j = jsg; j <= jeg; j++) {
      for (int k = ksg; k <= keg; k++) {
        gem[mags_b1][i][j][k] = gem[mag_bi][i][j][k] * (x1[i+1][j][k] - x1ctr[i][j][k]) / dx1[i][j][k] \
                            + gem[mag_bi][i+1][j][k] * (x1ctr[i][j][k] - x1[i][j][k]) / dx1[i][j][k];
        gem[mags_b2][i][j][k] = gem[mag_bj][i][j][k] * (x2[i][j+1][k] - x2ctr[i][j][k]) / dx2[i][j][k] \
                            + gem[mag_bj][i][j+1][k] * (x2ctr[i][j][k] - x2[i][j][k]) / dx2[i][j][k];
        gem[mags_b3][i][j][k] = gem[mag_bk][i][j][k] * (x3[i][j][k+1] - x3ctr[i][j][k]) / dx3[i][j][k] \
                            + gem[mag_bk][i][j][k+1] * (x3ctr[i][j][k] - x3[i][j][k]) / dx3[i][j][k];
      }
    }
  }

  divB_max = 0;
  #pragma omp parallel for collapse(3) schedule(static)
  for (int i = is; i <= ie; i++) {
    for (int j = js; j <= je; j++) {
      for (int k = ks; k <= ke; k++) {
        double div_b;
        div_b = ((gem[mag_bi][i+1][j][k]+gem[B0AfaceI_b1][i+1][j][k])*geo[face_idir][i+1][j][k]
                - (gem[mag_bi][i][j][k]+gem[B0AfaceI_b1][i][j][k])*geo[face_idir][i][j][k])/geo[vol_center][i][j][k] +
                ((gem[mag_bj][i][j+1][k]+gem[B0AfaceJ_b2][i][j+1][k])*geo[face_jdir][i][j+1][k]
                - (gem[mag_bj][i][j][k]+gem[B0AfaceJ_b2][i][j][k])*geo[face_jdir][i][j][k])/geo[vol_center][i][j][k] +
                ((gem[mag_bk][i][j][k+1]+gem[B0AfaceK_b3][i][j][k+1])*geo[face_kdir][i][j][k+1]
                - (gem[mag_bk][i][j][k]+gem[B0AfaceK_b3][i][j][k])*geo[face_kdir][i][j][k])/geo[vol_center][i][j][k];

        if (fabs(div_b) > divB_max) {
          divB_max = fabs(div_b);
        }

      }
    }
  }

  #pragma omp parallel for collapse(3) schedule(static)
  for (int i = isg; i <= ieg; i++) {
    for (int j = jsg; j <= jeg; j++) {
      for (int k = ksg; k <= keg; k++) {
        gem[magtot_b1][i][j][k] = gem[mags_b1][i][j][k] + gem[B0_b1][i][j][k];
        gem[magtot_b2][i][j][k] = gem[mags_b2][i][j][k] + gem[B0_b2][i][j][k];
        gem[magtot_b3][i][j][k] = gem[mags_b3][i][j][k] + gem[B0_b3][i][j][k];
      }
    }
  }

  // update conservative variables
  #pragma omp parallel for collapse(4) schedule(static)
  for (int s = 0; s < NS1; s++) {
    for (int i = is; i <= ie; i++) {
      for (int j = js; j <= je; j++) {
        for (int k = ks; k <= ke; k++) {
          // TODO: gasc_* are internal variables, don't save
          gas[s][gasc_rho][i][j][k] = gas[s][gas_rho][i][j][k];
          gas[s][gasc_rhov1][i][j][k] =
              gas[s][gas_rho][i][j][k] * gas[s][gas_v1][i][j][k];
          gas[s][gasc_rhov2][i][j][k] =
              gas[s][gas_rho][i][j][k] * gas[s][gas_v2][i][j][k];
          gas[s][gasc_rhov3][i][j][k] =
              gas[s][gas_rho][i][j][k] * gas[s][gas_v3][i][j][k];
          gas[s][gasc_eng][i][j][k] =
              0.5 * gas[s][gas_rho][i][j][k] *
                  (gas[s][gas_v1][i][j][k] * gas[s][gas_v1][i][j][k] +
                   gas[s][gas_v2][i][j][k] * gas[s][gas_v2][i][j][k] +
                   gas[s][gas_v3][i][j][k] * gas[s][gas_v3][i][j][k]) +
              gas[s][gas_p][i][j][k] / (gamma_val - 1.0);
          gas[s][gasc_entropy][i][j][k] =
              gas[s][gas_p_S][i][j][k] /
              pow(gas[s][gas_rho][i][j][k], gamma_val-1.0);
        }
      }
    }
  }
}

void get_AB2_variables() {
  // get the half step variables
  #pragma omp parallel for collapse(4) schedule(static)
  for (int s = 0; s < NS1; s++) {
    for (int i = isg; i <= ieg; i++) {
      for (int j = jsg; j <= jeg; j++) {
        for (int k = ksg; k <= keg; k++) {
          gas[s][gas_rho_h][i][j][k] = gas[s][gas_rho][i][j][k]
              + 0.5 * dt / dt0 * (gas[s][gas_rho][i][j][k] - gas[s][gas_rho_p][i][j][k]);
          gas[s][gas_v1_h][i][j][k] = gas[s][gas_v1][i][j][k]
              + 0.5 * dt / dt0 * (gas[s][gas_v1][i][j][k] - gas[s][gas_v1_p][i][j][k]);
          gas[s][gas_v2_h][i][j][k] = gas[s][gas_v2][i][j][k]
              + 0.5 * dt / dt0 * (gas[s][gas_v2][i][j][k] - gas[s][gas_v2_p][i][j][k]);
          gas[s][gas_v3_h][i][j][k] = gas[s][gas_v3][i][j][k]
              + 0.5 * dt / dt0 * (gas[s][gas_v3][i][j][k] - gas[s][gas_v3_p][i][j][k]);
          gas[s][gas_p_h][i][j][k] = gas[s][gas_p][i][j][k]
              + 0.5 * dt / dt0 * (gas[s][gas_p][i][j][k] - gas[s][gas_p_p][i][j][k]);
          gas[s][gas_p_S_h][i][j][k] = gas[s][gas_p_S][i][j][k]
              + 0.5 * dt / dt0 * (gas[s][gas_p_S][i][j][k] - gas[s][gas_p_S_p][i][j][k]);
        }
      }
    }
  }

  #pragma omp parallel for collapse(3) schedule(static)
  for (int i = isg; i <= ieg; i++) {
    for (int j = jsg; j <= jeg; j++) {
      for (int k = ksg; k <= keg; k++) {
        gem[mags_b1_h][i][j][k] = gem[mags_b1][i][j][k]
            + 0.5 * dt / dt0 * (gem[mags_b1][i][j][k] - gem[mags_b1_p][i][j][k]);
        gem[mags_b2_h][i][j][k] = gem[mags_b2][i][j][k]
            + 0.5 * dt / dt0 * (gem[mags_b2][i][j][k] - gem[mags_b2_p][i][j][k]);
        gem[mags_b3_h][i][j][k] = gem[mags_b3][i][j][k]
            + 0.5 * dt / dt0 * (gem[mags_b3][i][j][k] - gem[mags_b3_p][i][j][k]);
      }
    }
  }

  #pragma omp parallel for collapse(3) schedule(static)
  for (int i = isg; i <= ieg+1; i++) {
    for (int j = jsg; j <= jeg; j++) {
      for (int k = ksg; k <= keg; k++) {
        gem[mag_bi_h][i][j][k] = gem[mag_bi][i][j][k]
            + 0.5 * dt / dt0 * (gem[mag_bi][i][j][k] - gem[mag_bi_p][i][j][k]);
      }
    }
  }
  #pragma omp parallel for collapse(3) schedule(static)
  for (int i = isg; i <= ieg; i++) {
    for (int j = jsg; j <= jeg+1; j++) {
      for (int k = ksg; k <= keg; k++) {
          gem[mag_bj_h][i][j][k] = gem[mag_bj][i][j][k]
              + 0.5 * dt / dt0 * (gem[mag_bj][i][j][k] - gem[mag_bj_p][i][j][k]);
      }
    }
  }
  #pragma omp parallel for collapse(3) schedule(static)
  for (int i = isg; i <= ieg; i++) {
    for (int j = jsg; j <= jeg; j++) {
      for (int k = ksg; k <= keg+1; k++) {
          gem[mag_bk_h][i][j][k] = gem[mag_bk][i][j][k]
              + 0.5 * dt / dt0 * (gem[mag_bk][i][j][k] - gem[mag_bk_p][i][j][k]);
      }
    }
  }

  // get pdm limiting variables
  int d_index = gas_rho_p - gas_rho; // index difference between gas and gas_p
  #pragma omp parallel for collapse(5) schedule(static)
  for (int s = 0; s < NS1; s++) {
    for (int m = gas_rho; m <= gas_p_S; m++) {
      for (int i = isg; i <= ieg; i++) {
        for (int j = jsg; j <= jeg; j++) {
          for (int k = ksg; k <= keg; k++) {
            gas[s][m+d_index][i][j][k] = gas[s][m][i][j][k];
          }
        }
      }
    }
  }

  d_index = mags_b1_p - mags_b1; // index difference between mags and mags_p
  #pragma omp parallel for collapse(4) schedule(static)
  for (int m = mags_b1; m <= mags_b3; m++) {
    for (int i = isg ; i <= ieg ; i++) {
      for (int j = jsg; j <= jeg ; j++) {
        for (int k = ksg; k <= keg ; k++) {
          gem[m+d_index][i][j][k] = gem[m][i][j][k];
        }
      }
    }
  }

  d_index = mag_bi_p - mag_bi; // index difference between mag and mag_p
  #pragma omp parallel for collapse(4) schedule(static)
  for (int m = mag_bi; m <= mag_bk; m++) {
    for (int i = isg ; i <= ieg +1 ; i++) {
      for (int j = jsg; j <= jeg +1 ; j++) {
        for (int k = ksg; k <= keg +1; k++) {
          gem[m+d_index][i][j][k] = gem[m][i][j][k];
        }
      }
    }
  }

  // update dt0
  dt0 = dt;
}

/**
 * @brief reconstruction using PDMU7 in i direction
 *
 * @param field  The 4D array of field values
 * @param m0     The starting index of the field
 * @param m1     The ending index of the field
 * @param ml     The starting index of the left field
 * @param mr     The starting index of the right field
 */
void reconstruct_3dv_i(double ****field, int pdm_index, int m0, int m1, int ml, int mr) {
  #pragma omp parallel for collapse(4) schedule(static)
  for (int m = m0; m <= m1; m++) {
    for (int i = is; i <= ie + 1; i++) {
      for (int j = jsg; j <= jeg; j++) {
        for (int k = ksg; k <= keg; k++) {
          int d_l = ml - m0;
          int d_r = mr - m0;
          double fh0[7]; // value of half step for interpolation
          double fp0[7]; // value of this step beginning for PDMU limiter
          double wgtPDMU[7]; // weight for PDMU7
          int offset0;
          offset0 = 0;
          wgtPDMU[0] = rec[w1l_first_direction][i][j][k];
          wgtPDMU[1] = rec[w2l_first_direction][i][j][k];
          wgtPDMU[2] = rec[w3l_first_direction][i][j][k];
          wgtPDMU[3] = rec[w4l_first_direction][i][j][k];
          wgtPDMU[4] = rec[w5l_first_direction][i][j][k];
          wgtPDMU[5] = rec[w6l_first_direction][i][j][k];
          wgtPDMU[6] = rec[w7l_first_direction][i][j][k];
          for (int offset = -4; offset <= 2; offset++) {
            fh0[offset0] = field[m][i+offset][j][k];
            fp0[offset0] = field[m+pdm_index][i+offset][j][k];
            offset0++;
          }
          // PDMU7(gas(i-4:i+2,j,k,m),vol(i-4:i+2,j,k,center),vol(i,j,k,dir))
          field[m + d_l][i][j][k] = PDMU7(fh0, fp0, wgtPDMU, geo[GF_face_idir][i-1][j][k], geo[GF_face_idir][i][j][k],
                                          geo[GF_vol_idir][i-1][j][k],dx1[i-1][j][k]);
          offset0 = 0;
          wgtPDMU[6] = rec[w1r_first_direction][i][j][k];
          wgtPDMU[5] = rec[w2r_first_direction][i][j][k];
          wgtPDMU[4] = rec[w3r_first_direction][i][j][k];
          wgtPDMU[3] = rec[w4r_first_direction][i][j][k];
          wgtPDMU[2] = rec[w5r_first_direction][i][j][k];
          wgtPDMU[1] = rec[w6r_first_direction][i][j][k];
          wgtPDMU[0] = rec[w7r_first_direction][i][j][k];
          for (int offset = +3; offset >= -3; offset--) {
            fh0[offset0] = field[m][i+offset][j][k];
            fp0[offset0] = field[m+pdm_index][i+offset][j][k];
            offset0++;
          }
          // PDMU7(gas(i+3:i-3:-1,j,k,m),vol(i+3:i-3:-1,j,k,center),vol(i,j,k,dir))
          field[m + d_r][i][j][k] = PDMU7(fh0, fp0, wgtPDMU, geo[GF_face_idir][i+1][j][k], geo[GF_face_idir][i][j][k],
                                          geo[GF_vol_idir][i][j][k],dx1[i][j][k]);
        }
      }
    }
  }
}
void reconstruct_3dv_j(double ****field, int pdm_index, int m0, int m1, int ml, int mr) {
  #pragma omp parallel for collapse(4) schedule(static)
  for (int m = m0; m <= m1; m++) {
    for (int i = isg; i <= ieg; i++) {
      for (int j = js; j <= je + 1; j++) {
        for (int k = ksg; k <= keg; k++) {
          int d_l = ml - m0;
          int d_r = mr - m0;
          double fh0[7]; // value of half step for interpolation
          double fp0[7]; // value of this step beginning for PDMU limiter
          double wgtPDMU[7]; // weight for PDMU7
          int offset0;
          offset0 = 0;
          wgtPDMU[0] = rec[w1l_second_direction][i][j][k];
          wgtPDMU[1] = rec[w2l_second_direction][i][j][k];
          wgtPDMU[2] = rec[w3l_second_direction][i][j][k];
          wgtPDMU[3] = rec[w4l_second_direction][i][j][k];
          wgtPDMU[4] = rec[w5l_second_direction][i][j][k];
          wgtPDMU[5] = rec[w6l_second_direction][i][j][k];
          wgtPDMU[6] = rec[w7l_second_direction][i][j][k];
          for (int offset = -4; offset <= 2; offset++) {
            fh0[offset0] = field[m][i][j+offset][k];
            fp0[offset0] = field[m+pdm_index][i][j+offset][k];
            offset0++;
          }
          // PDMU7(gas(i,j-4:j+2,k,m),vol(i,j-4:j+2,k,center),vol(i,j,k,dir))
          field[m + d_l][i][j][k] = PDMU7(fh0, fp0, wgtPDMU, geo[GF_face_jdir][i][j-1][k], geo[GF_face_jdir][i][j][k],
                                          geo[GF_vol_jdir][i][j-1][k],dx2[i][j-1][k]);
          offset0 = 0;
          wgtPDMU[6] = rec[w1r_second_direction][i][j][k];
          wgtPDMU[5] = rec[w2r_second_direction][i][j][k];
          wgtPDMU[4] = rec[w3r_second_direction][i][j][k];
          wgtPDMU[3] = rec[w4r_second_direction][i][j][k];
          wgtPDMU[2] = rec[w5r_second_direction][i][j][k];
          wgtPDMU[1] = rec[w6r_second_direction][i][j][k];
          wgtPDMU[0] = rec[w7r_second_direction][i][j][k];
          for (int offset = +3; offset >= -3; offset--) {
            fh0[offset0] = field[m][i][j+offset][k];
            fp0[offset0] = field[m+pdm_index][i][j+offset][k];
            offset0++;
          }
          // PDMU7(gas(i,j+3:j-3:-1,k,m),vol(i,j+3:j-3:-1,k,center),vol(i,j,k,dir))
          field[m + d_r][i][j][k] = PDMU7(fh0, fp0, wgtPDMU, geo[GF_face_jdir][i][j+1][k], geo[GF_face_jdir][i][j][k],
                                          geo[GF_vol_jdir][i][j][k],dx2[i][j][k]);
        }
      }
    }
  }
}
void reconstruct_3dv_k(double ****field, int pdm_index, int m0, int m1, int ml, int mr) {
  #pragma omp parallel for collapse(4) schedule(static)
  for (int m = m0; m <= m1; m++) {
    for (int i = isg; i <= ieg; i++) {
      for (int j = jsg; j <= jeg; j++) {
        for (int k = ks; k <= ke + 1; k++) {
          int d_l = ml - m0;
          int d_r = mr - m0;
          double fh0[7]; // value of half step for interpolation
          double fp0[7]; // value of this step beginning for PDMU limiter
          double wgtPDMU[7]; // weight for PDMU7
          int offset0;
          offset0 = 0;
          wgtPDMU[0] = rec[w1l_third_direction][i][j][k];
          wgtPDMU[1] = rec[w2l_third_direction][i][j][k];
          wgtPDMU[2] = rec[w3l_third_direction][i][j][k];
          wgtPDMU[3] = rec[w4l_third_direction][i][j][k];
          wgtPDMU[4] = rec[w5l_third_direction][i][j][k];
          wgtPDMU[5] = rec[w6l_third_direction][i][j][k];
          wgtPDMU[6] = rec[w7l_third_direction][i][j][k];
          for (int offset = -4; offset <= 2; offset++) {
            fh0[offset0] = field[m][i][j][k+offset];
            fp0[offset0] = field[m+pdm_index][i][j][k+offset];
            offset0++;
          }
          // PDMU7(gas(i,j,k-4:k+2,m),vol(i,j,k-4:k+2,center),vol(i,j,k,dir))
          field[m + d_l][i][j][k] = PDMU7(fh0, fp0, wgtPDMU, geo[GF_face_kdir][i][j][k-1], geo[GF_face_kdir][i][j][k],
                                          geo[GF_vol_kdir][i][j][k-1],dx3[i][j][k-1]);
          offset0 = 0;
          wgtPDMU[6] = rec[w1r_third_direction][i][j][k];
          wgtPDMU[5] = rec[w2r_third_direction][i][j][k];
          wgtPDMU[4] = rec[w3r_third_direction][i][j][k];
          wgtPDMU[3] = rec[w4r_third_direction][i][j][k];
          wgtPDMU[2] = rec[w5r_third_direction][i][j][k];
          wgtPDMU[1] = rec[w6r_third_direction][i][j][k];
          wgtPDMU[0] = rec[w7r_third_direction][i][j][k];
          for (int offset = +3; offset >= -3; offset--) {
            fh0[offset0] = field[m][i][j][k+offset];
            fp0[offset0] = field[m+pdm_index][i][j][k+offset];
            offset0++;
          }
          // PDMU7(gas(i,j,k+3:k-3:-1,m),vol(i,j,k+3:k-3:-1,center),vol(i,j,k,dir))
          field[m + d_r][i][j][k] = PDMU7(fh0, fp0, wgtPDMU, geo[GF_face_kdir][i][j][k+1], geo[GF_face_kdir][i][j][k],
                                          geo[GF_vol_kdir][i][j][k],dx3[i][j][k]);
        }
      }
    }
  }
}

void reconstruct_3dv_gem(dir_t dir) {
  switch (dir) {
    // TODO: magsl_b*, magsr_b* are intermediate variables
    case i_dir:
      reconstruct_3dv_i(gem, pdm_index_gem, mags_b1_h, mags_b3_h, magsl_b1, magsr_b1);
      break;
    case j_dir:
      reconstruct_3dv_j(gem, pdm_index_gem, mags_b1_h, mags_b3_h, magsl_b1, magsr_b1);
      break;
    case k_dir:
      reconstruct_3dv_k(gem, pdm_index_gem, mags_b1_h, mags_b3_h, magsl_b1, magsr_b1);
      break;
    default:
      break;
  }
}
void reconstruct_3dv_gas(int dir) {
  // TODO: NS or NS1 ?
  for (int s = 0; s < NS; s++) {
    switch (dir) {
      // TODO: gasl_*, gasr_* are intermediate variables
      case i_dir:
        reconstruct_3dv_i(gas[s], pdm_index_gas, gas_rho_h, gas_p_S_h, gasl_rho, gasr_rho);
        break;
      case j_dir:
        reconstruct_3dv_j(gas[s], pdm_index_gas, gas_rho_h, gas_p_S_h, gasl_rho, gasr_rho);
        break;
      case k_dir:
        reconstruct_3dv_k(gas[s], pdm_index_gas, gas_rho_h, gas_p_S_h, gasl_rho, gasr_rho);
        break;
      default:
        break;
    }
  }
}
void center2corner_ij(int pdm_index, int m0, int m_interp, int m_interp_l, int m_interp_r, int m_ava){
  // i ----> j
  // ---------------i---------------
  for (int s = 0; s < NS; s++) {
    reconstruct_3dv_i(gas[s], pdm_index, m0, m0, m_interp_l, m_interp_r);
  }
  #pragma omp parallel for collapse(4) schedule(static)
  for (int s = 0; s < NS; s++) {
    for (int i = is; i <= ie + 1; i++) {
      for (int j = jsg; j <= jeg; j++) {
        for (int k = ks; k <= ke; k++) {
          gas[s][m_interp][i][j][k] = 0.5 * (gas[s][m_interp_l][i][j][k] + gas[s][m_interp_r][i][j][k]);
        }
      }
    }
  }
  // ---------------j---------------
  for (int s = 0; s < NS; s++) {
    reconstruct_3dv_j(gas[s], 0, m_interp, m_interp, m_interp_l, m_interp_r);
  }
  #pragma omp parallel for collapse(4) schedule(static)
  for (int s = 0; s < NS; s++) {
    for (int i = is; i <= ie+1; i++){
      for (int j = js; j <= je+1; j++){
        for (int k = ks; k <= ke; k++) {
          gas[s][m_ava][i][j][k] = 0.5 * (gas[s][m_interp_l][i][j][k] + gas[s][m_interp_r][i][j][k]);
        }
      }
    }
  }
}
void center2corner_jk(int pdm_index, int m0, int m_interp, int m_interp_l, int m_interp_r, int m_ava){
  // j ----> k
  // ---------------j---------------
  for (int s = 0; s < NS; s++) {
    reconstruct_3dv_j(gas[s], pdm_index, m0, m0, m_interp_l, m_interp_r);
  }
  #pragma omp parallel for collapse(4) schedule(static)
  for (int s = 0; s < NS; s++) {
    for (int i = is; i <= ie; i++){
      for (int j = js; j <= je+1; j++){
        for (int k = ksg; k <= keg; k++) {
          gas[s][m_interp][i][j][k] = 0.5 * (gas[s][m_interp_l][i][j][k] + gas[s][m_interp_r][i][j][k]);
        }
      }
    }
  }
  // ---------------k---------------
  for (int s = 0; s < NS; s++) {
    reconstruct_3dv_k(gas[s], 0, m_interp, m_interp, m_interp_l, m_interp_r);
  }
  #pragma omp parallel for collapse(4) schedule(static)
  for (int s = 0; s < NS; s++) {
    for (int i = is; i <= ie; i++){
      for (int j = js; j <= je+1; j++){
        for (int k = ks; k <= ke+1; k++) {
          gas[s][m_ava][i][j][k] = 0.5 * (gas[s][m_interp_l][i][j][k] + gas[s][m_interp_r][i][j][k]);
        }
      }
    }
  }
}
void center2corner_ki(int pdm_index, int m0, int m_interp, int m_interp_l, int m_interp_r, int m_ava){
  // k ----> i
  // ---------------k---------------
  for (int s = 0; s < NS; s++) {
    reconstruct_3dv_k(gas[s], pdm_index, m0, m0, m_interp_l, m_interp_r);
  }
  #pragma omp parallel for collapse(4) schedule(static)
  for (int s = 0; s < NS; s++) {
    for (int i = isg; i <= ieg; i++){
      for (int j = js; j <= je; j++){
        for (int k = ks; k <= ke+1; k++) {
          gas[s][m_interp][i][j][k] = 0.5 * (gas[s][m_interp_l][i][j][k] + gas[s][m_interp_r][i][j][k]);
        }
      }
    }
  }
  // ---------------i---------------
  for (int s = 0; s < NS; s++) {
    reconstruct_3dv_i(gas[s], 0, m_interp, m_interp, m_interp_l, m_interp_r);
  }
  #pragma omp parallel for collapse(4) schedule(static)
  for (int s = 0; s < NS; s++) {
    for (int i = is; i <= ie+1; i++){
      for (int j = js; j <= je; j++){
        for (int k = ks; k <= ke+1; k++) {
          gas[s][m_ava][i][j][k] = 0.5 * (gas[s][m_interp_l][i][j][k] + gas[s][m_interp_r][i][j][k]);
        }
      }
    }
  }
}
void get_e_fields_i() {

  center2corner_jk(pdm_index_gas, gas_v2_h, v2_interp, v2_interp_l, v2_interp_r, v2_avg);
  center2corner_jk(pdm_index_gas, gas_v3_h, v3_interp, v3_interp_l, v3_interp_r, v3_avg);
  reconstruct_3dv_k(gem, pdm_index_gem, mag_bj_h, mag_bj_h, bjl, bjr);
  reconstruct_3dv_j(gem, pdm_index_gem, mag_bk_h, mag_bk_h, bkl, bkr);

  // TODO: it's possible to just use a variable instead of an array for
  // those intermediate variables
  // gas[NS1-1][v2_interp][i][j][k] -> v2_interp

  int s = NS1 - 1;
  #pragma omp parallel for collapse(3) schedule(static)
  for (int i = is; i <= ie; i++) {
    for (int j = js; j <= je + 1; j++) {
      for (int k = ks; k <= ke + 1; k++) {
        // TODO: bi_avg, bj_avg, bk_avg are intermediate variables
        gem[bj_avg][i][j][k] = 0.5 * (gem[bjl][i][j][k] + gem[bjr][i][j][k]) +
                               gem[B0AedgeI_b2][i][j][k];
        gem[bk_avg][i][j][k] = 0.5 * (gem[bkl][i][j][k] + gem[bkr][i][j][k]) +
                               gem[B0AedgeI_b3][i][j][k];
        // TODO: efield_ei is intermediate variable
        gem[efield_ei][i][j][k] =
            -(gas[s][v2_avg][i][j][k] * gem[bk_avg][i][j][k] -
              gas[s][v3_avg][i][j][k] * gem[bj_avg][i][j][k]);
        // TODO: dvzz is intermediate variable
        gem[dvzz][i][j][k] = sqrt(
            (gem[bj_avg][i][j][k] * gem[bj_avg][i][j][k] +
             gem[bk_avg][i][j][k] * gem[bk_avg][i][j][k]) /
            (0.25 *
             (gas[s][gas_rho_h][i][j - 1][k - 1] + gas[s][gas_rho_h][i][j - 1][k] +
              gas[s][gas_rho_h][i][j][k - 1] + gas[s][gas_rho_h][i][j][k])));

        gem[dvzz][i][j][k] =
            gem[dvzz][i][j][k] * CA /
                sqrt(gem[dvzz][i][j][k] * gem[dvzz][i][j][k] + CA * CA) +
            sqrt(gas[s][v2_avg][i][j][k] * gas[s][v2_avg][i][j][k] +
                 gas[s][v3_avg][i][j][k] * gas[s][v3_avg][i][j][k]);

        gem[efield_ei][i][j][k] += 0.5 * gem[dvzz][i][j][k] *
                                   (gem[bjl][i][j][k] + gem[bkr][i][j][k] -
                                    gem[bjr][i][j][k] - gem[bkl][i][j][k]);
      }
    }
  }
}
void get_e_fields_j() {

  center2corner_ki(pdm_index_gas, gas_v3_h, v3_interp, v3_interp_l, v3_interp_r, v3_avg);
  center2corner_ki(pdm_index_gas, gas_v1_h, v1_interp, v1_interp_l, v1_interp_r, v1_avg);
  reconstruct_3dv_i(gem, pdm_index_gem, mag_bk_h, mag_bk_h, bkl, bkr);
  reconstruct_3dv_k(gem, pdm_index_gem, mag_bi_h, mag_bi_h, bil, bir);

  int s = NS1 - 1;
  #pragma omp parallel for collapse(3) schedule(static)
  for (int i = is; i <= ie + 1; i++) {
    for (int j = js; j <= je; j++) {
      for (int k = ks; k <= ke + 1; k++) {
        gem[bk_avg][i][j][k] = 0.5 * (gem[bkl][i][j][k] + gem[bkr][i][j][k]) +
            gem[B0AedgeJ_b3][i][j][k];
        gem[bi_avg][i][j][k] = 0.5 * (gem[bil][i][j][k] + gem[bir][i][j][k]) +
            gem[B0AedgeJ_b1][i][j][k];
        gem[efield_ej][i][j][k] =
            -(gas[s][v3_avg][i][j][k] * gem[bi_avg][i][j][k] -
                gas[s][v1_avg][i][j][k] * gem[bk_avg][i][j][k]);
        gem[dvzz][i][j][k] = sqrt(
            (gem[bk_avg][i][j][k] * gem[bk_avg][i][j][k] +
                gem[bi_avg][i][j][k] * gem[bi_avg][i][j][k]) /
                (0.25 *
                    (gas[s][gas_rho_h][i - 1][j][k - 1] + gas[s][gas_rho_h][i - 1][j][k] +
                        gas[s][gas_rho_h][i][j][k - 1] + gas[s][gas_rho_h][i][j][k])));

        gem[dvzz][i][j][k] =
            gem[dvzz][i][j][k] * CA /
                sqrt(gem[dvzz][i][j][k] * gem[dvzz][i][j][k] + CA * CA) +
                sqrt(gas[s][v1_avg][i][j][k] * gas[s][v1_avg][i][j][k] +
                    gas[s][v3_avg][i][j][k] * gas[s][v3_avg][i][j][k]);

        gem[efield_ej][i][j][k] += 0.5 * gem[dvzz][i][j][k] *
            (gem[bkl][i][j][k] + gem[bir][i][j][k] -
                gem[bkr][i][j][k] - gem[bil][i][j][k]);
      }
    }
  }
}
void get_e_fields_k() {

  center2corner_ij(pdm_index_gas, gas_v1_h, v1_interp, v1_interp_l, v1_interp_r, v1_avg);
  center2corner_ij(pdm_index_gas, gas_v2_h, v2_interp, v2_interp_l, v2_interp_r, v2_avg);
  reconstruct_3dv_j(gem, pdm_index_gem, mag_bi_h, mag_bi_h, bil, bir);
  reconstruct_3dv_i(gem, pdm_index_gem, mag_bj_h, mag_bj_h, bjl, bjr);

  int s = NS1 - 1;
  #pragma omp parallel for collapse(3) schedule(static)
  for (int i = is; i <= ie + 1; i++) {
    for (int j = js; j <= je + 1; j++) {
      for (int k = ks; k <= ke; k++) {
        gem[bi_avg][i][j][k] = 0.5 * (gem[bil][i][j][k] + gem[bir][i][j][k]) +
            gem[B0AedgeK_b1][i][j][k];
        gem[bj_avg][i][j][k] = 0.5 * (gem[bjl][i][j][k] + gem[bjr][i][j][k]) +
            gem[B0AedgeK_b2][i][j][k];
        gem[efield_ek][i][j][k] =
            -(gas[s][v1_avg][i][j][k] * gem[bj_avg][i][j][k] -
                gas[s][v2_avg][i][j][k] * gem[bi_avg][i][j][k]);
        gem[dvzz][i][j][k] = sqrt(
            (gem[bi_avg][i][j][k] * gem[bi_avg][i][j][k] +
                gem[bj_avg][i][j][k] * gem[bj_avg][i][j][k]) /
                (0.25 *
                    (gas[s][gas_rho_h][i - 1][j - 1][k] + gas[s][gas_rho_h][i - 1][j][k] +
                        gas[s][gas_rho_h][i][j - 1][k] + gas[s][gas_rho_h][i][j][k])));

        gem[dvzz][i][j][k] =
            gem[dvzz][i][j][k] * CA /
                sqrt(gem[dvzz][i][j][k] * gem[dvzz][i][j][k] + CA * CA) +
                sqrt(gas[s][v2_avg][i][j][k] * gas[s][v2_avg][i][j][k] +
                    gas[s][v1_avg][i][j][k] * gas[s][v1_avg][i][j][k]);

        gem[efield_ek][i][j][k] += 0.5 * gem[dvzz][i][j][k] *
            (gem[bil][i][j][k] + gem[bjr][i][j][k] -
                gem[bir][i][j][k] - gem[bjl][i][j][k]);
      }
    }
  }
}
void get_e_fields() {
  // This code calculates the Electric field components on cell edges - Ei, Ej
  // and Ek using high order constrained transport (Yee grid) method, details
  // described in Lyon et al., [2004]. The electric field is estimated as:
  //                      E = - v_avg x B_avg + eta*J
  // where v_avg is an average velocity at cell edges and B_upwind is the
  // upwinded magnetic field at cell edges chosen based on the average velocity.
  // Note that the eta*J term here is not usual resistive MHD term; rather it's
  // only turned on when the limiter detects a discontinuity in the B fields.
  // The eta is set to be the local fast speed averaged around the cell edges.
  // This term is important in the regions where Alfven waves are present.

  // ALGORITHM (Use Ek as an example):
  // STEP 1: interpolate cell-centered vx, vy to cell-edge (vx_avg, vy_avg);
  // STEP 2: reconstruct Bi in y-direction to get bi_left and bi_right at cell
  // edges
  //         reconstruct Bj in x-direction to get bj_left and bj_right at cell
  //         edges
  // STEP 3: calculate B_avg (simplified)
  // STEP 4: compute the electric field as Ek = - v_avg x B_avg + eta*j
  // STEP 5: the reconstructed Bi and Bi form a current at cell edges:
  //         J = bi_left + bj_right - bi_right - bj_left; the diffusion term is
  //         added as v_diffusive * J, v_diffusive is chosen as the average
  //         local fast mode speed. Note that this diffusive term contains
  //         both numerical and Alflven resistivity.

  get_e_fields_i();
  get_e_fields_j();
  get_e_fields_k();
}

void get_fluid_flux_rusanov_i() {
  #pragma omp parallel for collapse(4) schedule(static)
  for (int s = 0; s < NS1; s++) {
    for (int i = is; i <= ie + 1; i++) {
      for (int j = js; j <= je; j++) {
        for (int k = ks; k <= ke; k++) {
          double Vsound_l,Vsound_r, Vsound_entropy_l, Vsound_entropy_r, VF, VF_entropy,
              rhov1_l,rhov1_r, rhov2_l, rhov2_r, rhov3_l, rhov3_r, eng_l, eng_r, S_l, S_r;
          double B_tot_sq_l, B_tot_sq_r, Va_l, Va_r, Va_eff_l, Va_eff_r;
          Vsound_l = sqrt(gamma_val * gas[s][gasl_p][i][j][k] /
                          gas[s][gasl_rho][i][j][k]);
          Vsound_r = sqrt(gamma_val * gas[s][gasr_p][i][j][k] /
                          gas[s][gasr_rho][i][j][k]);
          Vsound_entropy_l = sqrt(gamma_val * gas[s][gasl_p_S][i][j][k] /
                                  gas[s][gasl_rho][i][j][k]);
          Vsound_entropy_r = sqrt(gamma_val * gas[s][gasr_p_S][i][j][k] /
                                  gas[s][gasr_rho][i][j][k]);

          B_tot_sq_l =
              (gem[magsl_b1][i][j][k] + gem[B0AfaceI_b1][i][j][k]) *
                  (gem[magsl_b1][i][j][k] + gem[B0AfaceI_b1][i][j][k]) +
              (gem[magsl_b2][i][j][k] + gem[B0AfaceI_b2][i][j][k]) *
                  (gem[magsl_b2][i][j][k] + gem[B0AfaceI_b2][i][j][k]) +
              (gem[magsl_b3][i][j][k] + gem[B0AfaceI_b3][i][j][k]) *
                  (gem[magsl_b3][i][j][k] + gem[B0AfaceI_b3][i][j][k]);

          B_tot_sq_r =
              (gem[magsr_b1][i][j][k] + gem[B0AfaceI_b1][i][j][k]) *
                  (gem[magsr_b1][i][j][k] + gem[B0AfaceI_b1][i][j][k]) +
              (gem[magsr_b2][i][j][k] + gem[B0AfaceI_b2][i][j][k]) *
                  (gem[magsr_b2][i][j][k] + gem[B0AfaceI_b2][i][j][k]) +
              (gem[magsr_b3][i][j][k] + gem[B0AfaceI_b3][i][j][k]) *
                  (gem[magsr_b3][i][j][k] + gem[B0AfaceI_b3][i][j][k]);

          Va_l = sqrt(B_tot_sq_l / gas[s][gasl_rho][i][j][k]);
          Va_r = sqrt(B_tot_sq_r / gas[s][gasr_rho][i][j][k]);

          Va_eff_l = Va_l * CA / (sqrt(Va_l * Va_l + CA * CA));
          Va_eff_r = Va_r * CA / (sqrt(Va_r * Va_r + CA * CA));

          VF = fmax(fabs(gas[s][gasl_v1][i][j][k]) +
                        sqrt(Vsound_l * Vsound_l + Va_eff_l * Va_eff_l),
                    fabs(gas[s][gasr_v1][i][j][k]) +
                        sqrt(Vsound_r * Vsound_r + Va_eff_r * Va_eff_r));
          VF_entropy = fmax(fabs(gas[s][gasl_v1][i][j][k]) +
                        sqrt(Vsound_entropy_l * Vsound_entropy_l + Va_eff_l * Va_eff_l),
                    fabs(gas[s][gasr_v1][i][j][k]) +
                        sqrt(Vsound_entropy_r * Vsound_entropy_r + Va_eff_r * Va_eff_r));

          rhov1_l = gas[s][gasl_rho][i][j][k] * gas[s][gasl_v1][i][j][k];
          rhov1_r = gas[s][gasr_rho][i][j][k] * gas[s][gasr_v1][i][j][k];

          rhov2_l = gas[s][gasl_rho][i][j][k] * gas[s][gasl_v2][i][j][k];
          rhov2_r = gas[s][gasr_rho][i][j][k] * gas[s][gasr_v2][i][j][k];

          rhov3_l = gas[s][gasl_rho][i][j][k] * gas[s][gasl_v3][i][j][k];
          rhov3_r = gas[s][gasr_rho][i][j][k] * gas[s][gasr_v3][i][j][k];

          eng_l = 0.5 * gas[s][gasl_rho][i][j][k] *
                      (pow(gas[s][gasl_v1][i][j][k], 2) +
                       pow(gas[s][gasl_v2][i][j][k], 2) +
                       pow(gas[s][gasl_v3][i][j][k], 2)) +
                  gas[s][gasl_p][i][j][k] / (gamma_val - 1.0);
          eng_r = 0.5 * gas[s][gasr_rho][i][j][k] *
                      (pow(gas[s][gasr_v1][i][j][k], 2) +
                       pow(gas[s][gasr_v2][i][j][k], 2) +
                       pow(gas[s][gasr_v3][i][j][k], 2)) +
                  gas[s][gasr_p][i][j][k] / (gamma_val - 1.0);

          S_l = gas[s][gasl_p_S][i][j][k] / pow(gas[s][gasl_rho][i][j][k], gamma_val-1.0);
          S_r = gas[s][gasr_p_S][i][j][k] / pow(gas[s][gasr_rho][i][j][k], gamma_val-1.0);

          // TODO: gascf_* are intermediate variables
          gas[s][gascf_rho][i][j][k] =
              0.5 *
              (rhov1_l + rhov1_r -
               VF * (gas[s][gasr_rho][i][j][k] - gas[s][gasl_rho][i][j][k]));
          gas[s][gascf_rhov1][i][j][k] =
              0.5 *
              (rhov1_l * gas[s][gasl_v1][i][j][k] + gas[s][gasl_p][i][j][k] +
               rhov1_r * gas[s][gasr_v1][i][j][k] + gas[s][gasr_p][i][j][k] -
               VF * (rhov1_r - rhov1_l));
          gas[s][gascf_rhov2][i][j][k] =
              0.5 *
              (rhov2_l * gas[s][gasl_v1][i][j][k] +
               rhov2_r * gas[s][gasr_v1][i][j][k] - VF * (rhov2_r - rhov2_l));
          gas[s][gascf_rhov3][i][j][k] =
              0.5 *
              (rhov3_l * gas[s][gasl_v1][i][j][k] +
               rhov3_r * gas[s][gasr_v1][i][j][k] - VF * (rhov3_r - rhov3_l));
          gas[s][gascf_eng][i][j][k] =
              0.5 * (eng_l * gas[s][gasl_v1][i][j][k] +
                     gas[s][gasl_p][i][j][k] * gas[s][gasl_v1][i][j][k] +
                     eng_r * gas[s][gasr_v1][i][j][k] +
                     gas[s][gasr_p][i][j][k] * gas[s][gasr_v1][i][j][k] -
                     VF * (eng_r - eng_l));
          gas[s][gascf_entropy][i][j][k] =
              0.5 * (S_l * gas[s][gasl_v1][i][j][k] +
                     S_r * gas[s][gasr_v1][i][j][k] -
                     VF_entropy * (S_r - S_l) );
        }
      }
    }
  }
}

void get_fluid_flux_rusanov_j() {
  #pragma omp parallel for collapse(4) schedule(static)
  for (int s = 0; s < NS1; s++) {
    for (int i = is; i <= ie; i++) {
      for (int j = js; j <= je + 1; j++) {
        for (int k = ks; k <= ke; k++) {
          double Vsound_l,Vsound_r, Vsound_entropy_l, Vsound_entropy_r, VF, VF_entropy,
              rhov1_l,rhov1_r, rhov2_l, rhov2_r, rhov3_l, rhov3_r, eng_l, eng_r, S_l, S_r;

          double B_tot_sq_l, B_tot_sq_r, Va_l, Va_r, Va_eff_l, Va_eff_r;
          Vsound_l = sqrt(gamma_val * gas[s][gasl_p][i][j][k] /
                          gas[s][gasl_rho][i][j][k]);
          Vsound_r = sqrt(gamma_val * gas[s][gasr_p][i][j][k] /
                          gas[s][gasr_rho][i][j][k]);
          Vsound_entropy_l = sqrt(gamma_val * gas[s][gasl_p_S][i][j][k] /
                                  gas[s][gasl_rho][i][j][k]);
          Vsound_entropy_r = sqrt(gamma_val * gas[s][gasr_p_S][i][j][k] /
                                  gas[s][gasr_rho][i][j][k]);

          B_tot_sq_l =
              (gem[magsl_b1][i][j][k] + gem[B0AfaceJ_b1][i][j][k]) *
                  (gem[magsl_b1][i][j][k] + gem[B0AfaceJ_b1][i][j][k]) +
              (gem[magsl_b2][i][j][k] + gem[B0AfaceJ_b2][i][j][k]) *
                  (gem[magsl_b2][i][j][k] + gem[B0AfaceJ_b2][i][j][k]) +
              (gem[magsl_b3][i][j][k] + gem[B0AfaceJ_b3][i][j][k]) *
                  (gem[magsl_b3][i][j][k] + gem[B0AfaceJ_b3][i][j][k]);

          B_tot_sq_r =
              (gem[magsr_b1][i][j][k] + gem[B0AfaceJ_b1][i][j][k]) *
                  (gem[magsr_b1][i][j][k] + gem[B0AfaceJ_b1][i][j][k]) +
              (gem[magsr_b2][i][j][k] + gem[B0AfaceJ_b2][i][j][k]) *
                  (gem[magsr_b2][i][j][k] + gem[B0AfaceJ_b2][i][j][k]) +
              (gem[magsr_b3][i][j][k] + gem[B0AfaceJ_b3][i][j][k]) *
                  (gem[magsr_b3][i][j][k] + gem[B0AfaceJ_b3][i][j][k]);

          Va_l = sqrt(B_tot_sq_l / gas[s][gasl_rho][i][j][k]);
          Va_r = sqrt(B_tot_sq_r / gas[s][gasr_rho][i][j][k]);

          Va_eff_l = Va_l * CA / (sqrt(Va_l * Va_l + CA * CA));
          Va_eff_r = Va_r * CA / (sqrt(Va_r * Va_r + CA * CA));

          VF = fmax(fabs(gas[s][gasl_v2][i][j][k]) +
                        sqrt(Vsound_l * Vsound_l + Va_eff_l * Va_eff_l),
                    fabs(gas[s][gasr_v2][i][j][k]) +
                        sqrt(Vsound_r * Vsound_r + Va_eff_r * Va_eff_r));
          VF_entropy = fmax(fabs(gas[s][gasl_v2][i][j][k]) +
                        sqrt(Vsound_entropy_l * Vsound_entropy_l + Va_eff_l * Va_eff_l),
                    fabs(gas[s][gasr_v2][i][j][k]) +
                        sqrt(Vsound_entropy_r * Vsound_entropy_r + Va_eff_r * Va_eff_r));

          rhov1_l = gas[s][gasl_rho][i][j][k] * gas[s][gasl_v1][i][j][k];
          rhov1_r = gas[s][gasr_rho][i][j][k] * gas[s][gasr_v1][i][j][k];

          rhov2_l = gas[s][gasl_rho][i][j][k] * gas[s][gasl_v2][i][j][k];
          rhov2_r = gas[s][gasr_rho][i][j][k] * gas[s][gasr_v2][i][j][k];

          rhov3_l = gas[s][gasl_rho][i][j][k] * gas[s][gasl_v3][i][j][k];
          rhov3_r = gas[s][gasr_rho][i][j][k] * gas[s][gasr_v3][i][j][k];

          eng_l = 0.5 * gas[s][gasl_rho][i][j][k] *
                      (pow(gas[s][gasl_v1][i][j][k], 2) +
                       pow(gas[s][gasl_v2][i][j][k], 2) +
                       pow(gas[s][gasl_v3][i][j][k], 2)) +
                  gas[s][gasl_p][i][j][k] / (gamma_val - 1);
          eng_r = 0.5 * gas[s][gasr_rho][i][j][k] *
                      (pow(gas[s][gasr_v1][i][j][k], 2) +
                       pow(gas[s][gasr_v2][i][j][k], 2) +
                       pow(gas[s][gasr_v3][i][j][k], 2)) +
                  gas[s][gasr_p][i][j][k] / (gamma_val - 1);

          S_l = gas[s][gasl_p_S][i][j][k] / pow(gas[s][gasl_rho][i][j][k], gamma_val-1.0);
          S_r = gas[s][gasr_p_S][i][j][k] / pow(gas[s][gasr_rho][i][j][k], gamma_val-1.0);

          gas[s][gascf_rho][i][j][k] =
              0.5 *
              (rhov2_l + rhov2_r -
               VF * (gas[s][gasr_rho][i][j][k] - gas[s][gasl_rho][i][j][k]));
          gas[s][gascf_rhov1][i][j][k] =
              0.5 *
              (rhov1_l * gas[s][gasl_v2][i][j][k] +
               rhov1_r * gas[s][gasr_v2][i][j][k] - VF * (rhov1_r - rhov1_l));
          gas[s][gascf_rhov2][i][j][k] =
              0.5 *
              (rhov2_l * gas[s][gasl_v2][i][j][k] + gas[s][gasl_p][i][j][k] +
               rhov2_r * gas[s][gasr_v2][i][j][k] + gas[s][gasr_p][i][j][k] -
               VF * (rhov2_r - rhov2_l));
          gas[s][gascf_rhov3][i][j][k] =
              0.5 *
              (rhov3_l * gas[s][gasl_v2][i][j][k] +
               rhov3_r * gas[s][gasr_v2][i][j][k] - VF * (rhov3_r - rhov3_l));
          gas[s][gascf_eng][i][j][k] =
              0.5 * (eng_l * gas[s][gasl_v2][i][j][k] +
                     gas[s][gasl_p][i][j][k] * gas[s][gasl_v2][i][j][k] +
                     eng_r * gas[s][gasr_v2][i][j][k] +
                     gas[s][gasr_p][i][j][k] * gas[s][gasr_v2][i][j][k] -
                     VF * (eng_r - eng_l));
          gas[s][gascf_entropy][i][j][k] =
              0.5 * (S_l * gas[s][gasl_v2][i][j][k] +
                     S_r * gas[s][gasr_v2][i][j][k] -
                     VF_entropy * (S_r - S_l) );
        }
      }
    }
  }
}

void get_fluid_flux_rusanov_k() {
  #pragma omp parallel for collapse(4) schedule(static)
  for (int s = 0; s < NS1; s++) {
    for (int i = is; i <= ie; i++) {
      for (int j = js; j <= je; j++) {
        for (int k = ks; k <= ke + 1; k++) {
          double Vsound_l,Vsound_r, Vsound_entropy_l, Vsound_entropy_r, VF, VF_entropy,
              rhov1_l,rhov1_r, rhov2_l, rhov2_r, rhov3_l, rhov3_r, eng_l, eng_r, S_l, S_r;
          double B_tot_sq_l, B_tot_sq_r, Va_l, Va_r, Va_eff_l, Va_eff_r;
          Vsound_l = sqrt(gamma_val * gas[s][gasl_p][i][j][k] /
                          gas[s][gasl_rho][i][j][k]);
          Vsound_r = sqrt(gamma_val * gas[s][gasr_p][i][j][k] /
                          gas[s][gasr_rho][i][j][k]);
          Vsound_entropy_l = sqrt(gamma_val * gas[s][gasl_p_S][i][j][k] /
                                  gas[s][gasl_rho][i][j][k]);
          Vsound_entropy_r = sqrt(gamma_val * gas[s][gasr_p_S][i][j][k] /
                                  gas[s][gasr_rho][i][j][k]);

          B_tot_sq_l =
              (gem[magsl_b1][i][j][k] + gem[B0AfaceK_b1][i][j][k]) *
                  (gem[magsl_b1][i][j][k] + gem[B0AfaceK_b1][i][j][k]) +
              (gem[magsl_b2][i][j][k] + gem[B0AfaceK_b2][i][j][k]) *
                  (gem[magsl_b2][i][j][k] + gem[B0AfaceK_b2][i][j][k]) +
              (gem[magsl_b3][i][j][k] + gem[B0AfaceK_b3][i][j][k]) *
                  (gem[magsl_b3][i][j][k] + gem[B0AfaceK_b3][i][j][k]);

          B_tot_sq_r =
              (gem[magsr_b1][i][j][k] + gem[B0AfaceK_b1][i][j][k]) *
                  (gem[magsr_b1][i][j][k] + gem[B0AfaceK_b1][i][j][k]) +
              (gem[magsr_b2][i][j][k] + gem[B0AfaceK_b2][i][j][k]) *
                  (gem[magsr_b2][i][j][k] + gem[B0AfaceK_b2][i][j][k]) +
              (gem[magsr_b3][i][j][k] + gem[B0AfaceK_b3][i][j][k]) *
                  (gem[magsr_b3][i][j][k] + gem[B0AfaceK_b3][i][j][k]);

          Va_l = sqrt(B_tot_sq_l / gas[s][gasl_rho][i][j][k]);
          Va_r = sqrt(B_tot_sq_r / gas[s][gasr_rho][i][j][k]);

          Va_eff_l = Va_l * CA / (sqrt(Va_l * Va_l + CA * CA));
          Va_eff_r = Va_r * CA / (sqrt(Va_r * Va_r + CA * CA));

          VF = fmax(fabs(gas[s][gasl_v3][i][j][k]) +
                        sqrt(Vsound_l * Vsound_l + Va_eff_l * Va_eff_l),
                    fabs(gas[s][gasr_v3][i][j][k]) +
                        sqrt(Vsound_r * Vsound_r + Va_eff_r * Va_eff_r));
          VF_entropy = fmax(fabs(gas[s][gasl_v3][i][j][k]) +
                        sqrt(Vsound_entropy_l * Vsound_entropy_l + Va_eff_l * Va_eff_l),
                    fabs(gas[s][gasr_v3][i][j][k]) +
                        sqrt(Vsound_entropy_r * Vsound_entropy_r + Va_eff_r * Va_eff_r));

          rhov1_l = gas[s][gasl_rho][i][j][k] * gas[s][gasl_v1][i][j][k];
          rhov1_r = gas[s][gasr_rho][i][j][k] * gas[s][gasr_v1][i][j][k];

          rhov2_l = gas[s][gasl_rho][i][j][k] * gas[s][gasl_v2][i][j][k];
          rhov2_r = gas[s][gasr_rho][i][j][k] * gas[s][gasr_v2][i][j][k];

          rhov3_l = gas[s][gasl_rho][i][j][k] * gas[s][gasl_v3][i][j][k];
          rhov3_r = gas[s][gasr_rho][i][j][k] * gas[s][gasr_v3][i][j][k];

          eng_l = 0.5 * gas[s][gasl_rho][i][j][k] *
                      (pow(gas[s][gasl_v1][i][j][k], 2) +
                       pow(gas[s][gasl_v2][i][j][k], 2) +
                       pow(gas[s][gasl_v3][i][j][k], 2)) +
                  gas[s][gasl_p][i][j][k] / (gamma_val - 1);
          eng_r = 0.5 * gas[s][gasr_rho][i][j][k] *
                      (pow(gas[s][gasr_v1][i][j][k], 2) +
                       pow(gas[s][gasr_v2][i][j][k], 2) +
                       pow(gas[s][gasr_v3][i][j][k], 2)) +
                  gas[s][gasr_p][i][j][k] / (gamma_val - 1);

          S_l = gas[s][gasl_p_S][i][j][k] / pow(gas[s][gasl_rho][i][j][k], gamma_val-1.0);
          S_r = gas[s][gasr_p_S][i][j][k] / pow(gas[s][gasr_rho][i][j][k], gamma_val-1.0);

          gas[s][gascf_rho][i][j][k] =
              0.5 *
              (rhov3_l + rhov3_r -
               VF * (gas[s][gasr_rho][i][j][k] - gas[s][gasl_rho][i][j][k]));

          gas[s][gascf_rhov1][i][j][k] =
              0.5 *
              (rhov1_l * gas[s][gasl_v3][i][j][k] +
               rhov1_r * gas[s][gasr_v3][i][j][k] - VF * (rhov1_r - rhov1_l));
          gas[s][gascf_rhov2][i][j][k] =
              0.5 *
              (rhov2_l * gas[s][gasl_v3][i][j][k] +
               rhov2_r * gas[s][gasr_v3][i][j][k] - VF * (rhov2_r - rhov2_l));
          gas[s][gascf_rhov3][i][j][k] =
              0.5 *
              (rhov3_l * gas[s][gasl_v3][i][j][k] + gas[s][gasl_p][i][j][k] +
               rhov3_r * gas[s][gasr_v3][i][j][k] + gas[s][gasr_p][i][j][k] -
               VF * (rhov3_r - rhov3_l));
          gas[s][gascf_eng][i][j][k] =
              0.5 * (eng_l * gas[s][gasl_v3][i][j][k] +
                     gas[s][gasl_p][i][j][k] * gas[s][gasl_v3][i][j][k] +
                     eng_r * gas[s][gasr_v3][i][j][k] +
                     gas[s][gasr_p][i][j][k] * gas[s][gasr_v3][i][j][k] -
                     VF * (eng_r - eng_l));
          gas[s][gascf_entropy][i][j][k] =
              0.5 * (S_l * gas[s][gasl_v3][i][j][k] +
                     S_r * gas[s][gasr_v3][i][j][k] -
                     VF_entropy * (S_r - S_l) );
        }
      }
    }
  }
}

void get_fluid_flux_rusanov(dir_t dir) {
  switch (dir) {
    case i_dir:
      get_fluid_flux_rusanov_i();
      break;
    case j_dir:
      get_fluid_flux_rusanov_j();
      break;
    case k_dir:
      get_fluid_flux_rusanov_k();
      break;
    default:
      break;
  }
}

void get_magnetic_stress_rusanov_i() {
  double bsq_l, bsq_r, B_tot_sq_l, B_tot_sq_r, B_sq, Bstress_1_p, Bstress_2_p,
      Bstress_3_p, Bstress_1_n, Bstress_2_n, Bstress_3_n;
  int i, j, k;
  int s = NS1 - 1;
  for (i = is; i <= ie + 1; i++) {
    for (j = js; j <= je; j++) {
      for (k = ks; k <= ke; k++) {
        bsq_l = (gem[magsl_b1][i][j][k] * gem[magsl_b1][i][j][k] +
                 gem[magsl_b2][i][j][k] * gem[magsl_b2][i][j][k] +
                 gem[magsl_b3][i][j][k] * gem[magsl_b3][i][j][k]) +
                2 * (gem[magsl_b1][i][j][k] * gem[B0AfaceI_b1][i][j][k] +
                     gem[magsl_b2][i][j][k] * gem[B0AfaceI_b2][i][j][k] +
                     gem[magsl_b3][i][j][k] * gem[B0AfaceI_b3][i][j][k]);

        bsq_r = (gem[magsr_b1][i][j][k] * gem[magsr_b1][i][j][k] +
                 gem[magsr_b2][i][j][k] * gem[magsr_b2][i][j][k] +
                 gem[magsr_b3][i][j][k] * gem[magsr_b3][i][j][k]) +
                2 * (gem[magsr_b1][i][j][k] * gem[B0AfaceI_b1][i][j][k] +
                     gem[magsr_b2][i][j][k] * gem[B0AfaceI_b2][i][j][k] +
                     gem[magsr_b3][i][j][k] * gem[B0AfaceI_b3][i][j][k]);

        Bstress_1_p =
            0.5 * (0.5 * bsq_l -
                   (gem[magsl_b1][i][j][k] + gem[B0AfaceI_b1][i][j][k]) *
                       gem[mag_bi][i][j][k] -
                   gem[magsl_b1][i][j][k] * gem[B0AfaceI_bn][i][j][k]);
        Bstress_2_p =
            0.5 * (-(gem[magsl_b2][i][j][k] + gem[B0AfaceI_b2][i][j][k]) *
                       gem[mag_bi][i][j][k] -
                   gem[magsl_b2][i][j][k] * gem[B0AfaceI_bn][i][j][k]);
        Bstress_3_p =
            0.5 * (-(gem[magsl_b3][i][j][k] + gem[B0AfaceI_b3][i][j][k]) *
                       gem[mag_bi][i][j][k] -
                   gem[magsl_b3][i][j][k] * gem[B0AfaceI_bn][i][j][k]);

        Bstress_1_n =
            0.5 * (0.5 * bsq_r -
                   (gem[magsr_b1][i][j][k] + gem[B0AfaceI_b1][i][j][k]) *
                       gem[mag_bi][i][j][k] -
                   gem[magsr_b1][i][j][k] * gem[B0AfaceI_bn][i][j][k]);
        Bstress_2_n =
            0.5 * (-(gem[magsr_b2][i][j][k] + gem[B0AfaceI_b2][i][j][k]) *
                       gem[mag_bi][i][j][k] -
                   gem[magsr_b2][i][j][k] * gem[B0AfaceI_bn][i][j][k]);
        Bstress_3_n =
            0.5 * (-(gem[magsr_b3][i][j][k] + gem[B0AfaceI_b3][i][j][k]) *
                       gem[mag_bi][i][j][k] -
                   gem[magsr_b3][i][j][k] * gem[B0AfaceI_bn][i][j][k]);

        gem[magsf_x1dir][i][j][k] = Bstress_1_p + Bstress_1_n;
        gem[magsf_x2dir][i][j][k] = Bstress_2_p + Bstress_2_n;
        gem[magsf_x3dir][i][j][k] = Bstress_3_p + Bstress_3_n;

        B_tot_sq_l = (gem[magsl_b1][i][j][k] + gem[B0AfaceI_b1][i][j][k]) *
                         (gem[magsl_b1][i][j][k] + gem[B0AfaceI_b1][i][j][k]) +
                     (gem[magsl_b2][i][j][k] + gem[B0AfaceI_b2][i][j][k]) *
                         (gem[magsl_b2][i][j][k] + gem[B0AfaceI_b2][i][j][k]) +
                     (gem[magsl_b3][i][j][k] + gem[B0AfaceI_b3][i][j][k]) *
                         (gem[magsl_b3][i][j][k] + gem[B0AfaceI_b3][i][j][k]);

        B_tot_sq_r = (gem[magsr_b1][i][j][k] + gem[B0AfaceI_b1][i][j][k]) *
                         (gem[magsr_b1][i][j][k] + gem[B0AfaceI_b1][i][j][k]) +
                     (gem[magsr_b2][i][j][k] + gem[B0AfaceI_b2][i][j][k]) *
                         (gem[magsr_b2][i][j][k] + gem[B0AfaceI_b2][i][j][k]) +
                     (gem[magsr_b3][i][j][k] + gem[B0AfaceI_b3][i][j][k]) *
                         (gem[magsr_b3][i][j][k] + gem[B0AfaceI_b3][i][j][k]);

        B_sq = (B_tot_sq_l + B_tot_sq_r) / 4 / CA;
        gem[magsf_x1dir][i][j][k] =
            gem[magsf_x1dir][i][j][k] -
            B_sq * (gas[s][gasr_v1][i][j][k] - gas[s][gasl_v1][i][j][k]);
        gem[magsf_x2dir][i][j][k] =
            gem[magsf_x2dir][i][j][k] -
            B_sq * (gas[s][gasr_v2][i][j][k] - gas[s][gasl_v2][i][j][k]);
        gem[magsf_x3dir][i][j][k] =
            gem[magsf_x3dir][i][j][k] -
            B_sq * (gas[s][gasr_v3][i][j][k] - gas[s][gasl_v3][i][j][k]);
      }
    }
  }
}

void get_magnetic_stress_rusanov_j() {
  double bsq_l, bsq_r, B_tot_sq_l, B_tot_sq_r, B_sq, Bstress_1_p, Bstress_2_p,
      Bstress_3_p, Bstress_1_n, Bstress_2_n, Bstress_3_n;
  int i, j, k;
  int s = NS1 - 1;
  for (i = is; i <= ie; i++) {
    for (j = js; j <= je + 1; j++) {
      for (k = ks; k <= ke; k++) {
        bsq_l = (gem[magsl_b1][i][j][k] * gem[magsl_b1][i][j][k] +
                 gem[magsl_b2][i][j][k] * gem[magsl_b2][i][j][k] +
                 gem[magsl_b3][i][j][k] * gem[magsl_b3][i][j][k]) +
                2 * (gem[magsl_b1][i][j][k] * gem[B0AfaceJ_b1][i][j][k] +
                     gem[magsl_b2][i][j][k] * gem[B0AfaceJ_b2][i][j][k] +
                     gem[magsl_b3][i][j][k] * gem[B0AfaceJ_b3][i][j][k]);

        bsq_r = (gem[magsr_b1][i][j][k] * gem[magsr_b1][i][j][k] +
                 gem[magsr_b2][i][j][k] * gem[magsr_b2][i][j][k] +
                 gem[magsr_b3][i][j][k] * gem[magsr_b3][i][j][k]) +
                2 * (gem[magsr_b1][i][j][k] * gem[B0AfaceJ_b1][i][j][k] +
                     gem[magsr_b2][i][j][k] * gem[B0AfaceJ_b2][i][j][k] +
                     gem[magsr_b3][i][j][k] * gem[B0AfaceJ_b3][i][j][k]);

        Bstress_1_p =
            0.5 * (-(gem[magsl_b1][i][j][k] + gem[B0AfaceJ_b1][i][j][k]) *
                       gem[mag_bj][i][j][k] -
                   gem[magsl_b1][i][j][k] * gem[B0AfaceJ_bn][i][j][k]);
        Bstress_2_p =
            0.5 * (0.5 * bsq_l -
                   (gem[magsl_b2][i][j][k] + gem[B0AfaceJ_b2][i][j][k]) *
                       gem[mag_bj][i][j][k] -
                   gem[magsl_b2][i][j][k] * gem[B0AfaceJ_bn][i][j][k]);
        Bstress_3_p =
            0.5 * (-(gem[magsl_b3][i][j][k] + gem[B0AfaceJ_b3][i][j][k]) *
                       gem[mag_bj][i][j][k] -
                   gem[magsl_b3][i][j][k] * gem[B0AfaceJ_bn][i][j][k]);

        Bstress_1_n =
            0.5 * (-(gem[magsr_b1][i][j][k] + gem[B0AfaceJ_b1][i][j][k]) *
                       gem[mag_bj][i][j][k] -
                   gem[magsr_b1][i][j][k] * gem[B0AfaceJ_bn][i][j][k]);
        Bstress_2_n =
            0.5 * (0.5 * bsq_r -
                   (gem[magsr_b2][i][j][k] + gem[B0AfaceJ_b2][i][j][k]) *
                       gem[mag_bj][i][j][k] -
                   gem[magsr_b2][i][j][k] * gem[B0AfaceJ_bn][i][j][k]);
        Bstress_3_n =
            0.5 * (-(gem[magsr_b3][i][j][k] + gem[B0AfaceJ_b3][i][j][k]) *
                       gem[mag_bj][i][j][k] -
                   gem[magsr_b3][i][j][k] * gem[B0AfaceJ_bn][i][j][k]);

        gem[magsf_x1dir][i][j][k] = Bstress_1_p + Bstress_1_n;
        gem[magsf_x2dir][i][j][k] = Bstress_2_p + Bstress_2_n;
        gem[magsf_x3dir][i][j][k] = Bstress_3_p + Bstress_3_n;

        B_tot_sq_l = (gem[magsl_b1][i][j][k] + gem[B0AfaceJ_b1][i][j][k]) *
                         (gem[magsl_b1][i][j][k] + gem[B0AfaceJ_b1][i][j][k]) +
                     (gem[magsl_b2][i][j][k] + gem[B0AfaceJ_b2][i][j][k]) *
                         (gem[magsl_b2][i][j][k] + gem[B0AfaceJ_b2][i][j][k]) +
                     (gem[magsl_b3][i][j][k] + gem[B0AfaceJ_b3][i][j][k]) *
                         (gem[magsl_b3][i][j][k] + gem[B0AfaceJ_b3][i][j][k]);

        B_tot_sq_r = (gem[magsr_b1][i][j][k] + gem[B0AfaceJ_b1][i][j][k]) *
                         (gem[magsr_b1][i][j][k] + gem[B0AfaceJ_b1][i][j][k]) +
                     (gem[magsr_b2][i][j][k] + gem[B0AfaceJ_b2][i][j][k]) *
                         (gem[magsr_b2][i][j][k] + gem[B0AfaceJ_b2][i][j][k]) +
                     (gem[magsr_b3][i][j][k] + gem[B0AfaceJ_b3][i][j][k]) *
                         (gem[magsr_b3][i][j][k] + gem[B0AfaceJ_b3][i][j][k]);

        B_sq = (B_tot_sq_l + B_tot_sq_r) / 4 / CA;
        gem[magsf_x1dir][i][j][k] =
            gem[magsf_x1dir][i][j][k] -
            B_sq * (gas[s][gasr_v1][i][j][k] - gas[s][gasl_v1][i][j][k]);
        gem[magsf_x2dir][i][j][k] =
            gem[magsf_x2dir][i][j][k] -
            B_sq * (gas[s][gasr_v2][i][j][k] - gas[s][gasl_v2][i][j][k]);
        gem[magsf_x3dir][i][j][k] =
            gem[magsf_x3dir][i][j][k] -
            B_sq * (gas[s][gasr_v3][i][j][k] - gas[s][gasl_v3][i][j][k]);
      }
    }
  }
}

void get_magnetic_stress_rusanov_k() {
  double bsq_l, bsq_r, B_tot_sq_l, B_tot_sq_r, B_sq, Bstress_1_p, Bstress_2_p,
      Bstress_3_p, Bstress_1_n, Bstress_2_n, Bstress_3_n;
  int i, j, k;
  int s = NS1 - 1;
  for (i = is; i <= ie; i++) {
    for (j = js; j <= je; j++) {
      for (k = ks; k <= ke + 1; k++) {
        bsq_l = (gem[magsl_b1][i][j][k] * gem[magsl_b1][i][j][k] +
                 gem[magsl_b2][i][j][k] * gem[magsl_b2][i][j][k] +
                 gem[magsl_b3][i][j][k] * gem[magsl_b3][i][j][k]) +
                2 * (gem[magsl_b1][i][j][k] * gem[B0AfaceK_b1][i][j][k] +
                     gem[magsl_b2][i][j][k] * gem[B0AfaceK_b2][i][j][k] +
                     gem[magsl_b3][i][j][k] * gem[B0AfaceK_b3][i][j][k]);
        bsq_r = (gem[magsr_b1][i][j][k] * gem[magsr_b1][i][j][k] +
                 gem[magsr_b2][i][j][k] * gem[magsr_b2][i][j][k] +
                 gem[magsr_b3][i][j][k] * gem[magsr_b3][i][j][k]) +
                2 * (gem[magsr_b1][i][j][k] * gem[B0AfaceK_b1][i][j][k] +
                     gem[magsr_b2][i][j][k] * gem[B0AfaceK_b2][i][j][k] +
                     gem[magsr_b3][i][j][k] * gem[B0AfaceK_b3][i][j][k]);

        Bstress_1_p =
            0.5 * (-(gem[magsl_b1][i][j][k] + gem[B0AfaceK_b1][i][j][k]) *
                       gem[mag_bk][i][j][k] -
                   gem[magsl_b1][i][j][k] * gem[B0AfaceK_bn][i][j][k]);
        Bstress_2_p =
            0.5 * (-(gem[magsl_b2][i][j][k] + gem[B0AfaceK_b2][i][j][k]) *
                       gem[mag_bk][i][j][k] -
                   gem[magsl_b2][i][j][k] * gem[B0AfaceK_bn][i][j][k]);
        Bstress_3_p =
            0.5 * (0.5 * bsq_l -
                   (gem[magsl_b3][i][j][k] + gem[B0AfaceK_b3][i][j][k]) *
                       gem[mag_bk][i][j][k] -
                   gem[magsl_b3][i][j][k] * gem[B0AfaceK_bn][i][j][k]);

        Bstress_1_n =
            0.5 * (-(gem[magsr_b1][i][j][k] + gem[B0AfaceK_b1][i][j][k]) *
                       gem[mag_bk][i][j][k] -
                   gem[magsr_b1][i][j][k] * gem[B0AfaceK_bn][i][j][k]);
        Bstress_2_n =
            0.5 * (-(gem[magsr_b2][i][j][k] + gem[B0AfaceK_b2][i][j][k]) *
                       gem[mag_bk][i][j][k] -
                   gem[magsr_b2][i][j][k] * gem[B0AfaceK_bn][i][j][k]);
        Bstress_3_n =
            0.5 * (0.5 * bsq_r -
                   (gem[magsr_b3][i][j][k] + gem[B0AfaceK_b3][i][j][k]) *
                       gem[mag_bk][i][j][k] -
                   gem[magsr_b3][i][j][k] * gem[B0AfaceK_bn][i][j][k]);

        gem[magsf_x1dir][i][j][k] = Bstress_1_p + Bstress_1_n;
        gem[magsf_x2dir][i][j][k] = Bstress_2_p + Bstress_2_n;
        gem[magsf_x3dir][i][j][k] = Bstress_3_p + Bstress_3_n;

        B_tot_sq_l = (gem[magsl_b1][i][j][k] + gem[B0AfaceK_b1][i][j][k]) *
                         (gem[magsl_b1][i][j][k] + gem[B0AfaceK_b1][i][j][k]) +
                     (gem[magsl_b2][i][j][k] + gem[B0AfaceK_b2][i][j][k]) *
                         (gem[magsl_b2][i][j][k] + gem[B0AfaceK_b2][i][j][k]) +
                     (gem[magsl_b3][i][j][k] + gem[B0AfaceK_b3][i][j][k]) *
                         (gem[magsl_b3][i][j][k] + gem[B0AfaceK_b3][i][j][k]);

        B_tot_sq_r = (gem[magsr_b1][i][j][k] + gem[B0AfaceK_b1][i][j][k]) *
                         (gem[magsr_b1][i][j][k] + gem[B0AfaceK_b1][i][j][k]) +
                     (gem[magsr_b2][i][j][k] + gem[B0AfaceK_b2][i][j][k]) *
                         (gem[magsr_b2][i][j][k] + gem[B0AfaceK_b2][i][j][k]) +
                     (gem[magsr_b3][i][j][k] + gem[B0AfaceK_b3][i][j][k]) *
                         (gem[magsr_b3][i][j][k] + gem[B0AfaceK_b3][i][j][k]);

        B_sq = (B_tot_sq_l + B_tot_sq_r) / 4 / CA;
        gem[magsf_x1dir][i][j][k] =
            gem[magsf_x1dir][i][j][k] -
            B_sq * (gas[s][gasr_v1][i][j][k] - gas[s][gasl_v1][i][j][k]);
        gem[magsf_x2dir][i][j][k] =
            gem[magsf_x2dir][i][j][k] -
            B_sq * (gas[s][gasr_v2][i][j][k] - gas[s][gasl_v2][i][j][k]);
        gem[magsf_x3dir][i][j][k] =
            gem[magsf_x3dir][i][j][k] -
            B_sq * (gas[s][gasr_v3][i][j][k] - gas[s][gasl_v3][i][j][k]);
      }
    }
  }
}
void get_magnetic_stress_rusanov(dir_t dir) {
  switch (dir) {
    case i_dir:
      get_magnetic_stress_rusanov_i();
      break;
    case j_dir:
      get_magnetic_stress_rusanov_j();
      break;
    case k_dir:
      get_magnetic_stress_rusanov_k();
      break;
    default:
      break;
  }
}

void g2int(int dir, double wt2[12][12], double eta[12][12], double psi[12][12],
           double x0, double x1, double x2, double x3, double y0, double y1,
           double y2, double y3, double z0, double z1, double z2, double z3,
           double *f1, double *f2, double *f3, double *fn, double *fsq,
           double *fn1, double *fn2, double *fn3) {
  // Initialize sums
  double sum_fx = 0, sum_fy = 0, sum_fz = 0, sum_fnorm = 0, sum_fsq = 0,
         sum_fnx = 0, sum_fny = 0, sum_fnz = 0;
  double total_area = 0;
  double xn = 0, yn = 0, zn = 0;

  switch (dir) {
    case i_dir:
      xn = 1;
      yn = 0;
      zn = 0;
      break;

    case j_dir:
      xn = 0;
      yn = 1;
      zn = 0;
      break;

    case k_dir:
      xn = 0;
      yn = 0;
      zn = 1;
      break;
  }

  for (int i = 0; i < 12; i++) {
    for (int j = 0; j < 12; j++) {
      // this line equals to x0 + dx(1).*eta + dx(2).*psi + dx(3).*eta.*psi in
      // MATLAB here show a form for orthogonal grid, much easy
      double x = x0 + psi[i][j] * (x2 - x0) + eta[i][j] * (x1 - x0);
      double y = y0 + psi[i][j] * (y2 - y0) + eta[i][j] * (y1 - y0);
      double z = z0 + psi[i][j] * (z2 - z0) + eta[i][j] * (z1 - z0);

      double bx = B10(x, y, z);
      double by = B20(x, y, z);
      double bz = B30(x, y, z);

      // double xeta = (1-psi[i][j])*(x1-x0) + psi[i][j]*(x3-x2);
      // double yeta = (1-psi[i][j])*(y1-y0) + psi[i][j]*(y3-y2);
      // double zeta = (1-psi[i][j])*(z1-z0) + psi[i][j]*(z3-z2);

      // yeta = dy(1) + dy(3).*psi;
      // zeta = dz(1) + dz(3).*psi;

      // xpsi = dx(2) + dx(3).*eta;
      // ypsi = dy(2) + dy(3).*eta;
      // zpsi = dz(2) + dz(3).*eta;

      //      double area = wt2[i][j];  // Assume this is the differential area
      //      element sum_fx += bx * area; sum_fy += by * area; sum_fz += bz *
      //      area; double bdotn = (bx * xn + by * yn + bz * zn); sum_fnorm +=
      //          bdotn * area;  // Modify according to actual computation of
      //          fnorm
      //      sum_fsq += (bx * bx + by * by + bz * bz) * area;
      //      sum_fnx += bx * bdotn * area;  // Modify according to actual
      //      computation sum_fny += by * bdotn * area; sum_fnz += bz * bdotn *
      //      area; total_area += area;

      double area =
          0.25 * wt2[i][j];  // Assume this is the differential area element
      sum_fx += bx * area;
      sum_fy += by * area;
      sum_fz += bz * area;
      double bdotn = (bx * xn + by * yn + bz * zn);
      sum_fnorm +=
          bdotn * area;  // Modify according to actual computation of fnorm
      sum_fsq += (bx * bx + by * by + bz * bz) * area;
      sum_fnx += bx * bdotn * area;  // Modify according to actual computation
      sum_fny += by * bdotn * area;
      sum_fnz += bz * bdotn * area;
      total_area += area;
    }
  }

  // Output results
  *f1 = sum_fx / total_area;
  *f2 = sum_fy / total_area;
  *f3 = sum_fz / total_area;
  *fn = sum_fnorm / total_area;
  *fsq = sum_fsq / total_area;
  *fn1 = sum_fnx / total_area;
  *fn2 = sum_fny / total_area;
  *fn3 = sum_fnz / total_area;
}

// Here g2int_ortho and g3int_ortho are only for background magnetic field calculation
void g2int_ortho(double (*B10)(double, double, double),double (*B20)(double, double, double),double (*B30)(double, double, double),\
            double (*H1)(double, double, double),double (*H2)(double, double, double),double (*H3)(double, double, double),\
            dir_t dir, double wt2[12][12], double eta[12][12], double psi[12][12], double face_area, \
            double x1_0, double x1_1, double x1_2, double x1_3, double x2_0, double x2_1, double x2_2, double x2_3, double x3_0, double x3_1, double x3_2, double x3_3, \
            double *f1, double *f2, double *f3, double *fn, double *fsq, double *fn1, double *fn2, double *fn3) {


  // Initialize sums
  double sum_f1 = 0, sum_f2 = 0, sum_f3 = 0, sum_fnorm = 0, sum_fsq = 0, sum_fn1 = 0, sum_fn2 = 0, sum_fn3 = 0;
  double x1,x2,x3,dx1,dx2,dx3,b1,b2,b3,h1,h2,h3;

  if (fabs(face_area) >= 1e-10){

    switch (dir)
    {
      case i_dir:
        for (int i = 0; i < 12; i++) {
          for (int j = 0; j < 12; j++) {
            // this line equals to x0 + dx(1).*eta + dx(2).*psi + dx(3).*eta.*psi in MATLAB
            // here show a form for orthogonal grid, much easy
            x1 = x1_0 + psi[i][j]*(x1_2-x1_0) + eta[i][j]*(x1_1-x1_0);
            x2 = x2_0 + psi[i][j]*(x2_2-x2_0) + eta[i][j]*(x2_1-x2_0);
            x3 = x3_0 + psi[i][j]*(x3_2-x3_0) + eta[i][j]*(x3_1-x3_0);
            dx1 = fmax(fmax(fabs(x1_1-x1_0), fabs(x1_2-x1_0)), fabs(x1_3-x1_0));
            dx2 = fmax(fmax(fabs(x2_1-x2_0), fabs(x2_2-x2_0)), fabs(x2_3-x2_0));
            dx3 = fmax(fmax(fabs(x3_1-x3_0), fabs(x3_2-x3_0)), fabs(x3_3-x3_0));

            b1 = B10(x1, x2, x3);
            b2 = B20(x1, x2, x3);
            b3 = B30(x1, x2, x3);
            h1 = H1(x1, x2, x3);
            h2 = H2(x1, x2, x3);
            h3 = H3(x1, x2, x3);

            double bn = b1;
            double bsq = b1 * b1 + b2 * b2 + b3 * b3;
            double bn1 = b1 * b1;
            double bn2 = b1 * b2;
            double bn3 = b1 * b3;

            sum_f1 += 0.25*dx2*dx3*wt2[i][j]*b1*h2*h3;
            sum_f2 += 0.25*dx2*dx3*wt2[i][j]*b2*h2*h3;
            sum_f3 += 0.25*dx2*dx3*wt2[i][j]*b3*h2*h3;
            sum_fnorm += 0.25*dx2*dx3*wt2[i][j]*b1*h2*h3;
            sum_fsq += 0.25*dx2*dx3*wt2[i][j]*bsq*h2*h3;
            sum_fn1 += 0.25*dx2*dx3*wt2[i][j]*bn1*h2*h3;
            sum_fn2 += 0.25*dx2*dx3*wt2[i][j]*bn2*h2*h3;
            sum_fn3 += 0.25*dx2*dx3*wt2[i][j]*bn3*h2*h3;
          }
        }
        // Output results
        *f1 = sum_f1 / face_area;
        *f2 = sum_f2 / face_area;
        *f3 = sum_f3 / face_area;
        *fn = sum_fnorm / face_area;
        *fsq = sum_fsq / face_area;
        *fn1 = sum_fn1 / face_area;
        *fn2 = sum_fn2 / face_area;
        *fn3 = sum_fn3 / face_area;
        break;

      case j_dir:
        for (int i = 0; i < 12; i++) {
          for (int j = 0; j < 12; j++) {
            // this line equals to x0 + dx(1).*eta + dx(2).*psi + dx(3).*eta.*psi in MATLAB
            // here show a form for orthogonal grid, much easy
            x1 = x1_0 + psi[i][j]*(x1_2-x1_0) + eta[i][j]*(x1_1-x1_0);
            x2 = x2_0 + psi[i][j]*(x2_2-x2_0) + eta[i][j]*(x2_1-x2_0);
            x3 = x3_0 + psi[i][j]*(x3_2-x3_0) + eta[i][j]*(x3_1-x3_0);
            dx1 = fmax(fmax(fabs(x1_1-x1_0), fabs(x1_2-x1_0)), fabs(x1_3-x1_0));
            dx2 = fmax(fmax(fabs(x2_1-x2_0), fabs(x2_2-x2_0)), fabs(x2_3-x2_0));
            dx3 = fmax(fmax(fabs(x3_1-x3_0), fabs(x3_2-x3_0)), fabs(x3_3-x3_0));

            b1 = B10(x1, x2, x3);
            b2 = B20(x1, x2, x3);
            b3 = B30(x1, x2, x3);
            h1 = H1(x1, x2, x3);
            h2 = H2(x1, x2, x3);
            h3 = H3(x1, x2, x3);

            double bn = b2;
            double bsq = b1 * b1 + b2 * b2 + b3 * b3;
            double bn1 = b2 * b1;
            double bn2 = b2 * b2;
            double bn3 = b2 * b3;

            sum_f1 += 0.25*dx1*dx3*wt2[i][j]*b1*h1*h3;
            sum_f2 += 0.25*dx1*dx3*wt2[i][j]*b2*h1*h3;
            sum_f3 += 0.25*dx1*dx3*wt2[i][j]*b3*h1*h3;
            sum_fnorm += 0.25*dx1*dx3*wt2[i][j]*bn*h1*h3;
            sum_fsq += 0.25*dx1*dx3*wt2[i][j]*bsq*h1*h3;
            sum_fn1 += 0.25*dx1*dx3*wt2[i][j]*bn1*h1*h3;
            sum_fn2 += 0.25*dx1*dx3*wt2[i][j]*bn2*h1*h3;
            sum_fn3 += 0.25*dx1*dx3*wt2[i][j]*bn3*h1*h3;
          }
        }
        // Output results
        *f1 = sum_f1 / face_area;
        *f2 = sum_f2 / face_area;
        *f3 = sum_f3 / face_area;
        *fn = sum_fnorm / face_area;
        *fsq = sum_fsq / face_area;
        *fn1 = sum_fn1 / face_area;
        *fn2 = sum_fn2 / face_area;
        *fn3 = sum_fn3 / face_area;
        break;
      case k_dir:
        for (int i = 0; i < 12; i++) {
          for (int j = 0; j < 12; j++) {
            // this line equals to x0 + dx(1).*eta + dx(2).*psi + dx(3).*eta.*psi in MATLAB
            // here show a form for orthogonal grid, much easy
            x1 = x1_0 + psi[i][j]*(x1_2-x1_0) + eta[i][j]*(x1_1-x1_0);
            x2 = x2_0 + psi[i][j]*(x2_2-x2_0) + eta[i][j]*(x2_1-x2_0);
            x3 = x3_0 + psi[i][j]*(x3_2-x3_0) + eta[i][j]*(x3_1-x3_0);
            dx1 = fmax(fmax(fabs(x1_1-x1_0), fabs(x1_2-x1_0)), fabs(x1_3-x1_0));
            dx2 = fmax(fmax(fabs(x2_1-x2_0), fabs(x2_2-x2_0)), fabs(x2_3-x2_0));
            dx3 = fmax(fmax(fabs(x3_1-x3_0), fabs(x3_2-x3_0)), fabs(x3_3-x3_0));

            b1 = B10(x1, x2, x3);
            b2 = B20(x1, x2, x3);
            b3 = B30(x1, x2, x3);
            h1 = H1(x1, x2, x3);
            h2 = H2(x1, x2, x3);
            h3 = H3(x1, x2, x3);

            double bn = b3;
            double bsq = b1 * b1 + b2 * b2 + b3 * b3;
            double bn1 = b3 * b1;
            double bn2 = b3 * b2;
            double bn3 = b3 * b3;

            sum_f1 += 0.25*dx1*dx2*wt2[i][j]*b1*h1*h2;
            sum_f2 += 0.25*dx1*dx2*wt2[i][j]*b2*h1*h2;
            sum_f3 += 0.25*dx1*dx2*wt2[i][j]*b3*h1*h2;
            sum_fnorm += 0.25*dx1*dx2*wt2[i][j]*bn*h1*h2;
            sum_fsq += 0.25*dx1*dx2*wt2[i][j]*bsq*h1*h2;
            sum_fn1 += 0.25*dx1*dx2*wt2[i][j]*bn1*h1*h2;
            sum_fn2 += 0.25*dx1*dx2*wt2[i][j]*bn2*h1*h2;
            sum_fn3 += 0.25*dx1*dx2*wt2[i][j]*bn3*h1*h2;
          }
        }
        // Output results
        *f1 = sum_f1 / face_area;
        *f2 = sum_f2 / face_area;
        *f3 = sum_f3 / face_area;
        *fn = sum_fnorm / face_area;
        *fsq = sum_fsq / face_area;
        *fn1 = sum_fn1 / face_area;
        *fn2 = sum_fn2 / face_area;
        *fn3 = sum_fn3 / face_area;
        break;
    }
  } else{
    switch (dir){
      case i_dir:
        *f1 = GaussianLineIntegral(B10, x1_0,x2_0,x3_0,x1_3,x2_3,x3_3);
        *f2 = GaussianLineIntegral(B20, x1_0,x2_0,x3_0,x1_3,x2_3,x3_3);
        *f3 = GaussianLineIntegral(B30, x1_0,x2_0,x3_0,x1_3,x2_3,x3_3);
        *fn = *f1;
        *fsq = (*f1)*(*f1) + (*f2)*(*f2) + (*f3)*(*f3);
        *fn1 = (*f1)*(*f1);
        *fn2 = (*f1)*(*f2);
        *fn3 = (*f1)*(*f3);
        break;
      case j_dir:
        *f1 = GaussianLineIntegral(B10, x1_0,x2_0,x3_0,x1_3,x2_3,x3_3);
        *f2 = GaussianLineIntegral(B20, x1_0,x2_0,x3_0,x1_3,x2_3,x3_3);
        *f3 = GaussianLineIntegral(B30, x1_0,x2_0,x3_0,x1_3,x2_3,x3_3);
        *fn = *f2;
        *fsq = (*f1)*(*f1) + (*f2)*(*f2) + (*f3)*(*f3);
        *fn1 = (*f2)*(*f1);
        *fn2 = (*f2)*(*f2);
        *fn3 = (*f2)*(*f3);
        break;
      case k_dir:
        *f1 = GaussianLineIntegral(B10, x1_0,x2_0,x3_0,x1_3,x2_3,x3_3);
        *f2 = GaussianLineIntegral(B20, x1_0,x2_0,x3_0,x1_3,x2_3,x3_3);
        *f3 = GaussianLineIntegral(B30, x1_0,x2_0,x3_0,x1_3,x2_3,x3_3);
        *fn = *f3;
        *fsq = (*f1)*(*f1) + (*f2)*(*f2) + (*f3)*(*f3);
        *fn1 = (*f3)*(*f1);
        *fn2 = (*f3)*(*f2);
        *fn3 = (*f3)*(*f3);
        break;
    }
  }
}

void g3int_ortho(double (*Source1_B0)(double, double, double),double (*Source2_B0)(double, double, double),double (*Source3_B0)(double, double, double),\
            double (*H1)(double, double, double),double (*H2)(double, double, double),double (*H3)(double, double, double),\
            double wt3[12][12][12], double eta[12][12][12], double psi[12][12][12], double zeta[12][12][12], double volume, \
            double x1n, double x1p, double x2n, double x2p, double x3n, double x3p,\
            double *S1_B0, double *S2_B0, double *S3_B0) {


  // Initialize sums
  double sum_S1 = 0, sum_S2 = 0, sum_S3 = 0;
  double x1,x2,x3,dx1,dx2,dx3,S1,S2,S3,h1,h2,h3;
  dx1 = x1p-x1n;
  dx2 = x2p-x2n;
  dx3 = x3p-x3n;

  for (int i = 0; i < 12; i++) {
    for (int j = 0; j < 12; j++) {
      for (int k = 0; k < 12; k++) {
        x1 = x1n + eta[i][j][k]*dx1;
        x2 = x2n + psi[i][j][k]*dx2;
        x3 = x3n + zeta[i][j][k]*dx3;

        S1 = Source1_B0(x1, x2, x3);
        S2 = Source2_B0(x1, x2, x3);
        S3 = Source3_B0(x1, x2, x3);
        h1 = H1(x1, x2, x3);
        h2 = H2(x1, x2, x3);
        h3 = H3(x1, x2, x3);
        sum_S1 += 0.125*dx1*dx2*dx3*wt3[i][j][k]*S1*h1*h2*h3;
        sum_S2 += 0.125*dx1*dx2*dx3*wt3[i][j][k]*S2*h1*h2*h3;
        sum_S3 += 0.125*dx1*dx2*dx3*wt3[i][j][k]*S3*h1*h2*h3;
      }
    }
  }

  // Output results
  *S1_B0 = sum_S1 / volume;
  *S2_B0 = sum_S2 / volume;
  *S3_B0 = sum_S3 / volume;
}

void set_background_field() {
  // Gaussian positions and weights for 2-D Gaussian Quadrature (12x12 point)

  log_info("Initializing Background Fields...\n");

  // initialize the B0_ptr
  #pragma omp parallel for collapse(3) schedule(static)
  for (int i = 0; i <= ieg; i++) {
    for (int j = 0; j <= jeg; j++) {
      for (int k = 0; k <= keg; k++) {
        gem[B0_b1][i][j][k] = B10(x1ctr[i][j][k], x2ctr[i][j][k], x3ctr[i][j][k]);
        gem[B0_b2][i][j][k] = B20(x1ctr[i][j][k], x2ctr[i][j][k], x3ctr[i][j][k]);
        gem[B0_b3][i][j][k] = B30(x1ctr[i][j][k], x2ctr[i][j][k], x3ctr[i][j][k]);
      }
    }
  }

  double a[] = {0.1252334085,  0.3678314989,  0.5873179542,  0.7699026741,
                0.9041172563,  0.9815606342,  -0.1252334085, -0.3678314989,
                -0.5873179542, -0.7699026741, -0.9041172563, -0.9815606342};
  double wt[] = {0.2491470458, 0.2334925365, 0.2031674267, 0.1600783285,
                 0.1069393259, 0.0471753363, 0.2491470458, 0.2334925365,
                 0.2031674267, 0.1600783285, 0.1069393259, 0.0471753363};

  double eta[12][12], psi[12][12], wt2[12][12];
  for (int i = 0; i < 12; i++) {
    for (int j = 0; j < 12; j++) {
      eta[i][j] = 0.5 * (1.0 + a[i]);
      psi[i][j] = 0.5 * (1.0 + a[j]);
      wt2[i][j] = wt[i] * wt[j];
    }
  }

  // Gaussian quadrature on i-faces
  #pragma omp parallel for collapse(3) schedule(static)
  for (int i = 0; i <= ieg + 1; i++) {
    for (int j = 0; j <= jeg; j++) {
      for (int k = 0; k <= keg; k++) {
        g2int_ortho(B10, B20, B30, &H1, &H2, &H3,
              i_dir, wt2, eta, psi, geo[face_idir][i][j][k],
              x1[i][j][k], x1[i][j + 1][k], x1[i][j][k + 1], x1[i][j + 1][k + 1],
              x2[i][j][k], x2[i][j + 1][k], x2[i][j][k + 1], x2[i][j + 1][k + 1],
              x3[i][j][k], x3[i][j + 1][k], x3[i][j][k + 1], x3[i][j + 1][k + 1],
              &gem[B0AfaceI_b1][i][j][k],
              &gem[B0AfaceI_b2][i][j][k], &gem[B0AfaceI_b3][i][j][k],
              &gem[B0AfaceI_bn][i][j][k], &gem[B0AfaceI_bsq][i][j][k],
              &gem[B0AfaceI_bn1][i][j][k], &gem[B0AfaceI_bn2][i][j][k],
              &gem[B0AfaceI_bn3][i][j][k]);
      }
    }
  }
  // Gaussian quadrature on j-faces
  #pragma omp parallel for collapse(3) schedule(static)
  for (int i = 0; i <= ieg; i++) {
    for (int j = 0; j <= jeg + 1; j++) {
      for (int k = 0; k <= keg; k++) {
        g2int_ortho(B10, B20, B30, &H1, &H2, &H3,
              j_dir, wt2, eta, psi, geo[face_jdir][i][j][k],
              x1[i][j][k], x1[i][j][k + 1], x1[i + 1][j][k], x1[i + 1][j][k + 1],
              x2[i][j][k], x2[i][j][k + 1], x2[i + 1][j][k], x2[i + 1][j][k + 1],
              x3[i][j][k], x3[i][j][k + 1], x3[i + 1][j][k], x3[i + 1][j][k + 1],
              &gem[B0AfaceJ_b1][i][j][k],
              &gem[B0AfaceJ_b2][i][j][k], &gem[B0AfaceJ_b3][i][j][k],
              &gem[B0AfaceJ_bn][i][j][k], &gem[B0AfaceJ_bsq][i][j][k],
              &gem[B0AfaceJ_bn1][i][j][k], &gem[B0AfaceJ_bn2][i][j][k],
              &gem[B0AfaceJ_bn3][i][j][k]);
      }
    }
  }
  // Gaussian quadrature on k-faces
  #pragma omp parallel for collapse(3) schedule(static)
  for (int i = 0; i <= ieg; i++) {
    for (int j = 0; j <= jeg; j++) {
      for (int k = 0; k <= keg + 1; k++) {
        g2int_ortho(B10, B20, B30, &H1, &H2, &H3,
              k_dir, wt2, eta, psi, geo[face_kdir][i][j][k],
              x1[i][j][k], x1[i + 1][j][k], x1[i][j + 1][k], x1[i + 1][j + 1][k],
              x2[i][j][k], x2[i + 1][j][k], x2[i][j + 1][k], x2[i + 1][j + 1][k],
              x3[i][j][k], x3[i + 1][j][k], x3[i][j + 1][k], x3[i + 1][j + 1][k],
              &gem[B0AfaceK_b1][i][j][k],
              &gem[B0AfaceK_b2][i][j][k], &gem[B0AfaceK_b3][i][j][k],
              &gem[B0AfaceK_bn][i][j][k], &gem[B0AfaceK_bsq][i][j][k],
              &gem[B0AfaceK_bn1][i][j][k], &gem[B0AfaceK_bn2][i][j][k],
              &gem[B0AfaceK_bn3][i][j][k]);
      }
    }
  }

  // integrate the background b field on edges for electric field calculations
  // i-edge background B field, integrate B along i edge
  #pragma omp parallel for collapse(3) schedule(static)
  for (int i = is; i <= ie; i++) {
    for (int j = js; j <= je + 1; j++) {
      for (int k = ks; k <= ke + 1; k++) {
        gem[B0AedgeI_b1][i][j][k] = GaussianLineIntegral(
            *B10, x1[i + 1][j][k], x2[i + 1][j][k], x3[i + 1][j][k],
            x1[i][j][k], x2[i][j][k], x3[i][j][k]);
        gem[B0AedgeI_b2][i][j][k] = GaussianLineIntegral(
            *B20, x1[i + 1][j][k], x2[i + 1][j][k], x3[i + 1][j][k],
            x1[i][j][k], x2[i][j][k], x3[i][j][k]);
        gem[B0AedgeI_b3][i][j][k] = GaussianLineIntegral(
            *B30, x1[i + 1][j][k], x2[i + 1][j][k], x3[i + 1][j][k],
            x1[i][j][k], x2[i][j][k], x3[i][j][k]);
      }
    }
  }
  // integrate B along j edge
  #pragma omp parallel for collapse(3) schedule(static)
  for (int i = is; i <= ie + 1; i++) {
    for (int j = js; j <= je; j++) {
      for (int k = ks; k <= ke + 1; k++) {
        gem[B0AedgeJ_b1][i][j][k] = GaussianLineIntegral(
            *B10, x1[i][j + 1][k], x2[i][j + 1][k], x3[i][j + 1][k],
            x1[i][j][k], x2[i][j][k], x3[i][j][k]);
        gem[B0AedgeJ_b2][i][j][k] = GaussianLineIntegral(
            *B20, x1[i][j + 1][k], x2[i][j + 1][k], x3[i][j + 1][k],
            x1[i][j][k], x2[i][j][k], x3[i][j][k]);
        gem[B0AedgeJ_b3][i][j][k] = GaussianLineIntegral(
            *B30, x1[i][j + 1][k], x2[i][j + 1][k], x3[i][j + 1][k],
            x1[i][j][k], x2[i][j][k], x3[i][j][k]);
      }
    }
  }
  // integrate B along k edge
  #pragma omp parallel for collapse(3) schedule(static)
  for (int i = is; i <= ie + 1; i++) {
    for (int j = js; j <= je + 1; j++) {
      for (int k = ks; k <= ke; k++) {
        gem[B0AedgeK_b1][i][j][k] = GaussianLineIntegral(
            *B10, x1[i][j][k + 1], x2[i][j][k + 1], x3[i][j][k + 1],
            x1[i][j][k], x2[i][j][k], x3[i][j][k]);
        gem[B0AedgeK_b2][i][j][k] = GaussianLineIntegral(
            *B20, x1[i][j][k + 1], x2[i][j][k + 1], x3[i][j][k + 1],
            x1[i][j][k], x2[i][j][k], x3[i][j][k]);
        gem[B0AedgeK_b3][i][j][k] = GaussianLineIntegral(
            *B30, x1[i][j][k + 1], x2[i][j][k + 1], x3[i][j][k + 1],
            x1[i][j][k], x2[i][j][k], x3[i][j][k]);
      }
    }
  }

  double eta_3d[12][12][12], psi_3d[12][12][12], zeta_3d[12][12][12], wt3[12][12][12];
  #pragma omp parallel for collapse(3) schedule(static)
  for (int i = 0; i < 12; i++) {
    for (int j = 0; j < 12; j++) {
      for (int k = 0; k < 12; k++) {
        eta_3d[i][j][k] = 0.5 * (1.0 + a[i]);
        psi_3d[i][j][k] = 0.5 * (1.0 + a[j]);
        zeta_3d[i][j][k] = 0.5 * (1.0 + a[k]);
        wt3[i][j][k] = wt[i] * wt[j] * wt[k];
      }
    }
  }

  // calculates the Lorentz force from the background magnetic field
  #pragma omp parallel for collapse(3) schedule(static)
  for (int i = 0; i <= ieg; i++) {
    for (int j = 0; j <= jeg; j++) {
      for (int k = 0; k <= keg; k++) {
        double I_stress_X_1, I_stress_X_0, I_stress_Y_1, I_stress_Y_0, I_stress_Z_1,
            I_stress_Z_0, J_stress_X_1, J_stress_X_0, J_stress_Y_1, J_stress_Y_0,
            J_stress_Z_1, J_stress_Z_0, K_stress_X_1, K_stress_X_0, K_stress_Y_1,
            K_stress_Y_0, K_stress_Z_1, K_stress_Z_0;
        double S1_B0, S2_B0, S3_B0;
        g3int_ortho(Source1_B0, Source2_B0, Source3_B0, &H1, &H2, &H3, wt3, eta_3d, psi_3d, zeta_3d,
                    geo[vol_center][i][j][k],x1[i][j][k], x1[i+1][j][k], x2[i][j][k], x2[i][j+1][k], x3[i][j][k], x3[i][j][k+1],
                    &S1_B0, &S2_B0, &S3_B0);
        // x1-direction
        I_stress_X_1 = 0.5 * gem[B0AfaceI_bsq][i + 1][j][k] -
                       gem[B0AfaceI_bn1][i + 1][j][k];
        I_stress_X_0 =
            0.5 * gem[B0AfaceI_bsq][i][j][k] - gem[B0AfaceI_bn1][i][j][k];
        J_stress_X_1 = -gem[B0AfaceJ_bn1][i][j + 1][k];
        J_stress_X_0 = -gem[B0AfaceJ_bn1][i][j][k];
        K_stress_X_1 = -gem[B0AfaceK_bn1][i][j][k + 1];
        K_stress_X_0 = -gem[B0AfaceK_bn1][i][j][k];
        gem[dmagsf_G_x1dir][i][j][k] =
            (I_stress_X_1 * geo[face_idir][i + 1][j][k] -
             I_stress_X_0 * geo[face_idir][i][j][k] +
             J_stress_X_1 * geo[face_jdir][i][j + 1][k] -
             J_stress_X_0 * geo[face_jdir][i][j][k] +
             K_stress_X_1 * geo[face_kdir][i][j][k + 1] -
             K_stress_X_0 * geo[face_kdir][i][j][k]) /
            geo[vol_center][i][j][k]
            - S1_B0;

        // x2-direction
        I_stress_Y_1 = -gem[B0AfaceI_bn2][i + 1][j][k];
        I_stress_Y_0 = -gem[B0AfaceI_bn2][i][j][k];
        J_stress_Y_1 = 0.5 * gem[B0AfaceJ_bsq][i][j + 1][k] -
                       gem[B0AfaceJ_bn2][i][j + 1][k];
        J_stress_Y_0 =
            0.5 * gem[B0AfaceJ_bsq][i][j][k] - gem[B0AfaceJ_bn2][i][j][k];
        K_stress_Y_1 = -gem[B0AfaceK_bn2][i][j][k + 1];
        K_stress_Y_0 = -gem[B0AfaceK_bn2][i][j][k];
        gem[dmagsf_G_x2dir][i][j][k] =
            (I_stress_Y_1 * geo[face_idir][i + 1][j][k] -
             I_stress_Y_0 * geo[face_idir][i][j][k] +
             J_stress_Y_1 * geo[face_jdir][i][j + 1][k] -
             J_stress_Y_0 * geo[face_jdir][i][j][k] +
             K_stress_Y_1 * geo[face_kdir][i][j][k + 1] -
             K_stress_Y_0 * geo[face_kdir][i][j][k]) /
            geo[vol_center][i][j][k]
            - S2_B0;

        // x3-direction
        I_stress_Z_1 = -gem[B0AfaceI_bn3][i + 1][j][k];
        I_stress_Z_0 = -gem[B0AfaceI_bn3][i][j][k];
        J_stress_Z_1 = -gem[B0AfaceJ_bn3][i][j + 1][k];
        J_stress_Z_0 = -gem[B0AfaceJ_bn3][i][j][k];
        K_stress_Z_1 = 0.5 * gem[B0AfaceK_bsq][i][j][k + 1] -
                       gem[B0AfaceK_bn3][i][j][k + 1];
        K_stress_Z_0 =
            0.5 * gem[B0AfaceK_bsq][i][j][k] - gem[B0AfaceK_bn3][i][j][k];
        gem[dmagsf_G_x3dir][i][j][k] =
            (I_stress_Z_1 * geo[face_idir][i + 1][j][k] -
             I_stress_Z_0 * geo[face_idir][i][j][k] +
             J_stress_Z_1 * geo[face_jdir][i][j + 1][k] -
             J_stress_Z_0 * geo[face_jdir][i][j][k] +
             K_stress_Z_1 * geo[face_kdir][i][j][k + 1] -
             K_stress_Z_0 * geo[face_kdir][i][j][k]) /
            geo[vol_center][i][j][k]
            - S3_B0;
      }
    }
  }
}

void reset_rho() {
  for (int s = 0; s < NS1; s++) {
    for (int i = is; i <= ie; i++) {
      for (int j = js; j <= je; j++) {
        for (int k = ks; k <= ke; k++) {
          if (gas[s][gas_rho][i][j][k] < rho_floor) {
            log_error(
                "Resetting density to floor value at position: %d %d %d, "
                "original value: %f\n",
                x1c[i][j][k], x2c[i][j][k], x3c[i][j][k],
                gas[s][gas_rho][i][j][k]);
            gas[s][gas_rho][i][j][k] = rho_floor;
          }
        }
      }
    }
  }
}

void reset_p() {
  for (int s = 0; s < NS1; s++) {
    for (int i = is; i <= ie; i++) {
      for (int j = js; j <= je; j++) {
        for (int k = ks; k <= ke; k++) {
          if (gas[s][gas_p][i][j][k] < p_floor) {
            log_error(
                "Resetting pressure to floor value at position: %d %d %d, "
                "original value: %f\n",
                x1c[i][j][k], x2c[i][j][k], x3c[i][j][k],
                gas[s][gas_p][i][j][k]);
            gas[s][gas_p][i][j][k] = p_floor;
          }
        }
      }
    }
  }
}

void check_positivity() {
  #pragma omp parallel for collapse(4) schedule(static)
  for (int s = 0; s < NS1; s++) {
    for (int i = is; i <= ie; i++) {
      for (int j = js; j <= je; j++) {
        for (int k = ks; k <= ke; k++) {
          if (gas[s][gas_rho][i][j][k] < rho_floor) {
            gas[s][gas_rho][i][j][k] = rho_floor;
            gas[s][gas_p][i][j][k] = p_floor;
          }
          if (gas[s][gas_p][i][j][k] < p_floor) {
//             gas[s][gas_p][i][j][k] = fmax(p_floor, gas[s][gas_p_S][i][j][k]);
             gas[s][gas_p][i][j][k] = p_floor;
          }
        }
      }
    }
  }
}


// minmod limiter function
double minmod(double x, double y) {
  return 0.5*(fsign(x) + fsign(y))*fmin(fabs(x), fabs(y));
}

double DipoleL(double r, double theta, double phi) {
  double L = r / sin(theta) / sin(theta);
  return L;
}

void PPM(double fm2, double fm1, double f, double fp1, double fp2, double* vm, double* vp) {
  double dvm2 = fm1 - fm2;
  double dvm1 = f - fm1;
  double dvp1 = fp1 - f;
  double dvp2 = fp2 - fp1;

  double dvc = 0.5 * (dvm2 + dvm1);
  double SM = 2.0 * minmod(dvm2, dvm1);
  double Sm1 = minmod(dvc, SM);

  dvc = 0.5 * (dvm1 + dvp1);
  SM = 2.0 * minmod(dvm1, dvp1);
  double Sp1 = minmod(dvc, SM);

  dvc = 0.5 * (dvp2 + dvp1);
  SM = 2.0 * minmod(dvp2, dvp1);
  double Sp2 = minmod(dvc, SM);

  *vp = 0.5 * (f + fp1) - (Sp2 - Sp1) / 6.0;
  *vm = 0.5 * (f + fm1) - (Sp1 - Sm1) / 6.0;

  double ap = *vp - f;
  double am = *vm - f;

  if (ap * am >= 0.0) {
    ap = 0.0;
    am = 0.0;
  } else {
    if (fabs(ap) >= 2.0 * fabs(am)) {
      ap = -2.0 * am;
    }
    if (fabs(am) >= 2.0 * fabs(ap)) {
      am = -2.0 * ap;
    }
  }

  *vp = f + ap;
  *vm = f + am;
}

void FilterBaseMode(double* Q, int N, double* Q0)
{
  const double scl = 2.0 / N;
  const double dtheta = 2.0 * PI / N;

  double a0 = 0.0, a1 = 0.0, b1 = 0.0;

  for (int i = 0; i < N; i++) {
    double theta = dtheta * (i + 0.5);
    a0 += Q[i];
    a1 += Q[i] * cos(theta);
    b1 += Q[i] * sin(theta);
  }

  a1 *= scl;
  b1 *= scl;
  a0 *= 0.5 * scl;

  for (int i = 0; i < N; i++) {
    double theta = dtheta * (i + 0.5);
    Q0[i] = a0 + a1*cos(theta) + b1*sin(theta);
    Q[i] -= Q0[i];
  }
}

void RingBlockAvg(double* Q, int N, int Nchunk)
{
  int ncell = N / Nchunk;
  double* blk_avg = (double*)malloc(Nchunk * sizeof(double));

  for (int ib = 0; ib < Nchunk; ib++) {
    int ist = ib * ncell;
    double sum = 0.0;

    for (int ic = 0; ic < ncell; ic++) {
      sum += Q[ist + ic];
    }
    blk_avg[ib] = sum / ncell;
  }

  double* ext_avg = (double*)malloc((Nchunk + 4) * sizeof(double));
  ext_avg[0] = blk_avg[Nchunk - 2];
  ext_avg[1] = blk_avg[Nchunk - 1];

  for (int ib = 0; ib < Nchunk; ib++) {
    ext_avg[ib + 2] = blk_avg[ib];
  }

  ext_avg[Nchunk + 2 + 0] = blk_avg[0];
  ext_avg[Nchunk + 2 + 1] = blk_avg[1];

  for (int ib = 2; ib < Nchunk + 2; ib++) {
    int ig = (ib - 2) * ncell;

    double fL, fR;
    PPM(ext_avg[ib-2], ext_avg[ib-1], ext_avg[ib], ext_avg[ib+1], ext_avg[ib+2], &fL, &fR);

    double a = 3.0 * (fL + fR - 2.0 * ext_avg[ib]);
    double b = 2.0 * (3.0 * ext_avg[ib] - fR - 2.0 * fL);
    double c = fL;

    for (int ic = 0; ic < ncell; ic++) {
      double kl = (double)ic + 1.0;
      Q[ig + ic] = (a/3.0/ncell/ncell)*(3*kl*kl - 3*kl + 1) + 0.5*b*(2*kl - 1)/ncell + c;
    }
  }

  free(blk_avg);
  free(ext_avg);
}


