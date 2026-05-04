// This file contains general functions for spherical geometry, including:
// Lame coefficients
// Geometry calculation
// Calculation of  reconstruction weights
// Source term
// get dt when ring average
// Transformation between SPHERICAL and CARTESIAN

#include <dirent.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <omp.h>

#include "log.h"

#include "curvilinear.h"
#include "config.h"
#include "solver.h"
#include "setup_mpi.h"
#include "common.h"
#include "spsolver.h"
#include "lu_solver.h"
#include "spsolver.h"
#include "lu_solver.h"

double H1(double x1, double x2, double x3) {
  return 1.0;
}

double H2(double x1, double x2, double x3) {
  return x1;
}

double H3(double x1, double x2, double x3) {
  return x1 * sin(x2);
}

int set_geometry_arrays() {
  // x1c
  for (int i = isg; i <= ieg; i++) {
    for (int j = jsg; j <= jeg + 1; j++) {
      for (int k = ksg; k <= keg + 1; k++) {
        x1c[i][j][k] =
            0.5*(x1[i][j][k] + x1[i + 1][j][k]);
      }
    }
  }
  // x2c
  for (int i = isg; i <= ieg + 1; i++) {
    for (int j = jsg; j <= jeg; j++) {
      for (int k = ksg; k <= keg + 1; k++) {
        x2c[i][j][k] =
            0.5*(x2[i][j][k] + x2[i][j + 1][k]);
      }
    }
  }
  // x3c
  for (int i = isg; i <= ieg + 1; i++) {
    for (int j = jsg; j <= jeg + 1; j++) {
      for (int k = ksg; k <= keg; k++) {
        x3c[i][j][k] =
            0.5*(x3[i][j][k] + x3[i][j][k + 1]);
      }
    }
  }

  // dx1 dx2 dx3 and Edges ------
  // ni x (nj+1) x (nk+1)
  // dx1 and edge_idir
  for (int i = isg; i <= (ieg); i++) {
    for (int j = jsg; j <= (jeg + 1); j++) {
      for (int k = ksg; k <= (keg + 1); k++) {
        dx1[i][j][k] = x1[i + 1][j][k] - x1[i][j][k];
        geo[edge_idir][i][j][k] = dx1[i][j][k];
      }
    }
  }
  // (ni+1) x nj x (nk+1)
  // dx2 and edge_jdir
  for (int i = isg; i <= (ieg + 1); i++) {
    for (int j = jsg; j <= (jeg); j++) {
      for (int k = ksg; k <= (keg + 1); k++) {
        dx2[i][j][k] = x2[i][j + 1][k] - x2[i][j][k];
        geo[edge_jdir][i][j][k] = x1[i][j][k]*dx2[i][j][k];
      }
    }
  }
  // (ni+1) x (nj+1) x nk
  // dx3 and edge_kdir
  for (int i = isg; i <= (ieg + 1); i++) {
    for (int j = jsg; j <= (jeg + 1); j++) {
      for (int k = ksg; k <= (keg); k++) {
        dx3[i][j][k] = x3[i][j][k + 1] - x3[i][j][k];
        geo[edge_kdir][i][j][k] = x1[i][j][k]*sin(x2[i][j][k])*dx3[i][j][k];
      }
    }
  }

  // dxi dxj dxk
  for (int i = isg; i <= (ieg); i++) {
    for (int j = jsg; j <= (jeg); j++) {
      for (int k = ksg; k <= (keg); k++) {
        geo[dxi][i][j][k] = dx1[i][j][k];
        geo[dxj][i][j][k] = x1c[i][j][k]*dx2[i][j][k];
        geo[dxk][i][j][k] = x1c[i][j][k]*sin(x2c[i][j][k])*dx3[i][j][k];
      }
    }
  }

  // Use the Edges to define Faces
  // (ni+1) x nj x nk : face_idir
  for (int i = isg; i <= (ieg + 1); i++) {
    for (int j = jsg; j <= (jeg); j++) {
      for (int k = ksg; k <= (keg); k++) {
        geo[face_idir][i][j][k] =
            -x1[i][j][k] * x1[i][j][k] * (cos(x2[i][j + 1][k]) - cos(x2[i][j][k]))*dx3[i][j][k];
      }
    }
  }
  // ni x (nj+1) x nk : face_jdir
  for (int i = isg; i <= (ieg); i++) {
    for (int j = jsg; j <= (jeg + 1); j++) {
      for (int k = ksg; k <= (keg); k++) {
        geo[face_jdir][i][j][k] =
            0.5*(x1[i+1][j][k]*x1[i+1][j][k]-x1[i][j][k]*x1[i][j][k])*sin(x2[i][j][k])*dx3[i][j][k];
      }
    }
  }
  // ni x nj x (nk+1) : face_kdir
  for (int i = isg; i <= (ieg); i++) {
    for (int j = jsg; j <= (jeg); j++) {
      for (int k = ksg; k <= (keg + 1); k++) {
        geo[face_kdir][i][j][k] =
            0.5*(x1[i+1][j][k]*x1[i+1][j][k]-x1[i][j][k]*x1[i][j][k])*dx2[i][j][k];
      }
    }
  }
  // ni x nj x nk : vol_center
  for (int i = isg; i <= (ieg); i++) {
    for (int j = jsg; j <= (jeg); j++) {
      for (int k = ksg; k <= (keg); k++) {
        geo[vol_center][i][j][k] =
            fabs(1.0/3.0*(pow(x1[i+1][j][k],3)-pow(x1[i][j][k],3))*
                (cos(x2[i][j+1][k])-cos(x2[i][j][k]))*dx3[i][j][k] );
      }
    }
  }
  // (ni+1) x nj x nk : GF_face_idir
  for (int i = isg; i <= (ieg + 1); i++) {
    for (int j = jsg; j <= (jeg); j++) {
      for (int k = ksg; k <= (keg); k++) {
        geo[GF_face_idir][i][j][k] =
            x1[i][j][k] * x1[i][j][k];
      }
    }
  }
  // ni x (nj+1) x nk : GF_face_jdir
  for (int i = isg; i <= (ieg); i++) {
    for (int j = jsg; j <= (jeg + 1); j++) {
      for (int k = ksg; k <= (keg); k++) {
        geo[GF_face_jdir][i][j][k] =
            sin(x2[i][j][k]);
      }
    }
  }
  // ni x nj x (nk+1) : GF_face_kdir
  for (int i = isg; i <= (ieg); i++) {
    for (int j = jsg; j <= (jeg); j++) {
      for (int k = ksg; k <= (keg + 1); k++) {
        geo[GF_face_kdir][i][j][k] = 1.0;
      }
    }
  }
  // ni x nj x nk : GF_vol
  for (int i = isg; i <= (ieg); i++) {
    for (int j = jsg; j <= (jeg); j++) {
      for (int k = ksg; k <= (keg); k++) {
        geo[GF_vol_idir][i][j][k] =
            1.0/3.0*(pow(x1[i+1][j][k],3)-pow(x1[i][j][k],3));
        geo[GF_vol_jdir][i][j][k] =
            -(cos(x2[i][j+1][k])-cos(x2[i][j][k]));
        geo[GF_vol_kdir][i][j][k] = dx3[i][j][k];
      }
    }
  }

  // x1ctr x2ctr x3ctr
//  for (int i = isg; i <= ieg; i++) {
//    for (int j = jsg; j <= jeg; j++) {
//      for (int k = ksg; k <= keg; k++) {
//        x1ctr[i][j][k] =
//            x1c[i][j][k] + 2.0 * x1c[i][j][k] * dx1[i][j][k] * dx1[i][j][k]
//                / (12.0*x1c[i][j][k]*x1c[i][j][k]+dx1[i][j][k]*dx1[i][j][k]);
//        x2ctr[i][j][k] =
//            ( sin(x2[i][j+1][k]) - sin(x2[i][j][k]) - (x2[i][j+1][k]* cos(x2[i][j+1][k]) - x2[i][j][k]* cos(x2[i][j][k])) )
//                / (-cos(x2[i][j+1][k]) + cos(x2[i][j][k]));
//        x3ctr[i][j][k] = x3c[i][j][k];
//      }
//    }
//  }
  // x1ctr
  for (int i = isg; i <= ieg; i++) {
    for (int j = jsg; j <= jeg+1; j++) {
      for (int k = ksg; k <= keg+1; k++) {
        x1ctr[i][j][k] =
            x1c[i][j][k] + 2.0 * x1c[i][j][k] * dx1[i][j][k] * dx1[i][j][k]
                / (12.0*x1c[i][j][k]*x1c[i][j][k]+dx1[i][j][k]*dx1[i][j][k]);
      }
    }
  }

  // x2ctr
  for (int i = isg; i <= ieg+1; i++) {
    for (int j = jsg; j <= jeg; j++) {
      for (int k = ksg; k <= keg+1; k++) {
        x2ctr[i][j][k] =
            ( sin(x2[i][j+1][k]) - sin(x2[i][j][k]) - (x2[i][j+1][k]* cos(x2[i][j+1][k]) - x2[i][j][k]* cos(x2[i][j][k])) )
                / (-cos(x2[i][j+1][k]) + cos(x2[i][j][k]));
      }
    }
  }

  // x3ctr
  for (int i = isg; i <= ieg+1; i++) {
    for (int j = jsg; j <= jeg+1; j++) {
      for (int k = ksg; k <= keg; k++) {
        x3ctr[i][j][k] = x3c[i][j][k];
      }
    }
  }

  // (ni+1) x nj x nk : vol_idir
  for (int i = is; i <= ie + 1; i++) {
    for (int j = jsg; j <= jeg; j++) {
      for (int k = ksg; k <= keg; k++) {
        double sum = 0.0;
        for (int onface = -4; onface <= 3; onface++) {
          // Adjust fac8th index for C
          sum += geo[vol_center][i + onface][j][k] * fac8th[onface + 4];
        }
        geo[vol_idir][i][j][k] = sum;
      }
    }
  }
  // ni x (nj+1) x nk : vol_jdir
  for (int i = isg; i <= ieg; i++) {
    for (int j = js; j <= je + 1; j++) {
      for (int k = ksg; k <= keg; k++) {
        double sum = 0.0;
        for (int onface = -4; onface <= 3; onface++) {
          // Adjust fac8th index for C
          sum += geo[vol_center][i][j + onface][k] * fac8th[onface + 4];
        }
        geo[vol_jdir][i][j][k] = sum;
      }
    }
  }
  // ni x nj x (nk+1) : vol_kdir
  for (int i = isg; i <= ieg; i++) {
    for (int j = jsg; j <= jeg; j++) {
      for (int k = ks; k <= ke + 1; k++) {
        double sum = 0.0;
        for (int onface = -4; onface <= 3; onface++) {
          // Adjust fac8th index for C
          sum += geo[vol_center][i][j][k + onface] * fac8th[onface + 4];
        }
        geo[vol_kdir][i][j][k] = sum;
      }
    }
  }
  // (ni+1) x nj x nk : IfaceAedgeK
  for (int i = is; i <= ie + 1; i++) {
    for (int j = js; j <= je + 1; j++) {
      for (int k = ks; k <= ke; k++) {
        double sumI = 0.0, sumJ = 0.0;
        for (int onface = -4; onface <= 3; onface++) {
          sumI += geo[face_idir][i][j + onface][k] * fac8th[onface + 4];
          sumJ += geo[face_jdir][i + onface][j][k] * fac8th[onface + 4];
        }
        geo[IfaceAedgeK][i][j][k] = sumI;
        geo[JfaceAedgeK][i][j][k] = sumJ;
      }
    }
  }
  for (int i = is; i <= ie; i++) {
    for (int j = js; j <= je + 1; j++) {
      for (int k = ks; k <= ke + 1; k++) {
        double sumJ = 0.0, sumK = 0.0;
        for (int onface = -4; onface <= 3; onface++) {
          sumJ += geo[face_jdir][i][j][k + onface] * fac8th[onface + 4];
          sumK += geo[face_kdir][i][j + onface][k] * fac8th[onface + 4];
        }
        geo[JfaceAedgeI][i][j][k] = sumJ;
        geo[KfaceAedgeI][i][j][k] = sumK;
      }
    }
  }
  for (int i = is; i <= ie + 1; i++) {
    for (int j = js; j <= je; j++) {
      for (int k = ks; k <= ke + 1; k++) {
        double sumK = 0.0, sumI = 0.0;
        for (int onface = -4; onface <= 3; onface++) {
          sumK += geo[face_kdir][i + onface][j][k] * fac8th[onface + 4];
          sumI += geo[face_idir][i][j][k + onface] * fac8th[onface + 4];
        }
        geo[KfaceAedgeJ][i][j][k] = sumK;
        geo[IfaceAedgeJ][i][j][k] = sumI;
      }
    }
  }
  return 0;
}

// Calculate reconstruction weights for curvilinear coordinates
void getweights_r_nonuniform(){
  double x1_s[config.ni + 2*NO2+1 + 2*NO2];
  for (int i = isg+NO2; i <= ieg+NO2+1; i++) {
    x1_s[i] = x1[i-NO2][0][0];
  }
  for (int i = isg+NO2-1; i>=0; i--) {
    x1_s[i] = x1_s[i+1] - dx1[0][0][0];
  }
  for (int i = ieg+NO2+2; i <= ieg+2*NO2+1; i++) {
    x1_s[i] = x1_s[i-1] + dx1[ieg][0][0];
  }


  int Order = 5, m=2, iL, iR, signum;;
  double x1_m, x1_p, value;
  iL = -Order/2;
  iR = iL + Order - 1;

#pragma omp parallel private(x1_m, x1_p, value, signum)
  {
    double **A = lu_alloc_matrix(Order);
    double *b = lu_alloc_vector(Order);
    double *x = lu_alloc_vector(Order);
    int *perm = lu_alloc_permutation(Order);

    // left side
#pragma omp for schedule(dynamic)
    for (int idx = isg; idx <= ieg+1; idx++) {
      // create matrix A and vector b
      for (int i = iL; i <= iR; i++) {
        x1_m = x1_s[idx+NO2+i-1];
        x1_p = x1_s[idx+NO2+i+0];

        for (int j = 0; j <= Order-1; j++) {
          value = (double)(m+1)/(double)(j+m+1)*
              (pow(x1_p, j+m+1) - pow(x1_m, j+m+1))/
              (pow(x1_p, m+1) - pow(x1_m, m+1));
          A[j][i+abs(iL)] = value;
        }

        value = pow(x1_s[idx+NO2+0], i+abs(iL));
        b[i+abs(iL)] = value;
      }

      // LU decomposition and solve
      if (lu_decomp(A, perm, Order) == 0) {
        lu_solve(A, perm, b, x, Order);
      }else {
        for (int jj = jsg; jj <= jeg + 1; jj++) {
          for (int kk = ksg; kk <= keg + 1; kk++) {
            rec[w1l_first_direction][idx][jj][kk] = 0.0;
            rec[w2l_first_direction][idx][jj][kk] = 0.0;
            rec[w3l_first_direction][idx][jj][kk] = 0.0;
            rec[w4l_first_direction][idx][jj][kk] = 0.0;
            rec[w5l_first_direction][idx][jj][kk] = 0.0;
            rec[w6l_first_direction][idx][jj][kk] = 0.0;
            rec[w7l_first_direction][idx][jj][kk] = 0.0;
          }
        }
        continue; // skip to next idx
      }
      for (int jj = jsg; jj <= jeg+1; jj++) {
        for (int kk = ksg; kk <= keg+1; kk++) {
          rec[w1l_first_direction][idx][jj][kk] = 0.0;
          rec[w2l_first_direction][idx][jj][kk] = x[0];
          rec[w3l_first_direction][idx][jj][kk] = x[1];
          rec[w4l_first_direction][idx][jj][kk] = x[2];
          rec[w5l_first_direction][idx][jj][kk] = x[3];
          rec[w6l_first_direction][idx][jj][kk] = x[4];
          rec[w7l_first_direction][idx][jj][kk] = 0.0;
        }
      }
    }

#pragma omp for schedule(dynamic)
    for (int idx = isg; idx <= ieg+1; idx++) {
      for (int i = iL; i <= iR; i++) {
        x1_m = x1_s[idx+NO2+i-0];
        x1_p = x1_s[idx+NO2+i+1];
        for (int j = 0; j <= Order-1; j++) {
          value = (double)(m+1)/(double)(j+m+1)*(pow(x1_p, j+m+1) - pow(x1_m, j+m+1))/(pow(x1_p, m+1) - pow(x1_m, m+1));
          A[j][i+abs(iL)] = value;
        }

        value = pow(x1_s[idx+NO2+0], i+abs(iL));
        b[i+abs(iL)] = value;
      }

      // Solve Ax = b (LU decomposition)
      if (lu_decomp(A, perm, Order) == 0) {
        lu_solve(A, perm, b, x, Order);
      }else {
        for (int jj = jsg; jj <= jeg + 1; jj++) {
          for (int kk = ksg; kk <= keg + 1; kk++) {
            rec[w1r_first_direction][idx][jj][kk] = 0.0;
            rec[w2r_first_direction][idx][jj][kk] = 0.0;
            rec[w3r_first_direction][idx][jj][kk] = 0.0;
            rec[w4r_first_direction][idx][jj][kk] = 0.0;
            rec[w5r_first_direction][idx][jj][kk] = 0.0;
            rec[w6r_first_direction][idx][jj][kk] = 0.0;
            rec[w7r_first_direction][idx][jj][kk] = 0.0;
          }
        }
        continue; // skip to next idx
      }
      for (int jj = jsg; jj <= jeg+1; jj++) {
        for (int kk = ksg; kk <= keg+1; kk++) {
          rec[w1r_first_direction][idx][jj][kk] = 0.0;
          rec[w2r_first_direction][idx][jj][kk] = x[0];
          rec[w3r_first_direction][idx][jj][kk] = x[1];
          rec[w4r_first_direction][idx][jj][kk] = x[2];
          rec[w5r_first_direction][idx][jj][kk] = x[3];
          rec[w6r_first_direction][idx][jj][kk] = x[4];
          rec[w7r_first_direction][idx][jj][kk] = 0.0;
        }
      }
    }
    lu_free_matrix(A, Order);
    lu_free_vector(b);
    lu_free_vector(x);
    lu_free_permutation(perm);
  }
}

void getweights_theta(){

  double x2_s[config.nj + 2*NO2+1 + 2*NO2];
  for (int j = jsg+NO2; j <= jeg+NO2+1; j++) {
    x2_s[j] = x2[0][j-NO2][0];
  }
  for (int j = jsg+NO2-1; j>=0; j--) {
    x2_s[j] = x2_s[j+1] - dx2[0][0][0];
  }
  for (int j = jeg+NO2+2; j <= jeg+2*NO2+1; j++) {
    x2_s[j] = x2_s[j-1] + dx2[0][jeg][0];
  }

  int Order = 5, m=2, iL, iR, signum;
  double x2_m, x2_p, value, sum1=0;
  iL = -Order/2;
  iR = iL + Order - 1;

#pragma omp parallel private(x2_m, x2_p, value, sum1, signum)
  {
    double **A = lu_alloc_matrix(Order);
    double *b = lu_alloc_vector(Order);
    double *x = lu_alloc_vector(Order);
    int *perm = lu_alloc_permutation(Order);

    // left side
#pragma omp for schedule(dynamic)
    for (int idx = jsg; idx <= jeg+1; idx++) {
      for (int i = iL; i <= iR; i++) {
        x2_m = x2_s[idx+NO2+i-1];
        x2_p = x2_s[idx+NO2+i+0];
        sum1 = 0.0;
        for (int j = 0; j <= Order-1; j++) {
          for (int k = 0; k <= j; k++) {
            sum1 = sum1 + factorial(k)*nchoosek(j, k)*(pow(x2_m,j-k)*cos(x2_m+k*PI/2.0) - pow(x2_p,j-k)*cos(x2_p+k*PI/2.0));
          }
          value = sum1/(cos(x2_m)-cos(x2_p));
          sum1 = 0;
          A[j][i+abs(iL)] = value;
        }
        value = pow(x2_s[idx+NO2+0], i+abs(iL));
        b[i+abs(iL)] = value;
      }

      // Solve Ax = b (LU decomposition)
      if (lu_decomp(A, perm, Order) == 0) {
        lu_solve(A, perm, b, x, Order);
      }else {
        for (int jj = jsg; jj <= jeg + 1; jj++) {
          for (int kk = ksg; kk <= keg + 1; kk++) {
            rec[w1l_second_direction][isg][jj][kk] = 0.0;
            rec[w2l_second_direction][isg][jj][kk] = 0.0;
            rec[w3l_second_direction][isg][jj][kk] = 0.0;
            rec[w4l_second_direction][isg][jj][kk] = 0.0;
            rec[w5l_second_direction][isg][jj][kk] = 0.0;
            rec[w6l_second_direction][isg][jj][kk] = 0.0;
            rec[w7l_second_direction][isg][jj][kk] = 0.0;
          }
        }
      }
      for (int ii = isg; ii <= ieg+1; ii++) {
        for (int kk = ksg; kk <= keg+1; kk++) {
          rec[w1l_second_direction][ii][idx][kk] = 0.0;
          rec[w2l_second_direction][ii][idx][kk] = x[0];
          rec[w3l_second_direction][ii][idx][kk] = x[1];
          rec[w4l_second_direction][ii][idx][kk] = x[2];
          rec[w5l_second_direction][ii][idx][kk] = x[3];
          rec[w6l_second_direction][ii][idx][kk] = x[4];
          rec[w7l_second_direction][ii][idx][kk] = 0.0;
        }
      }
    }

    // right side
#pragma omp for schedule(dynamic)
    for (int idx = jsg; idx <= jeg; idx++) {
      for (int i = iL; i <= iR; i++) {
        x2_m = x2_s[idx + NO2 + i - 0];
        x2_p = x2_s[idx + NO2 + i + 1];
        sum1 = 0.0;
        for (int j = 0; j <= Order - 1; j++) {
          for (int k = 0; k <= j; k++) {
            sum1 = sum1 + factorial(k) * nchoosek(j, k) * (pow(x2_m, j - k) * cos(x2_m + k * PI / 2) - pow(x2_p, j - k) * cos(x2_p + k * PI / 2));
          }
          value = sum1 / (cos(x2_m) - cos(x2_p));
          sum1 = 0;
          A[j][i + abs(iL)] = value;
        }

        value = pow(x2_s[idx + NO2 + 0], i + abs(iL));
        b[i + abs(iL)] = value;
      }

      // Solve Ax = b (LU decomposition)
      if (lu_decomp(A, perm, Order) == 0) {
        lu_solve(A, perm, b, x, Order);
      } else {
        for (int jj = jsg; jj <= jeg + 1; jj++) {
          for (int kk = ksg; kk <= keg + 1; kk++) {
            rec[w1r_second_direction][isg][jj][kk] = 0.0;
            rec[w2r_second_direction][isg][jj][kk] = 0.0;
            rec[w3r_second_direction][isg][jj][kk] = 0.0;
            rec[w4r_second_direction][isg][jj][kk] = 0.0;
            rec[w5r_second_direction][isg][jj][kk] = 0.0;
            rec[w6r_second_direction][isg][jj][kk] = 0.0;
            rec[w7r_second_direction][isg][jj][kk] = 0.0;
          }
        }
      }
      for (int ii = isg; ii <= ieg+1; ii++) {
        for (int kk = ksg; kk <= keg+1; kk++) {
          rec[w1r_second_direction][ii][idx][kk] = 0.0;
          rec[w2r_second_direction][ii][idx][kk] = x[0];
          rec[w3r_second_direction][ii][idx][kk] = x[1];
          rec[w4r_second_direction][ii][idx][kk] = x[2];
          rec[w5r_second_direction][ii][idx][kk] = x[3];
          rec[w6r_second_direction][ii][idx][kk] = x[4];
          rec[w7r_second_direction][ii][idx][kk] = 0.0;
        }
      }
    }

    lu_free_matrix(A, Order);
    lu_free_vector(b);
    lu_free_vector(x);
    lu_free_permutation(perm);
  }
}

void getweights_phi(){
  for (int ii = isg; ii <= ieg+1; ii++) {
    for (int jj = jsg; jj <= jeg + 1; jj++) {
      for (int kk = ksg; kk <= keg + 1; kk++) {
        rec[w1l_third_direction][ii][jj][kk] = 0.0;
        rec[w2l_third_direction][ii][jj][kk] = 1.0 / 30.0;
        rec[w3l_third_direction][ii][jj][kk] = -13.0 / 60.0;
        rec[w4l_third_direction][ii][jj][kk] = 47.0 / 60.0;
        rec[w5l_third_direction][ii][jj][kk] = 9.0 / 20.0;
        rec[w6l_third_direction][ii][jj][kk] = -1.0 / 20.0;
        rec[w7l_third_direction][ii][jj][kk] = 0.0;

        rec[w1r_third_direction][ii][jj][kk] = 0.0;
        rec[w2r_third_direction][ii][jj][kk] = -1.0 / 20.0;
        rec[w3r_third_direction][ii][jj][kk] = 9.0 / 20.0;
        rec[w4r_third_direction][ii][jj][kk] = 47.0 / 60.0;
        rec[w5r_third_direction][ii][jj][kk] = -13.0 / 60.0;
        rec[w6r_third_direction][ii][jj][kk] = 1.0 / 30.0;
        rec[w7r_third_direction][ii][jj][kk] = 0.0;
      }
    }
  }
}

void getweights(){
  // r direction
  getweights_r_nonuniform();
  // theta direction
  getweights_theta();
  // phi direction
  getweights_phi();
}

void Source_Term_i() {
  double M22_hd, M33_hd, M22_b, M33_b, Bsq, S2_HD_1, S2_HD_2, S2_B_1, S2_B_2;
  #pragma omp parallel for collapse(4) private(M22_hd, M33_hd, S2_HD_1, S2_HD_2) schedule(static)
  for (int s = 0 ; s < NS; s++) {
    for (int i = is ; i <= ie ; i++) {
      for (int j = js; j <= je; j++) {
        for (int k = ks; k <= ke; k++) {
          // R source term
          M22_hd = gas[s][gas_p_h][i][j][k] + gas[s][gas_rho_h][i][j][k]* gas[s][gas_v2_h][i][j][k]*gas[s][gas_v2_h][i][j][k];
          M33_hd = gas[s][gas_p_h][i][j][k] + gas[s][gas_rho_h][i][j][k]* gas[s][gas_v3_h][i][j][k]*gas[s][gas_v3_h][i][j][k];

          gas[s][S1_HD][i][j][k] = 0.5*(M22_hd + M33_hd)*(geo[face_idir][i+1][j][k] - geo[face_idir][i][j][k])/geo[vol_center][i][j][k];

          // Theta source term
          S2_HD_1 = -(gas[s][gascf_rhov2][i][j][k]*geo[face_idir][i][j][k] + gas[s][gascf_rhov2][i+1][j][k]*geo[face_idir][i+1][j][k])*dx1[i][j][k]/x1c[i][j][k]/2.0/geo[vol_center][i][j][k];
          S2_HD_2 = (geo[face_jdir][i][j+1][k] - geo[face_jdir][i][j][k])/geo[vol_center][i][j][k]*M33_hd;
          gas[s][S2_HD][i][j][k] = S2_HD_1 + S2_HD_2;

          // Phi source term (part 1)
          gas[s][S3_HD][i][j][k] = -(gas[s][gascf_rhov3][i][j][k]*geo[face_idir][i][j][k] + gas[s][gascf_rhov3][i+1][j][k]*geo[face_idir][i+1][j][k])*dx1[i][j][k]/x1c[i][j][k]/2.0/geo[vol_center][i][j][k];
        }
      }
    }
  }

  #pragma omp parallel for collapse(3) private(M22_b, M33_b, Bsq, S2_B_1, S2_B_2) schedule(static)
  for (int i = is ; i <= ie ; i++) {
    for (int j = js; j <= je; j++) {
      for (int k = ks; k <= ke; k++) {
        // R source term
        Bsq = gem[mags_b1_h][i][j][k]*gem[mags_b1_h][i][j][k] + gem[mags_b2_h][i][j][k]*gem[mags_b2_h][i][j][k] + gem[mags_b3_h][i][j][k]*gem[mags_b3_h][i][j][k] +
            2*(gem[mags_b1_h][i][j][k]*gem[B0_b1][i][j][k] + gem[mags_b2_h][i][j][k]*gem[B0_b2][i][j][k] + gem[mags_b3_h][i][j][k]*gem[B0_b3][i][j][k]);
        M22_b = (0.5*Bsq - (gem[mags_b2_h][i][j][k]+gem[B0_b2][i][j][k])*(gem[mags_b2_h][i][j][k]) - gem[mags_b2_h][i][j][k]*gem[B0_b2][i][j][k]);
        M33_b = (0.5*Bsq - (gem[mags_b3_h][i][j][k]+gem[B0_b3][i][j][k])*(gem[mags_b3_h][i][j][k]) - gem[mags_b3_h][i][j][k]*gem[B0_b3][i][j][k]);
        gem[S1_B][i][j][k] = 0.5*(M22_b + M33_b)*(geo[face_idir][i+1][j][k] - geo[face_idir][i][j][k])/geo[vol_center][i][j][k];

        // Theta source term
        S2_B_1 = -(gem[magsf_x2dir][i][j][k]*geo[face_idir][i][j][k] + gem[magsf_x2dir][i+1][j][k]*geo[face_idir][i+1][j][k])*dx1[i][j][k]/x1c[i][j][k]/2.0/geo[vol_center][i][j][k];
        S2_B_2 = (geo[face_jdir][i][j+1][k] - geo[face_jdir][i][j][k])/geo[vol_center][i][j][k]*M33_b;
        gem[S2_B][i][j][k] = S2_B_1 + S2_B_2;


        // Phi source term (part 1)
        gem[S3_B][i][j][k] = -(gem[magsf_x3dir][i][j][k]*geo[face_idir][i][j][k] + gem[magsf_x3dir][i+1][j][k]*geo[face_idir][i+1][j][k])*dx1[i][j][k]/x1c[i][j][k]/2.0/geo[vol_center][i][j][k];

      }
    }
  }

}

void Source_Term_j() {
  #pragma omp parallel for collapse(4) schedule(static)
  for (int s = 0 ; s < NS; s++) {
    for (int i = is ; i <= ie ; i++) {
      for (int j = js; j <= je; j++) {
        for (int k = ks; k <= ke; k++) {
          // Phi source term (part 2)
          gas[s][S3_HD][i][j][k] += -(gas[s][gascf_rhov3][i][j][k]*geo[face_jdir][i][j][k] + gas[s][gascf_rhov3][i][j+1][k]*geo[face_jdir][i][j+1][k]) *
              (sin(x2[i][j+1][k]) - sin(x2[i][j][k]))/(sin(x2[i][j+1][k]) + sin(x2[i][j][k]))/geo[vol_center][i][j][k];
        }
      }
    }
  }

  #pragma omp parallel for collapse(3) schedule(static)
  for (int i = is ; i <= ie ; i++) {
    for (int j = js; j <= je; j++) {
      for (int k = ks; k <= ke; k++) {
        // Phi source term (part 2)
        gem[S3_B][i][j][k] += -(gem[magsf_x3dir][i][j][k]*geo[face_jdir][i][j][k] + gem[magsf_x3dir][i][j+1][k]*geo[face_jdir][i][j+1][k]) *
            (sin(x2[i][j+1][k]) - sin(x2[i][j][k]))/(sin(x2[i][j+1][k]) + sin(x2[i][j][k]))/geo[vol_center][i][j][k];

      }
    }
  }
}

void Source_Term_k() {
}

void Source_Term(dir_t dir){
  switch (dir) {
    case i_dir:
      Source_Term_i();
      break;
    case j_dir:
      Source_Term_j();
      break;
    case k_dir:
      Source_Term_k();
      break;
  }
}

double get_dt_ring() {
  double Vfluid, Vsound, B_tot_sq, Va, Va_eff, VCFL, dtCFL, minDtCFL = INFINITY;
  int i, j, k;

  // Assuming gas has dimensions [4][keg+1][jeg+1][ieg+1] and indexed as
  // gas[variable][k][j][i]
  //  #pragma omp parallel for collapse(4) schedule(static) reduction(min:minDtCFL)
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
              VCFL / geo[dxk_ring][i][j][k]);
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

void cartesian_to_curvilinear_coord(double x, double y, double z,
                                    double* x1, double* x2, double* x3) {
  *x1 = sqrt(x*x + y*y + z*z);
  if (*x1 == 0.0) {
    *x2 = 0.0;
    *x3 = 0.0;
    return;
  }
  *x2 = acos(z / *x1);
  *x3 = atan2(y, x);
  if (*x3 < 0.0) {
    *x3 += 2.0 * PI;
  }
}

void curvilinear_to_cartesian_coord(double x1, double x2, double x3,
                                    double* x, double* y, double* z) {
  *x = x1 * sin(x2) * cos(x3);
  *y = x1 * sin(x2) * sin(x3);
  *z = x1 * cos(x2);
}

void cartesian_to_curvilinear_vector(double x, double y, double z,
                                     double vx, double vy, double vz,
                                     double* v1, double* v2, double* v3) {
  double x1, x2, x3;
  cartesian_to_curvilinear_coord(x, y, z, &x1, &x2, &x3);

  double sin_theta = sin(x2);
  double cos_theta = cos(x2);
  double sin_phi = sin(x3);
  double cos_phi = cos(x3);

  *v1 = vx * sin_theta * cos_phi + vy * sin_theta * sin_phi + vz * cos_theta;
  *v2 = vx * cos_theta * cos_phi + vy * cos_theta * sin_phi - vz * sin_theta;
  *v3 = -vx * sin_phi + vy * cos_phi;
}

void curvilinear_to_cartesian_vector(double x1, double x2, double x3,
                                     double v1, double v2, double v3,
                                     double* vx, double* vy, double* vz) {
  double x, y, z;
  curvilinear_to_cartesian_coord(x1, x2, x3, &x, &y, &z);
  double sin_theta = sin(x2);
  double cos_theta = cos(x2);
  double sin_phi = sin(x3);
  double cos_phi = cos(x3);

  *vx = v1 * sin_theta * cos_phi + v2 * cos_theta * cos_phi - v3 * sin_phi;
  *vy = v1 * sin_theta * sin_phi + v2 * cos_theta * sin_phi + v3 * cos_phi;
  *vz = v1 * cos_theta - v2 * sin_theta;
}