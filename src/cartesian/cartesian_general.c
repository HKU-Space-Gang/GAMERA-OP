// This file contains general functions for cartesian coordinates, including:
// Lame coefficients
// Geometry calculation
// Calculation of reconstruction weights
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
  return 1.0;
}

double H3(double x1, double x2, double x3) {
  return 1.0;
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
        geo[edge_jdir][i][j][k] = dx2[i][j][k];
      }
    }
  }
  // (ni+1) x (nj+1) x nk
  // dx3 and edge_kdir
  for (int i = isg; i <= (ieg + 1); i++) {
    for (int j = jsg; j <= (jeg + 1); j++) {
      for (int k = ksg; k <= (keg); k++) {
        dx3[i][j][k] = x3[i][j][k + 1] - x3[i][j][k];
        geo[edge_kdir][i][j][k] = dx3[i][j][k];
      }
    }
  }

  // dxi dxj dxk
  for (int i = isg; i <= (ieg); i++) {
    for (int j = jsg; j <= (jeg); j++) {
      for (int k = ksg; k <= (keg); k++) {
        geo[dxi][i][j][k] = dx1[i][j][k];
        geo[dxj][i][j][k] = dx2[i][j][k];
        geo[dxk][i][j][k] = dx3[i][j][k];
      }
    }
  }

  // Use the Edges to define Faces
  // (ni+1) x nj x nk : face_idir
  for (int i = isg; i <= (ieg + 1); i++) {
    for (int j = jsg; j <= (jeg); j++) {
      for (int k = ksg; k <= (keg); k++) {
        geo[face_idir][i][j][k] =
            dx2[i][j][k]*dx3[i][j][k];
      }
    }
  }
  // ni x (nj+1) x nk : face_jdir
  for (int i = isg; i <= (ieg); i++) {
    for (int j = jsg; j <= (jeg + 1); j++) {
      for (int k = ksg; k <= (keg); k++) {
        geo[face_jdir][i][j][k] =
            dx1[i][j][k]*dx3[i][j][k];
      }
    }
  }
  // ni x nj x (nk+1) : face_kdir
  for (int i = isg; i <= (ieg); i++) {
    for (int j = jsg; j <= (jeg); j++) {
      for (int k = ksg; k <= (keg + 1); k++) {
        geo[face_kdir][i][j][k] =
            dx1[i][j][k]*dx2[i][j][k];
      }
    }
  }
  // ni x nj x nk : vol_center
  for (int i = isg; i <= (ieg); i++) {
    for (int j = jsg; j <= (jeg); j++) {
      for (int k = ksg; k <= (keg); k++) {
        geo[vol_center][i][j][k] =
            dx1[i][j][k]*dx2[i][j][k]*dx3[i][j][k];
      }
    }
  }
  // (ni+1) x nj x nk : GF_face_idir
  for (int i = isg; i <= (ieg + 1); i++) {
    for (int j = jsg; j <= (jeg); j++) {
      for (int k = ksg; k <= (keg); k++) {
        geo[GF_face_idir][i][j][k] =
            geo[face_idir][i][j][k];
      }
    }
  }
  // ni x (nj+1) x nk : GF_face_jdir
  for (int i = isg; i <= (ieg); i++) {
    for (int j = jsg; j <= (jeg + 1); j++) {
      for (int k = ksg; k <= (keg); k++) {
        geo[GF_face_jdir][i][j][k] =
            geo[face_jdir][i][j][k];
      }
    }
  }
  // ni x nj x (nk+1) : GF_face_kdir
  for (int i = isg; i <= (ieg); i++) {
    for (int j = jsg; j <= (jeg); j++) {
      for (int k = ksg; k <= (keg + 1); k++) {
        geo[GF_face_kdir][i][j][k] = geo[face_kdir][i][j][k];
      }
    }
  }
  // ni x nj x nk : GF_vol
  for (int i = isg; i <= (ieg); i++) {
    for (int j = jsg; j <= (jeg); j++) {
      for (int k = ksg; k <= (keg); k++) {
        geo[GF_vol_idir][i][j][k] =
            geo[vol_center][i][j][k];
        geo[GF_vol_jdir][i][j][k] =
            geo[vol_center][i][j][k];
        geo[GF_vol_kdir][i][j][k] = geo[vol_center][i][j][k];
      }
    }
  }

  // x1ctr x2ctr x3ctr
  for (int i = isg; i <= ieg; i++) {
    for (int j = jsg; j <= jeg; j++) {
      for (int k = ksg; k <= keg; k++) {
        x1ctr[i][j][k] =
            x1c[i][j][k];
        x2ctr[i][j][k] =
            x2c[i][j][k];
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
void getweights_x1(){
  for (int ii = is; ii <= ie + 1; ii++) {
    for (int jj = jsg; jj <= jeg + 1; jj++) {
      for (int kk = ksg; kk <= keg + 1; kk++) {
        // double facPDMU7[7] = {-1.0/140, 5.0/84, -101.0/420, 319.0/420, 107.0/210,
        // -19.0/210, 1.0/105};
        rec[w1l_first_direction][ii][jj][kk] = -1.0/140.0 * geo[vol_center][ii][jj][kk] / geo[vol_idir][ii][jj][kk];
        rec[w2l_first_direction][ii][jj][kk] = 5.0/84.0 * geo[vol_center][ii][jj][kk] / geo[vol_idir][ii][jj][kk];
        rec[w3l_first_direction][ii][jj][kk] = -101.0/420.0 * geo[vol_center][ii][jj][kk] / geo[vol_idir][ii][jj][kk];
        rec[w4l_first_direction][ii][jj][kk] = 319.0/420.0 * geo[vol_center][ii][jj][kk] / geo[vol_idir][ii][jj][kk];
        rec[w5l_first_direction][ii][jj][kk] = 107.0/210.0 * geo[vol_center][ii][jj][kk] / geo[vol_idir][ii][jj][kk];
        rec[w6l_first_direction][ii][jj][kk] = -19.0/210.0 * geo[vol_center][ii][jj][kk] / geo[vol_idir][ii][jj][kk];
        rec[w7l_first_direction][ii][jj][kk] = 1.0/105.0 * geo[vol_center][ii][jj][kk] / geo[vol_idir][ii][jj][kk];

        rec[w1r_first_direction][ii][jj][kk] = 1.0/105.0 * geo[vol_center][ii][jj][kk] / geo[vol_idir][ii][jj][kk];
        rec[w2r_first_direction][ii][jj][kk] = -19.0/210.0 * geo[vol_center][ii][jj][kk] / geo[vol_idir][ii][jj][kk];
        rec[w3r_first_direction][ii][jj][kk] = 107.0/210.0 * geo[vol_center][ii][jj][kk] / geo[vol_idir][ii][jj][kk];
        rec[w4r_first_direction][ii][jj][kk] = 319.0/420.0 * geo[vol_center][ii][jj][kk] / geo[vol_idir][ii][jj][kk];
        rec[w5r_first_direction][ii][jj][kk] = -101.0/420.0 * geo[vol_center][ii][jj][kk] / geo[vol_idir][ii][jj][kk];
        rec[w6r_first_direction][ii][jj][kk] = 5.0/84.0 * geo[vol_center][ii][jj][kk] / geo[vol_idir][ii][jj][kk];
        rec[w7r_first_direction][ii][jj][kk] = -1.0/140.0 * geo[vol_center][ii][jj][kk] / geo[vol_idir][ii][jj][kk];
      }
    }
  }
}

void getweights_x2(){
  for (int ii = isg; ii <= ieg+1; ii++) {
    for (int jj = js; jj <= je + 1; jj++) {
      for (int kk = ksg; kk <= keg + 1; kk++) {
        rec[w1l_second_direction][ii][jj][kk] = -1.0/140.0 * geo[vol_center][ii][jj][kk] / geo[vol_jdir][ii][jj][kk];
        rec[w2l_second_direction][ii][jj][kk] = 5.0/84.0 * geo[vol_center][ii][jj][kk] / geo[vol_jdir][ii][jj][kk];
        rec[w3l_second_direction][ii][jj][kk] = -101.0/420.0 * geo[vol_center][ii][jj][kk] / geo[vol_jdir][ii][jj][kk];
        rec[w4l_second_direction][ii][jj][kk] = 319.0/420.0 * geo[vol_center][ii][jj][kk] / geo[vol_jdir][ii][jj][kk];
        rec[w5l_second_direction][ii][jj][kk] = 107.0/210.0 * geo[vol_center][ii][jj][kk] / geo[vol_jdir][ii][jj][kk];
        rec[w6l_second_direction][ii][jj][kk] = -19.0/210.0 * geo[vol_center][ii][jj][kk] / geo[vol_jdir][ii][jj][kk];
        rec[w7l_second_direction][ii][jj][kk] = 1.0/105.0 * geo[vol_center][ii][jj][kk] / geo[vol_jdir][ii][jj][kk];

        rec[w1r_second_direction][ii][jj][kk] = 1.0/105.0 * geo[vol_center][ii][jj][kk] / geo[vol_jdir][ii][jj][kk];
        rec[w2r_second_direction][ii][jj][kk] = -19.0/210.0 * geo[vol_center][ii][jj][kk] / geo[vol_jdir][ii][jj][kk];
        rec[w3r_second_direction][ii][jj][kk] = 107.0/210.0 * geo[vol_center][ii][jj][kk] / geo[vol_jdir][ii][jj][kk];
        rec[w4r_second_direction][ii][jj][kk] = 319.0/420.0 * geo[vol_center][ii][jj][kk] / geo[vol_jdir][ii][jj][kk];
        rec[w5r_second_direction][ii][jj][kk] = -101.0/420.0 * geo[vol_center][ii][jj][kk] / geo[vol_jdir][ii][jj][kk];
        rec[w6r_second_direction][ii][jj][kk] = 5.0/84.0 * geo[vol_center][ii][jj][kk] / geo[vol_jdir][ii][jj][kk];
        rec[w7r_second_direction][ii][jj][kk] = -1.0/140.0 * geo[vol_center][ii][jj][kk] / geo[vol_jdir][ii][jj][kk];
      }
    }
  }
}

void getweights_x3(){
  for (int ii = isg; ii <= ieg+1; ii++) {
    for (int jj = jsg; jj <= jeg + 1; jj++) {
      for (int kk = ks; kk <= ke + 1; kk++) {
        rec[w1l_third_direction][ii][jj][kk] = -1.0/140.0 * geo[vol_center][ii][jj][kk] / geo[vol_kdir][ii][jj][kk];
        rec[w2l_third_direction][ii][jj][kk] = 5.0/84.0 * geo[vol_center][ii][jj][kk] / geo[vol_kdir][ii][jj][kk];
        rec[w3l_third_direction][ii][jj][kk] = -101.0/420.0 * geo[vol_center][ii][jj][kk] / geo[vol_kdir][ii][jj][kk];
        rec[w4l_third_direction][ii][jj][kk] = 319.0/420.0 * geo[vol_center][ii][jj][kk] / geo[vol_kdir][ii][jj][kk];
        rec[w5l_third_direction][ii][jj][kk] = 107.0/210.0 * geo[vol_center][ii][jj][kk] / geo[vol_kdir][ii][jj][kk];
        rec[w6l_third_direction][ii][jj][kk] = -19.0/210.0 * geo[vol_center][ii][jj][kk] / geo[vol_kdir][ii][jj][kk];
        rec[w7l_third_direction][ii][jj][kk] = 1.0/105.0 * geo[vol_center][ii][jj][kk] / geo[vol_kdir][ii][jj][kk];

        rec[w1r_third_direction][ii][jj][kk] = 1.0/105.0 * geo[vol_center][ii][jj][kk] / geo[vol_kdir][ii][jj][kk];
        rec[w2r_third_direction][ii][jj][kk] = -19.0/210.0 * geo[vol_center][ii][jj][kk] / geo[vol_kdir][ii][jj][kk];
        rec[w3r_third_direction][ii][jj][kk] = 107.0/210.0 * geo[vol_center][ii][jj][kk] / geo[vol_kdir][ii][jj][kk];
        rec[w4r_third_direction][ii][jj][kk] = 319.0/420.0 * geo[vol_center][ii][jj][kk] / geo[vol_kdir][ii][jj][kk];
        rec[w5r_third_direction][ii][jj][kk] = -101.0/420.0 * geo[vol_center][ii][jj][kk] / geo[vol_kdir][ii][jj][kk];
        rec[w6r_third_direction][ii][jj][kk] = 5.0/84.0 * geo[vol_center][ii][jj][kk] / geo[vol_kdir][ii][jj][kk];
        rec[w7r_third_direction][ii][jj][kk] = -1.0/140.0 * geo[vol_center][ii][jj][kk] / geo[vol_kdir][ii][jj][kk];
      }
    }
  }
}

void getweights(){
  // x1 direction
  getweights_x1();
  // x2 direction
  getweights_x2();
  // x3 direction
  getweights_x3();
}

void Source_Term_i() {
}

void Source_Term_j() {
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
  log_error("get_dt_ring cannot be used in cartesian coordinates!\n");
  return 0.0;
}

void cartesian_to_curvilinear_coord(double x, double y, double z,
                                    double* x1, double* x2, double* x3) {
  *x1 = x;
  *x2 = y;
  *x3 = z;
}

void curvilinear_to_cartesian_coord(double x1, double x2, double x3,
                                    double* x, double* y, double* z) {
  *x = x1;
  *y = x2;
  *z = x3;
}

void cartesian_to_curvilinear_vector(double x, double y, double z,
                                     double vx, double vy, double vz,
                                     double* v1, double* v2, double* v3) {
  *v1 = vx;
  *v2 = vy;
  *v3 = vz;
}

void curvilinear_to_cartesian_vector(double x1, double x2, double x3,
                                     double v1, double v2, double v3,
                                     double* vx, double* vy, double* vz) {
  *vx = v1;
  *vy = v2;
  *vz = v3;
}