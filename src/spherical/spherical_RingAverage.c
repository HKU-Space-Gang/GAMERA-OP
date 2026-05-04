// This file contains ring average related functions for spherical coordinates, including:
// Initialization, setting geo data, setting config.
// Tackle E and B in pole
// Gather and broadcast ring average arrays

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

void init_RingAverage_variables() {
  // Ring related arrays configuration
  Ring_Ni = config.NI;
  Ring_Nj = config.NJ;
  Ring_Nk = config.nk_global + 2 * NG + 1;
}

void set_RingAverage_geo_data() {
  if (config.proc_dims[2] == 1) {
#pragma omp parallel for collapse(4) schedule(static)
    for (int m = face_idir_ring; m <= face_kdir_ring; m++) {
      for (int i = isg; i <= ieg + 1; i++) {
        for (int j = jsg; j <= jeg + 1; j++) {
          for (int k = ksg; k <= keg + 1; k++) {
            int diff = face_idir - face_idir_ring;
            geo_ring[m][i][j][k] = geo[m + diff][i][j][k];
          }
        }
      }
    }
  } else {
    for (int m = face_idir_ring; m <= face_kdir_ring; m++) {
      int diff = face_idir - face_idir_ring;
      RingAverage_gather_data_3d_array(geo[m + diff], geo_ring[m]);
    }
  }
}


void set_RingAverage_config() {
  if (config.nk_global == 256) {
    int temp[] = {8, 16, 32, 32, 64, 64, 64, 64, 128, 128, 128, 128, 128, 128, 128, 128};
    memcpy(Nchunk, temp, sizeof(temp));
    NAVR = sizeof(temp) / sizeof(temp[0]);
  }
  else if (config.nk_global == 128) {
    int temp[] = {8, 16, 32, 32, 64, 64, 64, 64};
    memcpy(Nchunk, temp, sizeof(temp));
    NAVR = sizeof(temp) / sizeof(temp[0]);
  }
  else if (config.nk_global == 64) {
    int temp[] = {8, 16, 32, 32};
    memcpy(Nchunk, temp, sizeof(temp));
    NAVR = sizeof(temp) / sizeof(temp[0]);
  }
  else if (config.nk_global == 32) {
    int temp[] = {8, 16};
    memcpy(Nchunk, temp, sizeof(temp));
    NAVR = sizeof(temp) / sizeof(temp[0]);
  }
  else{
    printf("config.nk_global must be 32, 64, 128, or 256\n");
    exit(0);
  }

  for (int i = isg; i <= ieg; i++) {
    for (int j = jsg; j <= jeg; j++) {
      for (int k = ksg; k <= (keg ); k++) {
        geo[dxk_ring][i][j][k] = geo[dxk][i][j][k];
      }
    }
  }

  int j_global_min = proc_coords[1] * config.nj_global / config.proc_dims[1];
  int j_global_max = (proc_coords[1] + 1) * config.nj_global / config.proc_dims[1] - 1;
  if (NAVR > j_global_min) {
    doRingAverage_thisrank_positivepole = 1;
  }
  if (NAVR > config.nj_global - j_global_max - 1) {
    doRingAverage_thisrank_negativepole = 1;
  }
  if (NAVR > config.nj_global - j_global_max - 1 && NAVR > j_global_min) {
    log_info("This rank should deal with the positive and negative pole together.");
  }

  // Eliminate the unnecessary ring average
  if (fabs(x2min_global) > 1e-10) {
    doRingAverage_thisrank_positivepole = 0;
  }
  if (fabs(x2max_global-PI) > 1e-10) {
    doRingAverage_thisrank_negativepole = 0;
  }

  log_info("Rank %d: doRingAverage_thisrank_positivepole = %d, doRingAverage_thisrank_negativepole = %d",
           rank, doRingAverage_thisrank_positivepole, doRingAverage_thisrank_negativepole);

  // positive pole
  if (doRingAverage_thisrank_positivepole == 1) {
    int j_max_temp;
    if (j_global_max < NAVR) {
      j_max_temp = config.nj - 1;
    } else {
      j_max_temp = NAVR - j_global_min - 1;
    }
    for (int jring = 0; jring <= j_max_temp; jring++) {
      int Nm = config.nk_global / Nchunk[jring + j_global_min];
      log_info("Rank %d: positive pole jring=%d, Nchunk=%d, Nm=%d",
               rank, jring, Nchunk[jring + j_global_min], Nm);
      for (int i = is; i <= ie; i++) {
        for (int k = ks; k <= ke; k++) {
          geo[dxk_ring][i][js + jring][k] = Nm * geo[dxk][i][js + jring][k];
        }
      }
    }
  }
  // negative pole
  if (doRingAverage_thisrank_negativepole == 1) {
    int j_max_temp;
    if (NAVR > config.nj_global - j_global_min - 1) {
      j_max_temp = config.nj - 1;
    } else {
      j_max_temp = NAVR - (config.nj_global - j_global_max);
    }
    for (int jring = 0; jring <= j_max_temp; jring++) {
      int Nm = config.nk_global / Nchunk[jring + (config.nj_global - j_global_max - 1)];
      log_info("Rank %d: negative pole jring=%d, Nchunk=%d, Nm=%d",
               rank, jring, Nchunk[jring + (config.nj_global - j_global_max - 1)], Nm);
      for (int i = is; i <= ie; i++) {
        for (int k = ks; k <= ke; k++) {
          geo[dxk_ring][i][je - jring][k] = Nm * geo[dxk][i][je - jring][k];
        }
      }
    }
  }

}


// For spherical coordinates, special treatment is needed for E and B fields at poles
void tackle_Efield_pole() {
  int ks_ring = ks;
  int ke_ring = ks + config.nk_global - 1;
  gather_Etackle_RingAverage_arrays();
  if (proc_coords[2] == 0) {
    if (fabs(x2min_global) < 1e-10 && proc_coords[1] == 0) {
      #pragma omp parallel for schedule(static)
      for (int i = is; i <= ie; i++) {
        double Ei_mean1 = 0.0;
        for (int k = ks_ring; k <= ke_ring + 1; k++) {
          Ei_mean1 = Ei_mean1 + gem_ring[efield_tackle_ei_ring][i][js][k];
        }
        Ei_mean1 = Ei_mean1 / (ke_ring - ks_ring + 2);
        for (int k = ks_ring; k <= ke_ring + 1; k++) {
          gem_ring[efield_tackle_ei_ring][i][js][k] = Ei_mean1;
        }
      }

      // E_k at the pole
      #pragma omp parallel for collapse(2) schedule(static)
      for (int i = is; i <= ie + 1; i++) {
        for (int k = ks_ring; k <= ke_ring; k++) {
          gem_ring[efield_tackle_ek_ring][i][js][k] = 0.0;
        }
      }
    }
    if (fabs(x2max_global - PI) < 1e-10 && proc_coords[1] == config.proc_dims[1] - 1) {
      #pragma omp parallel for schedule(static)
      for (int i = is; i <= ie; i++) {
        double Ei_mean2 = 0.0;
        for (int k = ks_ring; k <= ke_ring + 1; k++) {
          Ei_mean2 = Ei_mean2 + gem_ring[efield_tackle_ei_ring][i][je + 1][k];
        }
        Ei_mean2 = Ei_mean2 / (ke_ring - ks_ring + 2);
        for (int k = ks_ring; k <= ke_ring + 1; k++) {
          gem_ring[efield_tackle_ei_ring][i][je + 1][k] = Ei_mean2;
        }
      }

      // E_k at the pole
      #pragma omp parallel for collapse(2) schedule(static)
      for (int i = is; i <= ie + 1; i++) {
        for (int k = ks_ring; k <= ke_ring; k++) {
          gem_ring[efield_tackle_ek_ring][i][je + 1][k] = 0.0;
        }
      }
    }
  }
  broadcast_Etackle_RingAverage_arrays();
}

void tackle_Magfield_pole() {
  int ks_ring = ks;
  int ke_ring = ks + config.nk_global - 1;
  gather_value_RingAverage_arrays();
  if (proc_coords[2] == 0) {
    if (fabs(x2min_global) < 1e-10 && proc_coords[1] == 0) {
      #pragma omp parallel for collapse(2) schedule(static)
      for (int i = is; i <= ie; i++) {
        for (int k = ks_ring; k <= ke_ring; k++) {
//          gem_ring[mag_bj_ring][i][js][k] = gem_ring[mag_bj_ring][i][js + 1][k];
          gem_ring[mag_bj_ring][i][js][k] = 0.5*(gem_ring[mag_bj_ring][i][js + 1][k] -
              gem_ring[mag_bj_ring][i][js + 1][NG + (k - NG + config.nk_global / 2) % config.nk_global]);
        }
      }
    }

    if (fabs(x2max_global - PI) < 1e-10 && proc_coords[1] == config.proc_dims[1] - 1) {
      #pragma omp parallel for collapse(2) schedule(static)
      for (int i = is; i <= ie; i++) {
        for (int k = ks_ring; k <= ke_ring; k++) {
//          gem_ring[mag_bj_ring][i][je + 1][k] = gem_ring[mag_bj_ring][i][je][k];
          gem_ring[mag_bj_ring][i][je + 1][k] = 0.5*(gem_ring[mag_bj_ring][i][je][k] -
              gem_ring[mag_bj_ring][i][je][NG + (k - NG + config.nk_global / 2) % config.nk_global]);
        }
      }
    }
  }
  broadcast_value_RingAverage_arrays();
}

void HydroRingAverage_base()
{
  int j_max_temp = 0;
  int j_global_min = proc_coords[1] * config.nj_global / config.proc_dims[1];
  int j_global_max = (proc_coords[1] + 1) * config.nj_global / config.proc_dims[1] - 1;
  // positive pole
  if (doRingAverage_thisrank_positivepole == 1) {
    if (j_global_max < NAVR) {
      j_max_temp = config.nj - 1;
    } else {
      j_max_temp = NAVR - j_global_min - 1;
    }
  }
  // negative pole
  if (doRingAverage_thisrank_negativepole == 1) {
    if (NAVR > config.nj_global - j_global_min - 1) {
      j_max_temp = config.nj - 1;
    } else {
      j_max_temp = NAVR - (config.nj_global - j_global_max);
    }
  }

  // positive pole
  if (doRingAverage_thisrank_positivepole == 1) {
#pragma omp parallel for collapse(2) schedule(static)
    for (int s = 0; s < NS; s++) {
      for (int i = 0; i < config.ni; i++) {
        int i_idx = is + i;
        int N = config.nk_global;
        int ks_ring = ks;
        int ke_ring = ks + config.nk_global - 1;
        double *Q_fluid = (double *) malloc(N * sizeof(double));
        double *base = (double *) malloc(N * sizeof(double));
        double *rem = (double *) malloc(N * sizeof(double));

        for (int j = 0; j <= j_max_temp; j++) {
          int j_idx = js + j;
          int j_chunk = j + j_global_min;

          // 1. Density (rho)
          for (int k = ks_ring; k <= ke_ring; k++) Q_fluid[k - ks_ring] = gas_ring[s][gas_rho_ring][i_idx][j_idx][k];
          RingBlockAvg(Q_fluid, N, Nchunk[j_chunk]);
          for (int k = ks_ring; k <= ke_ring; k++) gas_ring[s][gas_rho_ring][i_idx][j_idx][k] = Q_fluid[k - ks_ring];

          // 2. Velocity (v1, v2, v3)
          // v1 component
          for (int k = ks_ring; k <= ke_ring; k++) Q_fluid[k - ks_ring] = gas_ring[s][gas_v1_ring][i_idx][j_idx][k];
          FilterBaseMode(Q_fluid, N, base);
          memcpy(rem, Q_fluid, N * sizeof(double));
          RingBlockAvg(rem, N, Nchunk[j_chunk]);
          for (int k = 0; k < N; k++) Q_fluid[k] = rem[k] + base[k];
          for (int k = ks_ring; k <= ke_ring; k++) gas_ring[s][gas_v1_ring][i_idx][j_idx][k] = Q_fluid[k - ks_ring];

          // v2 component (same as v1)
          for (int k = ks_ring; k <= ke_ring; k++) Q_fluid[k - ks_ring] = gas_ring[s][gas_v2_ring][i_idx][j_idx][k];
          FilterBaseMode(Q_fluid, N, base);
          memcpy(rem, Q_fluid, N * sizeof(double));
          RingBlockAvg(rem, N, Nchunk[j_chunk]);
          for (int k = 0; k < N; k++) Q_fluid[k] = rem[k] + base[k];
          for (int k = ks_ring; k <= ke_ring; k++) gas_ring[s][gas_v2_ring][i_idx][j_idx][k] = Q_fluid[k - ks];

          // v3 component (same as v1)
          for (int k = ks_ring; k <= ke_ring; k++) Q_fluid[k - ks_ring] = gas_ring[s][gas_v3_ring][i_idx][j_idx][k];
          FilterBaseMode(Q_fluid, N, base);
          memcpy(rem, Q_fluid, N * sizeof(double));
          RingBlockAvg(rem, N, Nchunk[j_chunk]);
          for (int k = 0; k < N; k++) Q_fluid[k] = rem[k] + base[k];
          for (int k = ks_ring; k <= ke_ring; k++) gas_ring[s][gas_v3_ring][i_idx][j_idx][k] = Q_fluid[k - ks];

          // 3. Pressure (p)
          for (int k = ks_ring; k <= ke_ring; k++) Q_fluid[k - ks] = gas_ring[s][gas_p_ring][i_idx][j_idx][k];
          RingBlockAvg(Q_fluid, N, Nchunk[j_chunk]);
          for (int k = ks_ring; k <= ke_ring; k++) gas_ring[s][gas_p_ring][i_idx][j_idx][k] = Q_fluid[k - ks_ring];
        }

        free(Q_fluid);
        free(base);
        free(rem);
      }
    }
  }

  // negative pole
  if (doRingAverage_thisrank_negativepole == 1) {
#pragma omp parallel for collapse(2) schedule(static)
    for (int s = 0; s < NS; s++) {
      for (int i = 0; i < config.ni; i++) {
        int i_idx = is + i;
        int N = config.nk_global;
        int ks_ring = ks;
        int ke_ring = ks + config.nk_global - 1;
        double *Q_fluid = (double *) malloc(N * sizeof(double));
        double *base = (double *) malloc(N * sizeof(double));
        double *rem = (double *) malloc(N * sizeof(double));

        for (int j = 0; j <= j_max_temp; j++) {
          int j_rev = je - j;
          int j_chunk = j + (config.nj_global - j_global_max - 1);

          // Density (rho)
          for (int k = ks_ring; k <= ke_ring; k++) Q_fluid[k-ks_ring] = gas_ring[s][gas_rho_ring][i_idx][j_rev][k];
          RingBlockAvg(Q_fluid, N, Nchunk[j_chunk]);
          for (int k = ks_ring; k <= ke_ring; k++) gas_ring[s][gas_rho_ring][i_idx][j_rev][k] = Q_fluid[k-ks_ring];

          // Velocity (v1, v2, v3)
          // v1 component
          for (int k = ks_ring; k <= ke_ring; k++) Q_fluid[k-ks_ring] = gas_ring[s][gas_v1_ring][i_idx][j_rev][k];
          FilterBaseMode(Q_fluid, N, base);
          memcpy(rem, Q_fluid, N * sizeof(double));
          RingBlockAvg(rem, N, Nchunk[j_chunk]);
          for (int k = 0; k < N; k++) Q_fluid[k] = rem[k] + base[k];
          for (int k = ks_ring; k <= ke_ring; k++) gas_ring[s][gas_v1_ring][i_idx][j_rev][k] = Q_fluid[k-ks_ring];

          // v2 component (same as v1)
          for (int k = ks_ring; k <= ke_ring; k++) Q_fluid[k-ks_ring] = gas_ring[s][gas_v2_ring][i_idx][j_rev][k];
          FilterBaseMode(Q_fluid, N, base);
          memcpy(rem, Q_fluid, N * sizeof(double));
          RingBlockAvg(rem, N, Nchunk[j_chunk]);
          for (int k = 0; k < N; k++) Q_fluid[k] = rem[k] + base[k];
          for (int k = ks_ring; k <= ke_ring; k++) gas_ring[s][gas_v2_ring][i_idx][j_rev][k] = Q_fluid[k-ks_ring];

          // v3 component (same as v1)
          for (int k = ks_ring; k <= ke_ring; k++) Q_fluid[k-ks_ring] = gas_ring[s][gas_v3_ring][i_idx][j_rev][k];
          FilterBaseMode(Q_fluid, N, base);
          memcpy(rem, Q_fluid, N * sizeof(double));
          RingBlockAvg(rem, N, Nchunk[j_chunk]);
          for (int k = 0; k < N; k++) Q_fluid[k] = rem[k] + base[k];
          for (int k = ks_ring; k <= ke_ring; k++) gas_ring[s][gas_v3_ring][i_idx][j_rev][k] = Q_fluid[k-ks_ring];

          // 3. Pressure (p)
          for (int k = ks_ring; k <= ke_ring; k++) Q_fluid[k-ks_ring] = gas_ring[s][gas_p_ring][i_idx][j_rev][k];
          RingBlockAvg(Q_fluid, N, Nchunk[j_chunk]);
          for (int k = ks_ring; k <= ke_ring; k++) gas_ring[s][gas_p_ring][i_idx][j_rev][k] = Q_fluid[k-ks_ring];
        }

        free(Q_fluid);
        free(base);
        free(rem);
      }
    }
  }
}


void MagneticRingAverage_base()
{
  int ksg_ring = ksg;
  int keg_ring = Ring_Nk - 2;
#pragma omp parallel for collapse(3) schedule(static)
  for (int i = isg; i <= ieg + 1; i++) {
    for (int j = jsg; j <= jeg + 1; j++) {
      for (int k = ksg_ring; k <= keg_ring + 1; k++) {
        gem_ring[efield_ei_ring][i][j][k] = 0.0;
        gem_ring[efield_ej_ring][i][j][k] = 0.0;
        gem_ring[efield_ek_ring][i][j][k] = 0.0;
      }
    }
  }

#pragma omp parallel for collapse(3) schedule(static)
  for (int i = isg; i <= ieg+1 ; i++) {
    for (int j = jsg; j <= jeg; j++) {
      for (int k = ksg_ring; k <= keg_ring; k++) {
        gem_ring[mag_bi_flux_ring][i][j][k] = gem_ring[mag_bi_ring][i][j][k] * geo_ring[face_idir_ring][i][j][k];
      }
    }
  }

#pragma omp parallel for collapse(3) schedule(static)
  for (int i = isg; i <= ieg; i++) {
    for (int j = jsg; j <= jeg+1; j++) {
      for (int k = ksg_ring; k <= keg_ring; k++) {
        gem_ring[mag_bj_flux_ring][i][j][k] = gem_ring[mag_bj_ring][i][j][k] * geo_ring[face_jdir_ring][i][j][k];
      }
    }
  }

  int j_max_temp = 0;
  int j_global_min = proc_coords[1] * config.nj_global / config.proc_dims[1];
  int j_global_max = (proc_coords[1] + 1) * config.nj_global / config.proc_dims[1] - 1;
  // positive pole
  if (doRingAverage_thisrank_positivepole == 1) {
    if (j_global_max < NAVR) {
      j_max_temp = config.nj - 1;
    } else {
      j_max_temp = NAVR - j_global_min - 1;
    }
  }
  // negative pole
  if (doRingAverage_thisrank_negativepole == 1) {
    if (NAVR > config.nj_global - j_global_min - 1) {
      j_max_temp = config.nj - 1;
    } else {
      j_max_temp = NAVR - (config.nj_global - j_global_max);
    }
  }

  // positive pole
  if (doRingAverage_thisrank_positivepole == 1) {
    // b_j_flux
#pragma omp parallel for schedule(static)
    for (int i = 0; i < config.ni; i++) {
      int ic_idx = is + i;
      int N = config.nk_global;
      int ks_ring = ks;
      int ke_ring = ks + config.nk_global - 1;
      double *flux = (double *) malloc(N * sizeof(double));
      double *flux_diff = (double *) malloc(N * sizeof(double));
      double *base = (double *) malloc(N * sizeof(double));
      double *rem = (double *) malloc(N * sizeof(double));

      for (int j = 0; j <= j_max_temp; j++) {
        int jf_idx = js + j + 1;
        int j_chunk = j + j_global_min;
        for (int k = ks_ring; k <= ke_ring; k++) flux[k - ks_ring] = gem_ring[mag_bj_flux_ring][ic_idx][jf_idx][k];
        FilterBaseMode(flux, N, base);
        memcpy(rem, flux, N * sizeof(double));
        RingBlockAvg(rem, N, Nchunk[j_chunk]);
        for (int k = 0; k < N; k++) {
          flux_diff[k] = rem[k] - flux[k];
          flux[k] = rem[k] + base[k];
        }

        for (int blk = 0; blk < Nchunk[j_chunk]; blk++) {
          int n_cell = N / Nchunk[j_chunk];
          int start = blk * n_cell;
          double sum = 0.0;
          for (int m = 0; m < n_cell-1 ; m++) {
            sum += flux_diff[start + m];
            gem_ring[efield_ei_ring][ic_idx][jf_idx][start + ks + 1 + m] = -sum;
          }
        }

      }
      free(flux);
      free(flux_diff);
      free(base);
      free(rem);
    }

    // b_i_flux
#pragma omp parallel for schedule(static)
    for (int i = 0; i < config.ni + 1; i++) {
      int if_idx = is + i;
      int N = config.nk_global;
      int ks_ring = ks;
      int ke_ring = ks + config.nk_global - 1;
      double *flux = (double *) malloc(N * sizeof(double));
      double *flux_diff = (double *) malloc(N * sizeof(double));
      double *base = (double *) malloc(N * sizeof(double));
      double *rem = (double *) malloc(N * sizeof(double));
      for (int j = 0; j <= j_max_temp; j++) {
        int jc_idx = js + j;
        int j_chunk = j + j_global_min;

        for (int k = ks_ring; k <= ke_ring; k++) flux[k - ks_ring] = gem_ring[mag_bi_flux_ring][if_idx][jc_idx][k];
        FilterBaseMode(flux, N, base);
        memcpy(rem, flux, N * sizeof(double));
        RingBlockAvg(rem, N, Nchunk[j_chunk]);
        for (int k = 0; k < N; k++) {
          flux_diff[k] = rem[k] - flux[k];
          flux[k] = rem[k] + base[k];
        }
        for (int blk = 0; blk < Nchunk[j_chunk]; blk++) {
          int n_cell = N / Nchunk[j_chunk];
          int start = blk * n_cell;
          double sum = 0.0;
          for (int m = 0; m < n_cell-1 ; m++) {
            sum += flux_diff[start + m];
            gem_ring[efield_ej_ring][if_idx][jc_idx][start + ks + 1 + m] = sum;
          }
        }
      }
      free(flux);
      free(flux_diff);
      free(base);
      free(rem);
    }
  }

  // negative pole
  if (doRingAverage_thisrank_negativepole == 1) {
    // b_j_flux
#pragma omp parallel for schedule(static)
    for (int i = 0; i < config.ni; i++) {
      int ic_idx = is + i;
      int N = config.nk_global;
      int ks_ring = ks;
      int ke_ring = ks + config.nk_global - 1;
      double *flux = (double *) malloc(N * sizeof(double));
      double *flux_diff = (double *) malloc(N * sizeof(double));
      double *base = (double *) malloc(N * sizeof(double));
      double *rem = (double *) malloc(N * sizeof(double));

      for (int j = 0; j <= j_max_temp; j++) {
        int jf_rev = je - j;
        int j_chunk = j + (config.nj_global - j_global_max - 1);
        for (int k = ks_ring; k <= ke_ring; k++) flux[k - ks_ring] = gem_ring[mag_bj_flux_ring][ic_idx][jf_rev][k];
        FilterBaseMode(flux, N, base);
        memcpy(rem, flux, N * sizeof(double));
        RingBlockAvg(rem, N, Nchunk[j_chunk]);
        for (int k = 0; k < N; k++) {
          flux_diff[k] = rem[k] - flux[k];
          flux[k] = rem[k] + base[k];
        }
        for (int blk = 0; blk < Nchunk[j]; blk++) {
          int n_cell = N / Nchunk[j_chunk];
          int start = blk * n_cell;
          double sum = 0.0;
          for (int m = 0; m < n_cell-1 ; m++) {
            sum += flux_diff[start + m];
            gem_ring[efield_ei_ring][ic_idx][jf_rev][start + ks + 1 + m] = -sum;
          }
        }
      }
      free(flux);
      free(flux_diff);
      free(base);
      free(rem);
    }

    // b_i_flux
#pragma omp parallel for schedule(static)
    for (int i = 0; i < config.ni + 1; i++) {
      int if_idx = is + i;
      int N = config.nk_global;
      int ks_ring = ks;
      int ke_ring = ks + config.nk_global - 1;
      double *flux = (double *) malloc(N * sizeof(double));
      double *flux_diff = (double *) malloc(N * sizeof(double));
      double *base = (double *) malloc(N * sizeof(double));
      double *rem = (double *) malloc(N * sizeof(double));
      for (int j = 0; j <= j_max_temp; j++) {
        int jc_rev = je - j;
        int j_chunk = j + (config.nj_global - j_global_max - 1);
        for (int k = ks_ring; k <= ke_ring; k++) flux[k - ks_ring] = gem_ring[mag_bi_flux_ring][if_idx][jc_rev][k];
        FilterBaseMode(flux, N, base);
        memcpy(rem, flux, N * sizeof(double));
        RingBlockAvg(rem, N, Nchunk[j_chunk]);
        for (int k = 0; k < N; k++) {
          flux_diff[k] = rem[k] - flux[k];
          flux[k] = rem[k] + base[k];
        }
        for (int blk = 0; blk < Nchunk[j_chunk]; blk++) {
          int n_cell = N / Nchunk[j_chunk];
          int start = blk * n_cell;
          double sum = 0.0;
          for (int m = 0; m < n_cell-1 ; m++) {
            sum += flux_diff[start + m];
            gem_ring[efield_ej_ring][if_idx][jc_rev][start + ks + 1 + m] = sum;
          }
        }
      }
      free(flux);
      free(flux_diff);
      free(base);
      free(rem);
    }
  }

  int ks_ring = ks;
  int ke_ring = ks + config.nk_global - 1;
#pragma omp parallel for collapse(3) schedule(static)
  for (int i = is; i <= ie+1; i++) {
    for (int j = js; j <= je; j++) {
      for (int k = ks_ring; k <= ke_ring; k++) {
        double curlE = 1.0 * ((gem_ring[efield_ek_ring][i][j+1][k] - gem_ring[efield_ek_ring][i][j][k])
            - (gem_ring[efield_ej_ring][i][j][k+1] - gem_ring[efield_ej_ring][i][j][k]));
        gem_ring[mag_bi_ring][i][j][k] -= curlE / geo_ring[face_idir_ring][i][j][k];
      }
    }
  }

#pragma omp parallel for collapse(3) schedule(static)
  for (int i = is; i <= ie; i++) {
    for (int j = js; j <= je+1; j++) {
      for (int k = ks_ring; k <= ke_ring; k++) {
        double curlE = 1.0 * ((gem_ring[efield_ei_ring][i][j][k+1] - gem_ring[efield_ei_ring][i][j][k])
            - (gem_ring[efield_ek_ring][i+1][j][k] - gem_ring[efield_ek_ring][i][j][k]));
        gem_ring[mag_bj_ring][i][j][k] -= curlE / geo_ring[face_jdir_ring][i][j][k];
      }
    }
  }

#pragma omp parallel for collapse(3) schedule(static)
  for (int i = is; i <= ie; i++) {
    for (int j = js; j <= je; j++) {
      for (int k = ks_ring; k <= ke_ring+1; k++) {
        double curlE = 1.0 * ((gem_ring[efield_ej_ring][i+1][j][k] - gem_ring[efield_ej_ring][i][j][k])
            - (gem_ring[efield_ei_ring][i][j+1][k] - gem_ring[efield_ei_ring][i][j][k]));
        gem_ring[mag_bk_ring][i][j][k] -= curlE / geo_ring[face_kdir_ring][i][j][k];
      }
    }
  }

  if (fabs(x2min_global) < 1e-10 && proc_coords[1] == 0) {
#pragma omp parallel for collapse(2) schedule(static)
    for (int i = is; i <= ie; i++) {
      for (int k = ks_ring; k <= ke_ring; k++) {
        gem_ring[mag_bj_ring][i][js][k] = gem_ring[mag_bj_ring][i][js + 1][k];
      }
    }
  }

  if (fabs(x2max_global - PI) < 1e-10 && proc_coords[1] == config.proc_dims[1] - 1) {
#pragma omp parallel for collapse(2) schedule(static)
    for (int i = is; i <= ie; i++) {
      for (int k = ks_ring; k <= ke_ring; k++) {
        gem_ring[mag_bj_ring][i][je + 1][k] = gem_ring[mag_bj_ring][i][je][k];
      }
    }
  }

  // positive pole
  if (fabs(x2min_global) < 1e-10 && proc_coords[1] == 0) {
#pragma omp parallel for schedule(static)
    for (int i = is; i <= ie; i++) {
      double mean_loop_js = 0.0;
      for (int k = ks_ring; k <= ke_ring+1; k++) {
        mean_loop_js += gem_ring[mag_bk_ring][i][js][k];
      }
      mean_loop_js = mean_loop_js / (config.nk_global+1);
      for (int k = ks_ring; k <= ke_ring + 1; k++) {
        gem_ring[mag_bk_ring][i][js][k] -= mean_loop_js;
      }
    }
  }

  // negative pole
  if (fabs(x2max_global - PI) < 1e-10 && proc_coords[1] == config.proc_dims[1] - 1) {
#pragma omp parallel for schedule(static)
    for (int i = is; i <= ie; i++) {
      double mean_loop_je = 0.0;
      for (int k = ks_ring; k <= ke_ring+1; k++) {
        mean_loop_je += gem_ring[mag_bk_ring][i][je][k];
      }
      mean_loop_je = mean_loop_je / (config.nk_global+1);
      for (int k = ks_ring; k <= ke_ring + 1; k++) {
        gem_ring[mag_bk_ring][i][je][k] -= mean_loop_je;
      }
    }
  }
}

void HydroRingAverage(){
  if (proc_coords[2] == 0) {
    HydroRingAverage_base();
  }
}

void MagneticRingAverage() {
  if (proc_coords[2] == 0) {
    MagneticRingAverage_base();
  }
}

