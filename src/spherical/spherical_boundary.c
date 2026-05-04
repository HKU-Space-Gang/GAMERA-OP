// This file contains boundary related functions for spherical coordinates, including the pole boundary.

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

// pole boundary condition only for spherical coordinates
void gas_bc_pole_j_low() {
  gather_value_RingAverage_arrays();
  if (proc_coords[2] == 0) {
    if (nbr_low[1] == MPI_PROC_NULL) {
      // pole boundary condition for j lowest boundary
#pragma omp parallel for collapse(5) schedule(static)
      for (int s = 0; s < NS1; s++) {
        for (int f = 0; f < NF_gas_prim; f++) {
          for (int i = 0; i < config.NI; i++) {
            for (int j = 0; j < NG; j++) {
              for (int k = 0; k < Ring_Nk; k++) {
                int diff = gas_v1_ring - gas_v1;
                if (f == gas_v2 || f == gas_v3) {
                  gas_ring[s][f][i][j][k] =
                      -gas_ring[s][f + diff][i][2 * NG - 1 - j + gas_onface_j[f]][NG
                          + (k - NG + config.nk_global / 2) % config.nk_global];
                } else {
                  gas_ring[s][f][i][j][k] =
                      gas_ring[s][f + diff][i][2 * NG - 1 - j + gas_onface_j[f]][NG
                          + (k - NG + config.nk_global / 2) % config.nk_global];
                }
              }
            }
          }
        }
      }
    }
  }
  broadcast_value_RingAverage_arrays();
}

void gas_bc_pole_j_high() {
  gather_value_RingAverage_arrays();
  if (proc_coords[2] == 0) {
    if (nbr_high[1] == MPI_PROC_NULL) {
      // pole boundary condition for j lowest boundary
#pragma omp parallel for collapse(5) schedule(static)
      for (int s = 0; s < NS1; s++) {
        for (int f = 0; f < NF_gas_prim; f++) {
          for (int i = 0; i < config.NI; i++) {
            for (int j = 0; j < NG; j++) {
              for (int k = 0; k < Ring_Nk; k++) {
                int diff = gas_v1_ring - gas_v1;
                if (f == gas_v2 || f == gas_v3) {
                  gas_ring[s][f][i][config.NJ - 1 - NG + j + gas_onface_j[f]][k] =
                      -gas_ring[s][f + diff][i][config.NJ - 1 - NG - 1 - j][NG
                          + (k - NG + config.nk_global / 2) % config.nk_global];
                } else {
                  gas_ring[s][f][i][config.NJ - 1 - NG + j + gas_onface_j[f]][k] =
                      gas_ring[s][f + diff][i][config.NJ - 1 - NG - 1 - j][NG
                          + (k - NG + config.nk_global / 2) % config.nk_global];
                }
              }
            }
          }
        }
      }
    }
  }
  broadcast_value_RingAverage_arrays();
}

void gem_bc_pole_j_low() {
  gather_value_RingAverage_arrays();
  if (proc_coords[2] == 0) {
    if (nbr_low[1] == MPI_PROC_NULL) {
      // pole boundary condition for j lowest boundary
#pragma omp parallel for collapse(4) schedule(static)
      for (int f = 0; f < NF_gem_prim; f++) {
        for (int i = 0; i < config.NI; i++) {
          for (int j = 0; j < NG; j++) {
            for (int k = 0; k < Ring_Nk; k++) {
              int diff = mag_bi_ring - mag_bi;
              if (f == mag_bi || f == mags_b1) {
                gem_ring[f][i][j][k] = gem_ring[f + diff][i][2 * NG - 1 - j + gem_onface_j[f]][NG
                    + (k - NG + config.nk_global / 2) % config.nk_global];
              } else {
                gem_ring[f][i][j][k] = -gem_ring[f + diff][i][2 * NG - 1 - j + gem_onface_j[f]][NG
                    + (k - NG + config.nk_global / 2) % config.nk_global];
              }
            }
          }
        }
      }
    }
  }
  broadcast_value_RingAverage_arrays();
}
void gem_bc_pole_j_high() {
  gather_value_RingAverage_arrays();
  if (proc_coords[2] == 0) {
    if (nbr_high[1] == MPI_PROC_NULL) {
      // pole boundary condition for j highest boundary
#pragma omp parallel for collapse(4) schedule(static)
      for (int f = 0; f < NF_gem_prim; f++) {
        for (int i = 0; i < config.NI; i++) {
          for (int j = 0; j < NG; j++) {
            for (int k = 0; k < Ring_Nk; k++) {
              int diff = mag_bi_ring - mag_bi;
              if (f == mag_bi || f == mags_b1) {
                gem_ring[f][i][config.NJ - 1 - NG + j + gem_onface_j[f]][k] =
                    gem_ring[f + diff][i][config.NJ - 1 - NG - 1 - j][NG
                        + (k - NG + config.nk_global / 2) % config.nk_global];
              } else {
                gem_ring[f][i][config.NJ - 1 - NG + j + gem_onface_j[f]][k] =
                    -gem_ring[f + diff][i][config.NJ - 1 - NG - 1 - j][NG
                        + (k - NG + config.nk_global / 2) % config.nk_global];
              }
            }
          }
        }
      }
    }
  }
  broadcast_value_RingAverage_arrays();
}
