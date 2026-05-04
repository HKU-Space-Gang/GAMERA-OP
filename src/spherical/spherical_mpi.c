// This file contains MPI related functions for spherical geometry, including:
// gather and broadcast functions for ring average, and other MPI communication related to spherical geometry.

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


void gather_Etackle_RingAverage_arrays() {
  if (config.proc_dims[2] == 1) {
    // Etackle gem_ring array
#pragma omp parallel for collapse(4) schedule(static)
    for (int m = efield_tackle_ei_ring; m <= efield_tackle_ek_ring; m++) {
      for (int i = isg; i <= ieg + 1; i++) {
        for (int j = jsg; j <= jeg + 1; j++) {
          for (int k = ksg; k <= keg + 1; k++) {
            int diff = efield_ei - efield_tackle_ei_ring;
            gem_ring[m][i][j][k] = gem[m + diff][i][j][k];
          }
        }
      }
    }
  } else {
    // Etackle gem_ring array
    for (int m = efield_tackle_ei_ring; m <= efield_tackle_ek_ring; m++) {
      int diff = efield_ei - efield_tackle_ei_ring;
      RingAverage_gather_data_3d_array(gem[m + diff], gem_ring[m]);
    }
  }
}

void broadcast_Etackle_RingAverage_arrays() {
  if (config.proc_dims[2] == 1) {
    // Etackle gem_ring array
#pragma omp parallel for collapse(4) schedule(static)
    for (int m = efield_tackle_ei_ring; m <= efield_tackle_ek_ring; m++) {
      for (int i = isg; i <= ieg + 1; i++) {
        for (int j = jsg; j <= jeg + 1; j++) {
          for (int k = ksg; k <= keg + 1; k++) {
            int diff = efield_ei - efield_tackle_ei_ring;
            gem[m + diff][i][j][k] = gem_ring[m][i][j][k];
          }
        }
      }
    }
  } else {
    // Etackle gem_ring array
    for (int m = efield_tackle_ei_ring; m <= efield_tackle_ek_ring; m++) {
      int diff = efield_ei - efield_tackle_ei_ring;
      RingAverage_broadcast_data_3d_array(gem_ring[m], gem[m + diff]);
    }
  }
}


// Ring relaxed functions and variables
int RingAverage_gather_data_3d_array(double ***field_array_scr, double ***field_array_dst) {
  // This function is for ring average
  // The rank with proc_coords[2] == 0 gathers data from the corresponding ranks
  // field_array_scr's size is config.NI * config.NJ * config.NK
  // field_array_dst's size is config.NI * config.NJ * RingAverage_Nk
  int k_start_local, k_end_local, k_start_global;
  int local_size = config.NI * config.NJ * config.NK;
  // allocate local array to send data
  double *local_array = malloc(local_size * sizeof(double));
  if (local_array == NULL) {
    fprintf(stderr, "Memory allocation failed for local_array\n");
    return -1;
  }
  for (int i = 0; i < config.NI; i++) {
    for (int j = 0; j < config.NJ; j++) {
      for (int k = 0; k < config.NK; k++) {
        int idx = i * config.NJ * config.NK + j * config.NK + k;
        local_array[idx] = field_array_scr[i][j][k];
      }
    }
  }
  // Rank with proc_coords[2] == 0 gathers data from the corresponding ranks
  if (proc_coords[2] == 0) {
    // when rank = 0
    k_end_local = (0 == config.proc_dims[2] - 1) ? config.NK-1 : ke;
    for (int i = 0; i < config.NI; i++) {
      for (int j = 0; j < config.NJ; j++) {
        for (int k = 0; k <= k_end_local; k++) {
          int src_idx = i * config.NJ * config.NK + j * config.NK + k;
          field_array_dst[i][j][k] = local_array[src_idx];
        }
      }
    }
    // Then receive and place data from other processes in k-direction
    for (int proc_k = 1; proc_k < config.proc_dims[2]; proc_k++) {
      double *recv_buffer = malloc(local_size * sizeof(double));
      if (recv_buffer == NULL) {
        fprintf(stderr, "Memory allocation failed for recv_buffer\n");
        free(local_array);
        return -1;
      }

      // Calculate source process rank
      int src_coords[3] = {proc_coords[0], proc_coords[1], proc_k};
      int src_rank;
      MPI_Cart_rank(comm_cart, src_coords, &src_rank);

      MPI_Recv(recv_buffer, local_size, MPI_DOUBLE,
               src_rank, 303, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      k_start_local = (proc_k == 0) ? 0 : ks;
      k_end_local = (proc_k == config.proc_dims[2] - 1) ? config.NK-1 : ke;
      k_start_global = (proc_k == 0) ? 0 : ks + proc_k * config.nk;
      for (int i = 0; i < config.NI; i++) {
        for (int j = 0; j < config.NJ; j++) {
          for (int k = k_start_local; k <= k_end_local; k++) {
            int src_idx = i * config.NJ * config.NK + j * config.NK + k;
            field_array_dst[i][j][k_start_global+k-k_start_local] = recv_buffer[src_idx];
          }
        }
      }
      free(recv_buffer);
    }
  } else {
    // For non-root processes, send data to root process (proc_coords[2] == 0)

    // Calculate destination process rank (root in k-direction)
    int dst_coords[3] = {proc_coords[0], proc_coords[1], 0};
    int dst_rank;
    MPI_Cart_rank(comm_cart, dst_coords, &dst_rank);

    // Send data to root process
    MPI_Send(local_array, local_size, MPI_DOUBLE, dst_rank, 303, MPI_COMM_WORLD);
  }

  free(local_array);
  return 0;
}

int RingAverage_broadcast_data_3d_array(double ***field_array_scr, double ***field_array_dst) {
  // This function is for ring average broadcast
  // The rank with proc_coords[2] == 0 broadcasts data to the corresponding ranks
  // field_array_scr's size is config.NI * config.NJ * RingAverage_Nk (only meaningful for root)
  // field_array_dst's size is config.NI * config.NJ * config.NK
  int local_size = config.NI * config.NJ * config.NK;

  // Rank with proc_coords[2] == 0 broadcasts data to the corresponding ranks
  if (proc_coords[2] == 0) {
    // First, prepare and send data to other processes in k-direction
    for (int proc_k = 0; proc_k < config.proc_dims[2]; proc_k++) {
      // Calculate the k range for this target process
      int k_start, k_end;
      k_start = proc_k * config.nk;
      k_end = proc_k * config.nk + config.NK - 1;

      double *scr_array = malloc(local_size * sizeof(double));
      for (int i = 0; i < config.NI; i++) {
        for (int j = 0; j < config.NJ; j++) {
          for (int k = k_start; k <= k_end; k++) {
            int src_idx = i * config.NJ * config.NK + (j - 0) * config.NK + (k - k_start);
            scr_array[src_idx] = field_array_scr[i][j][k];
          }
        }
      }

      // Calculate dest process rank
      int dest_coords[3] = {proc_coords[0], proc_coords[1], proc_k};
      int dest_rank;
      MPI_Cart_rank(comm_cart, dest_coords, &dest_rank);

      if (proc_k == 0) {
        // copy data to local array
        for (int i = 0; i < config.NI; i++) {
          for (int j = 0; j < config.NJ; j++) {
            for (int k = 0; k < config.NK; k++) {
              field_array_dst[i][j][k] = scr_array[i * config.NJ * config.NK + j * config.NK + k];
            }
          }
        }
      } else {
        MPI_Send(scr_array, local_size, MPI_DOUBLE, dest_rank, 304, MPI_COMM_WORLD);
      }
      free(scr_array);
    }

  } else {
    // For non-root processes, receive data from root process

    // Calculate source process rank (root in k-direction)
    int src_coords[3] = {proc_coords[0], proc_coords[1], 0};
    int src_rank;
    MPI_Cart_rank(comm_cart, src_coords, &src_rank);

    // receive data only when proc_coords[0] == 0
    double *dest_array = malloc(local_size * sizeof(double));
    MPI_Recv(dest_array, local_size, MPI_DOUBLE, src_rank, 304, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i = 0; i < config.NI; i++) {
      for (int j = 0; j < config.NJ; j++) {
        for (int k = 0; k < config.NK; k++) {
          field_array_dst[i][j][k] = dest_array[i * config.NJ * config.NK + j * config.NK + k];
        }
      }
    }
    free(dest_array);
  }

  return 0;
}

void gather_value_RingAverage_arrays() {
  if (config.proc_dims[2] == 1) {
    // gas_ring array
#pragma omp parallel for collapse(5) schedule(static)
    for (int s = 0; s < NS; s++) {
      for (int m = gas_rho_ring; m <= gas_p_S_ring; m++) {
        for (int i = isg; i <= ieg + 1; i++) {
          for (int j = jsg; j <= jeg + 1; j++) {
            for (int k = ksg; k <= keg + 1; k++) {
              int diff = gas_rho - gas_rho_ring;
              gas_ring[s][m][i][j][k] = gas[s][m + diff][i][j][k];
            }
          }
        }
      }
    }

    // gem_ring array
#pragma omp parallel for collapse(4) schedule(static)
    for (int m = mag_bi_ring; m <= mags_b3_ring; m++) {
      for (int i = isg; i <= ieg + 1; i++) {
        for (int j = jsg; j <= jeg + 1; j++) {
          for (int k = ksg; k <= keg + 1; k++) {
            int diff = mag_bi - mag_bi_ring;
            gem_ring[m][i][j][k] = gem[m + diff][i][j][k];
          }
        }
      }
    }
  } else {
    // gas_ring array
    for (int s = 0; s < NS; s++) {
      for (int m = gas_rho_ring; m <= gas_p_S_ring; m++) {
        int diff = gas_rho - gas_rho_ring;
        RingAverage_gather_data_3d_array(gas[s][m + diff], gas_ring[s][m]);
      }
    }

    // gem_ring array
    for (int m = mag_bi_ring; m <= mags_b3_ring; m++) {
      int diff = mag_bi - mag_bi_ring;
      RingAverage_gather_data_3d_array(gem[m + diff], gem_ring[m]);
    }
  }

}

void broadcast_value_RingAverage_arrays() {
  if (config.proc_dims[2] == 1) {
    // gas_ring array
#pragma omp parallel for collapse(5) schedule(static)
    for (int s = 0; s < NS; s++) {
      for (int m = gas_rho_ring; m <= gas_p_S_ring; m++) {
        for (int i = isg; i <= ieg + 1; i++) {
          for (int j = jsg; j <= jeg + 1; j++) {
            for (int k = ksg; k <= keg + 1; k++) {
              int diff = gas_rho - gas_rho_ring;
              gas[s][m + diff][i][j][k] = gas_ring[s][m][i][j][k];
            }
          }
        }
      }
    }

    // gem_ring array
#pragma omp parallel for collapse(4) schedule(static)
    for (int m = mag_bi_ring; m <= mags_b3_ring; m++) {
      for (int i = isg; i <= ieg + 1; i++) {
        for (int j = jsg; j <= jeg + 1; j++) {
          for (int k = ksg; k <= keg + 1; k++) {
            int diff = mag_bi - mag_bi_ring;
            gem[m + diff][i][j][k] = gem_ring[m][i][j][k];
          }
        }
      }
    }
  } else {
    // gas_ring array
    for (int s = 0; s < NS; s++) {
      for (int m = gas_rho_ring; m <= gas_p_S_ring; m++) {
        int diff = gas_rho - gas_rho_ring;
        RingAverage_broadcast_data_3d_array(gas_ring[s][m], gas[s][m + diff]);
      }
    }

    // gem_ring array
    for (int m = mag_bi_ring; m <= mags_b3_ring; m++) {
      int diff = mag_bi - mag_bi_ring;
      RingAverage_broadcast_data_3d_array(gem_ring[m], gem[m + diff]);
    }
  }
}
