#include "setup_mpi.h"

#include <stdlib.h>

#include "log.h"

#include "config.h"

// Rank and size of the MPI process
int rank, size;

// Cartesian communicator for grid topology
MPI_Comm comm_cart;

// Arrays to store neighboring ranks in each dimension (low and high)
int nbr_low[3], nbr_high[3];

// MPI datatypes for face communication in i, j, and k directions
MPI_Datatype i_face_type, j_face_type, k_face_type;

// MPI datatypes for sending and receiving data on each face in i, j, and k
// directions
MPI_Datatype i_face_send_low_type, i_face_recv_high_type, i_face_send_high_type,
    i_face_recv_low_type, j_face_send_low_type, j_face_recv_high_type,
    j_face_send_high_type, j_face_recv_low_type, k_face_send_low_type,
    k_face_recv_high_type, k_face_send_high_type, k_face_recv_low_type;

static int setup_cart_comm(int rank) {
  // periodic boundaries
  // 0: The dimension is not periodic
  // 1: The dimension is periodic
  // int periods[3] = {0, 0, 0};
  int periods[3] = {(problem_config.Boundary_i == BC_PERIODIC) ? 1 : 0,
                    (problem_config.Boundary_j == BC_PERIODIC) ? 1 : 0,
                    (problem_config.Boundary_k == BC_PERIODIC) ? 1 : 0};
  int reorder = 0;
  int status = MPI_Cart_create(MPI_COMM_WORLD, 3, config.proc_dims, periods,
                               reorder, &comm_cart);
  if (status != 0) {
    log_error("Failed to create cartisian communicator.");
    return status;
  }
  // get the coordinates of the current process
  status = MPI_Cart_coords(comm_cart, rank, 3, proc_coords);
  if (status != 0) {
    log_error("Failed to get the coordinates of the current process.");
    return status;
  }

  // get the rank of the neighbor processes
  for (int dim = 0; dim < 3; dim++) {
    status = MPI_Cart_shift(comm_cart, dim, 1, &nbr_low[dim], &nbr_high[dim]);
    if (status != 0) {
      log_error("Failed to get neighbor ranks.");
      return status;
    }
  }
  return 0;
}

// Create MPI datatypes for boundary exchange, with 4d array continous in memory
static void create_mpi_types_indexed(int nf, int NI, int NJ, int NK, int ng,
                                     int *onface_i, int *onface_j,
                                     int *onface_k) {
  // Allocate maximum needed size
  int max_blocks =
      nf * NI * NJ;  // Maximum number of blocks needed (for k-faces)
  int *blocklens = (int *)malloc(max_blocks * sizeof(int));
  int *displs = (int *)malloc(max_blocks * sizeof(int));

  if (blocklens == NULL || displs == NULL) {
    log_error("Failed to allocate memory for MPI type arrays");
    return;
  }

  // For i-faces: contiguous in j,k dimensions
  for (int f = 0; f < nf; f++) {
    blocklens[f] = ng * NJ * NK;  // size of each block
    // Send from low active region (NG + onface)
    displs[f] = (NI * NJ * NK) * f + ((NG + onface_i[f]) * NJ * NK);
  }
  MPI_Type_indexed(nf, blocklens, displs, MPI_DOUBLE, &i_face_send_low_type);
  MPI_Type_commit(&i_face_send_low_type);

  // Receive into high ghost region
  for (int f = 0; f < nf; f++) {
    blocklens[f] = ng * NJ * NK;
    displs[f] = (NI * NJ * NK) * f + ((NI - NG - 1 + onface_i[f]) * NJ * NK);
  }
  MPI_Type_indexed(nf, blocklens, displs, MPI_DOUBLE, &i_face_recv_high_type);
  MPI_Type_commit(&i_face_recv_high_type);

  // Send from high active region
  for (int f = 0; f < nf; f++) {
    blocklens[f] = ng * NJ * NK;
    displs[f] = (NI * NJ * NK) * f + ((NI - 2 * NG - 1) * NJ * NK);
  }
  MPI_Type_indexed(nf, blocklens, displs, MPI_DOUBLE, &i_face_send_high_type);
  MPI_Type_commit(&i_face_send_high_type);

  // Receive into low ghost region
  for (int f = 0; f < nf; f++) {
    blocklens[f] = ng * NJ * NK;
    displs[f] = (NI * NJ * NK) * f;
  }
  MPI_Type_indexed(nf, blocklens, displs, MPI_DOUBLE, &i_face_recv_low_type);
  MPI_Type_commit(&i_face_recv_low_type);

  // For j-faces: contiguous in k dimension
  for (int f = 0; f < nf; f++) {
    for (int i = 0; i < NI; i++) {
      int idx = f * NI + i;
      blocklens[idx] = ng * NK;  // size of each block
      // Send from low active region
      displs[idx] =
          (NI * NJ * NK) * f + (i * NJ * NK) + ((NG + onface_j[f]) * NK);
    }
  }
  MPI_Type_indexed(nf * NI, blocklens, displs, MPI_DOUBLE,
                   &j_face_send_low_type);
  MPI_Type_commit(&j_face_send_low_type);

  // Receive into high ghost region
  for (int f = 0; f < nf; f++) {
    for (int i = 0; i < NI; i++) {
      int idx = f * NI + i;
      blocklens[idx] = ng * NK;
      displs[idx] = (NI * NJ * NK) * f + (i * NJ * NK) +
                    ((NJ - NG - 1 + onface_j[f]) * NK);
    }
  }
  MPI_Type_indexed(nf * NI, blocklens, displs, MPI_DOUBLE,
                   &j_face_recv_high_type);
  MPI_Type_commit(&j_face_recv_high_type);

  // Send from high active region
  for (int f = 0; f < nf; f++) {
    for (int i = 0; i < NI; i++) {
      int idx = f * NI + i;
      blocklens[idx] = ng * NK;
      displs[idx] =
          (NI * NJ * NK) * f + (i * NJ * NK) + ((NJ - 2 * NG - 1) * NK);
    }
  }
  MPI_Type_indexed(nf * NI, blocklens, displs, MPI_DOUBLE,
                   &j_face_send_high_type);
  MPI_Type_commit(&j_face_send_high_type);

  // Receive into low ghost region
  for (int f = 0; f < nf; f++) {
    for (int i = 0; i < NI; i++) {
      int idx = f * NI + i;
      blocklens[idx] = ng * NK;
      displs[idx] = (NI * NJ * NK) * f + (i * NJ * NK);
    }
  }
  MPI_Type_indexed(nf * NI, blocklens, displs, MPI_DOUBLE,
                   &j_face_recv_low_type);
  MPI_Type_commit(&j_face_recv_low_type);

  // For k-faces: contiguous blocks
  for (int f = 0; f < nf; f++) {
    for (int i = 0; i < NI; i++) {
      for (int j = 0; j < NJ; j++) {
        int idx = (f * NI + i) * NJ + j;
        blocklens[idx] = ng;  // size of each block
        // Send from low active region
        displs[idx] =
            (NI * NJ * NK) * f + (i * NJ * NK) + (j * NK) + (NG + onface_k[f]);
      }
    }
  }
  MPI_Type_indexed(nf * NI * NJ, blocklens, displs, MPI_DOUBLE,
                   &k_face_send_low_type);
  MPI_Type_commit(&k_face_send_low_type);

  // Receive into high ghost region
  for (int f = 0; f < nf; f++) {
    for (int i = 0; i < NI; i++) {
      for (int j = 0; j < NJ; j++) {
        int idx = (f * NI + i) * NJ + j;
        blocklens[idx] = ng;
        displs[idx] = (NI * NJ * NK) * f + (i * NJ * NK) + (j * NK) +
                      (NK - NG - 1 + onface_k[f]);
      }
    }
  }
  MPI_Type_indexed(nf * NI * NJ, blocklens, displs, MPI_DOUBLE,
                   &k_face_recv_high_type);
  MPI_Type_commit(&k_face_recv_high_type);

  // Send from high active region
  for (int f = 0; f < nf; f++) {
    for (int i = 0; i < NI; i++) {
      for (int j = 0; j < NJ; j++) {
        int idx = (f * NI + i) * NJ + j;
        blocklens[idx] = ng;
        displs[idx] =
            (NI * NJ * NK) * f + (i * NJ * NK) + (j * NK) + (NK - 2 * NG - 1);
      }
    }
  }
  MPI_Type_indexed(nf * NI * NJ, blocklens, displs, MPI_DOUBLE,
                   &k_face_send_high_type);
  MPI_Type_commit(&k_face_send_high_type);

  // Receive into low ghost region
  for (int f = 0; f < nf; f++) {
    for (int i = 0; i < NI; i++) {
      for (int j = 0; j < NJ; j++) {
        int idx = (f * NI + i) * NJ + j;
        blocklens[idx] = ng;
        displs[idx] = (NI * NJ * NK) * f + (i * NJ * NK) + (j * NK);
      }
    }
  }
  MPI_Type_indexed(nf * NI * NJ, blocklens, displs, MPI_DOUBLE,
                   &k_face_recv_low_type);
  MPI_Type_commit(&k_face_recv_low_type);

  free(blocklens);
  free(displs);
}

static void free_mpi_types_indexed() {
  // Free MPI datatypes
  MPI_Type_free(&i_face_send_low_type);
  MPI_Type_free(&i_face_recv_high_type);
  MPI_Type_free(&i_face_send_high_type);
  MPI_Type_free(&i_face_recv_low_type);
  MPI_Type_free(&j_face_send_low_type);
  MPI_Type_free(&j_face_recv_high_type);
  MPI_Type_free(&j_face_send_high_type);
  MPI_Type_free(&j_face_recv_low_type);
  MPI_Type_free(&k_face_send_low_type);
  MPI_Type_free(&k_face_recv_high_type);
  MPI_Type_free(&k_face_send_high_type);
  MPI_Type_free(&k_face_recv_low_type);
}

static void create_mpi_types_vector(int nf, int NI, int NJ, int NK, int ng) {
  // For i-faces: contiguous in j,k dimensions
  MPI_Type_vector(nf, ng * NJ * NK, NI * NJ * NK, MPI_DOUBLE, &i_face_type);
  MPI_Type_commit(&i_face_type);

  // For j-faces: contiguous in k dimension for each nf step
  MPI_Type_vector(nf * NI, ng * NK, NJ * NK, MPI_DOUBLE, &j_face_type);
  MPI_Type_commit(&j_face_type);

  // For k-faces: contiguous blocks for each nf step
  MPI_Type_vector(nf * NI * NJ, ng, NK, MPI_DOUBLE, &k_face_type);
  MPI_Type_commit(&k_face_type);
}

// Free MPI datatypes
static void free_mpi_types_vector() {
  MPI_Type_free(&i_face_type);
  MPI_Type_free(&j_face_type);
  MPI_Type_free(&k_face_type);
}

int setup_comm_and_types(int rank) {
  int status = setup_cart_comm(rank);
  if (status != 0) {
    return status;
  }
  create_mpi_types_vector(NF_gas_prim, config.NI, config.NJ, config.NK, NG);

  create_mpi_types_indexed(NF_gem_prim, config.NI, config.NJ, config.NK, NG,
                           gem_onface_i, gem_onface_j, gem_onface_k);

  return 0;
}

void free_mpi_types() {
  free_mpi_types_vector();
  free_mpi_types_indexed();
  // Free the Cartesian communicator
  MPI_Comm_free(&comm_cart);
}
