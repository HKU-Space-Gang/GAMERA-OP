#include "log.h"

#include "config.h"
#include "setup_mpi.h"
#include "utils.h"

static int boundary_exchange_mpi_indexed(double *array, int *onface_i,
                                         int *onface_j, int *onface_k) {
  int status = 0;
  int tag_i0 = 10, tag_i1 = 11, tag_j0 = 20, tag_j1 = 21, tag_k0 = 30,
      tag_k1 = 31;

  // i-dimension exchange
  {
    // High ghost <- Low active
    status = MPI_Sendrecv(array, 1, i_face_send_low_type, nbr_low[0], tag_i0,
                          array, 1, i_face_recv_high_type, nbr_high[0], tag_i0,
                          comm_cart, MPI_STATUS_IGNORE);
    if (status != MPI_SUCCESS) {
      log_error("Error: MPI_Sendrecv failed in dim i, send low and recv high");
      return -1;
    }

    // High active -> Low ghost
    status = MPI_Sendrecv(array, 1, i_face_send_high_type, nbr_high[0], tag_i1,
                          array, 1, i_face_recv_low_type, nbr_low[0], tag_i1,
                          comm_cart, MPI_STATUS_IGNORE);
    if (status != MPI_SUCCESS) {
      log_error("Error: MPI_Sendrecv failed in dim i, send high and recv low");
      return -1;
    }
  }

  // j-dimension exchange
  {
    // High ghost <- Low active
    status = MPI_Sendrecv(array, 1, j_face_send_low_type, nbr_low[1], tag_j0,
                          array, 1, j_face_recv_high_type, nbr_high[1], tag_j0,
                          comm_cart, MPI_STATUS_IGNORE);
    if (status != MPI_SUCCESS) {
      log_error("Error: MPI_Sendrecv failed in dim j, send low and recv high");
      return -1;
    }

    // High active -> Low ghost
    status = MPI_Sendrecv(array, 1, j_face_send_high_type, nbr_high[1], tag_j1,
                          array, 1, j_face_recv_low_type, nbr_low[1], tag_j1,
                          comm_cart, MPI_STATUS_IGNORE);
    if (status != MPI_SUCCESS) {
      log_error("Error: MPI_Sendrecv failed in dim j, send high and recv low");
      return -1;
    }
  }

  // k-dimension exchange
  {
    // High ghost <- Low active
    status = MPI_Sendrecv(array, 1, k_face_send_low_type, nbr_low[2], tag_k0,
                          array, 1, k_face_recv_high_type, nbr_high[2], tag_k0,
                          comm_cart, MPI_STATUS_IGNORE);
    if (status != MPI_SUCCESS) {
      log_error("Error: MPI_Sendrecv failed in dim k, send low and recv high");
      return -1;
    }

    // High active -> Low ghost
    status = MPI_Sendrecv(array, 1, k_face_send_high_type, nbr_high[2], tag_k1,
                          array, 1, k_face_recv_low_type, nbr_low[2], tag_k1,
                          comm_cart, MPI_STATUS_IGNORE);
    if (status != MPI_SUCCESS) {
      log_error("Error: MPI_Sendrecv failed in dim k, send high and recv low");
      return -1;
    }
  }
  return 0;
}

static int boundary_exchange_mpi_vector(double *array, int NI, int NJ, int NK) {
  // This only work for arrays at cell centers due to mpi type vector requires
  // stride to be a single constant
  int status = 0;
  int tag_i0 = 10, tag_i1 = 11, tag_j0 = 20, tag_j1 = 21, tag_k0 = 30,
      tag_k1 = 31;
  // i-dimension exchange
  {
    // Calculate source and destination offsets
    int src_offset = (NG * NJ * NK);
    int dst_offset = ((NI - NG - 1) * NJ * NK);

    // High ghost <- Low active
    status = MPI_Sendrecv(array + src_offset, 1, i_face_type, nbr_low[0],
                          tag_i0, array + dst_offset, 1, i_face_type,
                          nbr_high[0], tag_i0, comm_cart, MPI_STATUS_IGNORE);
    if (status != MPI_SUCCESS) {
      log_error("Error: MPI_Sendrecv failed in dim i, send low and recv high");
      return -1;
    }

    //  High active -> Low ghost
    src_offset = ((NI - 1 - 2 * NG) * NJ * NK);
    dst_offset = 0;
    status = MPI_Sendrecv(array + src_offset, 1, i_face_type, nbr_high[0],
                          tag_i1, array + dst_offset, 1, i_face_type,
                          nbr_low[0], tag_i1, comm_cart, MPI_STATUS_IGNORE);
    if (status != MPI_SUCCESS) {
      log_error("Error: MPI_Sendrecv failed in dim i, send high and recv low");
      return -1;
    }
  }

  // j-dimension exchange
  {
    // Calculate source and destination offsets
    int src_offset = (NG * NK);
    int dst_offset = ((NJ - NG - 1) * NK);

    // High ghost <- Low active
    status = MPI_Sendrecv(array + src_offset, 1, j_face_type, nbr_low[1],
                          tag_j0, array + dst_offset, 1, j_face_type,
                          nbr_high[1], tag_j0, comm_cart, MPI_STATUS_IGNORE);
    if (status != MPI_SUCCESS) {
      log_error("Error: MPI_Sendrecv failed in dim j, send low and recv high");
      return -1;
    }

    // High active -> Low ghost
    src_offset = ((NJ - 1 - 2 * NG) * NK);
    dst_offset = 0;
    status = MPI_Sendrecv(array + src_offset, 1, j_face_type, nbr_high[1],
                          tag_j1, array + dst_offset, 1, j_face_type,
                          nbr_low[1], tag_j1, comm_cart, MPI_STATUS_IGNORE);
    if (status != MPI_SUCCESS) {
      log_error("Error: MPI_Sendrecv failed in dim j, send high and recv low");
      return -1;
    }
  }

  // k-dimension exchange
  {
    // Calculate source and destination offsets
    int src_offset = NG;
    int dst_offset = NK - NG - 1;

    // High ghost <- Low active
    status = MPI_Sendrecv(array + src_offset, 1, k_face_type, nbr_low[2],
                          tag_k0, array + dst_offset, 1, k_face_type,
                          nbr_high[2], tag_k0, comm_cart, MPI_STATUS_IGNORE);
    if (status != MPI_SUCCESS) {
      log_error("Error: MPI_Sendrecv failed in dim k, send low and recv high");
      return -1;
    }

    //  High active -> Low ghost
    src_offset = NK - 1 - 2 * NG;
    dst_offset = 0;
    status = MPI_Sendrecv(array + src_offset, 1, k_face_type, nbr_high[2],
                          tag_k1, array + dst_offset, 1, k_face_type,
                          nbr_low[2], tag_k1, comm_cart, MPI_STATUS_IGNORE);
    if (status != MPI_SUCCESS) {
      log_error("Error: MPI_Sendrecv failed in dim k, send high and recv low");
      return -1;
    }
  }
  return 0;
}

// sync boundary values with neighbor processes
int boundary_exchange_5d(double *****field_arrays, int ns, int NI, int NJ,
                         int NK) {
  int status = 0;
  int tag_i0 = 10, tag_i1 = 11, tag_j0 = 20, tag_j1 = 21, tag_k0 = 30,
      tag_k1 = 31;
  // Exchange in all dimensions for each species
  for (int s = 0; s < ns; s++) {
    double *array = &field_arrays[s][0][0][0][0];
    status = boundary_exchange_mpi_vector(array, NI, NJ, NK);
    if (status != 0) {
      log_error("Error: boundary_exchange_mpi_vector failed");
      return -1;
    }
  }

  return 0;
}

// sync boundary values with neighbor processes
int boundary_exchange_4d(double ****field_arrays, int *onface_i, int *onface_j,
                         int *onface_k) {
  int status = 0;
  // Get pointer to contiguous data starting at nf_start
  double *array = &field_arrays[0][0][0][0];

  // Exchange boundaries using MPI indexed types
  status = boundary_exchange_mpi_indexed(array, onface_i, onface_j, onface_k);
  if (status != 0) {
    log_error("Error: boundary_exchange_mpi_indexed failed");
    return -1;
  }

  return 0;
}
