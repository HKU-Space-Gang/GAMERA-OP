#ifndef SETUP_MPI_H
#define SETUP_MPI_H

#include <mpi.h>

// Rank and size of the MPI process
extern int rank, size;

// Cartesian communicator for grid topology
extern MPI_Comm comm_cart;

// Arrays to store neighboring ranks in each dimension (low and high)
extern int nbr_low[3], nbr_high[3];

// MPI datatypes for face communication in i, j, and k directions
extern MPI_Datatype i_face_type, j_face_type, k_face_type;

// MPI datatypes for sending and receiving data on each face in i, j, and k
// directions
extern MPI_Datatype i_face_send_low_type, i_face_recv_high_type, i_face_send_high_type,
    i_face_recv_low_type, j_face_send_low_type, j_face_recv_high_type,
    j_face_send_high_type, j_face_recv_low_type, k_face_send_low_type,
    k_face_recv_high_type, k_face_send_high_type, k_face_recv_low_type;

// Sets up the MPI communicator and datatypes for grid communication.
int setup_comm_and_types(int rank);

// Frees the MPI datatypes created for face communication.
void free_mpi_types();

#endif
