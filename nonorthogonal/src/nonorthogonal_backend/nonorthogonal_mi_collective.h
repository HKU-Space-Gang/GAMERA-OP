#ifndef GAMERA_NONORTHOGONAL_MI_COLLECTIVE_H
#define GAMERA_NONORTHOGONAL_MI_COLLECTIVE_H

#include <mpi.h>

/*
 * Complete one common collective before reporting a rank-local failure.
 * On success, global_maximum receives max(local_maximum) over communicator.
 * On failure, every rank returns -1 and global_maximum is left unchanged.
 */
int gamera_mi_collective_maximum_if_all_valid(
    MPI_Comm communicator, int local_failed, double local_maximum,
    double *global_maximum);

#endif
