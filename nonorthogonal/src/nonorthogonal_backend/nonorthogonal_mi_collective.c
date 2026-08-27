#include "nonorthogonal_mi_collective.h"

#include <math.h>
#include <stddef.h>

int gamera_mi_collective_maximum_if_all_valid(
    MPI_Comm communicator, int local_failed, double local_maximum,
    double *global_maximum) {
  const int invalid = local_failed != 0 || global_maximum == NULL ||
                      !isfinite(local_maximum);
  const double local[2] = {invalid ? 1.0 : 0.0,
                           invalid ? 0.0 : local_maximum};
  double global[2] = {0.0, 0.0};
  if (MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_MAX, communicator) !=
      MPI_SUCCESS) {
    return -1;
  }
  if (global[0] != 0.0 || !isfinite(global[1])) {
    return -1;
  }
  *global_maximum = global[1];
  return 0;
}
