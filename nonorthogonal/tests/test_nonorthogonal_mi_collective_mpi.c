#include "nonorthogonal_mi_collective.h"

#include <mpi.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define REQUIRE(condition)                                                     \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "rank %d requirement failed at %s:%d: %s\n", rank,    \
              __FILE__, __LINE__, #condition);                                 \
      (void)MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);                           \
      return EXIT_FAILURE;                                                     \
    }                                                                          \
  } while (0)

int main(int argc, char **argv) {
  int rank = -1;
  int size = 0;
  REQUIRE(MPI_Init(&argc, &argv) == MPI_SUCCESS);
  REQUIRE(MPI_Comm_rank(MPI_COMM_WORLD, &rank) == MPI_SUCCESS);
  REQUIRE(MPI_Comm_size(MPI_COMM_WORLD, &size) == MPI_SUCCESS);
  REQUIRE(size == 2);

  double global_maximum = -123.0;
  const int local_failed = rank == 1;
  REQUIRE(gamera_mi_collective_maximum_if_all_valid(
              MPI_COMM_WORLD, local_failed, 10.0 + (double)rank,
              &global_maximum) != 0);
  REQUIRE(global_maximum == -123.0);

  /* A collective immediately after the asymmetric failure proves that both
   * ranks completed the same failure reduction rather than stranding rank 0. */
  int survivors = 0;
  const int local_survivor = 1;
  REQUIRE(MPI_Allreduce(&local_survivor, &survivors, 1, MPI_INT, MPI_SUM,
                        MPI_COMM_WORLD) == MPI_SUCCESS);
  REQUIRE(survivors == size);

  global_maximum = -456.0;
  REQUIRE(gamera_mi_collective_maximum_if_all_valid(
              MPI_COMM_WORLD, 0, 20.0 + (double)rank,
              rank == 1 ? NULL : &global_maximum) != 0);
  if (rank == 0) {
    REQUIRE(global_maximum == -456.0);
  }

  global_maximum = -1.0;
  REQUIRE(gamera_mi_collective_maximum_if_all_valid(
              MPI_COMM_WORLD, 0, rank == 0 ? 2.5 : 7.25,
              &global_maximum) == 0);
  REQUIRE(fabs(global_maximum - 7.25) <= 1.0e-15);

  REQUIRE(MPI_Finalize() == MPI_SUCCESS);
  return EXIT_SUCCESS;
}
