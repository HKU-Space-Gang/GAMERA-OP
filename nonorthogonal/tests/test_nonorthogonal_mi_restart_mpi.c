#define _POSIX_C_SOURCE 200809L

#include "nonorthogonal_mi_restart.h"

#include <hdf5.h>
#include <mpi.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define REQUIRE(condition)                                                     \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "rank %d requirement failed at %s:%d: %s\n", rank,    \
              __FILE__, __LINE__, #condition);                                 \
      (void)MPI_Abort(MPI_COMM_WORLD, 1);                                      \
      return 1;                                                                \
    }                                                                          \
  } while (0)

static int create_empty_hdf5(const char *path) {
  hid_t file = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  return file < 0 || H5Fclose(file) < 0 ? -1 : 0;
}

static int add_unknown_link(const char *path) {
  hid_t file = H5Fopen(path, H5F_ACC_RDWR, H5P_DEFAULT);
  hid_t space = file < 0 ? -1 : H5Screate(H5S_SCALAR);
  hid_t dataset = space < 0 ? -1 : H5Dcreate2(
      file, "/mi/unknown", H5T_STD_U32LE, space, H5P_DEFAULT, H5P_DEFAULT,
      H5P_DEFAULT);
  uint32_t value = 1U;
  const int failed =
      file < 0 || space < 0 || dataset < 0 ||
      H5Dwrite(dataset, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT,
               &value) < 0;
  if (dataset >= 0) {
    H5Dclose(dataset);
  }
  if (space >= 0) {
    H5Sclose(space);
  }
  if (file >= 0) {
    H5Fclose(file);
  }
  return failed ? -1 : 0;
}

int main(int argc, char **argv) {
  int rank = -1;
  REQUIRE(MPI_Init(&argc, &argv) == MPI_SUCCESS);
  int size = 0;
  REQUIRE(MPI_Comm_rank(MPI_COMM_WORLD, &rank) == MPI_SUCCESS);
  REQUIRE(MPI_Comm_size(MPI_COMM_WORLD, &size) == MPI_SUCCESS);
  REQUIRE(size == 2);
  REQUIRE(argc == 2);

  char path[1024];
  REQUIRE(snprintf(path, sizeof(path), "%s/mi_restart_%ld.h5", argv[1],
                   (long)getpid()) > 0);
  gamera_mi_restart_state *state = NULL;
  if (rank == 0) {
    REQUIRE(gamera_mi_restart_create(3U, 4U, &state) == 0);
    state->solution_ready = 1U;
    state->snapshot_ready = 1U;
    state->next_update_s = 240.0;
    state->next_diagnostic_s = 240.0;
    state->update_count = 2U;
    state->snapshot_epoch_s = 120.0;
    state->snapshot_generation = 2U;
    state->maximum_cached_emf = 3.5e-4;
    state->maximum_prehold_difference = 8.0e-5;
    state->maximum_posthold_difference = 0.0;
    state->held_apply_count = 91U;
    state->precipitation_state_ready = 1U;
    state->dpb_initialized[0] = 1U;
    state->dpb_initialized[1] = 1U;
    state->dpb_radius_deg[0] = 20.5;
    state->dpb_radius_deg[1] = 21.25;
    for (size_t index = 0U; index < 24U; ++index) {
      state->potential_v[index] = -17000.0 + 137.0 * (double)index;
      state->pedersen_siemens[index] = 2.0 + 0.125 * (double)index;
      state->hall_siemens[index] = 1.0 + 0.0625 * (double)index;
    }
    REQUIRE(gamera_mi_restart_validate(state) == 0);
    REQUIRE(create_empty_hdf5(path) == 0);
    REQUIRE(gamera_mi_restart_append_hdf5(path, state) == 0);
    gamera_mi_restart_destroy(state);
    state = NULL;
    REQUIRE(gamera_mi_restart_read_hdf5(path, &state) == 0);
  }
  REQUIRE(gamera_mi_restart_broadcast(MPI_COMM_WORLD, 0, &state) == 0);
  REQUIRE(state != NULL && state->theta_points == 3U &&
          state->azimuth_points == 4U && state->update_count == 2U &&
          state->snapshot_generation == 2U &&
          state->held_apply_count == 91U && state->next_update_s == 240.0 &&
          state->precipitation_state_ready == 1U &&
          state->dpb_initialized[0] == 1U &&
          state->dpb_initialized[1] == 1U &&
          state->dpb_radius_deg[0] == 20.5 &&
          state->dpb_radius_deg[1] == 21.25 &&
          state->potential_v[0] == -17000.0 &&
          state->potential_v[23] == -17000.0 + 137.0 * 23.0 &&
          state->pedersen_siemens[23] == 2.0 + 0.125 * 23.0 &&
          state->hall_siemens[23] == 1.0 + 0.0625 * 23.0);
  unsigned char digest[GAMERA_MI_RESTART_DIGEST_BYTES];
  unsigned char root_digest[GAMERA_MI_RESTART_DIGEST_BYTES];
  REQUIRE(gamera_mi_restart_digest(state, digest) == 0);
  memcpy(root_digest, digest, sizeof(root_digest));
  REQUIRE(MPI_Bcast(root_digest, (int)sizeof(root_digest), MPI_BYTE, 0,
                    MPI_COMM_WORLD) == MPI_SUCCESS);
  REQUIRE(memcmp(root_digest, digest, sizeof(digest)) == 0);
  gamera_mi_restart_destroy(state);
  state = NULL;

  REQUIRE(MPI_Barrier(MPI_COMM_WORLD) == MPI_SUCCESS);
  if (rank == 0) {
    REQUIRE(add_unknown_link(path) == 0);
    REQUIRE(gamera_mi_restart_read_hdf5(path, &state) != 0);
    REQUIRE(state == NULL);
    REQUIRE(unlink(path) == 0);
  }
  REQUIRE(MPI_Barrier(MPI_COMM_WORLD) == MPI_SUCCESS);
  REQUIRE(MPI_Finalize() == MPI_SUCCESS);
  return 0;
}
