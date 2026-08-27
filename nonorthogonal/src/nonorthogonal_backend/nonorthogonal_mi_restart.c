#define _POSIX_C_SOURCE 200809L

#include "nonorthogonal_mi_restart.h"

#include "nonorthogonal_sha256.h"

#include <hdf5.h>

#include <float.h>
#include <fcntl.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

enum {
  MI_HEMISPHERE_COUNT = 2,
  MI_RESTART_MAX_THETA_POINTS = 16384,
  MI_RESTART_MAX_AZIMUTH_POINTS = 65536,
  MI_RESTART_MAX_POLAR_POINTS = 268435456
};

static int checked_product(size_t first, size_t second, size_t *product) {
  if (product == NULL || (first != 0U && second > SIZE_MAX / first)) {
    return 0;
  }
  *product = first * second;
  return 1;
}

static int potential_count(size_t theta_points, size_t azimuth_points,
                           size_t *count) {
  size_t polar = 0U;
  return theta_points >= 2U && azimuth_points >= 2U &&
         theta_points <= MI_RESTART_MAX_THETA_POINTS &&
         azimuth_points <= MI_RESTART_MAX_AZIMUTH_POINTS &&
         checked_product(theta_points, azimuth_points, &polar) &&
         polar <= MI_RESTART_MAX_POLAR_POINTS &&
         checked_product(MI_HEMISPHERE_COUNT, polar, count);
}

int gamera_mi_restart_create(size_t theta_points, size_t azimuth_points,
                             gamera_mi_restart_state **state) {
  if (state == NULL) {
    return -1;
  }
  *state = NULL;
  size_t count = 0U;
  if (!potential_count(theta_points, azimuth_points, &count) ||
      count > SIZE_MAX / sizeof(double)) {
    return -1;
  }
  gamera_mi_restart_state *candidate =
      (gamera_mi_restart_state *)calloc(1U, sizeof(*candidate));
  if (candidate == NULL) {
    return -1;
  }
  candidate->potential_v = (double *)calloc(count, sizeof(double));
  candidate->pedersen_siemens = (double *)calloc(count, sizeof(double));
  candidate->hall_siemens = (double *)calloc(count, sizeof(double));
  if (candidate->potential_v == NULL || candidate->pedersen_siemens == NULL ||
      candidate->hall_siemens == NULL) {
    free(candidate->hall_siemens);
    free(candidate->pedersen_siemens);
    free(candidate->potential_v);
    free(candidate);
    return -1;
  }
  candidate->schema_version = GAMERA_MI_RESTART_SCHEMA;
  candidate->theta_points = theta_points;
  candidate->azimuth_points = azimuth_points;
  *state = candidate;
  return 0;
}

void gamera_mi_restart_destroy(gamera_mi_restart_state *state) {
  if (state != NULL) {
    free(state->potential_v);
    free(state->pedersen_siemens);
    free(state->hall_siemens);
    memset(state, 0, sizeof(*state));
    free(state);
  }
}

int gamera_mi_restart_validate(const gamera_mi_restart_state *state) {
  size_t count = 0U;
  if (state == NULL || state->schema_version != GAMERA_MI_RESTART_SCHEMA ||
      state->solution_ready != 1U || state->snapshot_ready != 1U ||
      !potential_count(state->theta_points, state->azimuth_points, &count) ||
      state->potential_v == NULL || state->pedersen_siemens == NULL ||
      state->hall_siemens == NULL || state->update_count == 0U ||
      state->snapshot_generation != state->update_count ||
      !isfinite(state->next_update_s) || !(state->next_update_s > 0.0) ||
      !isfinite(state->next_diagnostic_s) ||
      state->next_diagnostic_s < 0.0 ||
      !isfinite(state->snapshot_epoch_s) || state->snapshot_epoch_s < 0.0 ||
      !isfinite(state->maximum_cached_emf) ||
      state->maximum_cached_emf < 0.0 ||
      !isfinite(state->maximum_prehold_difference) ||
      state->maximum_prehold_difference < 0.0 ||
      !isfinite(state->maximum_posthold_difference) ||
      state->maximum_posthold_difference < 0.0 ||
      state->precipitation_state_ready > 1U ||
      state->dpb_initialized[0] > 1U || state->dpb_initialized[1] > 1U ||
      !isfinite(state->dpb_radius_deg[0]) ||
      !isfinite(state->dpb_radius_deg[1]) ||
      state->dpb_radius_deg[0] < 0.0 || state->dpb_radius_deg[1] < 0.0 ||
      (state->precipitation_state_ready == 0U &&
       (state->dpb_initialized[0] != 0U ||
        state->dpb_initialized[1] != 0U || state->dpb_radius_deg[0] != 0.0 ||
        state->dpb_radius_deg[1] != 0.0))) {
    return -1;
  }
  for (size_t index = 0U; index < count; ++index) {
    if (!isfinite(state->potential_v[index]) ||
        !isfinite(state->pedersen_siemens[index]) ||
        !isfinite(state->hall_siemens[index]) ||
        state->pedersen_siemens[index] < 0.0 ||
        state->hall_siemens[index] < 0.0 ||
        (state->precipitation_state_ready == 0U &&
         (state->pedersen_siemens[index] != 0.0 ||
          state->hall_siemens[index] != 0.0))) {
      return -1;
    }
  }
  return 0;
}

static void digest_u32(gamera_no_sha256 *sha, uint32_t value) {
  unsigned char encoded[4];
  for (size_t byte = 0U; byte < sizeof(encoded); ++byte) {
    encoded[sizeof(encoded) - 1U - byte] = (unsigned char)(value & 0xffU);
    value >>= 8U;
  }
  gamera_no_sha256_update(sha, encoded, sizeof(encoded));
}

static void digest_u64(gamera_no_sha256 *sha, uint64_t value) {
  unsigned char encoded[8];
  for (size_t byte = 0U; byte < sizeof(encoded); ++byte) {
    encoded[sizeof(encoded) - 1U - byte] = (unsigned char)(value & 0xffU);
    value >>= 8U;
  }
  gamera_no_sha256_update(sha, encoded, sizeof(encoded));
}

static int digest_double(gamera_no_sha256 *sha, double value) {
  uint64_t bits = 0U;
  if (!isfinite(value) || sizeof(bits) != sizeof(value)) {
    return 0;
  }
  if (value == 0.0) {
    value = 0.0;
  }
  memcpy(&bits, &value, sizeof(bits));
  digest_u64(sha, bits);
  return 1;
}

int gamera_mi_restart_digest(
    const gamera_mi_restart_state *state,
    unsigned char digest[GAMERA_MI_RESTART_DIGEST_BYTES]) {
  static const unsigned char tag[] = "GAMERA_MI_RESTART_V2";
  size_t count = 0U;
  if (digest == NULL || gamera_mi_restart_validate(state) != 0 ||
      !potential_count(state->theta_points, state->azimuth_points, &count)) {
    return -1;
  }
  gamera_no_sha256 sha;
  gamera_no_sha256_init(&sha);
  gamera_no_sha256_update(&sha, tag, sizeof(tag) - 1U);
  digest_u32(&sha, state->schema_version);
  digest_u32(&sha, state->solution_ready);
  digest_u32(&sha, state->snapshot_ready);
  digest_u32(&sha, state->precipitation_state_ready);
  digest_u32(&sha, state->dpb_initialized[0]);
  digest_u32(&sha, state->dpb_initialized[1]);
  digest_u64(&sha, (uint64_t)state->theta_points);
  digest_u64(&sha, (uint64_t)state->azimuth_points);
  if (!digest_double(&sha, state->next_update_s) ||
      !digest_double(&sha, state->next_diagnostic_s) ||
      !digest_double(&sha, state->snapshot_epoch_s) ||
      !digest_double(&sha, state->maximum_cached_emf) ||
      !digest_double(&sha, state->maximum_prehold_difference) ||
      !digest_double(&sha, state->maximum_posthold_difference) ||
      !digest_double(&sha, state->dpb_radius_deg[0]) ||
      !digest_double(&sha, state->dpb_radius_deg[1])) {
    return -1;
  }
  digest_u64(&sha, state->update_count);
  digest_u64(&sha, state->snapshot_generation);
  digest_u64(&sha, state->held_apply_count);
  for (size_t index = 0U; index < count; ++index) {
    if (!digest_double(&sha, state->potential_v[index])) {
      return -1;
    }
  }
  for (size_t index = 0U; index < count; ++index) {
    if (!digest_double(&sha, state->pedersen_siemens[index]) ||
        !digest_double(&sha, state->hall_siemens[index])) {
      return -1;
    }
  }
  gamera_no_sha256_final(&sha, digest);
  return 0;
}

static int write_scalar(hid_t file, const char *path, hid_t file_type,
                        hid_t memory_type, const void *value) {
  hid_t space = H5Screate(H5S_SCALAR);
  hid_t dataset = space < 0 ? -1 : H5Dcreate2(
      file, path, file_type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  const int failed =
      space < 0 || dataset < 0 ||
      H5Dwrite(dataset, memory_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, value) < 0;
  if (dataset >= 0) {
    H5Dclose(dataset);
  }
  if (space >= 0) {
    H5Sclose(space);
  }
  return failed ? -1 : 0;
}

static int write_dataset(hid_t file, const char *path, hid_t file_type,
                         hid_t memory_type, int rank, const hsize_t *extent,
                         const void *value) {
  hid_t space = H5Screate_simple(rank, extent, NULL);
  hid_t dataset = space < 0 ? -1 : H5Dcreate2(
      file, path, file_type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  const int failed =
      space < 0 || dataset < 0 ||
      H5Dwrite(dataset, memory_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, value) < 0;
  if (dataset >= 0) {
    H5Dclose(dataset);
  }
  if (space >= 0) {
    H5Sclose(space);
  }
  return failed ? -1 : 0;
}

static int replace_scalar(hid_t file, const char *path, hid_t memory_type,
                          const void *value) {
  hid_t dataset = H5Dopen2(file, path, H5P_DEFAULT);
  const int failed =
      dataset < 0 ||
      H5Dwrite(dataset, memory_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, value) < 0;
  if (dataset >= 0) {
    H5Dclose(dataset);
  }
  return failed ? -1 : 0;
}

static int read_scalar(hid_t file, const char *path, hid_t memory_type,
                       void *value) {
  hid_t dataset = -1;
  H5E_BEGIN_TRY { dataset = H5Dopen2(file, path, H5P_DEFAULT); }
  H5E_END_TRY;
  hid_t space = dataset < 0 ? -1 : H5Dget_space(dataset);
  int failed = dataset < 0 || space < 0 ||
               H5Sget_simple_extent_type(space) != H5S_SCALAR;
  if (!failed &&
      H5Dread(dataset, memory_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, value) < 0) {
    failed = 1;
  }
  if (space >= 0) {
    H5Sclose(space);
  }
  if (dataset >= 0) {
    H5Dclose(dataset);
  }
  return failed ? -1 : 0;
}

static int read_dataset(hid_t file, const char *path, hid_t memory_type,
                        int expected_rank, const hsize_t *expected_extent,
                        void *value) {
  hid_t dataset = -1;
  H5E_BEGIN_TRY { dataset = H5Dopen2(file, path, H5P_DEFAULT); }
  H5E_END_TRY;
  hid_t space = dataset < 0 ? -1 : H5Dget_space(dataset);
  int failed = dataset < 0 || space < 0 ||
               H5Sget_simple_extent_ndims(space) != expected_rank;
  hsize_t extent[3] = {0U, 0U, 0U};
  if (!failed && H5Sget_simple_extent_dims(space, extent, NULL) < 0) {
    failed = 1;
  }
  for (int axis = 0; !failed && axis < expected_rank; ++axis) {
    failed = extent[axis] != expected_extent[axis];
  }
  if (!failed &&
      H5Dread(dataset, memory_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, value) < 0) {
    failed = 1;
  }
  if (space >= 0) {
    H5Sclose(space);
  }
  if (dataset >= 0) {
    H5Dclose(dataset);
  }
  return failed ? -1 : 0;
}

static int group_has_exact_links(hid_t file, const char *path,
                                 size_t expected_count,
                                 const char *const *expected_name) {
  hid_t group = -1;
  H5E_BEGIN_TRY { group = H5Gopen2(file, path, H5P_DEFAULT); }
  H5E_END_TRY;
  if (group < 0) {
    return 0;
  }
  H5G_info_t info;
  int valid = H5Gget_info(group, &info) >= 0 &&
              info.nlinks == (hsize_t)expected_count;
  for (size_t index = 0U; valid && index < expected_count; ++index) {
    valid = H5Lexists(group, expected_name[index], H5P_DEFAULT) > 0;
  }
  H5Gclose(group);
  return valid;
}

static int exact_mi_schema(hid_t file) {
  static const char *const name[] = {
      "schema_version", "complete", "payload_sha256", "theta_points",
      "azimuth_points", "solution_ready", "snapshot_ready",
      "next_update_s", "next_diagnostic_s", "update_count",
      "snapshot_epoch_s", "snapshot_generation", "maximum_cached_emf",
      "maximum_prehold_difference", "maximum_posthold_difference",
      "held_apply_count", "precipitation_state_ready", "dpb_initialized",
      "dpb_radius_deg", "potential_v", "pedersen_siemens", "hall_siemens"};
  return group_has_exact_links(file, "/mi", sizeof(name) / sizeof(name[0]),
                               name);
}

int gamera_mi_restart_append_hdf5(const char *path,
                                  const gamera_mi_restart_state *state) {
  if (path == NULL || gamera_mi_restart_validate(state) != 0) {
    return -1;
  }
  hid_t file = H5Fopen(path, H5F_ACC_RDWR, H5P_DEFAULT);
  if (file < 0 || H5Lexists(file, "/mi", H5P_DEFAULT) != 0) {
    if (file >= 0) {
      H5Fclose(file);
    }
    return -1;
  }
  hid_t group = H5Gcreate2(file, "/mi", H5P_DEFAULT, H5P_DEFAULT,
                           H5P_DEFAULT);
  if (group < 0) {
    H5Fclose(file);
    return -1;
  }
  H5Gclose(group);
  unsigned char digest[GAMERA_MI_RESTART_DIGEST_BYTES];
  uint32_t complete = 0U;
  uint64_t theta = (uint64_t)state->theta_points;
  uint64_t azimuth = (uint64_t)state->azimuth_points;
  const hsize_t digest_extent[1] = {GAMERA_MI_RESTART_DIGEST_BYTES};
  const hsize_t potential_extent[3] = {MI_HEMISPHERE_COUNT,
                                       (hsize_t)state->theta_points,
                                       (hsize_t)state->azimuth_points};
  const hsize_t hemisphere_extent[1] = {MI_HEMISPHERE_COUNT};
  int failed =
      gamera_mi_restart_digest(state, digest) != 0 ||
      write_scalar(file, "/mi/schema_version", H5T_STD_U32LE,
                   H5T_NATIVE_UINT32, &state->schema_version) != 0 ||
      write_scalar(file, "/mi/complete", H5T_STD_U32LE, H5T_NATIVE_UINT32,
                   &complete) != 0 ||
      write_dataset(file, "/mi/payload_sha256", H5T_STD_U8LE,
                    H5T_NATIVE_UCHAR, 1, digest_extent, digest) != 0 ||
      write_scalar(file, "/mi/theta_points", H5T_STD_U64LE,
                   H5T_NATIVE_UINT64, &theta) != 0 ||
      write_scalar(file, "/mi/azimuth_points", H5T_STD_U64LE,
                   H5T_NATIVE_UINT64, &azimuth) != 0 ||
      write_scalar(file, "/mi/solution_ready", H5T_STD_U32LE,
                   H5T_NATIVE_UINT32, &state->solution_ready) != 0 ||
      write_scalar(file, "/mi/snapshot_ready", H5T_STD_U32LE,
                   H5T_NATIVE_UINT32, &state->snapshot_ready) != 0 ||
      write_scalar(file, "/mi/next_update_s", H5T_IEEE_F64LE,
                   H5T_NATIVE_DOUBLE, &state->next_update_s) != 0 ||
      write_scalar(file, "/mi/next_diagnostic_s", H5T_IEEE_F64LE,
                   H5T_NATIVE_DOUBLE, &state->next_diagnostic_s) != 0 ||
      write_scalar(file, "/mi/update_count", H5T_STD_U64LE,
                   H5T_NATIVE_UINT64, &state->update_count) != 0 ||
      write_scalar(file, "/mi/snapshot_epoch_s", H5T_IEEE_F64LE,
                   H5T_NATIVE_DOUBLE, &state->snapshot_epoch_s) != 0 ||
      write_scalar(file, "/mi/snapshot_generation", H5T_STD_U64LE,
                   H5T_NATIVE_UINT64, &state->snapshot_generation) != 0 ||
      write_scalar(file, "/mi/maximum_cached_emf", H5T_IEEE_F64LE,
                   H5T_NATIVE_DOUBLE, &state->maximum_cached_emf) != 0 ||
      write_scalar(file, "/mi/maximum_prehold_difference", H5T_IEEE_F64LE,
                   H5T_NATIVE_DOUBLE,
                   &state->maximum_prehold_difference) != 0 ||
      write_scalar(file, "/mi/maximum_posthold_difference", H5T_IEEE_F64LE,
                   H5T_NATIVE_DOUBLE,
                   &state->maximum_posthold_difference) != 0 ||
      write_scalar(file, "/mi/held_apply_count", H5T_STD_U64LE,
                   H5T_NATIVE_UINT64, &state->held_apply_count) != 0 ||
      write_scalar(file, "/mi/precipitation_state_ready", H5T_STD_U32LE,
                   H5T_NATIVE_UINT32,
                   &state->precipitation_state_ready) != 0 ||
      write_dataset(file, "/mi/dpb_initialized", H5T_STD_U32LE,
                    H5T_NATIVE_UINT32, 1, hemisphere_extent,
                    state->dpb_initialized) != 0 ||
      write_dataset(file, "/mi/dpb_radius_deg", H5T_IEEE_F64LE,
                    H5T_NATIVE_DOUBLE, 1, hemisphere_extent,
                    state->dpb_radius_deg) != 0 ||
      write_dataset(file, "/mi/potential_v", H5T_IEEE_F64LE,
                    H5T_NATIVE_DOUBLE, 3, potential_extent,
                    state->potential_v) != 0 ||
      write_dataset(file, "/mi/pedersen_siemens", H5T_IEEE_F64LE,
                    H5T_NATIVE_DOUBLE, 3, potential_extent,
                    state->pedersen_siemens) != 0 ||
      write_dataset(file, "/mi/hall_siemens", H5T_IEEE_F64LE,
                    H5T_NATIVE_DOUBLE, 3, potential_extent,
                    state->hall_siemens) != 0;
  if (!failed) {
    complete = 1U;
    failed = replace_scalar(file, "/mi/complete", H5T_NATIVE_UINT32,
                            &complete) != 0 ||
             H5Fflush(file, H5F_SCOPE_GLOBAL) < 0;
  }
  const int close_failed = H5Fclose(file) < 0;
  int sync_failed = 0;
  if (!failed && !close_failed) {
    const int descriptor = open(path, O_RDONLY);
    if (descriptor < 0) {
      sync_failed = 1;
    } else {
      const int flush_failed = fsync(descriptor) != 0;
      const int descriptor_close_failed = close(descriptor) != 0;
      sync_failed = flush_failed || descriptor_close_failed;
    }
  }
  return failed || close_failed || sync_failed ? -1 : 0;
}

int gamera_mi_restart_hdf5_present(const char *path) {
  if (path == NULL) {
    return -1;
  }
  hid_t file = -1;
  H5E_BEGIN_TRY { file = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT); }
  H5E_END_TRY;
  if (file < 0) {
    return -1;
  }
  const htri_t present = H5Lexists(file, "/mi", H5P_DEFAULT);
  const int close_failed = H5Fclose(file) < 0;
  if (present < 0 || close_failed) {
    return -1;
  }
  return present > 0 ? 1 : 0;
}

int gamera_mi_restart_read_hdf5(const char *path,
                                gamera_mi_restart_state **state) {
  if (state == NULL) {
    return -1;
  }
  *state = NULL;
  if (path == NULL) {
    return -1;
  }
  hid_t file = -1;
  H5E_BEGIN_TRY { file = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT); }
  H5E_END_TRY;
  if (file < 0 || !exact_mi_schema(file)) {
    if (file >= 0) {
      H5Fclose(file);
    }
    return -1;
  }
  uint32_t schema = 0U;
  uint32_t complete = 0U;
  uint64_t theta = 0U;
  uint64_t azimuth = 0U;
  int failed =
      read_scalar(file, "/mi/schema_version", H5T_NATIVE_UINT32, &schema) !=
          0 ||
      read_scalar(file, "/mi/complete", H5T_NATIVE_UINT32, &complete) != 0 ||
      read_scalar(file, "/mi/theta_points", H5T_NATIVE_UINT64, &theta) != 0 ||
      read_scalar(file, "/mi/azimuth_points", H5T_NATIVE_UINT64,
                  &azimuth) != 0 ||
      schema != GAMERA_MI_RESTART_SCHEMA || complete != 1U ||
      theta > SIZE_MAX || azimuth > SIZE_MAX;
  gamera_mi_restart_state *candidate = NULL;
  if (!failed &&
      gamera_mi_restart_create((size_t)theta, (size_t)azimuth, &candidate) !=
          0) {
    failed = 1;
  }
  unsigned char stored[GAMERA_MI_RESTART_DIGEST_BYTES];
  unsigned char computed[GAMERA_MI_RESTART_DIGEST_BYTES];
  if (!failed) {
    const hsize_t digest_extent[1] = {GAMERA_MI_RESTART_DIGEST_BYTES};
    const hsize_t potential_extent[3] = {MI_HEMISPHERE_COUNT,
                                         (hsize_t)candidate->theta_points,
                                         (hsize_t)candidate->azimuth_points};
    const hsize_t hemisphere_extent[1] = {MI_HEMISPHERE_COUNT};
    failed =
        read_scalar(file, "/mi/solution_ready", H5T_NATIVE_UINT32,
                    &candidate->solution_ready) != 0 ||
        read_scalar(file, "/mi/snapshot_ready", H5T_NATIVE_UINT32,
                    &candidate->snapshot_ready) != 0 ||
        read_scalar(file, "/mi/next_update_s", H5T_NATIVE_DOUBLE,
                    &candidate->next_update_s) != 0 ||
        read_scalar(file, "/mi/next_diagnostic_s", H5T_NATIVE_DOUBLE,
                    &candidate->next_diagnostic_s) != 0 ||
        read_scalar(file, "/mi/update_count", H5T_NATIVE_UINT64,
                    &candidate->update_count) != 0 ||
        read_scalar(file, "/mi/snapshot_epoch_s", H5T_NATIVE_DOUBLE,
                    &candidate->snapshot_epoch_s) != 0 ||
        read_scalar(file, "/mi/snapshot_generation", H5T_NATIVE_UINT64,
                    &candidate->snapshot_generation) != 0 ||
        read_scalar(file, "/mi/maximum_cached_emf", H5T_NATIVE_DOUBLE,
                    &candidate->maximum_cached_emf) != 0 ||
        read_scalar(file, "/mi/maximum_prehold_difference",
                    H5T_NATIVE_DOUBLE,
                    &candidate->maximum_prehold_difference) != 0 ||
        read_scalar(file, "/mi/maximum_posthold_difference",
                    H5T_NATIVE_DOUBLE,
                    &candidate->maximum_posthold_difference) != 0 ||
        read_scalar(file, "/mi/held_apply_count", H5T_NATIVE_UINT64,
                    &candidate->held_apply_count) != 0 ||
        read_scalar(file, "/mi/precipitation_state_ready",
                    H5T_NATIVE_UINT32,
                    &candidate->precipitation_state_ready) != 0 ||
        read_dataset(file, "/mi/dpb_initialized", H5T_NATIVE_UINT32, 1,
                     hemisphere_extent, candidate->dpb_initialized) != 0 ||
        read_dataset(file, "/mi/dpb_radius_deg", H5T_NATIVE_DOUBLE, 1,
                     hemisphere_extent, candidate->dpb_radius_deg) != 0 ||
        read_dataset(file, "/mi/potential_v", H5T_NATIVE_DOUBLE, 3,
                     potential_extent, candidate->potential_v) != 0 ||
        read_dataset(file, "/mi/pedersen_siemens", H5T_NATIVE_DOUBLE, 3,
                     potential_extent, candidate->pedersen_siemens) != 0 ||
        read_dataset(file, "/mi/hall_siemens", H5T_NATIVE_DOUBLE, 3,
                     potential_extent, candidate->hall_siemens) != 0 ||
        read_dataset(file, "/mi/payload_sha256", H5T_NATIVE_UCHAR, 1,
                     digest_extent, stored) != 0 ||
        gamera_mi_restart_digest(candidate, computed) != 0 ||
        memcmp(stored, computed, sizeof(stored)) != 0;
  }
  if (H5Fclose(file) < 0) {
    failed = 1;
  }
  if (failed) {
    gamera_mi_restart_destroy(candidate);
    return -1;
  }
  *state = candidate;
  return 0;
}

int gamera_mi_restart_broadcast(MPI_Comm communicator, int root,
                                gamera_mi_restart_state **state) {
  int rank = -1;
  int size = 0;
  if (state == NULL || MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || root < 0 ||
      root >= size) {
    return -1;
  }
  uint64_t integer[11] = {0U};
  double real[8] = {0.0};
  int local_failed = rank == root && gamera_mi_restart_validate(*state) != 0;
  int global_failed = 0;
  if (MPI_Allreduce(&local_failed, &global_failed, 1, MPI_INT, MPI_MAX,
                    communicator) != MPI_SUCCESS ||
      global_failed) {
    return -1;
  }
  if (rank == root) {
    integer[0] = (*state)->schema_version;
    integer[1] = (*state)->solution_ready;
    integer[2] = (*state)->snapshot_ready;
    integer[3] = (uint64_t)(*state)->theta_points;
    integer[4] = (uint64_t)(*state)->azimuth_points;
    integer[5] = (*state)->update_count;
    integer[6] = (*state)->snapshot_generation;
    integer[7] = (*state)->held_apply_count;
    integer[8] = (*state)->precipitation_state_ready;
    integer[9] = (*state)->dpb_initialized[0];
    integer[10] = (*state)->dpb_initialized[1];
    real[0] = (*state)->next_update_s;
    real[1] = (*state)->next_diagnostic_s;
    real[2] = (*state)->snapshot_epoch_s;
    real[3] = (*state)->maximum_cached_emf;
    real[4] = (*state)->maximum_prehold_difference;
    real[5] = (*state)->maximum_posthold_difference;
    real[6] = (*state)->dpb_radius_deg[0];
    real[7] = (*state)->dpb_radius_deg[1];
  }
  if (MPI_Bcast(integer, 11, MPI_UINT64_T, root, communicator) != MPI_SUCCESS ||
      MPI_Bcast(real, 8, MPI_DOUBLE, root, communicator) != MPI_SUCCESS) {
    return -1;
  }
  if (rank != root) {
    if (integer[3] > SIZE_MAX || integer[4] > SIZE_MAX ||
        gamera_mi_restart_create((size_t)integer[3], (size_t)integer[4],
                                 state) != 0) {
      local_failed = 1;
    } else {
      (*state)->schema_version = (uint32_t)integer[0];
      (*state)->solution_ready = (uint32_t)integer[1];
      (*state)->snapshot_ready = (uint32_t)integer[2];
      (*state)->update_count = integer[5];
      (*state)->snapshot_generation = integer[6];
      (*state)->held_apply_count = integer[7];
      (*state)->precipitation_state_ready = (uint32_t)integer[8];
      (*state)->dpb_initialized[0] = (uint32_t)integer[9];
      (*state)->dpb_initialized[1] = (uint32_t)integer[10];
      (*state)->next_update_s = real[0];
      (*state)->next_diagnostic_s = real[1];
      (*state)->snapshot_epoch_s = real[2];
      (*state)->maximum_cached_emf = real[3];
      (*state)->maximum_prehold_difference = real[4];
      (*state)->maximum_posthold_difference = real[5];
      (*state)->dpb_radius_deg[0] = real[6];
      (*state)->dpb_radius_deg[1] = real[7];
    }
  }
  if (MPI_Allreduce(&local_failed, &global_failed, 1, MPI_INT, MPI_MAX,
                    communicator) != MPI_SUCCESS ||
      global_failed) {
    if (rank != root) {
      gamera_mi_restart_destroy(*state);
      *state = NULL;
    }
    return -1;
  }
  size_t count = 0U;
  if (!potential_count((*state)->theta_points, (*state)->azimuth_points,
                       &count) ||
      count > (size_t)INT_MAX ||
      MPI_Bcast((*state)->potential_v, (int)count, MPI_DOUBLE, root,
                communicator) != MPI_SUCCESS ||
      MPI_Bcast((*state)->pedersen_siemens, (int)count, MPI_DOUBLE, root,
                communicator) != MPI_SUCCESS ||
      MPI_Bcast((*state)->hall_siemens, (int)count, MPI_DOUBLE, root,
                communicator) != MPI_SUCCESS ||
      gamera_mi_restart_validate(*state) != 0) {
    if (rank != root) {
      gamera_mi_restart_destroy(*state);
      *state = NULL;
    }
    return -1;
  }
  return 0;
}
