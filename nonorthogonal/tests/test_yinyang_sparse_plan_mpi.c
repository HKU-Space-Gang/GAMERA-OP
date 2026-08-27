#include "nonorthogonal_yinyang_sparse_plan.h"

#include <math.h>
#include <mpi.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  int world_rank;
} pack_context;

static double expected_value(int donor_rank, size_t donor_cell, int channel) {
  return (double)donor_rank * 100000.0 + (double)donor_cell * 100.0 +
         (double)channel;
}

static int pack_values(const size_t *cells, size_t cell_count,
                       int channel_count, double *values, void *opaque) {
  const pack_context *context = (const pack_context *)opaque;
  if (context == NULL || (cell_count > 0 && (cells == NULL || values == NULL))) {
    return -1;
  }
  for (size_t item = 0; item < cell_count; ++item) {
    for (int channel = 0; channel < channel_count; ++channel) {
      values[item * (size_t)channel_count + (size_t)channel] =
          expected_value(context->world_rank, cells[item], channel);
    }
  }
  return 0;
}

static int pack_fail_on_rank_zero(const size_t *cells, size_t cell_count,
                                  int channel_count, double *values,
                                  void *opaque) {
  const pack_context *context = (const pack_context *)opaque;
  if (context != NULL && context->world_rank == 0) {
    return -1;
  }
  return pack_values(cells, cell_count, channel_count, values, opaque);
}

static int test_deduplicated_exchange(int rank, int size) {
  const int patch_size = size / 2;
  const int patch = rank / patch_size;
  const int patch_rank = rank % patch_size;
  const int allowed_begin = patch == 0 ? patch_size : 0;
  const int allowed_end = patch == 0 ? size : patch_size;
  const int donor_a = allowed_begin + patch_rank;
  const int donor_b =
      allowed_begin + ((patch_rank + 1) % patch_size);
  const gamera_no_sparse_reference references[6] = {
      {donor_a, 0}, {donor_a, 3}, {donor_a, 0},
      {donor_b, 5}, {donor_b, 7}, {donor_b, 5}};
  gamera_no_sparse_plan plan = GAMERA_NO_SPARSE_PLAN_INITIALIZER;
  int failed = gamera_no_sparse_plan_build(
                   &plan, MPI_COMM_WORLD, references, 6, 16, 3,
                   allowed_begin, allowed_end) != 0;
  pack_context context = {rank};
  if (!failed) {
    failed = gamera_no_sparse_plan_exchange(&plan, pack_values, &context) != 0;
  }
  for (size_t item = 0; !failed && item < 6; ++item) {
    for (int channel = 0; channel < 3; ++channel) {
      const double actual = gamera_no_sparse_plan_value(&plan, item, channel);
      const double expected = expected_value(
          references[item].donor_rank, references[item].donor_cell, channel);
      failed |= actual != expected;
    }
  }
  failed |= !plan.ready || plan.reference_count != 6 ||
            plan.receive_cell_count != 4U ||
            plan.call_count != 1 ||
            plan.useful_receive_bytes != 4U * 3U * sizeof(double);
  failed |= plan.reference_slot[0] != plan.reference_slot[2] ||
            plan.reference_slot[3] != plan.reference_slot[5];
  gamera_no_sparse_plan_destroy(&plan);
  return failed;
}

static int test_zero_degree_plan(int rank, int size) {
  const int patch_size = size / 2;
  const int patch = rank / patch_size;
  const int allowed_begin = patch == 0 ? patch_size : 0;
  const int allowed_end = patch == 0 ? size : patch_size;
  gamera_no_sparse_plan plan = GAMERA_NO_SPARSE_PLAN_INITIALIZER;
  int failed = gamera_no_sparse_plan_build(
                   &plan, MPI_COMM_WORLD, NULL, 0, 16, 2,
                   allowed_begin, allowed_end) != 0;
  pack_context context = {rank};
  if (!failed) {
    failed = gamera_no_sparse_plan_exchange(&plan, pack_values, &context) != 0;
  }
  failed |= !plan.ready || plan.indegree != 0 || plan.outdegree != 0 ||
            plan.receive_cell_count != 0 || plan.send_cell_count != 0 ||
            plan.call_count != 1 || plan.useful_send_bytes != 0 ||
            plan.useful_receive_bytes != 0;
  gamera_no_sparse_plan_destroy(&plan);
  return failed;
}

static int test_mixed_zero_degree(int rank, int size) {
  const int patch_size = size / 2;
  const int patch = rank / patch_size;
  const int allowed_begin = patch == 0 ? patch_size : 0;
  const int allowed_end = patch == 0 ? size : patch_size;
  gamera_no_sparse_reference reference = {patch_size, 11};
  const size_t reference_count = rank == 0 ? 1U : 0U;
  const gamera_no_sparse_reference *references =
      rank == 0 ? &reference : NULL;
  gamera_no_sparse_plan plan = GAMERA_NO_SPARSE_PLAN_INITIALIZER;
  int failed = gamera_no_sparse_plan_build(
                   &plan, MPI_COMM_WORLD, references, reference_count, 16, 2,
                   allowed_begin, allowed_end) != 0;
  pack_context context = {rank};
  if (!failed) {
    failed = gamera_no_sparse_plan_exchange(&plan, pack_values, &context) != 0;
  }
  if (!failed && rank == 0) {
    for (int channel = 0; channel < 2; ++channel) {
      failed |= gamera_no_sparse_plan_value(&plan, 0, channel) !=
                expected_value(patch_size, 11, channel);
    }
    failed |= plan.indegree != 1 || plan.receive_cell_count != 1;
  }
  if (!failed && rank == patch_size) {
    failed |= plan.outdegree != 1 || plan.send_cell_count != 1;
  }
  if (!failed && rank != 0 && rank != patch_size) {
    failed |= plan.indegree != 0 || plan.outdegree != 0;
  }
  gamera_no_sparse_plan_destroy(&plan);
  return failed;
}

static int test_collective_invalid_cell(int rank, int size) {
  const int patch_size = size / 2;
  const int patch = rank / patch_size;
  const int allowed_begin = patch == 0 ? patch_size : 0;
  const int allowed_end = patch == 0 ? size : patch_size;
  gamera_no_sparse_reference reference = {patch_size, 99};
  const size_t reference_count = rank == 0 ? 1U : 0U;
  gamera_no_sparse_plan plan = GAMERA_NO_SPARSE_PLAN_INITIALIZER;
  const int unexpectedly_succeeded =
      gamera_no_sparse_plan_build(
          &plan, MPI_COMM_WORLD, rank == 0 ? &reference : NULL,
          reference_count, 16, 2, allowed_begin, allowed_end) == 0;
  gamera_no_sparse_plan_destroy(&plan);
  return unexpectedly_succeeded;
}

static int test_collective_invalid_rank(int rank, int size) {
  const int patch_size = size / 2;
  const int patch = rank / patch_size;
  const int allowed_begin = patch == 0 ? patch_size : 0;
  const int allowed_end = patch == 0 ? size : patch_size;
  gamera_no_sparse_reference reference = {0, 1};
  const size_t reference_count = rank == 0 ? 1U : 0U;
  gamera_no_sparse_plan plan = GAMERA_NO_SPARSE_PLAN_INITIALIZER;
  const int unexpectedly_succeeded =
      gamera_no_sparse_plan_build(
          &plan, MPI_COMM_WORLD, rank == 0 ? &reference : NULL,
          reference_count, 16, 2, allowed_begin, allowed_end) == 0;
  gamera_no_sparse_plan_destroy(&plan);
  return unexpectedly_succeeded;
}

static int test_collective_invalid_metadata(int rank, int size) {
  const int patch_size = size / 2;
  const int patch = rank / patch_size;
  const int allowed_begin = patch == 0 ? patch_size : 0;
  const int allowed_end = patch == 0 ? size : patch_size;
  gamera_no_sparse_plan plan = GAMERA_NO_SPARSE_PLAN_INITIALIZER;
  const int channel_count = rank == 0 ? 0 : 2;
  const int unexpectedly_succeeded =
      gamera_no_sparse_plan_build(&plan, MPI_COMM_WORLD, NULL, 0, 16,
                                  channel_count, allowed_begin,
                                  allowed_end) == 0;
  gamera_no_sparse_plan_destroy(&plan);
  /* Destruction is intentionally idempotent after collective build failure. */
  gamera_no_sparse_plan_destroy(&plan);
  return unexpectedly_succeeded;
}

static int test_collective_pack_failure(int rank, int size) {
  const int patch_size = size / 2;
  const int patch = rank / patch_size;
  const int allowed_begin = patch == 0 ? patch_size : 0;
  const int allowed_end = patch == 0 ? size : patch_size;
  gamera_no_sparse_plan plan = GAMERA_NO_SPARSE_PLAN_INITIALIZER;
  int failed = gamera_no_sparse_plan_build(
                   &plan, MPI_COMM_WORLD, NULL, 0, 16, 2,
                   allowed_begin, allowed_end) != 0;
  pack_context context = {rank};
  if (!failed) {
    const int exchange_status = gamera_no_sparse_plan_exchange(
        &plan, pack_fail_on_rank_zero, &context);
    failed = exchange_status == 0 || plan.call_count != 0;
  }
  gamera_no_sparse_plan_destroy(&plan);
  return failed;
}

static int test_statistics_saturation(int rank, int size) {
  const int patch_size = size / 2;
  const int patch = rank / patch_size;
  const int allowed_begin = patch == 0 ? patch_size : 0;
  const int allowed_end = patch == 0 ? size : patch_size;
  const int donor = allowed_begin + rank % patch_size;
  const gamera_no_sparse_reference reference = {donor, 1};
  gamera_no_sparse_plan plan = GAMERA_NO_SPARSE_PLAN_INITIALIZER;
  int failed = gamera_no_sparse_plan_build(
                   &plan, MPI_COMM_WORLD, &reference, 1, 16, 2,
                   allowed_begin, allowed_end) != 0;
  pack_context context = {rank};
  if (!failed) {
    plan.call_count = UINT64_MAX;
    plan.useful_send_bytes = UINT64_MAX;
    plan.useful_receive_bytes = UINT64_MAX;
    failed = gamera_no_sparse_plan_exchange(&plan, pack_values, &context) != 0;
  }
  failed |= plan.call_count != UINT64_MAX ||
            plan.useful_send_bytes != UINT64_MAX ||
            plan.useful_receive_bytes != UINT64_MAX ||
            !plan.statistics_saturated;
  gamera_no_sparse_plan_destroy(&plan);
  return failed;
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
  int rank = -1;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (argc == 2 && strcmp(argv[1], "--abort-invalid-exchange") == 0) {
    const int patch_size = size / 2;
    const int patch = rank / patch_size;
    const int allowed_begin = patch == 0 ? patch_size : 0;
    const int allowed_end = patch == 0 ? size : patch_size;
    gamera_no_sparse_plan plan = GAMERA_NO_SPARSE_PLAN_INITIALIZER;
    if (gamera_no_sparse_plan_build(&plan, MPI_COMM_WORLD, NULL, 0, 16, 2,
                                    allowed_begin, allowed_end) == 0) {
      pack_context context = {rank};
      if (rank == 0) {
        plan.ready = 0;
      }
      (void)gamera_no_sparse_plan_exchange(&plan, pack_values, &context);
    }
    MPI_Finalize();
    return EXIT_SUCCESS;
  }
  if (argc == 2 && strcmp(argv[1], "--abort-null-build") == 0) {
    const int patch_size = size / 2;
    const int patch = rank / patch_size;
    const int allowed_begin = patch == 0 ? patch_size : 0;
    const int allowed_end = patch == 0 ? size : patch_size;
    gamera_no_sparse_plan plan = GAMERA_NO_SPARSE_PLAN_INITIALIZER;
    (void)gamera_no_sparse_plan_build(
        rank == 0 ? NULL : &plan, MPI_COMM_WORLD, NULL, 0, 16, 2,
        allowed_begin, allowed_end);
    MPI_Finalize();
    return EXIT_SUCCESS;
  }
  if (argc == 2 && strcmp(argv[1], "--abort-reused-build") == 0) {
    const int patch_size = size / 2;
    const int patch = rank / patch_size;
    const int allowed_begin = patch == 0 ? patch_size : 0;
    const int allowed_end = patch == 0 ? size : patch_size;
    gamera_no_sparse_plan plan = GAMERA_NO_SPARSE_PLAN_INITIALIZER;
    if (gamera_no_sparse_plan_build(&plan, MPI_COMM_WORLD, NULL, 0, 16, 2,
                                    allowed_begin, allowed_end) == 0) {
      if (rank == 0) {
        (void)gamera_no_sparse_plan_build(
            &plan, MPI_COMM_WORLD, NULL, 0, 16, 2, allowed_begin,
            allowed_end);
      } else {
        pack_context context = {rank};
        (void)gamera_no_sparse_plan_exchange(&plan, pack_values, &context);
      }
    }
    MPI_Finalize();
    return EXIT_SUCCESS;
  }
  int failed = size < 2 || size % 2 != 0;
  if (!failed) {
    failed |= test_deduplicated_exchange(rank, size);
    failed |= test_zero_degree_plan(rank, size);
    failed |= test_mixed_zero_degree(rank, size);
    failed |= test_collective_invalid_cell(rank, size);
    failed |= test_collective_invalid_rank(rank, size);
    failed |= test_collective_invalid_metadata(rank, size);
    failed |= test_collective_pack_failure(rank, size);
    failed |= test_statistics_saturation(rank, size);
  }
  int any_failed = 0;
  MPI_Allreduce(&failed, &any_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (rank == 0) {
    printf("sparse-plan MPI test: %s (ranks=%d)\n",
           any_failed ? "FAIL" : "PASS", size);
  }
  MPI_Finalize();
  return any_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
