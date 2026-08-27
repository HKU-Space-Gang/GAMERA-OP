#include "nonorthogonal_yinyang_sparse_plan.h"

#ifdef GAMERA_YINYANG_SPARSE_OVERSET

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#if MPI_VERSION < 3
#error "GAMERA_YINYANG_SPARSE_OVERSET requires MPI-3 distributed graph collectives"
#endif

typedef struct {
  int rank;
  size_t cell;
  size_t original;
} sparse_record;

static int sparse_record_compare(const void *left_pointer,
                                 const void *right_pointer) {
  const sparse_record *left = (const sparse_record *)left_pointer;
  const sparse_record *right = (const sparse_record *)right_pointer;
  if (left->rank != right->rank) {
    return left->rank < right->rank ? -1 : 1;
  }
  if (left->cell != right->cell) {
    return left->cell < right->cell ? -1 : 1;
  }
  if (left->original != right->original) {
    return left->original < right->original ? -1 : 1;
  }
  return 0;
}

static int sparse_collective_mpi_or_abort(MPI_Comm world_comm,
                                          int mpi_status) {
  if (mpi_status == MPI_SUCCESS) {
    return 0;
  }
  /*
   * MPI cannot promise that a failed collective returned on every peer or
   * left its communicator usable.  Returning into rank-divergent cleanup can
   * therefore deadlock or free a communicator twice.  Make that failure order
   * explicit and terminal; normal/local failures still use consensus gates.
   */
  MPI_Comm abort_comm =
      world_comm == MPI_COMM_NULL ? MPI_COMM_WORLD : world_comm;
  MPI_Abort(abort_comm, mpi_status);
  return -1;
}

static int sparse_programmer_error_or_abort(MPI_Comm world_comm) {
  return sparse_collective_mpi_or_abort(world_comm, MPI_ERR_ARG);
}

static int sparse_consensus(MPI_Comm communicator, int local_failed) {
  int any_failed = 0;
  const int mpi_status = MPI_Allreduce(&local_failed, &any_failed, 1, MPI_INT,
                                       MPI_MAX, communicator);
  if (sparse_collective_mpi_or_abort(communicator, mpi_status) != 0) {
    return -1;
  }
  return any_failed ? -1 : 0;
}

static int sparse_comm_free_or_abort(MPI_Comm world_comm,
                                     MPI_Comm *communicator) {
  if (communicator == NULL || *communicator == MPI_COMM_NULL) {
    return 0;
  }
  return sparse_collective_mpi_or_abort(world_comm,
                                        MPI_Comm_free(communicator));
}

static void sparse_accumulate_u64(uint64_t increment, uint64_t *total,
                                  int *saturated) {
  if (total == NULL || saturated == NULL) {
    return;
  }
  if (*total > UINT64_MAX - increment) {
    *total = UINT64_MAX;
    *saturated = 1;
  } else {
    *total += increment;
  }
}

static int sparse_multiply_size(size_t left, size_t right, size_t *result) {
  if (result == NULL || (right != 0 && left > SIZE_MAX / right)) {
    return -1;
  }
  *result = left * right;
  return 0;
}

static void *sparse_calloc(size_t count, size_t element_size) {
  const size_t nonzero_count = count == 0 ? 1 : count;
  if (element_size != 0 && nonzero_count > SIZE_MAX / element_size) {
    return NULL;
  }
  return calloc(nonzero_count, element_size);
}

static int sparse_rank_position(const int *ranks, int count, int target) {
  for (int index = 0; index < count; ++index) {
    if (ranks[index] == target) {
      return index;
    }
  }
  return -1;
}

static size_t sparse_unique_rank_begin(const int *unique_ranks,
                                       size_t unique_count, int target) {
  size_t lower = 0;
  size_t upper = unique_count;
  while (lower < upper) {
    const size_t middle = lower + (upper - lower) / 2;
    if (unique_ranks[middle] < target) {
      lower = middle + 1;
    } else {
      upper = middle;
    }
  }
  return lower;
}

static size_t sparse_unique_rank_end(const int *unique_ranks,
                                     size_t unique_count, int target) {
  size_t lower = sparse_unique_rank_begin(unique_ranks, unique_count, target);
  size_t upper = unique_count;
  while (lower < upper) {
    const size_t middle = lower + (upper - lower) / 2;
    if (unique_ranks[middle] <= target) {
      lower = middle + 1;
    } else {
      upper = middle;
    }
  }
  return lower;
}

static int sparse_accumulate_int_displacements(const int *counts, int count,
                                                int *displacements,
                                                size_t *total) {
  size_t used = 0;
  if (counts == NULL || displacements == NULL || total == NULL || count < 0) {
    return -1;
  }
  for (int index = 0; index < count; ++index) {
    if (counts[index] < 0 || used > (size_t)INT_MAX) {
      return -1;
    }
    displacements[index] = (int)used;
    if ((size_t)counts[index] > (size_t)INT_MAX - used) {
      return -1;
    }
    used += (size_t)counts[index];
  }
  *total = used;
  return 0;
}

static void sparse_release_arrays(gamera_no_sparse_plan *plan) {
  if (plan == NULL) {
    return;
  }
  free(plan->reference_slot);
  free(plan->send_cells);
  free(plan->send_values);
  free(plan->receive_values);
  free(plan->send_counts);
  free(plan->send_displacements);
  free(plan->receive_counts);
  free(plan->receive_displacements);
  plan->reference_slot = NULL;
  plan->send_cells = NULL;
  plan->send_values = NULL;
  plan->receive_values = NULL;
  plan->send_counts = NULL;
  plan->send_displacements = NULL;
  plan->receive_counts = NULL;
  plan->receive_displacements = NULL;
}

int gamera_no_sparse_plan_build(
    gamera_no_sparse_plan *plan, MPI_Comm world_comm,
    const gamera_no_sparse_reference *references, size_t reference_count,
    size_t local_active_cells, int channel_count, int allowed_rank_begin,
    int allowed_rank_end) {
  MPI_Comm request_comm = MPI_COMM_NULL;
  MPI_Comm value_comm = MPI_COMM_NULL;
  int request_comm_ready = 0;
  int value_comm_ready = 0;
  sparse_record *records = NULL;
  int *unique_ranks = NULL;
  uint64_t *unique_cells = NULL;
  size_t *reference_unique = NULL;
  size_t *unique_final_slot = NULL;
  int *donor_ranks = NULL;
  int *request_sources = NULL;
  int *request_destinations = NULL;
  int *request_send_counts = NULL;
  int *request_send_displacements = NULL;
  int *request_receive_counts = NULL;
  int *request_receive_displacements = NULL;
  uint64_t *request_send_cells = NULL;
  uint64_t *request_receive_cells = NULL;
  int *value_sources = NULL;
  int *value_destinations = NULL;
  int local_failed = 0;
  int world_rank = -1;
  int world_size = 0;
  size_t unique_count = 0;
  size_t donor_rank_count_size = 0;
  int donor_rank_count = 0;
  int request_indegree = 0;
  int request_outdegree = 0;
  int value_indegree = 0;
  int value_outdegree = 0;
  size_t request_receive_total = 0;

  if (plan == NULL || world_comm == MPI_COMM_NULL) {
    return sparse_programmer_error_or_abort(world_comm);
  }
  if (plan->ready || plan->value_comm != MPI_COMM_NULL) {
    /* Rebuilding a live plan would discard a collective communicator handle
     * during failure cleanup.  Treat this as a lifetime invariant violation,
     * not as a recoverable metadata error. */
    return sparse_programmer_error_or_abort(world_comm);
  }
  if (sparse_collective_mpi_or_abort(
          world_comm, MPI_Comm_rank(world_comm, &world_rank)) != 0 ||
      sparse_collective_mpi_or_abort(
          world_comm, MPI_Comm_size(world_comm, &world_size)) != 0) {
    return -1;
  }
  local_failed |= channel_count <= 0 || local_active_cells == 0;
  local_failed |= reference_count > 0 && references == NULL;
  local_failed |= allowed_rank_begin < 0 ||
                  allowed_rank_begin >= allowed_rank_end ||
                  allowed_rank_end > world_size;
  if (sparse_consensus(world_comm, local_failed) != 0) {
    goto fail;
  }

  records = (sparse_record *)sparse_calloc(reference_count, sizeof(*records));
  unique_ranks =
      (int *)sparse_calloc(reference_count, sizeof(*unique_ranks));
  unique_cells =
      (uint64_t *)sparse_calloc(reference_count, sizeof(*unique_cells));
  reference_unique =
      (size_t *)sparse_calloc(reference_count, sizeof(*reference_unique));
  unique_final_slot =
      (size_t *)sparse_calloc(reference_count, sizeof(*unique_final_slot));
  local_failed = records == NULL || unique_ranks == NULL ||
                 unique_cells == NULL || reference_unique == NULL ||
                 unique_final_slot == NULL;
  for (size_t item = 0; !local_failed && item < reference_count; ++item) {
    const int donor_rank = references[item].donor_rank;
    const size_t donor_cell = references[item].donor_cell;
    local_failed |= donor_rank < allowed_rank_begin ||
                    donor_rank >= allowed_rank_end || donor_rank == world_rank;
#if SIZE_MAX > UINT64_MAX
    local_failed |= donor_cell > (size_t)UINT64_MAX;
#endif
    records[item].rank = donor_rank;
    records[item].cell = donor_cell;
    records[item].original = item;
  }
  if (sparse_consensus(world_comm, local_failed) != 0) {
    goto fail;
  }

  qsort(records, reference_count, sizeof(*records), sparse_record_compare);
  for (size_t item = 0; item < reference_count; ++item) {
    if (item == 0 || records[item].rank != records[item - 1].rank ||
        records[item].cell != records[item - 1].cell) {
      unique_ranks[unique_count] = records[item].rank;
      unique_cells[unique_count] = (uint64_t)records[item].cell;
      ++unique_count;
    }
    reference_unique[records[item].original] = unique_count - 1;
  }
  for (size_t item = 0; item < unique_count; ++item) {
    if (item == 0 || unique_ranks[item] != unique_ranks[item - 1]) {
      ++donor_rank_count_size;
    }
  }
  local_failed = donor_rank_count_size > (size_t)INT_MAX;
  donor_rank_count = local_failed ? 0 : (int)donor_rank_count_size;
  donor_ranks = (int *)sparse_calloc(donor_rank_count_size,
                                     sizeof(*donor_ranks));
  local_failed |= donor_ranks == NULL;
  size_t donor_index = 0;
  for (size_t item = 0; !local_failed && item < unique_count; ++item) {
    if (item == 0 || unique_ranks[item] != unique_ranks[item - 1]) {
      donor_ranks[donor_index++] = unique_ranks[item];
    }
  }
  local_failed |= donor_index != donor_rank_count_size;
  if (sparse_consensus(world_comm, local_failed) != 0) {
    goto fail;
  }

  {
    const int source = world_rank;
    const int degree = donor_rank_count;
    const int dummy_rank = world_rank;
    const int *destinations =
        degree > 0 ? donor_ranks : &dummy_rank;
    const int mpi_status = MPI_Dist_graph_create(
        world_comm, 1, &source, &degree, destinations, MPI_UNWEIGHTED,
        MPI_INFO_NULL, 0, &request_comm);
    if (sparse_collective_mpi_or_abort(world_comm, mpi_status) != 0) {
      goto fail;
    }
  }
  request_comm_ready = 1;

  {
    int weighted = 0;
    const int mpi_status = MPI_Dist_graph_neighbors_count(
        request_comm, &request_indegree, &request_outdegree, &weighted);
    if (sparse_collective_mpi_or_abort(world_comm, mpi_status) != 0) {
      goto fail;
    }
    local_failed |= request_indegree < 0 || request_outdegree < 0 ||
                    request_outdegree != donor_rank_count;
  }
  request_sources =
      (int *)sparse_calloc((size_t)request_indegree, sizeof(*request_sources));
  request_destinations = (int *)sparse_calloc(
      (size_t)request_outdegree, sizeof(*request_destinations));
  request_send_counts = (int *)sparse_calloc(
      (size_t)request_outdegree, sizeof(*request_send_counts));
  request_send_displacements = (int *)sparse_calloc(
      (size_t)request_outdegree, sizeof(*request_send_displacements));
  request_receive_counts = (int *)sparse_calloc(
      (size_t)request_indegree, sizeof(*request_receive_counts));
  request_receive_displacements = (int *)sparse_calloc(
      (size_t)request_indegree, sizeof(*request_receive_displacements));
  local_failed |= request_sources == NULL || request_destinations == NULL ||
                  request_send_counts == NULL ||
                  request_send_displacements == NULL ||
                  request_receive_counts == NULL ||
                  request_receive_displacements == NULL;
  if (!local_failed) {
    int dummy_rank = world_rank;
    const int mpi_status = MPI_Dist_graph_neighbors(
        request_comm, request_indegree,
        request_indegree > 0 ? request_sources : &dummy_rank,
        MPI_UNWEIGHTED, request_outdegree,
        request_outdegree > 0 ? request_destinations : &dummy_rank,
        MPI_UNWEIGHTED);
    if (sparse_collective_mpi_or_abort(world_comm, mpi_status) != 0) {
      goto fail;
    }
  }
  for (int item = 0; !local_failed && item < request_outdegree; ++item) {
    const int donor = request_destinations[item];
    const size_t begin =
        sparse_unique_rank_begin(unique_ranks, unique_count, donor);
    const size_t end = sparse_unique_rank_end(unique_ranks, unique_count, donor);
    local_failed |= begin == end || end - begin > (size_t)INT_MAX ||
                    sparse_rank_position(donor_ranks, donor_rank_count,
                                         donor) < 0;
    if (!local_failed) {
      request_send_counts[item] = (int)(end - begin);
    }
  }
  size_t request_send_total = 0;
  if (!local_failed) {
    local_failed = sparse_accumulate_int_displacements(
                       request_send_counts, request_outdegree,
                       request_send_displacements, &request_send_total) != 0 ||
                   request_send_total != unique_count;
  }
  request_send_cells = (uint64_t *)sparse_calloc(
      request_send_total, sizeof(*request_send_cells));
  local_failed |= request_send_cells == NULL;
  for (int item = 0; !local_failed && item < request_outdegree; ++item) {
    const int donor = request_destinations[item];
    const size_t begin =
        sparse_unique_rank_begin(unique_ranks, unique_count, donor);
    for (int offset = 0; offset < request_send_counts[item]; ++offset) {
      request_send_cells[(size_t)request_send_displacements[item] +
                         (size_t)offset] = unique_cells[begin + (size_t)offset];
    }
  }
  if (sparse_consensus(world_comm, local_failed) != 0) {
    goto fail;
  }

  if (sparse_collective_mpi_or_abort(
          world_comm,
          MPI_Neighbor_alltoall(request_send_counts, 1, MPI_INT,
                                request_receive_counts, 1, MPI_INT,
                                request_comm)) != 0) {
    goto fail;
  }
  local_failed = sparse_accumulate_int_displacements(
                     request_receive_counts, request_indegree,
                     request_receive_displacements,
                     &request_receive_total) != 0;
  request_receive_cells = (uint64_t *)sparse_calloc(
      request_receive_total, sizeof(*request_receive_cells));
  local_failed |= request_receive_cells == NULL;
  if (sparse_consensus(world_comm, local_failed) != 0) {
    goto fail;
  }

  if (sparse_collective_mpi_or_abort(
          world_comm,
          MPI_Neighbor_alltoallv(
              request_send_cells, request_send_counts,
              request_send_displacements, MPI_UINT64_T,
              request_receive_cells, request_receive_counts,
              request_receive_displacements, MPI_UINT64_T,
              request_comm)) != 0) {
    goto fail;
  }
  for (int neighbor = 0; neighbor < request_indegree; ++neighbor) {
    uint64_t previous = 0;
    for (int item = 0; item < request_receive_counts[neighbor]; ++item) {
      const uint64_t cell = request_receive_cells
          [(size_t)request_receive_displacements[neighbor] + (size_t)item];
      local_failed |= cell >= (uint64_t)local_active_cells ||
                      (item > 0 && cell <= previous);
      previous = cell;
    }
  }
  if (sparse_consensus(world_comm, local_failed) != 0) {
    goto fail;
  }

  {
    const int dummy_rank = world_rank;
    const int *sources = request_outdegree > 0 ? request_destinations
                                                : &dummy_rank;
    const int *destinations = request_indegree > 0 ? request_sources
                                                    : &dummy_rank;
    const int mpi_status = MPI_Dist_graph_create_adjacent(
        world_comm, request_outdegree, sources, MPI_UNWEIGHTED,
        request_indegree, destinations, MPI_UNWEIGHTED, MPI_INFO_NULL, 0,
        &value_comm);
    if (sparse_collective_mpi_or_abort(world_comm, mpi_status) != 0) {
      goto fail;
    }
  }
  value_comm_ready = 1;

  {
    int weighted = 0;
    const int mpi_status = MPI_Dist_graph_neighbors_count(
        value_comm, &value_indegree, &value_outdegree, &weighted);
    if (sparse_collective_mpi_or_abort(world_comm, mpi_status) != 0) {
      goto fail;
    }
    local_failed |= value_indegree != request_outdegree ||
                    value_outdegree != request_indegree;
  }
  value_sources =
      (int *)sparse_calloc((size_t)value_indegree, sizeof(*value_sources));
  value_destinations = (int *)sparse_calloc((size_t)value_outdegree,
                                             sizeof(*value_destinations));
  plan->send_counts =
      (int *)sparse_calloc((size_t)value_outdegree, sizeof(*plan->send_counts));
  plan->send_displacements = (int *)sparse_calloc(
      (size_t)value_outdegree, sizeof(*plan->send_displacements));
  plan->receive_counts = (int *)sparse_calloc(
      (size_t)value_indegree, sizeof(*plan->receive_counts));
  plan->receive_displacements = (int *)sparse_calloc(
      (size_t)value_indegree, sizeof(*plan->receive_displacements));
  local_failed |= value_sources == NULL || value_destinations == NULL ||
                  plan->send_counts == NULL ||
                  plan->send_displacements == NULL ||
                  plan->receive_counts == NULL ||
                  plan->receive_displacements == NULL;
  if (!local_failed) {
    int dummy_rank = world_rank;
    const int mpi_status = MPI_Dist_graph_neighbors(
        value_comm, value_indegree,
        value_indegree > 0 ? value_sources : &dummy_rank,
        MPI_UNWEIGHTED, value_outdegree,
        value_outdegree > 0 ? value_destinations : &dummy_rank,
        MPI_UNWEIGHTED);
    if (sparse_collective_mpi_or_abort(world_comm, mpi_status) != 0) {
      goto fail;
    }
  }

  plan->send_cell_count = request_receive_total;
  plan->receive_cell_count = unique_count;
  plan->send_cells =
      (size_t *)sparse_calloc(plan->send_cell_count, sizeof(*plan->send_cells));
  plan->reference_slot = (size_t *)sparse_calloc(reference_count,
                                                  sizeof(*plan->reference_slot));
  local_failed |= plan->send_cells == NULL || plan->reference_slot == NULL;

  size_t send_cells_used = 0;
  for (int neighbor = 0; !local_failed && neighbor < value_outdegree;
       ++neighbor) {
    const int request_index = sparse_rank_position(
        request_sources, request_indegree, value_destinations[neighbor]);
    local_failed |= request_index < 0;
    if (local_failed) {
      break;
    }
    const size_t cell_count = (size_t)request_receive_counts[request_index];
    local_failed |= cell_count > (size_t)INT_MAX / (size_t)channel_count;
    for (size_t item = 0; !local_failed && item < cell_count; ++item) {
      const uint64_t cell = request_receive_cells
          [(size_t)request_receive_displacements[request_index] + item];
#if SIZE_MAX < UINT64_MAX
      local_failed |= cell > (uint64_t)SIZE_MAX;
#endif
      if (!local_failed) {
        plan->send_cells[send_cells_used + item] = (size_t)cell;
      }
    }
    if (!local_failed) {
      plan->send_counts[neighbor] = (int)(cell_count *
                                                  (size_t)channel_count);
      send_cells_used += cell_count;
    }
  }
  local_failed |= send_cells_used != plan->send_cell_count;

  size_t receive_cells_used = 0;
  for (int neighbor = 0; !local_failed && neighbor < value_indegree;
       ++neighbor) {
    const int donor = value_sources[neighbor];
    const int request_index = sparse_rank_position(
        request_destinations, request_outdegree, donor);
    const size_t begin =
        sparse_unique_rank_begin(unique_ranks, unique_count, donor);
    const size_t end = sparse_unique_rank_end(unique_ranks, unique_count, donor);
    local_failed |= request_index < 0 || begin == end ||
                    (size_t)request_send_counts[request_index] != end - begin ||
                    end - begin > (size_t)INT_MAX / (size_t)channel_count;
    if (local_failed) {
      break;
    }
    plan->receive_counts[neighbor] =
        (int)((end - begin) * (size_t)channel_count);
    for (size_t item = begin; item < end; ++item) {
      unique_final_slot[item] = receive_cells_used + item - begin;
    }
    receive_cells_used += end - begin;
  }
  local_failed |= receive_cells_used != plan->receive_cell_count;
  for (size_t item = 0; !local_failed && item < reference_count; ++item) {
    const size_t unique = reference_unique[item];
    local_failed |= unique >= unique_count ||
                    unique_final_slot[unique] >= plan->receive_cell_count;
    if (!local_failed) {
      plan->reference_slot[item] = unique_final_slot[unique];
    }
  }

  size_t send_value_count = 0;
  size_t receive_value_count = 0;
  if (!local_failed) {
    local_failed = sparse_accumulate_int_displacements(
                       plan->send_counts, value_outdegree,
                       plan->send_displacements, &send_value_count) != 0 ||
                   sparse_accumulate_int_displacements(
                       plan->receive_counts, value_indegree,
                       plan->receive_displacements,
                       &receive_value_count) != 0;
  }
  size_t expected_send_values = 0;
  size_t expected_receive_values = 0;
  if (!local_failed) {
    local_failed = sparse_multiply_size(plan->send_cell_count,
                                        (size_t)channel_count,
                                        &expected_send_values) != 0 ||
                   sparse_multiply_size(plan->receive_cell_count,
                                        (size_t)channel_count,
                                        &expected_receive_values) != 0 ||
                   send_value_count != expected_send_values ||
                   receive_value_count != expected_receive_values;
  }
  if (!local_failed) {
    local_failed = send_value_count > UINT64_MAX / sizeof(double) ||
                   receive_value_count > UINT64_MAX / sizeof(double);
  }
  plan->send_values =
      (double *)sparse_calloc(send_value_count, sizeof(*plan->send_values));
  plan->receive_values = (double *)sparse_calloc(
      receive_value_count, sizeof(*plan->receive_values));
  local_failed |= plan->send_values == NULL || plan->receive_values == NULL;

  if (sparse_consensus(world_comm, local_failed) != 0) {
    goto fail;
  }

  if (sparse_comm_free_or_abort(world_comm, &request_comm) != 0) {
    goto fail;
  }
  request_comm_ready = 0;
  plan->world_comm = world_comm;
  plan->value_comm = value_comm;
  value_comm = MPI_COMM_NULL;
  value_comm_ready = 0;
  plan->ready = 1;
  plan->channel_count = channel_count;
  plan->indegree = value_indegree;
  plan->outdegree = value_outdegree;
  plan->local_active_cells = local_active_cells;
  plan->reference_count = reference_count;
  plan->send_bytes_per_call =
      (uint64_t)send_value_count * (uint64_t)sizeof(double);
  plan->receive_bytes_per_call =
      (uint64_t)receive_value_count * (uint64_t)sizeof(double);

  free(records);
  free(unique_ranks);
  free(unique_cells);
  free(reference_unique);
  free(unique_final_slot);
  free(donor_ranks);
  free(request_sources);
  free(request_destinations);
  free(request_send_counts);
  free(request_send_displacements);
  free(request_receive_counts);
  free(request_receive_displacements);
  free(request_send_cells);
  free(request_receive_cells);
  free(value_sources);
  free(value_destinations);
  return 0;

fail:
  if (value_comm_ready) {
    sparse_comm_free_or_abort(world_comm, &value_comm);
    value_comm_ready = 0;
  }
  if (request_comm_ready) {
    sparse_comm_free_or_abort(world_comm, &request_comm);
    request_comm_ready = 0;
  }
  free(records);
  free(unique_ranks);
  free(unique_cells);
  free(reference_unique);
  free(unique_final_slot);
  free(donor_ranks);
  free(request_sources);
  free(request_destinations);
  free(request_send_counts);
  free(request_send_displacements);
  free(request_receive_counts);
  free(request_receive_displacements);
  free(request_send_cells);
  free(request_receive_cells);
  free(value_sources);
  free(value_destinations);
  sparse_release_arrays(plan);
  plan->world_comm = MPI_COMM_NULL;
  plan->value_comm = MPI_COMM_NULL;
  plan->ready = 0;
  return -1;
}

int gamera_no_sparse_plan_exchange(gamera_no_sparse_plan *plan,
                                   gamera_no_sparse_pack_fn pack,
                                   void *context) {
  /* Plan validity and all byte/count bounds are collective init invariants. */
  if (plan == NULL || !plan->ready || pack == NULL ||
      plan->world_comm == MPI_COMM_NULL || plan->value_comm == MPI_COMM_NULL) {
    const MPI_Comm world_comm =
        plan != NULL ? plan->world_comm : MPI_COMM_WORLD;
    return sparse_programmer_error_or_abort(world_comm);
  }

  const double pack_start = MPI_Wtime();
  const int local_failed =
      pack(plan->send_cells, plan->send_cell_count, plan->channel_count,
           plan->send_values, context) != 0;
  plan->last_pack_seconds = MPI_Wtime() - pack_start;
  if (sparse_consensus(plan->world_comm, local_failed) != 0) {
    return -1;
  }

  const double collective_start = MPI_Wtime();
  const int mpi_status = MPI_Neighbor_alltoallv(
      plan->send_values, plan->send_counts, plan->send_displacements,
      MPI_DOUBLE, plan->receive_values, plan->receive_counts,
      plan->receive_displacements, MPI_DOUBLE, plan->value_comm);
  plan->last_collective_seconds = MPI_Wtime() - collective_start;
  if (sparse_collective_mpi_or_abort(plan->world_comm, mpi_status) != 0) {
    return -1;
  }

  if (plan->call_count == UINT64_MAX) {
    plan->statistics_saturated = 1;
  } else {
    ++plan->call_count;
  }
  sparse_accumulate_u64(plan->send_bytes_per_call,
                        &plan->useful_send_bytes,
                        &plan->statistics_saturated);
  sparse_accumulate_u64(plan->receive_bytes_per_call,
                        &plan->useful_receive_bytes,
                        &plan->statistics_saturated);
  return 0;
}

void gamera_no_sparse_plan_destroy(gamera_no_sparse_plan *plan) {
  if (plan == NULL) {
    return;
  }
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  if (initialized) {
    MPI_Finalized(&finalized);
  }
  if (initialized && !finalized && plan->value_comm != MPI_COMM_NULL) {
    sparse_comm_free_or_abort(plan->world_comm, &plan->value_comm);
  }
  sparse_release_arrays(plan);
  memset(plan, 0, sizeof(*plan));
  plan->world_comm = MPI_COMM_NULL;
  plan->value_comm = MPI_COMM_NULL;
}

#else

/* Keep the globbed translation unit valid in default-OFF builds. */
typedef int gamera_no_sparse_plan_disabled_translation_unit;

#endif
