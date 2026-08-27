#ifndef NONORTHOGONAL_YINYANG_SPARSE_PLAN_H
#define NONORTHOGONAL_YINYANG_SPARSE_PLAN_H

#include <mpi.h>
#include <stddef.h>
#include <stdint.h>

typedef struct {
  int donor_rank;
  size_t donor_cell;
} gamera_no_sparse_reference;

typedef int (*gamera_no_sparse_pack_fn)(const size_t *local_cells,
                                        size_t cell_count,
                                        int channel_count,
                                        double *cell_major_values,
                                        void *context);

typedef struct {
  MPI_Comm world_comm;
  MPI_Comm value_comm;
  int ready;
  int channel_count;
  int indegree;
  int outdegree;
  size_t local_active_cells;
  size_t reference_count;
  size_t receive_cell_count;
  size_t send_cell_count;
  size_t *reference_slot;
  size_t *send_cells;
  double *send_values;
  double *receive_values;
  int *send_counts;
  int *send_displacements;
  int *receive_counts;
  int *receive_displacements;
  uint64_t call_count;
  uint64_t send_bytes_per_call;
  uint64_t receive_bytes_per_call;
  uint64_t useful_send_bytes;
  uint64_t useful_receive_bytes;
  int statistics_saturated;
  double last_pack_seconds;
  double last_collective_seconds;
} gamera_no_sparse_plan;

#define GAMERA_NO_SPARSE_PLAN_INITIALIZER \
  { .world_comm = MPI_COMM_NULL, .value_comm = MPI_COMM_NULL }

/*
 * Build a directed donor-value plan collectively on world_comm.  Every rank
 * must call build/exchange/destroy in the same lifetime order, starting from
 * a pristine GAMERA_NO_SPARSE_PLAN_INITIALIZER object.  Rebuilding a live
 * plan is a fatal collective-lifetime programmer error.  Every donor rank in
 * references must lie in [allowed_rank_begin, allowed_rank_end), and
 * reference_slot retains the input reference order.
 */
int gamera_no_sparse_plan_build(
    gamera_no_sparse_plan *plan, MPI_Comm world_comm,
    const gamera_no_sparse_reference *references, size_t reference_count,
    size_t local_active_cells, int channel_count, int allowed_rank_begin,
    int allowed_rank_end);

/*
 * Pack only requested cells, then perform one reverse-graph neighbor exchange.
 * A successfully built plan is a collective lifetime invariant: every rank
 * must call exchange in the same order with a non-NULL callback.  Callback
 * failures are synchronized before the neighbor collective.  MPI collective
 * failures are fatal because a portable rank-local recovery order does not
 * exist after a partially failed collective.
 */
int gamera_no_sparse_plan_exchange(gamera_no_sparse_plan *plan,
                                   gamera_no_sparse_pack_fn pack,
                                   void *context);

static inline double gamera_no_sparse_plan_value(
    const gamera_no_sparse_plan *plan, size_t reference_index, int channel) {
  return plan->receive_values[plan->reference_slot[reference_index] *
                                  (size_t)plan->channel_count +
                              (size_t)channel];
}

/* Collective when value_comm is live; call before MPI_Finalize on all ranks. */
void gamera_no_sparse_plan_destroy(gamera_no_sparse_plan *plan);

#endif
