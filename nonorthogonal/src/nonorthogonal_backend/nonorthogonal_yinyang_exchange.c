#include "nonorthogonal_yinyang_exchange.h"

#include "nonorthogonal_legacy_adapter.h"
#include "nonorthogonal_operators.h"
#include "nonorthogonal_state.h"
#include "nonorthogonal_step.h"
#include "nonorthogonal_yinyang.h"
#ifdef GAMERA_YINYANG_SPARSE_OVERSET
#include "nonorthogonal_yinyang_sparse_plan.h"
#endif

#include "config.h"
#include "log.h"
#include "problem.h"
#include "setup_mpi.h"

#include <math.h>
#include <float.h>
#include <mpi.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#if defined(GAMERA_YINYANG_HDIV_OPTIMIZED) ||                                \
    defined(GAMERA_YINYANG_HDIV_PROFILE)
#include <omp.h>
#endif

typedef struct {
  int receiver[3];
  int donor_rank[8];
  size_t donor_cell[8];
  double weight[8];
} receptor_t;

typedef struct {
  int direction;
  size_t receiver_edge;
  double blend_weight;
  int donor_rank[8];
  size_t donor_cell[8];
  double weight[8];
} edge_receptor_t;

typedef struct {
  int direction;
  size_t receiver_face;
  double blend_weight;
  int donor_rank[8];
  size_t donor_cell[8];
  double weight[8];
} face_receptor_t;

static receptor_t *receptors;
static size_t receptor_count;
#ifndef GAMERA_YINYANG_SPARSE_OVERSET
static double *send_buffer;
static double *gather_buffer;
#endif
static size_t local_active_cells;
#ifndef GAMERA_YINYANG_SPARSE_OVERSET
static int send_count;
#endif
static int exchange_ready;
static edge_receptor_t *edge_receptors;
static size_t edge_receptor_count;
#ifndef GAMERA_YINYANG_SPARSE_OVERSET
static double *emf_send_buffer;
static double *emf_gather_buffer;
static int emf_send_count;
#endif
static int emf_exchange_ready;
static face_receptor_t *face_receptors;
static size_t face_receptor_count;
static size_t *magnetic_recovery_cells;
static size_t magnetic_recovery_cell_count;
#ifndef GAMERA_YINYANG_SPARSE_OVERSET
static double *magnetic_send_buffer;
static double *magnetic_gather_buffer;
static int magnetic_send_count;
#endif
static int magnetic_exchange_ready;
#ifdef GAMERA_YINYANG_SPARSE_OVERSET
static gamera_no_sparse_plan fluid_full_plan =
    GAMERA_NO_SPARSE_PLAN_INITIALIZER;
static gamera_no_sparse_plan fluid_ghost_plan =
    GAMERA_NO_SPARSE_PLAN_INITIALIZER;
static gamera_no_sparse_plan emf_plan = GAMERA_NO_SPARSE_PLAN_INITIALIZER;
static gamera_no_sparse_plan magnetic_ghost_plan =
    GAMERA_NO_SPARSE_PLAN_INITIALIZER;
static size_t *fluid_ghost_reference_base;
#ifdef GAMERA_YINYANG_SPARSE_OVERSET_PROFILE
enum {
  GAMERA_NO_SPARSE_PROFILE_PACK = 0,
  GAMERA_NO_SPARSE_PROFILE_NEIGHBOR = 1,
  GAMERA_NO_SPARSE_PROFILE_RECONSTRUCT = 2,
  GAMERA_NO_SPARSE_PROFILE_TOTAL = 3,
  GAMERA_NO_SPARSE_PROFILE_PHASE_COUNT = 4
};
typedef struct {
  const char *kind;
  const char *phase;
  const char *mode;
  uint64_t samples;
  uint64_t send_bytes;
  uint64_t receive_bytes;
  uint64_t neighbor_collectives;
  uint64_t world_consensus_calls;
  int saturated;
  double sum[GAMERA_NO_SPARSE_PROFILE_PHASE_COUNT];
  double maximum[GAMERA_NO_SPARSE_PROFILE_PHASE_COUNT];
} sparse_profile_accumulator;
static sparse_profile_accumulator fluid_full_profile;
static sparse_profile_accumulator fluid_ghost_profile;
static sparse_profile_accumulator emf_profile;
static sparse_profile_accumulator magnetic_ghost_profile;
#endif
#endif
#ifdef GAMERA_YINYANG_HDIV_RECONCILE
static face_receptor_t *active_face_receptors;
static size_t active_face_receptor_count;
#ifdef GAMERA_YINYANG_SPARSE_OVERSET
static gamera_no_sparse_plan magnetic_active_both_plan =
    GAMERA_NO_SPARSE_PLAN_INITIALIZER;
static gamera_no_sparse_plan magnetic_active_current_plan =
    GAMERA_NO_SPARSE_PLAN_INITIALIZER;
#ifdef GAMERA_YINYANG_SPARSE_OVERSET_PROFILE
static sparse_profile_accumulator magnetic_active_both_profile;
static sparse_profile_accumulator magnetic_active_current_profile;
#endif
#endif
static double *hdiv_target_face_flux[2][GAMERA_NO_DIM];
static double *hdiv_send_buffer;
#ifndef GAMERA_YINYANG_HDIV_DISTRIBUTED
static double *hdiv_gather_buffer;
#endif
static int hdiv_send_count;
#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED
static size_t hdiv_local_face_extent[GAMERA_NO_DIM][GAMERA_NO_DIM];
static double *hdiv_local_face_base[2][GAMERA_NO_DIM];
static double *hdiv_local_face_flux[2][GAMERA_NO_DIM];
static double *hdiv_local_face_mobility[GAMERA_NO_DIM];
static size_t hdiv_local_cell_extent[GAMERA_NO_DIM];
static double *hdiv_face_halo_send;
static double *hdiv_face_halo_receive;
static size_t hdiv_face_halo_capacity;
static double *hdiv_cell_halo_send;
static double *hdiv_cell_halo_receive;
static size_t hdiv_cell_halo_capacity;
static int hdiv_distributed_layout_verified;
#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED_VERIFY
static int hdiv_distributed_verify_ready;
static int hdiv_distributed_verify_sentinel_ready;
static int hdiv_distributed_verify_old_audit_done;
static unsigned int hdiv_distributed_verify_projection_mask;
static double hdiv_distributed_verify_max_face[2];
static double hdiv_distributed_verify_max_before[2];
static double hdiv_distributed_verify_tolerance[2];
static double *hdiv_distributed_verify_old_face[GAMERA_NO_DIM];
#endif
#else
static size_t hdiv_global_face_extent[GAMERA_NO_DIM][GAMERA_NO_DIM];
static double *hdiv_global_face_base[2][GAMERA_NO_DIM];
static double *hdiv_global_face_flux[2][GAMERA_NO_DIM];
static double *hdiv_global_face_mobility[GAMERA_NO_DIM];
#endif
static double *hdiv_lambda;
static double *hdiv_residual;
static double *hdiv_preconditioned;
static double *hdiv_search;
static double *hdiv_operator_search;
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
static double *hdiv_cell_diagonal_cache;
static int hdiv_cell_diagonal_cache_ready;
static double *hdiv_dot_block_partial;
static size_t hdiv_dot_block_count;
#endif
#ifdef GAMERA_YINYANG_HDIV_PROFILE
static unsigned long hdiv_profile_calls;
enum {
  GAMERA_NO_HDIV_PROFILE_TOTAL = 0,
  GAMERA_NO_HDIV_PROFILE_PREPARE,
  GAMERA_NO_HDIV_PROFILE_GATHER,
  GAMERA_NO_HDIV_PROFILE_ASSEMBLE,
  GAMERA_NO_HDIV_PROFILE_PROJECT,
  GAMERA_NO_HDIV_PROFILE_RECOVER,
  GAMERA_NO_HDIV_PROFILE_PHASE_COUNT
};
static double hdiv_profile_sum[GAMERA_NO_HDIV_PROFILE_PHASE_COUNT];
static double hdiv_profile_maximum[GAMERA_NO_HDIV_PROFILE_PHASE_COUNT];
#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED
enum {
  GAMERA_NO_HDIV_DETAIL_CELL_HALO = 0,
  GAMERA_NO_HDIV_DETAIL_FACE_HALO,
  GAMERA_NO_HDIV_DETAIL_REDUCTION,
  GAMERA_NO_HDIV_DETAIL_OPERATOR,
  GAMERA_NO_HDIV_DETAIL_VECTOR,
  GAMERA_NO_HDIV_DETAIL_CORRECTION,
  GAMERA_NO_HDIV_DETAIL_COUNT
};
static int hdiv_profile_active;
static double hdiv_profile_detail_current[GAMERA_NO_HDIV_DETAIL_COUNT];
static unsigned long hdiv_profile_patch_allreduce_current;
static unsigned long hdiv_profile_world_allreduce_current;
static unsigned long hdiv_profile_math_allreduce_current;
static unsigned long hdiv_profile_cell_sendrecv_current;
static unsigned long hdiv_profile_face_sendrecv_current;
static unsigned long hdiv_profile_layout_sendrecv_current;
static unsigned long long hdiv_profile_cell_send_bytes_current;
static unsigned long long hdiv_profile_face_send_bytes_current;
static unsigned long hdiv_profile_pcg_iterations_current[2];
#endif
#endif
static int hdiv_exchange_ready;
static int hdiv_history_ready;
static double hdiv_last_time = NAN;
#endif

enum { GAMERA_NO_YINYANG_TIME_LEVELS = 2 };

int gamera_no_yinyang_receptor_count;
int gamera_no_yinyang_active_receptor_count;
double gamera_no_yinyang_max_donor_extrapolation;
int gamera_no_yinyang_edge_receptor_count;
double gamera_no_yinyang_max_emf_donor_extrapolation;
int gamera_no_yinyang_magnetic_face_receptor_count;
double gamera_no_yinyang_max_magnetic_donor_extrapolation;
int gamera_no_yinyang_active_magnetic_face_receptor_count;
double gamera_no_yinyang_hdiv_max_before;
double gamera_no_yinyang_hdiv_max_after;
double gamera_no_yinyang_hdiv_max_correction;
int gamera_no_yinyang_hdiv_max_iterations;
double gamera_no_yinyang_hdiv_min_weight = DBL_MAX;
double gamera_no_yinyang_hdiv_max_weight;
double gamera_no_yinyang_hdiv_min_metric_cosine = 1.0;

static void yinyang_exchange_destroy_local(void);

static gamera_no_vec3 vector_cross(gamera_no_vec3 left,
                                   gamera_no_vec3 right) {
  return (gamera_no_vec3){{
      left.value[1] * right.value[2] - left.value[2] * right.value[1],
      left.value[2] * right.value[0] - left.value[0] * right.value[2],
      left.value[0] * right.value[1] - left.value[1] * right.value[0]}};
}

static double vector_dot(gamera_no_vec3 left, gamera_no_vec3 right) {
  double result = 0.0;
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    result += left.value[component] * right.value[component];
  }
  return result;
}

static size_t local_active_index(int i, int j, int k) {
  return ((size_t)i * (size_t)config.nj + (size_t)j) *
             (size_t)config.nk +
         (size_t)k;
}

/*
 * A rank-local failure must never let one process return while its peers enter
 * the next world collective.  Besides hanging an error path, that can obscure
 * the original allocation/geometry failure behind an MPI timeout.  Call this
 * gate between every local-only phase and the following collective phase.
 */
static int world_consensus_failure(int local_failed, const char *stage) {
  int any_failed = 0;
#if defined(GAMERA_YINYANG_HDIV_PROFILE) &&                                \
    defined(GAMERA_YINYANG_HDIV_DISTRIBUTED)
  const double reduction_start =
      hdiv_profile_active ? omp_get_wtime() : 0.0;
#endif
  const int mpi_status =
      MPI_Allreduce(&local_failed, &any_failed, 1, MPI_INT, MPI_MAX,
                    MPI_COMM_WORLD);
#if defined(GAMERA_YINYANG_HDIV_PROFILE) &&                                \
    defined(GAMERA_YINYANG_HDIV_DISTRIBUTED)
  if (hdiv_profile_active) {
    hdiv_profile_detail_current[GAMERA_NO_HDIV_DETAIL_REDUCTION] +=
        omp_get_wtime() - reduction_start;
    ++hdiv_profile_world_allreduce_current;
  }
#endif
  if (mpi_status != MPI_SUCCESS) {
    /* A partially failed collective has no portable rank-local recovery or
     * cleanup order.  With MPI_ERRORS_RETURN, make the failure terminal just
     * as the sparse-plan graph collectives do. */
    MPI_Abort(MPI_COMM_WORLD, mpi_status);
    return -1;
  }
  if (any_failed) {
    if (rank == 0) {
      log_error("Yin-Yang exchange failed collectively during %s", stage);
    }
    return -1;
  }
  return 0;
}

static int owner_rank_and_cell(int donor_patch, const int global[3],
                               int *world_rank, size_t *cell) {
  int coords[3];
  int local[3];
  const int local_shape[3] = {config.ni, config.nj, config.nk};
  const int global_shape[3] = {config.ni_global, config.nj_global,
                               config.nk_global};
  for (int axis = 0; axis < 3; ++axis) {
    if (global[axis] < 0 || global[axis] >= global_shape[axis]) {
      return -1;
    }
    coords[axis] = global[axis] / local_shape[axis];
    local[axis] = global[axis] % local_shape[axis];
  }
  int donor_patch_rank = -1;
  if (MPI_Cart_rank(comm_cart, coords, &donor_patch_rank) != MPI_SUCCESS) {
    return -1;
  }
  *world_rank = donor_patch * patch_size + donor_patch_rank;
  *cell = local_active_index(local[0], local[1], local[2]);
  return 0;
}

static void interpolation_logical(double logical, int count, int *base,
                                  double *fraction, double *outside) {
  *outside = 0.0;
  if (logical < 0.0) {
    *outside = -logical;
    *base = 0;
    *fraction = 0.0;
  } else if (logical > (double)(count - 1)) {
    *outside = logical - (double)(count - 1);
    *base = count - 2;
    *fraction = 1.0;
  } else {
    *base = (int)floor(logical);
    if (*base >= count - 1) {
      *base = count - 2;
      *fraction = 1.0;
    } else {
      *fraction = logical - (double)*base;
    }
  }
}

static int interpolation_axis(int axis, double coordinate, double lower,
                              double upper, int count, int *base,
                              double *fraction, double *outside) {
  double logical;
#ifdef GAMERA_NONORTHOGONAL_HAS_RADIAL_MAP
  if (axis == GAMERA_NO_I) {
    if (problem_nonorthogonal_radial_logical(coordinate, &logical) != 0) {
      return -1;
    }
  } else
#endif
  {
    const double spacing = (upper - lower) / (double)count;
    logical = (coordinate - lower) / spacing - 0.5;
  }
  interpolation_logical(logical, count, base, fraction, outside);
  return 0;
}

static int append_receptor(size_t *used, size_t capacity, int i, int j,
                           int k, gamera_no_vec3 point) {
  if (*used >= capacity) {
    return -1;
  }
  const int donor_patch = 1 - patch_id;
  double logical[3];
  if (gamera_no_yinyang_global_to_logical(
          donor_patch, point, &logical[0], &logical[1], &logical[2]) != 0) {
    return -1;
  }
  const double lower[3] = {x1min_global, x2min_global, x3min_global};
  const double upper[3] = {x1max_global, x2max_global, x3max_global};
  const int global_shape[3] = {config.ni_global, config.nj_global,
                               config.nk_global};
  int base[3];
  double fraction[3];
  double receptor_extrapolation = 0.0;
  int invalid_extrapolation = 0;
  for (int axis = 0; axis < 3; ++axis) {
    double outside;
    if (interpolation_axis(axis, logical[axis], lower[axis], upper[axis],
                           global_shape[axis], &base[axis], &fraction[axis],
                           &outside) != 0) {
      return -1;
    }
    gamera_no_yinyang_max_donor_extrapolation =
        fmax(gamera_no_yinyang_max_donor_extrapolation, outside);
    receptor_extrapolation = fmax(receptor_extrapolation, outside);
    double allowed = 1.0e-2;
#ifdef GAMERA_NONORTHOGONAL_HAS_RADIAL_MAP
    /* A finite-volume centroid is slightly exterior to the parametric
     * midpoint in the widest stretched cells.  The paired patch has the same
     * shell; a small one-sided radial clamp is therefore geometrical, while
     * angular extrapolation remains subject to the original strict limit. */
    if (axis == GAMERA_NO_I) {
      allowed = 5.0e-2;
    }
#endif
    invalid_extrapolation |= outside > allowed;
  }
  /*
   * Curved-cell volume centroids from the two patches need not map to exactly
   * the same radial cell-center locus.  Permit only a tiny sub-cell clamp;
   * geometrical coverage failures remain fatal.
   */
  if (invalid_extrapolation) {
    return -1;
  }

  receptor_t *receptor = &receptors[*used];
  receptor->receiver[0] = i;
  receptor->receiver[1] = j;
  receptor->receiver[2] = k;
  int corner = 0;
  for (int di = 0; di <= 1; ++di) {
    for (int dj = 0; dj <= 1; ++dj) {
      for (int dk = 0; dk <= 1; ++dk, ++corner) {
        const int donor_global[3] = {base[0] + di, base[1] + dj,
                                     base[2] + dk};
        if (owner_rank_and_cell(donor_patch, donor_global,
                                &receptor->donor_rank[corner],
                                &receptor->donor_cell[corner]) != 0) {
          return -1;
        }
        const double wi = di == 0 ? 1.0 - fraction[0] : fraction[0];
        const double wj = dj == 0 ? 1.0 - fraction[1] : fraction[1];
        const double wk = dk == 0 ? 1.0 - fraction[2] : fraction[2];
        receptor->weight[corner] = wi * wj * wk;
      }
    }
  }
  ++*used;
  return 0;
}

static int append_edge_receptor(size_t *used, size_t capacity, int direction,
                                size_t receiver_edge, double blend_weight,
                                gamera_no_vec3 point) {
  if (*used >= capacity) {
    return -1;
  }
  const int donor_patch = 1 - patch_id;
  double logical[3];
  if (gamera_no_yinyang_global_to_logical(
          donor_patch, point, &logical[0], &logical[1], &logical[2]) != 0) {
    return -1;
  }
  const double lower[3] = {x1min_global, x2min_global, x3min_global};
  const double upper[3] = {x1max_global, x2max_global, x3max_global};
  const int global_shape[3] = {config.ni_global, config.nj_global,
                               config.nk_global};
  int base[3];
  double fraction[3];
  double receptor_extrapolation = 0.0;
  for (int axis = 0; axis < 3; ++axis) {
    double outside;
    if (interpolation_axis(axis, logical[axis], lower[axis], upper[axis],
                           global_shape[axis], &base[axis], &fraction[axis],
                           &outside) != 0) {
      return -1;
    }
    receptor_extrapolation = fmax(receptor_extrapolation, outside);
    gamera_no_yinyang_max_emf_donor_extrapolation =
        fmax(gamera_no_yinyang_max_emf_donor_extrapolation, outside);
  }
  /*
   * CT edges lie on cell faces in the two transverse directions.  A global
   * domain edge is therefore exactly one half cell outside the donor
   * cell-center lattice; a one-sided clamp of at most 0.5 cell is expected.
   */
  if (receptor_extrapolation > 0.525) {
    return -1;
  }

  edge_receptor_t *receptor = &edge_receptors[*used];
  receptor->direction = direction;
  receptor->receiver_edge = receiver_edge;
  receptor->blend_weight = blend_weight;
  int corner = 0;
  for (int di = 0; di <= 1; ++di) {
    for (int dj = 0; dj <= 1; ++dj) {
      for (int dk = 0; dk <= 1; ++dk, ++corner) {
        const int donor_global[3] = {base[0] + di, base[1] + dj,
                                     base[2] + dk};
        if (owner_rank_and_cell(donor_patch, donor_global,
                                &receptor->donor_rank[corner],
                                &receptor->donor_cell[corner]) != 0) {
          return -1;
        }
        const double wi = di == 0 ? 1.0 - fraction[0] : fraction[0];
        const double wj = dj == 0 ? 1.0 - fraction[1] : fraction[1];
        const double wk = dk == 0 ? 1.0 - fraction[2] : fraction[2];
        receptor->weight[corner] = wi * wj * wk;
      }
    }
  }
  ++*used;
  return 0;
}

static int append_face_receptor_to(face_receptor_t *list, size_t *used,
                                   size_t capacity, int direction,
                                   size_t receiver_face, double blend_weight,
                                   gamera_no_vec3 point) {
  if (*used >= capacity) {
    return -1;
  }
  const int donor_patch = 1 - patch_id;
  double logical[3];
  if (gamera_no_yinyang_global_to_logical(
          donor_patch, point, &logical[0], &logical[1], &logical[2]) != 0) {
    return -1;
  }
  const double lower[3] = {x1min_global, x2min_global, x3min_global};
  const double upper[3] = {x1max_global, x2max_global, x3max_global};
  const int global_shape[3] = {config.ni_global, config.nj_global,
                               config.nk_global};
  int base[3];
  double fraction[3];
  double receptor_extrapolation = 0.0;
  for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
    double outside;
    if (interpolation_axis(axis, logical[axis], lower[axis], upper[axis],
                           global_shape[axis], &base[axis], &fraction[axis],
                           &outside) != 0) {
      return -1;
    }
    receptor_extrapolation = fmax(receptor_extrapolation, outside);
    gamera_no_yinyang_max_magnetic_donor_extrapolation =
        fmax(gamera_no_yinyang_max_magnetic_donor_extrapolation, outside);
  }
  if (receptor_extrapolation > 0.525) {
    return -1;
  }
  face_receptor_t *receptor = &list[*used];
  receptor->direction = direction;
  receptor->receiver_face = receiver_face;
  receptor->blend_weight = blend_weight;
  int corner = 0;
  for (int di = 0; di <= 1; ++di) {
    for (int dj = 0; dj <= 1; ++dj) {
      for (int dk = 0; dk <= 1; ++dk, ++corner) {
        const int donor_global[3] = {base[0] + di, base[1] + dj,
                                     base[2] + dk};
        if (owner_rank_and_cell(donor_patch, donor_global,
                                &receptor->donor_rank[corner],
                                &receptor->donor_cell[corner]) != 0) {
          return -1;
        }
        const double wi = di == 0 ? 1.0 - fraction[0] : fraction[0];
        const double wj = dj == 0 ? 1.0 - fraction[1] : fraction[1];
        const double wk = dk == 0 ? 1.0 - fraction[2] : fraction[2];
        receptor->weight[corner] = wi * wj * wk;
      }
    }
  }
  ++*used;
  return 0;
}

static gamera_no_yinyang_angular_domain current_angular_domain(void) {
  const gamera_no_yinyang_angular_domain domain = {
      x2min_global, x2max_global, x3min_global, x3max_global,
      (size_t)config.nj_global, (size_t)config.nk_global, 2.0e-2};
  return domain;
}

static int angular_patch_margin(int patch, gamera_no_vec3 point, int *valid,
                                double *margin) {
  const gamera_no_yinyang_angular_domain domain = current_angular_domain();
  return gamera_no_yinyang_angular_margin(patch, point, &domain, valid,
                                          margin);
}

int gamera_no_yinyang_physical_owner(gamera_no_vec3 point,
                                     int *owner_patch, int *in_overlap,
                                     double margin_cells[2]) {
  const gamera_no_yinyang_angular_domain domain = current_angular_domain();
  return gamera_no_yinyang_composite_owner(
      point, &domain, owner_patch, in_overlap, margin_cells);
}

/*
 * In the geometrical overlap, only the patch whose point is farther from its
 * logical angular boundary is an owner.  The other patch is an active
 * receptor.  Updating these receptor cells after every full step prevents two
 * independent shock histories from developing in the redundant overlap.
 */
#ifndef GAMERA_YINYANG_MFE_INTERFACE
static int active_cell_is_receptor(gamera_no_vec3 point, int *result) {
  int owner;
  int overlap;
  if (result == NULL ||
      gamera_no_yinyang_physical_owner(point, &owner, &overlap, NULL) != 0) {
    return -1;
  }
  *result = overlap && owner != patch_id;
  return 0;
}
#endif

/*
 * Blend both patches' Cartesian electric fields through the geometrical
 * overlap.  At equal margins both use the same 50/50 state; the transfer
 * approaches the donor at the receiver fringe and the local value deep in
 * the owner.  This partition-of-unity avoids a curl-producing hard switch.
 */
static int active_edge_transfer(gamera_no_vec3 point, int *result,
                                double *weight) {
  const int donor_patch = 1 - patch_id;
  int self_valid;
  int donor_valid;
  double self_margin;
  double donor_margin;
  if (result == NULL || weight == NULL ||
      angular_patch_margin(patch_id, point, &self_valid, &self_margin) != 0 ||
      angular_patch_margin(donor_patch, point, &donor_valid, &donor_margin) !=
          0 ||
      !self_valid) {
    return -1;
  }
  if (!donor_valid) {
    *result = 0;
    *weight = 0.0;
    return 0;
  }
  const double blend_width_cells = 4.0;
  *weight = 0.5 *
            (1.0 + tanh((donor_margin - self_margin) / blend_width_cells));
  *result = *weight > 1.0e-6;
  return 0;
}

#ifdef GAMERA_YINYANG_HDIV_RECONCILE
/*
 * Form a smooth donor target throughout the geometrical overlap.  The target
 * follows the same complementary tanh partition used for edge EMF exchange;
 * a wide owner-side cutoff only avoids negligible far-tail transfers.  The
 * subsequent H(div) projection chooses its mobility independently.
 */
static int active_magnetic_face_transfer(gamera_no_vec3 point, int *target,
                                         double *blend_weight) {
  const int donor_patch = 1 - patch_id;
  int self_valid;
  int donor_valid;
  double self_margin;
  double donor_margin;
  if (target == NULL || blend_weight == NULL ||
      angular_patch_margin(patch_id, point, &self_valid, &self_margin) != 0 ||
      angular_patch_margin(donor_patch, point, &donor_valid, &donor_margin) !=
          0 ||
      !self_valid) {
    return -1;
  }
  if (!donor_valid) {
    *target = 0;
    *blend_weight = 0.0;
    return 0;
  }
  *blend_weight = 0.5 *
                  (1.0 + tanh((donor_margin - self_margin) / 4.0));
  if (self_margin - donor_margin <= 16.0) {
    *target = *blend_weight > 1.0e-6;
  } else {
    *target = 0;
    *blend_weight = 0.0;
  }
  return 0;
}
#endif

static int initialize_exchange(void) {
  if (exchange_ready) {
    return 0;
  }
  if (patch_count != 2 || (patch_id != GAMERA_NO_YIN_PATCH &&
                           patch_id != GAMERA_NO_YANG_PATCH) ||
      config.ni_global < 2 || config.nj_global < 2 ||
      config.nk_global < 2) {
    return -1;
  }
  gamera_no_grid *grid = gamera_no_legacy_grid();
  if (grid == NULL) {
    return -1;
  }
  const size_t capacity = (size_t)config.ni * grid->cell_extent[1] *
                          grid->cell_extent[2];
  receptors = (receptor_t *)calloc(capacity, sizeof(*receptors));
  local_active_cells =
      (size_t)config.ni * (size_t)config.nj * (size_t)config.nk;
  if (receptors == NULL || local_active_cells == 0
#ifndef GAMERA_YINYANG_SPARSE_OVERSET
      ||
      local_active_cells >
          (size_t)INT32_MAX /
              (GAMERA_NO_YINYANG_TIME_LEVELS * GAMERA_NO_FLUX_COUNT)
#endif
  ) {
    yinyang_exchange_destroy_local();
    return -1;
  }
#ifndef GAMERA_YINYANG_SPARSE_OVERSET
  send_count =
      (int)(local_active_cells * GAMERA_NO_FLUX_COUNT *
            GAMERA_NO_YINYANG_TIME_LEVELS);
  send_buffer = (double *)malloc((size_t)send_count * sizeof(*send_buffer));
  gather_buffer =
      (double *)malloc((size_t)send_count * (size_t)size *
                       sizeof(*gather_buffer));
  if (send_buffer == NULL || gather_buffer == NULL) {
    yinyang_exchange_destroy_local();
    return -1;
  }
#endif

  receptor_count = 0;
  size_t active_receptor_count = 0;
  gamera_no_yinyang_max_donor_extrapolation = 0.0;
  for (int i = is; i <= ie; ++i) {
    for (size_t j = 0; j < grid->cell_extent[1]; ++j) {
      const int global_j =
          proc_coords[1] * config.nj + (int)j - js;
      for (size_t k = 0; k < grid->cell_extent[2]; ++k) {
        const int global_k =
            proc_coords[2] * config.nk + (int)k - ks;
        const size_t cell = gamera_no_index3(
            grid->cell_extent, (size_t)i, j, k);
        const int angular_active =
            global_j >= 0 && global_j < config.nj_global && global_k >= 0 &&
            global_k < config.nk_global;
#ifdef GAMERA_YINYANG_MFE_INTERFACE
        /*
         * MFE-style overset coupling evolves every active cell on its local
         * patch.  Only angular ghosts are donor receptors; this keeps the
         * fluid predictor on the same patch history as its CT face fluxes.
         */
        if (angular_active) {
          continue;
        }
#else
        int active_receptor = 0;
        if (angular_active) {
          if (active_cell_is_receptor(grid->cell[cell].centroid,
                                      &active_receptor) != 0) {
            log_error("Yin-Yang ownership lookup failed for patch %d local "
                      "cell (%d,%zu,%zu)",
                      patch_id, i, j, k);
            yinyang_exchange_destroy_local();
            return -1;
          }
          if (!active_receptor) {
            continue;
          }
        }
#endif
        if (append_receptor(&receptor_count, capacity, i, (int)j, (int)k,
                            grid->cell[cell].centroid) != 0) {
          log_error("Yin-Yang donor lookup failed for patch %d local "
                    "cell (%d,%zu,%zu), max extrapolation %.6e cells",
                    patch_id, i, j, k,
                    gamera_no_yinyang_max_donor_extrapolation);
          yinyang_exchange_destroy_local();
          return -1;
        }
        if (angular_active) {
          ++active_receptor_count;
        }
      }
    }
  }
  gamera_no_yinyang_receptor_count = (int)receptor_count;
  gamera_no_yinyang_active_receptor_count = (int)active_receptor_count;
  exchange_ready = 1;
  log_info("Yin-Yang patch %d prepared %zu HD receptor cells (%zu active "
           "overlap, %zu angular ghost); maximum donor extrapolation %.3e "
           "cells",
           patch_id, receptor_count, active_receptor_count,
           receptor_count - active_receptor_count,
           gamera_no_yinyang_max_donor_extrapolation);
  return 0;
}

static gamera_no_vec3 edge_midpoint(const gamera_no_grid *grid,
                                    int direction, size_t i, size_t j,
                                    size_t k) {
  size_t upper[3] = {i, j, k};
  ++upper[direction];
  const gamera_no_vec3 start =
      grid->vertex[gamera_no_index3(grid->vertex_extent, i, j, k)];
  const gamera_no_vec3 end = grid->vertex[gamera_no_index3(
      grid->vertex_extent, upper[0], upper[1], upper[2])];
  gamera_no_vec3 midpoint;
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    midpoint.value[component] =
        0.5 * (start.value[component] + end.value[component]);
  }
  return midpoint;
}

static int initialize_emf_exchange(
    const gamera_no_grid *grid, const size_t active_lower[3],
    const size_t active_upper[3]) {
  if (emf_exchange_ready) {
    return 0;
  }
  if (grid == NULL || active_lower == NULL || active_upper == NULL ||
      patch_count != 2 || config.ni_global < 2 || config.nj_global < 2 ||
      config.nk_global < 2) {
    return -1;
  }
  size_t capacity = 0;
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    size_t extent[3];
    for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
      extent[axis] = active_upper[axis] - active_lower[axis] +
                     (axis == direction ? 0U : 1U);
    }
    capacity += extent[0] * extent[1] * extent[2];
  }
  edge_receptors =
      (edge_receptor_t *)calloc(capacity, sizeof(*edge_receptors));
  local_active_cells =
      (size_t)config.ni * (size_t)config.nj * (size_t)config.nk;
  if (edge_receptors == NULL || local_active_cells == 0
#ifndef GAMERA_YINYANG_SPARSE_OVERSET
      || local_active_cells > (size_t)INT32_MAX / GAMERA_NO_DIM
#endif
  ) {
    yinyang_exchange_destroy_local();
    return -1;
  }
#ifndef GAMERA_YINYANG_SPARSE_OVERSET
  emf_send_count = (int)(local_active_cells * GAMERA_NO_DIM);
  emf_send_buffer =
      (double *)malloc((size_t)emf_send_count * sizeof(*emf_send_buffer));
  emf_gather_buffer = (double *)malloc((size_t)emf_send_count * (size_t)size *
                                      sizeof(*emf_gather_buffer));
  if (emf_send_buffer == NULL || emf_gather_buffer == NULL) {
    yinyang_exchange_destroy_local();
    return -1;
  }
#endif

  edge_receptor_count = 0;
  gamera_no_yinyang_max_emf_donor_extrapolation = 0.0;
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const int transverse1 = (direction + 1) % GAMERA_NO_DIM;
    const int transverse2 = (direction + 2) % GAMERA_NO_DIM;
    size_t upper[3] = {active_upper[0], active_upper[1], active_upper[2]};
    ++upper[transverse1];
    ++upper[transverse2];
    for (size_t i = active_lower[0]; i < upper[0]; ++i) {
      for (size_t j = active_lower[1]; j < upper[1]; ++j) {
        for (size_t k = active_lower[2]; k < upper[2]; ++k) {
          const gamera_no_vec3 point =
              edge_midpoint(grid, direction, i, j, k);
          int is_receptor = 0;
          double blend_weight = 0.0;
          if (active_edge_transfer(point, &is_receptor, &blend_weight) != 0) {
            log_error("Yin-Yang EMF ownership lookup failed for patch %d "
                      "edge %d (%zu,%zu,%zu)",
                      patch_id, direction, i, j, k);
            yinyang_exchange_destroy_local();
            return -1;
          }
          if (!is_receptor) {
            continue;
          }
          const size_t edge = gamera_no_index3(
              grid->edge[direction].extent, i, j, k);
          if (grid->edge[direction].valid[edge] == 0U ||
              append_edge_receptor(&edge_receptor_count, capacity, direction,
                                   edge, blend_weight, point) != 0) {
            log_error("Yin-Yang EMF donor lookup failed for patch %d edge "
                      "%d (%zu,%zu,%zu), max extrapolation %.6e cells",
                      patch_id, direction, i, j, k,
                      gamera_no_yinyang_max_emf_donor_extrapolation);
            yinyang_exchange_destroy_local();
            return -1;
          }
        }
      }
    }
  }
  gamera_no_yinyang_edge_receptor_count = (int)edge_receptor_count;
  emf_exchange_ready = 1;
  log_info("Yin-Yang patch %d prepared %zu CT receptor edges; maximum "
           "donor extrapolation %.3e cells",
           patch_id, edge_receptor_count,
           gamera_no_yinyang_max_emf_donor_extrapolation);
  return 0;
}

static int reconstruct_cell_electric(
    const gamera_no_storage *storage, const gamera_no_grid *grid, size_t i,
    size_t j, size_t k, gamera_no_vec3 *electric) {
  gamera_no_vec3 row[GAMERA_NO_DIM] = {{{0.0}}};
  double projection[GAMERA_NO_DIM] = {0.0, 0.0, 0.0};
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    int transverse[2];
    int used = 0;
    for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
      if (axis != direction) {
        transverse[used++] = axis;
      }
    }
    for (int corner = 0; corner < 4; ++corner) {
      size_t coordinate[3] = {i, j, k};
      coordinate[transverse[0]] += (size_t)(corner & 1);
      coordinate[transverse[1]] += (size_t)((corner >> 1) & 1);
      const size_t edge = gamera_no_index3(
          grid->edge[direction].extent, coordinate[0], coordinate[1],
          coordinate[2]);
      const gamera_no_edge_geometry *geometry =
          &grid->edge[direction].value[edge];
      if (grid->edge[direction].valid[edge] == 0U ||
          !(geometry->length > 0.0)) {
        return -1;
      }
      projection[direction] +=
          0.25 * storage->edge_emf[direction][edge] / geometry->length;
      for (int component = 0; component < GAMERA_NO_DIM; ++component) {
        row[direction].value[component] +=
            0.25 * geometry->normal.value[component];
      }
    }
  }
  const gamera_no_vec3 cross12 = vector_cross(row[1], row[2]);
  const gamera_no_vec3 cross20 = vector_cross(row[2], row[0]);
  const gamera_no_vec3 cross01 = vector_cross(row[0], row[1]);
  const double determinant = vector_dot(row[0], cross12);
  if (!isfinite(determinant) || fabs(determinant) < 1.0e-12) {
    return -1;
  }
  for (int component = 0; component < GAMERA_NO_DIM; ++component) {
    electric->value[component] =
        (projection[0] * cross12.value[component] +
         projection[1] * cross20.value[component] +
         projection[2] * cross01.value[component]) /
        determinant;
    if (!isfinite(electric->value[component])) {
      return -1;
    }
  }
  return 0;
}

#ifndef GAMERA_YINYANG_SPARSE_OVERSET
static int pack_active_electric(const gamera_no_storage *storage,
                                const gamera_no_grid *grid) {
  for (int i = 0; i < config.ni; ++i) {
    for (int j = 0; j < config.nj; ++j) {
      for (int k = 0; k < config.nk; ++k) {
        gamera_no_vec3 electric;
        if (reconstruct_cell_electric(storage, grid, (size_t)(is + i),
                                      (size_t)(js + j), (size_t)(ks + k),
                                      &electric) != 0) {
          log_error("Yin-Yang electric-field reconstruction failed on patch "
                    "%d active cell (%d,%d,%d)",
                    patch_id, is + i, js + j, ks + k);
          return -1;
        }
        const size_t local = local_active_index(i, j, k);
        for (int component = 0; component < GAMERA_NO_DIM; ++component) {
          emf_send_buffer[(size_t)component * local_active_cells + local] =
              electric.value[component];
        }
      }
    }
  }
  return 0;
}
#endif

static int angular_face_coordinate_is_ghost(int face_direction, int axis,
                                            size_t coordinate) {
  const int local_count = axis == GAMERA_NO_J ? config.nj : config.nk;
  const int global_count =
      axis == GAMERA_NO_J ? config.nj_global : config.nk_global;
  const int global =
      proc_coords[axis] * local_count + (int)coordinate - NG;
  if (face_direction == axis) {
    return global < 0 || global > global_count;
  }
  return global < 0 || global >= global_count;
}

static int initialize_magnetic_exchange(const gamera_no_grid *grid) {
  if (magnetic_exchange_ready) {
    return 0;
  }
  if (grid == NULL || patch_count != 2) {
    return -1;
  }
  size_t capacity = 0;
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    capacity += gamera_no_element_count3(grid->face[direction].extent);
  }
  face_receptors =
      (face_receptor_t *)calloc(capacity, sizeof(*face_receptors));
  local_active_cells =
      (size_t)config.ni * (size_t)config.nj * (size_t)config.nk;
  if (face_receptors == NULL || local_active_cells == 0
#ifndef GAMERA_YINYANG_SPARSE_OVERSET
      || local_active_cells >
             (size_t)INT32_MAX / (2U * GAMERA_NO_DIM)
#endif
  ) {
    yinyang_exchange_destroy_local();
    return -1;
  }
#ifndef GAMERA_YINYANG_SPARSE_OVERSET
  magnetic_send_count =
      (int)(local_active_cells * 2U * GAMERA_NO_DIM);
  magnetic_send_buffer = (double *)malloc(
      (size_t)magnetic_send_count * sizeof(*magnetic_send_buffer));
  magnetic_gather_buffer = (double *)malloc(
      (size_t)magnetic_send_count * (size_t)size *
      sizeof(*magnetic_gather_buffer));
  if (magnetic_send_buffer == NULL || magnetic_gather_buffer == NULL) {
    yinyang_exchange_destroy_local();
    return -1;
  }
#endif

  face_receptor_count = 0;
  gamera_no_yinyang_max_magnetic_donor_extrapolation = 0.0;
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    for (size_t i = 0; i < grid->face[direction].extent[0]; ++i) {
      const int radial_active =
          direction == GAMERA_NO_I
              ? (i >= (size_t)is && i <= (size_t)(ie + 1))
              : (i >= (size_t)is && i <= (size_t)ie);
      if (!radial_active) {
        continue;
      }
      for (size_t j = 0; j < grid->face[direction].extent[1]; ++j) {
        for (size_t k = 0; k < grid->face[direction].extent[2]; ++k) {
          const int angular_ghost =
              angular_face_coordinate_is_ghost(direction, GAMERA_NO_J, j) ||
              angular_face_coordinate_is_ghost(direction, GAMERA_NO_K, k);
          if (!angular_ghost) {
            continue;
          }
          const size_t face = gamera_no_index3(
              grid->face[direction].extent, i, j, k);
          const gamera_no_vec3 point =
              grid->face[direction].value[face].centroid;
          if (append_face_receptor_to(face_receptors, &face_receptor_count,
                                      capacity, direction, face, 1.0,
                                      point) != 0) {
            log_error("Yin-Yang magnetic donor lookup failed for patch %d "
                      "face %d (%zu,%zu,%zu), max extrapolation %.6e cells",
                      patch_id, direction, i, j, k,
                      gamera_no_yinyang_max_magnetic_donor_extrapolation);
            yinyang_exchange_destroy_local();
            return -1;
          }
        }
      }
    }
  }
  gamera_no_yinyang_magnetic_face_receptor_count =
      (int)face_receptor_count;

  /* A ghost-face update can only change the cell-centred magnetic field in
   * cells sharing that face.  Cache that small, unique set once instead of
   * recovering B over the entire local block after every Yin-Yang exchange. */
  const size_t cell_count = gamera_no_element_count3(grid->cell_extent);
  unsigned char *marked = (unsigned char *)calloc(cell_count, sizeof(*marked));
  magnetic_recovery_cells =
      (size_t *)malloc(cell_count * sizeof(*magnetic_recovery_cells));
  if (marked == NULL || magnetic_recovery_cells == NULL) {
    free(marked);
    yinyang_exchange_destroy_local();
    return -1;
  }
  magnetic_recovery_cell_count = 0;
  for (size_t item = 0; item < face_receptor_count; ++item) {
    const face_receptor_t *receptor = &face_receptors[item];
    const size_t *extent = grid->face[receptor->direction].extent;
    const size_t plane = extent[1] * extent[2];
    const size_t i = receptor->receiver_face / plane;
    const size_t remainder = receptor->receiver_face % plane;
    const size_t j = remainder / extent[2];
    const size_t k = remainder % extent[2];
    const size_t face_coordinate[3] = {i, j, k};
    for (int side = 0; side < 2; ++side) {
      size_t cell_coordinate[3] = {i, j, k};
      if (side == 0) {
        if (face_coordinate[receptor->direction] == 0) {
          continue;
        }
        cell_coordinate[receptor->direction]--;
      }
      if (cell_coordinate[receptor->direction] >=
          grid->cell_extent[receptor->direction]) {
        continue;
      }
      const size_t cell = gamera_no_index3(
          grid->cell_extent, cell_coordinate[0], cell_coordinate[1],
          cell_coordinate[2]);
      if (!marked[cell]) {
        marked[cell] = 1;
        magnetic_recovery_cells[magnetic_recovery_cell_count++] = cell;
      }
    }
  }
  free(marked);
  magnetic_exchange_ready = 1;
  log_info("Yin-Yang patch %d prepared %zu angular-ghost magnetic faces; "
           "maximum donor extrapolation %.3e cells",
           patch_id, face_receptor_count,
           gamera_no_yinyang_max_magnetic_donor_extrapolation);
  return 0;
}

static int recover_magnetic_receptor_cells(
    const gamera_no_grid *grid, const double *const face_flux[3],
    gamera_no_vec3 *cell_magnetic) {
  int failed = 0;
#pragma omp parallel for reduction(| : failed) schedule(static)
  for (size_t item = 0; item < magnetic_recovery_cell_count; ++item) {
    const size_t cell = magnetic_recovery_cells[item];
    const size_t plane = grid->cell_extent[1] * grid->cell_extent[2];
    const size_t i = cell / plane;
    const size_t remainder = cell % plane;
    const size_t j = remainder / grid->cell_extent[2];
    const size_t k = remainder % grid->cell_extent[2];
    double local_flux[3][2];
    local_flux[GAMERA_NO_I][GAMERA_NO_LOWER] =
        face_flux[GAMERA_NO_I][gamera_no_index3(
            grid->face[GAMERA_NO_I].extent, i, j, k)];
    local_flux[GAMERA_NO_I][GAMERA_NO_UPPER] =
        face_flux[GAMERA_NO_I][gamera_no_index3(
            grid->face[GAMERA_NO_I].extent, i + 1, j, k)];
    local_flux[GAMERA_NO_J][GAMERA_NO_LOWER] =
        face_flux[GAMERA_NO_J][gamera_no_index3(
            grid->face[GAMERA_NO_J].extent, i, j, k)];
    local_flux[GAMERA_NO_J][GAMERA_NO_UPPER] =
        face_flux[GAMERA_NO_J][gamera_no_index3(
            grid->face[GAMERA_NO_J].extent, i, j + 1, k)];
    local_flux[GAMERA_NO_K][GAMERA_NO_LOWER] =
        face_flux[GAMERA_NO_K][gamera_no_index3(
            grid->face[GAMERA_NO_K].extent, i, j, k)];
    local_flux[GAMERA_NO_K][GAMERA_NO_UPPER] =
        face_flux[GAMERA_NO_K][gamera_no_index3(
            grid->face[GAMERA_NO_K].extent, i, j, k + 1)];
    if (gamera_no_flux_to_cell_field(&grid->cell[cell], local_flux,
                                     &cell_magnetic[cell], NULL) != 0) {
      failed = 1;
    }
  }
  return failed ? -1 : 0;
}

#ifdef GAMERA_YINYANG_SPARSE_OVERSET
typedef struct {
  const gamera_no_storage *storage;
  const gamera_no_grid *grid;
} sparse_state_context;

static int receptor_is_angular_active(const receptor_t *receptor) {
  if (receptor == NULL) {
    return 0;
  }
  const int global_j =
      proc_coords[1] * config.nj + receptor->receiver[1] - js;
  const int global_k =
      proc_coords[2] * config.nk + receptor->receiver[2] - ks;
  return global_j >= 0 && global_j < config.nj_global && global_k >= 0 &&
         global_k < config.nk_global;
}

static int sparse_reference_extent(size_t receptors_count, size_t *extent) {
  if (extent == NULL || receptors_count > SIZE_MAX / 8U) {
    return -1;
  }
  *extent = receptors_count * 8U;
  return 0;
}

static int sparse_donor_rank_range(int *begin, int *end) {
  if (begin == NULL || end == NULL || patch_count != 2 || patch_size <= 0 ||
      size != 2 * patch_size ||
      (patch_id != GAMERA_NO_YIN_PATCH && patch_id != GAMERA_NO_YANG_PATCH)) {
    return -1;
  }
  *begin = (1 - patch_id) * patch_size;
  *end = *begin + patch_size;
  return 0;
}

static void sparse_log_plan(const char *name,
                            const gamera_no_sparse_plan *plan) {
  if (name == NULL || plan == NULL || !plan->ready) {
    return;
  }
  log_info("Yin-Yang sparse plan rank=%d patch=%d kind=%s refs=%zu "
           "recv_unique=%zu requested_send_slots=%zu indegree=%d outdegree=%d "
           "channels=%d",
           rank, patch_id, name, plan->reference_count,
           plan->receive_cell_count, plan->send_cell_count, plan->indegree,
           plan->outdegree, plan->channel_count);
}

#ifdef GAMERA_YINYANG_SPARSE_OVERSET_PROFILE
static void sparse_profile_add_u64(uint64_t increment, uint64_t *total,
                                   int *saturated) {
  if (*total > UINT64_MAX - increment) {
    *total = UINT64_MAX;
    *saturated = 1;
  } else {
    *total += increment;
  }
}

static void sparse_profile_record(
    sparse_profile_accumulator *profile, gamera_no_sparse_plan *plan,
    const char *kind, const char *phase, const char *mode,
    double reconstruction_seconds, double total_seconds) {
  if (profile == NULL || plan == NULL || kind == NULL || phase == NULL ||
      mode == NULL || !plan->ready || plan->call_count == 0) {
    return;
  }
  const int warmup = plan->call_count == 1;
  const size_t receptor_count = plan->reference_count / 8U;
  log_info(
      "YINYANG_SPARSE_PROFILE "
      "{\"schema\":\"gamera.yinyang.sparse.v1\",\"record\":\"call\","
      "\"rank\":%d,\"patch\":%d,\"kind\":\"%s\",\"phase\":\"%s\","
      "\"mode\":\"%s\",\"call_index\":%llu,\"warmup\":%s,"
      "\"receptor_count\":%zu,\"corner_references\":%zu,"
      "\"unique_receive_cells\":%zu,\"send_cell_slots\":%zu,"
      "\"indegree\":%d,\"outdegree\":%d,\"channels\":%d,"
      "\"useful_send_bytes\":%llu,\"useful_receive_bytes\":%llu,"
      "\"pack_seconds\":%.17g,\"neighbor_seconds\":%.17g,"
      "\"reconstruction_seconds\":%.17g,\"total_seconds\":%.17g,"
      "\"neighbor_collectives\":1,\"world_consensus_calls\":2,"
      "\"statistics_saturated\":%s}",
      rank, patch_id, kind, phase, mode,
      (unsigned long long)plan->call_count, warmup ? "true" : "false",
      receptor_count, plan->reference_count, plan->receive_cell_count,
      plan->send_cell_count, plan->indegree, plan->outdegree,
      plan->channel_count, (unsigned long long)plan->send_bytes_per_call,
      (unsigned long long)plan->receive_bytes_per_call,
      plan->last_pack_seconds, plan->last_collective_seconds,
      reconstruction_seconds, total_seconds,
      plan->statistics_saturated ? "true" : "false");
  if (warmup) {
    return;
  }
  profile->kind = kind;
  profile->phase = phase;
  profile->mode = mode;
  sparse_profile_add_u64(1, &profile->samples, &profile->saturated);
  sparse_profile_add_u64(plan->send_bytes_per_call, &profile->send_bytes,
                         &profile->saturated);
  sparse_profile_add_u64(plan->receive_bytes_per_call,
                         &profile->receive_bytes, &profile->saturated);
  sparse_profile_add_u64(1, &profile->neighbor_collectives,
                         &profile->saturated);
  sparse_profile_add_u64(2, &profile->world_consensus_calls,
                         &profile->saturated);
  const double sample[GAMERA_NO_SPARSE_PROFILE_PHASE_COUNT] = {
      plan->last_pack_seconds, plan->last_collective_seconds,
      reconstruction_seconds, total_seconds};
  for (int item = 0; item < GAMERA_NO_SPARSE_PROFILE_PHASE_COUNT; ++item) {
    profile->sum[item] += sample[item];
    if (sample[item] > profile->maximum[item]) {
      profile->maximum[item] = sample[item];
    }
  }
}

static void sparse_profile_summary(
    const sparse_profile_accumulator *profile,
    const gamera_no_sparse_plan *plan) {
  if (profile == NULL || plan == NULL || profile->samples == 0 ||
      profile->kind == NULL || profile->phase == NULL ||
      profile->mode == NULL) {
    return;
  }
  const double inverse = 1.0 / (double)profile->samples;
  log_info(
      "YINYANG_SPARSE_PROFILE "
      "{\"schema\":\"gamera.yinyang.sparse.v1\","
      "\"record\":\"summary\",\"rank\":%d,\"patch\":%d,"
      "\"kind\":\"%s\",\"phase\":\"%s\",\"mode\":\"%s\","
      "\"samples\":%llu,\"warmup_excluded\":1,"
      "\"pack_mean_seconds\":%.17g,\"pack_max_seconds\":%.17g,"
      "\"neighbor_mean_seconds\":%.17g,\"neighbor_max_seconds\":%.17g,"
      "\"reconstruction_mean_seconds\":%.17g,"
      "\"reconstruction_max_seconds\":%.17g,"
      "\"total_mean_seconds\":%.17g,\"total_max_seconds\":%.17g,"
      "\"useful_send_bytes\":%llu,\"useful_receive_bytes\":%llu,"
      "\"neighbor_collectives\":%llu,\"world_consensus_calls\":%llu,"
      "\"statistics_saturated\":%s}",
      rank, patch_id, profile->kind, profile->phase, profile->mode,
      (unsigned long long)profile->samples,
      profile->sum[GAMERA_NO_SPARSE_PROFILE_PACK] * inverse,
      profile->maximum[GAMERA_NO_SPARSE_PROFILE_PACK],
      profile->sum[GAMERA_NO_SPARSE_PROFILE_NEIGHBOR] * inverse,
      profile->maximum[GAMERA_NO_SPARSE_PROFILE_NEIGHBOR],
      profile->sum[GAMERA_NO_SPARSE_PROFILE_RECONSTRUCT] * inverse,
      profile->maximum[GAMERA_NO_SPARSE_PROFILE_RECONSTRUCT],
      profile->sum[GAMERA_NO_SPARSE_PROFILE_TOTAL] * inverse,
      profile->maximum[GAMERA_NO_SPARSE_PROFILE_TOTAL],
      (unsigned long long)profile->send_bytes,
      (unsigned long long)profile->receive_bytes,
      (unsigned long long)profile->neighbor_collectives,
      (unsigned long long)profile->world_consensus_calls,
      (profile->saturated || plan->statistics_saturated) ? "true" : "false");
}
#endif

static int sparse_local_cell_coordinates(size_t local, int *i, int *j,
                                         int *k) {
  const size_t nj = (size_t)config.nj;
  const size_t nk = (size_t)config.nk;
  if (i == NULL || j == NULL || k == NULL || nj == 0 || nk == 0 ||
      nj > SIZE_MAX / nk || local >= local_active_cells) {
    return -1;
  }
  const size_t jk = nj * nk;
  const size_t qi = local / jk;
  const size_t remainder = local % jk;
  const size_t qj = remainder / nk;
  const size_t qk = remainder % nk;
  if (qi >= (size_t)config.ni || qj >= nj || qk >= nk ||
      qi > (size_t)INT32_MAX || qj > (size_t)INT32_MAX ||
      qk > (size_t)INT32_MAX) {
    return -1;
  }
  *i = (int)qi;
  *j = (int)qj;
  *k = (int)qk;
  return 0;
}

static int sparse_pack_fluid(const size_t *local_cells, size_t cell_count,
                             int channel_count, double *values,
                             void *opaque) {
  const sparse_state_context *context =
      (const sparse_state_context *)opaque;
  if (context == NULL || context->storage == NULL || context->grid == NULL ||
      channel_count !=
          GAMERA_NO_YINYANG_TIME_LEVELS * GAMERA_NO_FLUX_COUNT ||
      (cell_count > 0 && (local_cells == NULL || values == NULL))) {
    return -1;
  }
  for (size_t item = 0; item < cell_count; ++item) {
    int i;
    int j;
    int k;
    if (sparse_local_cell_coordinates(local_cells[item], &i, &j, &k) != 0) {
      return -1;
    }
    const size_t cell = gamera_no_index3(
        context->grid->cell_extent, (size_t)(is + i), (size_t)(js + j),
        (size_t)(ks + k));
#ifdef GAMERA_YINYANG_MFE_INTERFACE
    gamera_no_primitive primitive[GAMERA_NO_YINYANG_TIME_LEVELS];
    const double *conserved[GAMERA_NO_YINYANG_TIME_LEVELS] = {
        &context->storage->conserved[cell * GAMERA_NO_FLUX_COUNT],
        &context->storage->old_conserved[cell * GAMERA_NO_FLUX_COUNT]};
    for (int time_level = 0;
         time_level < GAMERA_NO_YINYANG_TIME_LEVELS; ++time_level) {
      if (gamera_no_conserved_to_primitive(
              conserved[time_level], gamma_val, rho_floor, p_floor,
              &primitive[time_level]) != 0) {
        return -1;
      }
    }
    const double level[GAMERA_NO_YINYANG_TIME_LEVELS]
                      [GAMERA_NO_FLUX_COUNT] = {
        {primitive[0].density, primitive[0].velocity.value[0],
         primitive[0].velocity.value[1], primitive[0].velocity.value[2],
         primitive[0].pressure},
        {primitive[1].density, primitive[1].velocity.value[0],
         primitive[1].velocity.value[1], primitive[1].velocity.value[2],
         primitive[1].pressure}};
#endif
    for (int time_level = 0;
         time_level < GAMERA_NO_YINYANG_TIME_LEVELS; ++time_level) {
      for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
        const size_t channel =
            (size_t)time_level * GAMERA_NO_FLUX_COUNT + (size_t)variable;
#ifdef GAMERA_YINYANG_MFE_INTERFACE
        values[item * (size_t)channel_count + channel] =
            level[time_level][variable];
#else
        const double *state = time_level == 0
                                  ? context->storage->conserved
                                  : context->storage->old_conserved;
        values[item * (size_t)channel_count + channel] =
            state[cell * GAMERA_NO_FLUX_COUNT + (size_t)variable];
#endif
      }
    }
  }
  return 0;
}

static int sparse_pack_electric(const size_t *local_cells, size_t cell_count,
                                int channel_count, double *values,
                                void *opaque) {
  const sparse_state_context *context =
      (const sparse_state_context *)opaque;
  if (context == NULL || context->storage == NULL || context->grid == NULL ||
      channel_count != GAMERA_NO_DIM ||
      (cell_count > 0 && (local_cells == NULL || values == NULL))) {
    return -1;
  }
  for (size_t item = 0; item < cell_count; ++item) {
    int i;
    int j;
    int k;
    gamera_no_vec3 electric;
    if (sparse_local_cell_coordinates(local_cells[item], &i, &j, &k) != 0 ||
        reconstruct_cell_electric(context->storage, context->grid,
                                  (size_t)(is + i), (size_t)(js + j),
                                  (size_t)(ks + k), &electric) != 0) {
      return -1;
    }
    for (int component = 0; component < GAMERA_NO_DIM; ++component) {
      values[item * (size_t)channel_count + (size_t)component] =
          electric.value[component];
    }
  }
  return 0;
}

static int sparse_pack_magnetic(const size_t *local_cells,
                                size_t cell_count, int channel_count,
                                double *values, void *opaque) {
  const sparse_state_context *context =
      (const sparse_state_context *)opaque;
  if (context == NULL || context->storage == NULL || context->grid == NULL ||
      (channel_count != GAMERA_NO_DIM &&
       channel_count != GAMERA_NO_YINYANG_TIME_LEVELS * GAMERA_NO_DIM) ||
      (cell_count > 0 && (local_cells == NULL || values == NULL))) {
    return -1;
  }
  const int level_count = channel_count / GAMERA_NO_DIM;
  for (size_t item = 0; item < cell_count; ++item) {
    int i;
    int j;
    int k;
    if (sparse_local_cell_coordinates(local_cells[item], &i, &j, &k) != 0) {
      return -1;
    }
    const size_t cell = gamera_no_index3(
        context->grid->cell_extent, (size_t)(is + i), (size_t)(js + j),
        (size_t)(ks + k));
    const gamera_no_vec3 level[GAMERA_NO_YINYANG_TIME_LEVELS] = {
        context->storage->cell_magnetic[cell],
        context->storage->old_cell_magnetic[cell]};
    for (int time_level = 0; time_level < level_count; ++time_level) {
      for (int component = 0; component < GAMERA_NO_DIM; ++component) {
        const size_t channel =
            (size_t)time_level * GAMERA_NO_DIM + (size_t)component;
        const double value = level[time_level].value[component];
        if (!isfinite(value)) {
          return -1;
        }
        values[item * (size_t)channel_count + channel] = value;
      }
    }
  }
  return 0;
}

static int initialize_sparse_fluid_plans(void) {
  if (fluid_full_plan.ready && fluid_ghost_plan.ready &&
      fluid_ghost_reference_base != NULL) {
    return 0;
  }
  if (fluid_full_plan.ready || fluid_ghost_plan.ready ||
      fluid_ghost_reference_base != NULL) {
    return -1;
  }
  size_t full_extent = 0;
  size_t ghost_receptor_count = 0;
  int local_failed =
      sparse_reference_extent(receptor_count, &full_extent) != 0;
  for (size_t item = 0; !local_failed && item < receptor_count; ++item) {
    if (!receptor_is_angular_active(&receptors[item])) {
      ++ghost_receptor_count;
    }
  }
  size_t ghost_extent = 0;
  local_failed |=
      sparse_reference_extent(ghost_receptor_count, &ghost_extent) != 0;
  gamera_no_sparse_reference *full =
      (gamera_no_sparse_reference *)calloc(full_extent == 0 ? 1 : full_extent,
                                            sizeof(*full));
  gamera_no_sparse_reference *ghost =
      (gamera_no_sparse_reference *)calloc(
          ghost_extent == 0 ? 1 : ghost_extent, sizeof(*ghost));
  fluid_ghost_reference_base = (size_t *)malloc(
      (receptor_count == 0 ? 1 : receptor_count) *
      sizeof(*fluid_ghost_reference_base));
  local_failed |= full == NULL || ghost == NULL ||
                  fluid_ghost_reference_base == NULL;
  size_t ghost_used = 0;
  for (size_t item = 0; !local_failed && item < receptor_count; ++item) {
    fluid_ghost_reference_base[item] = SIZE_MAX;
    const int is_ghost = !receptor_is_angular_active(&receptors[item]);
    if (is_ghost) {
      fluid_ghost_reference_base[item] = ghost_used;
    }
    for (int corner = 0; corner < 8; ++corner) {
      const size_t full_index = item * 8U + (size_t)corner;
      full[full_index].donor_rank = receptors[item].donor_rank[corner];
      full[full_index].donor_cell = receptors[item].donor_cell[corner];
      if (is_ghost) {
        ghost[ghost_used].donor_rank = receptors[item].donor_rank[corner];
        ghost[ghost_used].donor_cell = receptors[item].donor_cell[corner];
        ++ghost_used;
      }
    }
  }
  local_failed |= ghost_used != ghost_extent;
  int donor_begin = 0;
  int donor_end = 0;
  local_failed |= sparse_donor_rank_range(&donor_begin, &donor_end) != 0;
  if (world_consensus_failure(local_failed,
                              "sparse fluid request construction") != 0) {
    free(full);
    free(ghost);
    return -1;
  }
  int status = gamera_no_sparse_plan_build(
      &fluid_full_plan, MPI_COMM_WORLD, full, full_extent,
      local_active_cells,
      GAMERA_NO_YINYANG_TIME_LEVELS * GAMERA_NO_FLUX_COUNT, donor_begin,
      donor_end);
  if (status == 0) {
    status = gamera_no_sparse_plan_build(
        &fluid_ghost_plan, MPI_COMM_WORLD, ghost, ghost_extent,
        local_active_cells,
        GAMERA_NO_YINYANG_TIME_LEVELS * GAMERA_NO_FLUX_COUNT, donor_begin,
        donor_end);
  }
  free(full);
  free(ghost);
  if (status != 0) {
    return -1;
  }
  sparse_log_plan("fluid_full", &fluid_full_plan);
  sparse_log_plan("fluid_ghost", &fluid_ghost_plan);
  return 0;
}

static int initialize_sparse_edge_plan(void) {
  if (emf_plan.ready) {
    return 0;
  }
  size_t extent = 0;
  int local_failed =
      sparse_reference_extent(edge_receptor_count, &extent) != 0;
  gamera_no_sparse_reference *references =
      (gamera_no_sparse_reference *)calloc(extent == 0 ? 1 : extent,
                                            sizeof(*references));
  local_failed |= references == NULL;
  for (size_t item = 0; !local_failed && item < edge_receptor_count; ++item) {
    for (int corner = 0; corner < 8; ++corner) {
      const size_t index = item * 8U + (size_t)corner;
      references[index].donor_rank = edge_receptors[item].donor_rank[corner];
      references[index].donor_cell = edge_receptors[item].donor_cell[corner];
    }
  }
  int donor_begin = 0;
  int donor_end = 0;
  local_failed |= sparse_donor_rank_range(&donor_begin, &donor_end) != 0;
  if (world_consensus_failure(local_failed,
                              "sparse EMF request construction") != 0) {
    free(references);
    return -1;
  }
  const int status = gamera_no_sparse_plan_build(
      &emf_plan, MPI_COMM_WORLD, references, extent, local_active_cells,
      GAMERA_NO_DIM, donor_begin, donor_end);
  free(references);
  if (status != 0) {
    return -1;
  }
  sparse_log_plan("emf", &emf_plan);
  return 0;
}

static int initialize_sparse_magnetic_ghost_plan(void) {
  if (magnetic_ghost_plan.ready) {
    return 0;
  }
  size_t extent = 0;
  int local_failed =
      sparse_reference_extent(face_receptor_count, &extent) != 0;
  gamera_no_sparse_reference *references =
      (gamera_no_sparse_reference *)calloc(extent == 0 ? 1 : extent,
                                            sizeof(*references));
  local_failed |= references == NULL;
  for (size_t item = 0; !local_failed && item < face_receptor_count; ++item) {
    for (int corner = 0; corner < 8; ++corner) {
      const size_t index = item * 8U + (size_t)corner;
      references[index].donor_rank = face_receptors[item].donor_rank[corner];
      references[index].donor_cell = face_receptors[item].donor_cell[corner];
    }
  }
  int donor_begin = 0;
  int donor_end = 0;
  local_failed |= sparse_donor_rank_range(&donor_begin, &donor_end) != 0;
  if (world_consensus_failure(
          local_failed, "sparse magnetic-ghost request construction") != 0) {
    free(references);
    return -1;
  }
  const int status = gamera_no_sparse_plan_build(
      &magnetic_ghost_plan, MPI_COMM_WORLD, references, extent,
      local_active_cells, GAMERA_NO_YINYANG_TIME_LEVELS * GAMERA_NO_DIM,
      donor_begin, donor_end);
  free(references);
  if (status != 0) {
    return -1;
  }
  sparse_log_plan("magnetic_ghost", &magnetic_ghost_plan);
  return 0;
}

#ifdef GAMERA_YINYANG_HDIV_RECONCILE
static int initialize_sparse_active_magnetic_plans(void) {
  if (magnetic_active_both_plan.ready &&
      magnetic_active_current_plan.ready) {
    return 0;
  }
  if (magnetic_active_both_plan.ready ||
      magnetic_active_current_plan.ready) {
    return -1;
  }
  size_t extent = 0;
  int local_failed =
      sparse_reference_extent(active_face_receptor_count, &extent) != 0;
  gamera_no_sparse_reference *references =
      (gamera_no_sparse_reference *)calloc(extent == 0 ? 1 : extent,
                                            sizeof(*references));
  local_failed |= references == NULL;
  for (size_t item = 0; !local_failed && item < active_face_receptor_count;
       ++item) {
    for (int corner = 0; corner < 8; ++corner) {
      const size_t index = item * 8U + (size_t)corner;
      references[index].donor_rank =
          active_face_receptors[item].donor_rank[corner];
      references[index].donor_cell =
          active_face_receptors[item].donor_cell[corner];
    }
  }
  int donor_begin = 0;
  int donor_end = 0;
  local_failed |= sparse_donor_rank_range(&donor_begin, &donor_end) != 0;
  if (world_consensus_failure(
          local_failed, "sparse active-magnetic request construction") != 0) {
    free(references);
    return -1;
  }
  int status = gamera_no_sparse_plan_build(
      &magnetic_active_both_plan, MPI_COMM_WORLD, references, extent,
      local_active_cells, GAMERA_NO_YINYANG_TIME_LEVELS * GAMERA_NO_DIM,
      donor_begin, donor_end);
  if (status == 0) {
    status = gamera_no_sparse_plan_build(
        &magnetic_active_current_plan, MPI_COMM_WORLD, references, extent,
        local_active_cells, GAMERA_NO_DIM, donor_begin, donor_end);
  }
  free(references);
  if (status != 0) {
    return -1;
  }
  sparse_log_plan("magnetic_active_both", &magnetic_active_both_plan);
  sparse_log_plan("magnetic_active_current", &magnetic_active_current_plan);
  return 0;
}
#endif
#endif

#ifdef GAMERA_YINYANG_HDIV_RECONCILE
enum {
  GAMERA_NO_HDIV_FACE_SIDES = 2,
  GAMERA_NO_HDIV_TARGET_CHANNELS =
      GAMERA_NO_YINYANG_TIME_LEVELS * GAMERA_NO_DIM *
      GAMERA_NO_HDIV_FACE_SIDES,
  GAMERA_NO_HDIV_BASE_CHANNELS =
      GAMERA_NO_YINYANG_TIME_LEVELS * GAMERA_NO_DIM *
      GAMERA_NO_HDIV_FACE_SIDES,
  GAMERA_NO_HDIV_WEIGHT_CHANNELS =
      GAMERA_NO_DIM * GAMERA_NO_HDIV_FACE_SIDES,
  GAMERA_NO_HDIV_CHANNELS =
      GAMERA_NO_HDIV_TARGET_CHANNELS + GAMERA_NO_HDIV_BASE_CHANNELS +
      GAMERA_NO_HDIV_WEIGHT_CHANNELS
};

#ifndef GAMERA_YINYANG_HDIV_DISTRIBUTED
static size_t hdiv_global_cell_index(int i, int j, int k) {
  return ((size_t)i * (size_t)config.nj_global + (size_t)j) *
             (size_t)config.nk_global +
         (size_t)k;
}

static size_t hdiv_global_face_index(int direction, int i, int j, int k) {
  return gamera_no_index3(hdiv_global_face_extent[direction], (size_t)i,
                          (size_t)j, (size_t)k);
}

static double hdiv_cell_diagonal(int i, int j, int k);
#else
static size_t hdiv_local_cell_index(size_t i, size_t j, size_t k) {
  return gamera_no_index3(hdiv_local_cell_extent, i, j, k);
}

static size_t hdiv_local_face_index(int direction, size_t i, size_t j,
                                    size_t k) {
  return gamera_no_index3(hdiv_local_face_extent[direction], i, j, k);
}

static size_t hdiv_active_cell_index(int i, int j, int k) {
  return hdiv_local_cell_index((size_t)(is + i), (size_t)(js + j),
                               (size_t)(ks + k));
}
#endif

static size_t hdiv_target_channel(int time_level, int direction, int side) {
  return ((size_t)time_level * GAMERA_NO_DIM + (size_t)direction) *
             GAMERA_NO_HDIV_FACE_SIDES +
         (size_t)side;
}

static size_t hdiv_base_channel(int time_level, int direction, int side) {
  return (size_t)GAMERA_NO_HDIV_TARGET_CHANNELS +
         ((size_t)time_level * GAMERA_NO_DIM + (size_t)direction) *
             GAMERA_NO_HDIV_FACE_SIDES +
         (size_t)side;
}

static size_t hdiv_weight_channel(int direction, int side) {
  return (size_t)GAMERA_NO_HDIV_TARGET_CHANNELS +
         (size_t)GAMERA_NO_HDIV_BASE_CHANNELS +
         (size_t)direction * GAMERA_NO_HDIV_FACE_SIDES + (size_t)side;
}

/*
 * Diagonal finite-volume Hodge star for an integrated face flux.  With
 * w_f=A_f/d_f, minimizing sum(delta-Phi_f^2/w_f) approximates minimizing
 * the magnetic correction energy integral |delta B|^2 dV.  This matters on
 * the strongly stretched Earth grid: a unit-flux correction on a small
 * inner face must not be treated like the same correction at 100+ RE.
 *
 * Face normals point in the positive logical direction.  Interior distances
 * use adjacent cell centroids; an angular overset boundary uses the
 * one-sided center-to-face distance because lambda=0 is imposed at that
 * boundary.  The two physical radial-normal face planes are immobile.
 */
static int hdiv_metric_face_weight(const gamera_no_grid *grid, int direction,
                                   int side,
                                   const size_t cell_coordinate[3],
                                   const int global_cell[3], double *weight) {
  const int global_count[3] = {config.ni_global, config.nj_global,
                               config.nk_global};
  const int normal_coordinate = global_cell[direction] + side;
  if (grid == NULL || weight == NULL || direction < 0 ||
      direction >= GAMERA_NO_DIM ||
      (side != GAMERA_NO_LOWER && side != GAMERA_NO_UPPER) ||
      normal_coordinate < 0 ||
      normal_coordinate > global_count[direction]) {
    return -1;
  }
  if (direction == GAMERA_NO_I &&
      (normal_coordinate == 0 ||
       normal_coordinate == global_count[direction])) {
    *weight = 0.0;
    return 0;
  }

  size_t face_coordinate[3] = {cell_coordinate[0], cell_coordinate[1],
                               cell_coordinate[2]};
  face_coordinate[direction] += (size_t)side;
  const size_t face = gamera_no_index3(
      grid->face[direction].extent, face_coordinate[0], face_coordinate[1],
      face_coordinate[2]);
  const gamera_no_face_geometry *geometry =
      &grid->face[direction].value[face];

  gamera_no_vec3 displacement = {{0.0, 0.0, 0.0}};
  if (normal_coordinate == 0 ||
      normal_coordinate == global_count[direction]) {
    const size_t cell = gamera_no_index3(
        grid->cell_extent, cell_coordinate[0], cell_coordinate[1],
        cell_coordinate[2]);
    for (int component = 0; component < GAMERA_NO_DIM; ++component) {
      displacement.value[component] =
          normal_coordinate == 0
              ? grid->cell[cell].centroid.value[component] -
                    geometry->centroid.value[component]
              : geometry->centroid.value[component] -
                    grid->cell[cell].centroid.value[component];
    }
  } else {
    size_t lower[3] = {cell_coordinate[0], cell_coordinate[1],
                       cell_coordinate[2]};
    size_t upper[3] = {cell_coordinate[0], cell_coordinate[1],
                       cell_coordinate[2]};
    if (side == GAMERA_NO_LOWER) {
      --lower[direction];
    } else {
      ++upper[direction];
    }
    const size_t lower_cell = gamera_no_index3(
        grid->cell_extent, lower[0], lower[1], lower[2]);
    const size_t upper_cell = gamera_no_index3(
        grid->cell_extent, upper[0], upper[1], upper[2]);
    for (int component = 0; component < GAMERA_NO_DIM; ++component) {
      displacement.value[component] =
          grid->cell[upper_cell].centroid.value[component] -
          grid->cell[lower_cell].centroid.value[component];
    }
  }

  const double distance = vector_dot(displacement, geometry->normal);
  const double displacement_norm = sqrt(vector_dot(displacement, displacement));
  const double scale = fmax(displacement_norm, DBL_MIN);
  const double cosine = distance / scale;
  if (!isfinite(geometry->area) || geometry->area <= 0.0 ||
      !isfinite(distance) ||
      distance <= 128.0 * DBL_EPSILON * scale || !isfinite(cosine) ||
      cosine <= 128.0 * DBL_EPSILON) {
    log_error("Invalid H(div) metric on patch %d face d=%d side=%d "
              "global=(%d,%d,%d): area=%.17g distance=%.17g cos=%.17g",
              patch_id, direction, side, global_cell[0], global_cell[1],
              global_cell[2], geometry->area, distance, cosine);
    return -1;
  }
  *weight = geometry->area / distance;
  if (!isfinite(*weight) || *weight <= 0.0) {
    return -1;
  }
  gamera_no_yinyang_hdiv_min_metric_cosine =
      fmin(gamera_no_yinyang_hdiv_min_metric_cosine, cosine);
  return 0;
}

#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
/*
 * The grid is stationary, so its finite-volume Hodge weights never change
 * during a run.  Keep them in the persistent wire-buffer channels instead of
 * rebuilding six centroid/normal metrics for every active cell on every
 * timestep.  This initialization remains serial because
 * hdiv_metric_face_weight also records a diagnostic minimum cosine.
 */
static int cache_hdiv_local_metric_weights(const gamera_no_grid *grid) {
  for (int i = 0; i < config.ni; ++i) {
    for (int j = 0; j < config.nj; ++j) {
      for (int k = 0; k < config.nk; ++k) {
        const size_t local = local_active_index(i, j, k);
        const size_t lower[3] = {(size_t)(is + i), (size_t)(js + j),
                                 (size_t)(ks + k)};
        const int global[3] = {proc_coords[0] * config.ni + i,
                               proc_coords[1] * config.nj + j,
                               proc_coords[2] * config.nk + k};
        for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
          for (int side = 0; side < GAMERA_NO_HDIV_FACE_SIDES; ++side) {
            double metric_weight = 0.0;
            if (hdiv_metric_face_weight(grid, direction, side, lower, global,
                                        &metric_weight) != 0) {
              return -1;
            }
            hdiv_send_buffer[hdiv_weight_channel(direction, side) *
                                 local_active_cells +
                             local] = metric_weight;
          }
        }
      }
    }
  }
  return 0;
}
#endif

static int hdiv_face_touches_protected_radial_shell(int direction,
                                                     const int global[3]) {
  if (direction == GAMERA_NO_I) {
    return global[GAMERA_NO_I] <= 1 ||
           global[GAMERA_NO_I] >= config.ni_global - 1;
  }
  return global[GAMERA_NO_I] == 0 ||
         global[GAMERA_NO_I] == config.ni_global - 1;
}

static int initialize_hdiv_exchange(const gamera_no_grid *grid) {
  if (hdiv_exchange_ready) {
    return 0;
  }
  if (grid == NULL || patch_count != 2 || config.ni_global < 3 ||
      config.nj_global < 2 || config.nk_global < 2) {
    return -1;
  }
#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED
  /*
   * The local layout is dimension-generic, but the radial shared-face metric
   * path has not yet passed its independent qualification.  Reject that
   * capability explicitly instead of silently relying on a topology.
   */
  if (config.proc_dims[GAMERA_NO_I] != 1) {
    if (patch_rank == 0) {
      log_error("Distributed Yin-Yang H(div) currently requires "
                "proc_dims_i=1; requested %d",
                config.proc_dims[GAMERA_NO_I]);
    }
    return -1;
  }
#endif

  size_t active_capacity = 0;
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    size_t count[3] = {(size_t)config.ni, (size_t)config.nj,
                       (size_t)config.nk};
    ++count[direction];
    active_capacity += count[0] * count[1] * count[2];
  }
  active_face_receptors = (face_receptor_t *)calloc(
      active_capacity, sizeof(*active_face_receptors));
  if (active_face_receptors == NULL || local_active_cells == 0 ||
      local_active_cells > (size_t)INT32_MAX / GAMERA_NO_HDIV_CHANNELS) {
    yinyang_exchange_destroy_local();
    return -1;
  }
  hdiv_send_count = (int)(local_active_cells * GAMERA_NO_HDIV_CHANNELS);
  hdiv_send_buffer =
      (double *)malloc((size_t)hdiv_send_count * sizeof(*hdiv_send_buffer));
#ifndef GAMERA_YINYANG_HDIV_DISTRIBUTED
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
  const size_t hdiv_gather_ranks = (size_t)patch_size;
#else
  const size_t hdiv_gather_ranks = (size_t)size;
#endif
  hdiv_gather_buffer = (double *)malloc(
      (size_t)hdiv_send_count * hdiv_gather_ranks *
      sizeof(*hdiv_gather_buffer));
  if (hdiv_send_buffer == NULL || hdiv_gather_buffer == NULL) {
#else
  if (hdiv_send_buffer == NULL) {
#endif
    yinyang_exchange_destroy_local();
    return -1;
  }

  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t local_face_count =
        gamera_no_element_count3(grid->face[direction].extent);
    for (int time_level = 0; time_level < GAMERA_NO_YINYANG_TIME_LEVELS;
         ++time_level) {
      hdiv_target_face_flux[time_level][direction] = (double *)malloc(
          local_face_count * sizeof(*hdiv_target_face_flux[time_level]
                                                     [direction]));
    }
#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED
    for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
      hdiv_local_face_extent[direction][axis] =
          grid->face[direction].extent[axis];
    }
    for (int time_level = 0;
         time_level < GAMERA_NO_YINYANG_TIME_LEVELS; ++time_level) {
      hdiv_local_face_base[time_level][direction] =
          (double *)calloc(local_face_count,
                           sizeof(*hdiv_local_face_base[time_level]
                                                        [direction]));
      hdiv_local_face_flux[time_level][direction] =
          (double *)calloc(local_face_count,
                           sizeof(*hdiv_local_face_flux[time_level]
                                                        [direction]));
    }
    hdiv_local_face_mobility[direction] =
        (double *)calloc(local_face_count,
                         sizeof(*hdiv_local_face_mobility[direction]));
    if (hdiv_target_face_flux[0][direction] == NULL ||
        hdiv_target_face_flux[1][direction] == NULL ||
        hdiv_local_face_base[0][direction] == NULL ||
        hdiv_local_face_base[1][direction] == NULL ||
        hdiv_local_face_flux[0][direction] == NULL ||
        hdiv_local_face_flux[1][direction] == NULL ||
        hdiv_local_face_mobility[direction] == NULL) {
      yinyang_exchange_destroy_local();
      return -1;
    }
#else
    for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
      hdiv_global_face_extent[direction][axis] =
          (size_t)(axis == GAMERA_NO_I
                       ? config.ni_global
                       : (axis == GAMERA_NO_J ? config.nj_global
                                              : config.nk_global)) +
          (axis == direction ? 1U : 0U);
    }
    const size_t global_face_count = gamera_no_element_count3(
        hdiv_global_face_extent[direction]);
    for (int time_level = 0; time_level < GAMERA_NO_YINYANG_TIME_LEVELS;
         ++time_level) {
      hdiv_global_face_base[time_level][direction] = (double *)malloc(
          global_face_count *
          sizeof(*hdiv_global_face_base[time_level][direction]));
      hdiv_global_face_flux[time_level][direction] = (double *)malloc(
          global_face_count *
          sizeof(*hdiv_global_face_flux[time_level][direction]));
    }
    hdiv_global_face_mobility[direction] = (double *)malloc(
        global_face_count * sizeof(*hdiv_global_face_mobility[direction]));
    if (hdiv_target_face_flux[0][direction] == NULL ||
        hdiv_target_face_flux[1][direction] == NULL ||
        hdiv_global_face_base[0][direction] == NULL ||
        hdiv_global_face_base[1][direction] == NULL ||
        hdiv_global_face_flux[0][direction] == NULL ||
        hdiv_global_face_flux[1][direction] == NULL ||
        hdiv_global_face_mobility[direction] == NULL) {
      yinyang_exchange_destroy_local();
      return -1;
    }
#endif
  }

#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED
  for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
    hdiv_local_cell_extent[axis] = grid->cell_extent[axis];
  }
  const size_t hdiv_cell_count =
      gamera_no_element_count3(hdiv_local_cell_extent);

  hdiv_face_halo_capacity = 0;
  const int local_cells[3] = {config.ni, config.nj, config.nk};
  const int local_ghost_lower[3] = {is, js, ks};
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
      const size_t active =
          (size_t)local_cells[axis] + (axis == direction ? 1U : 0U);
      const size_t lower = (size_t)local_ghost_lower[axis];
      if (lower + active > hdiv_local_face_extent[direction][axis]) {
        yinyang_exchange_destroy_local();
        return -1;
      }
      const size_t upper =
          hdiv_local_face_extent[direction][axis] - lower - active;
      const size_t width_to_low =
          upper + (axis == direction ? 1U : 0U);
      const size_t width_to_high = lower;
      if (local_cells[axis] <= 0 || width_to_low > active ||
          width_to_high > (size_t)local_cells[axis]) {
        log_error("Distributed H(div) face halo is wider than the local "
                  "active block on patch %d axis=%d direction=%d: "
                  "cells=%d low_width=%zu high_width=%zu",
                  patch_id, axis, direction, local_cells[axis], width_to_low,
                  width_to_high);
        yinyang_exchange_destroy_local();
        return -1;
      }
      size_t plane = 1;
      for (int transverse = 0; transverse < GAMERA_NO_DIM; ++transverse) {
        if (transverse != axis) {
          plane *= hdiv_local_face_extent[direction][transverse];
        }
      }
      const size_t count =
          plane * (width_to_low > width_to_high ? width_to_low
                                                : width_to_high);
      hdiv_face_halo_capacity =
          hdiv_face_halo_capacity > count ? hdiv_face_halo_capacity : count;
    }
  }
  hdiv_cell_halo_capacity = 0;
  for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
    if (local_ghost_lower[axis] < 1 ||
        (size_t)local_ghost_lower[axis] + (size_t)local_cells[axis] >=
            hdiv_local_cell_extent[axis]) {
      log_error("Distributed H(div) requires one cell halo on both sides "
                "of axis %d (extent=%zu lower=%d cells=%d)",
                axis, hdiv_local_cell_extent[axis], local_ghost_lower[axis],
                local_cells[axis]);
      yinyang_exchange_destroy_local();
      return -1;
    }
    size_t plane = 1;
    for (int transverse = 0; transverse < GAMERA_NO_DIM; ++transverse) {
      if (transverse != axis) {
        plane *= (size_t)local_cells[transverse];
      }
    }
    hdiv_cell_halo_capacity =
        hdiv_cell_halo_capacity > plane ? hdiv_cell_halo_capacity : plane;
  }
  hdiv_face_halo_send =
      (double *)malloc(hdiv_face_halo_capacity * sizeof(*hdiv_face_halo_send));
  hdiv_face_halo_receive = (double *)malloc(
      hdiv_face_halo_capacity * sizeof(*hdiv_face_halo_receive));
  hdiv_cell_halo_send =
      (double *)malloc(hdiv_cell_halo_capacity * sizeof(*hdiv_cell_halo_send));
  hdiv_cell_halo_receive = (double *)malloc(
      hdiv_cell_halo_capacity * sizeof(*hdiv_cell_halo_receive));
  /*
   * Query the predefined attribute on MPI_COMM_WORLD.  Some MPI libraries do
   * not copy it onto a derived Cartesian communicator even though the same
   * implementation-wide tag limit applies there.
   */
  int *tag_upper = NULL;
  int tag_attribute_present = 0;
  if (MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &tag_upper,
                        &tag_attribute_present) != MPI_SUCCESS ||
      !tag_attribute_present || tag_upper == NULL || *tag_upper < 7705) {
    log_error("Distributed H(div) MPI tag limit is unavailable or below "
              "the required tag 7705");
    yinyang_exchange_destroy_local();
    return -1;
  }
#else
  const size_t hdiv_cell_count =
      (size_t)config.ni_global * (size_t)config.nj_global *
      (size_t)config.nk_global;
#endif
  hdiv_lambda = (double *)calloc(hdiv_cell_count, sizeof(*hdiv_lambda));
  hdiv_residual =
      (double *)malloc(hdiv_cell_count * sizeof(*hdiv_residual));
  hdiv_preconditioned =
      (double *)malloc(hdiv_cell_count * sizeof(*hdiv_preconditioned));
  hdiv_search =
      (double *)malloc(hdiv_cell_count * sizeof(*hdiv_search));
  hdiv_operator_search =
      (double *)malloc(hdiv_cell_count * sizeof(*hdiv_operator_search));
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
  hdiv_cell_diagonal_cache =
      (double *)malloc(hdiv_cell_count *
                       sizeof(*hdiv_cell_diagonal_cache));
  /*
   * Use one deterministic dot-product block per global radial slab.  The
   * partition therefore depends only on the runtime mesh, never on the
   * OpenMP team size.  Each block is accumulated in its original flat-index
   * order and the block partials are combined in radial order below.
   */
  hdiv_dot_block_count =
#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED
      (size_t)config.ni;
#else
      (size_t)config.ni_global;
#endif
  if (hdiv_dot_block_count > hdiv_cell_count) {
    hdiv_dot_block_count = hdiv_cell_count;
  }
  hdiv_dot_block_partial =
      (double *)malloc(hdiv_dot_block_count *
                       sizeof(*hdiv_dot_block_partial));
#endif
  if (hdiv_lambda == NULL || hdiv_residual == NULL ||
      hdiv_preconditioned == NULL || hdiv_search == NULL ||
      hdiv_operator_search == NULL
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
      || hdiv_cell_diagonal_cache == NULL || hdiv_dot_block_count == 0 ||
      hdiv_dot_block_partial == NULL
#endif
#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED
      || hdiv_face_halo_capacity == 0 || hdiv_cell_halo_capacity == 0 ||
      hdiv_face_halo_send == NULL || hdiv_face_halo_receive == NULL ||
      hdiv_cell_halo_send == NULL || hdiv_cell_halo_receive == NULL
#endif
  ) {
    yinyang_exchange_destroy_local();
    return -1;
  }

  active_face_receptor_count = 0;
  const size_t local_count[3] = {(size_t)config.ni, (size_t)config.nj,
                                 (size_t)config.nk};
  const size_t local_lower[3] = {(size_t)is, (size_t)js, (size_t)ks};
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    size_t count[3] = {local_count[0], local_count[1], local_count[2]};
    ++count[direction];
    for (size_t qi = 0; qi < count[0]; ++qi) {
      for (size_t qj = 0; qj < count[1]; ++qj) {
        for (size_t qk = 0; qk < count[2]; ++qk) {
          const size_t coordinate[3] = {local_lower[0] + qi,
                                        local_lower[1] + qj,
                                        local_lower[2] + qk};
          const int global[3] = {
              proc_coords[0] * config.ni + (int)qi,
              proc_coords[1] * config.nj + (int)qj,
              proc_coords[2] * config.nk + (int)qk};
          const size_t face = gamera_no_index3(
              grid->face[direction].extent, coordinate[0], coordinate[1],
              coordinate[2]);
          int target = 0;
          double blend_weight = 0.0;
          if (!hdiv_face_touches_protected_radial_shell(direction, global) &&
              active_magnetic_face_transfer(
                  grid->face[direction].value[face].centroid, &target,
                  &blend_weight) != 0) {
            log_error("Yin-Yang H(div) ownership lookup failed for patch %d "
                      "face %d (%zu,%zu,%zu)",
                      patch_id, direction, coordinate[0], coordinate[1],
                      coordinate[2]);
            yinyang_exchange_destroy_local();
            return -1;
          }
          if (target &&
              append_face_receptor_to(
                  active_face_receptors, &active_face_receptor_count,
                  active_capacity, direction, face, blend_weight,
                  grid->face[direction].value[face].centroid) != 0) {
            log_error("Yin-Yang H(div) donor lookup failed for patch %d face "
                      "%d (%zu,%zu,%zu)",
                      patch_id, direction, coordinate[0], coordinate[1],
                      coordinate[2]);
            yinyang_exchange_destroy_local();
            return -1;
          }
        }
      }
    }
  }
  gamera_no_yinyang_active_magnetic_face_receptor_count =
      (int)active_face_receptor_count;
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
  if (cache_hdiv_local_metric_weights(grid) != 0) {
    yinyang_exchange_destroy_local();
    return -1;
  }
  hdiv_cell_diagonal_cache_ready = 0;
#endif
  hdiv_exchange_ready = 1;
  log_info("Yin-Yang patch %d prepared %zu active magnetic face targets for "
           "H(div) reconciliation",
           patch_id, active_face_receptor_count);
  return 0;
}

static int prepare_hdiv_face_targets(const gamera_no_storage *storage,
                                     const gamera_no_grid *grid
#ifdef GAMERA_YINYANG_SPARSE_OVERSET
                                     ,
                                     const gamera_no_sparse_plan *donor_plan
#endif
                                     ) {
  const int level_count =
#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED
      hdiv_history_ready ? 1 : GAMERA_NO_YINYANG_TIME_LEVELS;
#else
      GAMERA_NO_YINYANG_TIME_LEVELS;
#endif
#ifdef GAMERA_YINYANG_SPARSE_OVERSET
  if (donor_plan == NULL || !donor_plan->ready ||
      donor_plan->channel_count != level_count * GAMERA_NO_DIM ||
      donor_plan->reference_count != active_face_receptor_count * 8U) {
    return -1;
  }
#endif
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t face_count =
        gamera_no_element_count3(grid->face[direction].extent);
    memset(hdiv_target_face_flux[0][direction], 0,
           face_count * sizeof(*hdiv_target_face_flux[0][direction]));
    if (level_count > 1) {
      memset(hdiv_target_face_flux[1][direction], 0,
             face_count * sizeof(*hdiv_target_face_flux[1][direction]));
    }
  }
  for (size_t item = 0; item < active_face_receptor_count; ++item) {
    const face_receptor_t *receptor = &active_face_receptors[item];
    gamera_no_vec3 magnetic[GAMERA_NO_YINYANG_TIME_LEVELS] = {{{0.0}}};
    for (int time_level = 0; time_level < level_count; ++time_level) {
      for (int component = 0; component < GAMERA_NO_DIM; ++component) {
        const size_t channel =
            (size_t)time_level * GAMERA_NO_DIM + (size_t)component;
        for (int corner = 0; corner < 8; ++corner) {
#ifdef GAMERA_YINYANG_SPARSE_OVERSET
          const double donor_value = gamera_no_sparse_plan_value(
              donor_plan, item * 8U + (size_t)corner, (int)channel);
#else
          const size_t offset =
              (size_t)receptor->donor_rank[corner] *
                  (size_t)magnetic_send_count +
              channel * local_active_cells + receptor->donor_cell[corner];
          const double donor_value = magnetic_gather_buffer[offset];
#endif
          magnetic[time_level].value[component] +=
              receptor->weight[corner] * donor_value;
        }
      }
    }
    const gamera_no_face_geometry *geometry =
        &grid->face[receptor->direction].value[receptor->receiver_face];
    for (int time_level = 0; time_level < level_count; ++time_level) {
      const double donor_flux =
          vector_dot(magnetic[time_level], geometry->area_vector);
      const double self_flux =
          time_level == 0
              ? storage->face_flux[receptor->direction]
                                  [receptor->receiver_face]
              : storage->old_face_flux[receptor->direction]
                                      [receptor->receiver_face];
      const double increment =
          receptor->blend_weight * (donor_flux - self_flux);
      if (!isfinite(increment)) {
        return -1;
      }
      hdiv_target_face_flux[time_level][receptor->direction]
                           [receptor->receiver_face] = increment;
    }
  }
  return 0;
}

static int pack_hdiv_face_targets(const gamera_no_storage *storage,
                                  const gamera_no_grid *grid) {
  if (storage == NULL || grid == NULL) {
    return -1;
  }
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
#pragma omp parallel for collapse(3) schedule(static)
#endif
  for (int i = 0; i < config.ni; ++i) {
    for (int j = 0; j < config.nj; ++j) {
      for (int k = 0; k < config.nk; ++k) {
        const size_t local = local_active_index(i, j, k);
        const size_t lower[3] = {(size_t)(is + i), (size_t)(js + j),
                                 (size_t)(ks + k)};
#ifndef GAMERA_YINYANG_HDIV_OPTIMIZED
        const int global[3] = {proc_coords[0] * config.ni + i,
                               proc_coords[1] * config.nj + j,
                               proc_coords[2] * config.nk + k};
#endif
        for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
          for (int side = 0; side < GAMERA_NO_HDIV_FACE_SIDES; ++side) {
            size_t coordinate[3] = {lower[0], lower[1], lower[2]};
            coordinate[direction] += (size_t)side;
            const size_t face = gamera_no_index3(
                grid->face[direction].extent, coordinate[0], coordinate[1],
                coordinate[2]);
            const int level_count =
#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED
                hdiv_history_ready ? 1 : GAMERA_NO_YINYANG_TIME_LEVELS;
#else
                GAMERA_NO_YINYANG_TIME_LEVELS;
#endif
            for (int time_level = 0; time_level < level_count;
                 ++time_level) {
              const size_t channel =
                  hdiv_target_channel(time_level, direction, side);
              hdiv_send_buffer[channel * local_active_cells + local] =
                  hdiv_target_face_flux[time_level][direction][face];
              const size_t base_channel =
                  hdiv_base_channel(time_level, direction, side);
              hdiv_send_buffer[base_channel * local_active_cells + local] =
                  time_level == 0 ? storage->face_flux[direction][face]
                                  : storage->old_face_flux[direction][face];
            }
#ifndef GAMERA_YINYANG_HDIV_OPTIMIZED
            double metric_weight = 0.0;
            if (hdiv_metric_face_weight(grid, direction, side, lower, global,
                                        &metric_weight) != 0) {
              return -1;
            }
            hdiv_send_buffer[hdiv_weight_channel(direction, side) *
                                 local_active_cells +
                                 local] = metric_weight;
#endif
          }
        }
      }
    }
  }
  return 0;
}

#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED
static int hdiv_patch_consensus_failure(int local_failed) {
  int patch_failed = 0;
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  const double reduction_start = omp_get_wtime();
#endif
  if (MPI_Allreduce(&local_failed, &patch_failed, 1, MPI_INT, MPI_MAX,
                    comm_cart) != MPI_SUCCESS) {
    return -1;
  }
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  hdiv_profile_detail_current[GAMERA_NO_HDIV_DETAIL_REDUCTION] +=
      omp_get_wtime() - reduction_start;
  ++hdiv_profile_patch_allreduce_current;
#endif
  return patch_failed ? -1 : 0;
}

static size_t hdiv_pack_face_slab(const double *field, int direction,
                                  int axis, size_t begin, size_t width,
                                  double *buffer) {
  size_t packed = 0;
  const size_t *extent = hdiv_local_face_extent[direction];
  size_t slab_begin[3] = {0, 0, 0};
  size_t slab_end[3] = {extent[0], extent[1], extent[2]};
  slab_begin[axis] = begin;
  slab_end[axis] = begin + width;
  for (size_t i = slab_begin[0]; i < slab_end[0]; ++i) {
    for (size_t j = slab_begin[1]; j < slab_end[1]; ++j) {
      for (size_t k = slab_begin[2]; k < slab_end[2]; ++k) {
        buffer[packed++] =
            field[hdiv_local_face_index(direction, i, j, k)];
      }
    }
  }
  return packed;
}

static size_t hdiv_unpack_face_slab(double *field, int direction, int axis,
                                    size_t begin, size_t width,
                                    const double *buffer) {
  size_t unpacked = 0;
  const size_t *extent = hdiv_local_face_extent[direction];
  size_t slab_begin[3] = {0, 0, 0};
  size_t slab_end[3] = {extent[0], extent[1], extent[2]};
  slab_begin[axis] = begin;
  slab_end[axis] = begin + width;
  for (size_t i = slab_begin[0]; i < slab_end[0]; ++i) {
    for (size_t j = slab_begin[1]; j < slab_end[1]; ++j) {
      for (size_t k = slab_begin[2]; k < slab_end[2]; ++k) {
        field[hdiv_local_face_index(direction, i, j, k)] =
            buffer[unpacked++];
      }
    }
  }
  return unpacked;
}

/*
 * Synchronize canonical active faces and every in-patch magnetic ghost.
 * Sending the lower active slab toward nbr_low includes the shared active
 * face when the exchange axis is the face normal.  The opposite transfer
 * excludes that duplicate and fills only lower ghosts.
 */
static int hdiv_verify_static_layout_collective(void) {
  const int local_cells[3] = {config.ni, config.nj, config.nk};
  const int local_lower[3] = {is, js, ks};
  int local_failed = 0;
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
      const size_t lower = (size_t)local_lower[axis];
      const size_t cells = (size_t)local_cells[axis];
      const size_t active = cells + (axis == direction ? 1U : 0U);
      const size_t extent = hdiv_local_face_extent[direction][axis];
      const size_t upper_ghost =
          lower + active <= extent ? extent - lower - active : 0;
      const size_t width_to_low =
          upper_ghost + (axis == direction ? 1U : 0U);
      const size_t width_to_high = lower;
      size_t plane = 1;
      for (int transverse = 0; transverse < GAMERA_NO_DIM; ++transverse) {
        if (transverse != axis) {
          plane *= hdiv_local_face_extent[direction][transverse];
        }
      }
      local_failed |= lower + active > extent || width_to_low > active ||
                      width_to_high > cells ||
                      plane * width_to_low > hdiv_face_halo_capacity ||
                      plane * width_to_high > hdiv_face_halo_capacity ||
                      plane * width_to_low > INT32_MAX ||
                      plane * width_to_high > INT32_MAX;

      /*
       * Count matching is topology, not iteration state.  Check each peer
       * pairing once so the hot halo path can execute a fixed Sendrecv
       * schedule without a preflight collective.
       */
      const unsigned long long count_to_low =
          (unsigned long long)(plane * width_to_low);
      const unsigned long long count_to_high =
          (unsigned long long)(plane * width_to_high);
      unsigned long long count_from_high = count_to_low;
      unsigned long long count_from_low = count_to_high;
      const int tag_base = 7500 + direction * 2 * GAMERA_NO_DIM + 2 * axis;
      local_failed |=
          MPI_Sendrecv(&count_to_low, 1, MPI_UNSIGNED_LONG_LONG,
                       nbr_low[axis], tag_base, &count_from_high, 1,
                       MPI_UNSIGNED_LONG_LONG, nbr_high[axis], tag_base,
                       comm_cart, MPI_STATUS_IGNORE) != MPI_SUCCESS;
#ifdef GAMERA_YINYANG_HDIV_PROFILE
      ++hdiv_profile_layout_sendrecv_current;
#endif
      local_failed |=
          MPI_Sendrecv(&count_to_high, 1, MPI_UNSIGNED_LONG_LONG,
                       nbr_high[axis], tag_base + 1, &count_from_low, 1,
                       MPI_UNSIGNED_LONG_LONG, nbr_low[axis], tag_base + 1,
                       comm_cart, MPI_STATUS_IGNORE) != MPI_SUCCESS;
#ifdef GAMERA_YINYANG_HDIV_PROFILE
      ++hdiv_profile_layout_sendrecv_current;
#endif
      if ((nbr_high[axis] != MPI_PROC_NULL &&
           count_from_high != count_to_low) ||
          (nbr_low[axis] != MPI_PROC_NULL &&
           count_from_low != count_to_high)) {
        local_failed = 1;
      }
    }
  }
  for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
    size_t plane = 1;
    for (int transverse = 0; transverse < GAMERA_NO_DIM; ++transverse) {
      if (transverse != axis) {
        plane *= (size_t)local_cells[transverse];
      }
    }
    local_failed |= local_cells[axis] <= 0 || local_lower[axis] < 1 ||
                    (size_t)local_lower[axis] +
                            (size_t)local_cells[axis] >=
                        hdiv_local_cell_extent[axis] ||
                    plane > hdiv_cell_halo_capacity || plane > INT32_MAX;
    const unsigned long long plane_count = (unsigned long long)plane;
    unsigned long long count_from_high = plane_count;
    unsigned long long count_from_low = plane_count;
    const int tag_base = 7540 + 2 * axis;
    local_failed |=
        MPI_Sendrecv(&plane_count, 1, MPI_UNSIGNED_LONG_LONG, nbr_low[axis],
                     tag_base, &count_from_high, 1, MPI_UNSIGNED_LONG_LONG,
                     nbr_high[axis], tag_base, comm_cart,
                     MPI_STATUS_IGNORE) != MPI_SUCCESS;
#ifdef GAMERA_YINYANG_HDIV_PROFILE
    ++hdiv_profile_layout_sendrecv_current;
#endif
    local_failed |=
        MPI_Sendrecv(&plane_count, 1, MPI_UNSIGNED_LONG_LONG, nbr_high[axis],
                     tag_base + 1, &count_from_low, 1,
                     MPI_UNSIGNED_LONG_LONG, nbr_low[axis], tag_base + 1,
                     comm_cart, MPI_STATUS_IGNORE) != MPI_SUCCESS;
#ifdef GAMERA_YINYANG_HDIV_PROFILE
    ++hdiv_profile_layout_sendrecv_current;
#endif
    if ((nbr_high[axis] != MPI_PROC_NULL &&
         count_from_high != plane_count) ||
        (nbr_low[axis] != MPI_PROC_NULL && count_from_low != plane_count)) {
      local_failed = 1;
    }
  }
  if (hdiv_patch_consensus_failure(local_failed) != 0) {
    return -1;
  }
  hdiv_distributed_layout_verified = 1;
  return 0;
}

static int hdiv_exchange_face_field(double *field, int direction) {
  const int local_cells[3] = {config.ni, config.nj, config.nk};
  const int local_lower[3] = {is, js, ks};
  if (!hdiv_distributed_layout_verified || field == NULL || direction < 0 ||
      direction >= GAMERA_NO_DIM) {
    return -1;
  }
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  const double halo_start = omp_get_wtime();
#endif
  int local_failed = 0;
  for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
    const size_t lower = (size_t)local_lower[axis];
    const size_t cells = (size_t)local_cells[axis];
    const size_t active = cells + (axis == direction ? 1U : 0U);
    const size_t extent = hdiv_local_face_extent[direction][axis];
    const size_t upper_ghost =
        lower + active <= extent ? extent - lower - active : 0;
    const size_t width_to_low =
        upper_ghost + (axis == direction ? 1U : 0U);
    const size_t width_to_high = lower;
    const int tag_base = 7600 + direction * 2 * GAMERA_NO_DIM + 2 * axis;

    if (nbr_low[axis] == MPI_PROC_NULL &&
        nbr_high[axis] == MPI_PROC_NULL) {
      continue;
    }

    if (width_to_low > 0) {
      size_t expected_count = width_to_low;
      for (int transverse = 0; transverse < GAMERA_NO_DIM; ++transverse) {
        if (transverse != axis) {
          expected_count *= hdiv_local_face_extent[direction][transverse];
        }
      }
      const size_t send_count = hdiv_pack_face_slab(
          field, direction, axis, lower, width_to_low, hdiv_face_halo_send);
      local_failed |= send_count != expected_count;
      const int status = MPI_Sendrecv(
          hdiv_face_halo_send, (int)expected_count, MPI_DOUBLE, nbr_low[axis],
          tag_base, hdiv_face_halo_receive, (int)expected_count, MPI_DOUBLE,
          nbr_high[axis], tag_base, comm_cart, MPI_STATUS_IGNORE);
#ifdef GAMERA_YINYANG_HDIV_PROFILE
      ++hdiv_profile_face_sendrecv_current;
      if (nbr_low[axis] != MPI_PROC_NULL) {
        hdiv_profile_face_send_bytes_current +=
            (unsigned long long)expected_count * sizeof(double);
      }
#endif
      if (status != MPI_SUCCESS) {
        local_failed = 1;
      } else if (nbr_high[axis] != MPI_PROC_NULL &&
                 hdiv_unpack_face_slab(field, direction, axis,
                                       lower + cells, width_to_low,
                                       hdiv_face_halo_receive) !=
                     expected_count) {
        local_failed = 1;
      }
    }

    if (width_to_high > 0) {
      size_t expected_count = width_to_high;
      for (int transverse = 0; transverse < GAMERA_NO_DIM; ++transverse) {
        if (transverse != axis) {
          expected_count *= hdiv_local_face_extent[direction][transverse];
        }
      }
      const size_t send_begin = lower + cells - width_to_high;
      const size_t send_count = hdiv_pack_face_slab(
          field, direction, axis, send_begin, width_to_high,
          hdiv_face_halo_send);
      local_failed |= send_count != expected_count;
      const int status = MPI_Sendrecv(
          hdiv_face_halo_send, (int)expected_count, MPI_DOUBLE,
          nbr_high[axis], tag_base + 1, hdiv_face_halo_receive,
          (int)expected_count, MPI_DOUBLE, nbr_low[axis], tag_base + 1,
          comm_cart, MPI_STATUS_IGNORE);
#ifdef GAMERA_YINYANG_HDIV_PROFILE
      ++hdiv_profile_face_sendrecv_current;
      if (nbr_high[axis] != MPI_PROC_NULL) {
        hdiv_profile_face_send_bytes_current +=
            (unsigned long long)expected_count * sizeof(double);
      }
#endif
      if (status != MPI_SUCCESS) {
        local_failed = 1;
      } else if (nbr_low[axis] != MPI_PROC_NULL &&
                 hdiv_unpack_face_slab(field, direction, axis,
                                       lower - width_to_high,
                                       width_to_high,
                                       hdiv_face_halo_receive) !=
                     expected_count) {
        local_failed = 1;
      }
    }
  }
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  hdiv_profile_detail_current[GAMERA_NO_HDIV_DETAIL_FACE_HALO] +=
      omp_get_wtime() - halo_start;
#endif
  return local_failed ? -1 : 0;
}

static size_t hdiv_pack_cell_plane(const double *field, int axis,
                                   size_t plane, double *buffer) {
  size_t packed = 0;
  const size_t lower[3] = {(size_t)is, (size_t)js, (size_t)ks};
  const size_t count[3] = {(size_t)config.ni, (size_t)config.nj,
                           (size_t)config.nk};
  size_t begin[3] = {lower[0], lower[1], lower[2]};
  size_t end[3] = {lower[0] + count[0], lower[1] + count[1],
                   lower[2] + count[2]};
  begin[axis] = plane;
  end[axis] = plane + 1U;
  for (size_t i = begin[0]; i < end[0]; ++i) {
    for (size_t j = begin[1]; j < end[1]; ++j) {
      for (size_t k = begin[2]; k < end[2]; ++k) {
        buffer[packed++] = field[hdiv_local_cell_index(i, j, k)];
      }
    }
  }
  return packed;
}

static size_t hdiv_unpack_cell_plane(double *field, int axis, size_t plane,
                                     const double *buffer) {
  size_t unpacked = 0;
  const size_t lower[3] = {(size_t)is, (size_t)js, (size_t)ks};
  const size_t count[3] = {(size_t)config.ni, (size_t)config.nj,
                           (size_t)config.nk};
  size_t begin[3] = {lower[0], lower[1], lower[2]};
  size_t end[3] = {lower[0] + count[0], lower[1] + count[1],
                   lower[2] + count[2]};
  begin[axis] = plane;
  end[axis] = plane + 1U;
  for (size_t i = begin[0]; i < end[0]; ++i) {
    for (size_t j = begin[1]; j < end[1]; ++j) {
      for (size_t k = begin[2]; k < end[2]; ++k) {
        field[hdiv_local_cell_index(i, j, k)] = buffer[unpacked++];
      }
    }
  }
  return unpacked;
}

static void hdiv_zero_cell_plane(double *field, int axis, size_t plane) {
  const size_t lower[3] = {(size_t)is, (size_t)js, (size_t)ks};
  const size_t count[3] = {(size_t)config.ni, (size_t)config.nj,
                           (size_t)config.nk};
  size_t begin[3] = {lower[0], lower[1], lower[2]};
  size_t end[3] = {lower[0] + count[0], lower[1] + count[1],
                   lower[2] + count[2]};
  begin[axis] = plane;
  end[axis] = plane + 1U;
  for (size_t i = begin[0]; i < end[0]; ++i) {
    for (size_t j = begin[1]; j < end[1]; ++j) {
      for (size_t k = begin[2]; k < end[2]; ++k) {
        field[hdiv_local_cell_index(i, j, k)] = 0.0;
      }
    }
  }
}

static int hdiv_exchange_cell_halo(double *field) {
  const size_t lower[3] = {(size_t)is, (size_t)js, (size_t)ks};
  const size_t count[3] = {(size_t)config.ni, (size_t)config.nj,
                           (size_t)config.nk};
  if (!hdiv_distributed_layout_verified || field == NULL) {
    return -1;
  }
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  const double halo_start = omp_get_wtime();
#endif
  int local_failed = 0;
  for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
    const int tag_base = 7700 + 2 * axis;
    size_t expected_count = 1;
    for (int transverse = 0; transverse < GAMERA_NO_DIM; ++transverse) {
      if (transverse != axis) {
        expected_count *= count[transverse];
      }
    }
    if (nbr_low[axis] == MPI_PROC_NULL &&
        nbr_high[axis] == MPI_PROC_NULL) {
      hdiv_zero_cell_plane(field, axis, lower[axis] - 1U);
      hdiv_zero_cell_plane(field, axis, lower[axis] + count[axis]);
      continue;
    }
    const size_t plane_count =
        hdiv_pack_cell_plane(field, axis, lower[axis], hdiv_cell_halo_send);
    local_failed |= plane_count != expected_count;
    const int lower_status = MPI_Sendrecv(
        hdiv_cell_halo_send, (int)expected_count, MPI_DOUBLE, nbr_low[axis],
        tag_base, hdiv_cell_halo_receive, (int)expected_count, MPI_DOUBLE,
        nbr_high[axis], tag_base, comm_cart, MPI_STATUS_IGNORE);
#ifdef GAMERA_YINYANG_HDIV_PROFILE
    ++hdiv_profile_cell_sendrecv_current;
    if (nbr_low[axis] != MPI_PROC_NULL) {
      hdiv_profile_cell_send_bytes_current +=
          (unsigned long long)expected_count * sizeof(double);
    }
#endif
    if (lower_status != MPI_SUCCESS) {
      local_failed = 1;
    } else if (nbr_high[axis] != MPI_PROC_NULL) {
      if (hdiv_unpack_cell_plane(field, axis, lower[axis] + count[axis],
                                 hdiv_cell_halo_receive) != expected_count) {
        local_failed = 1;
      }
    } else {
      hdiv_zero_cell_plane(field, axis, lower[axis] + count[axis]);
    }

    const size_t high_plane = lower[axis] + count[axis] - 1U;
    const size_t high_count =
        hdiv_pack_cell_plane(field, axis, high_plane, hdiv_cell_halo_send);
    local_failed |= high_count != expected_count;
    const int upper_status = MPI_Sendrecv(
        hdiv_cell_halo_send, (int)expected_count, MPI_DOUBLE, nbr_high[axis],
        tag_base + 1, hdiv_cell_halo_receive, (int)expected_count, MPI_DOUBLE,
        nbr_low[axis], tag_base + 1, comm_cart, MPI_STATUS_IGNORE);
#ifdef GAMERA_YINYANG_HDIV_PROFILE
    ++hdiv_profile_cell_sendrecv_current;
    if (nbr_high[axis] != MPI_PROC_NULL) {
      hdiv_profile_cell_send_bytes_current +=
          (unsigned long long)expected_count * sizeof(double);
    }
#endif
    if (upper_status != MPI_SUCCESS) {
      local_failed = 1;
    } else if (nbr_low[axis] != MPI_PROC_NULL) {
      if (hdiv_unpack_cell_plane(field, axis, lower[axis] - 1U,
                                 hdiv_cell_halo_receive) != expected_count) {
        local_failed = 1;
      }
    } else {
      hdiv_zero_cell_plane(field, axis, lower[axis] - 1U);
    }
  }
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  hdiv_profile_detail_current[GAMERA_NO_HDIV_DETAIL_CELL_HALO] +=
      omp_get_wtime() - halo_start;
#endif
  return local_failed ? -1 : 0;
}

static double hdiv_local_cell_net_flux(int time_level, int i, int j, int k) {
  size_t lower[3] = {(size_t)(is + i), (size_t)(js + j),
                     (size_t)(ks + k)};
  double net = 0.0;
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    size_t upper[3] = {lower[0], lower[1], lower[2]};
    ++upper[direction];
    net += hdiv_local_face_flux[time_level][direction]
                               [hdiv_local_face_index(
                                   direction, upper[0], upper[1], upper[2])] -
           hdiv_local_face_flux[time_level][direction]
                               [hdiv_local_face_index(
                                   direction, lower[0], lower[1], lower[2])];
  }
  return net;
}

static double hdiv_local_cell_diagonal(int i, int j, int k) {
  size_t lower[3] = {(size_t)(is + i), (size_t)(js + j),
                     (size_t)(ks + k)};
  double diagonal = 0.0;
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    size_t upper[3] = {lower[0], lower[1], lower[2]};
    ++upper[direction];
    diagonal += hdiv_local_face_mobility[direction][hdiv_local_face_index(
        direction, lower[0], lower[1], lower[2])];
    diagonal += hdiv_local_face_mobility[direction][hdiv_local_face_index(
        direction, upper[0], upper[1], upper[2])];
  }
  return diagonal;
}

static int assemble_hdiv_local_faces(int include_metric_weights) {
  const int include_old = !hdiv_history_ready;
  const int global_count[3] = {config.ni_global, config.nj_global,
                               config.nk_global};
  if (!hdiv_distributed_layout_verified &&
      hdiv_verify_static_layout_collective() != 0) {
    return -1;
  }
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t face_count =
        gamera_no_element_count3(hdiv_local_face_extent[direction]);
    memset(hdiv_local_face_base[0][direction], 0,
           face_count * sizeof(*hdiv_local_face_base[0][direction]));
    memset(hdiv_local_face_flux[0][direction], 0,
           face_count * sizeof(*hdiv_local_face_flux[0][direction]));
    if (include_old) {
      memset(hdiv_local_face_base[1][direction], 0,
             face_count * sizeof(*hdiv_local_face_base[1][direction]));
      memset(hdiv_local_face_flux[1][direction], 0,
             face_count * sizeof(*hdiv_local_face_flux[1][direction]));
    }
    if (include_metric_weights) {
      memset(hdiv_local_face_mobility[direction], 0,
             face_count * sizeof(*hdiv_local_face_mobility[direction]));
    }
  }

#pragma omp parallel for collapse(3) schedule(static)
  for (int i = 0; i < config.ni; ++i) {
    for (int j = 0; j < config.nj; ++j) {
      for (int k = 0; k < config.nk; ++k) {
        const size_t local = local_active_index(i, j, k);
        const size_t cell_lower[3] = {(size_t)(is + i), (size_t)(js + j),
                                      (size_t)(ks + k)};
        const int global[3] = {proc_coords[0] * config.ni + i,
                               proc_coords[1] * config.nj + j,
                               proc_coords[2] * config.nk + k};
        for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
          size_t face_coordinate[3] = {cell_lower[0], cell_lower[1],
                                       cell_lower[2]};
          size_t face = hdiv_local_face_index(
              direction, face_coordinate[0], face_coordinate[1],
              face_coordinate[2]);
          for (int time_level = 0; time_level < (include_old ? 2 : 1);
               ++time_level) {
            hdiv_local_face_base[time_level][direction][face] =
                hdiv_send_buffer[hdiv_base_channel(
                                     time_level, direction, GAMERA_NO_LOWER) *
                                     local_active_cells +
                                 local];
            hdiv_local_face_flux[time_level][direction][face] =
                hdiv_send_buffer[hdiv_target_channel(
                                     time_level, direction, GAMERA_NO_LOWER) *
                                     local_active_cells +
                                 local];
          }
          if (include_metric_weights) {
            hdiv_local_face_mobility[direction][face] =
                hdiv_send_buffer[hdiv_weight_channel(
                                     direction, GAMERA_NO_LOWER) *
                                     local_active_cells +
                                 local];
          }

          if (global[direction] == global_count[direction] - 1) {
            ++face_coordinate[direction];
            face = hdiv_local_face_index(
                direction, face_coordinate[0], face_coordinate[1],
                face_coordinate[2]);
            for (int time_level = 0; time_level < (include_old ? 2 : 1);
                 ++time_level) {
              hdiv_local_face_base[time_level][direction][face] =
                  hdiv_send_buffer[hdiv_base_channel(
                                       time_level, direction,
                                       GAMERA_NO_UPPER) *
                                       local_active_cells +
                                   local];
              hdiv_local_face_flux[time_level][direction][face] =
                  hdiv_send_buffer[hdiv_target_channel(
                                       time_level, direction,
                                       GAMERA_NO_UPPER) *
                                       local_active_cells +
                                   local];
            }
            if (include_metric_weights) {
              hdiv_local_face_mobility[direction][face] =
                  hdiv_send_buffer[hdiv_weight_channel(
                                       direction, GAMERA_NO_UPPER) *
                                       local_active_cells +
                                   local];
            }
          }
        }
      }
    }
  }

  int local_failed = 0;
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    local_failed |=
        hdiv_exchange_face_field(hdiv_local_face_base[0][direction],
                                 direction) != 0;
    local_failed |=
        hdiv_exchange_face_field(hdiv_local_face_flux[0][direction],
                                 direction) != 0;
    if (include_old) {
      local_failed |=
          hdiv_exchange_face_field(hdiv_local_face_base[1][direction],
                                   direction) != 0;
      local_failed |=
          hdiv_exchange_face_field(hdiv_local_face_flux[1][direction],
                                   direction) != 0;
    }
    if (include_metric_weights) {
      local_failed |=
          hdiv_exchange_face_field(hdiv_local_face_mobility[direction],
                                   direction) != 0;
    }
  }
  if (hdiv_patch_consensus_failure(local_failed) != 0) {
    return -1;
  }

  double local_min_weight = DBL_MAX;
  double local_max_weight = 0.0;
  for (int i = 0; i < config.ni; ++i) {
    for (int j = 0; j < config.nj; ++j) {
      for (int k = 0; k < config.nk; ++k) {
        const size_t cell_lower[3] = {(size_t)(is + i), (size_t)(js + j),
                                      (size_t)(ks + k)};
        const int global[3] = {proc_coords[0] * config.ni + i,
                               proc_coords[1] * config.nj + j,
                               proc_coords[2] * config.nk + k};
        for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
          const size_t face = hdiv_local_face_index(
              direction, cell_lower[0], cell_lower[1], cell_lower[2]);
          const double weight = hdiv_local_face_mobility[direction][face];
          const int locked_radial_face =
              direction == GAMERA_NO_I && global[GAMERA_NO_I] == 0;
          if (include_metric_weights &&
              (!isfinite(weight) || weight < 0.0 ||
               (locked_radial_face ? weight != 0.0 : weight == 0.0))) {
            local_failed = 1;
          }
          if (locked_radial_face &&
              (hdiv_local_face_flux[0][direction][face] != 0.0 ||
               (include_old &&
                hdiv_local_face_flux[1][direction][face] != 0.0))) {
            local_failed = 1;
          }
          if (include_metric_weights && weight > 0.0) {
            local_min_weight = fmin(local_min_weight, weight);
            local_max_weight = fmax(local_max_weight, weight);
          }
          if (global[direction] == global_count[direction] - 1) {
            size_t upper[3] = {cell_lower[0], cell_lower[1], cell_lower[2]};
            ++upper[direction];
            const size_t upper_face = hdiv_local_face_index(
                direction, upper[0], upper[1], upper[2]);
            const double upper_weight =
                hdiv_local_face_mobility[direction][upper_face];
            const int upper_locked =
                direction == GAMERA_NO_I &&
                global_count[GAMERA_NO_I] == global[GAMERA_NO_I] + 1;
            if (include_metric_weights &&
                (!isfinite(upper_weight) || upper_weight < 0.0 ||
                 (upper_locked ? upper_weight != 0.0
                               : upper_weight == 0.0))) {
              local_failed = 1;
            }
            if (upper_locked &&
                (hdiv_local_face_flux[0][direction][upper_face] != 0.0 ||
                 (include_old &&
                  hdiv_local_face_flux[1][direction][upper_face] != 0.0))) {
              local_failed = 1;
            }
            if (include_metric_weights && upper_weight > 0.0) {
              local_min_weight = fmin(local_min_weight, upper_weight);
              local_max_weight = fmax(local_max_weight, upper_weight);
            }
          }
        }
      }
    }
  }
  if (hdiv_patch_consensus_failure(local_failed) != 0) {
    return -1;
  }

  if (include_metric_weights) {
    double patch_min_weight = DBL_MAX;
    double patch_max_weight = 0.0;
    double patch_min_cosine = 1.0;
#ifdef GAMERA_YINYANG_HDIV_PROFILE
    const double reduction_start = omp_get_wtime();
#endif
    local_failed =
        MPI_Allreduce(&local_min_weight, &patch_min_weight, 1, MPI_DOUBLE,
                      MPI_MIN, comm_cart) != MPI_SUCCESS;
    local_failed |=
        MPI_Allreduce(&local_max_weight, &patch_max_weight, 1, MPI_DOUBLE,
                      MPI_MAX, comm_cart) != MPI_SUCCESS;
    local_failed |=
        MPI_Allreduce(&gamera_no_yinyang_hdiv_min_metric_cosine,
                      &patch_min_cosine, 1, MPI_DOUBLE, MPI_MIN,
                      comm_cart) != MPI_SUCCESS;
#ifdef GAMERA_YINYANG_HDIV_PROFILE
    hdiv_profile_detail_current[GAMERA_NO_HDIV_DETAIL_REDUCTION] +=
        omp_get_wtime() - reduction_start;
    hdiv_profile_patch_allreduce_current += 3;
#endif
    if (hdiv_patch_consensus_failure(local_failed) != 0) {
      return -1;
    }
    gamera_no_yinyang_hdiv_min_weight =
        fmin(gamera_no_yinyang_hdiv_min_weight, patch_min_weight);
    gamera_no_yinyang_hdiv_max_weight =
        fmax(gamera_no_yinyang_hdiv_max_weight, patch_max_weight);
    gamera_no_yinyang_hdiv_min_metric_cosine = patch_min_cosine;

#pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < config.ni; ++i) {
      for (int j = 0; j < config.nj; ++j) {
        for (int k = 0; k < config.nk; ++k) {
          hdiv_cell_diagonal_cache[hdiv_active_cell_index(i, j, k)] =
              hdiv_local_cell_diagonal(i, j, k);
        }
      }
    }
    hdiv_cell_diagonal_cache_ready = 1;
  }
  return 0;
}

static int apply_hdiv_local_operator(double *input, double *output) {
  const int halo_failed = hdiv_exchange_cell_halo(input) != 0;
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  const double operator_start = omp_get_wtime();
#endif
#pragma omp parallel for collapse(3) schedule(static)
  for (int i = 0; i < config.ni; ++i) {
    for (int j = 0; j < config.nj; ++j) {
      for (int k = 0; k < config.nk; ++k) {
        const size_t coordinate[3] = {(size_t)(is + i), (size_t)(js + j),
                                      (size_t)(ks + k)};
        const size_t cell = hdiv_active_cell_index(i, j, k);
        double value = 0.0;
        for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
          size_t lower_cell[3] = {coordinate[0], coordinate[1],
                                  coordinate[2]};
          size_t upper_cell[3] = {coordinate[0], coordinate[1],
                                  coordinate[2]};
          --lower_cell[direction];
          ++upper_cell[direction];
          size_t upper_face[3] = {coordinate[0], coordinate[1],
                                  coordinate[2]};
          ++upper_face[direction];
          const double lower_value = input[hdiv_local_cell_index(
              lower_cell[0], lower_cell[1], lower_cell[2])];
          const double upper_value = input[hdiv_local_cell_index(
              upper_cell[0], upper_cell[1], upper_cell[2])];
          const double lower_weight =
              hdiv_local_face_mobility[direction][hdiv_local_face_index(
                  direction, coordinate[0], coordinate[1], coordinate[2])];
          const double upper_weight =
              hdiv_local_face_mobility[direction][hdiv_local_face_index(
                  direction, upper_face[0], upper_face[1], upper_face[2])];
          value += lower_weight * (input[cell] - lower_value);
          value += upper_weight * (input[cell] - upper_value);
        }
        output[cell] = value;
      }
    }
  }
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  hdiv_profile_detail_current[GAMERA_NO_HDIV_DETAIL_OPERATOR] +=
      omp_get_wtime() - operator_start;
#endif
  return halo_failed ? -1 : 0;
}

static double hdiv_local_dot(const double *left, const double *right) {
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  const double vector_start = omp_get_wtime();
#endif
#pragma omp parallel for schedule(static)
  for (int i = 0; i < config.ni; ++i) {
    double partial = 0.0;
    for (int j = 0; j < config.nj; ++j) {
      for (int k = 0; k < config.nk; ++k) {
        const size_t cell = hdiv_active_cell_index(i, j, k);
        partial += left[cell] * right[cell];
      }
    }
    hdiv_dot_block_partial[i] = partial;
  }
  double sum = 0.0;
  for (int i = 0; i < config.ni; ++i) {
    sum += hdiv_dot_block_partial[i];
  }
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  hdiv_profile_detail_current[GAMERA_NO_HDIV_DETAIL_VECTOR] +=
      omp_get_wtime() - vector_start;
#endif
  return sum;
}

static int hdiv_patch_sum(double local, double *global) {
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  const double reduction_start = omp_get_wtime();
#endif
  const int status =
      MPI_Allreduce(&local, global, 1, MPI_DOUBLE, MPI_SUM, comm_cart);
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  hdiv_profile_detail_current[GAMERA_NO_HDIV_DETAIL_REDUCTION] +=
      omp_get_wtime() - reduction_start;
  ++hdiv_profile_patch_allreduce_current;
  ++hdiv_profile_math_allreduce_current;
#endif
  return status == MPI_SUCCESS ? 0 : -1;
}

static int hdiv_patch_max(double local, double *global) {
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  const double reduction_start = omp_get_wtime();
#endif
  const int status =
      MPI_Allreduce(&local, global, 1, MPI_DOUBLE, MPI_MAX, comm_cart);
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  hdiv_profile_detail_current[GAMERA_NO_HDIV_DETAIL_REDUCTION] +=
      omp_get_wtime() - reduction_start;
  ++hdiv_profile_patch_allreduce_current;
  ++hdiv_profile_math_allreduce_current;
#endif
  return status == MPI_SUCCESS ? 0 : -1;
}

/* Piggyback a locally observed halo/MPI failure on an existing mathematical
 * reduction.  The first element is the unchanged PCG scalar; the second is
 * only an aligned failure gate and does not alter the Hodge algebra. */
static int hdiv_patch_sum_with_failure(double local, int local_failed,
                                       double *global) {
  const double input[2] = {local, local_failed ? 1.0 : 0.0};
  double output[2] = {0.0, 0.0};
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  const double reduction_start = omp_get_wtime();
#endif
  const int status =
      MPI_Allreduce(input, output, 2, MPI_DOUBLE, MPI_SUM, comm_cart);
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  hdiv_profile_detail_current[GAMERA_NO_HDIV_DETAIL_REDUCTION] +=
      omp_get_wtime() - reduction_start;
  ++hdiv_profile_patch_allreduce_current;
  ++hdiv_profile_math_allreduce_current;
#endif
  if (status != MPI_SUCCESS) {
    return -1;
  }
  *global = output[0];
  return output[1] == 0.0 ? 0 : -1;
}

static int hdiv_patch_max_with_failure(double local, int local_failed,
                                       double *global) {
  const double input[2] = {local, local_failed ? 1.0 : 0.0};
  double output[2] = {0.0, 0.0};
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  const double reduction_start = omp_get_wtime();
#endif
  const int status =
      MPI_Allreduce(input, output, 2, MPI_DOUBLE, MPI_MAX, comm_cart);
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  hdiv_profile_detail_current[GAMERA_NO_HDIV_DETAIL_REDUCTION] +=
      omp_get_wtime() - reduction_start;
  ++hdiv_profile_patch_allreduce_current;
  ++hdiv_profile_math_allreduce_current;
#endif
  if (status != MPI_SUCCESS) {
    return -1;
  }
  *global = output[0];
  return output[1] == 0.0 ? 0 : -1;
}

#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED_VERIFY
static int hdiv_verify_same_bits(double left, double right) {
  return memcmp(&left, &right, sizeof(left)) == 0;
}

static int hdiv_verify_patch_rank_and_local(const int global[3],
                                            int *owner_patch_rank,
                                            size_t *owner_local) {
  const int local_shape[3] = {config.ni, config.nj, config.nk};
  const int global_shape[3] = {config.ni_global, config.nj_global,
                               config.nk_global};
  int coords[3];
  int local[3];
  for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
    if (global[axis] < 0 || global[axis] >= global_shape[axis] ||
        local_shape[axis] <= 0) {
      return -1;
    }
    coords[axis] = global[axis] / local_shape[axis];
    local[axis] = global[axis] % local_shape[axis];
    if (coords[axis] < 0 || coords[axis] >= config.proc_dims[axis]) {
      return -1;
    }
  }
  if (MPI_Cart_rank(comm_cart, coords, owner_patch_rank) != MPI_SUCCESS) {
    return -1;
  }
  *owner_local = local_active_index(local[0], local[1], local[2]);
  return 0;
}

static int hdiv_verify_gather_cell_channel(const double *gather,
                                           size_t rank_stride,
                                           const int global_cell[3],
                                           size_t channel, double *value) {
  int owner_patch_rank = -1;
  size_t owner_local = 0;
  if (gather == NULL || value == NULL ||
      hdiv_verify_patch_rank_and_local(global_cell, &owner_patch_rank,
                                       &owner_local) != 0) {
    return -1;
  }
  *value = gather[(size_t)owner_patch_rank * rank_stride +
                  channel * local_active_cells + owner_local];
  return 0;
}

static int hdiv_verify_gather_face_channel(const double *gather,
                                           size_t rank_stride, int direction,
                                           const int global_face[3],
                                           size_t lower_channel,
                                           size_t upper_channel,
                                           double *value) {
  const int global_shape[3] = {config.ni_global, config.nj_global,
                               config.nk_global};
  int owner_cell[3] = {global_face[0], global_face[1], global_face[2]};
  size_t channel = lower_channel;
  for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
    const int upper = global_shape[axis] + (axis == direction ? 1 : 0);
    if (global_face[axis] < 0 || global_face[axis] >= upper) {
      return -1;
    }
  }
  if (global_face[direction] == global_shape[direction]) {
    owner_cell[direction] = global_shape[direction] - 1;
    channel = upper_channel;
  }
  return hdiv_verify_gather_cell_channel(gather, rank_stride, owner_cell,
                                         channel, value);
}

static int hdiv_verify_snapshot_face_value(const double *gather,
                                           int time_level, int direction,
                                           const int global_face[3],
                                           int family, double *value) {
  size_t lower_channel = 0;
  size_t upper_channel = 0;
  if (family == 0) {
    lower_channel =
        hdiv_target_channel(time_level, direction, GAMERA_NO_LOWER);
    upper_channel =
        hdiv_target_channel(time_level, direction, GAMERA_NO_UPPER);
  } else if (family == 1) {
    lower_channel = hdiv_base_channel(time_level, direction, GAMERA_NO_LOWER);
    upper_channel = hdiv_base_channel(time_level, direction, GAMERA_NO_UPPER);
  } else {
    lower_channel = hdiv_weight_channel(direction, GAMERA_NO_LOWER);
    upper_channel = hdiv_weight_channel(direction, GAMERA_NO_UPPER);
  }
  return hdiv_verify_gather_face_channel(
      gather, (size_t)hdiv_send_count, direction, global_face, lower_channel,
      upper_channel, value);
}

static double hdiv_verify_lambda_formula(const int global_cell[3]) {
  return ldexp((double)(global_cell[0] + 1), -2) +
         ldexp((double)(global_cell[1] + 1), -6) +
         ldexp((double)(global_cell[2] + 1), -10);
}

static int hdiv_verify_global_lambda(const double *lambda_gather,
                                     const int global_cell[3],
                                     double *value) {
  if (global_cell[0] < 0 || global_cell[0] >= config.ni_global ||
      global_cell[1] < 0 || global_cell[1] >= config.nj_global ||
      global_cell[2] < 0 || global_cell[2] >= config.nk_global) {
    *value = 0.0;
    return 0;
  }
  return hdiv_verify_gather_cell_channel(lambda_gather, local_active_cells,
                                         global_cell, 0, value);
}

static int hdiv_verify_sentinel_ownership(void) {
  const size_t gather_count = (size_t)patch_size * (size_t)hdiv_send_count;
  double *sentinel_local =
      (double *)malloc((size_t)hdiv_send_count * sizeof(*sentinel_local));
  double *sentinel_gather =
      (double *)malloc(gather_count * sizeof(*sentinel_gather));
  size_t maximum_face_count = 0;
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t face_count =
        gamera_no_element_count3(hdiv_local_face_extent[direction]);
    maximum_face_count = maximum_face_count > face_count ? maximum_face_count
                                                         : face_count;
  }
  double *face = (double *)malloc(maximum_face_count * sizeof(*face));
  int local_failed = sentinel_local == NULL || sentinel_gather == NULL ||
                     face == NULL;
  if (hdiv_patch_consensus_failure(local_failed) != 0) {
    free(sentinel_local);
    free(sentinel_gather);
    free(face);
    return -1;
  }

  for (int channel = 0; channel < GAMERA_NO_HDIV_CHANNELS; ++channel) {
    for (size_t local = 0; local < local_active_cells; ++local) {
      /* Every wire owner/channel gets a distinct, exactly representable bit
       * pattern; in particular negative-rank upper and positive-rank lower
       * copies of an MPI interface can never alias accidentally. */
      sentinel_local[(size_t)channel * local_active_cells + local] =
          1.0 + (double)((size_t)patch_rank *
                             (size_t)GAMERA_NO_HDIV_CHANNELS *
                             local_active_cells +
                         (size_t)channel * local_active_cells + local);
    }
  }
  local_failed =
      MPI_Allgather(sentinel_local, hdiv_send_count, MPI_DOUBLE,
                    sentinel_gather, hdiv_send_count, MPI_DOUBLE,
                    comm_cart) != MPI_SUCCESS;
  if (hdiv_patch_consensus_failure(local_failed) != 0) {
    free(sentinel_local);
    free(sentinel_gather);
    free(face);
    return -1;
  }

  unsigned long long local_conflicts = 0;
  const int global_shape[3] = {config.ni_global, config.nj_global,
                               config.nk_global};
  const int local_shape[3] = {config.ni, config.nj, config.nk};
  const int local_lower[3] = {is, js, ks};
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t face_count =
        gamera_no_element_count3(hdiv_local_face_extent[direction]);
    for (int family = 0; family < 5; ++family) {
      const int time_level = family & 1;
      size_t lower_channel = 0;
      size_t upper_channel = 0;
      if (family < 2) {
        lower_channel =
            hdiv_target_channel(time_level, direction, GAMERA_NO_LOWER);
        upper_channel =
            hdiv_target_channel(time_level, direction, GAMERA_NO_UPPER);
      } else if (family < 4) {
        lower_channel =
            hdiv_base_channel(time_level, direction, GAMERA_NO_LOWER);
        upper_channel =
            hdiv_base_channel(time_level, direction, GAMERA_NO_UPPER);
      } else {
        lower_channel = hdiv_weight_channel(direction, GAMERA_NO_LOWER);
        upper_channel = hdiv_weight_channel(direction, GAMERA_NO_UPPER);
      }
      memset(face, 0, face_count * sizeof(*face));
      for (int i = 0; i < config.ni; ++i) {
        for (int j = 0; j < config.nj; ++j) {
          for (int k = 0; k < config.nk; ++k) {
            const int global_cell[3] = {
                proc_coords[0] * config.ni + i,
                proc_coords[1] * config.nj + j,
                proc_coords[2] * config.nk + k};
            const size_t local = local_active_index(i, j, k);
            size_t q[3] = {(size_t)(local_lower[0] + i),
                           (size_t)(local_lower[1] + j),
                           (size_t)(local_lower[2] + k)};
            face[hdiv_local_face_index(direction, q[0], q[1], q[2])] =
                sentinel_local[lower_channel * local_active_cells + local];
            if (global_cell[direction] == global_shape[direction] - 1) {
              ++q[direction];
              face[hdiv_local_face_index(direction, q[0], q[1], q[2])] =
                  sentinel_local[upper_channel * local_active_cells + local];
            }
          }
        }
      }
      local_failed |= hdiv_exchange_face_field(face, direction) != 0;

      const size_t *extent = hdiv_local_face_extent[direction];
      for (size_t qi = 0; qi < extent[0]; ++qi) {
        for (size_t qj = 0; qj < extent[1]; ++qj) {
          for (size_t qk = 0; qk < extent[2]; ++qk) {
            const size_t q[3] = {qi, qj, qk};
            const int global_face[3] = {
                proc_coords[0] * config.ni + (int)qi - local_lower[0],
                proc_coords[1] * config.nj + (int)qj - local_lower[1],
                proc_coords[2] * config.nk + (int)qk - local_lower[2]};
            int inside = 1;
            for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
              inside &= global_face[axis] >= 0 &&
                        global_face[axis] <
                            global_shape[axis] +
                                (axis == direction ? 1 : 0);
            }
            if (!inside) {
              continue;
            }
            double expected = 0.0;
            if (hdiv_verify_gather_face_channel(
                    sentinel_gather, (size_t)hdiv_send_count, direction,
                    global_face, lower_channel, upper_channel, &expected) !=
                    0 ||
                !hdiv_verify_same_bits(
                    face[hdiv_local_face_index(direction, qi, qj, qk)],
                    expected)) {
              local_failed = 1;
            }

            const int normal = global_face[direction];
            if (config.proc_dims[direction] > 1 && normal > 0 &&
                normal < global_shape[direction] &&
                normal % local_shape[direction] == 0) {
              int losing_cell[3] = {global_face[0], global_face[1],
                                    global_face[2]};
              --losing_cell[direction];
              double losing = 0.0;
              if (hdiv_verify_gather_cell_channel(
                      sentinel_gather, (size_t)hdiv_send_count, losing_cell,
                      upper_channel, &losing) != 0 ||
                  hdiv_verify_same_bits(expected, losing)) {
                local_failed = 1;
              } else {
                ++local_conflicts;
              }
            }
          }
        }
      }
    }
  }

  unsigned long long patch_conflicts = 0;
  local_failed |= MPI_Allreduce(&local_conflicts, &patch_conflicts, 1,
                                MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                                comm_cart) != MPI_SUCCESS;
  int needs_conflict = 0;
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    needs_conflict |= config.proc_dims[direction] > 1;
  }
  local_failed |= needs_conflict && patch_conflicts == 0;
  const int status = hdiv_patch_consensus_failure(local_failed);
  if (status == 0 && patch_rank == 0) {
    log_info("Distributed H(div) canonical sentinel oracle passed on patch "
             "%d with %llu conflicting upper/lower interface copies",
             patch_id, patch_conflicts);
  }
  free(sentinel_local);
  free(sentinel_gather);
  free(face);
  return status;
}

static int hdiv_verify_distributed_snapshot(int include_metric_weights) {
  (void)include_metric_weights;
  const int time_level_count =
      hdiv_history_ready ? 1 : GAMERA_NO_YINYANG_TIME_LEVELS;
  if (!hdiv_distributed_verify_sentinel_ready) {
    if (hdiv_verify_sentinel_ownership() != 0) {
      return -1;
    }
    hdiv_distributed_verify_sentinel_ready = 1;
  }
  hdiv_distributed_verify_projection_mask &=
      time_level_count == 1 ? ~1U : ~3U;

  const size_t gather_count = (size_t)patch_size * (size_t)hdiv_send_count;
  double *snapshot = (double *)malloc(gather_count * sizeof(*snapshot));
  int local_failed = snapshot == NULL;
  if (hdiv_patch_consensus_failure(local_failed) != 0) {
    free(snapshot);
    return -1;
  }
  local_failed = MPI_Allgather(hdiv_send_buffer, hdiv_send_count, MPI_DOUBLE,
                               snapshot, hdiv_send_count, MPI_DOUBLE,
                               comm_cart) != MPI_SUCCESS;
  if (hdiv_patch_consensus_failure(local_failed) != 0) {
    free(snapshot);
    return -1;
  }

  const int local_lower[3] = {is, js, ks};
  const int global_shape[3] = {config.ni_global, config.nj_global,
                               config.nk_global};
  double local_max_face[2] = {0.0, 0.0};
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t *extent = hdiv_local_face_extent[direction];
    for (size_t qi = 0; qi < extent[0]; ++qi) {
      for (size_t qj = 0; qj < extent[1]; ++qj) {
        for (size_t qk = 0; qk < extent[2]; ++qk) {
          const size_t q[3] = {qi, qj, qk};
          const int global_face[3] = {
              proc_coords[0] * config.ni + (int)qi - local_lower[0],
              proc_coords[1] * config.nj + (int)qj - local_lower[1],
              proc_coords[2] * config.nk + (int)qk - local_lower[2]};
          int inside = 1;
          for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
            inside &= global_face[axis] >= 0 &&
                      global_face[axis] <
                          global_shape[axis] + (axis == direction ? 1 : 0);
          }
          if (!inside) {
            continue;
          }
          const size_t face =
              hdiv_local_face_index(direction, q[0], q[1], q[2]);
          double expected_weight = 0.0;
          if (hdiv_verify_snapshot_face_value(snapshot, 0, direction,
                                              global_face, 2,
                                              &expected_weight) != 0 ||
              !hdiv_verify_same_bits(
                  hdiv_local_face_mobility[direction][face],
                  expected_weight)) {
            local_failed = 1;
          }
          for (int time_level = 0; time_level < time_level_count;
               ++time_level) {
            double expected_flux = 0.0;
            double expected_base = 0.0;
            if (hdiv_verify_snapshot_face_value(
                    snapshot, time_level, direction, global_face, 0,
                    &expected_flux) != 0 ||
                hdiv_verify_snapshot_face_value(
                    snapshot, time_level, direction, global_face, 1,
                    &expected_base) != 0 ||
                !hdiv_verify_same_bits(
                    hdiv_local_face_flux[time_level][direction][face],
                    expected_flux) ||
                !hdiv_verify_same_bits(
                    hdiv_local_face_base[time_level][direction][face],
                    expected_base)) {
              local_failed = 1;
            }
            local_max_face[time_level] =
                fmax(local_max_face[time_level], fabs(expected_flux));
            local_max_face[time_level] =
                fmax(local_max_face[time_level], fabs(expected_base));
          }
          if (direction == GAMERA_NO_I &&
              (global_face[GAMERA_NO_I] == 0 ||
               global_face[GAMERA_NO_I] == config.ni_global)) {
            const double positive_zero = 0.0;
            local_failed |=
                !hdiv_verify_same_bits(expected_weight, positive_zero);
            for (int time_level = 0; time_level < time_level_count;
                 ++time_level) {
              local_failed |= !hdiv_verify_same_bits(
                  hdiv_local_face_flux[time_level][direction][face],
                  positive_zero);
            }
          }
        }
      }
    }
  }

  double patch_max_face[2] = {0.0, 0.0};
  for (int time_level = 0; time_level < time_level_count; ++time_level) {
    local_failed |=
        MPI_Allreduce(&local_max_face[time_level], &patch_max_face[time_level],
                      1, MPI_DOUBLE, MPI_MAX, comm_cart) != MPI_SUCCESS;
  }

  double local_max_before[2] = {0.0, 0.0};
  const size_t cell_storage_count =
      gamera_no_element_count3(hdiv_local_cell_extent);
  memset(hdiv_lambda, 0, cell_storage_count * sizeof(*hdiv_lambda));
  for (int i = 0; i < config.ni; ++i) {
    for (int j = 0; j < config.nj; ++j) {
      for (int k = 0; k < config.nk; ++k) {
        const int global_cell[3] = {proc_coords[0] * config.ni + i,
                                    proc_coords[1] * config.nj + j,
                                    proc_coords[2] * config.nk + k};
        double expected_diagonal = 0.0;
        double expected_rhs[2] = {0.0, 0.0};
        for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
          int upper_face[3] = {global_cell[0], global_cell[1],
                               global_cell[2]};
          ++upper_face[direction];
          double lower_weight = 0.0;
          double upper_weight = 0.0;
          hdiv_verify_snapshot_face_value(snapshot, 0, direction, global_cell,
                                          2, &lower_weight);
          hdiv_verify_snapshot_face_value(snapshot, 0, direction, upper_face,
                                          2, &upper_weight);
          expected_diagonal += lower_weight;
          expected_diagonal += upper_weight;
          for (int time_level = 0; time_level < time_level_count;
               ++time_level) {
            double lower_flux = 0.0;
            double upper_flux = 0.0;
            hdiv_verify_snapshot_face_value(snapshot, time_level, direction,
                                            global_cell, 0, &lower_flux);
            hdiv_verify_snapshot_face_value(snapshot, time_level, direction,
                                            upper_face, 0, &upper_flux);
            expected_rhs[time_level] += upper_flux - lower_flux;
          }
        }
        const size_t cell = hdiv_active_cell_index(i, j, k);
        local_failed |= !hdiv_verify_same_bits(
            hdiv_cell_diagonal_cache[cell], expected_diagonal);
        for (int time_level = 0; time_level < time_level_count;
             ++time_level) {
          local_failed |= !hdiv_verify_same_bits(
              hdiv_local_cell_net_flux(time_level, i, j, k),
              expected_rhs[time_level]);
          local_max_before[time_level] =
              fmax(local_max_before[time_level], fabs(expected_rhs[time_level]));
        }
        hdiv_lambda[cell] = hdiv_verify_lambda_formula(global_cell);
      }
    }
  }

  const int operator_failed =
      apply_hdiv_local_operator(hdiv_lambda, hdiv_operator_search) != 0;
  local_failed |= operator_failed;
  for (int i = 0; i < config.ni; ++i) {
    for (int j = 0; j < config.nj; ++j) {
      for (int k = 0; k < config.nk; ++k) {
        const int global_cell[3] = {proc_coords[0] * config.ni + i,
                                    proc_coords[1] * config.nj + j,
                                    proc_coords[2] * config.nk + k};
        const double center = hdiv_verify_lambda_formula(global_cell);
        double expected_operator = 0.0;
        for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
          int lower_cell[3] = {global_cell[0], global_cell[1],
                               global_cell[2]};
          int upper_cell[3] = {global_cell[0], global_cell[1],
                               global_cell[2]};
          int upper_face[3] = {global_cell[0], global_cell[1],
                               global_cell[2]};
          --lower_cell[direction];
          ++upper_cell[direction];
          ++upper_face[direction];
          const double lower_value =
              lower_cell[direction] >= 0
                  ? hdiv_verify_lambda_formula(lower_cell)
                  : 0.0;
          const double upper_value =
              upper_cell[direction] < global_shape[direction]
                  ? hdiv_verify_lambda_formula(upper_cell)
                  : 0.0;
          double lower_weight = 0.0;
          double upper_weight = 0.0;
          hdiv_verify_snapshot_face_value(snapshot, 0, direction, global_cell,
                                          2, &lower_weight);
          hdiv_verify_snapshot_face_value(snapshot, 0, direction, upper_face,
                                          2, &upper_weight);
          expected_operator += lower_weight * (center - lower_value);
          expected_operator += upper_weight * (center - upper_value);
        }
        local_failed |= !hdiv_verify_same_bits(
            hdiv_operator_search[hdiv_active_cell_index(i, j, k)],
            expected_operator);
      }
    }
  }

  /* The deterministic lambda also independently checks the face correction
   * stencil against the global Dirichlet interpretation. */
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t count[3] = {
        (size_t)config.ni + (direction == GAMERA_NO_I ? 1U : 0U),
        (size_t)config.nj + (direction == GAMERA_NO_J ? 1U : 0U),
        (size_t)config.nk + (direction == GAMERA_NO_K ? 1U : 0U)};
    for (size_t i = 0; i < count[0]; ++i) {
      for (size_t j = 0; j < count[1]; ++j) {
        for (size_t k = 0; k < count[2]; ++k) {
          const int global_face[3] = {
              proc_coords[0] * config.ni + (int)i,
              proc_coords[1] * config.nj + (int)j,
              proc_coords[2] * config.nk + (int)k};
          int lower_cell[3] = {global_face[0], global_face[1],
                               global_face[2]};
          --lower_cell[direction];
          const double lambda_lower =
              lower_cell[direction] >= 0
                  ? hdiv_verify_lambda_formula(lower_cell)
                  : 0.0;
          const double lambda_upper =
              global_face[direction] < global_shape[direction]
                  ? hdiv_verify_lambda_formula(global_face)
                  : 0.0;
          double expected_weight = 0.0;
          hdiv_verify_snapshot_face_value(snapshot, 0, direction, global_face,
                                          2, &expected_weight);
          const double expected_correction =
              -expected_weight * (lambda_lower - lambda_upper);
          size_t local_face[3] = {(size_t)local_lower[0] + i,
                                  (size_t)local_lower[1] + j,
                                  (size_t)local_lower[2] + k};
          size_t local_lower_cell[3] = {local_face[0], local_face[1],
                                        local_face[2]};
          --local_lower_cell[direction];
          const size_t face = hdiv_local_face_index(
              direction, local_face[0], local_face[1], local_face[2]);
          const double local_correction =
              -hdiv_local_face_mobility[direction][face] *
              (hdiv_lambda[hdiv_local_cell_index(
                   local_lower_cell[0], local_lower_cell[1],
                   local_lower_cell[2])] -
               hdiv_lambda[hdiv_local_cell_index(
                   local_face[0], local_face[1], local_face[2])]);
          local_failed |= !hdiv_verify_same_bits(local_correction,
                                                 expected_correction);
        }
      }
    }
  }

  const double positive_zero = 0.0;
  for (int axis = GAMERA_NO_J; axis <= GAMERA_NO_K; ++axis) {
    if (proc_coords[axis] == 0) {
      size_t begin[3] = {(size_t)is, (size_t)js, (size_t)ks};
      size_t end[3] = {begin[0] + (size_t)config.ni,
                       begin[1] + (size_t)config.nj,
                       begin[2] + (size_t)config.nk};
      begin[axis] = (size_t)local_lower[axis] - 1U;
      end[axis] = begin[axis] + 1U;
      for (size_t i = begin[0]; i < end[0]; ++i)
        for (size_t j = begin[1]; j < end[1]; ++j)
          for (size_t k = begin[2]; k < end[2]; ++k)
            local_failed |= !hdiv_verify_same_bits(
                hdiv_lambda[hdiv_local_cell_index(i, j, k)], positive_zero);
    }
    if (proc_coords[axis] == config.proc_dims[axis] - 1) {
      size_t begin[3] = {(size_t)is, (size_t)js, (size_t)ks};
      size_t end[3] = {begin[0] + (size_t)config.ni,
                       begin[1] + (size_t)config.nj,
                       begin[2] + (size_t)config.nk};
      begin[axis] = (size_t)local_lower[axis] +
                    (size_t)(axis == GAMERA_NO_I
                                 ? config.ni
                                 : (axis == GAMERA_NO_J ? config.nj
                                                        : config.nk));
      end[axis] = begin[axis] + 1U;
      for (size_t i = begin[0]; i < end[0]; ++i)
        for (size_t j = begin[1]; j < end[1]; ++j)
          for (size_t k = begin[2]; k < end[2]; ++k)
            local_failed |= !hdiv_verify_same_bits(
                hdiv_lambda[hdiv_local_cell_index(i, j, k)], positive_zero);
    }
  }

  double patch_max_before[2] = {0.0, 0.0};
  for (int time_level = 0; time_level < time_level_count; ++time_level) {
    local_failed |= MPI_Allreduce(&local_max_before[time_level],
                                  &patch_max_before[time_level], 1,
                                  MPI_DOUBLE, MPI_MAX,
                                  comm_cart) != MPI_SUCCESS;
    hdiv_distributed_verify_max_face[time_level] =
        patch_max_face[time_level];
    hdiv_distributed_verify_max_before[time_level] =
        patch_max_before[time_level];
    hdiv_distributed_verify_tolerance[time_level] =
        1.0e-14 * patch_max_before[time_level] +
        64.0 * DBL_EPSILON *
            fmax(patch_max_face[time_level], DBL_MIN);
  }

  const int status = hdiv_patch_consensus_failure(local_failed);
  if (status == 0) {
    hdiv_distributed_verify_ready = 1;
    if (patch_rank == 0) {
      log_info("Distributed H(div) %s snapshot oracle passed on patch %d: "
               "canonical faces/RHS/diag/A/correction/Dirichlet/radial-lock",
               time_level_count == 1 ? "current-only" : "current+old",
               patch_id);
    }
  }
  free(snapshot);
  return status;
}

static int hdiv_verify_projection_result(int time_level,
                                         double reduced_max_correction,
                                         double reduced_max_after) {
  const unsigned int bit = 1U << (unsigned int)time_level;
  if (!hdiv_distributed_verify_ready ||
      (hdiv_distributed_verify_projection_mask & bit) != 0U) {
    return 0;
  }
  const size_t snapshot_count =
      (size_t)patch_size * (size_t)hdiv_send_count;
  const size_t lambda_count = (size_t)patch_size * local_active_cells;
  double *snapshot = (double *)malloc(snapshot_count * sizeof(*snapshot));
  double *lambda_local =
      (double *)malloc(local_active_cells * sizeof(*lambda_local));
  double *lambda_gather =
      (double *)malloc(lambda_count * sizeof(*lambda_gather));
  int local_failed = snapshot == NULL || lambda_local == NULL ||
                     lambda_gather == NULL;
  if (hdiv_patch_consensus_failure(local_failed) != 0) {
    free(snapshot);
    free(lambda_local);
    free(lambda_gather);
    return -1;
  }
  for (int i = 0; i < config.ni; ++i)
    for (int j = 0; j < config.nj; ++j)
      for (int k = 0; k < config.nk; ++k)
        lambda_local[local_active_index(i, j, k)] =
            hdiv_lambda[hdiv_active_cell_index(i, j, k)];
  local_failed =
      MPI_Allgather(hdiv_send_buffer, hdiv_send_count, MPI_DOUBLE, snapshot,
                    hdiv_send_count, MPI_DOUBLE, comm_cart) != MPI_SUCCESS;
  local_failed |=
      MPI_Allgather(lambda_local, (int)local_active_cells, MPI_DOUBLE,
                    lambda_gather, (int)local_active_cells, MPI_DOUBLE,
                    comm_cart) != MPI_SUCCESS;
  if (hdiv_patch_consensus_failure(local_failed) != 0) {
    free(snapshot);
    free(lambda_local);
    free(lambda_gather);
    return -1;
  }

  const int global_shape[3] = {config.ni_global, config.nj_global,
                               config.nk_global};
  const int local_lower[3] = {is, js, ks};
  double expected_max_correction = 0.0;
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t *extent = hdiv_local_face_extent[direction];
    for (size_t qi = 0; qi < extent[0]; ++qi) {
      for (size_t qj = 0; qj < extent[1]; ++qj) {
        for (size_t qk = 0; qk < extent[2]; ++qk) {
          const int global_face[3] = {
              proc_coords[0] * config.ni + (int)qi - local_lower[0],
              proc_coords[1] * config.nj + (int)qj - local_lower[1],
              proc_coords[2] * config.nk + (int)qk - local_lower[2]};
          int inside = 1;
          for (int axis = 0; axis < GAMERA_NO_DIM; ++axis)
            inside &= global_face[axis] >= 0 &&
                      global_face[axis] <
                          global_shape[axis] + (axis == direction ? 1 : 0);
          if (!inside) continue;
          int lower_cell[3] = {global_face[0], global_face[1],
                               global_face[2]};
          --lower_cell[direction];
          double lambda_lower = 0.0;
          double lambda_upper = 0.0;
          double weight = 0.0;
          double initial_flux = 0.0;
          hdiv_verify_global_lambda(lambda_gather, lower_cell, &lambda_lower);
          hdiv_verify_global_lambda(lambda_gather, global_face, &lambda_upper);
          hdiv_verify_snapshot_face_value(snapshot, 0, direction, global_face,
                                          2, &weight);
          hdiv_verify_snapshot_face_value(snapshot, time_level, direction,
                                          global_face, 0, &initial_flux);
          const double correction =
              -weight * (lambda_lower - lambda_upper);
          const double expected_flux = initial_flux + correction;
          const size_t face =
              hdiv_local_face_index(direction, qi, qj, qk);
          local_failed |= !hdiv_verify_same_bits(
              hdiv_local_face_flux[time_level][direction][face],
              expected_flux);
          expected_max_correction =
              fmax(expected_max_correction, fabs(correction));
        }
      }
    }
  }
  double oracle_max_correction = 0.0;
  local_failed |= MPI_Allreduce(&expected_max_correction,
                                &oracle_max_correction, 1, MPI_DOUBLE,
                                MPI_MAX, comm_cart) != MPI_SUCCESS;
  const int correction_max_mismatch =
      !hdiv_verify_same_bits(oracle_max_correction, reduced_max_correction);
  local_failed |= correction_max_mismatch;

  double local_oracle_after = 0.0;
  for (int i = 0; i < config.ni; ++i) {
    for (int j = 0; j < config.nj; ++j) {
      for (int k = 0; k < config.nk; ++k) {
        if (hdiv_cell_diagonal_cache[hdiv_active_cell_index(i, j, k)] > 0.0)
          local_oracle_after =
              fmax(local_oracle_after,
                   fabs(hdiv_local_cell_net_flux(time_level, i, j, k)));
      }
    }
  }
  double oracle_max_after = 0.0;
  local_failed |= MPI_Allreduce(&local_oracle_after, &oracle_max_after, 1,
                                MPI_DOUBLE, MPI_MAX,
                                comm_cart) != MPI_SUCCESS;
  const int after_max_mismatch =
      !hdiv_verify_same_bits(oracle_max_after, reduced_max_after);
  local_failed |= after_max_mismatch;
  if (patch_rank == 0 && (correction_max_mismatch || after_max_mismatch)) {
    log_error("Distributed H(div) projection reduction oracle mismatch "
              "patch=%d level=%d correction actual=%.17g expected=%.17g "
              "after actual=%.17g expected=%.17g",
              patch_id, time_level, reduced_max_correction,
              oracle_max_correction, reduced_max_after, oracle_max_after);
  }
  const int status = hdiv_patch_consensus_failure(local_failed);
  if (status == 0) {
    hdiv_distributed_verify_projection_mask |= bit;
    if (patch_rank == 0) {
      log_info("Distributed H(div) projected-face oracle passed on patch %d "
               "level=%d correction=%.6e after=%.6e",
               patch_id, time_level, reduced_max_correction,
               reduced_max_after);
    }
  }
  free(snapshot);
  free(lambda_local);
  free(lambda_gather);
  return status;
}

static int hdiv_verify_capture_old_storage(const gamera_no_storage *storage,
                                           const gamera_no_grid *grid) {
  int local_failed = storage == NULL || grid == NULL;
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t face_count =
        gamera_no_element_count3(grid->face[direction].extent);
    free(hdiv_distributed_verify_old_face[direction]);
    hdiv_distributed_verify_old_face[direction] =
        (double *)malloc(face_count *
                         sizeof(*hdiv_distributed_verify_old_face[direction]));
    local_failed |= hdiv_distributed_verify_old_face[direction] == NULL;
  }
  if (hdiv_patch_consensus_failure(local_failed) != 0) {
    for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
      free(hdiv_distributed_verify_old_face[direction]);
      hdiv_distributed_verify_old_face[direction] = NULL;
    }
    return -1;
  }
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t face_count =
        gamera_no_element_count3(grid->face[direction].extent);
    memcpy(hdiv_distributed_verify_old_face[direction],
           storage->old_face_flux[direction],
           face_count * sizeof(*storage->old_face_flux[direction]));
  }
  return 0;
}

static int hdiv_verify_old_storage_unchanged(
    const gamera_no_storage *storage, const gamera_no_grid *grid) {
  int local_failed = storage == NULL || grid == NULL;
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t face_count =
        gamera_no_element_count3(grid->face[direction].extent);
    local_failed |= hdiv_distributed_verify_old_face[direction] == NULL;
    if (hdiv_distributed_verify_old_face[direction] != NULL &&
        memcmp(hdiv_distributed_verify_old_face[direction],
               storage->old_face_flux[direction],
               face_count * sizeof(*storage->old_face_flux[direction])) != 0) {
      local_failed = 1;
    }
  }
  const int status = hdiv_patch_consensus_failure(local_failed);
  if (status == 0 && patch_rank == 0) {
    log_info("Distributed H(div) current-only audit kept old CT faces bitwise "
             "unchanged on patch %d",
             patch_id);
  }
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    free(hdiv_distributed_verify_old_face[direction]);
    hdiv_distributed_verify_old_face[direction] = NULL;
  }
  if (status == 0) {
    hdiv_distributed_verify_old_audit_done = 1;
  }
  return status;
}
#endif

static int project_hdiv_level(int time_level) {
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  double vector_start = 0.0;
#endif
  double local_max_face_flux = 0.0;
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t lower[3] = {(size_t)is, (size_t)js, (size_t)ks};
    const size_t count[3] = {
        (size_t)config.ni + (direction == GAMERA_NO_I ? 1U : 0U),
        (size_t)config.nj + (direction == GAMERA_NO_J ? 1U : 0U),
        (size_t)config.nk + (direction == GAMERA_NO_K ? 1U : 0U)};
#pragma omp parallel for collapse(3) schedule(static) reduction(max : local_max_face_flux)
    for (size_t i = 0; i < count[0]; ++i) {
      for (size_t j = 0; j < count[1]; ++j) {
        for (size_t k = 0; k < count[2]; ++k) {
          const size_t face = hdiv_local_face_index(
              direction, lower[0] + i, lower[1] + j, lower[2] + k);
          local_max_face_flux =
              fmax(local_max_face_flux,
                   fabs(hdiv_local_face_flux[time_level][direction][face]));
          local_max_face_flux =
              fmax(local_max_face_flux,
                   fabs(hdiv_local_face_base[time_level][direction][face]));
        }
      }
    }
  }
  double max_face_flux = 0.0;
  if (hdiv_patch_max(local_max_face_flux, &max_face_flux) != 0) {
    return -1;
  }
#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED_VERIFY
  if (hdiv_distributed_verify_ready &&
      (hdiv_distributed_verify_projection_mask &
       (1U << (unsigned int)time_level)) == 0U &&
      !hdiv_verify_same_bits(
          max_face_flux, hdiv_distributed_verify_max_face[time_level])) {
    if (patch_rank == 0) {
      log_error("Distributed H(div) max-face oracle mismatch patch=%d "
                "level=%d actual=%.17g expected=%.17g",
                patch_id, time_level, max_face_flux,
                hdiv_distributed_verify_max_face[time_level]);
    }
    return -1;
  }
#endif

  const size_t cell_storage_count =
      gamera_no_element_count3(hdiv_local_cell_extent);
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  vector_start = omp_get_wtime();
#endif
  memset(hdiv_lambda, 0, cell_storage_count * sizeof(*hdiv_lambda));
  double local_max_before = 0.0;
#pragma omp parallel for collapse(3) schedule(static) reduction(max : local_max_before)
  for (int i = 0; i < config.ni; ++i) {
    for (int j = 0; j < config.nj; ++j) {
      for (int k = 0; k < config.nk; ++k) {
        const size_t cell = hdiv_active_cell_index(i, j, k);
        const double diagonal = hdiv_cell_diagonal_cache[cell];
        if (diagonal > 0.0) {
          hdiv_residual[cell] = hdiv_local_cell_net_flux(time_level, i, j, k);
          local_max_before =
              fmax(local_max_before, fabs(hdiv_residual[cell]));
          hdiv_preconditioned[cell] = hdiv_residual[cell] / diagonal;
        } else {
          hdiv_residual[cell] = 0.0;
          hdiv_preconditioned[cell] = 0.0;
        }
        hdiv_search[cell] = hdiv_preconditioned[cell];
      }
    }
  }
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  hdiv_profile_detail_current[GAMERA_NO_HDIV_DETAIL_VECTOR] +=
      omp_get_wtime() - vector_start;
#endif
  double max_before = 0.0;
  if (hdiv_patch_max(local_max_before, &max_before) != 0) {
    return -1;
  }
#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED_VERIFY
  if (hdiv_distributed_verify_ready &&
      (hdiv_distributed_verify_projection_mask &
       (1U << (unsigned int)time_level)) == 0U &&
      !hdiv_verify_same_bits(
          max_before, hdiv_distributed_verify_max_before[time_level])) {
    if (patch_rank == 0) {
      log_error("Distributed H(div) max-before oracle mismatch patch=%d "
                "level=%d actual=%.17g expected=%.17g",
                patch_id, time_level, max_before,
                hdiv_distributed_verify_max_before[time_level]);
    }
    return -1;
  }
#endif
  gamera_no_yinyang_hdiv_max_before =
      fmax(gamera_no_yinyang_hdiv_max_before, max_before);
  if (max_before <=
      64.0 * DBL_EPSILON * fmax(max_face_flux, DBL_MIN)) {
    return 0;
  }

  double rho = 0.0;
  if (hdiv_patch_sum(hdiv_local_dot(hdiv_residual, hdiv_preconditioned),
                     &rho) != 0 ||
      !isfinite(rho) || rho < 0.0) {
    return -1;
  }
  const double tolerance =
      1.0e-14 * max_before +
      64.0 * DBL_EPSILON * fmax(max_face_flux, DBL_MIN);
#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED_VERIFY
  if (hdiv_distributed_verify_ready &&
      (hdiv_distributed_verify_projection_mask &
       (1U << (unsigned int)time_level)) == 0U &&
      !hdiv_verify_same_bits(
          tolerance, hdiv_distributed_verify_tolerance[time_level])) {
    if (patch_rank == 0) {
      log_error("Distributed H(div) tolerance oracle mismatch patch=%d "
                "level=%d actual=%.17g expected=%.17g",
                patch_id, time_level, tolerance,
                hdiv_distributed_verify_tolerance[time_level]);
    }
    return -1;
  }
#endif
  const int maximum_iterations = 4000;
  int iterations = 0;
  int converged = 0;
  while (iterations < maximum_iterations) {
    const int operator_failed =
        apply_hdiv_local_operator(hdiv_search, hdiv_operator_search) != 0;
    double denominator = 0.0;
    if (hdiv_patch_sum_with_failure(
            hdiv_local_dot(hdiv_search, hdiv_operator_search),
            operator_failed, &denominator) != 0 ||
        !isfinite(denominator) || denominator <= 0.0) {
      return -1;
    }
    const double alpha = rho / denominator;
    double local_residual_maximum = 0.0;
#ifdef GAMERA_YINYANG_HDIV_PROFILE
    vector_start = omp_get_wtime();
#endif
#pragma omp parallel for collapse(3) schedule(static) reduction(max : local_residual_maximum)
    for (int i = 0; i < config.ni; ++i) {
      for (int j = 0; j < config.nj; ++j) {
        for (int k = 0; k < config.nk; ++k) {
          const size_t cell = hdiv_active_cell_index(i, j, k);
          hdiv_lambda[cell] += alpha * hdiv_search[cell];
          hdiv_residual[cell] -= alpha * hdiv_operator_search[cell];
          local_residual_maximum =
              fmax(local_residual_maximum, fabs(hdiv_residual[cell]));
        }
      }
    }
#ifdef GAMERA_YINYANG_HDIV_PROFILE
    hdiv_profile_detail_current[GAMERA_NO_HDIV_DETAIL_VECTOR] +=
        omp_get_wtime() - vector_start;
#endif
    ++iterations;
    double residual_maximum = 0.0;
    if (hdiv_patch_max(local_residual_maximum, &residual_maximum) != 0) {
      return -1;
    }
    if (residual_maximum <= tolerance) {
      converged = 1;
      break;
    }

#ifdef GAMERA_YINYANG_HDIV_PROFILE
    vector_start = omp_get_wtime();
#endif
#pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < config.ni; ++i) {
      for (int j = 0; j < config.nj; ++j) {
        for (int k = 0; k < config.nk; ++k) {
          const size_t cell = hdiv_active_cell_index(i, j, k);
          const double diagonal = hdiv_cell_diagonal_cache[cell];
          hdiv_preconditioned[cell] =
              diagonal > 0.0 ? hdiv_residual[cell] / diagonal : 0.0;
        }
      }
    }
#ifdef GAMERA_YINYANG_HDIV_PROFILE
    hdiv_profile_detail_current[GAMERA_NO_HDIV_DETAIL_VECTOR] +=
        omp_get_wtime() - vector_start;
#endif
    double next_rho = 0.0;
    if (hdiv_patch_sum(hdiv_local_dot(hdiv_residual, hdiv_preconditioned),
                       &next_rho) != 0 ||
        !isfinite(next_rho) || next_rho < 0.0 || rho == 0.0) {
      return -1;
    }
    const double beta = next_rho / rho;
#ifdef GAMERA_YINYANG_HDIV_PROFILE
    vector_start = omp_get_wtime();
#endif
#pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < config.ni; ++i) {
      for (int j = 0; j < config.nj; ++j) {
        for (int k = 0; k < config.nk; ++k) {
          const size_t cell = hdiv_active_cell_index(i, j, k);
          hdiv_search[cell] =
              hdiv_preconditioned[cell] + beta * hdiv_search[cell];
        }
      }
    }
#ifdef GAMERA_YINYANG_HDIV_PROFILE
    hdiv_profile_detail_current[GAMERA_NO_HDIV_DETAIL_VECTOR] +=
        omp_get_wtime() - vector_start;
#endif
    rho = next_rho;
  }
  if (!converged) {
    return -1;
  }
  gamera_no_yinyang_hdiv_max_iterations =
      gamera_no_yinyang_hdiv_max_iterations > iterations
          ? gamera_no_yinyang_hdiv_max_iterations
          : iterations;
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  hdiv_profile_pcg_iterations_current[time_level] =
      (unsigned long)iterations;
#endif

  int correction_failed = hdiv_exchange_cell_halo(hdiv_lambda) != 0;
  double local_max_correction = 0.0;
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t lower[3] = {(size_t)is, (size_t)js, (size_t)ks};
    const size_t count[3] = {
        (size_t)config.ni + (direction == GAMERA_NO_I ? 1U : 0U),
        (size_t)config.nj + (direction == GAMERA_NO_J ? 1U : 0U),
        (size_t)config.nk + (direction == GAMERA_NO_K ? 1U : 0U)};
#ifdef GAMERA_YINYANG_HDIV_PROFILE
    const double correction_start = omp_get_wtime();
#endif
#pragma omp parallel for collapse(3) schedule(static) reduction(max : local_max_correction)
    for (size_t i = 0; i < count[0]; ++i) {
      for (size_t j = 0; j < count[1]; ++j) {
        for (size_t k = 0; k < count[2]; ++k) {
          const size_t face_coordinate[3] = {
              lower[0] + i, lower[1] + j, lower[2] + k};
          size_t lower_cell[3] = {face_coordinate[0], face_coordinate[1],
                                  face_coordinate[2]};
          --lower_cell[direction];
          const size_t face = hdiv_local_face_index(
              direction, face_coordinate[0], face_coordinate[1],
              face_coordinate[2]);
          const double correction =
              -hdiv_local_face_mobility[direction][face] *
              (hdiv_lambda[hdiv_local_cell_index(
                   lower_cell[0], lower_cell[1], lower_cell[2])] -
               hdiv_lambda[hdiv_local_cell_index(
                   face_coordinate[0], face_coordinate[1],
                   face_coordinate[2])]);
          hdiv_local_face_flux[time_level][direction][face] += correction;
          local_max_correction =
              fmax(local_max_correction, fabs(correction));
        }
      }
    }
#ifdef GAMERA_YINYANG_HDIV_PROFILE
    hdiv_profile_detail_current[GAMERA_NO_HDIV_DETAIL_CORRECTION] +=
        omp_get_wtime() - correction_start;
#endif
    correction_failed |=
        hdiv_exchange_face_field(
            hdiv_local_face_flux[time_level][direction], direction) != 0;
  }
  double max_correction = 0.0;
  if (hdiv_patch_max_with_failure(local_max_correction, correction_failed,
                                  &max_correction) != 0) {
    return -1;
  }
  gamera_no_yinyang_hdiv_max_correction =
      fmax(gamera_no_yinyang_hdiv_max_correction, max_correction);

  double local_max_after = 0.0;
#pragma omp parallel for collapse(3) schedule(static) reduction(max : local_max_after)
  for (int i = 0; i < config.ni; ++i) {
    for (int j = 0; j < config.nj; ++j) {
      for (int k = 0; k < config.nk; ++k) {
        const size_t cell = hdiv_active_cell_index(i, j, k);
        if (hdiv_cell_diagonal_cache[cell] > 0.0) {
          local_max_after =
              fmax(local_max_after,
                   fabs(hdiv_local_cell_net_flux(time_level, i, j, k)));
        }
      }
    }
  }
  double max_after = 0.0;
  if (hdiv_patch_max(local_max_after, &max_after) != 0) {
    return -1;
  }
  gamera_no_yinyang_hdiv_max_after =
      fmax(gamera_no_yinyang_hdiv_max_after, max_after);
#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED_VERIFY
  if (hdiv_verify_projection_result(time_level, max_correction, max_after) !=
      0) {
    return -1;
  }
#endif
  return max_after <=
                 fmax(16.0 * tolerance,
                      1024.0 * DBL_EPSILON *
                          fmax(max_face_flux, DBL_MIN))
             ? 0
             : -1;
}

static void scatter_hdiv_local_faces(gamera_no_storage *storage,
                                     const gamera_no_grid *grid) {
  const int local_lower[3] = {is, js, ks};
  const int global_cell_count[3] = {config.ni_global, config.nj_global,
                                    config.nk_global};
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t *extent = grid->face[direction].extent;
#pragma omp parallel for collapse(3) schedule(static)
    for (size_t qi = 0; qi < extent[0]; ++qi) {
      for (size_t qj = 0; qj < extent[1]; ++qj) {
        for (size_t qk = 0; qk < extent[2]; ++qk) {
          const size_t coordinate[3] = {qi, qj, qk};
          const int global[3] = {
              proc_coords[0] * config.ni + (int)qi - local_lower[0],
              proc_coords[1] * config.nj + (int)qj - local_lower[1],
              proc_coords[2] * config.nk + (int)qk - local_lower[2]};
          int inside_patch = 1;
          for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
            const int upper =
                global_cell_count[axis] + (axis == direction ? 1 : 0);
            inside_patch &= global[axis] >= 0 && global[axis] < upper;
          }
          if (!inside_patch) {
            continue;
          }
          const size_t face =
              hdiv_local_face_index(direction, qi, qj, qk);
          storage->face_flux[direction][face] =
              hdiv_local_face_base[0][direction][face] +
              hdiv_local_face_flux[0][direction][face];
          if (!hdiv_history_ready) {
            storage->old_face_flux[direction][face] =
                hdiv_local_face_base[1][direction][face] +
                hdiv_local_face_flux[1][direction][face];
          }
        }
      }
    }
  }
}
#else
static int assemble_hdiv_global_faces(int include_metric_weights,
                                      size_t gather_count) {
  const int global_count[3] = {config.ni_global, config.nj_global,
                               config.nk_global};
  for (int patch_rank = 0; patch_rank < patch_size; ++patch_rank) {
    int coords[3];
    if (MPI_Cart_coords(comm_cart, patch_rank, GAMERA_NO_DIM, coords) !=
        MPI_SUCCESS) {
      return -1;
    }
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
    const size_t rank_offset =
        (size_t)patch_rank * gather_count;
#else
    const int world_rank = patch_id * patch_size + patch_rank;
    const size_t rank_offset =
        (size_t)world_rank * (size_t)hdiv_send_count;
#endif
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
#pragma omp parallel for collapse(3) schedule(static)
#endif
    for (int i = 0; i < config.ni; ++i) {
      for (int j = 0; j < config.nj; ++j) {
        for (int k = 0; k < config.nk; ++k) {
          const int global[3] = {coords[0] * config.ni + i,
                                 coords[1] * config.nj + j,
                                 coords[2] * config.nk + k};
          const size_t local = local_active_index(i, j, k);
          for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
            int face_coordinate[3] = {global[0], global[1], global[2]};
            size_t face = hdiv_global_face_index(
                direction, face_coordinate[0], face_coordinate[1],
                face_coordinate[2]);
            for (int time_level = 0;
                 time_level < GAMERA_NO_YINYANG_TIME_LEVELS; ++time_level) {
              const size_t channel =
                  hdiv_target_channel(time_level, direction, GAMERA_NO_LOWER);
              const size_t base_channel =
                  hdiv_base_channel(time_level, direction, GAMERA_NO_LOWER);
              hdiv_global_face_base[time_level][direction][face] =
                  hdiv_gather_buffer[rank_offset +
                                      base_channel * local_active_cells +
                                      local];
              hdiv_global_face_flux[time_level][direction][face] =
                  hdiv_gather_buffer[rank_offset +
                                      channel * local_active_cells + local];
            }
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
            if (include_metric_weights) {
#endif
              const size_t weight_channel =
                  hdiv_weight_channel(direction, GAMERA_NO_LOWER);
              hdiv_global_face_mobility[direction][face] =
                  hdiv_gather_buffer[rank_offset +
                                      weight_channel * local_active_cells +
                                      local];
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
            }
#endif
            if (global[direction] == global_count[direction] - 1) {
              ++face_coordinate[direction];
              face = hdiv_global_face_index(
                  direction, face_coordinate[0], face_coordinate[1],
                  face_coordinate[2]);
              for (int time_level = 0;
                   time_level < GAMERA_NO_YINYANG_TIME_LEVELS; ++time_level) {
                const size_t channel = hdiv_target_channel(
                    time_level, direction, GAMERA_NO_UPPER);
                const size_t base_channel = hdiv_base_channel(
                    time_level, direction, GAMERA_NO_UPPER);
                hdiv_global_face_base[time_level][direction][face] =
                    hdiv_gather_buffer[rank_offset +
                                        base_channel * local_active_cells +
                                        local];
                hdiv_global_face_flux[time_level][direction][face] =
                    hdiv_gather_buffer[rank_offset +
                                        channel * local_active_cells + local];
              }
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
              if (include_metric_weights) {
#endif
                const size_t weight_channel =
                    hdiv_weight_channel(direction, GAMERA_NO_UPPER);
                hdiv_global_face_mobility[direction][face] =
                    hdiv_gather_buffer[rank_offset +
                                        weight_channel * local_active_cells +
                                        local];
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
              }
#endif
            }
          }
        }
      }
    }
  }
  /* The same canonical owner used for face flux also supplies its Hodge. */
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    for (int i = 0; i < (int)hdiv_global_face_extent[direction][0]; ++i) {
      for (int j = 0; j < (int)hdiv_global_face_extent[direction][1]; ++j) {
        for (int k = 0; k < (int)hdiv_global_face_extent[direction][2]; ++k) {
          const size_t face = hdiv_global_face_index(direction, i, j, k);
          const double weight = hdiv_global_face_mobility[direction][face];
          const int locked_radial_face =
              direction == GAMERA_NO_I &&
              (i == 0 || i == config.ni_global);
          if (include_metric_weights &&
              (!isfinite(weight) || weight < 0.0 ||
               (locked_radial_face ? weight != 0.0 : weight == 0.0))) {
            log_error("Invalid assembled H(div) weight on patch %d face "
                      "d=%d global=(%d,%d,%d): %.17g",
                      patch_id, direction, i, j, k, weight);
            return -1;
          }
          if (locked_radial_face &&
              (hdiv_global_face_flux[0][direction][face] != 0.0 ||
               hdiv_global_face_flux[1][direction][face] != 0.0)) {
            log_error("H(div) donor target attempted to change locked radial "
                      "face on patch %d global=(%d,%d,%d)",
                      patch_id, i, j, k);
            return -1;
          }
          if (include_metric_weights && weight > 0.0) {
            gamera_no_yinyang_hdiv_min_weight =
                fmin(gamera_no_yinyang_hdiv_min_weight, weight);
            gamera_no_yinyang_hdiv_max_weight =
                fmax(gamera_no_yinyang_hdiv_max_weight, weight);
          }
        }
      }
    }
  }
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
  if (!hdiv_cell_diagonal_cache_ready) {
#pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < config.ni_global; ++i) {
      for (int j = 0; j < config.nj_global; ++j) {
        for (int k = 0; k < config.nk_global; ++k) {
          const size_t cell = hdiv_global_cell_index(i, j, k);
          hdiv_cell_diagonal_cache[cell] = hdiv_cell_diagonal(i, j, k);
        }
      }
    }
    hdiv_cell_diagonal_cache_ready = 1;
  }
#endif
  return 0;
}

static double hdiv_cell_net_flux(int time_level, int i, int j, int k) {
  const int coordinate[3] = {i, j, k};
  double net = 0.0;
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    int upper[3] = {coordinate[0], coordinate[1], coordinate[2]};
    ++upper[direction];
    const size_t lower_face = hdiv_global_face_index(
        direction, coordinate[0], coordinate[1], coordinate[2]);
    const size_t upper_face = hdiv_global_face_index(
        direction, upper[0], upper[1], upper[2]);
    net += hdiv_global_face_flux[time_level][direction][upper_face] -
           hdiv_global_face_flux[time_level][direction][lower_face];
  }
  return net;
}

static double hdiv_cell_diagonal(int i, int j, int k) {
  const int coordinate[3] = {i, j, k};
  double diagonal = 0.0;
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    int upper[3] = {coordinate[0], coordinate[1], coordinate[2]};
    ++upper[direction];
    diagonal += hdiv_global_face_mobility[direction][hdiv_global_face_index(
        direction, coordinate[0], coordinate[1], coordinate[2])];
    diagonal += hdiv_global_face_mobility[direction][hdiv_global_face_index(
        direction, upper[0], upper[1], upper[2])];
  }
  return diagonal;
}

static void apply_hdiv_operator(const double *input, double *output) {
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
#pragma omp parallel for collapse(3) schedule(static)
#endif
  for (int i = 0; i < config.ni_global; ++i) {
    for (int j = 0; j < config.nj_global; ++j) {
      for (int k = 0; k < config.nk_global; ++k) {
        const int coordinate[3] = {i, j, k};
        const size_t cell = hdiv_global_cell_index(i, j, k);
        double value = 0.0;
        for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
          int lower[3] = {coordinate[0], coordinate[1], coordinate[2]};
          int upper[3] = {coordinate[0], coordinate[1], coordinate[2]};
          --lower[direction];
          ++upper[direction];
          const size_t lower_face = hdiv_global_face_index(
              direction, coordinate[0], coordinate[1], coordinate[2]);
          const size_t upper_face = hdiv_global_face_index(
              direction, upper[0], upper[1], upper[2]);
          const double lower_value =
              coordinate[direction] > 0
                  ? input[hdiv_global_cell_index(lower[0], lower[1], lower[2])]
                  : 0.0;
          const double upper_value =
              coordinate[direction] + 1 <
                      (direction == GAMERA_NO_I
                           ? config.ni_global
                           : (direction == GAMERA_NO_J ? config.nj_global
                                                      : config.nk_global))
                  ? input[hdiv_global_cell_index(upper[0], upper[1], upper[2])]
                  : 0.0;
          value += hdiv_global_face_mobility[direction][lower_face] *
                   (input[cell] - lower_value);
          value += hdiv_global_face_mobility[direction][upper_face] *
                   (input[cell] - upper_value);
        }
        output[cell] = value;
      }
    }
  }
}

static double hdiv_dot(const double *left, const double *right,
                       size_t count) {
  double sum = 0.0;
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
  /*
   * A native OpenMP reduction can merge partial sums in a team-dependent
   * order and consequently move the PCG trajectory by a few ulps.  Preserve
   * the legacy serial expression when only one thread is requested.  With a
   * larger team, evaluate runtime-mesh-derived fixed blocks in parallel and
   * merge them in a fixed order, so repeated runs and different team sizes
   * follow the same candidate trajectory.
   */
  if (omp_get_max_threads() > 1 && hdiv_dot_block_partial != NULL &&
      hdiv_dot_block_count > 1) {
    const size_t block_count =
        hdiv_dot_block_count < count ? hdiv_dot_block_count : count;
#pragma omp parallel for schedule(static)
    for (size_t block = 0; block < block_count; ++block) {
      const size_t begin = block * count / block_count;
      const size_t end = (block + 1U) * count / block_count;
      double partial = 0.0;
      for (size_t item = begin; item < end; ++item) {
        partial += left[item] * right[item];
      }
      hdiv_dot_block_partial[block] = partial;
    }
    for (size_t block = 0; block < block_count; ++block) {
      sum += hdiv_dot_block_partial[block];
    }
    return sum;
  }
#endif
  for (size_t item = 0; item < count; ++item) {
    sum += left[item] * right[item];
  }
  return sum;
}

static int project_hdiv_level(int time_level) {
  const size_t cell_count =
      (size_t)config.ni_global * (size_t)config.nj_global *
      (size_t)config.nk_global;
  double max_before = 0.0;
  double max_face_flux = 0.0;
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t face_count = gamera_no_element_count3(
        hdiv_global_face_extent[direction]);
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
#pragma omp parallel for schedule(static) reduction(max : max_face_flux)
#endif
    for (size_t face = 0; face < face_count; ++face) {
      max_face_flux = fmax(
          max_face_flux,
          fabs(hdiv_global_face_flux[time_level][direction][face]));
      max_face_flux = fmax(
          max_face_flux,
          fabs(hdiv_global_face_base[time_level][direction][face]));
    }
  }
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
#pragma omp parallel for collapse(3) schedule(static) reduction(max : max_before)
#endif
  for (int i = 0; i < config.ni_global; ++i) {
    for (int j = 0; j < config.nj_global; ++j) {
      for (int k = 0; k < config.nk_global; ++k) {
        const size_t cell = hdiv_global_cell_index(i, j, k);
        hdiv_lambda[cell] = 0.0;
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
        const double diagonal = hdiv_cell_diagonal_cache[cell];
        if (diagonal > 0.0) {
          hdiv_residual[cell] = hdiv_cell_net_flux(time_level, i, j, k);
          max_before = fmax(max_before, fabs(hdiv_residual[cell]));
          hdiv_preconditioned[cell] = hdiv_residual[cell] / diagonal;
#else
        const double diagonal = hdiv_cell_diagonal(i, j, k);
        if (diagonal > 0.0) {
          hdiv_residual[cell] = hdiv_cell_net_flux(time_level, i, j, k);
          max_before = fmax(max_before, fabs(hdiv_residual[cell]));
          hdiv_preconditioned[cell] = hdiv_residual[cell] / diagonal;
#endif
        } else {
          hdiv_preconditioned[cell] = 0.0;
          hdiv_residual[cell] = 0.0;
        }
        hdiv_search[cell] = hdiv_preconditioned[cell];
      }
    }
  }
  /*
   * Project only the donor-induced increment.  The residual field can carry
   * a nonzero, time-invariant discrete divergence which exactly cancels a
   * cut/background field.  Projecting base+increment to zero would destroy
   * that cancellation even though the total field was initially solenoidal.
   */
  gamera_no_yinyang_hdiv_max_before =
      fmax(gamera_no_yinyang_hdiv_max_before, max_before);
  if (max_before <=
      64.0 * DBL_EPSILON * fmax(max_face_flux, DBL_MIN)) {
    return 0;
  }

  double rho = hdiv_dot(hdiv_residual, hdiv_preconditioned, cell_count);
  if (!isfinite(rho) || rho < 0.0) {
    return -1;
  }
  const double tolerance =
      1.0e-14 * max_before +
      64.0 * DBL_EPSILON * fmax(max_face_flux, DBL_MIN);
  int iterations = 0;
  int converged = 0;
  const int maximum_iterations = 4000;
  while (iterations < maximum_iterations) {
    apply_hdiv_operator(hdiv_search, hdiv_operator_search);
    const double denominator =
        hdiv_dot(hdiv_search, hdiv_operator_search, cell_count);
    if (!isfinite(denominator) || denominator <= 0.0) {
      return -1;
    }
    const double alpha = rho / denominator;
    double residual_maximum = 0.0;
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
#pragma omp parallel for schedule(static) reduction(max : residual_maximum)
#endif
    for (size_t cell = 0; cell < cell_count; ++cell) {
      hdiv_lambda[cell] += alpha * hdiv_search[cell];
      hdiv_residual[cell] -= alpha * hdiv_operator_search[cell];
      residual_maximum =
          fmax(residual_maximum, fabs(hdiv_residual[cell]));
    }
    ++iterations;
    if (residual_maximum <= tolerance) {
      converged = 1;
      break;
    }
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
#pragma omp parallel for schedule(static)
    for (size_t cell = 0; cell < cell_count; ++cell) {
      const double diagonal = hdiv_cell_diagonal_cache[cell];
      hdiv_preconditioned[cell] =
          diagonal > 0.0 ? hdiv_residual[cell] / diagonal : 0.0;
    }
#else
    for (int i = 0; i < config.ni_global; ++i) {
      for (int j = 0; j < config.nj_global; ++j) {
        for (int k = 0; k < config.nk_global; ++k) {
          const size_t cell = hdiv_global_cell_index(i, j, k);
          const double diagonal = hdiv_cell_diagonal(i, j, k);
          hdiv_preconditioned[cell] =
              diagonal > 0.0 ? hdiv_residual[cell] / diagonal : 0.0;
        }
      }
    }
#endif
    const double next_rho =
        hdiv_dot(hdiv_residual, hdiv_preconditioned, cell_count);
    if (!isfinite(next_rho) || next_rho < 0.0 || rho == 0.0) {
      return -1;
    }
    const double beta = next_rho / rho;
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
#pragma omp parallel for schedule(static)
#endif
    for (size_t cell = 0; cell < cell_count; ++cell) {
      hdiv_search[cell] =
          hdiv_preconditioned[cell] + beta * hdiv_search[cell];
    }
    rho = next_rho;
  }
  if (!converged) {
    return -1;
  }
  gamera_no_yinyang_hdiv_max_iterations =
      gamera_no_yinyang_hdiv_max_iterations > iterations
          ? gamera_no_yinyang_hdiv_max_iterations
          : iterations;

  double max_correction = 0.0;
  const int global_count[3] = {config.ni_global, config.nj_global,
                               config.nk_global};
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
#pragma omp parallel for collapse(3) schedule(static) reduction(max : max_correction)
#endif
    for (int i = 0; i < (int)hdiv_global_face_extent[direction][0]; ++i) {
      for (int j = 0; j < (int)hdiv_global_face_extent[direction][1]; ++j) {
        for (int k = 0; k < (int)hdiv_global_face_extent[direction][2]; ++k) {
          const int coordinate[3] = {i, j, k};
          const int normal = coordinate[direction];
          double lambda_lower = 0.0;
          double lambda_upper = 0.0;
          if (normal > 0) {
            int lower_cell[3] = {i, j, k};
            --lower_cell[direction];
            lambda_lower = hdiv_lambda[hdiv_global_cell_index(
                lower_cell[0], lower_cell[1], lower_cell[2])];
          }
          if (normal < global_count[direction]) {
            int upper_cell[3] = {i, j, k};
            lambda_upper = hdiv_lambda[hdiv_global_cell_index(
                upper_cell[0], upper_cell[1], upper_cell[2])];
          }
          const size_t face =
              hdiv_global_face_index(direction, i, j, k);
          const double correction =
              -hdiv_global_face_mobility[direction][face] *
              (lambda_lower - lambda_upper);
          hdiv_global_face_flux[time_level][direction][face] += correction;
          max_correction = fmax(max_correction, fabs(correction));
        }
      }
    }
  }
  gamera_no_yinyang_hdiv_max_correction =
      fmax(gamera_no_yinyang_hdiv_max_correction, max_correction);

  double max_after = 0.0;
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
#pragma omp parallel for collapse(3) schedule(static) reduction(max : max_after)
#endif
  for (int i = 0; i < config.ni_global; ++i) {
    for (int j = 0; j < config.nj_global; ++j) {
      for (int k = 0; k < config.nk_global; ++k) {
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
        const size_t cell = hdiv_global_cell_index(i, j, k);
        if (hdiv_cell_diagonal_cache[cell] > 0.0) {
#else
        if (hdiv_cell_diagonal(i, j, k) > 0.0) {
#endif
          max_after = fmax(max_after,
                           fabs(hdiv_cell_net_flux(time_level, i, j, k)));
        }
      }
    }
  }
  gamera_no_yinyang_hdiv_max_after =
      fmax(gamera_no_yinyang_hdiv_max_after, max_after);
  return max_after <=
                 fmax(16.0 * tolerance,
                      1024.0 * DBL_EPSILON *
                          fmax(max_face_flux, DBL_MIN))
             ? 0
             : -1;
}

static void scatter_hdiv_global_faces(gamera_no_storage *storage,
                                      const gamera_no_grid *grid) {
  const int local_lower[3] = {is, js, ks};
  const int global_cell_count[3] = {config.ni_global, config.nj_global,
                                    config.nk_global};
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t *extent = grid->face[direction].extent;
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
#pragma omp parallel for collapse(3) schedule(static)
#endif
    for (size_t qi = 0; qi < extent[0]; ++qi) {
      for (size_t qj = 0; qj < extent[1]; ++qj) {
        for (size_t qk = 0; qk < extent[2]; ++qk) {
          const size_t coordinate[3] = {qi, qj, qk};
          const int global[3] = {
              proc_coords[0] * config.ni + (int)qi - local_lower[0],
              proc_coords[1] * config.nj + (int)qj - local_lower[1],
              proc_coords[2] * config.nk + (int)qk - local_lower[2]};
          int inside_patch = 1;
          for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
            const int upper =
                global_cell_count[axis] + (axis == direction ? 1 : 0);
            inside_patch &= global[axis] >= 0 && global[axis] < upper;
          }
          if (!inside_patch) {
            continue;
          }
          const size_t local_face = gamera_no_index3(
              grid->face[direction].extent, coordinate[0], coordinate[1],
              coordinate[2]);
          const size_t global_face = hdiv_global_face_index(
              direction, global[0], global[1], global[2]);
          storage->face_flux[direction][local_face] =
              hdiv_global_face_base[0][direction][global_face] +
              hdiv_global_face_flux[0][direction][global_face];
          if (!hdiv_history_ready) {
            storage->old_face_flux[direction][local_face] =
                hdiv_global_face_base[1][direction][global_face] +
                hdiv_global_face_flux[1][direction][global_face];
          }
        }
      }
    }
  }
}
#endif

static int reconcile_active_magnetic_faces(gamera_no_storage *storage,
                                           const gamera_no_grid *grid) {
  const int first_reconciliation = !hdiv_history_ready;
  const double time_tolerance =
      32.0 * DBL_EPSILON * fmax(1.0, fabs(time_sim));
  if (first_reconciliation && read_restart &&
      restart_has_yinyang_hdiv == GAMERA_YINYANG_HDIV_SCHEME_VERSION) {
    /*
     * This checkpoint already contains canonical/projected current and old CT
     * levels.  Reapplying the nonlinear donor target here would add an extra
     * operation that the uninterrupted trajectory never takes.
     */
    hdiv_history_ready = 1;
    hdiv_last_time = time_sim;
    log_info("Yin-Yang patch %d resumed canonical H(div) CT history at "
             "t=%.15g without reprojecting the checkpoint",
             patch_id, time_sim);
    return 0;
  }
  if (isfinite(hdiv_last_time) &&
      fabs(time_sim - hdiv_last_time) <= time_tolerance) {
    return 0;
  }
#if defined(GAMERA_YINYANG_HDIV_DISTRIBUTED_VERIFY) &&                     \
    defined(GAMERA_YINYANG_HDIV_DISTRIBUTED)
  const int verify_current_only_old_storage =
      hdiv_history_ready && !hdiv_distributed_verify_old_audit_done;
#endif
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  const double profile_total_start = omp_get_wtime();
  double profile_phase_start = profile_total_start;
  double profile_prepare = 0.0;
  double profile_gather = 0.0;
  double profile_assemble = 0.0;
  double profile_project = 0.0;
  double profile_scatter_recover = 0.0;
#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED
  memset(hdiv_profile_detail_current, 0,
         sizeof(hdiv_profile_detail_current));
  hdiv_profile_patch_allreduce_current = 0;
  hdiv_profile_world_allreduce_current = 0;
  hdiv_profile_math_allreduce_current = 0;
  hdiv_profile_cell_sendrecv_current = 0;
  hdiv_profile_face_sendrecv_current = 0;
  hdiv_profile_layout_sendrecv_current = 0;
  hdiv_profile_cell_send_bytes_current = 0;
  hdiv_profile_face_send_bytes_current = 0;
  memset(hdiv_profile_pcg_iterations_current, 0,
         sizeof(hdiv_profile_pcg_iterations_current));
  hdiv_profile_active = 1;
#endif
#endif
#ifdef GAMERA_YINYANG_SPARSE_OVERSET
  int local_failed = 0;
  if (!hdiv_exchange_ready || !magnetic_active_both_plan.ready ||
      !magnetic_active_current_plan.ready) {
    local_failed = initialize_hdiv_exchange(grid) != 0;
    if (world_consensus_failure(
            local_failed, "sparse H(div) receptor initialization") != 0) {
      return -1;
    }
    if (initialize_sparse_active_magnetic_plans() != 0) {
      return -1;
    }
  }
  gamera_no_sparse_plan *active_donor_plan =
#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED
      hdiv_history_ready ? &magnetic_active_current_plan
                         : &magnetic_active_both_plan;
#else
      &magnetic_active_both_plan;
#endif
  sparse_state_context sparse_context = {storage, grid};
#ifdef GAMERA_YINYANG_SPARSE_OVERSET_PROFILE
  const double sparse_profile_total_start = MPI_Wtime();
#endif
  local_failed = gamera_no_sparse_plan_exchange(
                     active_donor_plan, sparse_pack_magnetic,
                     &sparse_context) != 0;
#ifdef GAMERA_YINYANG_SPARSE_OVERSET_PROFILE
  const double sparse_profile_reconstruction_start = MPI_Wtime();
#endif
#else
  int local_failed = initialize_hdiv_exchange(grid) != 0;
#endif
  if (!local_failed) {
    local_failed = prepare_hdiv_face_targets(
                       storage, grid
#ifdef GAMERA_YINYANG_SPARSE_OVERSET
                       , active_donor_plan
#endif
                       ) != 0;
  }
  if (!local_failed) {
    local_failed = pack_hdiv_face_targets(storage, grid) != 0;
  }
  if (world_consensus_failure(local_failed, "H(div) target preparation") !=
      0) {
    return -1;
  }
#ifdef GAMERA_YINYANG_SPARSE_OVERSET_PROFILE
  sparse_profile_record(
      hdiv_history_ready ? &magnetic_active_current_profile
                         : &magnetic_active_both_profile,
      active_donor_plan, "magnetic", "pre_hdiv",
      hdiv_history_ready ? "current" : "current_old",
      MPI_Wtime() - sparse_profile_reconstruction_start,
      MPI_Wtime() - sparse_profile_total_start);
#endif
#if defined(GAMERA_YINYANG_HDIV_DISTRIBUTED_VERIFY) &&                     \
    defined(GAMERA_YINYANG_HDIV_DISTRIBUTED)
  if (verify_current_only_old_storage) {
    local_failed = hdiv_verify_capture_old_storage(storage, grid) != 0;
    if (world_consensus_failure(
            local_failed, "distributed H(div) old-level capture") != 0) {
      return -1;
    }
  }
#endif
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  profile_prepare = omp_get_wtime() - profile_phase_start;
  profile_phase_start = omp_get_wtime();
#endif
  int include_metric_weights = 1;
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
  include_metric_weights = !hdiv_cell_diagonal_cache_ready;
#endif
#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED
  local_failed = assemble_hdiv_local_faces(include_metric_weights) != 0;
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  profile_gather = omp_get_wtime() - profile_phase_start;
  profile_phase_start = omp_get_wtime();
#endif
  if (world_consensus_failure(local_failed,
                              "distributed H(div) face assembly") != 0) {
    return -1;
  }
#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED_VERIFY
  local_failed =
      hdiv_verify_distributed_snapshot(include_metric_weights) != 0;
  if (world_consensus_failure(local_failed,
                              "distributed H(div) verification oracle") !=
      0) {
    return -1;
  }
#endif
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  profile_assemble = omp_get_wtime() - profile_phase_start;
  profile_phase_start = omp_get_wtime();
#endif
#else
  int hdiv_collective_count = hdiv_send_count;
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
  if (!include_metric_weights) {
    hdiv_collective_count =
        (int)(local_active_cells *
              (GAMERA_NO_HDIV_TARGET_CHANNELS +
               GAMERA_NO_HDIV_BASE_CHANNELS));
  }
#endif
  const int hdiv_gather_status =
      MPI_Allgather(hdiv_send_buffer, hdiv_collective_count, MPI_DOUBLE,
                    hdiv_gather_buffer, hdiv_collective_count, MPI_DOUBLE,
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
                    comm_cart
#else
                    MPI_COMM_WORLD
#endif
                    );
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
  /*
   * The optimized gather is patch-local, but all later failure handling and
   * magnetic exchange collectives are world-wide.  Every rank must therefore
   * pass through a world consensus even when only one patch reports a gather
   * error; returning directly here could strand the other patch in its next
   * MPI_COMM_WORLD collective.
   */
  local_failed = hdiv_gather_status != MPI_SUCCESS;
  if (world_consensus_failure(local_failed, "H(div) face gather") != 0) {
    return -1;
  }
#else
  if (hdiv_gather_status != MPI_SUCCESS) {
    return -1;
  }
#endif
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  profile_gather = omp_get_wtime() - profile_phase_start;
  profile_phase_start = omp_get_wtime();
#endif

  local_failed =
      assemble_hdiv_global_faces(include_metric_weights,
                                 (size_t)hdiv_collective_count) != 0;
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  profile_assemble = omp_get_wtime() - profile_phase_start;
  profile_phase_start = omp_get_wtime();
#endif
#endif
  if (!local_failed) {
    local_failed = project_hdiv_level(0) != 0;
  }
  if (!local_failed && !hdiv_history_ready) {
    local_failed = project_hdiv_level(1) != 0;
  }
  if (world_consensus_failure(local_failed, "H(div) projection") != 0) {
    return -1;
  }
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  profile_project = omp_get_wtime() - profile_phase_start;
  profile_phase_start = omp_get_wtime();
#endif
#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED
  scatter_hdiv_local_faces(storage, grid);
#else
  scatter_hdiv_global_faces(storage, grid);
#endif
#if defined(GAMERA_YINYANG_HDIV_DISTRIBUTED_VERIFY) &&                     \
    defined(GAMERA_YINYANG_HDIV_DISTRIBUTED)
  if (verify_current_only_old_storage) {
    local_failed = hdiv_verify_old_storage_unchanged(storage, grid) != 0;
    if (world_consensus_failure(
            local_failed, "distributed H(div) old-level bit audit") != 0) {
      return -1;
    }
  }
#endif
  const double *current_flux[3] = {
      storage->face_flux[GAMERA_NO_I], storage->face_flux[GAMERA_NO_J],
      storage->face_flux[GAMERA_NO_K]};
  const double *old_flux[3] = {
      storage->old_face_flux[GAMERA_NO_I], storage->old_face_flux[GAMERA_NO_J],
      storage->old_face_flux[GAMERA_NO_K]};
  local_failed = gamera_no_recover_magnetic_field(
                     grid, current_flux, storage->cell_magnetic) != 0;
  if (!local_failed && !hdiv_history_ready) {
    local_failed = gamera_no_recover_magnetic_field(
                       grid, old_flux, storage->old_cell_magnetic) != 0;
  }
  if (world_consensus_failure(local_failed, "H(div) field recovery") != 0) {
    return -1;
  }
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  profile_scatter_recover = omp_get_wtime() - profile_phase_start;
  ++hdiv_profile_calls;
  if (hdiv_profile_calls > 1) {
    const double sample[GAMERA_NO_HDIV_PROFILE_PHASE_COUNT] = {
        omp_get_wtime() - profile_total_start, profile_prepare,
        profile_gather, profile_assemble, profile_project,
        profile_scatter_recover};
    for (int phase = 0; phase < GAMERA_NO_HDIV_PROFILE_PHASE_COUNT; ++phase) {
      hdiv_profile_sum[phase] += sample[phase];
      hdiv_profile_maximum[phase] =
          fmax(hdiv_profile_maximum[phase], sample[phase]);
    }
  }
#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED
  log_info("H(div) distributed profile rank=%d patch=%d call=%lu "
           "cell_halo_s=%.9g face_halo_s=%.9g reduction_s=%.9g "
           "operator_s=%.9g vector_s=%.9g correction_s=%.9g "
           "patch_allreduce=%lu world_allreduce=%lu "
           "math_allreduce=%lu collective_total=%lu "
           "cell_sendrecv=%lu face_sendrecv=%lu layout_sendrecv=%lu "
           "cell_send_bytes=%llu face_send_bytes=%llu "
           "pcg_iter_current=%lu pcg_iter_old=%lu",
           rank, patch_id, hdiv_profile_calls,
           hdiv_profile_detail_current[GAMERA_NO_HDIV_DETAIL_CELL_HALO],
           hdiv_profile_detail_current[GAMERA_NO_HDIV_DETAIL_FACE_HALO],
           hdiv_profile_detail_current[GAMERA_NO_HDIV_DETAIL_REDUCTION],
           hdiv_profile_detail_current[GAMERA_NO_HDIV_DETAIL_OPERATOR],
           hdiv_profile_detail_current[GAMERA_NO_HDIV_DETAIL_VECTOR],
           hdiv_profile_detail_current[GAMERA_NO_HDIV_DETAIL_CORRECTION],
           hdiv_profile_patch_allreduce_current,
           hdiv_profile_world_allreduce_current,
           hdiv_profile_math_allreduce_current,
           hdiv_profile_patch_allreduce_current +
               hdiv_profile_world_allreduce_current,
           hdiv_profile_cell_sendrecv_current,
           hdiv_profile_face_sendrecv_current,
           hdiv_profile_layout_sendrecv_current,
           hdiv_profile_cell_send_bytes_current,
           hdiv_profile_face_send_bytes_current,
           hdiv_profile_pcg_iterations_current[0],
           hdiv_profile_pcg_iterations_current[1]);
  hdiv_profile_active = 0;
#endif
#endif
  if (first_reconciliation) {
    log_info("Yin-Yang patch %d metric H(div) v%d initial projection: mobile "
             "max |D delta-Phi*| %.6e -> %.6e, max correction %.6e, "
             "iterations %d, Hodge=[%.6e,%.6e], min cos=%.6e",
             patch_id, GAMERA_YINYANG_HDIV_SCHEME_VERSION,
             gamera_no_yinyang_hdiv_max_before,
             gamera_no_yinyang_hdiv_max_after,
             gamera_no_yinyang_hdiv_max_correction,
             gamera_no_yinyang_hdiv_max_iterations,
             gamera_no_yinyang_hdiv_min_weight,
             gamera_no_yinyang_hdiv_max_weight,
             gamera_no_yinyang_hdiv_min_metric_cosine);
  }
  hdiv_history_ready = 1;
  hdiv_last_time = time_sim;
  return 0;
}
#endif

#ifndef GAMERA_YINYANG_SPARSE_OVERSET
static int pack_active_magnetic(const gamera_no_storage *storage,
                                const gamera_no_grid *grid) {
  if (storage == NULL || grid == NULL) {
    return -1;
  }
  for (int i = 0; i < config.ni; ++i) {
    for (int j = 0; j < config.nj; ++j) {
      for (int k = 0; k < config.nk; ++k) {
        const size_t local = local_active_index(i, j, k);
        const size_t cell = gamera_no_index3(
            grid->cell_extent, (size_t)(is + i), (size_t)(js + j),
            (size_t)(ks + k));
        const gamera_no_vec3 level[2] = {storage->cell_magnetic[cell],
                                         storage->old_cell_magnetic[cell]};
        for (int time_level = 0; time_level < 2; ++time_level) {
          for (int component = 0; component < GAMERA_NO_DIM; ++component) {
            const size_t channel =
                (size_t)time_level * GAMERA_NO_DIM + (size_t)component;
            const double value = level[time_level].value[component];
            if (!isfinite(value)) {
              return -1;
            }
            magnetic_send_buffer[channel * local_active_cells + local] =
                value;
          }
        }
      }
    }
  }
  return 0;
}

static int pack_active_fluid_state(void) {
  gamera_no_grid *grid = gamera_no_legacy_grid();
  gamera_no_storage *storage = gamera_no_legacy_storage();
  if (grid == NULL || storage == NULL) {
    return -1;
  }
  for (int i = 0; i < config.ni; ++i) {
    for (int j = 0; j < config.nj; ++j) {
      for (int k = 0; k < config.nk; ++k) {
        const size_t local = local_active_index(i, j, k);
        const size_t cell = gamera_no_index3(
            grid->cell_extent, (size_t)(is + i), (size_t)(js + j),
            (size_t)(ks + k));
#ifdef GAMERA_YINYANG_MFE_INTERFACE
        /*
         * These five wire channels deliberately carry primitive
         * [rho,vx,vy,vz,p], rather than the similarly sized conserved vector.
         * Interpolating primitives matches the reference Yin-Yang MFE fringe
         * treatment and avoids momentum/energy mixing across the interface.
         */
        gamera_no_primitive primitive[GAMERA_NO_YINYANG_TIME_LEVELS];
        const double *conserved[GAMERA_NO_YINYANG_TIME_LEVELS] = {
            &storage->conserved[cell * GAMERA_NO_FLUX_COUNT],
            &storage->old_conserved[cell * GAMERA_NO_FLUX_COUNT]};
        for (int time_level = 0;
             time_level < GAMERA_NO_YINYANG_TIME_LEVELS; ++time_level) {
          if (gamera_no_conserved_to_primitive(
                  conserved[time_level], gamma_val, rho_floor, p_floor,
                  &primitive[time_level]) != 0) {
            return -1;
          }
        }
        const double level[GAMERA_NO_YINYANG_TIME_LEVELS]
                          [GAMERA_NO_FLUX_COUNT] = {
            {primitive[0].density, primitive[0].velocity.value[0],
             primitive[0].velocity.value[1], primitive[0].velocity.value[2],
             primitive[0].pressure},
            {primitive[1].density, primitive[1].velocity.value[0],
             primitive[1].velocity.value[1], primitive[1].velocity.value[2],
             primitive[1].pressure}};
#else
        const double *level[GAMERA_NO_YINYANG_TIME_LEVELS] = {
            storage->conserved, storage->old_conserved};
#endif
        for (int time_level = 0;
             time_level < GAMERA_NO_YINYANG_TIME_LEVELS; ++time_level) {
          for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT;
               ++variable) {
            const size_t channel =
                (size_t)time_level * GAMERA_NO_FLUX_COUNT +
                (size_t)variable;
#ifdef GAMERA_YINYANG_MFE_INTERFACE
            send_buffer[channel * local_active_cells + local] =
                level[time_level][variable];
#else
            send_buffer[channel * local_active_cells + local] =
                level[time_level][cell * GAMERA_NO_FLUX_COUNT +
                                  (size_t)variable];
#endif
          }
        }
      }
    }
  }
  return 0;
}
#endif

#ifndef GAMERA_YINYANG_SPARSE_OVERSET
static int receptor_is_angular_active(const receptor_t *receptor) {
  if (receptor == NULL) {
    return 0;
  }
  const int global_j =
      proc_coords[1] * config.nj + receptor->receiver[1] - js;
  const int global_k =
      proc_coords[2] * config.nk + receptor->receiver[2] - ks;
  return global_j >= 0 && global_j < config.nj_global && global_k >= 0 &&
         global_k < config.nk_global;
}
#endif

static int exchange_hd(int include_active_receptors) {
  gamera_no_storage *storage = gamera_no_legacy_storage();
  int local_failed = storage == NULL;
#ifdef GAMERA_YINYANG_SPARSE_OVERSET
  if (!fluid_full_plan.ready || !fluid_ghost_plan.ready ||
      fluid_ghost_reference_base == NULL) {
    if (!local_failed) {
      local_failed = initialize_exchange() != 0;
    }
    if (world_consensus_failure(
            local_failed, "sparse fluid receptor initialization") != 0) {
      return -1;
    }
    if (initialize_sparse_fluid_plans() != 0) {
      return -1;
    }
  }
  gamera_no_grid *grid = gamera_no_legacy_grid();
  sparse_state_context sparse_context = {storage, grid};
  gamera_no_sparse_plan *fluid_plan =
      include_active_receptors ? &fluid_full_plan : &fluid_ghost_plan;
#ifdef GAMERA_YINYANG_SPARSE_OVERSET_PROFILE
  const double sparse_profile_total_start = MPI_Wtime();
#endif
  if (gamera_no_sparse_plan_exchange(fluid_plan, sparse_pack_fluid,
                                     &sparse_context) != 0) {
    return -1;
  }
#ifdef GAMERA_YINYANG_SPARSE_OVERSET_PROFILE
  const double sparse_profile_reconstruction_start = MPI_Wtime();
#endif
#else
  if (!local_failed) {
    local_failed = initialize_exchange() != 0;
  }
  if (!local_failed) {
    local_failed = pack_active_fluid_state() != 0;
  }
  if (world_consensus_failure(local_failed, "fluid donor packing") != 0) {
    return -1;
  }
  if (MPI_Allgather(send_buffer, send_count, MPI_DOUBLE, gather_buffer,
                    send_count, MPI_DOUBLE, MPI_COMM_WORLD) != MPI_SUCCESS) {
    return -1;
  }
#endif

  int reconstruction_failed = 0;
#pragma omp parallel for reduction(| : reconstruction_failed) schedule(static)
  for (size_t item = 0; item < receptor_count; ++item) {
    const receptor_t *receptor = &receptors[item];
    if (!include_active_receptors && receptor_is_angular_active(receptor)) {
      continue;
    }
#ifdef GAMERA_YINYANG_SPARSE_OVERSET
    if (!include_active_receptors &&
        fluid_ghost_reference_base[item] == SIZE_MAX) {
      reconstruction_failed = 1;
      continue;
    }
#endif
    double interpolated[GAMERA_NO_YINYANG_TIME_LEVELS]
                       [GAMERA_NO_FLUX_COUNT] = {{0.0}};
#ifdef GAMERA_YINYANG_SPARSE_OVERSET
    /* The sparse receive buffer stores all channels for a donor cell
     * contiguously.  Walk corners first so each slot lookup and cache line is
     * shared by all channels.  Each channel still accumulates corners in the
     * original 0..7 order, preserving the floating-point result. */
    for (int corner = 0; corner < 8; ++corner) {
      const size_t reference =
          include_active_receptors
              ? item * 8U + (size_t)corner
              : fluid_ghost_reference_base[item] + (size_t)corner;
      const size_t slot = fluid_plan->reference_slot[reference];
      const double *donor =
          &fluid_plan->receive_values[slot * (size_t)fluid_plan->channel_count];
      const double weight = receptor->weight[corner];
      for (int time_level = 0;
           time_level < GAMERA_NO_YINYANG_TIME_LEVELS; ++time_level) {
        for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
          const size_t channel =
              (size_t)time_level * GAMERA_NO_FLUX_COUNT + (size_t)variable;
          interpolated[time_level][variable] += weight * donor[channel];
        }
      }
    }
#else
    for (int time_level = 0;
         time_level < GAMERA_NO_YINYANG_TIME_LEVELS; ++time_level) {
      for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
        const size_t channel =
            (size_t)time_level * GAMERA_NO_FLUX_COUNT + (size_t)variable;
        for (int corner = 0; corner < 8; ++corner) {
          const size_t offset =
              (size_t)receptor->donor_rank[corner] * (size_t)send_count +
              channel * local_active_cells + receptor->donor_cell[corner];
          const double donor_value = gather_buffer[offset];
          interpolated[time_level][variable] +=
              receptor->weight[corner] * donor_value;
        }
      }
    }
#endif
    gamera_no_primitive primitive;
    gamera_no_primitive old_primitive;
#ifdef GAMERA_YINYANG_MFE_INTERFACE
    /* Decode the primitive wire channels, then convert both time levels back
     * to the conservative representation used by the GAMERA predictor. */
    primitive.density = interpolated[0][GAMERA_NO_FLUX_DENSITY];
    primitive.velocity.value[0] =
        interpolated[0][GAMERA_NO_FLUX_MOMENTUM_X];
    primitive.velocity.value[1] =
        interpolated[0][GAMERA_NO_FLUX_MOMENTUM_Y];
    primitive.velocity.value[2] =
        interpolated[0][GAMERA_NO_FLUX_MOMENTUM_Z];
    primitive.pressure = interpolated[0][GAMERA_NO_FLUX_ENERGY];
    old_primitive.density = interpolated[1][GAMERA_NO_FLUX_DENSITY];
    old_primitive.velocity.value[0] =
        interpolated[1][GAMERA_NO_FLUX_MOMENTUM_X];
    old_primitive.velocity.value[1] =
        interpolated[1][GAMERA_NO_FLUX_MOMENTUM_Y];
    old_primitive.velocity.value[2] =
        interpolated[1][GAMERA_NO_FLUX_MOMENTUM_Z];
    old_primitive.pressure = interpolated[1][GAMERA_NO_FLUX_ENERGY];
    double converted[GAMERA_NO_YINYANG_TIME_LEVELS][GAMERA_NO_FLUX_COUNT];
    if (gamera_no_primitive_to_conserved(&primitive, gamma_val, rho_floor,
                                         p_floor, converted[0]) != 0 ||
        gamera_no_primitive_to_conserved(&old_primitive, gamma_val, rho_floor,
                                         p_floor, converted[1]) != 0) {
      reconstruction_failed = 1;
      continue;
    }
#else
    if (gamera_no_conserved_to_primitive(
            interpolated[0], gamma_val, rho_floor, p_floor, &primitive) != 0 ||
        gamera_no_conserved_to_primitive(interpolated[1], gamma_val,
                                         rho_floor, p_floor,
                                         &old_primitive) != 0) {
      reconstruction_failed = 1;
      continue;
    }
#endif
    const int i = receptor->receiver[0];
    const int j = receptor->receiver[1];
    const int k = receptor->receiver[2];
    gas[0][gas_rho][i][j][k] = primitive.density;
    gas[0][gas_v1][i][j][k] = primitive.velocity.value[0];
    gas[0][gas_v2][i][j][k] = primitive.velocity.value[1];
    gas[0][gas_v3][i][j][k] = primitive.velocity.value[2];
    gas[0][gas_p][i][j][k] = primitive.pressure;
    gas[0][gas_p_S][i][j][k] = primitive.pressure;
    const size_t cell = gamera_no_index3(
        storage->cell_extent, (size_t)i, (size_t)j, (size_t)k);
    for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
#ifdef GAMERA_YINYANG_MFE_INTERFACE
      storage->conserved[cell * GAMERA_NO_FLUX_COUNT + (size_t)variable] =
          converted[0][variable];
      storage->old_conserved[cell * GAMERA_NO_FLUX_COUNT + (size_t)variable] =
          converted[1][variable];
#else
      storage->old_conserved[cell * GAMERA_NO_FLUX_COUNT +
                             (size_t)variable] = interpolated[1][variable];
#endif
    }
  }
  local_failed |= reconstruction_failed;
  const int reconstruction_status = world_consensus_failure(
      local_failed, "fluid donor reconstruction");
#ifdef GAMERA_YINYANG_SPARSE_OVERSET_PROFILE
  if (reconstruction_status == 0) {
    sparse_profile_record(
        include_active_receptors ? &fluid_full_profile
                                 : &fluid_ghost_profile,
        fluid_plan, "fluid", include_active_receptors ? "full" : "ghost",
        "current_old", MPI_Wtime() - sparse_profile_reconstruction_start,
        MPI_Wtime() - sparse_profile_total_start);
  }
#endif
  return reconstruction_status;
}

int gamera_no_yinyang_exchange_hd(void) { return exchange_hd(1); }

int gamera_no_yinyang_exchange_hd_ghosts_only(void) {
  return exchange_hd(0);
}

int gamera_no_yinyang_sync_edge_emf(
    gamera_no_storage *storage, const gamera_no_grid *grid,
    const size_t active_lower[3], const size_t active_upper[3],
    void *context) {
  (void)context;
  int local_failed = storage == NULL || grid == NULL;
#ifdef GAMERA_YINYANG_SPARSE_OVERSET
  if (!emf_plan.ready) {
    if (!local_failed &&
        initialize_emf_exchange(grid, active_lower, active_upper) != 0) {
      log_error("Yin-Yang CT receptor-edge initialization failed on patch %d",
                patch_id);
      local_failed = 1;
    }
    if (world_consensus_failure(
            local_failed, "sparse EMF receptor initialization") != 0) {
      return -1;
    }
    if (initialize_sparse_edge_plan() != 0) {
      return -1;
    }
  }
  sparse_state_context sparse_context = {storage, grid};
#ifdef GAMERA_YINYANG_SPARSE_OVERSET_PROFILE
  const double sparse_profile_total_start = MPI_Wtime();
#endif
  if (gamera_no_sparse_plan_exchange(&emf_plan, sparse_pack_electric,
                                     &sparse_context) != 0) {
    return -1;
  }
#ifdef GAMERA_YINYANG_SPARSE_OVERSET_PROFILE
  const double sparse_profile_reconstruction_start = MPI_Wtime();
#endif
#else
  if (!local_failed &&
      initialize_emf_exchange(grid, active_lower, active_upper) != 0) {
    log_error("Yin-Yang CT receptor-edge initialization failed on patch %d",
              patch_id);
    local_failed = 1;
  }
  if (!local_failed) {
    local_failed = pack_active_electric(storage, grid) != 0;
  }
  if (world_consensus_failure(local_failed, "edge-EMF donor packing") != 0) {
    return -1;
  }
  if (MPI_Allgather(emf_send_buffer, emf_send_count, MPI_DOUBLE,
                    emf_gather_buffer, emf_send_count, MPI_DOUBLE,
                    MPI_COMM_WORLD) != MPI_SUCCESS) {
    return -1;
  }
#endif

  for (size_t item = 0; item < edge_receptor_count; ++item) {
    const edge_receptor_t *receptor = &edge_receptors[item];
    gamera_no_vec3 electric = {{0.0, 0.0, 0.0}};
    for (int component = 0; component < GAMERA_NO_DIM; ++component) {
      for (int corner = 0; corner < 8; ++corner) {
#ifdef GAMERA_YINYANG_SPARSE_OVERSET
        const double donor_value = gamera_no_sparse_plan_value(
            &emf_plan, item * 8U + (size_t)corner, component);
#else
        const size_t offset =
            (size_t)receptor->donor_rank[corner] * (size_t)emf_send_count +
            (size_t)component * local_active_cells +
            receptor->donor_cell[corner];
        const double donor_value = emf_gather_buffer[offset];
#endif
        electric.value[component] +=
            receptor->weight[corner] * donor_value;
      }
    }
    const gamera_no_edge_geometry *geometry =
        &grid->edge[receptor->direction].value[receptor->receiver_edge];
    const double line_integral =
        geometry->length * vector_dot(electric, geometry->normal);
    if (!isfinite(line_integral)) {
      local_failed = 1;
      break;
    }
    double *receiver =
        &storage->edge_emf[receptor->direction][receptor->receiver_edge];
    *receiver = (1.0 - receptor->blend_weight) * (*receiver) +
                receptor->blend_weight * line_integral;
  }
  const int reconstruction_status = world_consensus_failure(
      local_failed, "edge-EMF donor reconstruction");
#ifdef GAMERA_YINYANG_SPARSE_OVERSET_PROFILE
  if (reconstruction_status == 0) {
    sparse_profile_record(
        &emf_profile, &emf_plan, "emf", "edge", "current",
        MPI_Wtime() - sparse_profile_reconstruction_start,
        MPI_Wtime() - sparse_profile_total_start);
  }
#endif
  return reconstruction_status;
}

int gamera_no_yinyang_exchange_magnetic_ghosts(void) {
  gamera_no_grid *grid = gamera_no_legacy_grid();
  gamera_no_storage *storage = gamera_no_legacy_storage();
  int local_failed = grid == NULL || storage == NULL;
#ifndef GAMERA_YINYANG_SPARSE_OVERSET
  if (world_consensus_failure(local_failed, "magnetic storage lookup") != 0) {
    return -1;
  }

  local_failed = initialize_magnetic_exchange(grid) != 0;
#else
  if (!magnetic_exchange_ready) {
    if (!local_failed) {
      local_failed = initialize_magnetic_exchange(grid) != 0;
    }
    if (world_consensus_failure(
            local_failed, "sparse magnetic receptor initialization") != 0) {
      return -1;
    }
  }
#endif
#ifdef GAMERA_YINYANG_SPARSE_OVERSET
#ifdef GAMERA_YINYANG_HDIV_RECONCILE
  if (reconcile_active_magnetic_faces(storage, grid) != 0) {
    return -1;
  }
#endif
  if (!magnetic_ghost_plan.ready) {
    if (initialize_sparse_magnetic_ghost_plan() != 0) {
      return -1;
    }
  }
  sparse_state_context sparse_context = {storage, grid};
#ifdef GAMERA_YINYANG_SPARSE_OVERSET_PROFILE
  const double sparse_profile_total_start = MPI_Wtime();
#endif
  if (gamera_no_sparse_plan_exchange(&magnetic_ghost_plan,
                                     sparse_pack_magnetic,
                                     &sparse_context) != 0) {
    return -1;
  }
#ifdef GAMERA_YINYANG_SPARSE_OVERSET_PROFILE
  const double sparse_profile_reconstruction_start = MPI_Wtime();
#endif
#else
  if (!local_failed) {
    local_failed = pack_active_magnetic(storage, grid) != 0;
  }
  if (world_consensus_failure(local_failed, "magnetic donor packing") != 0) {
    return -1;
  }
  if (MPI_Allgather(magnetic_send_buffer, magnetic_send_count, MPI_DOUBLE,
                    magnetic_gather_buffer, magnetic_send_count, MPI_DOUBLE,
                    MPI_COMM_WORLD) != MPI_SUCCESS) {
    return -1;
  }
#ifdef GAMERA_YINYANG_HDIV_RECONCILE
  /*
   * The first gather supplies donor active B for the face targets.  Reconcile
   * active fluxes, recover B, then gather again so angular ghosts see the
   * post-projection state rather than the stale CT branch.
   */
  if (reconcile_active_magnetic_faces(storage, grid) != 0) {
    return -1;
  }
  local_failed = pack_active_magnetic(storage, grid) != 0;
  if (world_consensus_failure(local_failed,
                              "post-H(div) magnetic donor packing") != 0) {
    return -1;
  }
  if (MPI_Allgather(magnetic_send_buffer, magnetic_send_count, MPI_DOUBLE,
                    magnetic_gather_buffer, magnetic_send_count, MPI_DOUBLE,
                    MPI_COMM_WORLD) != MPI_SUCCESS) {
    return -1;
  }
#endif
#endif

#pragma omp parallel for schedule(static)
  for (size_t item = 0; item < face_receptor_count; ++item) {
    const face_receptor_t *receptor = &face_receptors[item];
    gamera_no_vec3 magnetic[2] = {{{0.0}}};
#ifdef GAMERA_YINYANG_SPARSE_OVERSET
    /* As above, traverse the channel-contiguous donor slot once per corner
     * while retaining the original corner summation order per component. */
    for (int corner = 0; corner < 8; ++corner) {
      const size_t reference = item * 8U + (size_t)corner;
      const size_t slot = magnetic_ghost_plan.reference_slot[reference];
      const double *donor = &magnetic_ghost_plan.receive_values[
          slot * (size_t)magnetic_ghost_plan.channel_count];
      const double weight = receptor->weight[corner];
      for (int time_level = 0; time_level < 2; ++time_level) {
        for (int component = 0; component < GAMERA_NO_DIM; ++component) {
          const size_t channel =
              (size_t)time_level * GAMERA_NO_DIM + (size_t)component;
          magnetic[time_level].value[component] += weight * donor[channel];
        }
      }
    }
#else
    for (int time_level = 0; time_level < 2; ++time_level) {
      for (int component = 0; component < GAMERA_NO_DIM; ++component) {
        const size_t channel =
            (size_t)time_level * GAMERA_NO_DIM + (size_t)component;
        for (int corner = 0; corner < 8; ++corner) {
          const size_t offset =
              (size_t)receptor->donor_rank[corner] *
                  (size_t)magnetic_send_count +
              channel * local_active_cells + receptor->donor_cell[corner];
          const double donor_value = magnetic_gather_buffer[offset];
          magnetic[time_level].value[component] +=
              receptor->weight[corner] * donor_value;
        }
      }
    }
#endif
    const gamera_no_face_geometry *geometry =
        &grid->face[receptor->direction].value[receptor->receiver_face];
    storage->face_flux[receptor->direction][receptor->receiver_face] =
        vector_dot(magnetic[0], geometry->area_vector);
    storage->old_face_flux[receptor->direction][receptor->receiver_face] =
        vector_dot(magnetic[1], geometry->area_vector);
  }
  const double *current_flux[3] = {
      storage->face_flux[GAMERA_NO_I], storage->face_flux[GAMERA_NO_J],
      storage->face_flux[GAMERA_NO_K]};
  const double *old_flux[3] = {
      storage->old_face_flux[GAMERA_NO_I], storage->old_face_flux[GAMERA_NO_J],
      storage->old_face_flux[GAMERA_NO_K]};
  local_failed = recover_magnetic_receptor_cells(
                     grid, current_flux, storage->cell_magnetic) != 0;
  if (!local_failed) {
    local_failed = recover_magnetic_receptor_cells(
                       grid, old_flux, storage->old_cell_magnetic) != 0;
  }
  const int reconstruction_status = world_consensus_failure(
      local_failed, "magnetic ghost field recovery");
#ifdef GAMERA_YINYANG_SPARSE_OVERSET_PROFILE
  if (reconstruction_status == 0) {
    sparse_profile_record(
        &magnetic_ghost_profile, &magnetic_ghost_plan, "magnetic",
#ifdef GAMERA_YINYANG_HDIV_RECONCILE
        "post_hdiv",
#else
        "ghost",
#endif
        "current_old",
        MPI_Wtime() - sparse_profile_reconstruction_start,
        MPI_Wtime() - sparse_profile_total_start);
  }
#endif
  return reconstruction_status;
}

static void yinyang_exchange_destroy_local(void) {
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  if (hdiv_profile_calls > 1) {
    const double inverse_samples = 1.0 / (double)(hdiv_profile_calls - 1);
    log_info("H(div) profile summary rank=%d patch=%d samples=%lu "
             "warmup_excluded=1 total_mean=%.9g total_max=%.9g "
             "prepare_mean=%.9g prepare_max=%.9g "
             "gather_mean=%.9g gather_max=%.9g "
             "assemble_mean=%.9g assemble_max=%.9g "
             "project_mean=%.9g project_max=%.9g "
             "scatter_recover_mean=%.9g scatter_recover_max=%.9g",
             rank, patch_id, hdiv_profile_calls - 1,
             hdiv_profile_sum[GAMERA_NO_HDIV_PROFILE_TOTAL] * inverse_samples,
             hdiv_profile_maximum[GAMERA_NO_HDIV_PROFILE_TOTAL],
             hdiv_profile_sum[GAMERA_NO_HDIV_PROFILE_PREPARE] *
                 inverse_samples,
             hdiv_profile_maximum[GAMERA_NO_HDIV_PROFILE_PREPARE],
             hdiv_profile_sum[GAMERA_NO_HDIV_PROFILE_GATHER] * inverse_samples,
             hdiv_profile_maximum[GAMERA_NO_HDIV_PROFILE_GATHER],
             hdiv_profile_sum[GAMERA_NO_HDIV_PROFILE_ASSEMBLE] *
                 inverse_samples,
             hdiv_profile_maximum[GAMERA_NO_HDIV_PROFILE_ASSEMBLE],
             hdiv_profile_sum[GAMERA_NO_HDIV_PROFILE_PROJECT] * inverse_samples,
             hdiv_profile_maximum[GAMERA_NO_HDIV_PROFILE_PROJECT],
             hdiv_profile_sum[GAMERA_NO_HDIV_PROFILE_RECOVER] * inverse_samples,
             hdiv_profile_maximum[GAMERA_NO_HDIV_PROFILE_RECOVER]);
  }
#endif
  free(receptors);
#ifndef GAMERA_YINYANG_SPARSE_OVERSET
  free(send_buffer);
  free(gather_buffer);
#else
  free(fluid_ghost_reference_base);
#endif
  free(edge_receptors);
#ifndef GAMERA_YINYANG_SPARSE_OVERSET
  free(emf_send_buffer);
  free(emf_gather_buffer);
#endif
  free(face_receptors);
  free(magnetic_recovery_cells);
#ifndef GAMERA_YINYANG_SPARSE_OVERSET
  free(magnetic_send_buffer);
  free(magnetic_gather_buffer);
#endif
#ifdef GAMERA_YINYANG_HDIV_RECONCILE
  free(active_face_receptors);
  free(hdiv_send_buffer);
#ifndef GAMERA_YINYANG_HDIV_DISTRIBUTED
  free(hdiv_gather_buffer);
#else
  free(hdiv_face_halo_send);
  free(hdiv_face_halo_receive);
  free(hdiv_cell_halo_send);
  free(hdiv_cell_halo_receive);
#endif
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    for (int time_level = 0;
         time_level < GAMERA_NO_YINYANG_TIME_LEVELS; ++time_level) {
      free(hdiv_target_face_flux[time_level][direction]);
#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED
      free(hdiv_local_face_base[time_level][direction]);
      free(hdiv_local_face_flux[time_level][direction]);
      hdiv_local_face_base[time_level][direction] = NULL;
      hdiv_local_face_flux[time_level][direction] = NULL;
#else
      free(hdiv_global_face_base[time_level][direction]);
      free(hdiv_global_face_flux[time_level][direction]);
#endif
      hdiv_target_face_flux[time_level][direction] = NULL;
#ifndef GAMERA_YINYANG_HDIV_DISTRIBUTED
      hdiv_global_face_base[time_level][direction] = NULL;
      hdiv_global_face_flux[time_level][direction] = NULL;
#endif
    }
#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED
    free(hdiv_local_face_mobility[direction]);
    hdiv_local_face_mobility[direction] = NULL;
#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED_VERIFY
    free(hdiv_distributed_verify_old_face[direction]);
    hdiv_distributed_verify_old_face[direction] = NULL;
#endif
#else
    free(hdiv_global_face_mobility[direction]);
    hdiv_global_face_mobility[direction] = NULL;
#endif
  }
  free(hdiv_lambda);
  free(hdiv_residual);
  free(hdiv_preconditioned);
  free(hdiv_search);
  free(hdiv_operator_search);
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
  free(hdiv_cell_diagonal_cache);
  free(hdiv_dot_block_partial);
#endif
#endif
  receptors = NULL;
#ifndef GAMERA_YINYANG_SPARSE_OVERSET
  send_buffer = NULL;
  gather_buffer = NULL;
#else
  fluid_ghost_reference_base = NULL;
#endif
  edge_receptors = NULL;
#ifndef GAMERA_YINYANG_SPARSE_OVERSET
  emf_send_buffer = NULL;
  emf_gather_buffer = NULL;
#endif
  face_receptors = NULL;
  magnetic_recovery_cells = NULL;
#ifndef GAMERA_YINYANG_SPARSE_OVERSET
  magnetic_send_buffer = NULL;
  magnetic_gather_buffer = NULL;
#endif
#ifdef GAMERA_YINYANG_HDIV_RECONCILE
  active_face_receptors = NULL;
  hdiv_send_buffer = NULL;
#ifndef GAMERA_YINYANG_HDIV_DISTRIBUTED
  hdiv_gather_buffer = NULL;
#else
  hdiv_face_halo_send = NULL;
  hdiv_face_halo_receive = NULL;
  hdiv_face_halo_capacity = 0;
  hdiv_cell_halo_send = NULL;
  hdiv_cell_halo_receive = NULL;
  hdiv_cell_halo_capacity = 0;
  hdiv_distributed_layout_verified = 0;
#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED_VERIFY
  hdiv_distributed_verify_ready = 0;
  hdiv_distributed_verify_sentinel_ready = 0;
  hdiv_distributed_verify_old_audit_done = 0;
  hdiv_distributed_verify_projection_mask = 0U;
  memset(hdiv_distributed_verify_max_face, 0,
         sizeof(hdiv_distributed_verify_max_face));
  memset(hdiv_distributed_verify_max_before, 0,
         sizeof(hdiv_distributed_verify_max_before));
  memset(hdiv_distributed_verify_tolerance, 0,
         sizeof(hdiv_distributed_verify_tolerance));
#endif
  memset(hdiv_local_face_extent, 0, sizeof(hdiv_local_face_extent));
  memset(hdiv_local_cell_extent, 0, sizeof(hdiv_local_cell_extent));
#endif
  hdiv_lambda = NULL;
  hdiv_residual = NULL;
  hdiv_preconditioned = NULL;
  hdiv_search = NULL;
  hdiv_operator_search = NULL;
#ifdef GAMERA_YINYANG_HDIV_OPTIMIZED
  hdiv_cell_diagonal_cache = NULL;
  hdiv_cell_diagonal_cache_ready = 0;
  hdiv_dot_block_partial = NULL;
  hdiv_dot_block_count = 0;
#endif
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  hdiv_profile_calls = 0;
  memset(hdiv_profile_sum, 0, sizeof(hdiv_profile_sum));
  memset(hdiv_profile_maximum, 0, sizeof(hdiv_profile_maximum));
#ifdef GAMERA_YINYANG_HDIV_DISTRIBUTED
  hdiv_profile_active = 0;
  memset(hdiv_profile_detail_current, 0,
         sizeof(hdiv_profile_detail_current));
  hdiv_profile_patch_allreduce_current = 0;
  hdiv_profile_world_allreduce_current = 0;
  hdiv_profile_math_allreduce_current = 0;
  hdiv_profile_cell_sendrecv_current = 0;
  hdiv_profile_face_sendrecv_current = 0;
  hdiv_profile_layout_sendrecv_current = 0;
  hdiv_profile_cell_send_bytes_current = 0;
  hdiv_profile_face_send_bytes_current = 0;
  memset(hdiv_profile_pcg_iterations_current, 0,
         sizeof(hdiv_profile_pcg_iterations_current));
#endif
#endif
#endif
  receptor_count = 0;
  local_active_cells = 0;
#ifndef GAMERA_YINYANG_SPARSE_OVERSET
  send_count = 0;
#endif
  exchange_ready = 0;
  edge_receptor_count = 0;
#ifndef GAMERA_YINYANG_SPARSE_OVERSET
  emf_send_count = 0;
#endif
  emf_exchange_ready = 0;
  face_receptor_count = 0;
  magnetic_recovery_cell_count = 0;
#ifndef GAMERA_YINYANG_SPARSE_OVERSET
  magnetic_send_count = 0;
#endif
  magnetic_exchange_ready = 0;
#if defined(GAMERA_YINYANG_SPARSE_OVERSET) &&                              \
    defined(GAMERA_YINYANG_SPARSE_OVERSET_PROFILE)
  memset(&fluid_full_profile, 0, sizeof(fluid_full_profile));
  memset(&fluid_ghost_profile, 0, sizeof(fluid_ghost_profile));
  memset(&emf_profile, 0, sizeof(emf_profile));
  memset(&magnetic_ghost_profile, 0, sizeof(magnetic_ghost_profile));
#endif
#ifdef GAMERA_YINYANG_HDIV_RECONCILE
  active_face_receptor_count = 0;
  hdiv_send_count = 0;
  hdiv_exchange_ready = 0;
  hdiv_history_ready = 0;
  hdiv_last_time = NAN;
#if defined(GAMERA_YINYANG_SPARSE_OVERSET) &&                              \
    defined(GAMERA_YINYANG_SPARSE_OVERSET_PROFILE)
  memset(&magnetic_active_both_profile, 0,
         sizeof(magnetic_active_both_profile));
  memset(&magnetic_active_current_profile, 0,
         sizeof(magnetic_active_current_profile));
#endif
#endif
  gamera_no_yinyang_receptor_count = 0;
  gamera_no_yinyang_active_receptor_count = 0;
  gamera_no_yinyang_max_donor_extrapolation = 0.0;
  gamera_no_yinyang_edge_receptor_count = 0;
  gamera_no_yinyang_max_emf_donor_extrapolation = 0.0;
  gamera_no_yinyang_magnetic_face_receptor_count = 0;
  gamera_no_yinyang_max_magnetic_donor_extrapolation = 0.0;
  gamera_no_yinyang_active_magnetic_face_receptor_count = 0;
  gamera_no_yinyang_hdiv_max_before = 0.0;
  gamera_no_yinyang_hdiv_max_after = 0.0;
  gamera_no_yinyang_hdiv_max_correction = 0.0;
  gamera_no_yinyang_hdiv_max_iterations = 0;
  gamera_no_yinyang_hdiv_min_weight = DBL_MAX;
  gamera_no_yinyang_hdiv_max_weight = 0.0;
  gamera_no_yinyang_hdiv_min_metric_cosine = 1.0;
}

void gamera_no_yinyang_exchange_destroy(void) {
#ifdef GAMERA_YINYANG_SPARSE_OVERSET
#ifdef GAMERA_YINYANG_SPARSE_OVERSET_PROFILE
  sparse_profile_summary(&fluid_full_profile, &fluid_full_plan);
  sparse_profile_summary(&fluid_ghost_profile, &fluid_ghost_plan);
  sparse_profile_summary(&emf_profile, &emf_plan);
  sparse_profile_summary(&magnetic_ghost_profile, &magnetic_ghost_plan);
#ifdef GAMERA_YINYANG_HDIV_RECONCILE
  sparse_profile_summary(&magnetic_active_both_profile,
                         &magnetic_active_both_plan);
  sparse_profile_summary(&magnetic_active_current_profile,
                         &magnetic_active_current_plan);
#endif
#endif
  gamera_no_sparse_plan_destroy(&fluid_full_plan);
  gamera_no_sparse_plan_destroy(&fluid_ghost_plan);
  gamera_no_sparse_plan_destroy(&emf_plan);
  gamera_no_sparse_plan_destroy(&magnetic_ghost_plan);
#ifdef GAMERA_YINYANG_HDIV_RECONCILE
  gamera_no_sparse_plan_destroy(&magnetic_active_both_plan);
  gamera_no_sparse_plan_destroy(&magnetic_active_current_plan);
#endif
#endif
  yinyang_exchange_destroy_local();
}
