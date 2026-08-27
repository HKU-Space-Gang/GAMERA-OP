#include "nonorthogonal_driver.h"

#include "nonorthogonal_advance.h"
#include "nonorthogonal_background.h"
#include "nonorthogonal_legacy_adapter.h"
#include "nonorthogonal_step.h"
#ifdef GAMERA_NONORTHOGONAL_MATCHED_HPLUS_FEEDBACK
#include "nonorthogonal_inner_magnetosphere_runtime.h"
#endif
#ifdef GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC
#include "nonorthogonal_inner_magnetosphere_live_driver.h"
#include "nonorthogonal_inner_magnetosphere_online_runtime.h"
#endif
#ifdef GAMERA_MI_COUPLING
#include "nonorthogonal_mi_coupling.h"
#include "nonorthogonal_mi_restart.h"
#endif

#include "config.h"
#include "analysis_io.h"
#include "log.h"
#include "problem.h"
#include "setup_mpi.h"
#include "utils.h"
#ifdef GAMERA_YINYANG_BACKEND
#include "nonorthogonal_yinyang_exchange.h"
#endif

#include <math.h>
#include <mpi.h>
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#if defined(GAMERA_YINYANG_HDIV_PROFILE) || defined(GAMERA_BENCH_PROFILE_JSONL)
#include <omp.h>
#endif
#ifdef GAMERA_BENCH_PROFILE_JSONL
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct gamera_bench_step_sample {
  double step_excluding_checkpoint_io;
  double mhd_excluding_mi_checkpoint_io;
  double mi_in_step;
  int mi_updated;
} gamera_bench_step_sample;

enum { GAMERA_BENCH_MAX_STEPS = 4096 };

static int gamera_bench_env_steps(const char *name, int *result) {
  const char *text = getenv(name);
  if (text == NULL || text[0] == '\0' || result == NULL) {
    return -1;
  }
  errno = 0;
  char *end = NULL;
  const unsigned long value = strtoul(text, &end, 10);
  if (errno != 0 || end == text || *end != '\0' || value > (unsigned long)INT_MAX) {
    return -1;
  }
  *result = (int)value;
  return 0;
}

static int gamera_bench_safe_text(const char *text, int require_sha256) {
  if (text == NULL || text[0] == '\0') {
    return 0;
  }
  size_t length = 0;
  for (; text[length] != '\0'; ++length) {
    const unsigned char value = (unsigned char)text[length];
    const int safe = (value >= 'a' && value <= 'z') ||
                     (value >= 'A' && value <= 'Z') ||
                     (value >= '0' && value <= '9') || value == '_' ||
                     value == '-' || value == '.';
    if (!safe || length >= 255U) {
      return 0;
    }
  }
  if (require_sha256 && length != 64U) {
    return 0;
  }
  return 1;
}

static int gamera_bench_consensus(int local_error) {
  int global_error = 0;
  if (MPI_Allreduce(&local_error, &global_error, 1, MPI_INT, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS) {
    MPI_Abort(MPI_COMM_WORLD, 91);
    return 1;
  }
  return global_error != 0;
}

static int gamera_bench_write_profile(const gamera_bench_step_sample *samples,
                                      int total_steps, int warmup_steps,
                                      int timed_steps, double initial_mi_s,
                                      const char *run_id,
                                      const char *manifest_sha256,
                                      const char *plan_sha256,
                                      const char *binary_sha256,
                                      const char *build_manifest_sha256,
                                      int omp_team_threads,
                                      double final_time_code,
                                      double final_time_physical_s,
                                      double final_max_divB_total,
                                      int force_failure) {
  if (samples == NULL || total_steps <= 0 || warmup_steps < 0 ||
      timed_steps <= 0 || warmup_steps + timed_steps != total_steps ||
      !isfinite(final_time_code) || !isfinite(final_time_physical_s) ||
      !isfinite(final_max_divB_total) || final_time_code < 0.0 ||
      final_time_physical_s < 0.0 || final_max_divB_total < 0.0) {
    return -1;
  }
  int world_size = 0;
  if (MPI_Comm_size(MPI_COMM_WORLD, &world_size) != MPI_SUCCESS ||
      rank < 0 || rank >= world_size) {
    return -1;
  }
  char path[64];
  const int path_length =
      snprintf(path, sizeof(path), "profile.rank%08d.jsonl", rank);
  if (path_length <= 0 || (size_t)path_length >= sizeof(path)) {
    return -1;
  }
  FILE *stream = force_failure ? NULL : fopen(path, "w");
  if (stream == NULL) {
    return -1;
  }
  int status = 0;
  int timed_mi_updates = 0;
  for (int index = warmup_steps; index < total_steps; ++index) {
    timed_mi_updates += samples[index].mi_updated;
  }
  if (fprintf(stream,
              "{\"record\":\"run\",\"rank\":%d,\"world_size\":%d,"
              "\"warmup_steps\":%d,\"timed_steps\":%d,\"total_steps\":%d,"
              "\"primary_section\":\"step_excluding_checkpoint_io\","
              "\"initial_mi_s\":%.17g,\"timed_mi_updates\":%d,"
              "\"run_id\":\"%s\",\"manifest_sha256\":\"%s\","
              "\"plan_sha256\":\"%s\",\"binary_sha256\":\"%s\","
              "\"build_manifest_sha256\":\"%s\","
              "\"omp_max_threads\":%d,\"omp_dynamic\":%d,"
              "\"omp_team_threads\":%d,"
              "\"final_completed_steps\":%d,"
              "\"final_time_code\":%.17g,"
              "\"final_time_physical_s\":%.17g,"
              "\"final_max_divB_total\":%.17g,"
              "\"patch_id\":%d,\"patch_count\":%d,\"patch_rank\":%d,"
              "\"patch_size\":%d,\"proc_dims\":[%d,%d,%d],"
              "\"proc_coords\":[%d,%d,%d],\"global_grid\":[%d,%d,%d],"
              "\"local_grid\":[%d,%d,%d]}\n",
              rank, world_size, warmup_steps, timed_steps, total_steps,
              initial_mi_s, timed_mi_updates, run_id, manifest_sha256,
              plan_sha256, binary_sha256, build_manifest_sha256,
              omp_get_max_threads(), omp_get_dynamic(), omp_team_threads,
              total_steps, final_time_code, final_time_physical_s,
              final_max_divB_total,
              patch_id, patch_count, patch_rank, patch_size,
              config.proc_dims[0], config.proc_dims[1],
              config.proc_dims[2], proc_coords[0], proc_coords[1],
              proc_coords[2], config.ni_global, config.nj_global,
              config.nk_global, config.ni, config.nj, config.nk) < 0) {
    status = -1;
  }
  for (int index = 0; status == 0 && index < total_steps; ++index) {
    if (fprintf(stream,
                "{\"record\":\"step\",\"rank\":%d,\"step\":%d,"
                "\"sections\":{\"step_excluding_checkpoint_io\":%.17g,"
                "\"mhd_excluding_mi_checkpoint_io\":%.17g,"
                "\"mi_in_step\":%.17g,\"checkpoint_io_in_timed_section\":0.0}}\n",
                rank, index + 1,
                samples[index].step_excluding_checkpoint_io,
                samples[index].mhd_excluding_mi_checkpoint_io,
                samples[index].mi_in_step) < 0) {
      status = -1;
    }
  }
  if (fclose(stream) != 0) {
    status = -1;
  }
  return status;
}
#endif

static const size_t active_lower[3] = {NG, NG, NG};
static void active_upper_values(size_t upper[3]);

#ifdef GAMERA_NONORTHOGONAL_HAS_FLUID_FLUX_BOUNDARY
static int synchronize_fluid_flux(
    gamera_no_storage *storage, const gamera_no_grid *grid,
    const size_t lower[3], const size_t upper[3], void *context) {
  (void)context;
  return problem_nonorthogonal_fluid_flux_boundary(storage, grid, lower,
                                                    upper);
}
#endif

#ifdef GAMERA_NONORTHOGONAL_HAS_BACKGROUND_FIELD
static gamera_no_background_data driver_background;
static int driver_background_ready;
static double total_divB_max;

static int update_total_divergence_diagnostic(void) {
  if (!driver_background_ready) {
    return -1;
  }
  const gamera_no_grid *grid = gamera_no_legacy_grid();
  const gamera_no_storage *storage = gamera_no_legacy_storage();
  if (grid == NULL || storage == NULL) {
    return -1;
  }
  size_t upper[3];
  active_upper_values(upper);
  double local_maximum = 0.0;
  for (size_t i = active_lower[0]; i < upper[0]; ++i) {
    for (size_t j = active_lower[1]; j < upper[1]; ++j) {
      for (size_t k = active_lower[2]; k < upper[2]; ++k) {
        double net_flux = 0.0;
        for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
          size_t upper_face_coordinate[3] = {i, j, k};
          ++upper_face_coordinate[direction];
          const size_t lower_face = gamera_no_index3(
              grid->face[direction].extent, i, j, k);
          const size_t upper_face = gamera_no_index3(
              grid->face[direction].extent, upper_face_coordinate[0],
              upper_face_coordinate[1], upper_face_coordinate[2]);
          net_flux +=
              storage->face_flux[direction][upper_face] +
              driver_background.face_flux[direction][upper_face] -
              storage->face_flux[direction][lower_face] -
              driver_background.face_flux[direction][lower_face];
        }
        local_maximum = fmax(local_maximum, fabs(net_flux));
      }
    }
  }
  return MPI_Allreduce(&local_maximum, &total_divB_max, 1, MPI_DOUBLE,
                       MPI_MAX, MPI_COMM_WORLD) == MPI_SUCCESS
             ? 0
             : -1;
}

const gamera_no_background_field *gamera_no_driver_background_field(void) {
  return driver_background_ready ? &driver_background.field : NULL;
}
#endif

int gamera_no_driver_prepare(void) {
#ifdef GAMERA_NONORTHOGONAL_HAS_BACKGROUND_FIELD
  gamera_no_grid *grid = gamera_no_legacy_grid();
  if (grid == NULL || driver_background_ready ||
      gamera_no_background_create(
          grid, problem_nonorthogonal_background_field, NULL,
          &driver_background) != 0) {
    return -1;
  }
#ifdef GAMERA_NONORTHOGONAL_HAS_BACKGROUND_ADJUST
  if (problem_nonorthogonal_adjust_background(grid, &driver_background) !=
      0) {
    gamera_no_background_destroy(&driver_background);
    return -1;
  }
#endif
  driver_background_ready = 1;

#ifdef GAMERA_NONORTHOGONAL_HAS_BACKGROUND_RESIDUAL_INIT
  /* Checkpoints already store the evolved residual CT flux B1.  Repeating the
   * total-to-residual conversion on restart subtracts B0 a second time and
   * destroys div(B0+B1). */
  if (!read_restart &&
      problem_nonorthogonal_initialize_background_residual(
          grid, gamera_no_legacy_storage(), &driver_background) != 0) {
    gamera_no_driver_finalize();
    return -1;
  }
#endif

  if (update_total_divergence_diagnostic() != 0) {
    gamera_no_driver_finalize();
    return -1;
  }

  double local_divergence = 0.0;
  double local_force = 0.0;
  size_t upper[3];
  active_upper_values(upper);
  for (size_t i = active_lower[0]; i < upper[0]; ++i) {
    for (size_t j = active_lower[1]; j < upper[1]; ++j) {
      for (size_t k = active_lower[2]; k < upper[2]; ++k) {
        const size_t cell = gamera_no_index3(grid->cell_extent, i, j, k);
        double net_flux = 0.0;
        for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
          size_t upper_face_coordinate[3] = {i, j, k};
          ++upper_face_coordinate[direction];
          const size_t lower_face = gamera_no_index3(
              grid->face[direction].extent, i, j, k);
          const size_t upper_face = gamera_no_index3(
              grid->face[direction].extent, upper_face_coordinate[0],
              upper_face_coordinate[1], upper_face_coordinate[2]);
          net_flux += driver_background.face_flux[direction][upper_face] -
                      driver_background.face_flux[direction][lower_face];
        }
        local_divergence = fmax(local_divergence, fabs(net_flux));
        double force_square = 0.0;
        for (int component = 0; component < GAMERA_NO_DIM; ++component) {
          const double value =
              driver_background.cell_force[cell].value[component];
          force_square += value * value;
        }
        local_force = fmax(local_force, sqrt(force_square));
      }
    }
  }
  double global_divergence = 0.0;
  double global_force = 0.0;
  if (MPI_Allreduce(&local_divergence, &global_divergence, 1, MPI_DOUBLE,
                    MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Allreduce(&local_force, &global_force, 1, MPI_DOUBLE, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS) {
    gamera_no_driver_finalize();
    return -1;
  }
  log_info("Prepared 12-point background field on patch %d: "
           "max integrated divB0 %.4e, max |dpB0| %.4e",
           patch_id, global_divergence, global_force);
#endif
#ifdef GAMERA_MI_COUPLING
  if (gamera_mi_coupling_prepare(gamera_no_legacy_grid(),
                                 &driver_background) != 0) {
    gamera_no_driver_finalize();
    return -1;
  }
#endif
#ifdef GAMERA_NONORTHOGONAL_MATCHED_HPLUS_FEEDBACK
  if (gamera_no_im_runtime_prepare() != 0) {
    gamera_no_driver_finalize();
    return -1;
  }
#endif
#ifdef GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC
  if (gamera_no_im_live_driver_prepare() != 0) {
    log_error("Unable to prepare the default-OFF diagnostic-only live H+ "
              "adapter; MHD feedback remains disabled");
    gamera_no_driver_finalize();
    return -1;
  }
#endif
  return 0;
}

void gamera_no_driver_finalize(void) {
#ifdef GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC
  /* The runtime borrows the producer callback/context.  Release that binding
   * before destroying the producer that backs it. */
  gamera_no_im_online_runtime_finalize();
  gamera_no_im_live_driver_finalize();
#endif
#ifdef GAMERA_NONORTHOGONAL_MATCHED_HPLUS_FEEDBACK
  gamera_no_im_runtime_finalize();
#endif
#ifdef GAMERA_MI_COUPLING
  gamera_mi_coupling_finalize();
#endif
#ifdef GAMERA_NONORTHOGONAL_HAS_BACKGROUND_FIELD
  gamera_no_background_destroy(&driver_background);
  driver_background_ready = 0;
#endif
}

#if defined(GAMERA_YINYANG_MHD_BACKEND) || \
    defined(GAMERA_NONORTHOGONAL_HAS_EDGE_EMF_BOUNDARY)
static int synchronize_edge_emf(
    gamera_no_storage *storage, const gamera_no_grid *grid,
    const size_t lower[3], const size_t upper[3], void *context) {
  (void)context;
#if defined(GAMERA_YINYANG_MHD_BACKEND) && \
    !defined(GAMERA_YINYANG_MFE_INTERFACE)
  if (gamera_no_yinyang_sync_edge_emf(storage, grid, lower, upper, NULL) !=
      0) {
    return -1;
  }
#endif
#ifdef GAMERA_NONORTHOGONAL_HAS_EDGE_EMF_BOUNDARY
  if (problem_nonorthogonal_edge_emf_boundary(storage, grid, lower, upper) !=
      0) {
    return -1;
  }
#endif
#ifdef GAMERA_MI_COUPLING
  if (gamera_mi_coupling_apply_held_emf(storage, grid, lower, upper) != 0) {
    return -1;
  }
#endif
  return 0;
}
#endif

static void active_upper_values(size_t upper[3]) {
  upper[0] = (size_t)ie + 1U;
  upper[1] = (size_t)je + 1U;
  upper[2] = (size_t)ke + 1U;
}

static int exchange_current_state(int include_active_yinyang_receptors) {
  (void)include_active_yinyang_receptors;
  if (gamera_no_legacy_export() != 0 ||
      boundary_exchange_4d(gem, gem_onface_i, gem_onface_j,
                           gem_onface_k) != 0 ||
      boundary_exchange_5d(gas, NS1, config.NI, config.NJ, config.NK) != 0) {
    return -1;
  }
  boundary_conditions();
#ifdef GAMERA_YINYANG_BACKEND
  const int yinyang_status = include_active_yinyang_receptors
                                 ? gamera_no_yinyang_exchange_hd()
                                 : gamera_no_yinyang_exchange_hd_ghosts_only();
  if (yinyang_status != 0) {
    return -1;
  }
  /* Refresh radial/corner physical halos after angular donor interpolation. */
  boundary_conditions();
#endif
  if (gamera_no_legacy_import_current() != 0) {
    return -1;
  }
#ifdef GAMERA_NONORTHOGONAL_HAS_MAGNETIC_BOUNDARY
  if (problem_nonorthogonal_magnetic_boundary(
          gamera_no_legacy_storage(), gamera_no_legacy_grid()) != 0) {
    return -1;
  }
#endif
#ifdef GAMERA_YINYANG_MHD_BACKEND
  if (gamera_no_yinyang_exchange_magnetic_ghosts() != 0) {
    return -1;
  }
#ifdef GAMERA_NONORTHOGONAL_HAS_MAGNETIC_BOUNDARY
  /* Angular overset interpolation may touch radial corner ghosts. */
  if (problem_nonorthogonal_magnetic_boundary(
          gamera_no_legacy_storage(), gamera_no_legacy_grid()) != 0) {
    return -1;
  }
#endif
#endif
  if (gamera_no_legacy_export() != 0) {
    return -1;
  }
  const double local_divergence = divB_max;
  if (MPI_Allreduce(&local_divergence, &divB_max, 1, MPI_DOUBLE, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS) {
    return -1;
  }
#ifdef GAMERA_NONORTHOGONAL_HAS_BACKGROUND_FIELD
  if (update_total_divergence_diagnostic() != 0) {
    return -1;
  }
#endif
  return 0;
}

static gamera_no_advance_options advance_options(double source_time) {
  gamera_no_sweep_options stress = {
      gamma_val, rho_floor, p_floor, PDMB, 0.25, 0.25, CA,
      true,      true,      false,   false};
  gamera_no_emf_options emf = {rho_floor, PDMB, 0.5, CA, CFL,
                                0.0,       false, false};
  gamera_no_update_options update = {
      gamma_val, rho_floor, p_floor, CA, true, false, false};
#ifdef GAMERA_NONORTHOGONAL_USE_BORIS
  stress.use_boris = true;
  emf.use_boris = true;
  update.use_boris = true;
#endif
#ifdef GAMERA_NONORTHOGONAL_HAS_BACKGROUND_FIELD
  stress.use_background = true;
  emf.use_background = true;
  update.use_background = true;
#endif
  gamera_no_advance_options result = {
      .stress = stress,
      .emf = emf,
      .update = update,
      .cell_source = NULL,
      .source_context = NULL,
      .indexed_cell_source = NULL,
      .indexed_source_context = NULL,
      .source_time = source_time,
      .fluid_flux_sync = NULL,
      .fluid_flux_context = NULL,
      .edge_emf_sync = NULL,
      .edge_emf_context = NULL};
#ifdef GAMERA_NONORTHOGONAL_HAS_CELL_SOURCE
  result.cell_source = problem_nonorthogonal_cell_source;
#endif
#ifdef GAMERA_NONORTHOGONAL_MATCHED_HPLUS_FEEDBACK
  if (gamera_no_im_runtime_source_active()) {
    result.indexed_cell_source = gamera_no_im_pressure_source_callback;
    result.indexed_source_context = gamera_no_im_runtime_source_cache();
  }
#endif
#ifdef GAMERA_NONORTHOGONAL_ONLINE_HPLUS_PRESSURE_FEEDBACK
  if (gamera_no_im_live_driver_source_active()) {
    result.indexed_cell_source =
        gamera_no_im_live_driver_source_callback;
    result.indexed_source_context =
        gamera_no_im_live_driver_source_context();
  }
#endif
#ifdef GAMERA_NONORTHOGONAL_HAS_FLUID_FLUX_BOUNDARY
  result.fluid_flux_sync = synchronize_fluid_flux;
#endif
#ifdef GAMERA_YINYANG_MHD_BACKEND
  result.edge_emf_sync = synchronize_edge_emf;
#elif defined(GAMERA_NONORTHOGONAL_HAS_EDGE_EMF_BOUNDARY)
  result.edge_emf_sync = synchronize_edge_emf;
#endif
  return result;
}

static int global_timestep(double *result) {
  gamera_no_grid *grid = gamera_no_legacy_grid();
  gamera_no_storage *storage = gamera_no_legacy_storage();
  if (grid == NULL || storage == NULL || result == NULL) {
    return -1;
  }
  size_t active_upper[3];
  active_upper_values(active_upper);
  gamera_no_timestep_options options = {
      gamma_val, rho_floor, p_floor, CFL, CA, true, false, false};
#ifdef GAMERA_NONORTHOGONAL_USE_BORIS
  options.use_boris = true;
#endif
#ifdef GAMERA_NONORTHOGONAL_HAS_BACKGROUND_FIELD
  options.use_background = true;
#endif
  double local;
  const gamera_no_vec3 *background_cell = NULL;
#ifdef GAMERA_NONORTHOGONAL_HAS_BACKGROUND_FIELD
  if (!driver_background_ready) {
    return -1;
  }
  background_cell = driver_background.cell_magnetic;
#endif
  if (gamera_no_local_timestep(grid, storage->conserved,
                               storage->cell_magnetic, active_lower,
                               active_upper, &options, background_cell,
                               &local) != 0) {
    return -1;
  }
  return MPI_Allreduce(&local, result, 1, MPI_DOUBLE, MPI_MIN,
                       MPI_COMM_WORLD) ==
                 MPI_SUCCESS
             ? 0
             : -1;
}

static int timed_analysis_output(void) {
  const double begin = MPI_Wtime();
  const int local_failed = dump_analysis_hdf5() != 0;
  const double local_seconds = MPI_Wtime() - begin;
  int global_failed = 0;
  double maximum_seconds = 0.0;
  if (MPI_Allreduce(&local_failed, &global_failed, 1, MPI_INT, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Allreduce(&local_seconds, &maximum_seconds, 1, MPI_DOUBLE, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS) {
    return -1;
  }
  if (rank == 0) {
    log_info("Compact analysis I/O max-rank wall time %.9g s",
             maximum_seconds);
  }
  return global_failed ? -1 : 0;
}

#ifdef GAMERA_MI_COUPLING
static int append_regular_mi_restart(void) {
  gamera_mi_restart_state *state = NULL;
  int local_failed = gamera_mi_coupling_export_restart(&state) != 0;
  int global_failed = 0;
  if (MPI_Allreduce(&local_failed, &global_failed, 1, MPI_INT, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS ||
      global_failed) {
    gamera_mi_restart_destroy(state);
    return -1;
  }
  if (rank == 0) {
    char path[256];
    int count;
    if (patch_count > 1) {
      count = snprintf(path, sizeof(path),
                       "%s_p%d_%02d-%02d-%02d_%06d.h5",
                       BASE_RESTART_NAME, patch_id, proc_coords[0],
                       proc_coords[1], proc_coords[2], hdf_seq_num);
    } else {
      count = snprintf(path, sizeof(path), out_file_pattern,
                       BASE_RESTART_NAME, proc_coords[0], proc_coords[1],
                       proc_coords[2], hdf_seq_num, "h5");
    }
    local_failed = count < 0 || (size_t)count >= sizeof(path) ||
                   gamera_mi_restart_append_hdf5(path, state) != 0;
  }
  gamera_mi_restart_destroy(state);
  if (MPI_Allreduce(&local_failed, &global_failed, 1, MPI_INT, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS ||
      global_failed) {
    return -1;
  }
  return 0;
}

static int restore_regular_mi_restart(void) {
#ifdef GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC
  if (im_online_hplus_config.enabled) {
    return 0;
  }
#endif
  if (!read_restart || !mi_config.enabled) {
    return 0;
  }
  gamera_mi_restart_state *state = NULL;
  int present = 0;
  int local_failed = 0;
  if (rank == 0) {
    present = restart_source_filename[0] == '\0'
                  ? -1
                  : gamera_mi_restart_hdf5_present(
                        restart_source_filename);
    local_failed = present < 0;
    if (!local_failed && present > 0) {
      local_failed = gamera_mi_restart_read_hdf5(
                         restart_source_filename, &state) != 0;
    }
  }
  int global_failed = 0;
  if (MPI_Bcast(&present, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Allreduce(&local_failed, &global_failed, 1, MPI_INT, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS ||
      global_failed) {
    gamera_mi_restart_destroy(state);
    return -1;
  }
  if (present == 0) {
    if (mi_config.electron_precipitation_enabled) {
      if (rank == 0) {
        log_error("Precipitation-enabled restart is missing the required "
                  "electrostatic M-I state in %s",
                  restart_source_filename);
      }
      return -1;
    }
    if (rank == 0) {
      log_warn("Legacy restart has no electrostatic M-I state; performing "
               "a fresh constant-conductance solve at restart time");
    }
    return 0;
  }
  if (gamera_mi_restart_broadcast(MPI_COMM_WORLD, 0, &state) != 0) {
    gamera_mi_restart_destroy(state);
    return -1;
  }
  local_failed = gamera_mi_coupling_restore_restart(state, time_sim) != 0;
  gamera_mi_restart_destroy(state);
  if (MPI_Allreduce(&local_failed, &global_failed, 1, MPI_INT, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS ||
      global_failed) {
    return -1;
  }
  if (rank == 0) {
    log_info("Restored electrostatic M-I potential, conductance, and hybrid "
             "DPB state from %s", restart_source_filename);
  }
  return 0;
}
#endif

static int timed_restart_output(void) {
  const double begin = MPI_Wtime();
  char original_directory[PATH_MAX];
  char checkpoint_directory[PATH_MAX];
  const int checkpoint_sequence = hdf_seq_num + 1;
  int local_failed = 0;
  int global_failed = 0;
  int use_online_hplus_manifest = 0;
#ifdef GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC
  use_online_hplus_manifest = im_online_hplus_config.enabled != 0;
#endif

  if (getcwd(original_directory, sizeof(original_directory)) == NULL ||
      snprintf(checkpoint_directory, sizeof(checkpoint_directory),
               "restart/checkpoint_%06d", checkpoint_sequence) < 0) {
    local_failed = 1;
  }
  if (rank == 0 && !local_failed) {
    if ((mkdir("restart", 0750) != 0 && errno != EEXIST) ||
        mkdir(checkpoint_directory, 0750) != 0) {
      log_error("Unable to create restart checkpoint directory %s: %s",
                checkpoint_directory, strerror(errno));
      local_failed = 1;
    }
  }
  if (MPI_Allreduce(&local_failed, &global_failed, 1, MPI_INT, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS || global_failed ||
      MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
    return -1;
  }
  if (chdir(checkpoint_directory) != 0) {
    log_error("Unable to enter restart checkpoint directory %s: %s",
              checkpoint_directory, strerror(errno));
    local_failed = 1;
  } else if (dump_data_hdf5(BASE_RESTART_NAME) != 0 ||
             dump_extra_data(BASE_RESTART_NAME) != 0
#ifdef GAMERA_NONORTHOGONAL_MATCHED_HPLUS_FEEDBACK
             || gamera_no_im_runtime_write_restart() != 0
#endif
  ) {
    local_failed = 1;
  }
  if (!local_failed && hdf_seq_num != checkpoint_sequence) {
    log_error("Restart rank sequence mismatch: wrote %d, expected %d",
              hdf_seq_num, checkpoint_sequence);
    local_failed = 1;
  }
  if (MPI_Allreduce(&local_failed, &global_failed, 1, MPI_INT, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS) {
    return -1;
  }
#ifdef GAMERA_MI_COUPLING
  if (!global_failed && !use_online_hplus_manifest &&
      append_regular_mi_restart() != 0) {
    local_failed = 1;
    global_failed = 1;
  }
#endif
  if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
    return -1;
  }
  if (!global_failed && use_online_hplus_manifest) {
    int world_size = 0;
    const int manifest_setup_failed =
        MPI_Comm_size(MPI_COMM_WORLD, &world_size) != MPI_SUCCESS ||
        world_size <= 0;
    int manifest_setup_failed_global = 0;
    if (MPI_Allreduce(&manifest_setup_failed, &manifest_setup_failed_global,
                      1, MPI_INT, MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS) {
      return -1;
    }
    local_failed = manifest_setup_failed_global;
#ifdef GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC
    if (!local_failed) {
      local_failed = gamera_no_im_live_driver_checkpoint_publish(
                         ".", (uint64_t)checkpoint_sequence, time_sim,
                         (uint64_t)world_size) != 0;
    }
#endif
  }
  if (rank == 0 && !global_failed && !use_online_hplus_manifest) {
    int world_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    FILE *manifest = fopen("manifest.json.tmp", "wx");
    if (manifest == NULL) {
      local_failed = 1;
    } else {
      int manifest_failed =
          fprintf(manifest,
                "{\n  \"schema_version\": 1,\n"
                "  \"complete\": true,\n"
                "  \"checkpoint_sequence\": %d,\n"
                "  \"time_code\": %.17g,\n"
                "  \"time_seconds\": %.17g,\n"
                "  \"rank_file_count\": %d\n}\n",
                checkpoint_sequence, time_sim,
                time_sim * norm_config.Time_Norm, world_size) < 0 ||
          fflush(manifest) != 0 || fsync(fileno(manifest)) != 0;
      if (fclose(manifest) != 0) {
        manifest_failed = 1;
      }
      manifest = NULL;
      if (manifest_failed ||
          rename("manifest.json.tmp", "manifest.json") != 0) {
        local_failed = 1;
      }
    }
    if (local_failed) {
      log_error("Unable to publish restart checkpoint manifest");
      if (manifest != NULL) {
        fclose(manifest);
      }
    }
  }
  if (chdir(original_directory) != 0) {
    log_error("Unable to return from restart checkpoint directory: %s",
              strerror(errno));
    local_failed = 1;
  }
  const double local_seconds = MPI_Wtime() - begin;
  double maximum_seconds = 0.0;
  if (MPI_Allreduce(&local_failed, &global_failed, 1, MPI_INT, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Allreduce(&local_seconds, &maximum_seconds, 1, MPI_DOUBLE, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS) {
    return -1;
  }
  if (rank == 0) {
    log_info("Full restart I/O max-rank wall time %.9g s", maximum_seconds);
  }
  return global_failed ? -1 : 0;
}

int solve_nonorthogonal(int maximum_steps) {
  gamera_no_grid *grid = gamera_no_legacy_grid();
  gamera_no_storage *storage = gamera_no_legacy_storage();
  const int initial_exchange_includes_active_receptors = !read_restart;
#ifdef GAMERA_BENCH_PROFILE_JSONL
  int bench_warmup_steps = 0;
  int bench_timed_steps = 0;
  const char *bench_run_id = getenv("GAMERA_BENCH_RUN_ID");
  const char *bench_manifest_sha256 =
      getenv("GAMERA_BENCH_MANIFEST_SHA256");
  const char *bench_plan_sha256 = getenv("GAMERA_BENCH_PLAN_SHA256");
  const char *bench_binary_sha256 = getenv("GAMERA_BENCH_BINARY_SHA256");
  const char *bench_build_manifest_sha256 =
      getenv("GAMERA_BENCH_BUILD_MANIFEST_SHA256");
  const char *bench_fault_kind = getenv("GAMERA_BENCH_FAULT_KIND");
  int bench_fault_rank = -1;
  int bench_expected_omp_threads = 0;
  const int bench_omp_max_threads = omp_get_max_threads();
  const int bench_omp_dynamic = omp_get_dynamic();
  int bench_omp_team_threads = 0;
#pragma omp parallel
  {
#pragma omp single
    { bench_omp_team_threads = omp_get_num_threads(); }
  }
  int bench_contract_error =
      gamera_bench_env_steps("GAMERA_BENCH_WARMUP_STEPS",
                             &bench_warmup_steps) != 0 ||
      gamera_bench_env_steps("GAMERA_BENCH_TIMED_STEPS",
                             &bench_timed_steps) != 0 ||
      gamera_bench_env_steps("OMP_NUM_THREADS",
                             &bench_expected_omp_threads) != 0 ||
      bench_timed_steps <= 0 ||
      bench_expected_omp_threads <= 0 || bench_omp_dynamic != 0 ||
      bench_omp_max_threads != bench_expected_omp_threads ||
      bench_omp_team_threads != bench_expected_omp_threads ||
      bench_warmup_steps > INT_MAX - bench_timed_steps ||
      maximum_steps != bench_warmup_steps + bench_timed_steps ||
      maximum_steps <= 0 || maximum_steps > GAMERA_BENCH_MAX_STEPS ||
      !gamera_bench_safe_text(bench_run_id, 0) ||
      !gamera_bench_safe_text(bench_manifest_sha256, 1) ||
      !gamera_bench_safe_text(bench_plan_sha256, 1) ||
      !gamera_bench_safe_text(bench_binary_sha256, 1) ||
      !gamera_bench_safe_text(bench_build_manifest_sha256, 1);
  if (bench_fault_kind != NULL) {
    bench_contract_error |=
        gamera_bench_env_steps("GAMERA_BENCH_FAULT_RANK",
                               &bench_fault_rank) != 0 ||
        (strcmp(bench_fault_kind, "env") != 0 &&
         strcmp(bench_fault_kind, "timer") != 0 &&
         strcmp(bench_fault_kind, "write") != 0);
  }
  if (bench_fault_rank == rank && bench_fault_kind != NULL &&
      strcmp(bench_fault_kind, "env") == 0) {
    bench_contract_error = 1;
  }
  if (gamera_bench_consensus(bench_contract_error)) {
    log_error("Invalid fixed-Nt benchmark contract: Nt=%d", maximum_steps);
    return -1;
  }
  int bench_sample_count = 0;
  double bench_initial_mi_s = 0.0;
  int bench_timer_error = 0;
  gamera_bench_step_sample bench_samples[GAMERA_BENCH_MAX_STEPS] = {{0}};
#endif
  if (grid == NULL || storage == NULL || maximum_steps < 0 ||
      exchange_current_state(initial_exchange_includes_active_receptors) !=
          0) {
    return -1;
  }
  if (rank == 0) {
    if (analysis_output_enabled) {
      log_info("Analysis I/O standard: one parallel HDF5 file per patch "
               "(analysis_p0/analysis_p1); legacy per-rank analysis is "
               "retired");
    } else {
      log_info("Routine analysis output is disabled; legacy per-rank "
               "analysis remains retired");
    }
  }
#if defined(GAMERA_MI_COUPLING) ||                                         \
    defined(GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC)
  int mi_updated = 0;
#endif
#ifdef GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC
  /* A strict restart restores M-I from the composite sidecar before the
   * initial maybe_update below.  This preserves the held potential/EMF and
   * prevents a duplicate restart-time solve or diagnostic publication. */
  if (read_restart &&
      gamera_no_im_live_driver_initialize(time_sim, 0,
                                          &driver_background) < 0) {
    log_error("Initial live-H+ restart restore failed before electrostatic "
              "M-I scheduling; no new pressure source was accepted");
    return -1;
  }
#endif
#ifdef GAMERA_MI_COUPLING
#ifdef GAMERA_BENCH_PROFILE_JSONL
  const double bench_initial_mi_start = omp_get_wtime();
#endif
  if (restore_regular_mi_restart() != 0) {
    log_error("Initial electrostatic M-I restart restore failed");
    return -1;
  }
  mi_updated = gamera_mi_coupling_maybe_update(time_sim);
  if (mi_updated < 0 ||
      (mi_updated > 0 &&
       exchange_current_state(0) != 0)) {
    log_error("Initial electrostatic M-I coupling update failed");
    return -1;
  }
#ifdef GAMERA_BENCH_PROFILE_JSONL
  bench_initial_mi_s = omp_get_wtime() - bench_initial_mi_start;
  if (!isfinite(bench_initial_mi_s) || bench_initial_mi_s < 0.0) {
    bench_initial_mi_s = 0.0;
    bench_timer_error = 1;
  }
#endif
#endif

#ifdef GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC
  /* The live diagnostic gate deliberately performs one t=0 collection after
   * the synchronized MHD halo and a fresh M-I solve.  That exact shell defines
   * the H+ reference volume and initial Maxwellian.  Thereafter the adapter
   * gathers only at canonical runtime due epochs. */
  if (!read_restart &&
      gamera_no_im_live_driver_initialize(time_sim, mi_updated,
                                          &driver_background) < 0) {
    log_error("Initial live-H+ t=0 bootstrap failed; no pressure source was "
              "published");
    return -1;
  }
#endif

  if (!read_restart && global_timestep(&dt0) != 0) {
    log_error("Unable to initialize the non-orthogonal AB2 history timestep");
    return -1;
  }
  if (!isfinite(dt0) || dt0 <= 0.0) {
    log_error("Invalid non-orthogonal AB2 history timestep dt0=%.17g%s", dt0,
              read_restart ? " in restart checkpoint" : "");
    return -1;
  }

  double last_restart_time = time_sim;
  double last_analysis_time = time_sim;
  double next_analysis_time = INFINITY;
  double next_restart_time = INFINITY;
  if (analysis_output_enabled) {
    if (initialize_analysis_sequence() != 0) {
      log_error("Unable to recover compact-analysis output sequence");
      return -1;
    }
    next_analysis_time =
        output_interval *
        (floor(time_sim / output_interval + 1.0e-12) + 1.0);
  }
  if (restart_interval > 0.0) {
    next_restart_time =
        restart_interval *
        (floor(time_sim / restart_interval + 1.0e-12) + 1.0);
  }
  if (!read_restart) {
    hdf_seq_num = -1;
    log_seq_num = 0;
    if (timed_restart_output() != 0 ||
        (analysis_output_enabled && timed_analysis_output() != 0)) {
      return -1;
    }
    last_restart_time = 0.0;
    if (analysis_output_enabled) {
      last_analysis_time = 0.0;
    }
  }

  time_sim_start = time_sim;
  size_t active_upper[3];
  active_upper_values(active_upper);
#ifdef GAMERA_NONORTHOGONAL_HAS_FLUID_FLUX_BOUNDARY
  unsigned long long wall_clamped_window[2] = {0ULL, 0ULL};
  double wall_max_window[26] = {0.0};
  gamera_no_vec3 wall_db_location_window[2] = {
      {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}};
#ifdef GAMERA_MI_COUPLING
  double next_uncoupled_wall_diagnostic_s =
      (floor(time_sim * norm_config.Time_Norm / 10.0) + 1.0) * 10.0;
#endif
#endif
#ifdef GAMERA_YINYANG_HDIV_PROFILE
  unsigned long profile_step_calls = 0;
  double profile_step_sum = 0.0;
  double profile_step_maximum = 0.0;
#endif
  for (int step = 1; step <= maximum_steps; ++step) {
    if (time_sim >= time_stop) {
      break;
    }
#ifdef GAMERA_YINYANG_HDIV_PROFILE
    const double profile_step_start = omp_get_wtime();
#endif
#ifdef GAMERA_BENCH_PROFILE_JSONL
    const double bench_step_start = omp_get_wtime();
    double bench_mi_s = 0.0;
#endif
    if (global_timestep(&dt) != 0 || !isfinite(dt) || dt <= 0.0) {
      return -1;
    }
#ifdef GAMERA_MI_COUPLING
    dt = gamera_mi_coupling_limit_timestep(time_sim, dt);
#endif
#ifdef GAMERA_NONORTHOGONAL_MATCHED_HPLUS_FEEDBACK
    double feedback_limited_dt = dt;
    if (gamera_no_im_runtime_limit_timestep(time_sim, dt,
                                            &feedback_limited_dt) != 0) {
      log_error("Unable to land timestep on matched-H+ coupling cadence");
      return -1;
    }
    dt = feedback_limited_dt;
#endif
#ifdef GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC
    if (gamera_no_im_online_runtime_enabled()) {
      double online_limited_dt = dt;
      if (gamera_no_im_online_runtime_limit_timestep(
              time_sim, dt, &online_limited_dt) != 0) {
        log_error("Unable to land timestep on diagnostic online-H+ cadence");
        return -1;
      }
      dt = online_limited_dt;
    }
#endif
    if (time_sim + dt > time_stop) {
      dt = time_stop - time_sim;
    }
    if (!(dt > 0.0)) {
      break;
    }
    const double predictor_ratio = 0.5 * dt / dt0;
    const gamera_no_advance_options options =
        advance_options(time_sim + 0.5 * dt);
    const gamera_no_background_field *background = NULL;
#ifdef GAMERA_NONORTHOGONAL_HAS_BACKGROUND_FIELD
    background = gamera_no_driver_background_field();
#endif
    if (gamera_no_advance(storage, grid, active_lower, active_upper,
                          predictor_ratio, dt, &options, background) != 0) {
      log_error("Non-orthogonal advance failed at step %d", step);
      return -1;
    }
#ifdef GAMERA_EARTH_DIPOLE_BACKGROUND
    size_t chilled_low_density = 0U;
    size_t chilled_high_sound = 0U;
    double chillout_maximum_sound = 0.0;
    double chillout_maximum_pressure = 0.0;
    if (gamera_no_apply_kaiju_chillout(
            storage, grid, active_lower, active_upper, dt, gamma_val,
            rho_floor, p_floor, CA, 1.0e-3,
            1.0e-2 / norm_config.u_Norm, &chilled_low_density,
            &chilled_high_sound, &chillout_maximum_sound,
            &chillout_maximum_pressure) != 0) {
      log_error("Kaiju ChillOut failed at step %d", step);
      return -1;
    }
    if (chilled_low_density > 0U || chilled_high_sound > 0U) {
      log_info("Kaiju ChillOut rank %d: low-density=%zu high-sound=%zu "
               "pre-cooling max cs=%.6g pressure=%.6g code",
               rank, chilled_low_density, chilled_high_sound,
               chillout_maximum_sound, chillout_maximum_pressure);
    }
#endif
    if (storage->nuclear_hogs_face_count > 0U) {
      log_info("Kaiju nuclear HOGS rank %d: faces=%zu max interface "
               "speed=%.6g code (threshold=%.6g)",
               rank, storage->nuclear_hogs_face_count,
               storage->nuclear_hogs_max_interface_speed, 1.5 * CA);
    }
#ifdef GAMERA_NONORTHOGONAL_HAS_FLUID_FLUX_BOUNDARY
    for (int hemisphere = 0; hemisphere < 2; ++hemisphere) {
      wall_clamped_window[hemisphere] +=
          (unsigned long long)
              storage->inner_wall_clamped_face_count[hemisphere];
      wall_max_window[hemisphere] =
          fmax(wall_max_window[hemisphere],
               storage->inner_wall_positive_mass_max[hemisphere]);
      wall_max_window[2 + hemisphere] =
          fmax(wall_max_window[2 + hemisphere],
               storage->inner_wall_positive_energy_max[hemisphere]);
      wall_max_window[4 + hemisphere] =
          fmax(wall_max_window[4 + hemisphere],
               storage->inner_wall_fluid_momentum_max[hemisphere]);
      wall_max_window[6 + hemisphere] =
          fmax(wall_max_window[6 + hemisphere],
               storage->inner_wall_maxwell_momentum_max[hemisphere]);
      wall_max_window[8 + hemisphere] =
          fmax(wall_max_window[8 + hemisphere],
               storage->inner_wall_density_max[hemisphere]);
      wall_max_window[10 + hemisphere] =
          fmax(wall_max_window[10 + hemisphere],
               storage->inner_wall_pressure_max[hemisphere]);
      wall_max_window[12 + hemisphere] =
          fmax(wall_max_window[12 + hemisphere],
               storage->inner_wall_pressure_gradient_max[hemisphere]);
      wall_max_window[14 + hemisphere] =
          fmax(wall_max_window[14 + hemisphere],
               storage->inner_wall_speed_max[hemisphere]);
      if (storage->inner_wall_residual_magnetic_max[hemisphere] >
          wall_max_window[16 + hemisphere]) {
        wall_max_window[16 + hemisphere] =
            storage->inner_wall_residual_magnetic_max[hemisphere];
        wall_db_location_window[hemisphere] =
            storage->inner_wall_residual_magnetic_location[hemisphere];
      }
      wall_max_window[18 + hemisphere] =
          fmax(wall_max_window[18 + hemisphere],
               storage->inner_wall_lowlat_pressure_max[hemisphere]);
      wall_max_window[20 + hemisphere] = fmax(
          wall_max_window[20 + hemisphere],
          storage->inner_wall_lowlat_pressure_gradient_max[hemisphere]);
      wall_max_window[22 + hemisphere] =
          fmax(wall_max_window[22 + hemisphere],
               storage->inner_wall_lowlat_speed_max[hemisphere]);
      wall_max_window[24 + hemisphere] = fmax(
          wall_max_window[24 + hemisphere],
          storage->inner_wall_lowlat_residual_magnetic_max[hemisphere]);
    }
#endif
    dt0 = dt;
    time_sim += dt;
#ifdef GAMERA_NONORTHOGONAL_MATCHED_HPLUS_FEEDBACK
    if (gamera_no_im_runtime_enabled()) {
      double canonical_time = time_sim;
      if (gamera_no_im_runtime_canonicalize_step_end(
              time_sim, time_stop, &canonical_time) != 0) {
        log_error("Unable to canonicalize matched-H+ cadence/stop landing");
        return -1;
      }
      time_sim = canonical_time;
    }
#endif
#ifdef GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC
    if (gamera_no_im_online_runtime_enabled()) {
      double canonical_time = time_sim;
      if (gamera_no_im_online_runtime_canonicalize_step_end(
              time_sim, time_stop, &canonical_time) != 0) {
        log_error("Unable to canonicalize diagnostic online-H+ "
                  "cadence/stop landing");
        return -1;
      }
      time_sim = canonical_time;
    }
#endif

    if (exchange_current_state(1) != 0) {
      log_error("Non-orthogonal halo exchange failed at step %d", step);
      return -1;
    }
#ifdef GAMERA_MI_COUPLING
#ifdef GAMERA_BENCH_PROFILE_JSONL
    const double bench_mi_start = omp_get_wtime();
#endif
    mi_updated = gamera_mi_coupling_maybe_update(time_sim);
    if (mi_updated < 0 ||
        (mi_updated > 0 && exchange_current_state(0) != 0)) {
      log_error("Electrostatic M-I coupling update failed at step %d", step);
      return -1;
    }
#ifdef GAMERA_NONORTHOGONAL_HAS_FLUID_FLUX_BOUNDARY
    const int uncoupled_wall_diagnostic_due =
        !mi_config.enabled &&
        time_sim * norm_config.Time_Norm + 1.0e-10 >=
            next_uncoupled_wall_diagnostic_s;
    if (mi_updated > 0 || uncoupled_wall_diagnostic_due) {
      unsigned long long global_count[2] = {0ULL, 0ULL};
      double global_maximum[26] = {0.0};
      if (MPI_Allreduce(wall_clamped_window, global_count, 2,
                        MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                        MPI_COMM_WORLD) != MPI_SUCCESS ||
          MPI_Allreduce(wall_max_window, global_maximum, 26, MPI_DOUBLE,
                        MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS) {
        return -1;
      }
      if (rank == 0) {
        log_info("Inner-wall flux window at t=%.6f s N/S: clamped="
                 "%llu/%llu max positive mass=%.6g/%.6g energy=%.6g/%.6g "
                 "max |fluid momentum|=%.6g/%.6g |Maxwell momentum|="
                 "%.6g/%.6g",
                 time_sim * norm_config.Time_Norm, global_count[0],
                 global_count[1], global_maximum[0],
                 global_maximum[1], global_maximum[2], global_maximum[3],
                 global_maximum[4], global_maximum[5], global_maximum[6],
                 global_maximum[7]);
        log_info("Inner active shell at t=%.6f s N/S: max density="
                 "%.6g/%.6g pressure="
                 "%.6g/%.6g |dP/dr|=%.6g/%.6g |v|=%.6g/%.6g "
                 "|dB|=%.6g/%.6g",
                 time_sim * norm_config.Time_Norm, global_maximum[8],
                 global_maximum[9], global_maximum[10],
                 global_maximum[11], global_maximum[12],
                 global_maximum[13], global_maximum[14],
                 global_maximum[15], global_maximum[16],
                 global_maximum[17]);
        log_info("Inner active shell low-|MLAT|<=15 deg at t=%.6f s N/S: "
                 "max pressure="
                 "%.6g/%.6g |dP/dr|=%.6g/%.6g |v|=%.6g/%.6g "
                 "|dB|=%.6g/%.6g",
                 time_sim * norm_config.Time_Norm, global_maximum[18],
                 global_maximum[19],
                 global_maximum[20], global_maximum[21],
                 global_maximum[22], global_maximum[23],
                 global_maximum[24], global_maximum[25]);
      }
      for (int hemisphere = 0; hemisphere < 2; ++hemisphere) {
        if (wall_max_window[16 + hemisphere] > 0.0) {
          const gamera_no_vec3 point = wall_db_location_window[hemisphere];
          const double radius = sqrt(point.value[0] * point.value[0] +
                                     point.value[1] * point.value[1] +
                                     point.value[2] * point.value[2]);
          const double latitude =
              asin(point.value[2] / radius) * 180.0 / PI;
          double longitude = atan2(point.value[1], point.value[0]);
          if (longitude < 0.0) {
            longitude += 2.0 * PI;
          }
          log_info("Inner-shell local |dB| maximum at t=%.6f s rank=%d "
                   "patch=%d hemisphere=%s value=%.6g xyz=(%.6g,%.6g,%.6g) "
                   "MLAT=%.6g MLT=%.6g",
                   time_sim * norm_config.Time_Norm, rank, patch_id,
                   hemisphere == 0 ? "N" : "S",
                   wall_max_window[16 + hemisphere], point.value[0],
                   point.value[1], point.value[2], latitude,
                   /* The ionosphere/polar-plot convention is phi=0 at
                    * midnight (00 MLT), increasing toward 06 MLT. */
                   fmod(longitude * 12.0 / PI, 24.0));
        }
      }
      for (int hemisphere = 0; hemisphere < 2; ++hemisphere) {
        wall_clamped_window[hemisphere] = 0ULL;
      }
      for (int diagnostic = 0; diagnostic < 26; ++diagnostic) {
        wall_max_window[diagnostic] = 0.0;
      }
      for (int hemisphere = 0; hemisphere < 2; ++hemisphere) {
        wall_db_location_window[hemisphere] =
            (gamera_no_vec3){{0.0, 0.0, 0.0}};
      }
      if (uncoupled_wall_diagnostic_due) {
        do {
          next_uncoupled_wall_diagnostic_s += 10.0;
        } while (time_sim * norm_config.Time_Norm + 1.0e-10 >=
                 next_uncoupled_wall_diagnostic_s);
      }
    }
#endif
#ifdef GAMERA_BENCH_PROFILE_JSONL
    bench_mi_s = omp_get_wtime() - bench_mi_start;
#endif
#endif

#ifdef GAMERA_NONORTHOGONAL_MATCHED_HPLUS_FEEDBACK
    if (gamera_no_im_runtime_maybe_refresh(time_sim) < 0) {
      log_error("Matched-H+ prescribed target/source refresh failed at step %d",
                step);
      return -1;
    }
#endif
#ifdef GAMERA_NONORTHOGONAL_ONLINE_HPLUS_DIAGNOSTIC
    if (gamera_no_im_live_driver_refresh(time_sim, mi_updated,
                                         &driver_background) < 0) {
      log_error("Live-H+ due-shell/feedback refresh failed at step %d; "
                "stopping before another source interval",
                step);
      return -1;
    }
#endif

#ifdef GAMERA_BENCH_PROFILE_JSONL
    double bench_step_excluding_checkpoint_io =
        omp_get_wtime() - bench_step_start;
    const int forced_timer_failure =
        bench_fault_rank == rank && bench_fault_kind != NULL &&
        strcmp(bench_fault_kind, "timer") == 0 && step == 1;
    const int invalid_timer =
        bench_sample_count >= maximum_steps ||
        !isfinite(bench_step_excluding_checkpoint_io) ||
        !isfinite(bench_mi_s) ||
        bench_step_excluding_checkpoint_io < 0.0 || bench_mi_s < 0.0 ||
        bench_mi_s > bench_step_excluding_checkpoint_io + 1.0e-9;
    bench_timer_error |= invalid_timer || forced_timer_failure;
    if (invalid_timer) {
      bench_step_excluding_checkpoint_io = 0.0;
      bench_mi_s = 0.0;
    }
    bench_samples[bench_sample_count++] = (gamera_bench_step_sample){
        bench_step_excluding_checkpoint_io,
        fmax(0.0, bench_step_excluding_checkpoint_io - bench_mi_s),
        bench_mi_s,
#ifdef GAMERA_MI_COUPLING
        mi_updated > 0
#else
        0
#endif
    };
#endif
#ifdef GAMERA_YINYANG_HDIV_PROFILE
    ++profile_step_calls;
    if (profile_step_calls > 1) {
      const double profile_step_sample = omp_get_wtime() - profile_step_start;
      profile_step_sum += profile_step_sample;
      profile_step_maximum =
          fmax(profile_step_maximum, profile_step_sample);
    }
#endif
    const double output_tolerance =
        2.0e-12 * fmax(1.0, fabs(time_sim));
    if (analysis_output_enabled) {
      if (time_sim + output_tolerance >= next_analysis_time) {
        if (timed_analysis_output() != 0) {
          return -1;
        }
        last_analysis_time = time_sim;
        do {
          next_analysis_time += output_interval;
        } while (time_sim + output_tolerance >= next_analysis_time);
      }
    }
    if (restart_interval > 0.0 &&
        time_sim + output_tolerance >= next_restart_time) {
      if (timed_restart_output() != 0) {
        return -1;
      }
      last_restart_time = time_sim;
      do {
        next_restart_time += restart_interval;
      } while (time_sim + output_tolerance >= next_restart_time);
    }
    if (time_sim >= 0.1 * output_interval * (double)log_seq_num) {
#ifdef GAMERA_NONORTHOGONAL_HAS_BACKGROUND_FIELD
      log_info("Non-orthogonal step %d rank %d time %.15f dt %.15f "
               "max_divB1 %.4e max_divB_total %.4e",
               step, rank, time_sim, dt, divB_max, total_divB_max);
#else
      log_info("Non-orthogonal step %d rank %d time %.15f dt %.15f "
               "max_divB %.4e",
               step, rank, time_sim, dt, divB_max);
#endif
      ++log_seq_num;
    }
  }

#ifdef GAMERA_YINYANG_HDIV_PROFILE
  if (profile_step_calls > 1) {
    log_info("Core-step profile summary rank=%d patch=%d samples=%lu "
             "warmup_excluded=1 io_excluded=1 mean=%.9g max=%.9g",
             rank, patch_id, profile_step_calls - 1,
             profile_step_sum / (double)(profile_step_calls - 1),
             profile_step_maximum);
  }
#endif

#ifdef GAMERA_BENCH_PROFILE_JSONL
  double bench_final_max_divB_total = 0.0;
  int bench_final_diagnostic_error = 0;
#ifdef GAMERA_NONORTHOGONAL_HAS_BACKGROUND_FIELD
  bench_final_diagnostic_error =
      !isfinite(total_divB_max) || total_divB_max < 0.0;
  bench_final_max_divB_total = total_divB_max;
#else
  bench_final_diagnostic_error = 1;
#endif
  bench_final_diagnostic_error |=
      bench_sample_count != maximum_steps || !isfinite(time_sim) ||
      !isfinite(time_sim * norm_config.Time_Norm);
  log_info("Formal benchmark final diagnostic rank=%d patch=%d "
           "completed_steps=%d time_code=%.17g time_s=%.17g "
           "final_max_divB_total=%.17g",
           rank, patch_id, bench_sample_count, time_sim,
           time_sim * norm_config.Time_Norm,
           bench_final_max_divB_total);

  int bench_final_error = bench_sample_count != maximum_steps ||
                          bench_timer_error ||
                          bench_final_diagnostic_error;
  bench_final_error |=
      gamera_bench_write_profile(
          bench_samples, maximum_steps, bench_warmup_steps,
          bench_timed_steps, bench_initial_mi_s, bench_run_id,
          bench_manifest_sha256, bench_plan_sha256, bench_binary_sha256,
          bench_build_manifest_sha256, bench_omp_team_threads,
          time_sim, time_sim * norm_config.Time_Norm,
          bench_final_max_divB_total,
          bench_fault_rank == rank && bench_fault_kind != NULL &&
              strcmp(bench_fault_kind, "write") == 0) != 0;
  if (gamera_bench_consensus(bench_final_error)) {
    log_error("Incomplete or invalid per-rank fixed-Nt benchmark profile");
    return -1;
  }
#endif

  if (analysis_output_enabled && last_analysis_time != time_sim &&
      timed_analysis_output() != 0) {
    return -1;
  }
  if (last_restart_time != time_sim && timed_restart_output() != 0) {
    return -1;
  }
  return 0;
}
