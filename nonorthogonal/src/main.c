#include <mpi.h>
#include <stdio.h>
#include <unistd.h>

#include "log.h"

#include "config.h"
#include "problem.h"
#include "setup_mpi.h"
#include "solver.h"
#include "curvilinear.h"
#ifdef GAMERA_NONORTHOGONAL_BACKEND
#include "nonorthogonal_driver.h"
#include "nonorthogonal_legacy_adapter.h"
#ifdef GAMERA_YINYANG_BACKEND
#include "nonorthogonal_yinyang_exchange.h"
#endif
#endif

// the main funcation that accepts mpi configs and calls the solver
int main(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <config_file>\n", argv[0]);
    return 1;
  }
  const char *config_file = argv[1];
  int exit_status = 0;

  double total_time;
  clock_t start_time, end_time;
  start_time = clock();
  // initialize mpi
  MPI_Init(&argc, &argv);
  // get rank and size
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  set_problem_config();

  if (initialize_config(rank, size, config_file) != 0) {
    log_error("Failed to initialize config.");
    return -1;
  }

  int status;
  status = setup_comm_and_types(rank);
  if (status != 0) {
    log_error("Failed to setup cartesian communicator or types.");
    return -1;
  }

  if (initialize_logging(BASE_LOG_NAME, rank) != 0) {
    log_error("Failed to initialize logging.");
    return -1;
  }

//  extra_operations_start();
  if (initialize_solver(BASE_DAT_NAME) != 0) {
    log_error("Failed to initialize solver.");
    exit_status = 1;
    goto done;
  }
#ifdef GAMERA_TIME_DEPENDENT_WIND
  if (problem_runtime_init() != 0) {
    log_error("Failed to initialize time-dependent solar-wind boundary");
    exit_status = 1;
    goto done;
  }
#endif

  // initialize reconstruction weights
#ifndef GAMERA_NONORTHOGONAL_BACKEND
  getweights();
#endif
  if (!read_restart) {
    // variables are not initialized yet, initialize them
    log_info("Initializing problem");
    problem_init();
#ifdef GAMERA_NONORTHOGONAL_BACKEND
    if (gamera_no_legacy_grid() == NULL ||
        gamera_no_legacy_storage() == NULL) {
      log_error("Failed to initialize non-orthogonal problem state");
      exit_status = 1;
      goto done;
    }
#endif
  }
#ifdef GAMERA_NONORTHOGONAL_BACKEND
  else if (gamera_no_legacy_adapter_create(true) != 0) {
    log_error("Failed to import non-orthogonal restart state");
    exit_status = 1;
    goto done;
  }
#else
// set constant background magnetic field related variables
  set_background_field();
#endif
#ifdef GAMERA_NONORTHOGONAL_BACKEND
  if (gamera_no_driver_prepare() != 0) {
    log_error("Failed to prepare non-orthogonal driver data");
    exit_status = 1;
    goto done;
  }
#endif
  // ---------------- end of mpi setup ----------------------------

  // ring average config
  if (doRingAverage == 1) {
    set_RingAverage_geo_data();
    set_RingAverage_config();
  }

  // the main loop
  dtout = 10.0;

#ifdef GAMERA_NONORTHOGONAL_BACKEND
  status = solve_nonorthogonal(Nt);
#else
  status = solve(Nt);
#endif
  if (status != 0) {
    log_error("Error: solver failed.");
    exit_status = 1;
  }

  log_info("End of simulation: rank %d, time %.15f, time simulated %.15f\n",
           rank, time_sim, time_sim - time_sim_start);

  end_time = clock();
  total_time = (double) (end_time - start_time) / CLOCKS_PER_SEC;

  log_info("Total running time = %f second, time per step = %f \n", total_time,
           total_time / Nt);

  done:
#ifdef GAMERA_TIME_DEPENDENT_WIND
  problem_runtime_finalize();
#endif
#ifdef GAMERA_NONORTHOGONAL_BACKEND
  gamera_no_driver_finalize();
#ifdef GAMERA_YINYANG_BACKEND
  gamera_no_yinyang_exchange_destroy();
#endif
  gamera_no_legacy_adapter_destroy();
#endif
  finalize_solver();
  // finalize mpi
  free_mpi_types();
  log_info("Rank %d finalized.\n", rank);
  fclose(log_file);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return exit_status;
}
