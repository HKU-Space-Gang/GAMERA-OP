#include <mpi.h>
#include <stdio.h>
#include <unistd.h>

#include "log.h"

#include "config.h"
#include "problem.h"
#include "setup_mpi.h"
#include "solver.h"
#include "curvilinear.h"

// the main funcation that accepts mpi configs and calls the solver
int main(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <config_file>\n", argv[0]);
    return 1;
  }
  const char *config_file = argv[1];

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
    goto done;
  }

  // initialize reconstruction weights
  getweights();
  if (!read_restart) {
    // variables are not initialized yet, initialize them
    log_info("Initializing problem");
    problem_init();
  }
// set constant background magnetic field related variables
  set_background_field();
  // ---------------- end of mpi setup ----------------------------

  // ring average config
  if (doRingAverage == 1) {
    set_RingAverage_geo_data();
    set_RingAverage_config();
  }

  // the main loop
  dtout = 10.0;

  status = solve(Nt);
  if (status != 0) {
    log_error("Error: solver failed.");
  }

  log_info("End of simulation: rank %d, time %.15f, time simulated %.15f\n",
           rank, time_sim, time_sim - time_sim_start);

  end_time = clock();
  total_time = (double) (end_time - start_time) / CLOCKS_PER_SEC;

  log_info("Total running time = %f second, time per step = %f \n", total_time,
           total_time / Nt);

  done:
  finalize_solver();
  // finalize mpi
  free_mpi_types();
  log_info("Rank %d finalized.\n", rank);
  fclose(log_file);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
