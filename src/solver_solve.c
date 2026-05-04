//This file implements the main MHD (Magnetohydrodynamics) solver which
//advances the solution in time using AB2 method. It
//handles both the gas dynamics and magnetic field evolution, including:
//Conservative variable updates
//Boundary exchanges
//Electric field calculations
//Time step restrictions
//Data output

#include <mpi.h>
#include <math.h>
#include <omp.h>

#include "log.h"

#include "config.h"
#include "problem.h"
#include "setup_mpi.h"
#include "solver.h"
#include "utils.h"
#include "curvilinear.h"

int solve(int nt) {

  // exchange ghost cells
  int status =
      boundary_exchange_4d(gem, gem_onface_i, gem_onface_j, gem_onface_k);
  if (status != 0) {
    log_error("Error: boundary exchange for gem failed.");
    return -1;
  }

  status = boundary_exchange_5d(gas, NS1, config.NI, config.NJ, config.NK);
  if (status != 0) {
    log_error("Error: boundary exchange for gas failed.");
    return -1;
  }

  // Apply boundary conditions
  boundary_conditions();

  get_derived_variables();

  double last_dump_time = time_sim;
  if (!read_restart) {
    // dump the initial state
    log_info("Dump the initial state.");
    hdf_seq_num = -1;
    log_seq_num = 0;
    dump_data_hdf5(BASE_DAT_NAME);
    dump_extra_data(BASE_EXT_NAME);
    last_dump_time = 0;
  }

  if (!read_restart) {
    // get pdm limiting variables
    int d_index = gas_rho_p - gas_rho; // index difference between gas and gas_p
#pragma omp parallel for collapse(5) schedule(static)
    for (int s = 0; s < NS1; s++) {
      for (int m = gas_rho; m <= gas_p_S; m++) {
        for (int i = isg; i <= ieg; i++) {
          for (int j = jsg; j <= jeg; j++) {
            for (int k = ksg; k <= keg; k++) {
              gas[s][m + d_index][i][j][k] = gas[s][m][i][j][k];
            }
          }
        }
      }
    }
    d_index = mags_b1_p - mags_b1; // index difference between mags and mags_p
#pragma omp parallel for collapse(4) schedule(static)
    for (int m = mags_b1; m <= mags_b3; m++) {
      for (int i = isg; i <= ieg; i++) {
        for (int j = jsg; j <= jeg; j++) {
          for (int k = ksg; k <= keg; k++) {
            gem[m + d_index][i][j][k] = gem[m][i][j][k];
          }
        }
      }
    }

    d_index = mag_bi_p - mag_bi; // index difference between mag and mag_p
#pragma omp parallel for collapse(4) schedule(static)
    for (int m = mag_bi; m <= mag_bk; m++) {
      for (int i = isg; i <= ieg + 1; i++) {
        for (int j = jsg; j <= jeg + 1; j++) {
          for (int k = ksg; k <= keg + 1; k++) {
            gem[m + d_index][i][j][k] = gem[m][i][j][k];
          }
        }
      }
    }

    // Initialize dt0
    dt0 = get_dt();
    double dt_min_1;
    MPI_Allreduce(&dt0, &dt_min_1, 1, MPI_DOUBLE, MPI_MIN, comm_cart);
    dt0 = dt_min_1;
  }

  time_sim_start = time_sim;
  for (int n = 1; n <= nt; n++) {
    dt = get_dt();
    double dt_min;
    MPI_Allreduce(&dt, &dt_min, 1, MPI_DOUBLE, MPI_MIN, comm_cart);
    dt = dt_min;
    time_sim += dt;

    if (dtout < dt) {
      log_error(
          "Rank=%d, Something wrong with dt, dt may be too large. dtout=%lf, "
          "dt=%lf\n",
          rank, dtout, dt);
    }

    // AB2
    get_AB2_variables();

    int diff = gasc0_rho - gasc_rho;
#pragma omp parallel for collapse(5) schedule(static)
    for (int s = 0; s < NS1; s++) {
      for (int m = gasc_rho; m <= gasc_entropy; m++) {
        for (int i = is; i <= ie; i++) {
          for (int j = js; j <= je; j++) {
            for (int k = ks; k <= ke; k++) {
              // TODO: gasc0_* are intermediate variables
              // Store current gas state for RK stages
              gas[s][m + diff][i][j][k] = gas[s][m][i][j][k];
            }
          }
        }
      }
    }

    // Save magnetic field state
    diff = mag0_bi - mag_bi;
#pragma omp parallel for collapse(4) schedule(static)
    for (int m = mag_bi; m <= mag_bk; m++) {
      for (int i = isg; i <= ieg + 1; i++) {
        for (int j = jsg; j <= jeg + 1; j++) {
          for (int k = ksg; k <= keg + 1; k++) {
            // Store current magnetic field state for RK stages
            gem[m + diff][i][j][k] = gem[m][i][j][k] * 1.0;
          }
        }
      }
    }

    //  save velocity for Boris...
#pragma omp parallel for collapse(4) schedule(static)
    for (int s = 0; s < NS1; s++) {
      for (int i = is; i <= ie; i++) {
        for (int j = js; j <= je; j++) {
          for (int k = ks; k <= ke; k++) {
            // Cache velocity components for later use in Boris algorithm
            gas[s][v1_0][i][j][k] = gas[s][gas_v1][i][j][k];
            gas[s][v2_0][i][j][k] = gas[s][gas_v2][i][j][k];
            gas[s][v3_0][i][j][k] = gas[s][gas_v3][i][j][k];
          }
        }
      }
    }

    // Calculate electric fields for current RK stage
    get_e_fields();

    for (dir_t direction = i_dir; direction <= k_dir; direction++) {
      // gas --> gasl, gasr
      reconstruct_3dv_gas(direction);
      reconstruct_3dv_gem(direction);

      get_fluid_flux_rusanov(direction);
      get_magnetic_stress_rusanov(direction);

      Source_Term(direction);

      // i
      int diff = dgascf_rho - gascf_rho;

#pragma omp parallel for collapse(5) schedule(static)
      for (int s = 0; s < NS1; s++) {
        for (int m = gascf_rho; m <= gascf_entropy; m++) {
          for (int i = is; i <= ie; i++) {
            for (int j = js; j <= je; j++) {
              for (int k = ks; k <= ke; k++) {

                int m_dgascf = m + diff;

                switch (direction) {
                  case i_dir:
                    // TODO: m_dgascf is an intermediate variable
                    gas[s][m_dgascf][i][j][k] =
                        (gas[s][m][i + 1][j][k] *
                            geo[face_idir][i + 1][j][k] -
                            gas[s][m][i][j][k] * geo[face_idir][i][j][k]) /
                            geo[vol_center][i][j][k];
                    break;

                  case j_dir:
                    gas[s][m_dgascf][i][j][k] +=
                        (gas[s][m][i][j + 1][k] *
                            geo[face_jdir][i][j + 1][k] -
                            gas[s][m][i][j][k] * geo[face_jdir][i][j][k]) /
                            geo[vol_center][i][j][k];
                    break;

                  case k_dir:
                    gas[s][m_dgascf][i][j][k] +=
                        (gas[s][m][i][j][k + 1] *
                            geo[face_kdir][i][j][k + 1] -
                            gas[s][m][i][j][k] * geo[face_kdir][i][j][k]) /
                            geo[vol_center][i][j][k];
                    break;
                }
              }
            }
          }
        }
      }

      // same as (idir:kdir)
#pragma omp parallel for collapse(4) schedule(static)
      for (int m = magsf_x1dir; m <= magsf_x3dir; m++) {
        for (int i = is; i <= ie; i++) {
          for (int j = js; j <= je; j++) {
            for (int k = ks; k <= ke; k++) {
              int diff = dmagsf_x1dir - magsf_x1dir;

              switch (direction) {
                case i_dir:
                  gem[m + diff][i][j][k] =
                      (gem[m][i + 1][j][k] * geo[face_idir][i + 1][j][k] -
                          gem[m][i][j][k] * geo[face_idir][i][j][k]) /
                          geo[vol_center][i][j][k];
                  break;

                case j_dir:
                  gem[m + diff][i][j][k] +=
                      (gem[m][i][j + 1][k] * geo[face_jdir][i][j + 1][k] -
                          gem[m][i][j][k] * geo[face_jdir][i][j][k]) /
                          geo[vol_center][i][j][k];
                  break;

                case k_dir:
                  gem[m + diff][i][j][k] +=
                      (gem[m][i][j][k + 1] * geo[face_kdir][i][j][k + 1] -
                          gem[m][i][j][k] * geo[face_kdir][i][j][k]) /
                          geo[vol_center][i][j][k];
                  break;
              }

            }
          }
        }
      }
    }  // end idir - kdir

    if (doRingAverage == 1) {
      tackle_Efield_pole();
    }

    // update magnetic field...
#pragma omp parallel for collapse(3) schedule(static)
    for (int i = is; i <= ie + 1; i++) {
      for (int j = js; j <= je; j++) {
        for (int k = ks; k <= ke; k++) {
          gem[mag_bi][i][j][k] -=
              dt *
                  ((gem[efield_ek][i][j + 1][k] * geo[edge_kdir][i][j + 1][k] -
                      gem[efield_ek][i][j][k] * geo[edge_kdir][i][j][k]) -
                      (gem[efield_ej][i][j][k + 1] * geo[edge_jdir][i][j][k + 1] -
                          gem[efield_ej][i][j][k] * geo[edge_jdir][i][j][k])) /
                  geo[face_idir][i][j][k];
        }
      }
    }

#pragma omp parallel for collapse(3) schedule(static)
    for (int i = is; i <= ie; i++) {
      for (int j = js; j <= je + 1; j++) {
        for (int k = ks; k <= ke; k++) {
          gem[mag_bj][i][j][k] -=
              dt *
                  ((gem[efield_ei][i][j][k + 1] * geo[edge_idir][i][j][k + 1] -
                      gem[efield_ei][i][j][k] * geo[edge_idir][i][j][k]) -
                      (gem[efield_ek][i + 1][j][k] * geo[edge_kdir][i + 1][j][k] -
                          gem[efield_ek][i][j][k] * geo[edge_kdir][i][j][k])) /
                  geo[face_jdir][i][j][k];
        }
      }
    }

#pragma omp parallel for collapse(3) schedule(static)
    for (int i = is; i <= ie; i++) {
      for (int j = js; j <= je; j++) {
        for (int k = ks; k <= ke + 1; k++) {
          gem[mag_bk][i][j][k] -=
              dt *
                  ((gem[efield_ej][i + 1][j][k] * geo[edge_jdir][i + 1][j][k] -
                      gem[efield_ej][i][j][k] * geo[edge_jdir][i][j][k]) -
                      (gem[efield_ei][i][j + 1][k] * geo[edge_idir][i][j + 1][k] -
                          gem[efield_ei][i][j][k] * geo[edge_idir][i][j][k])) /
                  geo[face_kdir][i][j][k];
        }
      }
    }

    if (doRingAverage == 1) {
      tackle_Magfield_pole();
    }

    // add the lorentz force terms from the background field B0 - not the
    // deltas are -dt*volume_integral./volume, only need active cells
#pragma omp parallel for collapse(4) schedule(static)
    for (int m = dmagsf_G_x1dir; m <= dmagsf_G_x3dir; m++) {
      for (int i = is; i <= ie; i++) {
        for (int j = js; j <= je; j++) {
          for (int k = ks; k <= ke; k++) {
            int diff = dmagsf_x1dir - dmagsf_G_x1dir;
            gem[m + diff][i][j][k] += gem[m][i][j][k];
          }
        }
      }
    }

#pragma omp parallel for collapse(4) schedule(static)
    for (int s = 0; s < NS1; s++) {
      for (int i = is; i <= ie; i++) {
        for (int j = js; j <= je; j++) {
          for (int k = ks; k <= ke; k++) {
            gas[s][gasc_rho][i][j][k] -= dt * gas[s][dgascf_rho][i][j][k];
            gas[s][gasc_rhov1][i][j][k] = gas[s][gasc_rhov1][i][j][k] -
                dt * gas[s][dgascf_rhov1][i][j][k] +
                dt * gas[s][S1_HD][i][j][k];
            gas[s][gasc_rhov2][i][j][k] = gas[s][gasc_rhov2][i][j][k] -
                dt * gas[s][dgascf_rhov2][i][j][k] +
                dt * gas[s][S2_HD][i][j][k];
            gas[s][gasc_rhov3][i][j][k] = gas[s][gasc_rhov3][i][j][k] -
                dt * gas[s][dgascf_rhov3][i][j][k] +
                dt * gas[s][S3_HD][i][j][k];
            gas[s][gasc_eng][i][j][k] -= dt * gas[s][dgascf_eng][i][j][k];
            gas[s][gasc_entropy][i][j][k] -= dt * gas[s][dgascf_entropy][i][j][k];
          }
        }
      }
    }

    // now, conservative --> primitive
#pragma omp parallel for collapse(4) schedule(static)
    for (int s = 0; s < NS1; s++) {
      for (int i = is; i <= ie; i++) {
        for (int j = js; j <= je; j++) {
          for (int k = ks; k <= ke; k++) {
            gas[s][gas_rho][i][j][k] = gas[s][gasc_rho][i][j][k];

            gas[s][gas_v1][i][j][k] =
                gas[s][gasc_rhov1][i][j][k] / gas[s][gasc_rho][i][j][k];
            gas[s][gas_v2][i][j][k] =
                gas[s][gasc_rhov2][i][j][k] / gas[s][gasc_rho][i][j][k];
            gas[s][gas_v3][i][j][k] =
                gas[s][gasc_rhov3][i][j][k] / gas[s][gasc_rho][i][j][k];

            gas[s][gas_p][i][j][k] = (gas[s][gasc_eng][i][j][k] -
                0.5 *
                    (gas[s][gas_v1][i][j][k] *
                        gas[s][gas_v1][i][j][k] +
                        gas[s][gas_v2][i][j][k] *
                            gas[s][gas_v2][i][j][k] +
                        gas[s][gas_v3][i][j][k] *
                            gas[s][gas_v3][i][j][k]) *
                    gas[s][gas_rho][i][j][k]) *
                (gamma_val - 1.0);
            gas[s][gas_p_S][i][j][k] = gas[s][gasc_entropy][i][j][k] * pow(gas[s][gasc_rho][i][j][k], gamma_val - 1.0);
          }
        }
      }
    }

//    // reset the rho and pressure
//    reset_rho();
//    reset_p();
    check_positivity();

    //  calculate the magtot in immediate step
#pragma omp parallel for collapse(3) schedule(static)
    for (int i = is; i <= ie; i++) {
      for (int j = js; j <= je; j++) {
        for (int k = ks; k <= ke; k++) {
          gem[mags_b1][i][j][k] = gem[mag_bi][i][j][k] * (x1[i + 1][j][k] - x1ctr[i][j][k]) / dx1[i][j][k]
              + gem[mag_bi][i + 1][j][k] * (x1ctr[i][j][k] - x1[i][j][k]) / dx1[i][j][k];
          gem[mags_b2][i][j][k] = gem[mag_bj][i][j][k] * (x2[i][j + 1][k] - x2ctr[i][j][k]) / dx2[i][j][k]
              + gem[mag_bj][i][j + 1][k] * (x2ctr[i][j][k] - x2[i][j][k]) / dx2[i][j][k];
          gem[mags_b3][i][j][k] = gem[mag_bk][i][j][k] * (x3[i][j][k + 1] - x3ctr[i][j][k]) / dx3[i][j][k]
              + gem[mag_bk][i][j][k + 1] * (x3ctr[i][j][k] - x3[i][j][k]) / dx3[i][j][k];
        }
      }
    }

#pragma omp parallel for collapse(3) schedule(static)
    for (int i = is; i <= ie; i++) {
      for (int j = js; j <= je; j++) {
        for (int k = ks; k <= ke; k++) {
          gem[magtot_b1][i][j][k] = 0.5 * (gem[mags_b1][i][j][k] + gem[mags_b1_p][i][j][k]) + gem[B0_b1][i][j][k];
          gem[magtot_b2][i][j][k] = 0.5 * (gem[mags_b2][i][j][k] + gem[mags_b2_p][i][j][k]) + gem[B0_b2][i][j][k];
          gem[magtot_b3][i][j][k] = 0.5 * (gem[mags_b3][i][j][k] + gem[mags_b3_p][i][j][k]) + gem[B0_b3][i][j][k];
        }
      }
    }

    int s = NS1 - 1;
#pragma omp parallel for collapse(3) schedule(static)
    for (int i = is; i <= ie; i++) {
      for (int j = js; j <= je; j++) {
        for (int k = ks; k <= ke; k++) {
          //  first estimate the Alfven constant for perp momentum
          gem[alf_ratio][i][j][k] =
              (gem[magtot_b1][i][j][k] * gem[magtot_b1][i][j][k] +
                  gem[magtot_b2][i][j][k] * gem[magtot_b2][i][j][k] +
                  gem[magtot_b3][i][j][k] * gem[magtot_b3][i][j][k]) /
                  gas[s][gas_rho][i][j][k] / (CA * CA);
          // use the updated rho at T(n+1) for the alfven correction
          gem[perp_ratio][i][j][k] = 1 / (1 + gem[alf_ratio][i][j][k]);
          // then update the momentum at t(N+1) with Alfven correction
          gem[dv_alf][i][j][k] =
              gem[alf_ratio][i][j][k] * (-dt * gas[s][dgascf_rho][i][j][k]);
          gas[s][v1_tmp][i][j][k] =
              gem[alf_ratio][i][j][k] *
                  (-dt * gas[s][dgascf_rhov1][i][j][k] + dt * gas[s][S1_HD][i][j][k] -
                      (-dt * gas[s][dgascf_rho][i][j][k]) * gas[s][v1_0][i][j][k]);
          gas[s][v2_tmp][i][j][k] =
              gem[alf_ratio][i][j][k] *
                  (-dt * gas[s][dgascf_rhov2][i][j][k] + dt * gas[s][S2_HD][i][j][k] -
                      (-dt * gas[s][dgascf_rho][i][j][k]) * gas[s][v2_0][i][j][k]);
          gas[s][v3_tmp][i][j][k] =
              gem[alf_ratio][i][j][k] *
                  (-dt * gas[s][dgascf_rhov3][i][j][k] + dt * gas[s][S3_HD][i][j][k] -
                      (-dt * gas[s][dgascf_rho][i][j][k]) * gas[s][v3_0][i][j][k]);
          gas[s][bdotv][i][j][k] =
              gem[magtot_b1][i][j][k] * gas[s][v1_tmp][i][j][k] +
                  gem[magtot_b2][i][j][k] * gas[s][v2_tmp][i][j][k] +
                  gem[magtot_b3][i][j][k] * gas[s][v3_tmp][i][j][k];
          gas[s][bdotv][i][j][k] =
              gas[s][bdotv][i][j][k] /
                  (gem[magtot_b1][i][j][k] * gem[magtot_b1][i][j][k] +
                      gem[magtot_b2][i][j][k] * gem[magtot_b2][i][j][k] +
                      gem[magtot_b3][i][j][k] * gem[magtot_b3][i][j][k] + 1e-10);
          // if CA>>1, then rhovx = rhovx + dpxF + dpxB
          gas[s][gasc_rhov1][i][j][k] = gas[s][gasc0_rhov1][i][j][k] +
              gem[perp_ratio][i][j][k] *
                  ((-dt * gas[s][dgascf_rhov1][i][j][k] + dt * gas[s][S1_HD][i][j][k]) +
                      (-dt * gem[dmagsf_x1dir][i][j][k] + dt * gem[S1_B][i][j][k]) +
                      gem[dv_alf][i][j][k] * gas[s][v1_0][i][j][k] +
                      gem[magtot_b1][i][j][k] * gas[s][bdotv][i][j][k]);
          // rhovy = rhovy + dpyF + dpyB
          gas[s][gasc_rhov2][i][j][k] = gas[s][gasc0_rhov2][i][j][k] +
              gem[perp_ratio][i][j][k] *
                  ((-dt * gas[s][dgascf_rhov2][i][j][k] + dt * gas[s][S2_HD][i][j][k]) +
                      (-dt * gem[dmagsf_x2dir][i][j][k] + dt * gem[S2_B][i][j][k]) +
                      gem[dv_alf][i][j][k] * gas[s][v2_0][i][j][k] +
                      gem[magtot_b2][i][j][k] * gas[s][bdotv][i][j][k]);
          // rhovz = rhovz + dpzF + dpzB
          gas[s][gasc_rhov3][i][j][k] = gas[s][gasc0_rhov3][i][j][k] +
              gem[perp_ratio][i][j][k] *
                  ((-dt * gas[s][dgascf_rhov3][i][j][k] + dt * gas[s][S3_HD][i][j][k]) +
                      (-dt * gem[dmagsf_x3dir][i][j][k] + dt * gem[S3_B][i][j][k]) +
                      gem[dv_alf][i][j][k] * gas[s][v3_0][i][j][k] +
                      gem[magtot_b3][i][j][k] * gas[s][bdotv][i][j][k]);

          gas[s][gasc_rhov1][i][j][k] += dt * gas[s][gasc_rho][i][j][k] * G01;
          gas[s][gasc_rhov2][i][j][k] += dt * gas[s][gasc_rho][i][j][k] * G02;
          gas[s][gasc_rhov3][i][j][k] += dt * gas[s][gasc_rho][i][j][k] * G03;


          // get velocities with magnetic stress - now bulk vx, vy, vz are
          // solved
          gas[s][gas_v1][i][j][k] =
              gas[s][gasc_rhov1][i][j][k] / gas[s][gasc_rho][i][j][k];
          gas[s][gas_v2][i][j][k] =
              gas[s][gasc_rhov2][i][j][k] / gas[s][gasc_rho][i][j][k];
          gas[s][gas_v3][i][j][k] =
              gas[s][gasc_rhov3][i][j][k] / gas[s][gasc_rho][i][j][k];
        }
      }
    }

    if (doRingAverage == 1) {
      gather_value_RingAverage_arrays();
      HydroRingAverage();
      MagneticRingAverage();
      broadcast_value_RingAverage_arrays();
    }

    status =
        boundary_exchange_4d(gem, gem_onface_i, gem_onface_j, gem_onface_k);
    if (status != 0) {
      log_error("Error: boundary exchange for gem failed.");
      return -1;
    }
    status = boundary_exchange_5d(gas, NS1, config.NI, config.NJ, config.NK);
    if (status != 0) {
      log_error("Error: boundary exchange for gas failed.");
      return -1;
    }

    // Apply boundary conditions
    boundary_conditions();

    get_derived_variables();


    // Output data:
    if (time_sim >= output_interval * (hdf_seq_num + 1)) {
      dump_data_hdf5(BASE_DAT_NAME);
      log_info("Output data at time %.15f\n", time_sim);
      last_dump_time = time_sim;
    }

    if (time_sim >= 0.1 * output_interval * (log_seq_num)) {
      log_info("Status at step %d: rank %d, Time %.15f, dt %.15f, max_divB for this rank is %.4e\n",
               n, rank, time_sim, dt, divB_max);
      log_seq_num = log_seq_num + 1;
    }

    if (time_sim > time_stop) {
      log_info("Step %d, Time %f, dt %f\n", n, time_sim, dt);
      log_info("This is the end\n");
      break;  // Exits the loop
    }
  }  // end do n=1:Nt

  if (last_dump_time != time_sim) {
    dump_data_hdf5(BASE_DAT_NAME);
  }

  return 0;
}
