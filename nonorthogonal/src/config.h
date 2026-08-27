#ifndef CONFIG_H
#define CONFIG_H

#include <stdbool.h>
#include <stdio.h>

// #define MULTIFLUID

// NG is the number of ghost points
#define NG 4
#define NO2 NG
#define NO (NO2 * 2)
#define PI 3.14159265358979323846

#ifdef MULTIFLUID
#define NS 2
#define NS1 (NS + 1)
#else
#define NS 1
#define NS1 NS
#endif

#define out_prefix_pattern "%s_%02d-%02d-%02d_"
#define out_file_pattern (out_prefix_pattern "%06d.%s")
#define BASE_DAT_NAME "mhd"
#define BASE_ANALYSIS_NAME "analysis"
#define BASE_RESTART_NAME "restart"
//#define BASE_EXT_NAME "extra"
#define BASE_EXT_NAME "mhd"
#define BASE_LOG_NAME "LOG"

// Boundary conditions
typedef enum {
  BC_NON_PERIODIC = 1,
  BC_PERIODIC
} Boundary_condition_t;

// ni, nj, nk are the number of active grid points in each dimension
typedef struct {
  int ni_global;
  int nj_global;
  int nk_global;
  // (active center size)
  int ni;
  int nj;
  int nk;
  // NI, NJ, NK are the grid size in each dimension(conner size)
  int NI;
  int NJ;
  int NK;
  int proc_dims[3];
  int proc_coords[3];
} config_t;

extern config_t config;
extern int *proc_coords;

// sequence number of the simulation
extern int hdf_seq_num;
extern int analysis_seq_num;
extern int log_seq_num;
extern bool read_restart;
#define GAMERA_RESTART_PATH_MAX 1024
/* Rank-local MHD restart file selected by solver initialization. Rank zero's
 * path also carries the collective electrostatic M-I restart group. */
extern char restart_source_filename[GAMERA_RESTART_PATH_MAX];
enum { GAMERA_YINYANG_HDIV_SCHEME_VERSION = 3 };
/* H(div) checkpoint scheme version.  Legacy files leave this zero; version 1
 * used a unit flux norm, version 2 projected the full residual field with an
 * A/d Hodge, and version 3 preserves background-field cancellation by
 * projecting only the Yin-Yang donor increment. */
extern int restart_has_yinyang_hdiv;

typedef struct {
  Boundary_condition_t Boundary_i;
  Boundary_condition_t Boundary_j;
  Boundary_condition_t Boundary_k;
} Problem_config_t;

typedef struct {
  double mu0;
  double x_Norm;
  double u_Norm;
  double Time_Norm;
  double rho_Norm;
  double p_Norm;
  double B_Norm;
  double Magnetic_moment_Norm;
  double Omega_Norm;
  double Sigma_Norm;
} Norm_t;
extern Norm_t norm_config;

#define GAMERA_WIND_PATH_MAX 1024
typedef struct {
  int enabled;
  char file[GAMERA_WIND_PATH_MAX];
  /* 0: already normalized code units; 1: Kaiju/space-weather physical units. */
  int physical_units;
  int linear_interpolation;
  double velocity_si_scale;
  double reference[3];
  double time_offset;
  int coefficient_override;
  double by_coefficient;
  double bz_coefficient;
  double bx_offset;
  /* -1: infer from ByC/BzC, 0: retain input Bx, 1: enforce Lyon relation. */
  int enforce_bx_relation;
} Wind_config_t;
extern Wind_config_t wind_config;

typedef struct {
  int enabled;
  int diagnostics_enabled;
  /* Fedder-style MHD-to-electron precipitation followed by Robinson
   * conductance.  This is independent of whether the variable-conductance
   * tensor is fed back into the electrostatic solve. */
  int electron_precipitation_enabled;
  int electron_precipitation_hall_model;
  /* Fedder-1995 thermal number-flux coefficient.  Keeping this explicit in
   * the run configuration makes calibrated precipitation reproducible. */
  double electron_precipitation_beta;
  double electron_precipitation_ramp_s;
  int electron_precipitation_conductance_feedback_enabled;
  double electron_precipitation_conductance_smoothing_s;
  double euv_f107;
  /* Unit need not be normalized in input.  This is the direction from Earth
   * toward the Sun in the model's native Cartesian coordinates. */
  double sunward_direction[3];
  int hybrid_dpb_enabled;
  double dpb_absolute_potential_v;
  double dpb_adaptive_potential_fraction;
  double dpb_adaptive_potential_floor_v;
  double dpb_potential_equatorward_offset_deg;
  double dpb_fac_equatorward_offset_deg;
  double dpb_fac_absolute_floor_a_m2;
  double dpb_minimum_boundary_latitude_deg;
  double dpb_maximum_candidate_latitude_deg;
  double dpb_transition_width_deg;
  double dpb_temporal_timescale_s;
  double dpb_maximum_slew_deg_per_s;
  double dpb_quiet_initial_nightside_latitude_deg;
  double dpb_oval_center_latitude_deg;
  double dpb_oval_center_mlt_h;
  double dpb_nightside_poleward_limit_deg;
  double minimum_pedersen_siemens;
  double minimum_hall_siemens;
  double maximum_hall_to_pedersen;
  double coupling_interval_s;
  double diagnostics_interval_s;
  double ionosphere_radius_re;
  double pedersen_siemens;
  double hall_siemens;
  int current_shell; /* zero-based active radial cell */
  /* Linear refinement of an automatically derived ionospheric grid relative
   * to the nominal Yin--Yang angular cell spacing.  Explicit longitude and
   * colatitude counts continue to take precedence. */
  double angular_oversampling;
  int longitude_count;
  int colatitude_count;
  int maximum_iterations;
  double relative_tolerance;
  double absolute_tolerance;
} MI_config_t;
extern MI_config_t mi_config;

#define GAMERA_IM_FEEDBACK_PATH_MAX 1024
typedef struct {
  int enabled;
  char provider[64];
  char prescribed_state_directory[GAMERA_IM_FEEDBACK_PATH_MAX];
  char prescribed_map_directory[GAMERA_IM_FEEDBACK_PATH_MAX];
  char provenance_sha256[65];
  int prescribed_end_generation;
  double coupling_interval_s;
  double activation_time_s;
  double relaxation_time_s;
  double matching_time_s;
  double maximum_fractional_pressure_change;
  double pressure_floor_pa;
  double ion_pressure_fraction;
  double taper_full_l;
  double taper_zero_l;
} IM_feedback_config_t;
extern IM_feedback_config_t im_feedback_config;

/*
 * Online H+ provider seam.  By itself it advances and exposes a diagnostic
 * kinetic target without changing MHD.  The separate, default-OFF live
 * pressure-feedback option may consume the same committed snapshots.  An
 * enabled provider must register an explicit deterministic mapped-shell
 * callback at runtime.
 */
typedef struct {
  int enabled;
  char provider[64];
  double coupling_interval_s;
  double activation_time_s;
  double cfl_fraction;
  double ion_pressure_fraction;
} IM_online_hplus_config_t;
extern IM_online_hplus_config_t im_online_hplus_config;

/*
 * Live pressure-only feedback is intentionally configured independently from
 * the retired prescribed-state verification path.  Its cadence, activation,
 * and H+ ion-pressure fraction are inherited from im_online_hplus_config so
 * one online provider cannot be advanced on two subtly different clocks.
 */
typedef struct {
  int enabled;
  double relaxation_time_s;
  double matching_time_s;
  double maximum_fractional_pressure_change;
  double pressure_floor_pa;
  double taper_full_l;
  double taper_zero_l;
} IM_online_hplus_feedback_config_t;
extern IM_online_hplus_feedback_config_t im_online_hplus_feedback_config;

extern int Nt;
extern double time_stop;
extern double output_interval;
/* Compact one-file-per-patch analysis output is the default.  Setting it to
 * zero disables routine analysis; it does not restore the retired per-rank
 * analysis stream.  restart_interval is expressed in code time; zero means
 * initial/final checkpoints only. */
extern int analysis_output_enabled;
extern double restart_interval;
//extern int restart;
extern config_t config;
extern FILE *log_file;
extern int pdm_index_gas;
extern int pdm_index_gem;

// Problem configuration
extern Problem_config_t problem_config;
// Domain boundaries
extern double x1min_global, x1max_global;
extern double x2min_global, x2max_global;
extern double x3min_global, x3max_global;
extern double x1min, x1max;
extern double x2min, x2max;
extern double x3min, x3max;

// Ringavarage option
extern int doRingAverage;

// 0 keeps the pole machinery but makes every azimuthal chunk one cell
// (identity/no coarsening); 1 uses the standard Ring Average chunk schedule.
extern int ringAverageMode;

// Define whether this rank is in the region of ringaverage
extern int doRingAverage_thisrank_positivepole;
extern int doRingAverage_thisrank_negativepole;

// Nchunck store number of chunk in ring average direction
extern int Nchunk[32];

// NAVR is the number of belts should be averaged
extern int NAVR;

// Physical constants
extern const double gamma_val;
extern const double CFL;
extern const double PDMB;
// extern const double time_stop;
extern double CA;
extern const double nlf_cs;
extern const double nlf_va;
// Set floor value for density and pressure
extern double rho_floor;
extern double p_floor;

// Gravitational field
extern double G01;
extern double G02;
extern double G03;

// Divergence of magnetic field
extern double divB_max;

// Factor arrays (Direct initialization)
extern const double fac8th[8];
extern const double facPDMU7[7];

// Grid indices
// is = isg+NO2, ie = is+Nx1 -1, isg = 0, ieg = ie + NO2,
extern int is, ie, isg, ieg;
extern int js, je, jsg, jeg;
extern int ks, ke, ksg, keg;

// Time step
extern double dtout, dt, dt0, dxx, dyy, dzz, time_sim, time_sim_start;

// standalone 3D arrays
extern double ***x1, ***x2, ***x3;
extern double ***x1c, ***x2c, ***x3c;
extern double ***x1ctr, ***x2ctr, ***x3ctr;
extern double ***dx1, ***dx2, ***dx3;


typedef enum {
  i_dir = 0,
  j_dir,
  k_dir,
} dir_t;
// 4D constant arrays
// constant geometry arrays
extern double ****geo;
enum {
  edge_idir = 0,
  edge_jdir,
  edge_kdir,
  dxi,
  dxj,
  dxk,
  dxi_ring,
  dxj_ring,
  dxk_ring,
  face_idir,
  face_jdir,
  face_kdir,
  vol_idir,
  vol_jdir,
  vol_kdir,
  vol_center,
  GF_face_idir,
  GF_face_jdir,
  GF_face_kdir,
  GF_vol_idir,
  GF_vol_jdir,
  GF_vol_kdir,
  IfaceAedgeJ,
  IfaceAedgeK,
  JfaceAedgeK,
  JfaceAedgeI,
  KfaceAedgeI,
  KfaceAedgeJ,
  NF_geo
};
extern int geo_onface_i[NF_geo];
extern int geo_onface_j[NF_geo];
extern int geo_onface_k[NF_geo];

// reconstruction related arrays
extern double ****rec;
enum {
  // left states
  // first direction
  w1l_first_direction = 0,
  w2l_first_direction,
  w3l_first_direction,
  w4l_first_direction,
  w5l_first_direction,
  w6l_first_direction,
  w7l_first_direction,
  // second direction
  w1l_second_direction,
  w2l_second_direction,
  w3l_second_direction,
  w4l_second_direction,
  w5l_second_direction,
  w6l_second_direction,
  w7l_second_direction,
  // third direction
  w1l_third_direction,
  w2l_third_direction,
  w3l_third_direction,
  w4l_third_direction,
  w5l_third_direction,
  w6l_third_direction,
  w7l_third_direction,
  // right states
  // first direction
  w1r_first_direction,
  w2r_first_direction,
  w3r_first_direction,
  w4r_first_direction,
  w5r_first_direction,
  w6r_first_direction,
  w7r_first_direction,
  // second direction
  w1r_second_direction,
  w2r_second_direction,
  w3r_second_direction,
  w4r_second_direction,
  w5r_second_direction,
  w6r_second_direction,
  w7r_second_direction,
  // third direction
  w1r_third_direction,
  w2r_third_direction,
  w3r_third_direction,
  w4r_third_direction,
  w5r_third_direction,
  w6r_third_direction,
  w7r_third_direction,
  NF_rec
};

// 4D variable arrays for geometry, electric and magnetic fields which need
// exchanging geometry, electric and magnetic fields
extern double ****gem;
enum {
  // mag
  mag_bi = 0,
  mag_bj,
  mag_bk,

  // mags: even boundary in i, j, k
  mags_b1,
  mags_b2,
  mags_b3,

  // mag_h to store the half step value of mag
  mag_bi_h,
  mag_bj_h,
  mag_bk_h,
  // mag_p to store the beginning value of mag for PDMU limiters
  mag_bi_p,
  mag_bj_p,
  mag_bk_p,

  // mags_h to store the half step value of mag
  mags_b1_h,
  mags_b2_h,
  mags_b3_h,
  // mags_p to store the beginning value of mag for PDMU limiters
  mags_b1_p,
  mags_b2_p,
  mags_b3_p,

  bi_avg,
  bj_avg,
  bk_avg,
  bil,
  bir,
  bjl,
  bjr,
  bkl,
  bkr,
  dvzz,

  // alf_ratio
  alf_ratio,
  // perp_ratio
  perp_ratio,
  // dv_alf
  dv_alf,

  // mag0
  mag0_bi,
  mag0_bj,
  mag0_bk,
  // magsl
  magsl_b1,
  magsl_b2,
  magsl_b3,
  // magsr
  magsr_b1,
  magsr_b2,
  magsr_b3,
  // magsf
  magsf_x1dir,
  magsf_x2dir,
  magsf_x3dir,
  // dmagsf
  dmagsf_x1dir,
  dmagsf_x2dir,
  dmagsf_x3dir,
  // efield
  efield_ei,
  efield_ej,
  efield_ek,
  // LA*
  LAi,
  LAj,
  LAk,
  // B0
  B0_b1,
  B0_b2,
  B0_b3,
  // B0AfaceI
  B0AfaceI_b1,
  B0AfaceI_b2,
  B0AfaceI_b3,
  B0AfaceI_bn,
  B0AfaceI_bsq,
  B0AfaceI_bn1,
  B0AfaceI_bn2,
  B0AfaceI_bn3,
  // B0AfaceJ
  B0AfaceJ_b1,
  B0AfaceJ_b2,
  B0AfaceJ_b3,
  B0AfaceJ_bn,
  B0AfaceJ_bsq,
  B0AfaceJ_bn1,
  B0AfaceJ_bn2,
  B0AfaceJ_bn3,
  // B0AfaceK
  B0AfaceK_b1,
  B0AfaceK_b2,
  B0AfaceK_b3,
  B0AfaceK_bn,
  B0AfaceK_bsq,
  B0AfaceK_bn1,
  B0AfaceK_bn2,
  B0AfaceK_bn3,
  // B0Aedge_I
  B0AedgeI_b1,
  B0AedgeI_b2,
  B0AedgeI_b3,
  // B0Aedge_J
  B0AedgeJ_b1,
  B0AedgeJ_b2,
  B0AedgeJ_b3,
  // B0Aedge_K
  B0AedgeK_b1,
  B0AedgeK_b2,
  B0AedgeK_b3,
  // dmagsf_G
  dmagsf_G_x1dir,
  dmagsf_G_x2dir,
  dmagsf_G_x3dir,
  // source term due to mag
  S1_B,
  S2_B,
  S3_B,
  // magtot
  magtot_b1,
  magtot_b2,
  magtot_b3,
  NF_gem
};

extern int gem_onface_i[NF_gem];
extern int gem_onface_j[NF_gem];
extern int gem_onface_k[NF_gem];
extern int NF_gem_prim;
extern int NF_gem_store;

// 5D arrays
// species related
extern double *****gas;
enum {
  // gas: even boundary in i, j, k
  gas_rho = 0,
  gas_v1,
  gas_v2,
  gas_v3,
  gas_p,
  gas_p_S,

  // gas_h to store the half step value of gas
  gas_rho_h,
  gas_v1_h,
  gas_v2_h,
  gas_v3_h,
  gas_p_h,
  gas_p_S_h,

  // gas_p to store the beginning value of gas for PDMU limiters
  gas_rho_p,
  gas_v1_p,
  gas_v2_p,
  gas_v3_p,
  gas_p_p,
  gas_p_S_p,

  v1_interp_l,
  v2_interp_l,
  v3_interp_l,
  v1_interp_r,
  v2_interp_r,
  v3_interp_r,
  v1_interp,
  v2_interp,
  v3_interp,
  v1_avg,
  v2_avg,
  v3_avg,

  // gasl
  gasl_rho,
  gasl_v1,
  gasl_v2,
  gasl_v3,
  gasl_p,
  gasl_p_S,
  // gasr
  gasr_rho,
  gasr_v1,
  gasr_v2,
  gasr_v3,
  gasr_p,
  gasr_p_S,
  // gasc
  gasc_rho,
  gasc_rhov1,
  gasc_rhov2,
  gasc_rhov3,
  gasc_eng,
  gasc_entropy,
  // gascf
  gascf_rho,
  gascf_rhov1,
  gascf_rhov2,
  gascf_rhov3,
  gascf_eng,
  gascf_entropy,
  // gasc0
  gasc0_rho,
  gasc0_rhov1,
  gasc0_rhov2,
  gasc0_rhov3,
  gasc0_eng,
  gasc0_entropy,
  // dgascf
  dgascf_rho,
  dgascf_rhov1,
  dgascf_rhov2,
  dgascf_rhov3,
  dgascf_eng,
  dgascf_entropy,
  // v1_tmp
  v1_tmp,
  v2_tmp,
  v3_tmp,
  // v1_0
  v1_0,
  v2_0,
  v3_0,
  // rhov1_gas
  rhov1_fluid,
  rhov2_fluid,
  rhov3_fluid,
  rhov1_h,
  rhov2_h,
  rhov3_h,
  // source term due to hydro
  S1_HD,
  S2_HD,
  S3_HD,
  // bdotv
  bdotv,
  NF_gas
};
extern int gas_onface_i[NF_gas];
extern int gas_onface_j[NF_gas];
extern int gas_onface_k[NF_gas];
extern int NF_gas_prim;
extern int NF_gas_store;

typedef enum { LOW_BOUNDARY = 0, HIGH_BOUNDARY = 1 } boundary_t;
/**
 * Reads the configuration from a file.
 *
 * @param filename The name of the file to read from.
 * @return 0 on success, -1 on failure.
 */
int read_config(const char *filename);

/**
 * Initialize the configuration of the solver.
 *
 * @param rank The rank of the process.
 * @param size The size of the communicator.
 * @param filename The name of the configuration file to read.
 * @return 0 on success, -1 on failure.
 */
int initialize_config(const int rank, const int size, const char *filename);

/**
 * Initializes grid indices for the solver.
 * This function sets up the grid indices (is, ie, isg, ieg, etc.) based on the global configuration.
 */
void initialize_grid_indices(void);

/**
 * Initializes logging for the solver.
 *
 * @param base_filename The base name of the log file.
 * @param rank The rank of the process.
 * @return 0 on success, -1 on failure.
 */
int initialize_logging(const char *base_filename, int rank);

#endif
