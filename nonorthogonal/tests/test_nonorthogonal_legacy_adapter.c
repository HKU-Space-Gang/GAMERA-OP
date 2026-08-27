#include "config.h"
#include "nonorthogonal_legacy_adapter.h"
#include "nonorthogonal_state.h"
#include "utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

config_t config;
double ***x1;
double ***x2;
double ***x3;
double ***x1ctr;
double ***x2ctr;
double ***x3ctr;
double ****geo;
double ****gem;
double *****gas;
const double gamma_val = 5.0 / 3.0;
double rho_floor = 1.0e-12;
double p_floor = 1.0e-12;
double divB_max;
int is;
int ie;
int isg;
int ieg;
int js;
int je;
int jsg;
int jeg;
int ks;
int ke;
int ksg;
int keg;

static int failures;

static void expect_near(const char *name, double actual, double expected,
                        double tolerance) {
  if (!isfinite(actual) || fabs(actual - expected) > tolerance) {
    fprintf(stderr,
            "FAIL %s: actual=%.17g expected=%.17g tolerance=%.3e\n", name,
            actual, expected, tolerance);
    ++failures;
  }
}

static double dot3(const double left[3], const double right[3]) {
  return left[0] * right[0] + left[1] * right[1] +
         left[2] * right[2];
}

static double norm3(const double vector[3]) {
  return sqrt(dot3(vector, vector));
}

static int allocate_legacy_arrays(void) {
  const size_t ni = (size_t)config.NI;
  const size_t nj = (size_t)config.NJ;
  const size_t nk = (size_t)config.NK;
  x1 = (double ***)alloc_3d_array(ni, nj, nk, sizeof(double));
  x2 = (double ***)alloc_3d_array(ni, nj, nk, sizeof(double));
  x3 = (double ***)alloc_3d_array(ni, nj, nk, sizeof(double));
  x1ctr = (double ***)alloc_3d_array(ni, nj, nk, sizeof(double));
  x2ctr = (double ***)alloc_3d_array(ni, nj, nk, sizeof(double));
  x3ctr = (double ***)alloc_3d_array(ni, nj, nk, sizeof(double));
  geo = (double ****)alloc_4d_array_contiguous(
      (size_t)NF_geo, ni, nj, nk, sizeof(double));
  gem = (double ****)alloc_4d_array_contiguous(
      (size_t)NF_gem, ni, nj, nk, sizeof(double));
  gas = (double *****)alloc_5d_array_with_4d_contiguous(
      (size_t)NS1, (size_t)NF_gas, ni, nj, nk, sizeof(double));
  return x1 == NULL || x2 == NULL || x3 == NULL || x1ctr == NULL ||
                 x2ctr == NULL || x3ctr == NULL || geo == NULL ||
                 gem == NULL || gas == NULL
             ? -1
             : 0;
}

static void free_legacy_arrays(void) {
  if (x1 != NULL) {
    free_3d_array((void ***)x1);
  }
  if (x2 != NULL) {
    free_3d_array((void ***)x2);
  }
  if (x3 != NULL) {
    free_3d_array((void ***)x3);
  }
  if (x1ctr != NULL) {
    free_3d_array((void ***)x1ctr);
  }
  if (x2ctr != NULL) {
    free_3d_array((void ***)x2ctr);
  }
  if (x3ctr != NULL) {
    free_3d_array((void ***)x3ctr);
  }
  if (geo != NULL) {
    free_4d_array_contiguous((void ****)geo, (size_t)NF_geo);
  }
  if (gem != NULL) {
    free_4d_array_contiguous((void ****)gem, (size_t)NF_gem);
  }
  if (gas != NULL) {
    free_5d_array_with_4d_contiguous((void *****)gas, (size_t)NS1,
                                     (size_t)NF_gas);
  }
}

static void initialize_uniform_legacy_state(void) {
  const double magnetic[3] = {0.38, -0.21, 0.16};
  const double area_vector[3][3] = {{1.0, -0.12, 0.0084},
                                     {0.0035, 1.0, -0.07},
                                     {-0.05, 0.006, 1.0}};
  const int face_slot[3] = {mag_bi, mag_bj, mag_bk};
  for (int i = 0; i < config.NI; ++i) {
    for (int j = 0; j < config.NJ; ++j) {
      for (int k = 0; k < config.NK; ++k) {
        x1[i][j][k] = (double)i + 0.12 * (double)j;
        x2[i][j][k] = (double)j + 0.07 * (double)k;
        x3[i][j][k] = (double)k + 0.05 * (double)i;
        gas[0][gas_rho][i][j][k] = 1.2;
        gas[0][gas_v1][i][j][k] = 0.23;
        gas[0][gas_v2][i][j][k] = -0.17;
        gas[0][gas_v3][i][j][k] = 0.09;
        gas[0][gas_p][i][j][k] = 0.74;
        gas[0][gas_p_S][i][j][k] = 0.74;
        for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
          gem[face_slot[direction]][i][j][k] =
              dot3(magnetic, area_vector[direction]) /
              norm3(area_vector[direction]);
        }
      }
    }
  }
}

int main(void) {
  config.NI = 3;
  config.NJ = 3;
  config.NK = 3;
  isg = jsg = ksg = 0;
  ieg = jeg = keg = 1;
  is = js = ks = 0;
  ie = je = ke = 1;
  if (allocate_legacy_arrays() != 0) {
    fprintf(stderr, "FAIL allocating legacy adapter arrays\n");
    free_legacy_arrays();
    return EXIT_FAILURE;
  }
  initialize_uniform_legacy_state();

  if (gamera_no_legacy_adapter_create(false) != 0) {
    fprintf(stderr, "FAIL creating legacy adapter\n");
    free_legacy_arrays();
    return EXIT_FAILURE;
  }
  gamera_no_storage *storage = gamera_no_legacy_storage();
  gamera_no_grid *grid = gamera_no_legacy_grid();
  if (storage == NULL || grid == NULL) {
    fprintf(stderr, "FAIL legacy adapter accessors\n");
    ++failures;
  } else {
    const size_t cell = gamera_no_index3(grid->cell_extent, 0, 0, 0);
    expect_near("adapter Bx", storage->cell_magnetic[cell].value[0], 0.38,
                8.0e-13);
    expect_near("adapter By", storage->cell_magnetic[cell].value[1], -0.21,
                8.0e-13);
    expect_near("adapter Bz", storage->cell_magnetic[cell].value[2], 0.16,
                8.0e-13);
    expect_near("adapter volume", geo[vol_center][0][0][0],
                grid->cell[cell].volume, 2.0e-14);

    gas[0][gas_rho][0][0][0] = 1.35;
    gas[0][gas_v1][0][0][0] = -0.11;
    gas[0][gas_p][0][0][0] = 0.66;
    if (gamera_no_legacy_import_current() != 0) {
      fprintf(stderr, "FAIL importing modified legacy state\n");
      ++failures;
    } else {
      gamera_no_primitive imported;
      (void)gamera_no_conserved_to_primitive(
          &storage->conserved[cell * GAMERA_NO_FLUX_COUNT], gamma_val,
          rho_floor, p_floor, &imported);
      expect_near("adapter imported density", imported.density, 1.35,
                  2.0e-14);
      expect_near("adapter imported velocity", imported.velocity.value[0],
                  -0.11, 2.0e-14);
      expect_near("adapter imported pressure", imported.pressure, 0.66,
                  2.0e-14);
    }
    if (gamera_no_legacy_export() != 0) {
      fprintf(stderr, "FAIL exporting modified adapter state\n");
      ++failures;
    }
    expect_near("adapter exported density", gas[0][gas_rho][0][0][0], 1.35,
                2.0e-14);
  }

  gamera_no_legacy_adapter_destroy();
  free_legacy_arrays();
  if (failures != 0) {
    fprintf(stderr, "%d legacy adapter checks failed\n", failures);
    return EXIT_FAILURE;
  }
  printf("All non-orthogonal legacy adapter checks passed.\n");
  return EXIT_SUCCESS;
}
