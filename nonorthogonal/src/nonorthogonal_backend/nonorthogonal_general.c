#include "curvilinear.h"

#include <float.h>

double H1(double x1_value, double x2_value, double x3_value) {
  (void)x1_value;
  (void)x2_value;
  (void)x3_value;
  return 1.0;
}

double H2(double x1_value, double x2_value, double x3_value) {
  (void)x1_value;
  (void)x2_value;
  (void)x3_value;
  return 1.0;
}

double H3(double x1_value, double x2_value, double x3_value) {
  (void)x1_value;
  (void)x2_value;
  (void)x3_value;
  return 1.0;
}

/* Exact metric arrays are exported when the non-orthogonal adapter is built. */
int set_geometry_arrays(void) { return 0; }
void getweights(void) {}
void Source_Term(dir_t direction) { (void)direction; }
double get_dt_ring(void) { return DBL_MAX; }

void cartesian_to_curvilinear_coord(double x, double y, double z, double *x1,
                                    double *x2, double *x3) {
  *x1 = x;
  *x2 = y;
  *x3 = z;
}

void curvilinear_to_cartesian_coord(double x1, double x2, double x3, double *x,
                                    double *y, double *z) {
  *x = x1;
  *y = x2;
  *z = x3;
}

void cartesian_to_curvilinear_vector(double x, double y, double z, double vx,
                                     double vy, double vz, double *v1,
                                     double *v2, double *v3) {
  (void)x;
  (void)y;
  (void)z;
  *v1 = vx;
  *v2 = vy;
  *v3 = vz;
}

void curvilinear_to_cartesian_vector(double x1, double x2, double x3,
                                     double v1, double v2, double v3,
                                     double *vx, double *vy, double *vz) {
  (void)x1;
  (void)x2;
  (void)x3;
  *vx = v1;
  *vy = v2;
  *vz = v3;
}
