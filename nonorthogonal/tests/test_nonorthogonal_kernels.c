#include "nonorthogonal_advance.h"
#include "nonorthogonal_geometry.h"
#include "nonorthogonal_flux.h"
#include "nonorthogonal_grid.h"
#include "nonorthogonal_operators.h"
#include "nonorthogonal_reconstruction.h"
#include "nonorthogonal_state.h"
#include "nonorthogonal_step.h"
#include "nonorthogonal_storage.h"
#include "nonorthogonal_sweep.h"
#include "nonorthogonal_yinyang.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0;

static void expect_near(const char *name, double actual, double expected,
                        double tolerance) {
  if (!isfinite(actual) || fabs(actual - expected) > tolerance) {
    fprintf(stderr, "FAIL %s: actual=%.17g expected=%.17g tolerance=%.3e\n",
            name, actual, expected, tolerance);
    ++failures;
  }
}

static double dot(gamera_no_vec3 a, gamera_no_vec3 b) {
  return a.value[0] * b.value[0] + a.value[1] * b.value[1] +
         a.value[2] * b.value[2];
}

static gamera_no_vec3 cross(gamera_no_vec3 a, gamera_no_vec3 b) {
  gamera_no_vec3 result = {
      {a.value[1] * b.value[2] - a.value[2] * b.value[1],
       a.value[2] * b.value[0] - a.value[0] * b.value[2],
       a.value[0] * b.value[1] - a.value[1] * b.value[0]}};
  return result;
}

static double norm(gamera_no_vec3 vector) { return sqrt(dot(vector, vector)); }

static gamera_no_vec3 affine_point(gamera_no_vec3 origin, gamera_no_vec3 a,
                                   gamera_no_vec3 b, gamera_no_vec3 c, int i,
                                   int j, int k) {
  gamera_no_vec3 point = origin;
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    point.value[d] += i * a.value[d] + j * b.value[d] + k * c.value[d];
  }
  return point;
}

static void fill_affine_cell(gamera_no_vec3 corners[8], gamera_no_vec3 origin,
                             gamera_no_vec3 a, gamera_no_vec3 b,
                             gamera_no_vec3 c) {
  const int logical_corner[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0},
                                    {0, 1, 0}, {0, 0, 1}, {1, 0, 1},
                                    {1, 1, 1}, {0, 1, 1}};
  for (int n = 0; n < 8; ++n) {
    corners[n] = affine_point(origin, a, b, c, logical_corner[n][0],
                              logical_corner[n][1], logical_corner[n][2]);
  }
}

static void test_unit_cube(void) {
  const gamera_no_vec3 origin = {{0.0, 0.0, 0.0}};
  const gamera_no_vec3 a = {{1.0, 0.0, 0.0}};
  const gamera_no_vec3 b = {{0.0, 1.0, 0.0}};
  const gamera_no_vec3 c = {{0.0, 0.0, 1.0}};
  gamera_no_vec3 corners[8];
  gamera_no_cell_geometry geometry;
  fill_affine_cell(corners, origin, a, b, c);

  if (gamera_no_compute_cell_geometry(corners, &geometry) != 0) {
    fprintf(stderr, "FAIL unit cube geometry returned an error\n");
    ++failures;
    return;
  }
  expect_near("cube volume", geometry.volume, 1.0, 2.0e-14);
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    expect_near("cube centroid", geometry.centroid.value[d], 0.5, 2.0e-14);
    expect_near("cube lower area", geometry.face[d][0].area, 1.0, 2.0e-14);
    expect_near("cube upper area", geometry.face[d][1].area, 1.0, 2.0e-14);
    expect_near("cube cfl length", geometry.cfl_length[d], 1.0, 3.0e-14);
    expect_near("cube lower normal", geometry.face[d][0].normal.value[d], 1.0,
                2.0e-14);
    expect_near("cube upper normal", geometry.face[d][1].normal.value[d], 1.0,
                2.0e-14);
  }
}

static void test_affine_skew_cell(void) {
  const gamera_no_vec3 origin = {{0.4, -0.7, 1.2}};
  const gamera_no_vec3 a = {{2.0, 0.2, 0.1}};
  const gamera_no_vec3 b = {{0.3, 1.5, 0.2}};
  const gamera_no_vec3 c = {{0.1, 0.4, 1.2}};
  gamera_no_vec3 corners[8];
  gamera_no_cell_geometry geometry;
  fill_affine_cell(corners, origin, a, b, c);

  if (gamera_no_compute_cell_geometry(corners, &geometry) != 0) {
    fprintf(stderr, "FAIL affine geometry returned an error\n");
    ++failures;
    return;
  }

  const double exact_volume = fabs(dot(a, cross(b, c)));
  const gamera_no_vec3 area_vector[3] = {cross(b, c), cross(c, a),
                                         cross(a, b)};
  expect_near("affine volume", geometry.volume, exact_volume, 5.0e-14);
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    const double expected_centroid =
        origin.value[d] + 0.5 * (a.value[d] + b.value[d] + c.value[d]);
    expect_near("affine centroid", geometry.centroid.value[d],
                expected_centroid, 5.0e-14);
  }
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    expect_near("affine face area", geometry.face[direction][0].area,
                norm(area_vector[direction]), 5.0e-14);
    for (int side = 0; side < 2; ++side) {
      for (int component = 0; component < GAMERA_NO_DIM; ++component) {
        expect_near("affine oriented area vector",
                    geometry.face[direction][side]
                        .area_vector.value[component],
                    area_vector[direction].value[component], 5.0e-14);
      }
    }
    expect_near("affine triad orthogonality",
                dot(geometry.face[direction][0].normal,
                    geometry.face[direction][0].tangent1),
                0.0, 2.0e-14);
  }

  const gamera_no_vec3 expected_field = {{0.8, -0.35, 1.1}};
  double face_flux[GAMERA_NO_DIM][2];
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const double flux = dot(expected_field, area_vector[direction]);
    face_flux[direction][0] = flux;
    face_flux[direction][1] = flux;
  }
  gamera_no_vec3 recovered;
  double net_flux = 1.0;
  if (gamera_no_flux_to_cell_field(&geometry, face_flux, &recovered,
                                   &net_flux) != 0) {
    fprintf(stderr, "FAIL affine flux-to-field returned an error\n");
    ++failures;
  } else {
    expect_near("affine net flux", net_flux, 0.0, 2.0e-14);
    for (int d = 0; d < GAMERA_NO_DIM; ++d) {
      expect_near("affine recovered field", recovered.value[d],
                  expected_field.value[d], 8.0e-14);
    }
  }
}

static void test_warped_cell(void) {
  const gamera_no_vec3 origin = {{-0.3, 0.2, 0.4}};
  const gamera_no_vec3 a = {{1.2, 0.1, 0.0}};
  const gamera_no_vec3 b = {{0.2, 0.9, 0.1}};
  const gamera_no_vec3 c = {{0.0, 0.3, 1.1}};
  gamera_no_vec3 corners[8];
  gamera_no_cell_geometry geometry;
  fill_affine_cell(corners, origin, a, b, c);
  corners[6].value[0] += 0.17;
  corners[6].value[1] -= 0.08;
  corners[6].value[2] += 0.13;

  if (gamera_no_compute_cell_geometry(corners, &geometry) != 0) {
    fprintf(stderr, "FAIL warped geometry returned an error\n");
    ++failures;
    return;
  }
  if (!(geometry.volume > 0.0) || !isfinite(geometry.volume)) {
    fprintf(stderr, "FAIL warped cell has invalid volume\n");
    ++failures;
  }
  /* Gold values emitted by the unmodified Fortran quadrature.F90 routines. */
  const double gold_centroid[3] = {0.423300568523973020,
                                   0.850640454819178382,
                                   1.01858480680142360};
  const double gold_area[3][2] = {
      {0.986711710683520926, 1.01362568054571667},
      {1.37262522197430092, 1.52933258246423787},
      {1.06681769764098000, 1.10775095336575724}};
  const double gold_face_centroid[3][2][3] = {
      {{-0.2, 0.8, 1.0},
       {1.04488343642762294, 0.887236279866936317, 1.02826354371433948}},
      {{0.3, 0.4, 0.95},
       {0.552354468297568291, 1.28309623908319370, 1.09494530078500385}},
      {{0.4, 0.7, 0.45},
       {0.437995002306707237, 0.990589094289441530, 1.58418667298518701}}};
  expect_near("warped Fortran volume", geometry.volume, 1.17164999999999897,
              8.0e-14);
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    expect_near("warped Fortran centroid", geometry.centroid.value[d],
                gold_centroid[d], 8.0e-14);
  }
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    for (int side = 0; side < 2; ++side) {
      const gamera_no_face_geometry *face = &geometry.face[direction][side];
      expect_near("warped Fortran face area", face->area,
                  gold_area[direction][side], 8.0e-14);
      for (int d = 0; d < GAMERA_NO_DIM; ++d) {
        expect_near("warped Fortran face centroid", face->centroid.value[d],
                    gold_face_centroid[direction][side][d], 8.0e-14);
      }
      expect_near("warped normal norm", norm(face->normal), 1.0, 2.0e-14);
      expect_near("warped tangent1 norm", norm(face->tangent1), 1.0,
                  2.0e-14);
      expect_near("warped tangent2 norm", norm(face->tangent2), 1.0,
                  2.0e-14);
      expect_near("warped normal-tangent1", dot(face->normal, face->tangent1),
                  0.0, 2.0e-14);
      expect_near("warped normal-tangent2", dot(face->normal, face->tangent2),
                  0.0, 2.0e-14);
    }
  }
}

static void test_edge_field_solve(void) {
  const double normal1[2] = {0.8, 0.6};
  const double normal2[2] = {-0.3, sqrt(0.91)};
  const double expected[2] = {1.25, -0.45};
  const double face1 = normal1[0] * expected[0] + normal1[1] * expected[1];
  const double face2 = normal2[0] * expected[0] + normal2[1] * expected[1];
  double recovered[2];
  if (gamera_no_solve_edge_field(normal1, normal2, face1, face2, recovered) !=
      0) {
    fprintf(stderr, "FAIL edge field solve returned an error\n");
    ++failures;
  } else {
    expect_near("edge field b1", recovered[0], expected[0], 2.0e-14);
    expect_near("edge field b2", recovered[1], expected[1], 2.0e-14);
  }

  const double singular[2] = {1.6, 1.2};
  if (gamera_no_solve_edge_field(normal1, singular, face1, face2, recovered) ==
      0) {
    fprintf(stderr, "FAIL singular edge transform was accepted\n");
    ++failures;
  }
}

static void test_edge_geometry_and_normal_interpolation(void) {
  const gamera_no_vec3 start = {{0.0, 0.0, 0.0}};
  const gamera_no_vec3 end = {{0.0, 0.0, 2.5}};
  const gamera_no_vec3 transverse_average = {{1.0, 0.0, 0.0}};
  gamera_no_edge_geometry edge;
  if (gamera_no_compute_edge_geometry(start, end, transverse_average, &edge) !=
      0) {
    fprintf(stderr, "FAIL edge geometry returned an error\n");
    ++failures;
    return;
  }
  expect_near("edge length", edge.length, 2.5, 2.0e-14);
  expect_near("edge normal norm", norm(edge.normal), 1.0, 2.0e-14);
  expect_near("edge tangent1 norm", norm(edge.tangent1), 1.0, 2.0e-14);
  expect_near("edge tangent2 norm", norm(edge.tangent2), 1.0, 2.0e-14);
  expect_near("edge N dot T1", dot(edge.normal, edge.tangent1), 0.0,
              2.0e-14);
  expect_near("edge N dot T2", dot(edge.normal, edge.tangent2), 0.0,
              2.0e-14);

  double area[8] = {0.7, 0.8, 0.95, 1.1, 1.25, 1.4, 1.6, 1.85};
  gamera_no_vec3 normal_x[8];
  gamera_no_vec3 normal_y[8];
  for (int n = 0; n < 8; ++n) {
    normal_x[n] = (gamera_no_vec3){{1.0, 0.0, 0.0}};
    normal_y[n] = (gamera_no_vec3){{0.0, 1.0, 0.0}};
  }
  double row_x[2];
  double row_y[2];
  if (gamera_no_interpolate_face_normal_to_edge(area, normal_x, &edge, row_x) !=
          0 ||
      gamera_no_interpolate_face_normal_to_edge(area, normal_y, &edge, row_y) !=
          0) {
    fprintf(stderr, "FAIL edge normal interpolation returned an error\n");
    ++failures;
    return;
  }
  /* For this frame T1=-y and T2=+x. */
  expect_near("edge projected x/T1", row_x[0], 0.0, 2.0e-14);
  expect_near("edge projected x/T2", row_x[1], 1.0, 2.0e-14);
  expect_near("edge projected y/T1", row_y[0], -1.0, 2.0e-14);
  expect_near("edge projected y/T2", row_y[1], 0.0, 2.0e-14);

  const double expected[2] = {-0.62, 1.18};
  const double bx = row_x[0] * expected[0] + row_x[1] * expected[1];
  const double by = row_y[0] * expected[0] + row_y[1] * expected[1];
  double recovered[2];
  if (gamera_no_solve_edge_field(row_x, row_y, bx, by, recovered) != 0) {
    fprintf(stderr, "FAIL interpolated edge transform solve returned an error\n");
    ++failures;
  } else {
    expect_near("interpolated edge field b1", recovered[0], expected[0],
                2.0e-14);
    expect_near("interpolated edge field b2", recovered[1], expected[1],
                2.0e-14);
  }
}

static void test_ct_divergence_identity(void) {
  const double emf_i[2][2] = {{0.13, -0.27}, {0.41, 0.09}};
  const double emf_j[2][2] = {{-0.22, 0.37}, {0.16, -0.31}};
  const double emf_k[2][2] = {{0.45, -0.11}, {-0.29, 0.34}};
  double increment[GAMERA_NO_DIM][2];
  gamera_no_ct_face_increments(emf_i, emf_j, emf_k, 0.07, increment);

  const double divergence_increment =
      increment[GAMERA_NO_I][GAMERA_NO_UPPER] -
      increment[GAMERA_NO_I][GAMERA_NO_LOWER] +
      increment[GAMERA_NO_J][GAMERA_NO_UPPER] -
      increment[GAMERA_NO_J][GAMERA_NO_LOWER] +
      increment[GAMERA_NO_K][GAMERA_NO_UPPER] -
      increment[GAMERA_NO_K][GAMERA_NO_LOWER];
  expect_near("CT divergence-of-curl", divergence_increment, 0.0, 2.0e-16);
}

static void test_edge_emf(void) {
  gamera_no_edge_geometry edge = {0};
  edge.length = 1.7;
  double emf;
  double diffusion;
  if (gamera_no_compute_edge_emf(0.33, -0.41, 0.52, 0.27, -0.18, 0.93,
                                 &edge, 0.5, true, 1.4, 0.3, 0.02, &emf,
                                 &diffusion) != 0) {
    fprintf(stderr, "FAIL edge EMF returned an error\n");
    ++failures;
    return;
  }
  expect_near("edge EMF diffusion", diffusion, 1.0836549182839057,
              3.0e-14);
  expect_near("edge EMF", emf, -0.6797092024974376, 3.0e-14);

  if (gamera_no_compute_edge_emf(0.33, -0.41, 0.52, 0.27, -0.18, 0.93,
                                 &edge, 0.5, true, 1.4, 0.3, 1.0, &emf,
                                 &diffusion) != 0) {
    fprintf(stderr, "FAIL capped edge EMF returned an error\n");
    ++failures;
    return;
  }
  expect_near("edge EMF CFL cap", diffusion, 0.51, 2.0e-14);
  expect_near("capped edge EMF", emf, -0.59194, 2.0e-14);
}

static void test_reconstruction(void) {
  const double constant[GAMERA_NO_RECON_STENCIL] =
      {2.75, 2.75, 2.75, 2.75, 2.75, 2.75, 2.75, 2.75};
  const double varying_volume[GAMERA_NO_RECON_STENCIL] =
      {0.73, 0.81, 0.94, 1.02, 1.17, 1.31, 1.46, 1.64};
  double left;
  double right;
  if (gamera_no_reconstruct_up7_pdm(varying_volume, constant, 1.0, &left,
                                    &right) != 0) {
    fprintf(stderr, "FAIL constant reconstruction returned an error\n");
    ++failures;
  } else {
    expect_near("constant reconstruction left", left, 2.75, 2.0e-14);
    expect_near("constant reconstruction right", right, 2.75, 2.0e-14);
  }

  double linear[GAMERA_NO_RECON_STENCIL];
  double unit_volume[GAMERA_NO_RECON_STENCIL];
  for (int n = 0; n < GAMERA_NO_RECON_STENCIL; ++n) {
    const double x = (double)n - 3.5;
    linear[n] = 1.4 - 0.3 * x;
    unit_volume[n] = 1.0;
  }
  if (gamera_no_reconstruct_up7_pdm(unit_volume, linear, 1.0, &left,
                                    &right) != 0) {
    fprintf(stderr, "FAIL linear reconstruction returned an error\n");
    ++failures;
  } else {
    expect_near("linear reconstruction left", left, 1.4, 2.0e-14);
    expect_near("linear reconstruction right", right, 1.4, 2.0e-14);
  }

  const double jump[GAMERA_NO_RECON_STENCIL] =
      {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
  if (gamera_no_reconstruct_up7_pdm(unit_volume, jump, 1.0, &left, &right) !=
      0) {
    fprintf(stderr, "FAIL jump reconstruction returned an error\n");
    ++failures;
  } else {
    if (left < 0.0 || left > 1.0 || right < 0.0 || right > 1.0) {
      fprintf(stderr, "FAIL PDM jump values escaped neighbor bounds\n");
      ++failures;
    }
  }

  expect_near("central8 constant", gamera_no_central8(constant), 2.75,
              2.0e-14);
  expect_near("central6 constant", gamera_no_central6(constant), 2.75,
              2.0e-14);

  const double gold_primitive[GAMERA_NO_RECON_STENCIL] =
      {0.31, 0.48, 0.72, 1.05, 1.29, 1.37, 1.21, 0.86};
  if (gamera_no_reconstruct_up7_pdm(varying_volume, gold_primitive, 1.0,
                                    &left, &right) != 0) {
    fprintf(stderr, "FAIL Fortran-gold reconstruction returned an error\n");
    ++failures;
  } else {
    expect_near("Fortran-gold reconstruction left", left,
                1.19178156891945175, 3.0e-14);
    expect_near("Fortran-gold reconstruction right", right,
                1.19391748061845759, 3.0e-14);
  }
}

static void test_fortran_flux_gold(void) {
  const gamera_no_primitive state[2] = {
      {1.2, {{0.7, -0.2, 0.15}}, 0.9},
      {0.85, {{-0.35, 0.4, -0.1}}, 0.65}};
  const gamera_no_vec3 magnetic[2] = {{{0.45, -0.3, 0.2}},
                                      {{0.38, -0.18, 0.27}}};
  gamera_no_face_geometry face = {0};
  face.normal = (gamera_no_vec3){{0.8, 0.36, 0.48}};
  face.tangent1 = (gamera_no_vec3){{0.6, -0.48, -0.64}};
  face.tangent2 = (gamera_no_vec3){{0.0, 0.8, -0.6}};

  gamera_no_fluid_flux fluid;
  if (gamera_no_kinetic_fluid_flux(state, 5.0 / 3.0, &face, &fluid) != 0) {
    fprintf(stderr, "FAIL kinetic fluid flux returned an error\n");
    ++failures;
    return;
  }
  const double gold_fluid[5] = {0.453096698490859873,
                                1.55440363356794276,
                                0.0572190315226311780,
                                0.665421441742305997,
                                0.995030259167600684};
  const double gold_jump[5] = {-0.35, -1.1375, 0.58, -0.265,
                               -0.5821875};
  const double gold_velocity_jump[3] = {-1.05, 0.6, -0.25};
  const double gold_normal_velocity[2] = {0.56, -0.184};
  for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
    expect_near("Fortran-gold fluid flux", fluid.conserved[variable],
                gold_fluid[variable], 5.0e-14);
    expect_near("Fortran-gold conserved jump",
                fluid.conserved_jump[variable], gold_jump[variable],
                5.0e-14);
  }
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    expect_near("Fortran-gold velocity jump", fluid.velocity_jump.value[d],
                gold_velocity_jump[d], 5.0e-14);
  }
  for (int side = 0; side < 2; ++side) {
    expect_near("Fortran-gold normal velocity", fluid.normal_velocity[side],
                gold_normal_velocity[side], 5.0e-14);
  }

  gamera_no_maxwell_flux maxwell;
  if (gamera_no_kinetic_maxwell_flux(
          state, magnetic, 0.12, &face, true,
          (gamera_no_vec3){{0.2, 0.05, -0.1}}, 0.08, true, 1.7,
          fluid.normal_velocity, &maxwell) != 0) {
    fprintf(stderr, "FAIL kinetic Maxwell flux returned an error\n");
    ++failures;
    return;
  }
  const double gold_maxwell[3] = {0.0636265312070995737,
                                  0.147454531285714979,
                                  0.0775737756422462954};
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    expect_near("Fortran-gold Maxwell flux", maxwell.momentum.value[d],
                gold_maxwell[d], 5.0e-14);
  }
  expect_near("Fortran-gold Alfven diffusion",
              maxwell.alfven_diffusion_speed, 0.612298816333177509,
              5.0e-14);
  expect_near("Fortran-gold magnetic pressure sum",
              maxwell.magnetic_pressure_sum, 0.438600000000000101,
              5.0e-14);

  if (gamera_no_apply_hogs(&fluid, &maxwell, 0.25, 0.25, true, 1.7) != 0) {
    fprintf(stderr, "FAIL HOGS returned an error\n");
    ++failures;
    return;
  }
  const double gold_fluid_hogs[5] = {0.506672844920012855,
                                     1.72852610946269003,
                                     -0.0315642968456795736,
                                     0.705986238324379012,
                                     1.08414843845109354};
  const double gold_maxwell_hogs[3] = {0.0880194355224904329,
                                       0.133515728819777330,
                                       0.0833816100030536494};
  for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
    expect_near("Fortran-gold fluid HOGS", fluid.conserved[variable],
                gold_fluid_hogs[variable], 7.0e-14);
  }
  for (int d = 0; d < GAMERA_NO_DIM; ++d) {
    expect_near("Fortran-gold Maxwell HOGS", maxwell.momentum.value[d],
                gold_maxwell_hogs[d], 7.0e-14);
  }
}

static void test_fortran_state_gold(void) {
  const double gamma = 5.0 / 3.0;
  const double density_floor = 1.0e-12;
  const double pressure_floor = 1.0e-12;
  const gamera_no_primitive old_primitive =
      {1.0, {{0.2, -0.1, 0.3}}, 0.7};
  const gamera_no_primitive current_primitive =
      {1.1, {{0.35, -0.08, 0.25}}, 0.8};
  double old_conserved[GAMERA_NO_FLUX_COUNT];
  double current_conserved[GAMERA_NO_FLUX_COUNT];
  double predicted[GAMERA_NO_FLUX_COUNT];
  if (gamera_no_primitive_to_conserved(
          &old_primitive, gamma, density_floor, pressure_floor,
          old_conserved) != 0 ||
      gamera_no_primitive_to_conserved(
          &current_primitive, gamma, density_floor, pressure_floor,
          current_conserved) != 0 ||
      gamera_no_predict_cell(old_conserved, current_conserved, 0.45, gamma,
                             density_floor, pressure_floor, predicted) != 0) {
    fprintf(stderr, "FAIL state predictor returned an error\n");
    ++failures;
    return;
  }
  const double gold_predictor[5] = {1.145,
                                    0.4780375,
                                    -0.081295,
                                    0.2604875,
                                    1.39980675375000008};
  for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
    expect_near("Fortran-gold predictor", predicted[variable],
                gold_predictor[variable], 5.0e-14);
  }

  if (gamera_no_apply_maxwell(predicted,
                              (gamera_no_vec3){{0.12, -0.07, 0.05}}, 0.03,
                              gamma, density_floor, pressure_floor) != 0) {
    fprintf(stderr, "FAIL Maxwell update returned an error\n");
    ++failures;
    return;
  }
  const double gold_maxwell[5] = {1.145,
                                  0.4816375,
                                  -0.083395,
                                  0.2619875,
                                  1.40180867143558952};
  for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
    expect_near("Fortran-gold Maxwell update", predicted[variable],
                gold_maxwell[variable], 5.0e-14);
  }

  const gamera_no_primitive boris_updated_primitive =
      {1.05, {{0.4, 0.1, -0.15}}, 0.75};
  const gamera_no_primitive boris_old_primitive =
      {1.0, {{0.3, 0.05, -0.1}}, 0.72};
  double boris_updated[5];
  double boris_old[5];
  const double hydro_rate[5] = {0.08, 0.12, -0.04, 0.03, 0.2};
  if (gamera_no_primitive_to_conserved(
          &boris_updated_primitive, gamma, density_floor, pressure_floor,
          boris_updated) != 0 ||
      gamera_no_primitive_to_conserved(
          &boris_old_primitive, gamma, density_floor, pressure_floor,
          boris_old) != 0 ||
      gamera_no_apply_boris(
          boris_updated, boris_old, (gamera_no_vec3){{0.6, -0.2, 0.3}},
          (gamera_no_vec3){{0.55, -0.25, 0.28}}, hydro_rate,
          (gamera_no_vec3){{0.15, -0.11, 0.09}}, 1.7, 0.03, gamma,
          density_floor, pressure_floor) != 0) {
    fprintf(stderr, "FAIL Boris update returned an error\n");
    ++failures;
    return;
  }
  const double gold_boris[5] = {1.05,
                                0.303087506797275918,
                                0.0476917251487173691,
                                -0.0989437181033525764,
                                1.17448876036969763};
  for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
    expect_near("Fortran-gold Boris update", boris_updated[variable],
                gold_boris[variable], 7.0e-14);
  }
}

static void test_nuclear_hogs(void) {
  gamera_no_fluid_flux fluid = {0};
  const gamera_no_primitive interface[2] = {
      {1.0, {{15.0, 0.0, 0.0}}, 1.0},
      {1.0, {{15.2, 0.0, 0.0}}, 1.0}};
  const double lower[GAMERA_NO_FLUX_COUNT] = {1.0, 2.0, 3.0, 4.0, 5.0};
  const double upper[GAMERA_NO_FLUX_COUNT] = {1.5, 3.0, 2.0, 4.5, 7.0};
  bool applied = false;
  if (gamera_no_apply_nuclear_hogs(&fluid, interface, lower, upper, true,
                                   10.0, &applied) != 0 || !applied) {
    fprintf(stderr, "FAIL Kaiju nuclear HOGS did not trigger\n");
    ++failures;
    return;
  }
  for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
    expect_near("Kaiju nuclear HOGS flux", fluid.conserved[variable],
                -10.0 * (upper[variable] - lower[variable]), 1.0e-15);
  }
  gamera_no_fluid_flux quiet = {0};
  gamera_no_primitive slow[2] = {interface[0], interface[1]};
  slow[0].velocity.value[0] = 14.9;
  slow[1].velocity.value[0] = 14.9;
  if (gamera_no_apply_nuclear_hogs(&quiet, slow, lower, upper, true, 10.0,
                                   &applied) != 0 || applied) {
    fprintf(stderr, "FAIL Kaiju nuclear HOGS sub-threshold behavior\n");
    ++failures;
  }
}

static void test_grid_and_storage_layout(void) {
  const size_t cells[3] = {2, 2, 2};
  const size_t vertices_extent[3] = {3, 3, 3};
  gamera_no_vec3 vertices[27];
  for (size_t i = 0; i < vertices_extent[0]; ++i) {
    for (size_t j = 0; j < vertices_extent[1]; ++j) {
      for (size_t k = 0; k < vertices_extent[2]; ++k) {
        gamera_no_vec3 point =
            {{(double)i + 0.1 * (double)j,
              (double)j + 0.07 * (double)k,
              (double)k + 0.04 * (double)i}};
        point.value[0] +=
            0.01 * (double)i * (double)j * (double)k;
        vertices[gamera_no_index3(vertices_extent, i, j, k)] = point;
      }
    }
  }

  gamera_no_grid grid;
  if (gamera_no_grid_create(cells, vertices, &grid) != 0) {
    fprintf(stderr, "FAIL staggered grid creation returned an error\n");
    ++failures;
    return;
  }
  if (gamera_no_element_count3(grid.cell_extent) != 8 ||
      gamera_no_element_count3(grid.face[GAMERA_NO_I].extent) != 12 ||
      gamera_no_element_count3(grid.face[GAMERA_NO_J].extent) != 12 ||
      gamera_no_element_count3(grid.face[GAMERA_NO_K].extent) != 12 ||
      gamera_no_element_count3(grid.edge[GAMERA_NO_I].extent) != 18 ||
      gamera_no_element_count3(grid.edge[GAMERA_NO_J].extent) != 18 ||
      gamera_no_element_count3(grid.edge[GAMERA_NO_K].extent) != 18) {
    fprintf(stderr, "FAIL staggered grid extents are inconsistent\n");
    ++failures;
  }

  const gamera_no_face_geometry *shared_face =
      &grid.face[GAMERA_NO_I]
           .value[gamera_no_index3(grid.face[GAMERA_NO_I].extent, 1, 0, 0)];
  const gamera_no_cell_geometry *left_cell =
      &grid.cell[gamera_no_index3(grid.cell_extent, 0, 0, 0)];
  const gamera_no_cell_geometry *right_cell =
      &grid.cell[gamera_no_index3(grid.cell_extent, 1, 0, 0)];
  expect_near("shared face from left cell", shared_face->area,
              left_cell->face[GAMERA_NO_I][GAMERA_NO_UPPER].area, 5.0e-14);
  expect_near("shared face from right cell", shared_face->area,
              right_cell->face[GAMERA_NO_I][GAMERA_NO_LOWER].area, 5.0e-14);

  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    size_t valid_count = 0;
    const size_t edge_count =
        gamera_no_element_count3(grid.edge[direction].extent);
    for (size_t n = 0; n < edge_count; ++n) {
      valid_count += grid.edge[direction].valid[n] != 0U ? 1U : 0U;
    }
    if (valid_count != 2) {
      fprintf(stderr, "FAIL interior edge count in direction %d: %zu\n",
              direction, valid_count);
      ++failures;
    }
  }
  gamera_no_grid_destroy(&grid);
  if (grid.vertex != NULL || grid.cell != NULL) {
    fprintf(stderr, "FAIL grid destroy did not clear ownership\n");
    ++failures;
  }

  gamera_no_storage storage;
  const size_t state_cells[3] = {3, 2, 4};
  if (gamera_no_storage_create(state_cells, &storage) != 0) {
    fprintf(stderr, "FAIL non-orthogonal storage creation returned an error\n");
    ++failures;
    return;
  }
  if (storage.conserved == NULL || storage.face_flux[GAMERA_NO_J] == NULL ||
      storage.edge_emf[GAMERA_NO_K] == NULL || storage.conserved[0] != 0.0 ||
      storage.face_flux[GAMERA_NO_I][0] != 0.0 ||
      storage.fluid_face_flux[GAMERA_NO_I] == NULL ||
      storage.maxwell_face_flux[GAMERA_NO_J] == NULL ||
      storage.hydro_rate == NULL || storage.maxwell_rate == NULL) {
    fprintf(stderr, "FAIL non-orthogonal storage initialization\n");
    ++failures;
  }
  gamera_no_storage_destroy(&storage);
  if (storage.conserved != NULL || storage.face_flux[0] != NULL) {
    fprintf(stderr, "FAIL storage destroy did not clear ownership\n");
    ++failures;
  }
}

static void test_uniform_face_sweep(void) {
  const size_t extent[3] = {9, 1, 1};
  const size_t vertex_extent[3] = {10, 2, 2};
  gamera_no_vec3 vertices[40];
  for (size_t i = 0; i < vertex_extent[0]; ++i) {
    for (size_t j = 0; j < vertex_extent[1]; ++j) {
      for (size_t k = 0; k < vertex_extent[2]; ++k) {
        vertices[gamera_no_index3(vertex_extent, i, j, k)] =
            (gamera_no_vec3){{(double)i + 0.12 * (double)j,
                              (double)j + 0.07 * (double)k,
                              (double)k + 0.05 * (double)i}};
      }
    }
  }
  gamera_no_grid grid = {0};
  gamera_no_storage storage = {0};
  if (gamera_no_grid_create(extent, vertices, &grid) != 0 ||
      gamera_no_storage_create(extent, &storage) != 0) {
    fprintf(stderr, "FAIL uniform sweep setup returned an error\n");
    ++failures;
    gamera_no_storage_destroy(&storage);
    gamera_no_grid_destroy(&grid);
    return;
  }

  const gamera_no_primitive primitive =
      {1.2, {{0.23, -0.17, 0.09}}, 0.74};
  const gamera_no_vec3 magnetic = {{0.38, -0.21, 0.16}};
  const size_t cell_count = gamera_no_element_count3(extent);
  for (size_t cell = 0; cell < cell_count; ++cell) {
    if (gamera_no_primitive_to_conserved(
            &primitive, 5.0 / 3.0, 1.0e-12, 1.0e-12,
            &storage.predicted_conserved[cell * GAMERA_NO_FLUX_COUNT]) != 0) {
      fprintf(stderr, "FAIL uniform sweep conserved conversion\n");
      ++failures;
      break;
    }
    storage.predicted_cell_magnetic[cell] = magnetic;
  }
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t face_count =
        gamera_no_element_count3(grid.face[direction].extent);
    for (size_t face = 0; face < face_count; ++face) {
      storage.predicted_face_flux[direction][face] =
          dot(magnetic, grid.face[direction].value[face].area_vector);
    }
  }

  const double *face_flux[3] = {
      storage.predicted_face_flux[GAMERA_NO_I],
      storage.predicted_face_flux[GAMERA_NO_J],
      storage.predicted_face_flux[GAMERA_NO_K]};
  const size_t active_lower[3] = {4, 0, 0};
  const size_t active_upper[3] = {5, 1, 1};
  const gamera_no_sweep_options options = {
      5.0 / 3.0, 1.0e-12, 1.0e-12, 1.0, 0.15, 0.2, 2.3,
      true,        true,      true,      false};
  if (gamera_no_sweep_face_fluxes(
          &grid, &storage, storage.predicted_conserved,
          storage.predicted_cell_magnetic, face_flux, GAMERA_NO_I,
          active_lower, active_upper, &options, NULL) != 0) {
    fprintf(stderr, "FAIL uniform non-orthogonal face sweep\n");
    ++failures;
  } else {
    const size_t lower_face =
        gamera_no_index3(grid.face[GAMERA_NO_I].extent, 4, 0, 0);
    const size_t upper_face =
        gamera_no_index3(grid.face[GAMERA_NO_I].extent, 5, 0, 0);
    for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
      expect_near(
          "uniform sweep fluid conservation",
          storage.fluid_face_flux[GAMERA_NO_I]
                                 [lower_face * GAMERA_NO_FLUX_COUNT +
                                  (size_t)variable],
          storage.fluid_face_flux[GAMERA_NO_I]
                                 [upper_face * GAMERA_NO_FLUX_COUNT +
                                  (size_t)variable],
          2.0e-13);
    }
    for (int component = 0; component < GAMERA_NO_DIM; ++component) {
      expect_near(
          "uniform sweep Maxwell conservation",
          storage.maxwell_face_flux[GAMERA_NO_I][lower_face].value[component],
          storage.maxwell_face_flux[GAMERA_NO_I][upper_face].value[component],
          2.0e-13);
    }
  }

  double local_dt = 0.0;
  const gamera_no_timestep_options timestep_options = {
      5.0 / 3.0, 1.0e-12, 1.0e-12, 0.3, 2.3, true, true, false};
  if (gamera_no_local_timestep(
          &grid, storage.predicted_conserved,
          storage.predicted_cell_magnetic, active_lower, active_upper,
          &timestep_options, NULL, &local_dt) != 0) {
    fprintf(stderr, "FAIL non-orthogonal local timestep\n");
    ++failures;
  } else {
    const size_t cell = gamera_no_index3(extent, 4, 0, 0);
    const double sound_speed =
        sqrt((5.0 / 3.0) * primitive.pressure / primitive.density);
    const double raw_alfven = norm(magnetic) / sqrt(primitive.density);
    const double corrected_alfven =
        2.3 * raw_alfven / hypot(2.3, raw_alfven);
    double length = grid.cell[cell].cfl_length[0];
    length = fmin(length, grid.cell[cell].cfl_length[1]);
    length = fmin(length, grid.cell[cell].cfl_length[2]);
    const double expected_dt =
        0.3 * length /
        (norm(primitive.velocity) + hypot(sound_speed, corrected_alfven));
    expect_near("non-orthogonal local timestep", local_dt, expected_dt,
                3.0e-15);
  }

  if (gamera_no_calculate_stress_rates(&grid, &storage, active_lower,
                                       active_upper, true, NULL) != 0) {
    fprintf(stderr, "FAIL stress divergence assembly\n");
    ++failures;
  } else {
    const size_t cell = gamera_no_index3(extent, 4, 0, 0);
    for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
      expect_near("uniform sweep hydro rate",
                  storage.hydro_rate[cell * GAMERA_NO_FLUX_COUNT +
                                     (size_t)variable],
                  0.0, 3.0e-13);
    }
    for (int component = 0; component < GAMERA_NO_DIM; ++component) {
      expect_near("uniform sweep Maxwell rate",
                  storage.maxwell_rate[cell].value[component], 0.0, 3.0e-13);
    }
  }

  gamera_no_storage_destroy(&storage);
  gamera_no_grid_destroy(&grid);
}

static void test_uniform_edge_emf_sweep(void) {
  const size_t extent[3] = {1, 9, 9};
  const size_t vertex_extent[3] = {2, 10, 10};
  gamera_no_vec3 vertices[200];
  for (size_t i = 0; i < vertex_extent[0]; ++i) {
    for (size_t j = 0; j < vertex_extent[1]; ++j) {
      for (size_t k = 0; k < vertex_extent[2]; ++k) {
        vertices[gamera_no_index3(vertex_extent, i, j, k)] =
            (gamera_no_vec3){{(double)i + 0.12 * (double)j,
                              (double)j + 0.07 * (double)k,
                              (double)k + 0.05 * (double)i}};
      }
    }
  }
  gamera_no_grid grid = {0};
  gamera_no_storage storage = {0};
  if (gamera_no_grid_create(extent, vertices, &grid) != 0 ||
      gamera_no_storage_create(extent, &storage) != 0) {
    fprintf(stderr, "FAIL uniform edge sweep setup returned an error\n");
    ++failures;
    gamera_no_storage_destroy(&storage);
    gamera_no_grid_destroy(&grid);
    return;
  }

  const gamera_no_primitive primitive =
      {1.2, {{0.23, -0.17, 0.09}}, 0.74};
  const gamera_no_vec3 magnetic = {{0.38, -0.21, 0.16}};
  const size_t cell_count = gamera_no_element_count3(extent);
  for (size_t cell = 0; cell < cell_count; ++cell) {
    (void)gamera_no_primitive_to_conserved(
        &primitive, 5.0 / 3.0, 1.0e-12, 1.0e-12,
        &storage.predicted_conserved[cell * GAMERA_NO_FLUX_COUNT]);
  }
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t face_count =
        gamera_no_element_count3(grid.face[direction].extent);
    for (size_t face = 0; face < face_count; ++face) {
      storage.predicted_face_flux[direction][face] =
          dot(magnetic, grid.face[direction].value[face].area_vector);
    }
  }
  const double *face_flux[3] = {
      storage.predicted_face_flux[GAMERA_NO_I],
      storage.predicted_face_flux[GAMERA_NO_J],
      storage.predicted_face_flux[GAMERA_NO_K]};
  const size_t active_lower[3] = {0, 4, 4};
  const size_t active_upper[3] = {1, 5, 5};
  const gamera_no_emf_options options = {1.0e-12, 1.0, 0.6, 2.3,
                                         0.8,     0.02, true, false};
  if (gamera_no_sweep_edge_emf(
          &grid, &storage, storage.predicted_conserved, face_flux,
          GAMERA_NO_I, active_lower, active_upper, &options, NULL) != 0) {
    fprintf(stderr, "FAIL uniform non-orthogonal edge EMF sweep\n");
    ++failures;
  } else {
    const gamera_no_vec3 velocity_cross_magnetic =
        cross(primitive.velocity, magnetic);
    for (size_t j = active_lower[1]; j <= active_upper[1]; ++j) {
      for (size_t k = active_lower[2]; k <= active_upper[2]; ++k) {
        const size_t edge =
            gamera_no_index3(grid.edge[GAMERA_NO_I].extent, 0, j, k);
        const gamera_no_edge_geometry *geometry =
            &grid.edge[GAMERA_NO_I].value[edge];
        const double expected =
            -dot(velocity_cross_magnetic, geometry->normal) *
            geometry->length;
        expect_near("uniform edge ideal EMF",
                    storage.edge_emf[GAMERA_NO_I][edge], expected, 8.0e-13);
      }
    }
  }
  gamera_no_storage_destroy(&storage);
  gamera_no_grid_destroy(&grid);
}

static void test_grid_level_ct_and_predictor(void) {
  const size_t extent[3] = {2, 2, 2};
  const size_t vertex_extent[3] = {3, 3, 3};
  gamera_no_vec3 vertices[27];
  for (size_t i = 0; i < vertex_extent[0]; ++i) {
    for (size_t j = 0; j < vertex_extent[1]; ++j) {
      for (size_t k = 0; k < vertex_extent[2]; ++k) {
        vertices[gamera_no_index3(vertex_extent, i, j, k)] =
            (gamera_no_vec3){{(double)i + 0.08 * (double)j,
                              (double)j + 0.05 * (double)k,
                              (double)k + 0.03 * (double)i}};
      }
    }
  }
  gamera_no_grid grid = {0};
  gamera_no_storage storage = {0};
  if (gamera_no_grid_create(extent, vertices, &grid) != 0 ||
      gamera_no_storage_create(extent, &storage) != 0) {
    fprintf(stderr, "FAIL grid-level CT setup returned an error\n");
    ++failures;
    gamera_no_grid_destroy(&grid);
    gamera_no_storage_destroy(&storage);
    return;
  }

  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t face_count =
        gamera_no_element_count3(storage.face_extent[direction]);
    const size_t edge_count =
        gamera_no_element_count3(storage.edge_extent[direction]);
    for (size_t n = 0; n < face_count; ++n) {
      storage.face_flux[direction][n] =
          0.01 * (double)(direction + 1) * (double)(n + 1);
      storage.old_face_flux[direction][n] =
          storage.face_flux[direction][n] - 0.001 * (double)(n + 1);
    }
    for (size_t n = 0; n < edge_count; ++n) {
      storage.edge_emf[direction][n] =
          0.003 * (double)(direction + 2) * (double)(n + 1) - 0.02;
    }
  }

  double old_divergence[8];
  for (size_t i = 0; i < extent[0]; ++i) {
    for (size_t j = 0; j < extent[1]; ++j) {
      for (size_t k = 0; k < extent[2]; ++k) {
        const size_t cell = gamera_no_index3(extent, i, j, k);
        old_divergence[cell] = gamera_no_cell_net_flux(&storage, i, j, k);
      }
    }
  }
  if (gamera_no_advance_ct(&storage, 0.07) != 0) {
    fprintf(stderr, "FAIL grid-level CT update returned an error\n");
    ++failures;
  } else {
    for (size_t i = 0; i < extent[0]; ++i) {
      for (size_t j = 0; j < extent[1]; ++j) {
        for (size_t k = 0; k < extent[2]; ++k) {
          const size_t cell = gamera_no_index3(extent, i, j, k);
          expect_near("grid CT divergence preservation",
                      gamera_no_cell_net_flux(&storage, i, j, k),
                      old_divergence[cell], 3.0e-16);
        }
      }
    }
  }

  const double *current_flux[3] = {storage.face_flux[GAMERA_NO_I],
                                   storage.face_flux[GAMERA_NO_J],
                                   storage.face_flux[GAMERA_NO_K]};
  if (gamera_no_recover_magnetic_field(&grid, current_flux,
                                       storage.cell_magnetic) != 0) {
    fprintf(stderr, "FAIL grid magnetic recovery returned an error\n");
    ++failures;
  }

  const gamera_no_primitive old_primitive =
      {1.0, {{0.1, -0.2, 0.3}}, 0.8};
  const gamera_no_primitive current_primitive =
      {1.04, {{0.13, -0.18, 0.27}}, 0.84};
  double old_cell[5];
  double current_cell[5];
  (void)gamera_no_primitive_to_conserved(&old_primitive, 5.0 / 3.0, 1.0e-12,
                                         1.0e-12, old_cell);
  (void)gamera_no_primitive_to_conserved(
      &current_primitive, 5.0 / 3.0, 1.0e-12, 1.0e-12, current_cell);
  const size_t cell_count = gamera_no_element_count3(extent);
  for (size_t cell = 0; cell < cell_count; ++cell) {
    memcpy(&storage.old_conserved[cell * 5], old_cell, sizeof(old_cell));
    memcpy(&storage.conserved[cell * 5], current_cell, sizeof(current_cell));
  }
  const double expected_predicted_face =
      storage.face_flux[GAMERA_NO_I][0] +
      0.5 * (storage.face_flux[GAMERA_NO_I][0] -
             storage.old_face_flux[GAMERA_NO_I][0]);
  if (gamera_no_predict_storage(&storage, &grid, 0.5, 5.0 / 3.0, 1.0e-12,
                                1.0e-12) != 0) {
    fprintf(stderr, "FAIL grid predictor returned an error\n");
    ++failures;
  } else {
    expect_near("grid predicted face flux",
                storage.predicted_face_flux[GAMERA_NO_I][0],
                expected_predicted_face, 2.0e-16);
    for (size_t cell = 0; cell < cell_count; ++cell) {
      if (!isfinite(storage.predicted_cell_magnetic[cell].value[0]) ||
          !isfinite(storage.predicted_conserved[cell * 5])) {
        fprintf(stderr, "FAIL grid predictor produced non-finite state\n");
        ++failures;
        break;
      }
    }
  }
  gamera_no_storage_destroy(&storage);
  gamera_no_grid_destroy(&grid);
}

typedef struct indexed_source_test_context {
  size_t expected_cell;
  int point_calls;
  int indexed_calls;
} indexed_source_test_context;

static int zero_point_source(
    gamera_no_vec3 point,
    const double predicted[GAMERA_NO_FLUX_COUNT], double time,
    void *context, double rate[GAMERA_NO_FLUX_COUNT]) {
  indexed_source_test_context *state = context;
  if (state == NULL || predicted == NULL || rate == NULL ||
      !isfinite(point.value[0]) || !isfinite(time)) {
    return -1;
  }
  ++state->point_calls;
  return 0;
}

static int indexed_energy_source(
    size_t cell, gamera_no_vec3 point,
    const double predicted[GAMERA_NO_FLUX_COUNT], double time,
    void *context, double rate[GAMERA_NO_FLUX_COUNT]) {
  indexed_source_test_context *state = context;
  if (state == NULL || predicted == NULL || rate == NULL ||
      cell != state->expected_cell || !isfinite(point.value[0]) ||
      !isfinite(time)) {
    return -1;
  }
  ++state->indexed_calls;
  rate[GAMERA_NO_FLUX_ENERGY] = 2.0;
  return 0;
}

static int indexed_nonfinite_source(
    size_t cell, gamera_no_vec3 point,
    const double predicted[GAMERA_NO_FLUX_COUNT], double time,
    void *context, double rate[GAMERA_NO_FLUX_COUNT]) {
  indexed_source_test_context *state = context;
  if (state == NULL || predicted == NULL || rate == NULL ||
      cell != state->expected_cell || !isfinite(point.value[0]) ||
      !isfinite(time)) {
    return -1;
  }
  ++state->indexed_calls;
  rate[GAMERA_NO_FLUX_DENSITY] = 7.0;
  rate[GAMERA_NO_FLUX_ENERGY] = NAN;
  return 0;
}

static void test_complete_uniform_advance(void) {
  const size_t extent[3] = {9, 9, 9};
  const size_t vertex_extent[3] = {10, 10, 10};
  gamera_no_vec3 vertices[1000];
  for (size_t i = 0; i < vertex_extent[0]; ++i) {
    for (size_t j = 0; j < vertex_extent[1]; ++j) {
      for (size_t k = 0; k < vertex_extent[2]; ++k) {
        vertices[gamera_no_index3(vertex_extent, i, j, k)] =
            (gamera_no_vec3){{(double)i + 0.12 * (double)j,
                              (double)j + 0.07 * (double)k,
                              (double)k + 0.05 * (double)i}};
      }
    }
  }
  gamera_no_grid grid = {0};
  gamera_no_storage storage = {0};
  if (gamera_no_grid_create(extent, vertices, &grid) != 0 ||
      gamera_no_storage_create(extent, &storage) != 0) {
    fprintf(stderr, "FAIL complete advance setup returned an error\n");
    ++failures;
    gamera_no_storage_destroy(&storage);
    gamera_no_grid_destroy(&grid);
    return;
  }

  const gamera_no_primitive primitive =
      {1.2, {{0.23, -0.17, 0.09}}, 0.74};
  const gamera_no_vec3 magnetic = {{0.38, -0.21, 0.16}};
  double conserved[GAMERA_NO_FLUX_COUNT];
  (void)gamera_no_primitive_to_conserved(&primitive, 5.0 / 3.0, 1.0e-12,
                                         1.0e-12, conserved);
  const size_t cell_count = gamera_no_element_count3(extent);
  for (size_t cell = 0; cell < cell_count; ++cell) {
    memcpy(&storage.conserved[cell * GAMERA_NO_FLUX_COUNT], conserved,
           sizeof(conserved));
    memcpy(&storage.old_conserved[cell * GAMERA_NO_FLUX_COUNT], conserved,
           sizeof(conserved));
    storage.cell_magnetic[cell] = magnetic;
    storage.old_cell_magnetic[cell] = magnetic;
  }
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    const size_t face_count =
        gamera_no_element_count3(grid.face[direction].extent);
    for (size_t face = 0; face < face_count; ++face) {
      const double flux =
          dot(magnetic, grid.face[direction].value[face].area_vector);
      storage.face_flux[direction][face] = flux;
      storage.old_face_flux[direction][face] = flux;
    }
  }

  const size_t active_lower[3] = {4, 4, 4};
  const size_t active_upper[3] = {5, 5, 5};
  const size_t active_cell = gamera_no_index3(extent, 4, 4, 4);
  indexed_source_test_context source_context = {active_cell, 0, 0};
  const gamera_no_advance_options options = {
      {5.0 / 3.0, 1.0e-12, 1.0e-12, 1.0, 0.15, 0.2, 2.3,
       true,        true,      true,      false},
      {1.0e-12, 1.0, 0.6, 2.3, 0.8, 0.01, true, false},
      {5.0 / 3.0, 1.0e-12, 1.0e-12, 2.3, true, true, false},
      zero_point_source, &source_context, indexed_energy_source,
      &source_context, 0.0};

  /* Reject a bad source row before adding any of its otherwise-finite fields. */
  indexed_source_test_context bad_source_context = {active_cell, 0, 0};
  gamera_no_advance_options bad_options = options;
  bad_options.cell_source = NULL;
  bad_options.source_context = NULL;
  bad_options.indexed_cell_source = indexed_nonfinite_source;
  bad_options.indexed_source_context = &bad_source_context;
  if (gamera_no_advance(&storage, &grid, active_lower, active_upper, 0.5,
                        0.01, &bad_options, NULL) == 0) {
    fprintf(stderr, "FAIL nonfinite indexed source was accepted\n");
    ++failures;
  }
  for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
    expect_near("nonfinite source is cell-transactional",
                storage.hydro_rate[active_cell * GAMERA_NO_FLUX_COUNT +
                                   (size_t)variable],
                0.0, 2.0e-12);
  }

  const double old_net_flux =
      gamera_no_cell_net_flux(&storage, 4, 4, 4);
  if (gamera_no_advance(&storage, &grid, active_lower, active_upper, 0.5,
                        0.01, &options, NULL) != 0) {
    fprintf(stderr, "FAIL complete non-orthogonal uniform advance\n");
    ++failures;
  } else {
    for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
      const double expected =
          conserved[variable] +
          (variable == GAMERA_NO_FLUX_ENERGY ? 0.02 : 0.0);
      expect_near("complete uniform conserved state",
                  storage.conserved[active_cell * GAMERA_NO_FLUX_COUNT +
                                    (size_t)variable],
                  expected, 2.0e-12);
      expect_near(
          "complete uniform hydro rate",
          storage.hydro_rate[active_cell * GAMERA_NO_FLUX_COUNT +
                             (size_t)variable],
          variable == GAMERA_NO_FLUX_ENERGY ? 2.0 : 0.0, 2.0e-12);
    }
    for (int component = 0; component < GAMERA_NO_DIM; ++component) {
      expect_near("complete uniform magnetic state",
                  storage.cell_magnetic[active_cell].value[component],
                  magnetic.value[component], 2.0e-12);
      expect_near("complete uniform Maxwell rate",
                  storage.maxwell_rate[active_cell].value[component], 0.0,
                  2.0e-12);
    }
    expect_near("complete uniform divergence",
                gamera_no_cell_net_flux(&storage, 4, 4, 4), old_net_flux,
                2.0e-13);
    if (source_context.point_calls != 1 ||
        source_context.indexed_calls != 1) {
      fprintf(stderr,
              "FAIL additive point/indexed source calls=%d/%d expected=1/1\n",
              source_context.point_calls, source_context.indexed_calls);
      ++failures;
    }
  }
  gamera_no_storage_destroy(&storage);
  gamera_no_grid_destroy(&grid);
}

static void test_inner_wall_mass_energy_flux_trap(void) {
  const size_t extent[3] = {10U, 10U, 10U};
  gamera_no_storage storage;
  if (gamera_no_storage_create(extent, &storage) != 0) {
    fprintf(stderr, "FAIL inner-wall flux-trap storage allocation\n");
    ++failures;
    return;
  }
  const size_t lower[3] = {4U, 4U, 4U};
  const size_t upper[3] = {6U, 6U, 6U};
  const size_t face_a =
      gamera_no_index3(storage.face_extent[GAMERA_NO_I], 4U, 4U, 4U);
  const size_t face_b =
      gamera_no_index3(storage.face_extent[GAMERA_NO_I], 4U, 4U, 5U);
  const size_t face_c =
      gamera_no_index3(storage.face_extent[GAMERA_NO_I], 4U, 5U, 4U);
  double *a = &storage.fluid_face_flux[GAMERA_NO_I]
                                      [face_a * GAMERA_NO_FLUX_COUNT];
  double *b = &storage.fluid_face_flux[GAMERA_NO_I]
                                      [face_b * GAMERA_NO_FLUX_COUNT];
  double *c = &storage.fluid_face_flux[GAMERA_NO_I]
                                      [face_c * GAMERA_NO_FLUX_COUNT];
  a[GAMERA_NO_FLUX_DENSITY] = 2.0;
  a[GAMERA_NO_FLUX_ENERGY] = 3.0;
  a[GAMERA_NO_FLUX_MOMENTUM_X] = 4.0;
  b[GAMERA_NO_FLUX_DENSITY] = 2.0;
  b[GAMERA_NO_FLUX_ENERGY] = -3.0;
  c[GAMERA_NO_FLUX_DENSITY] = -2.0;
  c[GAMERA_NO_FLUX_ENERGY] = 3.0;

  size_t count = 0U;
  if (gamera_no_trap_inner_outward_mass_energy_flux(
          &storage, lower, upper, &count) != 0) {
    fprintf(stderr, "FAIL inner-wall flux trap returned an error\n");
    ++failures;
  } else {
    expect_near("inner-wall injected mass", a[GAMERA_NO_FLUX_DENSITY], 0.0,
                0.0);
    expect_near("inner-wall injected energy", a[GAMERA_NO_FLUX_ENERGY], 0.0,
                0.0);
    expect_near("inner-wall momentum untouched",
                a[GAMERA_NO_FLUX_MOMENTUM_X], 4.0, 0.0);
    expect_near("inner-wall inward energy retained",
                b[GAMERA_NO_FLUX_ENERGY], -3.0, 0.0);
    expect_near("inner-wall inward mass retained",
                c[GAMERA_NO_FLUX_DENSITY], -2.0, 0.0);
    expect_near("inner-wall energy without outward mass retained",
                c[GAMERA_NO_FLUX_ENERGY], 3.0, 0.0);
    if (count != 2U) {
      fprintf(stderr,
              "FAIL inner-wall clamped face count: actual=%zu expected=2\n",
              count);
      ++failures;
    }
  }
  gamera_no_storage_destroy(&storage);
}

static void test_kaiju_chillout(void) {
  const size_t extent[3] = {1U, 1U, 1U};
  const size_t vertex_extent[3] = {2U, 2U, 2U};
  gamera_no_vec3 vertices[8];
  for (size_t i = 0; i < 2U; ++i) {
    for (size_t j = 0; j < 2U; ++j) {
      for (size_t k = 0; k < 2U; ++k) {
        vertices[gamera_no_index3(vertex_extent, i, j, k)] =
            (gamera_no_vec3){{2.0 + (double)i, (double)j, (double)k}};
      }
    }
  }
  gamera_no_grid grid = {0};
  gamera_no_storage storage = {0};
  if (gamera_no_grid_create(extent, vertices, &grid) != 0 ||
      gamera_no_storage_create(extent, &storage) != 0) {
    fprintf(stderr, "FAIL Kaiju ChillOut setup returned an error\n");
    ++failures;
    gamera_no_storage_destroy(&storage);
    gamera_no_grid_destroy(&grid);
    return;
  }
  const double gamma = 5.0 / 3.0;
  const double dt_test = 0.01;
  const gamera_no_primitive hot = {1.0, {{0.2, 0.0, 0.0}}, 600.0};
  (void)gamera_no_primitive_to_conserved(
      &hot, gamma, 1.0e-6, 1.0e-8, storage.conserved);
  size_t low_count = 0U;
  size_t hot_count = 0U;
  double max_cs = 0.0;
  double max_p = 0.0;
  const size_t lower[3] = {0U, 0U, 0U};
  const size_t upper[3] = {1U, 1U, 1U};
  if (gamera_no_apply_kaiju_chillout(
          &storage, &grid, lower, upper, dt_test, gamma, 1.0e-6,
          1.0e-8, 10.0, 1.0e-3, 1.0e-7, &low_count, &hot_count,
          &max_cs, &max_p) != 0) {
    fprintf(stderr, "FAIL Kaiju ChillOut returned an error\n");
    ++failures;
  } else {
    gamera_no_primitive cooled;
    (void)gamera_no_conserved_to_primitive(
        storage.conserved, gamma, 1.0e-6, 1.0e-8, &cooled);
    const double radius = sqrt(2.5 * 2.5 + 0.5 * 0.5 + 0.5 * 0.5);
    const double dipole_l = radius * radius * radius /
                            (2.5 * 2.5 + 0.5 * 0.5);
    const double sound = sqrt(gamma * hot.pressure / hot.density);
    const double target = hot.density * 15.0 * 15.0 / gamma;
    const double expected =
        hot.pressure - dt_test * sound / dipole_l * (hot.pressure - target);
    expect_near("Kaiju ChillOut hot pressure", cooled.pressure, expected,
                2.0e-12);
    if (low_count != 0U || hot_count != 1U) {
      fprintf(stderr,
              "FAIL Kaiju ChillOut counts low/hot=%zu/%zu expected=0/1\n",
              low_count, hot_count);
      ++failures;
    }
    expect_near("Kaiju ChillOut maximum sound", max_cs, sound, 2.0e-12);
    expect_near("Kaiju ChillOut maximum pressure", max_p, hot.pressure,
                2.0e-12);
  }

  const gamera_no_primitive tenuous =
      {1.0e-4, {{0.0, 0.0, 0.0}}, 1.0};
  (void)gamera_no_primitive_to_conserved(
      &tenuous, gamma, 1.0e-6, 1.0e-8, storage.conserved);
  if (gamera_no_apply_kaiju_chillout(
          &storage, &grid, lower, upper, dt_test, gamma, 1.0e-6,
          1.0e-8, 10.0, 1.0e-3, 1.0e-7, &low_count, &hot_count,
          &max_cs, &max_p) != 0) {
    fprintf(stderr, "FAIL Kaiju low-density ChillOut returned an error\n");
    ++failures;
  } else {
    gamera_no_primitive cooled;
    (void)gamera_no_conserved_to_primitive(
        storage.conserved, gamma, 1.0e-6, 1.0e-8, &cooled);
    expect_near("Kaiju low-density pressure", cooled.pressure, 1.0e-8,
                1.0e-15);
    if (low_count != 1U) {
      fprintf(stderr,
              "FAIL Kaiju low-density count=%zu expected=1\n", low_count);
      ++failures;
    }
  }
  gamera_no_storage_destroy(&storage);
  gamera_no_grid_destroy(&grid);
}

static void test_yinyang_composite_owner(void) {
  const gamera_no_yinyang_angular_domain domain = {
      0.5, 2.5, -2.0, 2.0, 64U, 72U, 2.0e-2};
  int owner = -1;
  int overlap = -1;
  double margin[2] = {NAN, NAN};
  if (gamera_no_yinyang_composite_owner(
          (gamera_no_vec3){{1.0, 0.0, 0.0}}, &domain, &owner, &overlap,
          margin) != 0 ||
      owner != GAMERA_NO_YIN_PATCH || overlap != 0) {
    fprintf(stderr, "FAIL Yin-only composite ownership\n");
    ++failures;
  }
  if (gamera_no_yinyang_composite_owner(
          (gamera_no_vec3){{0.0, 0.0, 1.0}}, &domain, &owner, &overlap,
          margin) != 0 ||
      owner != GAMERA_NO_YANG_PATCH || overlap != 0) {
    fprintf(stderr, "FAIL Yang-only composite ownership\n");
    ++failures;
  }
  if (gamera_no_yinyang_composite_owner(
          (gamera_no_vec3){{0.0, 1.0, 1.0}}, &domain, &owner, &overlap,
          margin) != 0 ||
      owner != GAMERA_NO_YIN_PATCH || overlap != 1) {
    fprintf(stderr, "FAIL tied-overlap composite ownership\n");
    ++failures;
  } else {
    expect_near("Yin/Yang tied overlap margin", margin[0], margin[1],
                2.0e-13);
  }
}

int main(void) {
  test_unit_cube();
  test_affine_skew_cell();
  test_warped_cell();
  test_edge_field_solve();
  test_edge_geometry_and_normal_interpolation();
  test_ct_divergence_identity();
  test_edge_emf();
  test_reconstruction();
  test_fortran_flux_gold();
  test_nuclear_hogs();
  test_fortran_state_gold();
  test_grid_and_storage_layout();
  test_grid_level_ct_and_predictor();
  test_uniform_face_sweep();
  test_uniform_edge_emf_sweep();
  test_inner_wall_mass_energy_flux_trap();
  test_kaiju_chillout();
  test_yinyang_composite_owner();
  test_complete_uniform_advance();

  if (failures != 0) {
    fprintf(stderr, "%d non-orthogonal kernel checks failed\n", failures);
    return EXIT_FAILURE;
  }
  printf("All non-orthogonal geometry/operator checks passed.\n");
  return EXIT_SUCCESS;
}
