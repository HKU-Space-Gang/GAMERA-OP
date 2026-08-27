#include "nonorthogonal_legacy_adapter.h"

#include "nonorthogonal_state.h"
#include "nonorthogonal_step.h"

#include "config.h"

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static gamera_no_grid adapter_grid;
static gamera_no_storage adapter_storage;
static bool adapter_ready;

static const int current_face_slot[3] = {mag_bi, mag_bj, mag_bk};
static const int old_face_slot[3] = {mag_bi_p, mag_bj_p, mag_bk_p};
static const int current_cell_field_slot[3] = {mags_b1, mags_b2, mags_b3};
static const int old_cell_field_slot[3] = {mags_b1_p, mags_b2_p, mags_b3_p};
static const int total_cell_field_slot[3] = {magtot_b1, magtot_b2,
                                             magtot_b3};
static const int edge_emf_slot[3] = {efield_ei, efield_ej, efield_ek};
static const int face_geometry_slot[3] = {face_idir, face_jdir, face_kdir};
static const int edge_geometry_slot[3] = {edge_idir, edge_jdir, edge_kdir};

static double vertex_edge_length(int direction, size_t i, size_t j,
                                 size_t k) {
  size_t end[3] = {i, j, k};
  ++end[direction];
  const double delta[3] = {x1[end[0]][end[1]][end[2]] - x1[i][j][k],
                           x2[end[0]][end[1]][end[2]] - x2[i][j][k],
                           x3[end[0]][end[1]][end[2]] - x3[i][j][k]};
  return sqrt(delta[0] * delta[0] + delta[1] * delta[1] +
              delta[2] * delta[2]);
}

static int build_geometry(void) {
  const size_t cell_extent[3] = {(size_t)config.NI - 1U,
                                 (size_t)config.NJ - 1U,
                                 (size_t)config.NK - 1U};
  const size_t vertex_extent[3] = {(size_t)config.NI, (size_t)config.NJ,
                                   (size_t)config.NK};
  const size_t vertex_count = gamera_no_element_count3(vertex_extent);
  if (vertex_count == 0 || vertex_count > SIZE_MAX / sizeof(gamera_no_vec3)) {
    return -1;
  }
  gamera_no_vec3 *vertices =
      (gamera_no_vec3 *)malloc(vertex_count * sizeof(*vertices));
  if (vertices == NULL) {
    return -1;
  }
  for (size_t i = 0; i < vertex_extent[0]; ++i) {
    for (size_t j = 0; j < vertex_extent[1]; ++j) {
      for (size_t k = 0; k < vertex_extent[2]; ++k) {
        const size_t vertex = gamera_no_index3(vertex_extent, i, j, k);
        vertices[vertex] =
            (gamera_no_vec3){{x1[i][j][k], x2[i][j][k], x3[i][j][k]}};
      }
    }
  }
  const int status =
      gamera_no_grid_create(cell_extent, vertices, &adapter_grid);
  free(vertices);
  return status;
}

static int import_face_state(double *destination[3],
                             const int legacy_slot[3]) {
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    for (size_t i = 0; i < adapter_grid.face[direction].extent[0]; ++i) {
      for (size_t j = 0; j < adapter_grid.face[direction].extent[1]; ++j) {
        for (size_t k = 0; k < adapter_grid.face[direction].extent[2]; ++k) {
          const size_t face = gamera_no_index3(
              adapter_grid.face[direction].extent, i, j, k);
          destination[direction][face] =
              gem[legacy_slot[direction]][i][j][k] *
              adapter_grid.face[direction].value[face].area;
        }
      }
    }
  }
  return 0;
}

static int import_primitive_state(double *destination, bool old_state) {
  const int offset = old_state ? gas_rho_p - gas_rho : 0;
  for (size_t i = 0; i < adapter_grid.cell_extent[0]; ++i) {
    for (size_t j = 0; j < adapter_grid.cell_extent[1]; ++j) {
      for (size_t k = 0; k < adapter_grid.cell_extent[2]; ++k) {
        const size_t cell = gamera_no_index3(adapter_grid.cell_extent, i, j, k);
        const gamera_no_primitive primitive = {
            gas[0][gas_rho + offset][i][j][k],
            {{gas[0][gas_v1 + offset][i][j][k],
              gas[0][gas_v2 + offset][i][j][k],
              gas[0][gas_v3 + offset][i][j][k]}},
            gas[0][gas_p + offset][i][j][k]};
        if (gamera_no_primitive_to_conserved(
                &primitive, gamma_val, rho_floor, p_floor,
                &destination[cell * GAMERA_NO_FLUX_COUNT]) != 0) {
          return -1;
        }
      }
    }
  }
  return 0;
}

static void export_geometry(void) {
  for (size_t i = 0; i < adapter_grid.cell_extent[0]; ++i) {
    for (size_t j = 0; j < adapter_grid.cell_extent[1]; ++j) {
      for (size_t k = 0; k < adapter_grid.cell_extent[2]; ++k) {
        const size_t cell = gamera_no_index3(adapter_grid.cell_extent, i, j, k);
        geo[vol_center][i][j][k] = adapter_grid.cell[cell].volume;
        x1ctr[i][j][k] = adapter_grid.cell[cell].centroid.value[0];
        x2ctr[i][j][k] = adapter_grid.cell[cell].centroid.value[1];
        x3ctr[i][j][k] = adapter_grid.cell[cell].centroid.value[2];
      }
    }
  }
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    for (size_t i = 0; i < adapter_grid.face[direction].extent[0]; ++i) {
      for (size_t j = 0; j < adapter_grid.face[direction].extent[1]; ++j) {
        for (size_t k = 0; k < adapter_grid.face[direction].extent[2]; ++k) {
          const size_t face = gamera_no_index3(
              adapter_grid.face[direction].extent, i, j, k);
          geo[face_geometry_slot[direction]][i][j][k] =
              adapter_grid.face[direction].value[face].area;
        }
      }
    }
    for (size_t i = 0; i < adapter_grid.edge[direction].extent[0]; ++i) {
      for (size_t j = 0; j < adapter_grid.edge[direction].extent[1]; ++j) {
        for (size_t k = 0; k < adapter_grid.edge[direction].extent[2]; ++k) {
          geo[edge_geometry_slot[direction]][i][j][k] =
              vertex_edge_length(direction, i, j, k);
        }
      }
    }
  }
}

static int allocate_adapter(void) {
  if (adapter_ready || config.NI < 2 || config.NJ < 2 || config.NK < 2) {
    return -1;
  }
  if (build_geometry() != 0 ||
      gamera_no_storage_create(adapter_grid.cell_extent, &adapter_storage) !=
          0) {
    gamera_no_legacy_adapter_destroy();
    return -1;
  }
  adapter_ready = true;
  export_geometry();
  return 0;
}

int gamera_no_legacy_adapter_create_empty(void) {
  return allocate_adapter();
}

int gamera_no_legacy_adapter_create(bool import_restart_history) {
  if (allocate_adapter() != 0) {
    return -1;
  }
  if (gamera_no_legacy_import_current() != 0) {
    gamera_no_legacy_adapter_destroy();
    return -1;
  }
  if (import_restart_history) {
    double *old_flux[3] = {adapter_storage.old_face_flux[GAMERA_NO_I],
                           adapter_storage.old_face_flux[GAMERA_NO_J],
                           adapter_storage.old_face_flux[GAMERA_NO_K]};
    if (import_primitive_state(adapter_storage.old_conserved, true) != 0 ||
        import_face_state(old_flux, old_face_slot) != 0) {
      gamera_no_legacy_adapter_destroy();
      return -1;
    }
    const double *old_flux_const[3] = {old_flux[GAMERA_NO_I],
                                       old_flux[GAMERA_NO_J],
                                       old_flux[GAMERA_NO_K]};
    if (gamera_no_recover_magnetic_field(
            &adapter_grid, old_flux_const,
            adapter_storage.old_cell_magnetic) != 0) {
      gamera_no_legacy_adapter_destroy();
      return -1;
    }
  } else if (gamera_no_save_current_as_old(&adapter_storage) != 0) {
    gamera_no_legacy_adapter_destroy();
    return -1;
  }
  return gamera_no_legacy_export();
}

void gamera_no_legacy_adapter_destroy(void) {
  gamera_no_storage_destroy(&adapter_storage);
  gamera_no_grid_destroy(&adapter_grid);
  adapter_ready = false;
}

int gamera_no_legacy_import_current(void) {
  if (!adapter_ready ||
      import_primitive_state(adapter_storage.conserved, false) != 0) {
    return -1;
  }
  double *current_flux[3] = {adapter_storage.face_flux[GAMERA_NO_I],
                             adapter_storage.face_flux[GAMERA_NO_J],
                             adapter_storage.face_flux[GAMERA_NO_K]};
  if (import_face_state(current_flux, current_face_slot) != 0) {
    return -1;
  }
  const double *current_flux_const[3] = {
      current_flux[GAMERA_NO_I], current_flux[GAMERA_NO_J],
      current_flux[GAMERA_NO_K]};
  return gamera_no_recover_magnetic_field(
      &adapter_grid, current_flux_const, adapter_storage.cell_magnetic);
}

static void export_cell_state(void) {
  for (size_t i = 0; i < adapter_grid.cell_extent[0]; ++i) {
    for (size_t j = 0; j < adapter_grid.cell_extent[1]; ++j) {
      for (size_t k = 0; k < adapter_grid.cell_extent[2]; ++k) {
        const size_t cell = gamera_no_index3(adapter_grid.cell_extent, i, j, k);
        gamera_no_primitive primitive;
        gamera_no_primitive old_primitive;
        (void)gamera_no_conserved_to_primitive(
            &adapter_storage.conserved[cell * GAMERA_NO_FLUX_COUNT],
            gamma_val, rho_floor, p_floor, &primitive);
        (void)gamera_no_conserved_to_primitive(
            &adapter_storage.old_conserved[cell * GAMERA_NO_FLUX_COUNT],
            gamma_val, rho_floor, p_floor, &old_primitive);

        gas[0][gas_rho][i][j][k] = primitive.density;
        gas[0][gas_v1][i][j][k] = primitive.velocity.value[0];
        gas[0][gas_v2][i][j][k] = primitive.velocity.value[1];
        gas[0][gas_v3][i][j][k] = primitive.velocity.value[2];
        gas[0][gas_p][i][j][k] = primitive.pressure;
        gas[0][gas_p_S][i][j][k] = primitive.pressure;
        gas[0][gas_rho_p][i][j][k] = old_primitive.density;
        gas[0][gas_v1_p][i][j][k] = old_primitive.velocity.value[0];
        gas[0][gas_v2_p][i][j][k] = old_primitive.velocity.value[1];
        gas[0][gas_v3_p][i][j][k] = old_primitive.velocity.value[2];
        gas[0][gas_p_p][i][j][k] = old_primitive.pressure;
        gas[0][gas_p_S_p][i][j][k] = old_primitive.pressure;

        for (int variable = 0; variable < GAMERA_NO_FLUX_COUNT; ++variable) {
          gas[0][gasc_rho + variable][i][j][k] =
              adapter_storage
                  .conserved[cell * GAMERA_NO_FLUX_COUNT + (size_t)variable];
        }
        for (int component = 0; component < GAMERA_NO_DIM; ++component) {
          const double field =
              adapter_storage.cell_magnetic[cell].value[component];
          gem[current_cell_field_slot[component]][i][j][k] = field;
          gem[old_cell_field_slot[component]][i][j][k] =
              adapter_storage.old_cell_magnetic[cell].value[component];
          gem[total_cell_field_slot[component]][i][j][k] = field;
        }
      }
    }
  }
}

static void export_face_and_edge_state(void) {
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    for (size_t i = 0; i < adapter_grid.face[direction].extent[0]; ++i) {
      for (size_t j = 0; j < adapter_grid.face[direction].extent[1]; ++j) {
        for (size_t k = 0; k < adapter_grid.face[direction].extent[2]; ++k) {
          const size_t face = gamera_no_index3(
              adapter_grid.face[direction].extent, i, j, k);
          const double inverse_area =
              1.0 / adapter_grid.face[direction].value[face].area;
          gem[current_face_slot[direction]][i][j][k] =
              adapter_storage.face_flux[direction][face] * inverse_area;
          gem[old_face_slot[direction]][i][j][k] =
              adapter_storage.old_face_flux[direction][face] * inverse_area;
        }
      }
    }
    for (size_t i = 0; i < adapter_grid.edge[direction].extent[0]; ++i) {
      for (size_t j = 0; j < adapter_grid.edge[direction].extent[1]; ++j) {
        for (size_t k = 0; k < adapter_grid.edge[direction].extent[2]; ++k) {
          const size_t edge = gamera_no_index3(
              adapter_grid.edge[direction].extent, i, j, k);
          if (adapter_grid.edge[direction].valid[edge] != 0U) {
            gem[edge_emf_slot[direction]][i][j][k] =
                adapter_storage.edge_emf[direction][edge] /
                adapter_grid.edge[direction].value[edge].length;
          }
        }
      }
    }
  }
}

int gamera_no_legacy_export(void) {
  if (!adapter_ready) {
    return -1;
  }
  export_cell_state();
  export_face_and_edge_state();
  divB_max = 0.0;
  for (int i = is; i <= ie; ++i) {
    for (int j = js; j <= je; ++j) {
      for (int k = ks; k <= ke; ++k) {
        const size_t cell = gamera_no_index3(
            adapter_grid.cell_extent, (size_t)i, (size_t)j, (size_t)k);
        const double divergence =
            fabs(gamera_no_cell_net_flux(&adapter_storage, (size_t)i,
                                         (size_t)j, (size_t)k) /
                 adapter_grid.cell[cell].volume);
        divB_max = fmax(divB_max, divergence);
      }
    }
  }
  return 0;
}

gamera_no_grid *gamera_no_legacy_grid(void) {
  return adapter_ready ? &adapter_grid : NULL;
}

gamera_no_storage *gamera_no_legacy_storage(void) {
  return adapter_ready ? &adapter_storage : NULL;
}
