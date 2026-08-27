#include "nonorthogonal_storage.h"

#include "nonorthogonal_grid.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int allocate_items(size_t count, size_t item_size, void **items) {
  if (items == NULL || count == 0 || item_size == 0 ||
      count > SIZE_MAX / item_size) {
    return -1;
  }
  *items = calloc(count, item_size);
  return *items == NULL ? -1 : 0;
}

int gamera_no_storage_create(const size_t cell_extent[3],
                             gamera_no_storage *storage) {
  if (cell_extent == NULL || storage == NULL) {
    return -1;
  }
  memset(storage, 0, sizeof(*storage));
  for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
    if (cell_extent[axis] == 0 || cell_extent[axis] == SIZE_MAX) {
      return -1;
    }
    storage->cell_extent[axis] = cell_extent[axis];
  }

  const size_t cell_count = gamera_no_element_count3(storage->cell_extent);
  if (cell_count == 0 || cell_count > SIZE_MAX / GAMERA_NO_FLUX_COUNT) {
    return -1;
  }
  const size_t conserved_count = cell_count * GAMERA_NO_FLUX_COUNT;
  if (allocate_items(conserved_count, sizeof(*storage->conserved),
                     (void **)&storage->conserved) != 0 ||
      allocate_items(conserved_count, sizeof(*storage->old_conserved),
                     (void **)&storage->old_conserved) != 0 ||
      allocate_items(conserved_count, sizeof(*storage->predicted_conserved),
                     (void **)&storage->predicted_conserved) != 0 ||
      allocate_items(cell_count, sizeof(*storage->cell_magnetic),
                     (void **)&storage->cell_magnetic) != 0 ||
      allocate_items(cell_count, sizeof(*storage->old_cell_magnetic),
                     (void **)&storage->old_cell_magnetic) != 0 ||
      allocate_items(cell_count, sizeof(*storage->predicted_cell_magnetic),
                     (void **)&storage->predicted_cell_magnetic) != 0 ||
      allocate_items(conserved_count, sizeof(*storage->hydro_rate),
                     (void **)&storage->hydro_rate) != 0 ||
      allocate_items(cell_count, sizeof(*storage->maxwell_rate),
                     (void **)&storage->maxwell_rate) != 0) {
    gamera_no_storage_destroy(storage);
    return -1;
  }

  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    for (int axis = 0; axis < GAMERA_NO_DIM; ++axis) {
      storage->face_extent[direction][axis] =
          cell_extent[axis] + (axis == direction ? 1U : 0U);
      storage->edge_extent[direction][axis] =
          cell_extent[axis] + (axis == direction ? 0U : 1U);
    }
    const size_t face_count =
        gamera_no_element_count3(storage->face_extent[direction]);
    const size_t edge_count =
        gamera_no_element_count3(storage->edge_extent[direction]);
    if (face_count == 0 || face_count > SIZE_MAX / GAMERA_NO_FLUX_COUNT ||
        allocate_items(face_count, sizeof(*storage->face_flux[direction]),
                       (void **)&storage->face_flux[direction]) != 0 ||
        allocate_items(face_count, sizeof(*storage->old_face_flux[direction]),
                       (void **)&storage->old_face_flux[direction]) != 0 ||
        allocate_items(
            face_count, sizeof(*storage->predicted_face_flux[direction]),
            (void **)&storage->predicted_face_flux[direction]) != 0 ||
        allocate_items(face_count * GAMERA_NO_FLUX_COUNT,
                       sizeof(*storage->fluid_face_flux[direction]),
                       (void **)&storage->fluid_face_flux[direction]) != 0 ||
        allocate_items(face_count,
                       sizeof(*storage->maxwell_face_flux[direction]),
                       (void **)&storage->maxwell_face_flux[direction]) != 0 ||
        allocate_items(edge_count, sizeof(*storage->edge_emf[direction]),
                       (void **)&storage->edge_emf[direction]) != 0) {
      gamera_no_storage_destroy(storage);
      return -1;
    }
  }
  return 0;
}

void gamera_no_storage_destroy(gamera_no_storage *storage) {
  if (storage == NULL) {
    return;
  }
  free(storage->conserved);
  free(storage->old_conserved);
  free(storage->predicted_conserved);
  free(storage->cell_magnetic);
  free(storage->old_cell_magnetic);
  free(storage->predicted_cell_magnetic);
  free(storage->hydro_rate);
  free(storage->maxwell_rate);
  for (int direction = 0; direction < GAMERA_NO_DIM; ++direction) {
    free(storage->face_flux[direction]);
    free(storage->old_face_flux[direction]);
    free(storage->predicted_face_flux[direction]);
    free(storage->edge_emf[direction]);
    free(storage->fluid_face_flux[direction]);
    free(storage->maxwell_face_flux[direction]);
  }
  memset(storage, 0, sizeof(*storage));
}
