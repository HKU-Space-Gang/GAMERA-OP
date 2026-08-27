#include <stdlib.h>

#include "log.h"

#include "utils.h"

/*
    Allocate a 2D array [nx][ny].
    Usage: array = (double **)alloc_2d_array(nx, ny, sizeof(double));
*/
void **alloc_2d_array(size_t nx, size_t ny, size_t size) {
  void **array;
  size_t i;

  if ((array = (void **)calloc(nx, sizeof(void *))) == NULL) {
    log_error("failed 1 alloc_2d (%zu)", nx);
    return NULL;
  }

  if ((array[0] = (void *)calloc(nx * ny, size)) == NULL) {
    log_error("failed 2 alloc_2d (%zu,%zu)", nx, ny);
    free((void *)array);
    return NULL;
  }

  for (i = 1; i < nx; i++) {
    array[i] = (void *)((unsigned char *)array[0] + i * ny * size);
  }

  return array;
}

/*
    Free a 2D array [nx][ny].
    Usage: free_2d_array((void **)array);
*/
void free_2d_array(void **array) {
  free(array[0]);
  free(array);
}

/*
    Allocate a 3D array [nx][ny][nz] using pre-allocated contiguous data.
    Usage: array = (double ***)alloc_3d_array_on_data(data, nx, ny, nz,
   sizeof(double));
*/
static void ***alloc_3d_array_on_data(void *data, size_t nx, size_t ny,
                                      size_t nz, size_t size) {
  void ***array;
  size_t i, j;

  if ((array = (void ***)calloc(nx, sizeof(void **))) == NULL) {
    log_error("failed 1 alloc_3d (%zu)", nx);
    return NULL;
  }

  if ((array[0] = (void **)calloc(nx * ny, sizeof(void *))) == NULL) {
    log_error("failed 2 alloc_3d (%zu)", nx * ny);
    free((void *)array);
    return NULL;
  }

  for (i = 1; i < nx; i++) {
    array[i] = (void **)((unsigned char *)array[0] + i * ny * sizeof(void *));
  }
  array[0][0] = data;

  for (j = 1; j < ny; j++) {
    array[0][j] = (void **)((unsigned char *)array[0][j - 1] + nz * size);
  }

  for (i = 1; i < nx; i++) {
    array[i][0] = (void **)((unsigned char *)array[i - 1][0] + ny * nz * size);
    for (j = 1; j < ny; j++) {
      array[i][j] = (void **)((unsigned char *)array[i][j - 1] + nz * size);
    }
  }

  return array;
}

/*
    Allocate a contiguous 4D array [nt][nx][ny][nz].
    Usage: array = (double ****)alloc_4d_array_contiguous(nt, nx, ny, nz,
   sizeof(double));
*/
void ****alloc_4d_array_contiguous(size_t nt, size_t nx, size_t ny, size_t nz,
                                   size_t size) {
  void ****array;
  if ((array = (void ****)calloc(nt, sizeof(void ***))) == NULL) {
    log_error("Failed to allocate memory for 4D array.");
    return NULL;
  }
  void *data = (void *)calloc(nt * nx * ny * nz * size, 1);
  if (data == NULL) {
    log_error("Failed to allocate memory for 4D array.");
    free(array);
    return NULL;
  }
  for (int t = 0; t < nt; t++) {
    array[t] = alloc_3d_array_on_data(data + t * nx * ny * nz * size, nx, ny,
                                      nz, size);
    if (array[t] == NULL) {
      free(data);
      for (int tt = 0; tt < t; tt++) {
        free(array[tt][0]);
        free(array[tt]);
      }
      free(array);
      return NULL;
    }
  }
  return array;
}

/*
    Free a contiguous 4D array [nt][nx][ny][nz].
    Usage: free_4d_array_contiguous((void ****)array, nt);
*/
void free_4d_array_contiguous(void ****array, size_t nt) {
  // array[0][0][0] is the address of the first element of the 1D data array
  free((void *)array[0][0][0]);
  for (int t = 0; t < nt; t++) {
    free(array[t][0]);
    free(array[t]);
  }
  free(array);
}

/*
    Allocate a 5D array [ns][nt][nx][ny][nz] using contiguous 4D blocks.
    Usage: array = (double *****)alloc_5d_array_with_4d_contiguous(ns, nt, nx,
   ny, nz, sizeof(double));
*/
void *****alloc_5d_array_with_4d_contiguous(size_t ns, size_t nt, size_t nx,
                                            size_t ny, size_t nz, size_t size) {
  void *****array;
  if ((array = (void *****)calloc(ns, sizeof(void ****))) == NULL) {
    log_error("Failed to allocate memory for 5D array.");
    return NULL;
  }
  for (int s = 0; s < ns; s++) {
    array[s] = alloc_4d_array_contiguous(nt, nx, ny, nz, size);
    if (array[s] == NULL) {
      for (int ss = 0; ss < s; ss++) {
        free_4d_array_contiguous(array[ss], nt);
      }
      free(array);
      return NULL;
    }
  }
  return array;
}

/*
    Free a 5D array [ns][nt][nx][ny][nz].
    Usage: free_5d_array_with_4d_contiguous((void *****)array, ns, nt);
*/
void free_5d_array_with_4d_contiguous(void *****array, size_t ns, size_t nt) {
  for (int s = 0; s < ns; s++) {
    free_4d_array_contiguous(array[s], nt);
  }
  free(array);
}


/*
    Allocate a 3D array [nx][ny][nz].
    Usage: array = (double ***)alloc_3d_array(nx, ny, nz, sizeof(double));
*/
void ***alloc_3d_array(size_t nx, size_t ny, size_t nz, size_t size) {
  void ***array;
  size_t i, j;

  if ((array = (void ***)calloc(nx, sizeof(void **))) == NULL) {
    log_error("failed 1 alloc_3d (%zu)", nx);
    return NULL;
  }

  if ((array[0] = (void **)calloc(nx * ny, sizeof(void *))) == NULL) {
    log_error("failed 2 alloc_3d (%zu)", nx * ny);
    free((void *)array);
    return NULL;
  }

  for (i = 1; i < nx; i++) {
    array[i] = (void **)((unsigned char *)array[0] + i * ny * sizeof(void *));
  }

  if ((array[0][0] = (void *)calloc(nx * ny * nz, size)) == NULL) {
    log_error("failed 3 alloc_3d (%zu,%zu,%zu,%zu)", nx, ny, ny, size);
    free((void *)array[0]);
    free((void *)array);
    return NULL;
  }

  for (j = 1; j < ny; j++) {
    array[0][j] = (void **)((unsigned char *)array[0][j - 1] + nz * size);
  }

  for (i = 1; i < nx; i++) {
    array[i][0] = (void **)((unsigned char *)array[i - 1][0] + ny * nz * size);
    for (j = 1; j < ny; j++) {
      array[i][j] = (void **)((unsigned char *)array[i][j - 1] + nz * size);
    }
  }

  return array;
}

/*
    Free a 3D array [nx][ny][nz].
    Usage: free_3d_array((void ***)array);
*/
void free_3d_array(void ***array) {
  free(array[0][0]);
  free(array[0]);
  free(array);
}

/*
    Allocate a 4D array [nt][nx][ny][nz].
    Usage: array = (double ****)alloc_4d_array(nt, nx, ny, nz, sizeof(double));
*/
void ****alloc_4d_array(size_t nt, size_t nx, size_t ny, size_t nz,
                        size_t size) {
  void ****array;

  if ((array = (void ****)calloc(nt, sizeof(void ***))) == NULL) {
    log_error("failed 1 alloc_4d (%zu,%zu,%zu,%zu)", nt, nx, ny, nz);
    return NULL;
  }

  for (int t = 0; t < nt; t++) {
    if ((array[t] = alloc_3d_array(nx, ny, nz, size)) == NULL) {
      for (int tt = 0; tt < t; tt++) {
        free_3d_array(array[tt]);
      }
    }
  }

  return array;
}

/*
    Free a 4D array [nt][nx][ny][nz].
    Usage: free_4d_array((void ****)array, size_t nt);
*/
void free_4d_array(void ****array, size_t nt) {
  for (int t = 0; t < nt; t++) {
    free_3d_array(array[t]);
  }
  free(array);
}

/*
    Allocate a 5D array [ns][nt][nx][ny][nz].
    Usage: array = (double *****)alloc_5d_array(ns, nt, nx, ny, nz,
   sizeof(double));
*/
void *****alloc_5d_array(size_t ns, size_t nt, size_t nx, size_t ny, size_t nz,
                         size_t size) {
  void *****array;

  if ((array = (void *****)calloc(ns, sizeof(void ****))) == NULL) {
    log_error("failed 1 alloc_5d (%zu,%zu,%zu,%zu,%zu)", ns, nt, nx, ny, nz);
    return NULL;
  }

  for (int s = 0; s < ns; s++) {
    if ((array[s] = alloc_4d_array(nt, nx, ny, nz, size)) == NULL) {
      for (int ss = 0; ss < s; ss++) {
        free_4d_array(array[ss], nt);
      }
    }
  }

  return array;
}

/*
    Free a 5D array [ns][nt][nx][ny][nz].
    Usage: free_5d_array((void *****)array, size_t ns, size_t nt);
*/
void free_5d_array(void *****array, size_t ns, size_t nt) {
  for (int s = 0; s < ns; s++) {
    free_4d_array(array[s], nt);
  }
  free(array);
}
