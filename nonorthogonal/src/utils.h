#ifndef UTILS_H
#define UTILS_H

#include <mpi.h>
#include <stddef.h>

// Reads the dimensions of a dataset from an HDF5 file.
int read_data_dims(char *filename, char *dataset_name, int *dims, int rank);

// Reads grid data from an HDF5 file.
int read_grid_data_hdf5(const char *filename, double ***x1_temp,
                        double ***x2_temp, double ***x3_temp);

int read_wind_data_hdf5(const char *filename, double *wind_time, double *wind_rho, double *wind_vx, double *wind_vy, double *wind_vz, double *wind_p,
                        double *wind_bx, double *wind_by, double *wind_bz) ;

// Generates an output filename based on a base filename, sequence number, and extension.
void get_out_filename(const char *base_filename, char *filename,
                      size_t filename_size, int seq_num, const char *extension);

int dump_check_data(const char *base_filename);

// Dumps data to an HDF5 file.
int dump_data_hdf5(const char *filename);

// Dumps additional data to an HDF5 file.
int dump_extra_data(const char *base_filename);

// Reads data from an HDF5 file.
int read_data_hdf5(const char *filename);

// Synchronizes the boundary values of a 4D field array with neighboring processes.
int boundary_exchange_4d(double ****field_arrays, int *onface_i, int *onface_j,
                         int *onface_k);

// Synchronizes the boundary values of a 5D field array with neighboring processes.
int boundary_exchange_5d(double *****field_arrays, int ns, int NI, int NJ,
                         int NK);

// Allocates memory for a 4D contiguous array.
void ****alloc_4d_array_contiguous(size_t nt, size_t nx, size_t ny, size_t nz,
                                   size_t size);

// Frees memory allocated by alloc_4d_array_contiguous.
void free_4d_array_contiguous(void ****array, size_t nt);

// Allocates memory for a 5D array with the last 4 dimensions contiguous.
void *****alloc_5d_array_with_4d_contiguous(size_t ns, size_t nt, size_t nx,
                                            size_t ny, size_t nz, size_t size);

// Frees memory allocated by alloc_5d_array_with_4d_contiguous.
void free_5d_array_with_4d_contiguous(void *****array, size_t ns, size_t nt);

void **alloc_2d_array(size_t nx, size_t ny, size_t size);

void free_2d_array(void **array);

// Allocates memory for a 3D array.
void ***alloc_3d_array(size_t nx, size_t ny, size_t nz, size_t size);

// Frees memory allocated by alloc_3d_array.
void free_3d_array(void ***array);

// Allocates memory for a 4D array.
void ****alloc_4d_array(size_t nt, size_t nx, size_t ny, size_t nz,
                        size_t size);

// Frees memory allocated by alloc_4d_array.
void free_4d_array(void ****array, size_t nt);

// Allocates memory for a 5D array.
void *****alloc_5d_array(size_t ns, size_t nt, size_t nx, size_t ny, size_t nz,
                         size_t size);

// Frees memory allocated by alloc_5d_array.
void free_5d_array(void *****array, size_t ns, size_t nt);

#endif
