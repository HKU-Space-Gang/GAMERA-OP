#include <H5Tpublic.h>
#include <hdf5.h>
#include <string.h>

#include "log.h"

#include "config.h"
#include "utils.h"
#include "curvilinear.h"

hid_t create_compressed_dcpl(int rank, const hsize_t *dims, int level) {
  if (rank < 1 || rank > 5) {
    log_error("Error: Unsupported rank %d for compression", rank);
    return -1;
  }

  hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
  if (dcpl_id < 0) return -1;

  hsize_t chunk_dims[5];  // 最多5维
  for (int i = 0; i < rank; ++i) {
    chunk_dims[i] = dims[i] < 64 ? dims[i] : 64;
  }

  if (H5Pset_chunk(dcpl_id, rank, chunk_dims) < 0) {
    H5Pclose(dcpl_id);
    return -1;
  }

  if (H5Pset_deflate(dcpl_id, level) < 0) {
    H5Pclose(dcpl_id);
    return -1;
  }

  return dcpl_id;
}

// write simple 1d data to hdf5 file
// return 0 if successful, -1 if failed
static herr_t write_simple_data_1d(hid_t file_id, const char *name,
                                   hid_t type_id, const void *data,
                                   hsize_t size) {
  herr_t status = 0;
  hsize_t dims[1] = {size};
  hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
  hid_t dcpl_id = create_compressed_dcpl(1, dims, 6);
//  hid_t dataset_id = H5Dcreate2(file_id, name, type_id, dataspace_id,
//                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t dataset_id = H5Dcreate2(file_id, name, type_id, dataspace_id,
                                H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  if (dataset_id < 0) {
    log_error("Error: Unable to create dataset %s", name);
    return -1;
  }
  status = H5Dwrite(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  if (status < 0) {
    log_error("Error: Unable to write dataset %s", name);
    H5Pclose(dcpl_id);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    return -1;
  }
  H5Pclose(dcpl_id);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  return 0;
}

static herr_t write_multiple_data_3d(hid_t file_id, const char *data_names[],
                                     hid_t type_id, double ***data[], int NI,
                                     int NJ, int NK, int num_datasets) {
  herr_t ret = 0;
  herr_t status = 0;
  hid_t dataset_id, dataspace_id, memspace_id;
  hsize_t mem_dims[3] = {NI, NJ, NK};

  // Create memory space and dataspace once, as they are the same for all
  // datasets
  memspace_id = H5Screate_simple(3, mem_dims, NULL);
  dataspace_id = H5Screate_simple(3, mem_dims, NULL);
  hid_t dcpl_id = create_compressed_dcpl(3, mem_dims, 6);

  for (int i = 0; i < num_datasets; i++) {
    log_info("Writing %s data", data_names[i]);
//    dataset_id = H5Dcreate2(file_id, data_names[i], type_id, dataspace_id,
//                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    dataset_id = H5Dcreate2(file_id, data_names[i], type_id, dataspace_id,
                                  H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    if (dataset_id < 0) {
      log_error("Error: Unable to create dataset for %s", data_names[i]);
      ret = -1;
      break;
    }

    status = H5Dwrite(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      &data[i][0][0][0]);
    if (status < 0) {
      log_error("Error: Unable to write data for %s", data_names[i]);
      ret = -1;
      H5Pclose(dcpl_id);
      H5Dclose(dataset_id);
      break;
    }

    H5Dclose(dataset_id);
  }

  // Close the dataspace and memory space after the loop
  H5Pclose(dcpl_id);
  H5Sclose(dataspace_id);
  H5Sclose(memspace_id);
  return ret;
}

static herr_t write_data_4d(hid_t file_id, const char *data_name, hid_t type_id,
                            double ****data, int NF, int NI, int NJ, int NK) {
  log_info("Writing %s data", data_name);
  herr_t ret = 0;
  herr_t status = 0;
  hid_t dataset_id, dataspace_id, memspace_id;
  log_info("Creating memory space");
  // Create memory space for one 4D field
  hsize_t mem_dims[4] = {NF, NI, NJ, NK};
  memspace_id = H5Screate_simple(4, mem_dims, NULL);

  log_info("Saving %s data", data_name);
  // Create the dataspace and dataset for field data
  hsize_t dims_4d[4] = {NF, NI, NJ, NK};
  dataspace_id = H5Screate_simple(4, dims_4d, NULL);
  hid_t dcpl_id = create_compressed_dcpl(4, mem_dims, 6);
//  dataset_id = H5Dcreate2(file_id, data_name, type_id, dataspace_id,
//                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dataset_id = H5Dcreate2(file_id, data_name, type_id, dataspace_id,
                          H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  if (dataset_id < 0) {
    log_error("Error: Unable to create dataset for %s", data_name);
    H5Pclose(dcpl_id);
    H5Sclose(dataspace_id);
    return -1;
  }

  hsize_t count[4] = {NF, NI, NJ, NK};
  hsize_t start[4] = {0, 0, 0, 0};
  status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count,
                               NULL);
  if (status < 0) {
    log_error("Error: Unable to select hyperslab for %s", data_name);
    ret = -1;
  }
  status = H5Dwrite(dataset_id, type_id, memspace_id, dataspace_id, H5P_DEFAULT,
                    &data[0][0][0][0]);
  if (status < 0) {
    log_error("Error: Unable to write data for %s", data_name);
    ret = -1;
  }
  log_info("H5Dwrite for field %s done", data_name);
  // Close and release resources
  H5Pclose(dcpl_id);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  return ret;
}

static herr_t write_data_5d(hid_t file_id, const char *data_name, hid_t type_id,
                            double *****data, int NF, int NI, int NJ, int NK) {
  log_info("Writing %s data", data_name);
  herr_t ret = 0;
  herr_t status = 0;
  hid_t dataset_id, dataspace_id, memspace_id;
  log_info("Creating memory space for %s", data_name);
  // Create memory space for one 4D field
  hsize_t mem_dims[4] = {NF, NI, NJ, NK};
  memspace_id = H5Screate_simple(4, mem_dims, NULL);

  log_info("Saving %s data", data_name);
  // Create the dataspace and dataset for field data
  hsize_t dims_5d[5] = {NS1, NF, NI, NJ, NK};
  dataspace_id = H5Screate_simple(5, dims_5d, NULL);
  hid_t dcpl_id = create_compressed_dcpl(5, dims_5d, 6);
//  dataset_id = H5Dcreate2(file_id, data_name, type_id, dataspace_id,
//                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dataset_id = H5Dcreate2(file_id, data_name, type_id, dataspace_id,
                          H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  if (dataset_id < 0) {
    log_error("Error: Unable to create dataset for %s", data_name);
    H5Pclose(dcpl_id);
    H5Sclose(dataspace_id);
    H5Sclose(memspace_id);
    return -1;
  }

  // Write data species by species
  hsize_t count_5d[5] = {1, NF, NI, NJ, NK};
  hsize_t start[5] = {0, 0, 0, 0, 0};
  for (int s = 0; s < NS1; s++) {
    start[0] = s;
    status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL,
                                 count_5d, NULL);
    if (status < 0) {
      log_error("Error: Unable to select hyperslab for %s[%d]", data_name, s);
      ret = -1;
      break;
    }
    status = H5Dwrite(dataset_id, type_id, memspace_id, dataspace_id,
                      H5P_DEFAULT, &data[s][0][0][0][0]);
    if (status < 0) {
      log_error("Error: Unable to write data for %s[%d]", data_name, s);
      ret = -1;
      break;
    }
  }
  // Close and release resources
  H5Pclose(dcpl_id);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  H5Sclose(memspace_id);
  return ret;
}

void get_out_filename(const char *base_filename, char *filename,
                      size_t filename_size, int seq_num,
                      const char *extension) {
  // Create filename with timestamp and coordinates
  //  snprintf(filename, sizeof(filename),
  //           "%s_%03d-%03d-%03d_%02d%02d%02d-%02d%02d%02d.h5", base_filename,
  //           proc_coords[0], proc_coords[1], proc_coords[2],
  //           timeinfo->tm_year - 100, timeinfo->tm_mon + 1, timeinfo->tm_mday,
  //           timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);

  snprintf(filename, filename_size, out_file_pattern, base_filename,
           proc_coords[0], proc_coords[1], proc_coords[2], seq_num, extension);
}


int dump_extra_data(const char *base_filename) {
  int status = 0;
  hid_t file_id;
  char filename[256];
  get_out_filename(base_filename, filename, sizeof(filename), hdf_seq_num,
                   "h5");
  log_info("Dumping extra hdf5 file: %s", filename);
  // Create a new HDF5 file
//  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
  if (file_id < 0) {
    log_error("Error: Unable to create HDF5 file %s", filename);
    return -1;
  }
  // Write geo data
  // --------------------------------------------------------------
  status = write_data_4d(file_id, "/geo", H5T_NATIVE_DOUBLE, geo, NF_geo,
                         config.NI, config.NJ, config.NK);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }

  if (doRingAverage == 1) {
    status = write_data_4d(file_id, "/geo_ring", H5T_NATIVE_DOUBLE, geo_ring, NF_geo_ring,
                           Ring_Ni, Ring_Nj, Ring_Nk);
    if (status < 0) {
      H5Fclose(file_id);
      return -1;
    }
  }

  // write rec data
  status = write_data_4d(file_id, "/rec", H5T_NATIVE_DOUBLE, rec, NF_rec,
                           config.NI, config.NJ, config.NK);
  if (status < 0) {
      H5Fclose(file_id);
      return -1;
  }

//  //write partial gem data
//  status = write_data_4d(file_id, "/gem_partial", H5T_NATIVE_DOUBLE, gem, NF_gem,
//                         config.NI, config.NJ, config.NK);
//  if (status < 0) {
//    H5Fclose(file_id);
//    return -1;
//  }

  //write partial gem data
  const char *data_names_background_field[] = {"/B0_i", "/B0_j", "/B0_k"};
  double ***data_background_field[] = {gem[B0AfaceI_b1], gem[B0AfaceJ_b2], gem[B0AfaceK_b3]};
  int num_datasets_background_field = 3;

  status = write_multiple_data_3d(file_id, data_names_background_field, H5T_NATIVE_DOUBLE,
                                  data_background_field, config.NI, config.NJ, config.NK,
                                  num_datasets_background_field);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }

  const char *data_names_background_field_B0[] = {"/B0_1", "/B0_2", "/B0_3"};
  double ***data_background_field_B0[] = {gem[B0_b1], gem[B0_b2], gem[B0_b3]};
  int num_datasets_background_field_B0 = 3;

  status = write_multiple_data_3d(file_id, data_names_background_field_B0, H5T_NATIVE_DOUBLE,
                                  data_background_field_B0, config.NI, config.NJ, config.NK,
                                  num_datasets_background_field_B0);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }

  // write center grid data
  const char *data_names_center[] = {"/x1c", "/x2c", "/x3c"};
  double ***data_center[] = {x1c, x2c, x3c};
  int num_datasets_center = 3;

  status = write_multiple_data_3d(file_id, data_names_center, H5T_NATIVE_DOUBLE,
                                  data_center, config.NI, config.NJ, config.NK,
                                  num_datasets_center);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }

  // write ctr grid data
  const char *data_names_ctr[] = {"/x1ctr", "/x2ctr", "/x3ctr"};
  double ***data_ctr[] = {x1ctr, x2ctr, x3ctr};
  int num_datasets_ctr = 3;

  status = write_multiple_data_3d(file_id, data_names_ctr, H5T_NATIVE_DOUBLE,
                                  data_ctr, config.NI, config.NJ, config.NK,
                                  num_datasets_ctr);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }

  // dx
  const char *data_names_dx[] = {"/dx1", "/dx2", "/dx3"};
  double ***data_dx[] = {dx1, dx2, dx3};
  int num_datasets_dx = 3;

  status = write_multiple_data_3d(file_id, data_names_dx, H5T_NATIVE_DOUBLE,
                                  data_dx, config.NI, config.NJ, config.NK,
                                  num_datasets_dx);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }

  H5Fclose(file_id);
  return 0;
}

// Dump data to HDF5 file
// return 0 if successful, -1 if failed
int dump_data_hdf5(const char *base_filename) {
  int ret = 0;
  char filename[256];
  hid_t file_id;
  herr_t status = 0;
  hdf_seq_num++;
  get_out_filename(base_filename, filename, sizeof(filename), hdf_seq_num,
                   "h5");
  log_info("Dumping hdf5 file: %s", filename);
  // Create a new HDF5 file
  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file_id < 0) {
    log_error("Error: Unable to create HDF5 file %s", filename);
    return -1;
  }

  // Save time as a dataset
  status = write_simple_data_1d(file_id, "/time_sim", H5T_NATIVE_DOUBLE,
                                &time_sim, 1);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  // Save log_seq_num
  status = write_simple_data_1d(file_id, "/log_seq_num", H5T_NATIVE_DOUBLE,
                                &log_seq_num, 1);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  // Save dt0
  status = write_simple_data_1d(file_id, "/dt0", H5T_NATIVE_DOUBLE,
                                &dt0, 1);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }

  // Create and write the config dataset
  status = write_simple_data_1d(file_id, "/config", H5T_NATIVE_INT, &config,
                                sizeof(config_t) / sizeof(int));
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  // Write the problem dataset
  status = write_simple_data_1d(
      file_id, "/problem_config", H5T_NATIVE_INT, &problem_config,
      sizeof(Problem_config_t) / sizeof(Boundary_condition_t));
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
//  // Write planet config
//  status = write_simple_data_1d(file_id, "/planet_config", H5T_NATIVE_INT,
//                                &planet_config, sizeof(Planet_config_t) / sizeof(double));
//  if (status < 0) {
//    H5Fclose(file_id);
//    return -1;
//  }
  // Create and write the norm_config dataset
  double *ptr = (double *)&norm_config;
  status = write_simple_data_1d(file_id, "/norm_config", H5T_NATIVE_DOUBLE, ptr,
                                sizeof(Norm_t) / sizeof(double));
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  // Create and write the proc_coords dataset
  status = write_simple_data_1d(file_id, "/proc_coords", H5T_NATIVE_INT,
                                proc_coords, 3);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  // write grid data
  status = write_simple_data_1d(file_id, "/x1min_global", H5T_NATIVE_DOUBLE,
                                &x1min_global, 1);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  status = write_simple_data_1d(file_id, "/x1max_global", H5T_NATIVE_DOUBLE,
                                &x1max_global, 1);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  status = write_simple_data_1d(file_id, "/x2min_global", H5T_NATIVE_DOUBLE,
                                &x2min_global, 1);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  status = write_simple_data_1d(file_id, "/x2max_global", H5T_NATIVE_DOUBLE,
                                &x2max_global, 1);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  status = write_simple_data_1d(file_id, "/x3min_global", H5T_NATIVE_DOUBLE,
                                &x3min_global, 1);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  status = write_simple_data_1d(file_id, "/x3max_global", H5T_NATIVE_DOUBLE,
                                &x3max_global, 1);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }

  const char *data_names[] = {"/x1", "/x2", "/x3"};
  double ***data[] = {x1, x2, x3};
  int num_datasets = 3;

  status =
      write_multiple_data_3d(file_id, data_names, H5T_NATIVE_DOUBLE, data,
                             config.NI, config.NJ, config.NK, num_datasets);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  // Write gem data
  // --------------------------------------------------------------
  status = write_data_4d(file_id, "/gem", H5T_NATIVE_DOUBLE, gem, NF_gem_store,
                         config.NI, config.NJ, config.NK);
  if (status < 0) {
    log_error("Error: Unable to write gem data to HDF5 file %s with code %d",
              filename, status);
    log_info("gem: %p, mag_bk: %d\n", gem, mag_bk);
    H5Fclose(file_id);
    return -1;
  }

//  status = write_data_4d(file_id, "/gem_ring", H5T_NATIVE_DOUBLE, gem_ring, NF_gem_ring,
//                            Ring_Ni, Ring_Nj, Ring_Nk);
//  if (status < 0) {
//    log_error("Error: Unable to write gem data to HDF5 file %s with code %d",
//              filename, status);
//    log_info("gem: %p, mag_bk: %d\n", gem, mag_bk);
//    H5Fclose(file_id);
//    return -1;
//  }

  // Write gas data
  status = write_data_5d(file_id, "/gas", H5T_NATIVE_DOUBLE, gas, NF_gas_store,
                         config.NI, config.NJ, config.NK);
  if (status < 0) {
    log_error("Error: Unable to write gas data to HDF5 file %s with code %d",
              filename, status);
    log_info("gas: %p, gas_p: %d\n", gas, gas_p);
    H5Fclose(file_id);
    return -1;
  }

//  status = write_data_5d(file_id, "/gas_ring", H5T_NATIVE_DOUBLE, gas_ring, NF_gas_ring,
//                         Ring_Ni, Ring_Nj, Ring_Nk);
//  if (status < 0) {
//    log_error("Error: Unable to write gas data to HDF5 file %s with code %d",
//              filename, status);
//    log_info("gas: %p, gas_p: %d\n", gas, gas_p);
//    H5Fclose(file_id);
//    return -1;
//  }

  // Close and release resources
  H5Fclose(file_id);
  return ret;
}

// read simple 1d data from hdf5 file
// return 0 if successful, -1 if failed
int read_simple_data_1d(hid_t file_id, const char *name, hid_t type_id,
                        void *data) {
  herr_t status = 0;
  hid_t dataset_id = H5Dopen2(file_id, name, H5P_DEFAULT);
  if (dataset_id < 0) {
    log_error("Error: Unable to open dataset %s", name);
    return -1;
  }

  status = H5Dread(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  if (status < 0) {
    log_error("Error: Unable to read dataset %s", name);
    H5Dclose(dataset_id);
    return -1;
  }

  H5Dclose(dataset_id);
  return 0;
}

static herr_t read_multiple_data_3d(hid_t file_id, const char *data_names[],
                                    hid_t type_id, double ***data[], int NI,
                                    int NJ, int NK, int num_datasets) {
  herr_t ret = 0;
  herr_t status = 0;
  hid_t dataset_id, dataspace_id, memspace_id;
  hsize_t mem_dims[3] = {NI, NJ, NK};

  // Create memory space once, as it is the same for all datasets
  memspace_id = H5Screate_simple(3, mem_dims, NULL);

  for (int i = 0; i < num_datasets; i++) {
    log_info("Reading %s data", data_names[i]);
    dataset_id = H5Dopen2(file_id, data_names[i], H5P_DEFAULT);
    if (dataset_id < 0) {
      log_error("Error: Unable to open dataset %s", data_names[i]);
      ret = -1;
      break;
    }

    // Get the dataspace for the current dataset
    dataspace_id = H5Dget_space(dataset_id);
    if (dataspace_id < 0) {
      log_error("Error: Unable to get dataspace for %s", data_names[i]);
      H5Dclose(dataset_id);
      ret = -1;
      break;
    }

    status = H5Dread(dataset_id, type_id, memspace_id, dataspace_id,
                     H5P_DEFAULT, &data[i][0][0][0]);
    if (status < 0) {
      log_error("Error: Unable to read data for %s", data_names[i]);
      ret = -1;
      H5Dclose(dataset_id);
      H5Sclose(dataspace_id);
      break;
    }

    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
  }

  // Close the memory space after the loop
  H5Sclose(memspace_id);
  return ret;
}

static herr_t read_data_4d(hid_t file_id, const char *data_name, hid_t type_id,
                           double ****data, int NF, int NI, int NJ, int NK) {
  log_info("Reading %s data", data_name);
  herr_t status = 0;
  hid_t dataset_id, dataspace_id, memspace_id;
  // Create memory space for one 4D field
  hsize_t mem_dims[4] = {NF, NI, NJ, NK};
  memspace_id = H5Screate_simple(4, mem_dims, NULL);
  // read gem data -------------------------------------------------------------
  dataset_id = H5Dopen2(file_id, data_name, H5P_DEFAULT);
  if (dataset_id < 0) {
    log_error("Error: Unable to open %s dataset", data_name);
    return -1;
  }
  dataspace_id = H5Dget_space(dataset_id);
  hsize_t count_4d[4] = {NF, NI, NJ, NK};
  hsize_t onface_4d[4] = {0, 0, 0, 0};
  status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, onface_4d, NULL,
                               count_4d, NULL);
  if (status < 0) {
    log_error("Error: Unable to select hyperslab for %s", data_name);
    return -1;
  }
  status = H5Dread(dataset_id, type_id, memspace_id, dataspace_id, H5P_DEFAULT,
                   &data[0][0][0][0]);
  if (status < 0) {
    log_error("Error: Unable to read data for %s", data_name);
    return -1;
  }
  H5Sclose(dataspace_id);
  H5Dclose(dataset_id);
  return 0;
}

static herr_t read_data_5d(hid_t file_id, const char *data_name, hid_t type_id,
                           double *****data, int NF, int NI, int NJ, int NK) {
  log_info("Reading %s data", data_name);
  herr_t status = 0;
  hid_t dataset_id, dataspace_id, memspace_id;
  // Create memory space for one 3D field
  hsize_t mem_dims[4] = {NF, NI, NJ, NK};
  memspace_id = H5Screate_simple(4, mem_dims, NULL);

  dataset_id = H5Dopen2(file_id, data_name, H5P_DEFAULT);
  if (dataset_id < 0) {
    log_error("Error: Unable to open %s dataset", data_name);
    H5Sclose(memspace_id);
    return -1;
  }
  dataspace_id = H5Dget_space(dataset_id);

  hsize_t count_5d[5] = {1, NF, NI, NJ, NK};
  hsize_t start[5] = {0, 0, 0, 0, 0};
  for (int s = 0; s < NS1; s++) {
    start[0] = s;
    status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL,
                                 count_5d, NULL);
    if (status < 0) {
      log_error("Error: Unable to select hyperslab for %s[%d]", data_name, s);
      H5Sclose(memspace_id);
      H5Sclose(dataspace_id);
      H5Dclose(dataset_id);
      return -1;
    }
    status = H5Dread(dataset_id, type_id, memspace_id, dataspace_id,
                     H5P_DEFAULT, &data[s][0][0][0][0]);
    if (status < 0) {
      log_error("Error: Unable to read data for %s[%d]", data_name, s);
      H5Sclose(memspace_id);
      H5Sclose(dataspace_id);
      H5Dclose(dataset_id);
      return -1;
    }
  }

  // Close resources
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset_id);
  return 0;
}

// Read data from HDF5 file
// return 0 if successful, -1 if failed
int read_data_hdf5(const char *filename) {
  log_info("Reading hdf5 file: %s", filename);

  hid_t file_id;
  herr_t status = 0;

  // Open the HDF5 file
  file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) {
    log_error("Error: Unable to open HDF5 file %s", filename);
    return -1;
  }
  // read time from dataset
  status =
      read_simple_data_1d(file_id, "/time_sim", H5T_NATIVE_DOUBLE, &time_sim);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  log_info("Read time: %.15f", time_sim);
  // read log_seq_num from dataset
  status =
      read_simple_data_1d(file_id, "/log_seq_num", H5T_NATIVE_DOUBLE, &log_seq_num);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  log_info("log_seq_num: %d", log_seq_num);
  // read dt0 from dataset
  status =
      read_simple_data_1d(file_id, "/dt0", H5T_NATIVE_DOUBLE, &dt0);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  log_info("dt0: %.15f", dt0);
  // Read the config data
  config_t config_read;
  status =
      read_simple_data_1d(file_id, "/config", H5T_NATIVE_INT, &config_read);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  // check if the config is the same as the one in the file
  if (memcmp(&config_read, &config, sizeof(config_t)) != 0) {
    log_error("Config is not the same as the one in the file.");
    H5Fclose(file_id);
    return -1;
  }

  // read problem config
  status = read_simple_data_1d(file_id, "/problem_config", H5T_NATIVE_INT,
                               &problem_config);
  // check if the config is the same as the one in the file
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }

  // read norm_config from the file
  Norm_t norm_config_read;
  double *ptr = (double *)&norm_config_read;
  status = read_simple_data_1d(file_id, "/norm_config", H5T_NATIVE_DOUBLE, ptr);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  if (memcmp(&norm_config_read, &norm_config, sizeof(Norm_t)) != 0) {
    log_error("Norm config is not the same as the one in the file.");
    H5Fclose(file_id);
    return -1;
  }

  // read grid data
  // -------------------------------------------------------------
  status =
      read_simple_data_1d(file_id, "/x1min_global", H5T_NATIVE_DOUBLE, &x1min_global);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  log_info("x1min_global: %.15f", x1min_global);
  status =
      read_simple_data_1d(file_id, "/x1max_global", H5T_NATIVE_DOUBLE, &x1max_global);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  log_info("x1max_global: %.15f", x1max_global);
  status =
      read_simple_data_1d(file_id, "/x2min_global", H5T_NATIVE_DOUBLE, &x2min_global);
  if (status < 0) {
      H5Fclose(file_id);
      return -1;
  }
  log_info("x2min_global: %.15f", x2min_global);
  status =
      read_simple_data_1d(file_id, "/x2max_global", H5T_NATIVE_DOUBLE, &x2max_global);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  log_info("x2max_global: %.15f", x2max_global);
  status =
      read_simple_data_1d(file_id, "/x3min_global", H5T_NATIVE_DOUBLE, &x3min_global);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  log_info("x3min_global: %.15f", x3min_global);
  status =
      read_simple_data_1d(file_id, "/x3max_global", H5T_NATIVE_DOUBLE, &x3max_global);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  log_info("x3max_global: %.15f", x3max_global);

  const char *data_names[] = {"/x1", "/x2", "/x3"};
  double ***data[] = {x1, x2, x3};
  int num_datasets = 3;

  status = read_multiple_data_3d(file_id, data_names, H5T_NATIVE_DOUBLE, data,
                                 config.NI, config.NJ, config.NK, num_datasets);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  // read gem data -------------------------------------------------------------
  status = read_data_4d(file_id, "/gem", H5T_NATIVE_DOUBLE, gem, NF_gem_store,
                        config.NI, config.NJ, config.NK);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  // Read gas data
  status = read_data_5d(file_id, "/gas", H5T_NATIVE_DOUBLE, gas, NF_gas_store,
                        config.NI, config.NJ, config.NK);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }

  // Close resources
  H5Fclose(file_id);
  return 0;
}

int read_data_dims(char *filename, char *dataset_name, int *dims, int rank) {

  if (!filename || !dataset_name || !dims || rank <= 0) {
    fprintf(stderr, "Invalid parameters to read_data_dims\n");
    return -1;
  }

  hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) {
    fprintf(stderr, "Cannot open HDF5 file: %s\n", filename);
    return -1;
  }

  hid_t dataset_id = H5Dopen(file_id, dataset_name, H5P_DEFAULT);
  if (dataset_id < 0) {
    fprintf(stderr, "Cannot open dataset: %s\n", dataset_name);
    H5Fclose(file_id);
    return -1;
  }

  hid_t dataspace_id = H5Dget_space(dataset_id);
  if (dataspace_id < 0) {
    fprintf(stderr, "Cannot get dataspace for dataset: %s\n", dataset_name);
    H5Dclose(dataset_id);
    H5Fclose(file_id);
    return -1;
  }

  int actual_rank = H5Sget_simple_extent_ndims(dataspace_id);
  if (actual_rank < 0) {
    fprintf(stderr, "Cannot get rank for dataset: %s\n", dataset_name);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    H5Fclose(file_id);
    return -1;
  }

  if (actual_rank != rank) {
    fprintf(stderr, "Dimension mismatch for dataset %s: expected rank=%d, actual rank=%d\n",
            dataset_name, rank, actual_rank);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    H5Fclose(file_id);
    return -1;
  }

  hsize_t h5_dims[H5S_MAX_RANK];
  int status = H5Sget_simple_extent_dims(dataspace_id, h5_dims, NULL);
  if (status < 0) {
    fprintf(stderr, "Cannot get dimensions for dataset: %s\n", dataset_name);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    H5Fclose(file_id);
    return -1;
  }

  for (int i = 0; i < rank; ++i) {
    dims[i] = (int)h5_dims[i];
  }

  H5Sclose(dataspace_id);
  H5Dclose(dataset_id);
  H5Fclose(file_id);

  return 0;
}

int read_grid_data_hdf5(const char *filename, double ***x1_temp,
                        double ***x2_temp, double ***x3_temp) {
  log_info("Reading hdf5 file: %s", filename);

  hid_t file_id;
  herr_t status;

  // Open the HDF5 file
  file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) {
    log_error("Error: Unable to open HDF5 file %s", filename);
    return -1;
  }

  // -------------------------------------------------------------
  const char *data_names[] = {"/x1grid", "/x2grid", "/x3grid"};
  double ***data[] = {x1_temp, x2_temp, x3_temp};
  int num_datasets = 3;

  status = read_multiple_data_3d(file_id, data_names, H5T_NATIVE_DOUBLE, data,
                                 config.nk_global + 2 * NG + 1,
                                 config.nj_global + 2 * NG + 1,
                                 config.ni_global + 2 * NG + 1, num_datasets);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }

  // Close resources
  H5Fclose(file_id);
  return 0;
}

int read_wind_data_hdf5(const char *filename, double *wind_time, double *wind_rho, double *wind_vx, double *wind_vy, double *wind_vz, double *wind_p,
                        double *wind_bx, double *wind_by, double *wind_bz) {
  log_info("Reading wind hdf5 file: %s", filename);

  hid_t file_id;
  herr_t status;

  // Open the HDF5 file
  file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) {
    log_error("Error: Unable to open HDF5 file %s", filename);
    return -1;
  }

  // -------------------------------------------------------------
  // read Time
  status = read_simple_data_1d(file_id, "/T", H5T_NATIVE_DOUBLE, wind_time);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  // read Gas
  status = read_simple_data_1d(file_id, "/D", H5T_NATIVE_DOUBLE, wind_rho);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  status = read_simple_data_1d(file_id, "/Vx", H5T_NATIVE_DOUBLE, wind_vx);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  status = read_simple_data_1d(file_id, "/Vy", H5T_NATIVE_DOUBLE, wind_vy);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  status = read_simple_data_1d(file_id, "/Vz", H5T_NATIVE_DOUBLE, wind_vz);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  status = read_simple_data_1d(file_id, "/P", H5T_NATIVE_DOUBLE, wind_p);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  // read B field
  status = read_simple_data_1d(file_id, "/Bx", H5T_NATIVE_DOUBLE, wind_bx);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  status = read_simple_data_1d(file_id, "/By", H5T_NATIVE_DOUBLE, wind_by);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }
  status = read_simple_data_1d(file_id, "/Bz", H5T_NATIVE_DOUBLE, wind_bz);
  if (status < 0) {
    H5Fclose(file_id);
    return -1;
  }

  // Close resources
  H5Fclose(file_id);
  return 0;
}
