#define H5_USE_16_API 1
#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "util.h"

// check return code of HDF5 call and print a warning in case of an error
#define HDF5_ERROR(fn_call) HDF5_ERROR_FUNC((fn_call), #fn_call, \
                                            __LINE__, __FILE__)
static hid_t HDF5_ERROR_FUNC(hid_t error_code, const char* fn_call,
                           const int line, const char *file)
{
  if (error_code < 0)
  {
    fprintf (stderr, "%s (%d): HDF5 call '%s' returned error code %d",
             __FILE__, __LINE__, fn_call, (int)error_code);
  }
  return error_code;
}

static void print_usage(const char *name)
{
  fprintf(stderr, "usage: %s nn na Rin Rout file"
      "where nn is the number of radial modes,\n"
      "na is the number of angular modes,\n"
      "Rin and Rout are the inner and outer extraction radii\n"
      "and file is the hdf5 with the cauchy data", name);
}

static void write_meta(hid_t fileid,
                   int nn,
                   int na,
                   double Rin,
                   double Rout,
                   int s);

int main(int argc, char **argv)
{
  if (argc != 6)
  {
    print_usage(argv[0]);
    return EXIT_FAILURE;
  }

  const int nn = atoi(argv[1]);
  const int na = atoi(argv[2]);
  const double Rin = atof(argv[3]);
  const double Rout = atof(argv[4]);
  const char *file = argv[5];

  hid_t       file_id;

  util_verify_hdf5_file(file);

  file_id = HDF5_ERROR(H5Fopen(file, H5F_ACC_RDWR, H5P_DEFAULT));

  fprintf(stderr, "Run Parameters\n"
                  "... nn = %d\n"
                  "... na = %d\n"
                  "... Rin = %e\n"
                  "... Rout = %e\n", nn, na, Rin, Rout);

  write_meta(file_id, nn, na, Rin, Rout, 0);

  /* Cleanup */
  HDF5_ERROR(H5Fclose(file_id));

  return EXIT_SUCCESS;
}


static void write_meta(hid_t fileid,
                   int nn,
                   int na,
                   double Rin,
                   double Rout,
                   int s)
{
  hid_t attribute_id, dataspace_id, group_id;
  int ds[2] = {nn, na};
  hsize_t oD2 = 2;
  hsize_t oD1 = 1;

  char name[]="/metadata";
  group_id= HDF5_ERROR(H5Gcreate(fileid, name, 0));
  dataspace_id =  HDF5_ERROR(H5Screate_simple(1, &oD2, NULL));
  attribute_id = HDF5_ERROR(H5Acreate(group_id, "dim", H5T_NATIVE_INT,
                                      dataspace_id, H5P_DEFAULT));
  HDF5_ERROR(H5Awrite(attribute_id, H5T_NATIVE_INT, ds));
  HDF5_ERROR(H5Aclose(attribute_id));
  HDF5_ERROR(H5Sclose(dataspace_id));

  dataspace_id = HDF5_ERROR(H5Screate_simple(1, &oD1, NULL));
  attribute_id = HDF5_ERROR(H5Acreate(group_id, "spin", H5T_NATIVE_INT,
                                      dataspace_id, H5P_DEFAULT));
  HDF5_ERROR(H5Awrite(attribute_id, H5T_NATIVE_INT, &s));
  HDF5_ERROR(H5Aclose(attribute_id));
  HDF5_ERROR(H5Sclose(dataspace_id));

  dataspace_id = HDF5_ERROR(H5Screate_simple(1, &oD1, NULL));
  attribute_id = HDF5_ERROR(H5Acreate(group_id, "Rin", H5T_NATIVE_DOUBLE,
                                      dataspace_id, H5P_DEFAULT));
  HDF5_ERROR(H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &Rin));
  HDF5_ERROR(H5Aclose(attribute_id));
  HDF5_ERROR(H5Sclose(dataspace_id));

  dataspace_id = HDF5_ERROR(H5Screate_simple(1, &oD1, NULL));
  attribute_id = HDF5_ERROR(H5Acreate(group_id, "Rout", H5T_NATIVE_DOUBLE,
                                      dataspace_id, H5P_DEFAULT));
  HDF5_ERROR(H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &Rout));
  HDF5_ERROR(H5Aclose(attribute_id));
  HDF5_ERROR(H5Sclose(dataspace_id));

  HDF5_ERROR(H5Gclose(group_id));
}
