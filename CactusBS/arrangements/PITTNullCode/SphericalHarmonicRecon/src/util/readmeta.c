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
  fprintf(stderr, "usage: %s file\n"
                  " where file is the hdf5 file with the cauchy data\n",
                   name);
  exit(-1);
}


static void read_meta(hid_t file_id,
                   int *nn,
                   int *na,
                   double *Rin,
                   double *Rout,
                   int *s)
{

  int dim[2];

  HDF5_ERROR(H5LTget_attribute_int(file_id,"/metadata/","spin", s));
  HDF5_ERROR(H5LTget_attribute_int(file_id,"/metadata/","dim", dim));
  HDF5_ERROR(H5LTget_attribute_double(file_id,"/metadata/","Rin", Rin));
  HDF5_ERROR(H5LTget_attribute_double(file_id,"/metadata/","Rout", Rout));

  *nn = dim[0];
  *na = dim[1];
}

int main(int argc, char **argv)
{
  int nn;
  int na;
  double Rin;
  double Rout;
  int spin;
  if (argc != 2)
  {
    print_usage(argv[0]);
  }

  const char *file = argv[1];

  util_verify_hdf5_file(file);

  hid_t       file_id;

  file_id = HDF5_ERROR(H5Fopen(file, H5F_ACC_RDONLY, H5P_DEFAULT));

  read_meta(file_id, &nn, &na, &Rin, &Rout, &spin);

  fprintf(stderr, "Run Parameters\n"
                  "... nn = %d\n"
                  "... na = %d\n"
                  "... Rin = %e\n"
                  "... Rout = %e\n", nn, na, Rin, Rout);

  /* Cleanup */
  HDF5_ERROR(H5Fclose(file_id));

  return EXIT_SUCCESS;
}
