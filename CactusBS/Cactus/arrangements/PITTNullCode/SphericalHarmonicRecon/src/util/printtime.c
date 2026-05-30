#define H5_USE_16_API 1
#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "util.h"

#define pow2(x_) ((x_)*(x_))

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
  fprintf(stderr, "usage: %s it file\n"
                  " where it is the iteration number"
                  " and file is the hdf5 file with the cauchy data\n",
                   name);
  exit(-1);
}

int main(int argc, char **argv)
{
  if (argc != 3)
  {
    print_usage(argv[0]);
  }

  const int it = atoi(argv[1]);
  const char *file = argv[2];

  assert(it >= 0);

  util_verify_hdf5_file(file);

  hid_t       file_id;

  file_id = HDF5_ERROR(H5Fopen(file, H5F_ACC_RDONLY, H5P_DEFAULT));

  double time;

  time = sqrt(-1);
  char name[1024];
  snprintf(name, sizeof name, "/%d/", it);
  HDF5_ERROR(H5LTget_attribute_double(file_id,name, "Time", &time));
  printf("%d: %e\n", it, time);
  HDF5_ERROR(H5Fclose(file_id));
  return 0;
}
