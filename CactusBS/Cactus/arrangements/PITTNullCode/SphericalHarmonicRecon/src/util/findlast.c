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
    fprintf (stderr, "%s (%d): HDF5 call '%s' returned error code %d\n",
             __FILE__, __LINE__, fn_call, (int)error_code);
    exit(1);
  }
  return error_code;
}

static void print_usage(const char *name)
{
  fprintf(stderr, "usage: %s hdf5file\n"
                  "  where hdf5file contains the cauchy data\n",
                  name);
  exit(1);
}

int main(int argc, char **argv)
{

  hid_t       file_id;
  herr_t      status;
  char name[1024];
  if (argc!=2)
  {
    print_usage(argv[0]);
  }
  const char *infile = argv[1];
  util_verify_hdf5_file(infile);

  file_id = HDF5_ERROR(H5Fopen(infile, H5F_ACC_RDONLY, H5P_DEFAULT));

  int delta = 1000;
  int last = 0;
  double time;

  while (delta)
  {
    snprintf(name, sizeof name, "/%d/", last+delta);
    H5E_BEGIN_TRY {
      status = H5LTget_attribute_double(file_id,name, "Time", &time);
    } H5E_END_TRY;
    if (status < 0)
    {
      delta /= 2;
      continue;
    }
    last += delta;
  };

  time = atof("NAN");
  snprintf(name, sizeof name, "/%d/", last);
  HDF5_ERROR(H5LTget_attribute_double(file_id,name, "Time", &time));
  printf("%d: %e\n", last, time);
  HDF5_ERROR(H5Fclose(file_id));
  return 0;
}
