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
  fprintf(stderr, "usage: %s nlev n l m comp file\n"
      "where nlev is the number of timelevels in the hdf5 file,\n"
      "  and n, l, m describe the mode,\n"
      "  comp is an integer (0-9) giving the component to output,\n"
      "  file is the hdf5 file.\n", name);
}

static void read_meta(hid_t file_id,
                   int *nn,
                   int *na,
                   double *Rin,
                   double *Rout,
                   int *s);

static void read_in_data_for_component(
                   const hid_t file_id,
                   double *rC,
                   double *iC,
                   const int nlev,
                   const int skip,
                   const int ncoefs,
                   const char *gfname);

int main(int argc, char **argv)
{
  if (argc != 7)
  {
    print_usage(argv[0]);
    return EXIT_FAILURE;
  }

  const int nlev  = atoi(argv[1]);
  const int n  = atoi(argv[2]);
  const int l  = atoi(argv[3]);
  const int m  = atoi(argv[4]);
  const int comp  = atoi(argv[5]);
  const char *name = argv[6];
  int err = EXIT_SUCCESS;

  const int skip = 1;
  hid_t       file_id;

  util_verify_hdf5_file(name);
  file_id = HDF5_ERROR(H5Fopen(name, H5F_ACC_RDONLY, H5P_DEFAULT));

  int nn;
  int na;
  double Rin;
  double Rout;
  int spin;

  read_meta(file_id, &nn, &na, &Rin, &Rout, &spin);


  fprintf(stderr, "Run Parameters\n"
                  "... nn = %d\n"
                  "... na = %d\n"
                  "... Rin = %e\n"
                  "... Rout = %e\n", nn, na, Rin, Rout);

  double *rC = calloc(nn*na*nlev, sizeof(double)); 
  double *iC = calloc(nn*na*nlev, sizeof(double)); 

  double *time = calloc(nlev, sizeof(double));

  const char *gfnames[] = {
                            "gxx",
                            "gxy",
                            "gxz",
                            "gyy",
                            "gyz",
                            "gzz",
                            "betax",
                            "betay",
                            "betaz",
                            "alp"
        };

  int i;

  const int lmax = (int)(sqrt(na)+0.1) - 1; 
  
  assert (rC);
  assert (iC);

  herr_t      status;

  if (n<0 || n >= nn)
  {
    fprintf (stderr, "n must be between 0 and %d\n", nn-1);
    err = EXIT_FAILURE;
    goto end;
  }

  if (l<0 || l > lmax)
  {
    fprintf(stderr, "l must be between 0 and %d\n",
      (int)lmax);
    err = EXIT_FAILURE;
    goto end;
  }

  if (m<-l || m > l)
  {
    fprintf(stderr, "m must be between -%d and %d\n", l, l);
    err = EXIT_FAILURE;
    goto end;
  }

  if (comp >= sizeof(gfnames)/sizeof(*gfnames) ||
      comp < 0)
  {
    fprintf(stderr, "comp parameter must be between 0 and %d\n",
        (int)(sizeof(gfnames)/sizeof(*gfnames) -1));
    err = EXIT_FAILURE;
    goto end;
  }

  fprintf(stderr, "Writing n=%d, l=%d, m=%d of %s\n", n, l,m, gfnames[comp]);
  for (i=0; i < nlev; i++)
  {
    time[i] = sqrt(-1);
    char name[1024];
    snprintf(name, sizeof name, "/%d/", i*skip);
    status = HDF5_ERROR(H5LTget_attribute_double(file_id,name, "Time", time+i));

    if (status)
    {
      fprintf(stderr, "problem with time at %d\n", i);
    }
  }


  fprintf(stderr, "Reading in %s\n", gfnames[comp]);
  read_in_data_for_component(file_id, rC, iC,
             nlev, skip, nn*na, gfnames[comp]);
  fprintf(stderr, "done Reading in %s\n", gfnames[comp]);

  for (i=0; i < nlev; i++)
  {
#define MODE (l*l+l+m + n*na)
    printf("%20.16e %20.16e %20.16e\n",
          time[i], rC[i*nn*na+MODE], iC[i*nn*na+MODE]);
  }

  end:
  /* Cleanup */
  HDF5_ERROR(H5Fclose(file_id));

  free(rC);
  free(iC);
  free(time);

  return err;
}


static void read_in_data_for_component(
                   const hid_t file_id,
                   double *rC,
                   double *iC,
                   const int nlev,
                   const int skip,
                   const int ncoefs,
                   const char *gfname)
{
  int i;
  herr_t      status;

  for (i=0; i < nlev; i++)
  {
    char name[1024];

    snprintf(name, sizeof(name), "/%d/%s/re", i*skip, gfname);
    status = HDF5_ERROR(H5LTread_dataset_double(file_id, name, rC+i*ncoefs));

    if (status)
    {
      fprintf(stderr, "problem with %s\n", name);
    }

    snprintf(name, sizeof(name), "/%d/%s/im", i, gfname);
    status = HDF5_ERROR(H5LTread_dataset_double(file_id, name, iC+i*ncoefs));
    if (status)
    {
      fprintf(stderr, "problem with %s\n", name);
    }
  }
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
