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
  fprintf(stderr, "usage: %s nlast Tscale file\n"
      "  where nlast is the last iteration,\n"
      "  Tscale is damping timescale (in physical units)\n"
      "  file is the hdf5 file with the cauchy data\n", name);
}

static void filter(const int n,
                   const double c, /* c = dt/Tscale */
                   const double *in,
                   double *out);

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

static void write_meta(hid_t fileid,
                   int nn,
                   int na,
                   double Rin,
                   double Rout,
                   int s);

static void write_out_data_for_component(
                   hid_t fileid,
                   double *rC,
                   double *iC,
                   int nlev,
                   int nn,
                   int na,
                   const char *name,
                   double *time);

int main(int argc, char **argv)
{
  if (argc != 4)
  {
    print_usage(argv[0]);
    return EXIT_FAILURE;
  }

  const int last = atoi(argv[1]);
  const double Tscale = atof(argv[2]);
  const char *infile = argv[3];
  const int skip = 1;
  const int nlev = last/skip + 1;
  hid_t       file_id;
  char *outfile = NULL;

  util_verify_hdf5_file(infile);

  util_verify_hdf5_file(infile);
  if (util_gen_out_hdf5_filename(infile, &outfile, "ft") != 0)
  {
    fprintf(stderr, "Error creating output file\n");
    exit(-1);
  }

  file_id = HDF5_ERROR(H5Fopen(infile, H5F_ACC_RDONLY, H5P_DEFAULT));

#ifdef OVERRIDE_PARAMETERS
  const int nn = 7;
  const int na = 49;
  const double Rin = 18.0;
  const double Rout = 22.0;
#else
  int nn;
  int na;
  double Rin;
  double Rout;
  int spin;

  read_meta(file_id, &nn, &na, &Rin, &Rout, &spin);

#endif

  fprintf(stderr, "Run Parameters\n"
                  "... nn = %d\n"
                  "... na = %d\n"
                  "... Rin = %e\n"
                  "... Rout = %e\n", nn, na, Rin, Rout);

  double *rC = calloc(nn*na*nlev, sizeof(double)); 
  double *iC = calloc(nn*na*nlev, sizeof(double)); 

  double *rCt = calloc(nn*na*nlev, sizeof(double)); 
  double *iCt = calloc(nn*na*nlev, sizeof(double)); 

  double *time = calloc(nlev, sizeof(double));

  double *in_r = calloc(nlev, sizeof(double));
  double *in_i = calloc(nlev, sizeof(double));

  double *out_r = calloc(nlev, sizeof(double));
  double *out_i = calloc(nlev, sizeof(double));


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

  int i, j, k;
  
  assert (rC);
  assert (iC);
  assert (rCt);
  assert (iCt);

  herr_t      status;
  hid_t fout_id;

  fout_id = HDF5_ERROR(H5Fcreate (outfile, H5F_ACC_TRUNC,
                       H5P_DEFAULT, H5P_DEFAULT));

  write_meta(fout_id, nn, na, Rin, Rout, 0);

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

  const double dt=time[1] - time[0];

  for (i=0; i < 10; i++)
  {
    fprintf(stderr, "Reading in %s\n", gfnames[i]);
    read_in_data_for_component(file_id, rC, iC,
               nlev, skip, nn*na, gfnames[i]);
    fprintf(stderr, "done Reading in %s\n", gfnames[i]);

    for (j=0; j < nn*na; j++)
    {
      fprintf(stderr, "transforming %d of %d \n", j, nn*na);
      for (k=0; k < nlev; k++)
      {
        in_r[k] = rC[k*(nn*na) + j];
        in_i[k] = iC[k*(nn*na) + j];
      }

      filter(nlev, dt/Tscale, in_r, out_r);
      filter(nlev, dt/Tscale, in_i, out_i);

      /* fftw transform is not normalized,
         so we need to renormalize here */
      for (k=0; k < nlev; k++)
      {
        rC[k*(nn*na) + j] = out_r[k];
        iC[k*(nn*na) + j] = out_i[k];
      }

    }

    write_out_data_for_component(fout_id, rC, iC,
         nlev, nn, na, gfnames[i], time);
  }

  /* Cleanup */
  HDF5_ERROR(H5Fclose(file_id));
  HDF5_ERROR(H5Fclose(fout_id));

  free(in_r);
  free(in_i);
  free(out_r);
  free(out_i);
  free(rC);
  free(iC);
  free(rCt);
  free(iCt);
  free(time);

  return EXIT_SUCCESS;
}

static void filter(const int n,
                   const double c, /* c = dt/Tscale */
                   const double *in,
                   double *out)
{
  const double k1 = 2.0 * c / (2.0 + c);
  const double k2 = (2.0 - c) / (2.0 + c);
  int i;


  out[0] = in[0];

  for (i=1; i < n; i++)
  {
    const double in_a = 0.5 * (in[i] + in[i-1]);
    out[i] = k1 * in_a + k2 * out[i-1];
  }
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

  group_id = HDF5_ERROR(H5Gcreate(fileid, name, 0));

  dataspace_id = HDF5_ERROR(H5Screate_simple(1, &oD2, NULL));
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

  dataspace_id =  HDF5_ERROR(H5Screate_simple(1, &oD1, NULL));
  attribute_id = H5Acreate(group_id, "Rout", H5T_NATIVE_DOUBLE,
                  dataspace_id, HDF5_ERROR(H5P_DEFAULT));
  HDF5_ERROR(H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &Rout));
  HDF5_ERROR(H5Aclose(attribute_id));
  HDF5_ERROR(H5Sclose(dataspace_id));

  HDF5_ERROR(H5Gclose(group_id));
}


static void write_out_data_for_component(
                   hid_t fileid,
                   double *rC,
                   double *iC,
                   int nlev,
                   int nn,
                   int na,
                   const char *name,
                   double *time)
{
  static int first = 1;
#define BUFFSIZE 1024
  char buff[BUFFSIZE];
  int i;
  hsize_t dims[2];
  herr_t  status;
  dims[0] = nn;
  dims[1] = na;

  hid_t group_id, dataset_id, attribute_id, dataspace_id;

  fprintf(stderr, "Writing %s\n", name);

  if (first)
  {
    first=0;

    for (i=0; i < nlev; i++)
    {
      hsize_t oneD = 1;
      snprintf(buff, BUFFSIZE-1, "/%d", i);
      group_id = HDF5_ERROR(H5Gcreate(fileid, buff, 0));
      dataspace_id =  HDF5_ERROR(H5Screate_simple(1, &oneD, NULL));
      attribute_id = H5Acreate(group_id, "Time", H5T_NATIVE_DOUBLE,
        dataspace_id, HDF5_ERROR(H5P_DEFAULT));
      HDF5_ERROR(H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, time+i));
      HDF5_ERROR(H5Aclose(attribute_id));
      HDF5_ERROR(H5Sclose(dataspace_id));
      HDF5_ERROR(H5Gclose(group_id));
    }
  }

  for (i=0; i < nlev; i++)
  {
    snprintf(buff, BUFFSIZE-1, "/%d/%s", i, name);
    group_id = HDF5_ERROR(H5Gcreate(fileid, buff, 0));
    HDF5_ERROR(H5Gclose(group_id));

    snprintf(buff, BUFFSIZE-1, "/%d/%s/re", i, name);
    dataspace_id =  HDF5_ERROR(H5Screate_simple(2, dims, NULL));
    dataset_id =  HDF5_ERROR(H5Dcreate(fileid, buff, H5T_NATIVE_DOUBLE,
                                       dataspace_id, HDF5_ERROR(H5P_DEFAULT)));
    HDF5_ERROR(H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                                 H5S_ALL, H5P_DEFAULT, rC + i*(na*nn)));
    HDF5_ERROR(H5Dclose(dataset_id));
    HDF5_ERROR(H5Sclose(dataspace_id));


    snprintf(buff, BUFFSIZE-1, "/%d/%s/im", i, name);
    dataspace_id =  HDF5_ERROR(H5Screate_simple(2, dims, NULL));
    dataset_id =  HDF5_ERROR(H5Dcreate(fileid, buff, H5T_NATIVE_DOUBLE,
                                       dataspace_id, H5P_DEFAULT));
    HDF5_ERROR(H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                                 H5S_ALL, H5P_DEFAULT, iC + i*(na*nn)));
    HDF5_ERROR(H5Dclose(dataset_id));
    HDF5_ERROR(H5Sclose(dataspace_id));

  }
}

static void read_meta(hid_t file_id,
                   int *nn,
                   int *na,
                   double *Rin,
                   double *Rout,
                   int *s)
{

  herr_t  status;
  int dim[2];

  HDF5_ERROR(H5LTget_attribute_int(file_id,"/metadata/","spin", s));
  HDF5_ERROR(H5LTget_attribute_int(file_id,"/metadata/","dim", dim));
  HDF5_ERROR(H5LTget_attribute_double(file_id,"/metadata/","Rin", Rin));
  HDF5_ERROR(H5LTget_attribute_double(file_id,"/metadata/","Rout", Rout));

  *nn = dim[0];
  *na = dim[1];
}
