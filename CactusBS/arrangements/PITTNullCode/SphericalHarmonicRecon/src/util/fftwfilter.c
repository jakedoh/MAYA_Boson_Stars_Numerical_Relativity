#define H5_USE_16_API 1
#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>
#include "util.h"

#define pow2(x_) ((x_)*(x_))
#define pow3(x_) ((x_)*(x_)*(x_))

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
    exit (1);
  }
  return error_code;
}

static void print_usage(const char *name)
{
  fprintf(stderr, "\n\nusage: %s nlast omegamax kdamp hdf5file\n"
      "where nlast is the last iteration,\n"
      "omegamax is the largest frequency (in physical units),\n"
      "kdamp is the damping length in fft k index,\n"
      "hdf5file is the filename of the hdf5 containing the data.\n"
      "  The transform is truncated by multipling g(k)\n"
      "   by 1/2[1-erf(-(k-kmax)/kdamp)], where \n"
      "   kmax is obtained from omegamax.\n"
      "This program performs an fft filtering on the data in the\n"
      "input hdf5file 'FILE.h5' and generates a new file\n"
      "'FILE_ft.hf' with the filtered data and its time derivative\n\n",
       name);
}

static void fft_filter(const int nmax,
                   const int kmax,
                   const double ikd,
                   fftw_complex out[]);

static void fft_time_deriv(const int nmax,
                           const double dt,
                           fftw_complex *out);
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
                   double *rCt,
                   double *iCt,
                   int nlev,
                   int nn,
                   int na,
                   const char *name,
                   double *time);

int main(int argc, char **argv)
{
  if (argc != 5)
  {
    print_usage(argv[0]);
    return EXIT_FAILURE;
  }

#define SMALL 1.0e-90
  const int last = atoi(argv[1]);
  const double omegamax = atof(argv[2]);
  const double kdamp = atof(argv[3]);
  const char *infile = argv[4];
  const int skip = 1;
  const int nlev = last/skip + 1;
  const double ikd = 1.0 / (kdamp + SMALL);
  hid_t       file_id;
  char *outfile = NULL;

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


  fftw_complex *in, *out;
  fftw_plan pf;
  fftw_plan pr;

  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nlev);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nlev);

  pf = fftw_plan_dft_1d(nlev, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  pr = fftw_plan_dft_1d(nlev, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);


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

  int kmax = (int)(omegamax * nlev *
                 ((time[nlev-1] - time[0])/
                 (2.0*M_PI*(double)(nlev - 1))));  

  assert(kmax>0);

  if (kmax > nlev/2)
  {
    kmax = nlev;
    fprintf(stderr, "WARNING, CANT FILTER, FREQ TO HIGH\n");
  }


  for (i=0; i < 10; i++)
  {
    fprintf(stderr, "Reading in %s\n", gfnames[i]);
    read_in_data_for_component(file_id, rC, iC,
               nlev, skip, nn*na, gfnames[i]);
    fprintf(stderr, "done Reading in %s\n", gfnames[i]);

    for (j=0; j < nn*na; j++)
    {
      fprintf(stderr, "transforming %d of %d \n", j, nn*na);

/* this bunch of nastiness adds a quadratic to f(t) such that 
   for g(t) = f(t) + a t + b t^2
   g(0) = g(T)
   g'(0) = g'(T)
  we then perform the FFT on g(t), reducing the Gibbs phenomenon
  f' is calculated using 1 sided derivative
*/


/* new version: Here we require that
        (g[n-1] - g[n-2])/dt * dt + g[n-1] = g[0]  i.e. if we extrapolate
                                                   to get g[n], then
                                                   g[n] = g[0] by
                                                   periodicity.
      and
        g[0] - g[n-1] = g[1] - g[0]      i.e. The two one sided derivatives
                                          at i=0 agree
*/
#define NEW_FILTER
#if defined (NEW_FILTER)
      const double idt = 1.0 / (time[1] - time[0]);
      const double idt2 = idt * idt;
      const int n = nlev;
      const double oon = 1.0 / (double) n;
      const double oonm2 = 1.0 / (double) (n-2);
#define f(i_) rC[(i_)*(nn*na) + j]
      const double ra = -0.5 * idt * oonm2 * oon * 
             ((-6 + n*(2 + n))*f(0) - (-2 + n*n)*f(1) - 
             (2 + (-2 + n)*n)*f(-2 + n) + (6 + (-4 + n)*n)*f(-1 + n));
      const double rb = -0.5 * idt2 * oonm2 *
             (-f(0) + f(1) + f(-2 + n) - f(-1 + n));
#undef f
#define f(i_) iC[(i_)*(nn*na) + j]
      const double ia = -0.5 * idt * oonm2 * oon * 
             ((-6 + n*(2 + n))*f(0) - (-2 + n*n)*f(1) - 
             (2 + (-2 + n)*n)*f(-2 + n) + (6 + (-4 + n)*n)*f(-1 + n));
      const double ib = -0.5 * idt2 * oonm2 *
             (-f(0) + f(1) + f(-2 + n) - f(-1 + n));

      const double rc = 0;
      const double ic = 0;

#undef f
#elif defined (OLD_FILTER)
  /* this one has a bug. It assumes that
        g[n-1] should be equal to g[0]. This is wrong.
     we actually should have g[n] = g[0] (where g[n]
     is the hypothetical data point beyond the last element
     of g[] */

      const double idt = 1.0 / (time[1] - time[0]);
      const double T = time[nlev-1] - time[0];
      const double r_l =rC[(0)*(nn*na) + j];
      const double r_u =rC[(nlev-1)*(nn*na) + j];
      const double i_l =iC[(0)*(nn*na) + j];
      const double i_u =iC[(nlev-1)*(nn*na) + j];

      const double rp_l =  idt*(rC[(1)*(nn*na) + j] - rC[(0)*(nn*na) + j]);
      const double rp_u =  idt*(rC[(nlev-1)*(nn*na) + j] - rC[(nlev-2)*(nn*na) + j]);
      const double ip_l =  idt*(iC[(1)*(nn*na) + j] - iC[(0)*(nn*na) + j]);
      const double ip_u =  idt*(iC[(nlev-1)*(nn*na) + j] - iC[(nlev-2)*(nn*na) + j]);

      const double rb = 0.5 * (rp_u - rp_l) / T;
      const double ib = 0.5 * (ip_u - ip_l) / T;
      const double ra = (r_u - r_l)/T - rb*T;
      const double ia = (i_u - i_l)/T - ib*T;
      const double rc = 0;
      const double ic = 0;

#else
      const double rb = 0;
      const double ib = 0;
      const double ra = 0;
      const double ia = 0;
      const double rc = 0;
      const double ic = 0;

#endif

      for (k=0; k < nlev; k++)
      {
        const double tt = time[k] - time[0];
        const double rpoly = ra * tt + rb * tt * tt + rc * tt * tt * tt;
        const double ipoly = ia * tt + ib * tt * tt + ic * tt * tt * tt;

        in[k][0]= rC[(k)*(nn*na) + j] - rpoly;
        in[k][1]= iC[(k)*(nn*na) + j] - ipoly;
      }

      fftw_execute(pf);

      fft_filter(nlev, kmax, ikd, out);

      fftw_execute(pr);

      /* fftw transform is not normalized,
         so we need to renormalize here */
      for (k=0; k < nlev; k++)
      {
        const double tt = time[k] - time[0];
        const double rpoly = ra * tt + rb * tt * tt + rc * tt * tt * tt;
        const double ipoly = ia * tt + ib * tt * tt + ic * tt * tt * tt;

        rC[(k)*(nn*na) + j] = in[k][0] / ((double)(nlev)) + rpoly;
        iC[(k)*(nn*na) + j] = in[k][1] / ((double)(nlev)) + ipoly;
      }

      fft_time_deriv(nlev, dt, out);

      fftw_execute(pr);

      /* fftw transform is not normalized,
         so we need to renormalize here */
      for (k=0; k < nlev; k++)
      {
        const double tt = time[k] - time[0];
        const double drpoly = ra + 2.0 * rb * tt + 3.0 * rc * tt * tt;
        const double dipoly = ia + 2.0 * ib * tt + 3.0 * ic * tt * tt;

        rCt[(k)*(nn*na) + j] = in[k][0] / ((double)(nlev)) + drpoly;
        iCt[(k)*(nn*na) + j] = in[k][1] / ((double)(nlev)) + dipoly;
      }


    }

    write_out_data_for_component(fout_id, rC, iC, rCt, iCt,
         nlev, nn, na, gfnames[i], time);

#undef DEBUG_TEST
#ifdef DEBUG_TEST
    if (i==0)
    {
      for (k=0; k < nlev; k++)
      {
#define L 2
#define M 2
#define MODE (L*L+L+M)
         printf("%20.16e %20.16e %20.16e %20.16e %20.16e\n",
          time[k], rC[k*nn*na+MODE], iC[k*nn*na+MODE],  rCt[k*nn*na+MODE], iCt[k*nn*na+MODE]);
      }
      fflush(stdout);
      exit(0);
    }
#endif
  }

  /* Cleanup */
  HDF5_ERROR(H5Fclose(file_id));
  HDF5_ERROR(H5Fclose(fout_id));

  fftw_destroy_plan(pf);
  fftw_destroy_plan(pr);
  fftw_free(in);
  fftw_free(out);
  free(rC);
  free(iC);
  free(rCt);
  free(iCt);
  free(time);

  return EXIT_SUCCESS;
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

static void fft_filter(const int nmax,
           const int kmax,
           const double ikd,
           fftw_complex out[])
{
  int i;

  for (i=0; i <= nmax/2; i++)
  {
    const double f = 0.5 * (1 + erf(-ikd * (double)(i-kmax)));
    out[i][0] *= f;
    out[i][1] *= f;
    out[nmax-i][0] *= f;
    out[nmax-i][1] *= f;

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
  group_id= HDF5_ERROR(H5Gcreate(fileid, name, 0));
  dataspace_id =  HDF5_ERROR(H5Screate_simple(1, &oD2, NULL));
  attribute_id = HDF5_ERROR(H5Acreate(group_id, "dim", H5T_NATIVE_INT,
                       dataspace_id, H5P_DEFAULT));
  HDF5_ERROR(H5Awrite(attribute_id, H5T_NATIVE_INT, ds));
  HDF5_ERROR(H5Aclose(attribute_id));
  HDF5_ERROR(H5Sclose(dataspace_id));

  dataspace_id =  HDF5_ERROR(H5Screate_simple(1, &oD1, NULL));
  attribute_id = HDF5_ERROR(H5Acreate(group_id, "spin", H5T_NATIVE_INT,
                        dataspace_id, H5P_DEFAULT));
  HDF5_ERROR(H5Awrite(attribute_id, H5T_NATIVE_INT, &s));
  HDF5_ERROR(H5Aclose(attribute_id));
  HDF5_ERROR(H5Sclose(dataspace_id));

  dataspace_id =  HDF5_ERROR(H5Screate_simple(1, &oD1, NULL));
  attribute_id = HDF5_ERROR(H5Acreate(group_id, "Rin", H5T_NATIVE_DOUBLE,
                  dataspace_id, H5P_DEFAULT));
  HDF5_ERROR(H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &Rin));
  HDF5_ERROR(H5Aclose(attribute_id));
  HDF5_ERROR(H5Sclose(dataspace_id));

  dataspace_id =  HDF5_ERROR(H5Screate_simple(1, &oD1, NULL));
  attribute_id = HDF5_ERROR(H5Acreate(group_id, "Rout", H5T_NATIVE_DOUBLE,
                  dataspace_id, H5P_DEFAULT));
  HDF5_ERROR(H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &Rout));
  HDF5_ERROR(H5Aclose(attribute_id));
  HDF5_ERROR(H5Sclose(dataspace_id));

  HDF5_ERROR(H5Gclose(group_id));
}


static void write_out_data_for_component(
                   hid_t fileid,
                   double *rC,
                   double *iC,
                   double *rCt,
                   double *iCt,
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
      attribute_id = HDF5_ERROR(H5Acreate(group_id, "Time", H5T_NATIVE_DOUBLE,
        dataspace_id, H5P_DEFAULT));
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
           dataspace_id, H5P_DEFAULT));
    HDF5_ERROR(H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, rC + (i)*(na*nn)));
    HDF5_ERROR(H5Dclose(dataset_id));
    HDF5_ERROR(H5Sclose(dataspace_id));

    snprintf(buff, BUFFSIZE-1, "/%d/%s/dt_re", i, name);
    dataspace_id =  HDF5_ERROR(H5Screate_simple(2, dims, NULL));
    dataset_id =  HDF5_ERROR(H5Dcreate(fileid, buff, H5T_NATIVE_DOUBLE,
           dataspace_id, H5P_DEFAULT));
    HDF5_ERROR(H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, rCt + (i)*(na*nn)));
    HDF5_ERROR(H5Dclose(dataset_id));
    HDF5_ERROR(H5Sclose(dataspace_id));



    snprintf(buff, BUFFSIZE-1, "/%d/%s/im", i, name);
    dataspace_id =  HDF5_ERROR(H5Screate_simple(2, dims, NULL));
    dataset_id =  HDF5_ERROR(H5Dcreate(fileid, buff, H5T_NATIVE_DOUBLE,
           dataspace_id, H5P_DEFAULT));
    HDF5_ERROR(H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, iC + (i)*(na*nn)));
    HDF5_ERROR(H5Dclose(dataset_id));
    HDF5_ERROR(H5Sclose(dataspace_id));

    snprintf(buff, BUFFSIZE-1, "/%d/%s/dt_im", i, name);
    dataspace_id =  HDF5_ERROR(H5Screate_simple(2, dims, NULL));
    dataset_id =  HDF5_ERROR(H5Dcreate(fileid, buff, H5T_NATIVE_DOUBLE,
           dataspace_id, H5P_DEFAULT));
    HDF5_ERROR(H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, iCt + (i)*(na*nn)));
    HDF5_ERROR(H5Dclose(dataset_id));
    HDF5_ERROR(H5Sclose(dataspace_id));
  }
}

static void fft_time_deriv(const int nmax,
                           const double dt,
                           fftw_complex *out)
{
  int i;

  const double wm = M_PI*2.0 / (nmax *dt);
  out[0][0] = 0;
  out[0][1] = 0;

  /* negative frequencies are stored in elements (n-i) */
  for (i=1; i<=nmax/2; i++)
  {
    const double omega = (double) i * wm;
    const double re = out[i][0];
    const double im = out[i][1];
    const double re2 = out[nmax-i][0];
    const double im2 = out[nmax-i][1];

    out[i][0] = -im * omega;
    out[i][1] = +re * omega;
    out[nmax-i][0] = +im2 * omega;
    out[nmax-i][1] = -re2 * omega;
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
