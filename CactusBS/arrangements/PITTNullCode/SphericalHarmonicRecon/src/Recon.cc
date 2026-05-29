#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "util_Table.h"
#include "util_ErrorCodes.h"
#include "Chebyshev.hh"
#include <stdio.h>
#include <assert.h>

#include <iostream>
#include <math.h>
#include <complex>
#include <mpi.h>

#define H5_USE_16_API 1
#include <hdf5.h>
#include <hdf5_hl.h>

static int hdf5_warn_level = CCTK_WARN_ABORT;

#define HDF5_ERROR(fn_call)                                                   \
        do {                                                                  \
          hid_t _error_code = fn_call;                                        \
                                                                              \
                                                                              \
          if (_error_code < 0)                                                \
          {                                                                   \
            CCTK_VWarn (hdf5_warn_level, __LINE__, __FILE__, CCTK_THORNSTRING,\
                        "HDF5 call '%s' returned error code %d",              \
                        #fn_call, (int)_error_code);                          \
          }                                                                   \
        } while (0)


#undef SPHERICALHARMONICRECON_DEBUG

#define BUFFSIZE 1024
#define NUM_TIMELEVELS 5

#define NUM_METRIC_COMPONENTS 10
#define GXX 0
#define GXY 1
#define GXZ 2
#define GYY 3
#define GYZ 4
#define GZZ 5
#define BETAX 6
#define BETAY 7
#define BETAZ 8
#define ALP 9

static const char *gfname[] = 
    { 
      "gxx",   /* GXX = 0 */
      "gxy",   /* GXY = 1 */
      "gxz",   /* GXZ = 2 */
      "gyy",   /* GYY = 3 */
      "gyz",   /* GYZ = 4 */
      "gzz",   /* GZZ = 5 */
      "betax", /* BETAX = 6 */
      "betay", /* BETAY = 7 */
      "betaz", /* BETAZ = 8 */
      "alp"    /* ALP = 9 */
    };

/* state info vars */
static int initialized = 0;
static int read_data = 0;

/*
   spin: the spinweight of the spherical harmonics (0) 
   Rin: The inner radius of the extraction world tube
   Rout: The outer radius of the extraction world  tube
   nn: Number of Chebyshev polynomials (nn = nmax - 1)
   na: number of angular coefficients (Sum (2l+1) )
   nl: numbero f l modes (nl = lmax + 1)
   time0: first dumptime in data file (usually t=0)
   dtime: time difference between dumps in file
*/

static int spin;
static CCTK_REAL Rin;
static CCTK_REAL Rout;

static int nn;
static int na;
static int nl;

static CCTK_REAL time0;
static CCTK_REAL dtime;

static CCTK_REAL U[CHEBYSHEV_NMAX+1];
static CCTK_REAL Ux[CHEBYSHEV_NMAX+1];
static CCTK_REAL Ur[CHEBYSHEV_NMAX+1];

static CCTK_REAL *retmp[NUM_TIMELEVELS];
static CCTK_REAL *imtmp[NUM_TIMELEVELS];
static CCTK_REAL *rertmp[NUM_TIMELEVELS];
static CCTK_REAL *imrtmp[NUM_TIMELEVELS];
static CCTK_REAL *rettmp[NUM_TIMELEVELS];
static CCTK_REAL *imttmp[NUM_TIMELEVELS];

static CCTK_REAL *ReC[NUM_METRIC_COMPONENTS];
static CCTK_REAL *ImC[NUM_METRIC_COMPONENTS];
static CCTK_REAL *ReCr[NUM_METRIC_COMPONENTS];
static CCTK_REAL *ImCr[NUM_METRIC_COMPONENTS];
static CCTK_REAL *ReCt[NUM_METRIC_COMPONENTS];
static CCTK_REAL *ImCt[NUM_METRIC_COMPONENTS];

static CCTK_REAL *reall = NULL;
static CCTK_REAL *imall = NULL;
static CCTK_REAL *reallt = NULL;
static CCTK_REAL *imallt = NULL;

static int MyProc = -1;
static MPI_Comm comm_world = MPI_COMM_NULL;

static char buff[BUFFSIZE];


/* local function prototypes */

static void UxtoUr(void);

static int angcoef(const hid_t file_id, const int comp,
           const int it0, const CCTK_REAL r_extract,
           double **retmp, double **imtmp,
           double **rettmp, double **imttmp,
           double **rertmp, double **imrtmp,
           const int read_dt);

static int allocate_coef(void);

static int constcoef_4th(const hid_t file_id, const int comp,
    const CCTK_REAL cctk_time, const CCTK_REAL r_extract,
    const int read_dt);

static int constcoef_2nd(const hid_t file_id, const int comp,
    const CCTK_REAL cctk_time, const CCTK_REAL r_extract,
    const int read_dt);

static void print_coefs(const int comp);

static CCTK_REAL GetReCoef(const int comp,
         const int l, const int m);

static CCTK_REAL GetImCoef(const int comp,
         const int l, const int m);


static CCTK_REAL GetReDrCoef(const int comp,
         const int l, const int m);

static CCTK_REAL GetImDrCoef(const int comp,
         const int l, const int m);

static CCTK_REAL GetReDtCoef(const int comp,
         const int l, const int m);

static CCTK_REAL GetImDtCoef(const int comp,
         const int l, const int m);



using namespace std;
using namespace recon_Chebyshev;

/* externally visible function */
extern "C"
{
  /* Cactus Scheduled functions */

  void SphericalHarmonicRecon_Startup(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_PARAMETERS;

    hid_t       file_id;
    herr_t      status;

    MyProc = CCTK_MyProc(cctkGH);

    if (CCTK_IsFunctionAliased("GetMPICommWorld"))
    {
      comm_world = *(const MPI_Comm *)GetMPICommWorld(cctkGH);
    }
    else
    {
      comm_world = MPI_COMM_WORLD;
    }

    if (!MyProc)
    {
      int dim[2];

      file_id = H5Fopen(metric_data_filename, H5F_ACC_RDONLY, H5P_DEFAULT);    
      if(file_id<0) {
        string message = string("cannot open the file ") +metric_data_filename;
        CCTK_WARN(CCTK_WARN_ABORT, message.c_str());
      }

      if (override_extraction_parameters)
      {
        CCTK_WARN(CCTK_WARN_ALERT,
             "WARNING: !!!!!!!!!!!!!!!!! \n"
             "overrriding extraction parameter in metric_decomp\n"
             "This is a really stupid thing to do... unless you have to.");

        spin = override_spin;
        dim[0] = override_nn;
        dim[1] = override_na;
        Rout = override_Rout;
        Rin = override_Rin;
      }
      else
      {
        HDF5_ERROR(status = H5LTget_attribute_int(file_id,"/metadata/","spin",&spin));
        HDF5_ERROR(status = H5LTget_attribute_int(file_id,"/metadata/","dim",dim));
        HDF5_ERROR(status = H5LTget_attribute_double(file_id,"/metadata/","Rin",&Rin));
        HDF5_ERROR(status = H5LTget_attribute_double(file_id,"/metadata/","Rout",&Rout));
      }

      nn = dim[0];
      na = dim[1];
      nl = static_cast<int>(round(sqrt(double(na+spin*spin))-abs(double(spin))));

      HDF5_ERROR(status = H5LTget_attribute_double(file_id,"/0/", "Time", &time0));
      HDF5_ERROR(status = H5LTget_attribute_double(file_id,"/1/", "Time", &dtime));
      dtime -= time0;


      HDF5_ERROR(status = H5Fclose(file_id));
    }

    const int root = 0; /* proc 0 sends the data to everyone else */

    int pack_int[] = {spin, nn, na, nl};
    double pack_double[] = {Rin, Rout, time0, dtime};

    MPI_Bcast (pack_int, sizeof pack_int / sizeof pack_int[0],
          MPI_INT, root, comm_world);
    MPI_Bcast (pack_double, sizeof pack_double / sizeof pack_double[0],
          MPI_DOUBLE, root, comm_world);

    spin = pack_int[0];
    nn = pack_int[1];
    na = pack_int[2];
    nl = pack_int[3];

    Rin = pack_double[0];
    Rout = pack_double[1];
    time0 = pack_double[2];
    dtime = pack_double[3];

    cout << "SphericalHarmonicRecon Parameters" << endl;
    cout << "# Spin = " << spin << endl;
    cout << "# Rin = " << Rin << endl;
    cout << "# Rout = " << Rout << endl;
    cout << "# NA = " << na << endl;
    cout << "# NL = " << nl << endl;
    cout << "# NN = " << nn << endl;
    cout << "# dt = " << dtime << endl;


    assert (r_extract > Rin && r_extract < Rout);
    assert(nn < CHEBYSHEV_NMAX);

    allocate_coef();


    const CCTK_REAL x = (2.0 * r_extract - (Rout + Rin))/(Rout - Rin);
    ChebyshevU(nn-1, x,  U, Ux);
    UxtoUr();

    initialized = 1;
  }

  void SphericalHarmonicRecon_ReadData(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    hid_t       file_id;
    herr_t      status;

    if (! initialized)
    {
      CCTK_WARN(CCTK_WARN_ABORT, "Schedule mismatch. \n"
         "      SphericalHarmonicRecon_ReadData must be called after\n"
         "      SphericalHarmonicRecon_Startup");
    }

    if (!MyProc)
    {
      file_id = H5Fopen(metric_data_filename, H5F_ACC_RDONLY, H5P_DEFAULT);    
      if(file_id<0) {
        string message = string("cannot open the file ") +metric_data_filename;
        CCTK_WARN(CCTK_WARN_ABORT, message.c_str());
      }

      if (order==2)
      {

        constcoef_2nd(file_id, GXX, cctk_time, r_extract, time_derivative_in_file);
        constcoef_2nd(file_id, GXY, cctk_time, r_extract, time_derivative_in_file);
        constcoef_2nd(file_id, GXZ, cctk_time, r_extract, time_derivative_in_file);
        constcoef_2nd(file_id, GYY, cctk_time, r_extract, time_derivative_in_file);
        constcoef_2nd(file_id, GYZ, cctk_time, r_extract, time_derivative_in_file);
        constcoef_2nd(file_id, GZZ, cctk_time, r_extract, time_derivative_in_file);
        constcoef_2nd(file_id, BETAX, cctk_time, r_extract, time_derivative_in_file);
        constcoef_2nd(file_id, BETAY, cctk_time, r_extract, time_derivative_in_file);
        constcoef_2nd(file_id, BETAZ, cctk_time, r_extract, time_derivative_in_file);
        constcoef_2nd(file_id, ALP, cctk_time, r_extract, time_derivative_in_file);
      }
      else if (order==4)
      {
        constcoef_4th(file_id, GXX, cctk_time, r_extract, time_derivative_in_file);
        constcoef_4th(file_id, GXY, cctk_time, r_extract, time_derivative_in_file);
        constcoef_4th(file_id, GXZ, cctk_time, r_extract, time_derivative_in_file);
        constcoef_4th(file_id, GYY, cctk_time, r_extract, time_derivative_in_file);
        constcoef_4th(file_id, GYZ, cctk_time, r_extract, time_derivative_in_file);
        constcoef_4th(file_id, GZZ, cctk_time, r_extract, time_derivative_in_file);
        constcoef_4th(file_id, BETAX, cctk_time, r_extract, time_derivative_in_file);
        constcoef_4th(file_id, BETAY, cctk_time, r_extract, time_derivative_in_file);
        constcoef_4th(file_id, BETAZ, cctk_time, r_extract, time_derivative_in_file);
        constcoef_4th(file_id, ALP, cctk_time, r_extract, time_derivative_in_file);
      }
      else
      {
        CCTK_WARN(CCTK_WARN_ABORT, "order not implemented");
      }

      HDF5_ERROR(status = H5Fclose(file_id));
    }


    for (int c=0; c < NUM_METRIC_COMPONENTS; c++)
    {
      const int root = 0; /* proc 0 sends the data */
      MPI_Bcast(ReC[c], na, MPI_DOUBLE, root, comm_world);
      MPI_Bcast(ImC[c], na, MPI_DOUBLE, root, comm_world);
      MPI_Bcast(ReCr[c], na, MPI_DOUBLE, root, comm_world);
      MPI_Bcast(ImCr[c], na, MPI_DOUBLE, root, comm_world);
      MPI_Bcast(ReCt[c], na, MPI_DOUBLE, root, comm_world);
      MPI_Bcast(ImCt[c], na, MPI_DOUBLE, 0, comm_world);
    }

#   ifdef SPHERICALHARMONICRECON_DEBUG
    cout << "TimeStep " << cctk_time << endl;
    print_coefs(GXX);
#   endif
    read_data = 1;
  }

  void SphericalHarmonicRecon_PostStep(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    read_data = 0;
  }

  /* Cactus Aliased Functions */

  CCTK_INT SphericalHarmonicRecon_GetParameters (CCTK_INT *lmax,
                                            CCTK_REAL *r_inner,
                                            CCTK_REAL *r_outer)
  {
    if (!initialized)
    {
      *lmax = -1;
      return -1;
    }
    *lmax = nl - 1; 
    *r_inner = Rin;
    *r_outer = Rout;

    return 0;
  }

  CCTK_INT SphericalHarmonicRecon_GetCurrentCoefs(
            CCTK_INT l, /* the ell mode */
            CCTK_INT m, /* the m mode */
            CCTK_REAL reC[], /* real part of (l,m) coefficient
                                of the 10 metric functions */
            CCTK_REAL imC[], /* ditto for the imaginary part */
            CCTK_REAL reCr[], /* real part of (l,m) coefficient
                                 of the r-derivative of the 10
                                 metric functions */
            CCTK_REAL imCr[], /* ditto for the imaginary part */
            CCTK_REAL reCt[], /* real part of (l,m) coefficient
                                 of the t-derivative of the 10
                                 metric functions */
            CCTK_REAL imCt[] /* ditto for the imaginary part */
   )
  {
    const int n = l*l + l + m;

    if (!read_data)
    {
      CCTK_WARN(CCTK_WARN_ABORT, "Error: data has not been read in");
      return -1;
    }

    for (int c = 0; c < NUM_METRIC_COMPONENTS; c++)
    {
      reC[c] = ReC[c][n]; 
      imC[c] = ImC[c][n]; 
      reCr[c] = ReCr[c][n]; 
      imCr[c] = ImCr[c][n]; 
      reCt[c] = ReCt[c][n]; 
      imCt[c] = ImCt[c][n]; 
    }

    return 0;
  }


}
 
static CCTK_REAL GetReCoef(const int comp, const int l, const int m)
{
  const int n = l*l + l + m;
  return ReC[comp][n];
}

static CCTK_REAL GetImCoef(const int comp, const int l, const int m)
{
  const int n = l*l + l + m;
  return ImC[comp][n];
}

static CCTK_REAL GetReDrCoef(const int comp, const int l, const int m)
{
  const int n = l*l + l + m;
  return ReCr[comp][n];
}

static CCTK_REAL GetImDrCoef(const int comp, const int l, const int m)
{
  const int n = l*l + l + m;
  return ImCr[comp][n];
}

static CCTK_REAL GetReDtCoef(const int comp, const int l, const int m)
{
  const int n = l*l + l + m;
  return ReCt[comp][n];
}

static CCTK_REAL GetImDtCoef(const int comp, const int l, const int m)
{
  const int n = l*l + l + m;
  return ImCt[comp][n];
}      

static void print_coefs(const int comp)
{
  for (int l=0; l < nl; l++)
  {
    for (int m=-l; m <= l; m++)
    {
      cout << "C(" << l << ", " << m << ") = ("
           << GetReDtCoef(comp, l, m) << ", "
           << GetImDtCoef(comp, l, m) << ")" << endl;
    }
  }
}

static void UxtoUr(void)
{
  const CCTK_REAL dxdr = 2.0 / (Rout - Rin);
  for (int n=0; n< nn; n++)
  {
    Ur[n] = Ux[n] * dxdr;
  }
}


static int constcoef_4th(const hid_t file_id, const int comp,
    const CCTK_REAL cctk_time, const CCTK_REAL r_extract,
    const int read_dt)
{
  const int it = static_cast<int>((cctk_time - time0)/dtime);
  const int it0 = it-2>0 ? it - 2: 0;

  const CCTK_REAL x0 = time0 + it0*dtime;
  const CCTK_REAL x = cctk_time;
  const CCTK_REAL dt = dtime;
  const CCTK_REAL dt2 = dt*dt;
  const CCTK_REAL dt3 = dt2*dt;
  const CCTK_REAL idt = 1.0 / dt;
  const CCTK_REAL idt4 = idt * idt * idt * idt;


#define INTERP_4th(f) \
  ((idt4*(-4*dt + x - x0)*(-3*dt + x - x0)*(-2*dt + x - x0)*\
      (-dt + x - x0)*f[0][j])/24. -\
   (idt4*(x - x0)*(-4*dt + x - x0)*(-3*dt + x - x0)*\
      (-2*dt + x - x0)*f[1][j])/6. +\
   (idt4*(x - x0)*(-4*dt + x - x0)*(-3*dt + x - x0)*\
      (-dt + x - x0)*f[2][j])/4. -\
   (idt4*(x - x0)*(-4*dt + x - x0)*(-2*dt + x - x0)*\
      (-dt + x - x0)*f[3][j])/6. +\
   (idt4*(x - x0)*(-3*dt + x - x0)*(-2*dt + x - x0)*\
      (-dt + x - x0)*f[4][j])/24.)

#define INTERPD_4th(f) \
    ((idt4*(dt3*(-25*f[0][j] + 48*f[1][j] -\
     36*f[2][j] + 16*f[3][j] - 3*f[4][j]) + \
     (x - x0)*(dt2*(35*f[0][j] - 104*f[1][j] +\
     114*f[2][j] - 56*f[3][j] + 11*f[4][j]) +\
     (x - x0)* (2*(x - x0)*(f[0][j] - 4*f[1][j] +\
     6*f[2][j] - 4*f[3][j] + f[4][j]) - 3*dt*\
     (5*f[0][j] - 18*f[1][j] + 24*f[2][j] - 14*f[3][j] +\
     3*f[4][j])))))/12. )

  angcoef(file_id, comp, it0, r_extract,
          retmp, imtmp,
          rettmp, imttmp,
          rertmp, imrtmp,
          read_dt);

  for (int j=0; j<na; j++)
  {
    ReC[comp][j] = INTERP_4th(retmp);
    ImC[comp][j] = INTERP_4th(imtmp);
    ReCr[comp][j] = INTERP_4th(rertmp);
    ImCr[comp][j] = INTERP_4th(imrtmp);
    if (read_dt)
    {
      ReCt[comp][j] = INTERP_4th(rettmp);
      ImCt[comp][j] = INTERP_4th(imttmp);
    }
    else
    {
      ReCt[comp][j] = INTERPD_4th(retmp);
      ImCt[comp][j] = INTERPD_4th(imtmp);
    }
  }

  return 0;
}

static int constcoef_2nd(const hid_t file_id, const int comp,
    const CCTK_REAL cctk_time, const CCTK_REAL r_extract,
    const int read_dt)
{
  const int it = static_cast<int>((cctk_time - time0)/dtime);
  const int it0 = it-2>0 ? it - 2: 0;

  const CCTK_REAL x0 = time0 + it0*dtime;
  const CCTK_REAL x = cctk_time;
  const CCTK_REAL dt = dtime;
  const CCTK_REAL dt2 = dt*dt;
  const CCTK_REAL idt = 1.0 / dt;
  const CCTK_REAL idt2 = idt * idt;
  const CCTK_REAL xm = x - x0;
  const CCTK_REAL xm2 = xm * xm;


#define INTERP_2nd(f) \
  ((idt2*(2*dt2*(3*f[1][j] - 3*f[2][j] + f[3][j]) + \
       xm2*(f[1][j] - 2*f[2][j] + f[3][j]) - \
       dt*xm*(5*f[1][j] - 8*f[2][j] + 3*f[3][j])))/2.)

#define INTERPD_2nd(f) \
  ((idt2*(dt*(-5*f[1][j] + 8*f[2][j] - 3*f[3][j]) +\
       2*xm*(f[1][j] - 2*f[2][j] + f[3][j])))/2.)

  angcoef(file_id, comp, it0, r_extract,
          retmp, imtmp,
          rettmp, imttmp,
          rertmp, imrtmp,
          read_dt);

  for (int j=0; j<na; j++)
  {
    ReC[comp][j] = INTERP_2nd(retmp);
    ImC[comp][j] = INTERP_2nd(imtmp);
    ReCr[comp][j] = INTERP_2nd(rertmp);
    ImCr[comp][j] = INTERP_2nd(imrtmp);
    if (read_dt)
    {
      ReCt[comp][j] = INTERP_2nd(rettmp);
      ImCt[comp][j] = INTERP_2nd(imttmp);
    }
    else
    {
      ReCt[comp][j] = INTERPD_2nd(retmp);
      ImCt[comp][j] = INTERPD_2nd(imtmp);
    }
  }

  return 0;
}

static int angcoef(const hid_t file_id, const int comp,
           const int it0, const CCTK_REAL r_extract,
           double **retmp, double **imtmp,
           double **rettmp, double **imttmp,
           double **rertmp, double **imrtmp,
           const int read_dt)
{
  herr_t      status;

  for (int t=0; t < NUM_TIMELEVELS; t++)
  {
    const int it = it0+t;

    snprintf(buff, BUFFSIZE-1, "/%d/%s/re", it, gfname[comp]);
    status = H5LTread_dataset_double(file_id, buff, reall);
    if (status)
    {
      cerr << "Could not read data from " << buff << endl;
      CCTK_WARN(CCTK_WARN_ABORT, "Dataset read failure");
    }

    snprintf(buff, BUFFSIZE-1, "/%d/%s/im", it, gfname[comp]);
    status = H5LTread_dataset_double(file_id, buff, imall);

    if (status)
    {
      cerr << "Could not read data from " << buff << endl;
      CCTK_WARN(CCTK_WARN_ABORT, "Dataset read failure");
    }

    if (read_dt)
    {
      snprintf(buff, BUFFSIZE-1, "/%d/%s/dt_re", it, gfname[comp]);
      status = H5LTread_dataset_double(file_id, buff, reallt);
      if (status)
      {
	cerr << "Could not read data from " << buff << endl;
        CCTK_WARN(CCTK_WARN_ABORT, "Dataset read failure");
      }

      snprintf(buff, BUFFSIZE-1, "/%d/%s/dt_im", it, gfname[comp]);
      status = H5LTread_dataset_double(file_id, buff, imallt);

      if (status)
      {
	cerr << "Could not read data from " << buff << endl;
        CCTK_WARN(CCTK_WARN_ABORT, "Dataset read failure");
      }

    }

    for (int j=0; j<na; j++)
    {
      retmp[t][j] = 0;
      imtmp[t][j] = 0;
      rertmp[t][j] = 0;
      imrtmp[t][j] = 0;
      rettmp[t][j] = 0;
      imttmp[t][j] = 0;
    }

    for (int n =0; n < nn; n++)
    {
      for (int j=0; j<na; j++)
      {
        retmp[t][j] += reall[n*na+j]*U[n];
        imtmp[t][j] += imall[n*na+j]*U[n];
        rertmp[t][j] += reall[n*na+j]*Ur[n];
        imrtmp[t][j] += imall[n*na+j]*Ur[n];
      }
    }

    if (read_dt)
    {
      for (int n =0; n < nn; n++)
      {
	for (int j=0; j<na; j++)
	{
	  rettmp[t][j] += reallt[n*na+j]*U[n];
	  imttmp[t][j] += imallt[n*na+j]*U[n];
	}
      }
    }
  }

  return 0;
}

static int allocate_coef(void)
{
  for (int t=0; t < NUM_TIMELEVELS; t++)
  {
    retmp[t] = new CCTK_REAL[na];
    imtmp[t] = new CCTK_REAL[na];
    rettmp[t] = new CCTK_REAL[na];
    imttmp[t] = new CCTK_REAL[na];
    rertmp[t] = new CCTK_REAL[na];
    imrtmp[t] = new CCTK_REAL[na];
  }

  for (int c=0; c < NUM_METRIC_COMPONENTS; c++)
  {
    ReC[c] = new CCTK_REAL[na];
    ImC[c] = new CCTK_REAL[na];
    ReCr[c] = new CCTK_REAL[na];
    ImCr[c] = new CCTK_REAL[na];
    ReCt[c] = new CCTK_REAL[na];
    ImCt[c] = new CCTK_REAL[na];
  }

  reall = new CCTK_REAL[nn*na];
  imall = new CCTK_REAL[nn*na];

  reallt = new CCTK_REAL[nn*na];
  imallt = new CCTK_REAL[nn*na];

  return 0;
}
