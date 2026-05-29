#include <iostream>
#include <complex>
#include <math.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "Legendre.hh"
#include "sYlm.hh"
#define BUFFSIZE 1024

using namespace std;

struct decomp_info
{
  int s;
  int nl;
  int nn;
  double Rin;
  double Rout;
  const double *re;
  const double *im;
};

struct output_info
{
  double xminus;
  double yminus;
  double xplus;
  double yplus;
  double z;
  int nx;
  int ny;
  int nmax;
  int lmax;
  int mmax;
};

static void dump_2d_data (const struct output_info *oinfo,
                          const struct decomp_info *dinfo);

int main(int argc, char **argv)
{
  if (argc < 13)
  {
    cerr << "Usage: " << argv[0] << " hdf5file iteration xm ym xp yp z nx ny nmax lmax mmax" << endl;
    return -1;
  }

  const char *filename = argv[1];
  FILE *file = fopen(filename, "rb");
  if (!file)
  {
    cerr << "File: " << filename << " does not exist !" << endl;
    return -2;
  }
  else
  {
    fclose(file);
    file = NULL;
  }

  const int iteration = atol(argv[2]);
  const double xm = atof(argv[3]);
  const double ym = atof(argv[4]);
  const double xp = atof(argv[5]);
  const double yp = atof(argv[6]);
  const double z = atof(argv[7]);
  const int nx = atol(argv[8]);
  const int ny= atol(argv[9]);
  const int nmax = atol(argv[10]);
  const int lmax = atol(argv[11]);
  const int mmax = atol(argv[12]);

  { /* check parameters */
  }

  
  
  char buff[BUFFSIZE];
  double Rin, Rout;
  int spin;
  int dim[2];

  hid_t       file_id, dataset_id;
  herr_t      status;

  file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  status = H5LTget_attribute_int(file_id,"/","spin",&spin);
  status = H5LTget_attribute_int(file_id,"/","dim",dim);
  status = H5LTget_attribute_double(file_id,"/","Rin",&Rin);
  status = H5LTget_attribute_double(file_id,"/","Rout",&Rout);

  /* the order of the dimensions has changed to C
     so dim[0] = nn
        dim[1] = na
  */
  const int nn = dim[0];
  const int na = dim[1];
  const int nl = static_cast<int>(round(sqrt(na+spin*spin)-abs(spin)));

  cout << "# Spin = " << spin << endl;
  cout << "# Rin = " << Rin << endl;
  cout << "# Rout = " << Rout << endl;
  cout << "# NA = " << na << endl;
  cout << "# NL = " << nl << endl;
  cout << "# NN = " << nn << endl;

  double time;
  snprintf(buff, BUFFSIZE-1, "/%d", iteration);
  status = H5LTget_attribute_double(file_id,buff, "Time", &time);
  cout << "# Time = " << time << endl;

  double *re = NULL;
  double *im = NULL;

  re = new  double[na*nn];
  im = new  double[na*nn];

  snprintf(buff, BUFFSIZE-1, "/%d/re", iteration);
  H5LTread_dataset_double(file_id, buff, re);

  snprintf(buff, BUFFSIZE-1, "/%d/im", iteration);
  H5LTread_dataset_double(file_id, buff, im);

  status = H5Fclose(file_id);

  struct decomp_info dinfo = {spin, nl, nn, Rin, Rout, re, im};
  struct output_info oinfo = {xm, ym, xp, yp, z, nx, ny, nmax, lmax, mmax};
  dump_2d_data(&oinfo, &dinfo);

  delete [] re;
  delete [] im;
  return 0;


}

static void dump_2d_data (const struct output_info *oinfo,
                          const struct decomp_info *dinfo)
{
  using namespace decomp_Legendre;
  using namespace decomp_sYlm;
  const double xm = oinfo->xminus;
  const double xp = oinfo->xplus;
  const double ym = oinfo->yminus;
  const double yp = oinfo->yplus;
  const double z = oinfo->z;
  const int nx = oinfo->nx;
  const int ny = oinfo->ny;
  const int nmax = oinfo->nmax;
  const int lmax = oinfo->lmax;
  const int mmax = oinfo->mmax;
  const double Rout = dinfo->Rout;
  const double Rin = dinfo->Rin;
  const int nl = dinfo->nl;
  const int nn = dinfo->nn;
  const int s = dinfo->s;
  const int na = nl*(nl + 2*abs(s));
  const double *re = dinfo->re;
  const double *im = dinfo->im;

  const double dx = (xp - xm) / (static_cast<double>(nx-1));
  const double dy = (yp - ym) / (static_cast<double>(ny-1));
  const double dXdr = 2.0 / (Rout - Rin);

  for (int j=0; j < ny; j++)
  {
    const double y = ym + j*dy;
    for (int i=0; i< nx; i++)
    {
      const double x = xm + i*dx;
      const double r = sqrt(x*x + y*y + z*z) ;
      const double mu = z/(r+1.0e-100);
      const double phi = atan2(y,x);
      const double Xl = -1.0 +(r-Rin)*dXdr;

      complex<double> val = 0;

      for (int n =0; n < nn; n++)
      {
        int sfun = 0;
        const double Pn = fabs(Xl)<1.0 ? LegendreP(n, Xl) : 0;
        if (n>nmax)
        {
          break;
        }

        for (int ll =0; ll < nl; ll++)
        {
          const int l = abs(s) + ll;
          for (int m=-l; m <=l; m++, sfun++)
          {
            const complex<double> coef (re[n*na + sfun], im[n*na + sfun]);
            const complex<double> sYlm = sYlm_mu(s,l,m,mu,phi);
            if (l > lmax || abs(m) > mmax)
            {
              continue;
            }
            val += coef*sYlm*Pn;
          }
        }
      }
      cout << x << "\t" << y << "\t"
           << r*val.real() << "\t"
           << r*val.imag() << "\n";
    }
    cout << "\n";
  }
}
