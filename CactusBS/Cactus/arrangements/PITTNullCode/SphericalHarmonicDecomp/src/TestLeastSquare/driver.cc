#include <iostream>
#include <complex>
#include <sys/time.h>
#include <cmath>
#include "sYlm.hh"
#include "decomp.hh"

#include "TestLeastSquare.h"

#ifdef USE_LEGENDRE
#  include "Legendre.hh"
#else
#  include "Chebyshev.hh"
#endif


#define Max(a_,b_) ((a_)>(b_)? (a_):(b_))
#define Min(a_,b_) ((a_)<(b_)? (a_):(b_))
#define ABS(x_) ((x_)>0 ? (x_) : (-(x_)))

using namespace decomp_sYlm;
#ifdef USE_LEGENDRE
using namespace decomp_Legendre;
#else
using namespace decomp_Chebyshev;
#endif
using namespace decomp_decompose;
using namespace std;

static void
    print_coeffs(int s, int nl, int nn, double *re, double *im);

static void print_coeffs2
    (int s, int nl, double *re, double *im);

static double get_walltime ()
  {
    // get the current time
    struct timeval tv;
    gettimeofday (& tv, 0);
    return tv.tv_sec + tv.tv_usec / double (1.0e6);
  }

int main(int argc, char ** argv)
{
  using namespace decomp_matrix_class;
  using namespace decomp_sYlm;
  using namespace decomp_decompose;
  using namespace std;

  const int s = -2;
  const int nl = 6;
  const int nn = 10;
  const int nmu = 21;
  const int nphi =43;
  const int nx = 20;
  double tbeg, tinvert, tpostdata, tpostdecomp;

  tbeg = get_walltime();
  const decomp_info *dinfo_p = initialize(s, nl, nn, nmu, nphi, nx);
  tinvert = get_walltime();

  const int nlmmodes = nl*(nl+2*ABS(s));

  double *re = new double[nmu*nphi*nx];
  double *im = new double[nmu*nphi*nx];

  double *ore = new double[nlmmodes*nx];
  double *oim = new double[nlmmodes*nx];

  double *ncolloc = new double[nx];
  double *mucolloc = new double[nmu*nphi];
  double *phicolloc = new double[nmu*nphi];
  double *lp = new double [nn*nx];
  dinfo_p->get_ncolloc(nx, ncolloc);
  dinfo_p->get_mucolloc(nmu*nphi, mucolloc);
  dinfo_p->get_phicolloc(nmu*nphi, phicolloc);

# ifdef USE_LEGENDRE
  cout << "Testing Legendre polynomial +sYlm decomp" << endl;
#else
  cout << "Testing ChebyshevU polynomial +sYlm decomp" << endl;
#endif


  for (int i=0; i < nn; i++)
  {
    for (int j=0; j < nx; j++)
    {
#ifdef USE_LEGENDRE
      lp[i*nx+j] = LegendreP(i,ncolloc[j]);
#else
      lp[i*nx+j] = ChebyshevU(i,ncolloc[j]);
#endif
    }
  }

  for (int r=0; r < nx; r++)
  {

    for (int a=0; a < nmu*nphi; a++)
    {
      const double mu = mucolloc[a];
      const double phi = phicolloc[a];

      complex<double> val = 0;
      for (int ll=0; ll < nl; ll++)
      {
        const int l = ABS(s) + ll;

	for (int m = -l; m < l+1; m++)
	{
          complex <double> ylm = sYlm_mu(s,l,m,mu,phi);
          for (int n=0; n < nn; n++)
          {
            complex<double> c(l,m);
            c *= static_cast <complex<double> > (lp[n*nx+r]*(n+1));
	    val += c*ylm;
	  }
	}
	re[r*nmu*nphi+a] = val.real();
	im[r*nmu*nphi+a] = val.imag();
      }
    }
  }
  tpostdata = get_walltime();

  decompose3D (dinfo_p, re, im, ore, oim);
  //decompose2D (dinfo_p, re, im);
  tpostdecomp = get_walltime();

  print_coeffs (s, nl, nn, ore, oim);
  //print_coeffs2 (s, nl, re, im);

  cout << "Timings:" << endl << "  Decomposition: " << tpostdecomp - tpostdata << endl <<
                                "  Generating data:  " << tpostdata - tinvert << endl <<
                                "  Inverting matrix:  " << tinvert - tbeg << endl;

  delete [] re;
  delete [] im;
  delete [] ore;
  delete [] oim;
  delete [] ncolloc;
  delete [] mucolloc;
  delete [] phicolloc;
  delete [] lp;
  delete dinfo_p;
}

static void print_coeffs
    (int s, int nl, int nn, double *re, double *im)
{
  const int nlmmodes = nl*(nl+2*ABS(s));
  cout.precision(15);
  for (int n=0; n < nn; n++)
  {
    int sf=0;
    for (int ll=0; ll < nl; ll++)
    {
      const int l = ABS(s) + ll;

      for (int m=-l; m < l+1; m++, sf++)
      {
	cout << "(n=" << n << ", l=" << l << ", m="<<m <<") = (" << re[sf+n*nlmmodes] << ", "<<
	  im[sf+n*nlmmodes] << ")" << endl;
      }
    }
  }
}

static void print_coeffs2
    (int s, int nl, double *re, double *im)
{
  int sfun = 0;
  for (int ll=0; ll < nl; ll++)
  {
    const int l = ABS(s) + ll;
    for (int m=-l; m <l+1; m++, sfun++)
    {
      cout << "(l=" << l << ", m="<<m <<") = (" << re[sfun] << ", "<<
	  im[sfun] << ")" << endl;
    }
  }
}
