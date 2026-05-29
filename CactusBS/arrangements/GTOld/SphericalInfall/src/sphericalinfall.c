#include <stdio.h>

#include "assert.h"
#include "math.h"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#define DEBUG

// scheduled routine
void SphericalInfall_Initial ( CCTK_ARGUMENTS );

// math helpers
static inline double pow2 (const double x)
{
  return x*x;
}
static inline double pow4 (const double x)
{
  const double y = x*x;
  return y*y;
}

#define M 1.0
void SphericalInfall_Initial ( CCTK_ARGUMENTS ) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL Dp, Ep, Vp; // Hawleys D, E and V^r
  CCTK_REAL drho, deps, dpress, dvelx, dvely, dvelz, dw_lorentz; // HydroBase's primitives
  CCTK_REAL Ur, Urb;    // u^r, u^{\bar r}
  CCTK_REAL Ux, Uy, Uz; // isotropic cartesion components
  CCTK_REAL dr_drbar;   // diff(r,rbar) where rbar is the isotropic radial coordinate
  CCTK_REAL rbar, rsw;  // isotropic and Schwarzschild radius
  CCTK_REAL xx,yy,zz;   // isotropic coordinates
  CCTK_REAL g11,g12,g13,g22,g23,g33; // spatial metric components
  CCTK_REAL beta1, beta2, beta3, alpha; // shift and lapse
  CCTK_REAL gtt, grr;   // Schwarzschild metric components

  for ( int k = 0 ; k < cctk_lsh[2] ; k++) {
    for ( int j = 0 ; j < cctk_lsh[1] ; j++) {
      for (int i = 0 ; i < cctk_lsh[0] ; i++) {

        const int ind = CCTK_GFINDEX3D (cctkGH, i, j, k);

        xx = x[ind];
        yy = y[ind];
        zz = z[ind];
        rbar = sqrt(pow2(xx)+pow2(yy)+pow2(zz));
	rbar = pow(pow(rbar,4)+pow(epsilon,4),0.25); // keep away from the puncture
        rsw  = rbar*pow2(1+M/(2*rbar));
        dr_drbar = 1-pow2(M)/(4*pow2(rbar));

        // the metric is always assumed to be Schwarzschild (essentially from Schutz 10.88)
        beta1 = beta2 = beta3 = 0;
        alpha = (1.-M/(2*rbar))/(1.+M/(2*rbar));
        g12 = g13 = g23 = 0.;
        g11 = g22 = g33 = pow4(1+M/(2*rbar));

        if(rbar > rmin)
        {
            // Schwarzschild metric components
            gtt = -(1-2*M/rsw);
            grr = -1/gtt;

            // hydro quantities (Eqs. 67a--67d in ApJ 277:296)
            Dp = initial_dens / ( pow2(rsw) * sqrt( (2*M/rsw) * (1-2*M/rsw) ) );
            Ep = initial_energy / ( pow(sqrt(2*M/rsw)*pow2(rsw), eos_gamma) * pow(1-2*M/rsw,(eos_gamma+1)/4.) );
            Vp = sqrt(2*M/rsw) * (1-2*M/rsw); // hmm, this is positive...

            // get Lorentz factor from -1 = u^mu u_mu, u^t = W/alpha, V^mu = u^mu/u^0
            // in particular Eq. 1a in ApJ 277:296 is wrong 
            // (eg. since V^mu V_mu = (..)^2 u^mu u_mu < 0)
            dw_lorentz = 1 / sqrt(1 - grr * pow2(Vp/alpha));
            assert(dw_lorentz >= 1.);

            // four velocity (ApJ 277:296 1b)
            Ur  = -dw_lorentz/alpha * Vp; // throw in a minus sign to get infalling matter
            Urb = Ur / dr_drbar;
            Ux = Urb * xx/rbar;
            Uy = Urb * yy/rbar;
            Uz = Urb * zz/rbar;
        } 
        else // inside EH, set everything to atmosphere
        {
            Dp = rho_min;
            Ep =  rho_min / initial_dens * initial_energy; // this is just a made up ratio so that E<<D
            Vp = 0.;
            dw_lorentz = 1.;

            Ur = Urb = Ux = Uy = Uz = 0.;
        }

        // set the primitives
        drho = Dp / dw_lorentz;
        deps = Ep / Dp;
        dpress = drho*deps*(eos_gamma-1.);
        dvelx = Ux/dw_lorentz + beta1/alpha;
        dvely = Uy/dw_lorentz + beta2/alpha;
        dvelz = Uz/dw_lorentz + beta3/alpha;

        // fill in the grid functions
        gxx[ind] = g11;
        gxy[ind] = g12;
        gxz[ind] = g13;
        gyy[ind] = g22;
        gyz[ind] = g23;
        gzz[ind] = g33;
        alp[ind] = alpha;
        betax[ind] = beta1;
        betay[ind] = beta2;
        betaz[ind] = beta3;

        rho[ind] = drho;
        velx[ind] = dvelx;
        vely[ind] = dvely;
        velz[ind] = dvelz;
        eps[ind] = deps;
        press[ind] = dpress;

	assert(drho > 0.);
	assert(deps > 0.);
	assert(dpress >= 0.);
	assert(dw_lorentz == dw_lorentz);
	assert(dw_lorentz >= 1.);
	assert(( g11*pow2(dvelx)+g22*pow2(dvely)+g33*pow2(dvelz) +
                          2.*g12*dvelx*dvely + 2.*g13*dvelx*dvelz + 
			  2.*g23*dvely*dvelz ) <= 1. );

      }
    }
  }
}
