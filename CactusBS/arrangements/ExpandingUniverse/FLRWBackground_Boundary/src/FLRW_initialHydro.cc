#include <math.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "FLRWBackground.h"


/*
Populates the hydro and spacetime GFs with a FLRW background.
*/
extern "C"
void initialHydro_FLRW(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  double rhocrit = 3.0*H0*H0/(8.0*PI);

  /* Precalculate the values that will populate the GFs */
  double eps0 = getEps( rhocrit , K_intrinsic , a0 , tol);
  double rho0 = rhocrit/(1+eps0);
  double press0 = eps0*rho0/3.0;
  double sx0;
  double sy0;
  double sz0;
  double dens0;
  double tau0;
  CCTK_REAL vlowx, vlowy, vlowz;
  CCTK_REAL D_h_w,det;

  /* Populate the GFs */
  for (int k = 0; k < cctk_lsh[2]; k++)
  {
    for (int j = 0; j < cctk_lsh[1]; j++)
    {
      for (int i = 0; i < cctk_lsh[0]; i++)
      {
        CCTK_INT ind = CCTK_GFINDEX3D(cctkGH, i , j , k);

        CCTK_INT indx = CCTK_GFINDEX3D(cctkGH, i , j , k );
        CCTK_INT indy = CCTK_GFINDEX3D(cctkGH, i , j , k );
        CCTK_INT indz = CCTK_GFINDEX3D(cctkGH, i , j , k );

        CCTK_REAL R2 = x[ind]*x[ind] + y[ind]*y[ind] + z[ind]*z[ind];
        CCTK_REAL R = pow(R2 , 0.5);
        
        if(FLRW_excision and R < FLRW_excisionRadius) {
          rho[ind] = TINY;
          eps[ind] = eps0;
          press[ind] = TINY*eps0;
          
        } else {
          rho[ind] = rho0;
          eps[ind] = eps0;
          press[ind] = press0;
        }

        
//        w_lorentz[ind] = 1.;
        velx[indx] = 0.;
        vely[indy] = 0.;
        velz[indz] = 0.;
        
        w_lorentz[ind] = 1/sqrt(1.0-(
                                  gxx[ind] * velx[ind] * velx[ind]+
                                  gyy[ind] * vely[ind] * vely[ind]+
                                  gzz[ind] * velz[ind] * velz[ind]+
                                2.0*gxy[ind] * velx[ind] * vely[ind]+
                                2.0*gxz[ind] * velx[ind] * velz[ind]+
                                2.0*gyz[ind] * vely[ind] * velz[ind]));
      SpatialDet(gxx[ind],gxy[ind],gxz[ind],gyy[ind],gyz[ind],gzz[ind],&det);


      dens[ind] = sqrt(det) * w_lorentz[ind] * rho[ind];
      D_h_w = dens[ind] * (1.0 + eps[ind] + press[ind]/rho[ind]) * w_lorentz[ind];

      vlowx = gxx[ind]*velx[ind] + gxy[ind]*vely[ind] + gxz[ind]*velz[ind];
      vlowy = gxy[ind]*velx[ind] + gyy[ind]*vely[ind] + gyz[ind]*velz[ind];
      vlowz = gxz[ind]*velx[ind] + gyz[ind]*vely[ind] + gzz[ind]*velz[ind];

      sx[ind] = D_h_w * vlowx;
      sy[ind] = D_h_w * vlowy;
      sz[ind] = D_h_w * vlowz;

      /* One could use D_h_w here, but it would introduce more error */
      tau[ind] = sqrt(det) * (rho[ind]  *w_lorentz[ind]*w_lorentz[ind] *
                                         (1.0 + eps[ind]) +
                              press[ind]*(w_lorentz[ind]*w_lorentz[ind]-1.0) ) -
                 dens[ind];


      }
    }
  }

}
