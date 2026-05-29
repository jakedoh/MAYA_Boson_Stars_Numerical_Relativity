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

  /* Populate the GFs */
  for (int k = 0; k < cctk_lsh[2]; k++)
  {
    for (int j = 0; j < cctk_lsh[1]; j++)
    {
      for (int i = 0; i < cctk_lsh[0]; i++)
      {
        CCTK_INT ind = CCTK_GFINDEX3D(cctkGH, i , j , k);

        CCTK_INT indx = CCTK_VECTGFINDEX3D(cctkGH, i , j , k , 0);
        CCTK_INT indy = CCTK_VECTGFINDEX3D(cctkGH, i , j , k , 1);
        CCTK_INT indz = CCTK_VECTGFINDEX3D(cctkGH, i , j , k , 2);

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

        vel[indx] = 0.;
        vel[indy] = 0.;
        vel[indz] = 0.;
        w_lorentz[ind] = 1.0;
      }
    }
  }

}
