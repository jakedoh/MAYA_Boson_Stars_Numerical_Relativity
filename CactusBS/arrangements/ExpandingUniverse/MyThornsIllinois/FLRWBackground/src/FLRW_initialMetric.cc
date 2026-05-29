
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
void initialMetric_FLRW(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* Precalculate the values that will populate the GFs */
  double g0 = a0*a0;
  double k0 = -H0;

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

        gxx[ind] = g0;
        gxy[ind] = 0.;
        gxz[ind] = 0.;
        gyy[ind] = g0;
        gyz[ind] = 0.;
        gzz[ind] = g0;

        kxx[ind] = k0;
        kxy[ind] = 0.;
        kxz[ind] = 0.;
        kyy[ind] = k0;
        kyz[ind] = 0.;
        kzz[ind] = k0;

        alp[ind] = a0;
        betax[ind] = 0.;
        betay[ind] = 0.;
        betaz[ind] = 0.;
      }
    }
  }

}
