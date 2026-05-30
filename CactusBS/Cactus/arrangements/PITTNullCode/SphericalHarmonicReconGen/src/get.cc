#include <stdio.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "vars.hh"


extern "C" 
{

/* Cactus Aliased Functions */

  CCTK_INT SphericalHarmonicReconGeneric_GetParameters (CCTK_INT *lmax,
                                            CCTK_REAL *r_inner,
                                            CCTK_REAL *r_outer)
  {
    // we don't have Rin or Rout, so we always return an error.
    // This routine is not called by any method within the Nullcode anyway
    return -1;
  }



  CCTK_INT SphericalHarmonicReconGeneric_GetCurrentCoefs(
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
    if (!SHR::read_data)
    {
      CCTK_WARN(CCTK_WARN_ABORT, "Error: data has not been read in");
      return -1;
    }

    double factor = 1;
    if (SHR::use_Condon_Shortley_phase_factor) 
       factor = pow(-1.0, m);

    for (int c = 0; c < NUM_METRIC_COMPONENTS; c++)
    {
      reC[c] = factor * SHR::C[c]->get_mode(l, m).real(); 
      imC[c] = factor * SHR::C[c]->get_mode(l, m).imag(); 
      reCr[c] = factor * SHR::Cr[c]->get_mode(l, m).real(); 
      imCr[c] = factor * SHR::Cr[c]->get_mode(l, m).imag();
      if (SHR::time_derivative_in_file) { 
         reCt[c] = factor * SHR::Ct[c]->get_mode(l, m).real(); 
         imCt[c] = factor * SHR::Ct[c]->get_mode(l, m).imag();
      } else {
         reCt[c] = factor * SHR::C[c]->get_mode_dt(l, m).real(); 
         imCt[c] = factor * SHR::C[c]->get_mode_dt(l, m).imag();
      }
    }

    return 0;
  }

}

