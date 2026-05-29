#include <cstdio>
#include <cmath>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_DefineThorn.h"

#include "carpet.hh"

namespace MassDistribution {

#ifdef HAVE_CARPET
  using namespace Carpet;
#endif

  // from CarpetLibs/defs.cc
  template <typename T>
  inline T ipow_helper (T x, unsigned int y)
  {
    T z = y&1 ? x : 1;
    while (y >>= 1)
    {
      x *= x;
      if (y & 1) z *= x;
    }
    return z;
  }


  extern "C"
    void MassDistribution_get_do_level_every(const cGH * cctkGH, CCTK_INT * do_every)
    {

      *do_every = 1;

#ifdef HAVE_CARPET

      if ( CCTK_IsThornActive ( "Carpet" ) ) {
        *do_every
          = ipow(mgfact, mglevel) * (maxtimereflevelfact / timereffacts.at(reflevel));
      } 
     
#endif

    }

  extern "C"
    CCTK_FCALL void CCTK_FNAME (MassDistribution_get_do_level_every) (CCTK_POINTER_TO_CONST * cctkGH, CCTK_INT * do_every)
    {
      MassDistribution_get_do_level_every((const cGH *)* cctkGH,  do_every);
    }

} // namespace MassDistribution
