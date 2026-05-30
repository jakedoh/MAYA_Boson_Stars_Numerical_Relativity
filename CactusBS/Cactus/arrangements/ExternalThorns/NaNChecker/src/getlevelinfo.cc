/* $Header$ */

#include <cstdio>
#include <cmath>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_DefineThorn.h"

#include "carpet.hh"

namespace NaNChecker {

#ifdef HAVE_CARPET
  using namespace Carpet;
#endif

  extern "C"
    void NaNChecker_levelinfo(const cGH * cctkGH,
                              CCTK_INT * my_level,
                              CCTK_INT * n_levels)
    {

      *my_level = 0;
      *n_levels = 1;

#ifdef HAVE_CARPET

      if ( CCTK_IsThornActive ( "Carpet" ) ) {
        *my_level = reflevel;
        *n_levels = reflevels;
      } 
     
#endif

    }

  extern "C"
    CCTK_FCALL void CCTK_FNAME (NaNChecker_levelinfo) (CCTK_POINTER_TO_CONST * cctkGH, CCTK_INT * my_level, CCTK_INT * n_levels)
    {
      NaNChecker_levelinfo((const cGH *)* cctkGH,  my_level, n_levels);
    }

} // namespace NaNChecker
