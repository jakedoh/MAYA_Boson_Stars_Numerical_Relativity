/* $Header$ */

#include <cstdio>
#include <cmath>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_DefineThorn.h"

#include "carpet.hh"

namespace Whisky {

#ifdef HAVE_CARPET
  using namespace Carpet;
#endif

  extern "C"
    int Whisky_AllowPastTimelevels(void)
    {
      int retval = 0;

#ifdef HAVE_CARPET

      if ( CCTK_IsThornActive ( "Carpet" ) ) {
        retval = do_allow_past_timelevels;
      } 
     
#endif

      return retval;

    }

  extern "C"
    CCTK_FCALL void CCTK_FNAME (Whisky_AllowPastTimelevels) (int *retval)
    {
      *retval = Whisky_AllowPastTimelevels();
    }

} // namespace Whisky
