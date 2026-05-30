/* $Header$ */

#include <cstdio>
#include <cmath>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_DefineThorn.h"

#include "carpet.hh"

namespace WeylScal4 {

#ifdef HAVE_CARPET
  using namespace Carpet;
#endif

  extern "C"
    CCTK_INT WeylScal4_rnd_iteration(CCTK_INT it)
    {
    
      CCTK_INT retval = it;

#ifdef HAVE_CARPET
      // integer time step size for this level
      CCTK_INT idt = maxtimereflevelfact/timereflevelfact;

      // round iteration upward to nearest multiple of timereflevelfact
      if ( CCTK_IsThornActive ( "Carpet" ) ) {
        if ( it % idt == 0 ) {
          retval = it;
        } else {
          retval = it - it % idt + idt;
        }
      } 
     
#endif

      return retval;
    }

  extern "C"
    CCTK_FCALL CCTK_INT CCTK_FNAME (WeylScal4_rnd_iteration) (CCTK_INT * pit)
    {
      return WeylScal4_rnd_iteration(*pit);
    }

} // namespace WeylScal4
