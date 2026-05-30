#include <cstdio>
#include <cmath>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#ifdef HAVE_CARPET
#include "carpet.hh"
#endif

namespace EHFinder {

#ifdef HAVE_CARPET
  using namespace Carpet;

    void loop_over_levels (CCTK_ARGUMENTS, CCTK_INT current_level)
    {
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS

      CCTK_INT ierr;
      int i,j,k,l,index;

      // Do the work.
      enter_level_mode ( cctkGH, current_level );

        ierr = CallScheduleGroup ( cctkGH, "Euler_ReInitializeEvolve" );
      leave_level_mode ( cctkGH );

      // Only count iterations on the finest level.
      if ( current_level == reflevels-1 ) {
        (*niter_reinit)++;
      }


      // If necessary go to the next level.
      if ( current_level < reflevels-1) {
        loop_over_levels (CCTK_PASS_CTOC,current_level+1);
      }  

      // If we are not on the coarsest level, it's time to take the next step.
      if ( current_level > 0 ) {

        enter_level_mode ( cctkGH, current_level );
          ierr = CallScheduleGroup ( cctkGH, "Euler_ReInitializeEvolve" );
          Restrict ( cctkGH );
        leave_level_mode ( cctkGH );

        // Again only count iterations on the finest level.
        if ( current_level == reflevels-1 ) {
          (*niter_reinit)++;
        }

        // If necessary go to the next level.
        if ( current_level < reflevels-1 ) {
          loop_over_levels (CCTK_PASS_CTOC,current_level+1);
        }

        // It's time to exit.
        return;
      }
  
    }
#endif

  extern "C"
    void EHFinder_PreReInitialize_Carpet(CCTK_ARGUMENTS)
    {
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS

#ifdef HAVE_CARPET
      my_level = reflevel;
      leave_level_mode ( cctkGH ); 

      // If we are called on the finest level call the EHFinder_PreReInitialize
      // schedule group otherwise make sure that no re-initialization routines
      // are executed.

      if ( my_level == reflevels - 1 ) {
        ierr = CallScheduleGroup ( cctkGH, "EHFinder_PreReInitialize" );
      } else {
        *re_init_control = 0;
      }
      enter_level_mode ( cctkGH, my_level );
#endif
    }

  extern "C"
    void EHFinder_ReInitialize_Wrapper(CCTK_ARGUMENTS)
    {
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS

#ifdef HAVE_CARPET
      // Get the current level.
      my_level = reflevel;

      // Leave level mode.
      leave_level_mode ( cctkGH ); 

      // If we were called on the finest level then loop recursively over all
      // levels.
      if ( my_level == reflevels - 1 ) {

        loop_over_levels (CCTK_PASS_CTOC,current_level);

        // Check if the re-initialization is completed.
        ierr = CallScheduleGroup ( cctkGH, "Euler_PostStep" );
      }

      // Go back to level mode.
      enter_level_mode ( cctkGH, my_level );

#endif
    }
} // namespace EHFinder
