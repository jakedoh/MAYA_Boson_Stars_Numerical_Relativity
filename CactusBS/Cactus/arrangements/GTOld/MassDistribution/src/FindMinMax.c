#include <assert.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "MassDistribution.h"

#define HUGE_VAL 1e60

void MassDistribution_InitMinMax(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int vol;
  int do_any_calc, my_calc_every;
 
  /* check if any of the volumes need computation */
  do_any_calc = 0;
  for (vol=0 ; vol<intvolumes; vol++)
  {
    my_calc_every = calc_every_vol[vol] >= 0 ? calc_every_vol[vol] : calc_every;
    assert(my_calc_every >= 0);
    HistogramActive[vol] = my_calc_every > 0 && 
                           cctk_iteration % my_calc_every == 0 &&
                           calc_after_vol[vol] <= cctk_iteration &&
                           cctk_iteration <= dont_calc_after_vol[vol];
    HistogramValid[vol] = 1;
    if (HistogramActive[vol] && HistogramValid[vol])
    {
      if(CCTK_Equals(bin_mode[vol], "automatic") || 
         CCTK_Equals(bin_mode[vol], "automatic-volume") ||
         CCTK_Equals(bin_mode[vol], "mass-fraction"))
      {
        /* initialize the persistent arrays for the local loops */
        if(bin_variable_minimum[vol] == bin_variable_maximum[vol])
        {
          HistogramMinimum[vol] = +HUGE_VAL;
          HistogramMaximum[vol] = -HUGE_VAL;
        }
        else
        {
          HistogramMinimum[vol] = bin_variable_minimum[vol];
          HistogramMaximum[vol] = bin_variable_maximum[vol];
        }
        
        do_any_calc = 1;
      }
      else if(CCTK_Equals(bin_mode[vol], "manual"))
      {
        /* always set min/max of manual mode bins even if there is nothing else to do */
        HistogramMinimum[vol] = bin_variable_minimum[vol]; 
        HistogramMaximum[vol] = bin_variable_maximum[vol]; 
      }
      else
      {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Unknown range mode '%s' for volume %d.", bin_mode[vol], vol);
        return; /* NOTREACHED */
      }
    }
  }

  if (do_any_calc) 
  {
    if (verbose >= 2)
      CCTK_VInfo(CCTK_THORNSTRING,"Finding min/max at iteration %d",cctk_iteration);

    /* kickstart the while loop in schedule.ccl */
    *have_work_to_do = 1;
  }
  else
  {
    *have_work_to_do = 0;
  }
}
/* find min/max in a single patch and processor*/
void MassDistribution_FindMinMax(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int do_every;
 
  /* first check if we should actually run in this refinement level
   * this works around a weakness in Carpet's global loop-local tags that would
   * call us on each and every existing level instead of only the ones a local
   * function would be called on */
  MassDistribution_get_do_level_every(cctkGH, &do_every);
  if(cctk_iteration % do_every != 0)
    return;
 
  MassDistribution_sumLocalMass(CCTK_PASS_CTOC, find_minmax);
}
     
/* reduce minima and maxima from all processors */
void MassDistribution_ReduceMinMax(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int vol, ierr;
  int num_auto_vols;
  /* note the order of the indices, we have packets of (min,max) for each volume */
  CCTK_REAL local_minmax_buffer[intvolumes][2];
  CCTK_REAL global_minmax_buffer[intvolumes][2];

  /* get the min/max for volumes that are determined automatically */
  num_auto_vols = 0;
  for(vol = 0 ; vol < intvolumes ; vol++)
  {
    if(HistogramActive[vol] && HistogramValid[vol])
    {
      if(CCTK_Equals(bin_mode[vol], "automatic") || 
         CCTK_Equals(bin_mode[vol], "automatic-volume") ||
         CCTK_Equals(bin_mode[vol], "mass-fraction"))
      {
          local_minmax_buffer[num_auto_vols][0] = -HistogramMinimum[vol]; 
          local_minmax_buffer[num_auto_vols][1] = +HistogramMaximum[vol]; 
          num_auto_vols += 1;
      }
      else if(CCTK_Equals(bin_mode[vol], "manual"))
      {
          /* do nothing */
      }
      else
      {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Unknown range mode '%s' for volume %d.", bin_mode[vol], vol);
        return; /* NOTREACHED */
      }
    }
  }

  /* do all reductions at once to reduce communication overhead */
  const int max_op = CCTK_ReductionHandle ("maximum");
  if (max_op < 0)
    CCTK_WARN (CCTK_WARN_ABORT, "error when obtaining 'maximum' reduction handle");
  ierr = CCTK_ReduceLocalArray1D (cctkGH, -1, max_op, local_minmax_buffer,
          global_minmax_buffer, 2*num_auto_vols, CCTK_VARIABLE_REAL);
  if (ierr)
    CCTK_WARN (CCTK_WARN_ABORT, "error performing 'maximum' reduction");

  /* unpack and merge min/max information into static array again */
  num_auto_vols = 0;
  for(vol = 0 ; vol < intvolumes ; vol++)
  {
    if(HistogramActive[vol] && HistogramValid[vol])
    {
      if(CCTK_Equals(bin_mode[vol], "automatic") || 
         CCTK_Equals(bin_mode[vol], "automatic-volume") ||
         CCTK_Equals(bin_mode[vol], "mass-fraction"))
      {
        HistogramMinimum[vol] = -global_minmax_buffer[num_auto_vols][0]; 
        HistogramMaximum[vol] = +global_minmax_buffer[num_auto_vols][1]; 
        /* make a bit of room at the end bins to be sure not to miss anything
         * due to roundoff or because many values are actually at the found
         * maximum */
        HistogramMinimum[vol] -= 0.001*(HistogramMaximum[vol]-HistogramMinimum[vol])/number_of_bins;
        HistogramMaximum[vol] += 0.001*(HistogramMaximum[vol]-HistogramMinimum[vol])/number_of_bins;
        num_auto_vols += 1;
      }
      else if(CCTK_Equals(bin_mode[vol], "manual"))
      {
        /* do nothing */
      }
      else
      {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Unknown range mode '%s' for volume %d.", bin_mode[vol], vol);
        return; /* NOTREACHED */
      }

      if (verbose >= 2) {
        CCTK_VInfo(CCTK_THORNSTRING,"Using range [%g,%g] for volume %d",
                   HistogramMinimum[vol], HistogramMaximum[vol], vol);
      }
    }
  }
}
