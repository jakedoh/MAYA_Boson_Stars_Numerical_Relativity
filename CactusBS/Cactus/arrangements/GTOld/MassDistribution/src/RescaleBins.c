#include <assert.h>
#include <float.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "MassDistribution.h"

// standard helper macros
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

void MassDistribution_RescaleBins(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL total_mass, partial_mass;
  /* the 99% range is from the left edge of first_bin to the right edge of
   * last_bin, centered around peak_bin*/ 
  CCTK_INT first_bin, last_bin, peak_bin;
  CCTK_INT more_work;

  more_work = 0;
  for(int vol = 0 ; vol < intvolumes ; vol++)
  {
    if(HistogramActive[vol] && HistogramValid[vol]) 
    {
      if(CCTK_Equals(bin_mode[vol], "mass-fraction"))
      {
        const CCTK_REAL *mass = &MassPerSpecificEnergy[vol*TOTAL_NUMBER_OF_BINS]; 

        /* sumLocalMasses stored the total mass for us */
        total_mass = mass[TOTAL_MASS_BIN];

        /* now find a range around the maximum such that it covers
         * mass_fraction[vol] of the mass */
        peak_bin = 0;
        for(int bin = 1 ; bin < number_of_bins ; bin++)
        {
          if(mass[bin] > mass[peak_bin])
            peak_bin = bin;
        }
        first_bin = last_bin = peak_bin;
        partial_mass = mass[peak_bin];
        while(partial_mass < mass_fraction[vol]*total_mass)
        {
          if(!(first_bin > 0 || last_bin < number_of_bins-1))
          {
            CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                       "Unexpected error: bin range cannot not grow [%g,%g], "
                       "in volume %d, current partial mass is %g, "
                       "total_mass is %g", 
                       HistogramMinimum[vol], HistogramMaximum[vol], vol,
                       partial_mass, total_mass);
            HistogramValid[vol] = 0;
            break; /* we need a further if outside to go to the next volume */
          }

          /* move first_bin to the left if it contains more mass than last_bin+1
             or if last_bin is the very last bin already
             same for the last_bin as well */
          if(first_bin > 0 && 
            (last_bin == number_of_bins-1 || mass[first_bin-1] >= mass[last_bin+1]))
          {
            first_bin -= 1;
            partial_mass += mass[first_bin];
          }
          else if(last_bin < number_of_bins-1 && 
            (first_bin == 0 || mass[first_bin-1] < mass[last_bin+1]))
          {
            last_bin += 1;
            partial_mass += mass[last_bin];
          }
          else
          {
            CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                       "Internal error: bin range did not grow [%d,%d], "
                       "candidates were {%g,%g}, current partial mass is %g, "
                       "total_mass is %g, volume is %d", first_bin, last_bin, 
                       mass[MAX(first_bin-1,0)], 
                       mass[MIN(last_bin+1, number_of_bins-1)],
                       partial_mass, total_mass, vol);
            return; /* NOTREACHED */
          }
          if(verbose >= 4)
          {
            CCTK_VInfo(CCTK_THORNSTRING,"Testing %g%% range [%g,%g] for volume %d, fractional mass: %g%%.",
                       mass_fraction[vol]*100.0,
                       HistogramLabels[vol * NUMBER_OF_LABELS + first_bin], 
                       HistogramLabels[vol * NUMBER_OF_LABELS + last_bin + 1], vol,
                       partial_mass/total_mass*100.);
          }
        }
        /* skip further processing if we could not find a valid range for some
         * reason */
        if(!HistogramValid[vol])
          continue;
        if(!(first_bin <= last_bin))
        {
          CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Scaling failed in volume %d, min/max [%g,%g], "
                     "first/last [%d,%d], total_mass is %g", 
                     vol, HistogramLabels[vol * NUMBER_OF_LABELS + first_bin], 
                     HistogramLabels[vol * NUMBER_OF_LABELS + last_bin + 1],
                     first_bin, last_bin, total_mass);
          HistogramValid[vol] = 0;
          continue;
        }

        if(last_bin-first_bin+1 < acceptable_fraction_of_bins[vol]*number_of_bins)
        {
          /* give some information of what we are doing */
          if(verbose >= 2 || (*have_work_to_do == max_number_of_iterations && verbose >= 1))
          {
            CCTK_VInfo(CCTK_THORNSTRING,"Rejecting %g%% range [%g,%g] for volume %d (only %g%% of bins); rescaling.",
                       mass_fraction[vol]*100.0,
                       HistogramLabels[vol * NUMBER_OF_LABELS + first_bin], 
                       HistogramLabels[vol * NUMBER_OF_LABELS + last_bin + 1], 
                       vol, (last_bin-first_bin+1.)/number_of_bins*100.0);
          }

          /* set new min/max values */
          HistogramMinimum[vol] = HistogramLabels[vol * NUMBER_OF_LABELS + first_bin];
          HistogramMaximum[vol] = HistogramLabels[vol * NUMBER_OF_LABELS + last_bin + 1];

          if(!(HistogramMaximum[vol] >= HistogramMinimum[vol]))
          {
            CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                       "Empty histogram in volume %d, min/max [%g,%g], "
                       "first/last [%d,%d], total_mass is %g", 
                       vol, HistogramLabels[vol * NUMBER_OF_LABELS + first_bin], 
                       HistogramLabels[vol * NUMBER_OF_LABELS + last_bin + 1],
                       first_bin, last_bin, total_mass);
            HistogramValid[vol] = 0;
            continue;
          }

          more_work = 1; /* set mark for scheduler to re-run the integration */
        }
        else
        {
          /* accept values, let user know */
          if(verbose >= 2)
          {
            CCTK_VInfo(CCTK_THORNSTRING,"Accepting %g%% range [%g,%g] for volume %d",
                       mass_fraction[vol]*100.0,
                       HistogramMinimum[vol], HistogramMaximum[vol], vol);
          }
        }
      } // bin_mode == "mass-fraction"
    } // HistogramActive && HistogramValid
  }

  if(more_work)
  {
    if(*have_work_to_do < max_number_of_iterations) /* try again */
    {
      *have_work_to_do += 1;
    }
    else /* give up */
    {
      *have_work_to_do = 0;
      /* warn and inform user about why we fail */
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Maximum number of iterations %d exceeded when rescaling the histogram.", max_number_of_iterations);
    }
  }
  else /* all done */
  {
    *have_work_to_do = 0;
  }
}
