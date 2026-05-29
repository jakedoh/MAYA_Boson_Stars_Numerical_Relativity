#include <assert.h>
#include <float.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "MassDistribution.h"

void MassDistribution_InitMasses(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int vol;

  /* check if any of the volumes need computation */
  *have_work_to_do = 0;
  for (vol=0 ; vol<intvolumes; vol++)
  {
    if(HistogramActive[vol] && HistogramValid[vol])
    {
      if (verbose >= 2)
        CCTK_VInfo(CCTK_THORNSTRING,"Calculating mass distribution mass at iteration %d",cctk_iteration);
      /* kickstart the while loop in schedule.ccl */
      *have_work_to_do = 1;
      break;
    }
  }
}

void MassDistribution_ZeroMasses(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int vol;
  enum bin_scaling_type scaling;

  /* initialize the persistent arrays (incl. total mass) for the local loops */
  for (int i=0; i<intvolumes*TOTAL_NUMBER_OF_BINS; i++) {
    MassPerSpecificEnergy[i] = 0.0;
  }

  for (vol=0 ; vol<intvolumes; vol++)
  {
    scaling = MassDistribution_resolve_keyword_parameter("bin_scaling", bin_scaling[vol], NUM_BIN_SCALINGS, bin_scalings);
    switch (scaling) /* linear! size of bin */
    {
      case bin_scaling_linear:
        for(int bin=0; bin<NUMBER_OF_LABELS; bin++)
          HistogramLabels[vol*NUMBER_OF_LABELS + bin] = HistogramMinimum[vol] + bin*(HistogramMaximum[vol]-HistogramMinimum[vol])/number_of_bins;
        break;
      case bin_scaling_logarithmic:
        for(int bin=0; bin<NUMBER_OF_LABELS; bin++)
          HistogramLabels[vol*NUMBER_OF_LABELS + bin] = HistogramMinimum[vol]*pow(HistogramMaximum[vol]/HistogramMinimum[vol], (double)bin/number_of_bins);
        break;
      default:
        assert(0 == "Internal error");
        break;
    }
  }
}


/* add up mass in a single patch and accumulate to processor total */
void MassDistribution_SumLocalRestMass(CCTK_ARGUMENTS)
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
 
  MassDistribution_sumLocalMass(CCTK_PASS_CTOC, integrate);
}
     
/* sum masses from all processors */
void MassDistribution_SumGlobalRestMass(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int ierr;
  CCTK_REAL *massbuffer_send, *massbuffer_recv;
  int num_active_volumes;
  

  /* some temp. storage */
  massbuffer_send = malloc(sizeof(CCTK_REAL) * intvolumes*TOTAL_NUMBER_OF_BINS);
  assert(massbuffer_send != NULL);
  massbuffer_recv = malloc(sizeof(CCTK_REAL) * intvolumes*TOTAL_NUMBER_OF_BINS);
  assert(massbuffer_recv != NULL);

  /* Get the summation operator */
  const int sum_op = CCTK_ReductionHandle ("sum");
  if (sum_op < 0)  CCTK_WARN (CCTK_WARN_ABORT, "error when obtaining 'sum' reduction handle");

  /* copy data from active histograms into transfer buffer */
  /* TODO: check if using something like ReduceLocalArrays would be more useful */
  num_active_volumes = 0; 
  for (int vol=0 ; vol<intvolumes; vol++) {
    if(HistogramActive[vol] && HistogramValid[vol]) 
    {
      for(int bin = 0 ; bin<TOTAL_NUMBER_OF_BINS ; bin++) {
        massbuffer_send[num_active_volumes*TOTAL_NUMBER_OF_BINS + bin] = 
          MassPerSpecificEnergy[vol*TOTAL_NUMBER_OF_BINS + bin]; 
      }
      num_active_volumes += 1;
    }
  }

  /* do all (incl. total mass) reductions at once to reduce communication overhead */
  /* BCast sum to all processors in case we need to rescale things */
  /* TODO: think about Bcasting just the rescaled min/max after rescaling */
  ierr = CCTK_ReduceLocalArray1D (cctkGH, -1, sum_op,
      massbuffer_send, massbuffer_recv,
      num_active_volumes*TOTAL_NUMBER_OF_BINS, CCTK_VARIABLE_REAL);
  assert(!ierr);
  if (ierr) CCTK_WARN (CCTK_WARN_ABORT, "error performing 'sum' reduction");

  /* copy back into Histogram and apply scalings */
  num_active_volumes = 0; 
  for(int vol = 0 ; vol < intvolumes ; vol++) {
    if(HistogramActive[vol] && HistogramValid[vol]) 
    {
      for(int bin = 0 ; bin<TOTAL_NUMBER_OF_BINS ; bin++) {
        MassPerSpecificEnergy[vol*TOTAL_NUMBER_OF_BINS + bin] = 
          massbuffer_recv[num_active_volumes*TOTAL_NUMBER_OF_BINS + bin];
      }
      /* total mass */
      MassPerSpecificEnergy[vol*TOTAL_NUMBER_OF_BINS + TOTAL_MASS_BIN] = 
        massbuffer_recv[num_active_volumes*TOTAL_NUMBER_OF_BINS + TOTAL_MASS_BIN];
      /* mass above and below threshold */
      MassPerSpecificEnergy[vol*TOTAL_NUMBER_OF_BINS + BELOW_MASS_BIN] = 
        massbuffer_recv[num_active_volumes*TOTAL_NUMBER_OF_BINS + BELOW_MASS_BIN];
      MassPerSpecificEnergy[vol*TOTAL_NUMBER_OF_BINS + ABOVE_MASS_BIN] = 
        massbuffer_recv[num_active_volumes*TOTAL_NUMBER_OF_BINS + ABOVE_MASS_BIN];
      num_active_volumes += 1;
    }
  }

  free(massbuffer_send);
  free(massbuffer_recv);
}
