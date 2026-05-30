#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "util_Table.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>

/* This routine should only be called at cctk_iteration=start_tracing_particles_iteration */
void InitializeRK4IterationCounterToOne(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(cctk_iteration!=start_tracing_particles_iteration) return;

  *RK4IterationCounter = 1;

  if(verbose==2) printf("particle_tracerET: Just set RK4IterationCounter to %d\n",*RK4IterationCounter);
}
