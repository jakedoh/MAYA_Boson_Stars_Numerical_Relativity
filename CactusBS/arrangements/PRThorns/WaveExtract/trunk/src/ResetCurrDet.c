 /*@@
   @file      ResetCurrDet.c
   @date      16 Jul 2003
   @author    Frank Herrmann
   @desc 
              Reset the current_detector value, this is needed for the while loop in
              the scheduler.
   @enddesc 
   @version $Id: ResetCurrDet.c 37 2008-02-14 03:30:05Z schnetter $
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(HerrmannCVS_WaveExtract_ResetCurrDet_c)

void WavExtr_ResetCurrDet(CCTK_ARGUMENTS);

 /*@@
   @routine    WavExtr_Reset_CurrDet
   @date       16 Jul 2003
   @author     Frank Herrmann
   @desc 
               Reset the current_detector value, this is needed for the while loop in
               the scheduler.
   @enddesc 
   @calls     
@@*/
void WavExtr_ResetCurrDet(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (verbose>0)
    CCTK_VInfo(CCTK_THORNSTRING,"Calling WaveExtract at time: %f",(double)cctkGH->cctk_time);

  if (verbose>2)
    CCTK_INFO("Reset the current_detector value");

  if (verbose >3)
    CCTK_VInfo(CCTK_THORNSTRING,"iteration %d detector %d out_every_det[] %d",cctk_iteration,*current_detector,out_every_det[*current_detector]);

  if (cctk_iteration != 0) {
    if (cctk_iteration % out_every_det[*current_detector] != 0)
    {
      *do_nothing=1;
    }
    else if (Cauchy_time_ID!=0) {
      if((cctk_iteration-1) % out_every_det[*current_detector] == 0) {
        *do_nothing=0;
      }
      else if ((cctk_iteration+1) % out_every_det[*current_detector] == 0) {
        *do_nothing=0;
      }
    }

    if (cctk_iteration % out_every_det[*current_detector] == 0) {
      *do_nothing=0;
    }
  }
  else
    *do_nothing = 0;

  if (*do_nothing == 1) {
    return;
  }

  *current_detector=maximum_detector_number;

}
