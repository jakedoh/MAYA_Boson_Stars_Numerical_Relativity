#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"


#include "cctk.h" 
#include "cctk_Arguments.h" 
#include "cctk_Parameters.h" 


// functions called
void ADM_EJP_ProjectSphere(CCTK_ARGUMENTS);
void ADM_EJP_Compute_EJP(CCTK_ARGUMENTS);


// this function is only necessary because the fucking schedule.ccl
// loop mechanism seems to be broken!
void ADM_EJP_loop_over_detectors (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (number_of_detectors==0) {
    CCTK_INFO("no detectors requested");
  }

  for (int i=0;i<number_of_detectors;i++)
  {
    if (out_det_ev[i]<=0) { 
      CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
                 "skip this detector: out_det_ev[%d]=%d <=0.",i,out_det_ev[i]);
      continue;
    }
    if (cctk_iteration % out_det_ev[i] != 0) {
      if (verbose>1) {
        CCTK_VInfo(CCTK_THORNSTRING,
                   "wrong iteration: cctk_iteraation=%d, out_det_ev[%d]=%d",
                   cctk_iteration,i,out_det_ev[i]);
      }
      continue;
    }

    *current_detector=i;

    if (verbose>1) {
      CCTK_VInfo(CCTK_THORNSTRING,"detector=%d",i);
    }

    if (verbose>2) {
      CCTK_INFO("calling ProjectSphere");
    }

    ADM_EJP_ProjectSphere(CCTK_PASS_CTOC);

    if (verbose>2) {
      CCTK_INFO("calling Compute_EJP");
    }

    ADM_EJP_Compute_EJP(CCTK_PASS_CTOC);

    if (verbose>2) {
      CCTK_INFO("done with this detector");
    }

    lastt[i]=cctk_time;
  }
}
