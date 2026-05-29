#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"


#include "cctk.h" 
#include "cctk_Arguments.h" 
#include "cctk_Parameters.h" 


// functions called
void scalarAnalysis_ProjectSphere(CCTK_ARGUMENTS);
void scalarAnalysis_Compute_EPRad(CCTK_ARGUMENTS);


int scalarAnalysis_setup_out_every_det(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (verbose>2) {
    CCTK_INFO("setup out_every per detector");
  }

  CCTK_INT ierr=0;
  char default_out_every[2048]; // FIXME XXX let's hope that we never need more
  sprintf(default_out_every,"%d",out_every);
  

  // initialize detectors which are not set to out_every
  for (int i=0;i<number_of_detectors;i++) {
    if (out_every_det[i]==-1) {
      char paramstr[2048]; // FIXME XXX let's hope that we never need more
      sprintf(paramstr,"out_every_det[%d]",i);

      ierr=CCTK_ParameterSet(paramstr,"scalarAnalysis",default_out_every);
      if (ierr<0) {
	CCTK_WARN(1,"parameter set for detector failed!");
	return -1;
      }
    }
    if (verbose>1) {
      CCTK_VInfo(CCTK_THORNSTRING,"out_every_det[%d]=%d",
		 i,out_every_det[i]);
    }
  }

  return 1;
}


// this function is only necessary because the fucking schedule.ccl
// loop mechanism seems to be broken!
void scalarAnalysis_loop_over_detectors (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (number_of_detectors==0) {
    CCTK_INFO("no detectors requested");
  }

  CCTK_INT ierr=scalarAnalysis_setup_out_every_det(CCTK_PASS_CTOC);
  if (ierr<0) {
    CCTK_WARN(1,"failed to setup detector out_every");
    return;
  }

  for (int i=0;i<number_of_detectors;i++) {
    *current_detector=i;

    if (verbose>1) {
      CCTK_VInfo(CCTK_THORNSTRING,"detector=%d",i);
    }

    // check if it is time to do something
    if (cctk_iteration % out_every_det[i] != 0) {
      if (verbose>1) {
	CCTK_VInfo(CCTK_THORNSTRING,
		   "wrong iteration: cctk_iteration=%d, out_every_det=%d, detector=%d",
		   cctk_iteration,out_every_det[i],i);
      }
      continue; // jump to next detector
    }
    else if (verbose>1) {
      CCTK_VInfo(CCTK_THORNSTRING,
		 "start computations at cctk_iteration=%d, out_every_det=%d, detector=%d",
		 cctk_iteration,out_every_det[i],i);
    }

    // setup theta and phi integrals
    *ntheta=maxntheta;
    *nphi=maxnphi;
    if (ntheta_det[i]>maxntheta) {
      CCTK_WARN(1,"ntheta_det larger than maxntheta, using maxntheta instead");
      *ntheta=maxntheta;
    }
    if (nphi_det[i]>maxnphi) {
      CCTK_WARN(1,"nphi_det larger than maxnphi, using maxnphi instead");
      *nphi=maxnphi;
    }

    if (ntheta_det[i]>0) {
      *ntheta=ntheta_det[i];
      if (verbose>1) {
	CCTK_VInfo(CCTK_THORNSTRING,"setting ntheta to %d",*ntheta);
      }
    }

    if (nphi_det[i]>0) {
      *nphi=nphi_det[i];
      if (verbose>1) {
	CCTK_VInfo(CCTK_THORNSTRING,"setting nphi to %d",*nphi);
      }
    }

    if (verbose>2) {
      CCTK_VInfo(CCTK_THORNSTRING,"ntheta=%d nphi=%d",*ntheta,*nphi);
    }

    if (verbose>2) {
      CCTK_INFO("calling ProjectSphere");
    }

    scalarAnalysis_ProjectSphere(CCTK_PASS_CTOC);

    if (verbose>2) {
      CCTK_INFO("calling Compute_EPRad");
    }

    scalarAnalysis_Compute_EPRad(CCTK_PASS_CTOC);

    if (verbose>2) {
      CCTK_INFO("done with this detector");
    }
  }
}
