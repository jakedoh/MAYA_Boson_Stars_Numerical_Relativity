#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void Ylm_ParamCheck(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  char valstr[3048]; // FIXME XXX let's hope that we never need more
  CCTK_INT ierr;



  CCTK_INFO("Checking parameters");
  
  if (ntheta != -1) {
    CCTK_VInfo(CCTK_THORNSTRING,"update maxntheta to: %d",ntheta);
    sprintf(valstr,"%d",ntheta);
    ierr=CCTK_ParameterSet("maxntheta","Ylm_Decomp",valstr);
    if (ierr<0) {
      CCTK_WARN(1,"unable to update maxntheta!");
    }
  }

  if (nphi != -1) {
    CCTK_VInfo(CCTK_THORNSTRING,"update maxnphi to: %d",nphi);
    sprintf(valstr,"%d",nphi);
    ierr=CCTK_ParameterSet("maxnphi","Ylm_Decomp",valstr);
    if (ierr<0) {
      CCTK_WARN(1,"unable to update maxnphi!");
    }
  }

}
