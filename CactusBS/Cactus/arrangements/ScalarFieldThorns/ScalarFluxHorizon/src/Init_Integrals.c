#include "cctk.h" 
#include "cctk_Arguments.h" 
#include "cctk_Parameters.h" 

void ScalarFluxHorizon_Init_Integrals (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;


  const CCTK_INT npoints=maxntheta*maxnphi*number_of_detectors;

  // init last time
  for (int i=0;i<number_of_detectors;i++) {
    lastt[i]=0.;
  }

  // files don't exist initially
  *ScalarFluxHorizon_files_exist=0;

  // initialize integral values
  for (int i=0;i<number_of_detectors;i++)
  {
    energy[i]=0.;
    Px[i]=0.;
    Py[i]=0.;
    Pz[i]=0.;
    Jz[i]=0.;
  }

  if (verbose>0) {
    CCTK_INFO("initialized integrals");
  }

  if (out_offset!=0) {
    CCTK_WARN(1,"out_offset is deprecated, assuming 0");
  }

  *do_nothing=0;
}
