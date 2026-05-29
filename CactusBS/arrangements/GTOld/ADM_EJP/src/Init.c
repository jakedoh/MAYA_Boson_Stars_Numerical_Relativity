#include "cctk.h" 
#include "cctk_Arguments.h" 
#include "cctk_Parameters.h" 

void ADM_EJP_Init (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;


  // files don't exist initially
  *adm_ejp_files_exist=0;

  // initialize integral values
  for (int i=0;i<number_of_detectors;i++)
  {
    energy[i]=0.;
    Px[i]=0.;
    Py[i]=0.;
    Pz[i]=0.;
    Jx[i]=0.;
    Jy[i]=0.;
    Jz[i]=0.;
    Erad[i]=0;
    Pxrad[i]=0.;
    Pyrad[i]=0.;
    Pzrad[i]=0.;
    Sxrad[i]=0.;
    Syrad[i]=0.;
    Szrad[i]=0.;

    lastt[i]=0.;

    if (out_every_det[i]==-1) {
      out_det_ev[i]=out_every;
    }
    else {
      out_det_ev[i]=out_every_det[i];
    }
  }

  if (verbose>0) {
    CCTK_INFO("initialization done");
  }

  if (out_offset !=0) {
    CCTK_WARN(1,"out_offset parameter is deprecated - automatically assume 0");
  }

  *do_nothing=0;
}
