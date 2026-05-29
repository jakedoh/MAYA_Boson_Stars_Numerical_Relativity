#include "stdio.h"

#include "cctk.h" 
#include "cctk_Arguments.h" 
#include "cctk_Parameters.h" 

void Psi4Analysis_Integrate_Psi4 (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (*do_nothing==1) {
    return;
  }

  if (*rdet<=0) {
    CCTK_VInfo(CCTK_THORNSTRING,"no 2D integral computations, rdet=%g current_detector=%d",
	       *rdet,*current_detector);
    return;
  }

  if (verbose>2) {
    CCTK_INFO("compute time integrals on 2D surface");
  }

  CCTK_REAL dt=cctk_time- lastt[*current_detector];

  if (verbose>3) {
    CCTK_VInfo(CCTK_THORNSTRING,
               "Integrate_Psi4.c: *current_detector=%d",
               *current_detector);
    CCTK_INT npoints=*ntheta* *nphi;
    CCTK_VInfo(CCTK_THORNSTRING,
               "Integrate_Psi4.c: npoints=%d dt=%g cctk_time=%g lastt=%g",
               npoints,dt,cctk_time,lastt[*current_detector]);
  }

  // loop over the sphere and update time integrals
  CCTK_INT ind2d=0;
  CCTK_INT ind3d=0;
  for (int i=0;i<*ntheta;i++) {
    for (int j=0;j<*nphi;j++) {
      ind2d=i+ maxntheta*j;
      ind3d=i+ maxntheta*j+ maxntheta*maxnphi*(*current_detector);

      int_psi4r_dt[ind3d]+=tetrad_factor*psi4_2d_re[ind2d]*dt;
      int_psi4i_dt[ind3d]+=tetrad_factor*psi4_2d_im[ind2d]*dt;

      int2_psi4r_dt2[ind3d]+=int_psi4r_dt[ind3d]*dt;
      int2_psi4i_dt2[ind3d]+=int_psi4i_dt[ind3d]*dt;
    }
  }
}
