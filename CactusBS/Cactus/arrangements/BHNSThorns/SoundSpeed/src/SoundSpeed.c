
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "EOS_Base.h"

#define SQR(x) ((x) * (x))

#if ET_HYDROBASE
#define velx (&vel[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely (&vel[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz (&vel[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#endif

#define INITVALUE (42)

/* External Functions */
void recalc_soundspeed (CCTK_ARGUMENTS);
void soundspeed_init (CCTK_ARGUMENTS);


void soundspeed_init (CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_LOOP3_ALL( soundspeed_initgfs, cctkGH, i,j,k )
  {
    CCTK_INT idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
    SoundSpeed[idx]=INITVALUE;
    MachNum[idx]=INITVALUE;
  } CCTK_ENDLOOP3_ALL( soundspeed_initgfs ); 

}

void recalc_soundspeed (CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (verbose > 0) {
     CCTK_INFO("Entering Speed of Sound Calculation");
  }

  /* NB: I should technically also check if eos_handle_var is a CCTK_INT or
   * not, but there are only so many things I am willing to check */
  CCTK_INT *eos_handle;
  eos_handle = (CCTK_INT *)CCTK_VarDataPtr(cctkGH, 0, eos_handle_var);
  if (NULL == eos_handle) {
     CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                "invalid eos_handle_var '%s'", eos_handle_var);
     return; /* NOTREACHED */
  }

  CCTK_LOOP3_ALL( soundspeed_setGFs, cctkGH, i,j,k ) {
//  #pragma omp parallel for firstprivate(eos_handle)
//  for (int k = 0; k < cctk_lsh[2]; k++)
//    for (int j = 0; j < cctk_lsh[1]; j++)
//      for (int i = 0; i < cctk_lsh[0]; i++)

      CCTK_REAL rhoL, epsL, velxL, velyL, velzL;
      CCTK_REAL pressL, dpdrho, dpdeps, cs2, enthalpy;
      CCTK_REAL gxxL, gxyL, gxzL, gyyL, gyzL, gzzL, v2;
      CCTK_REAL soundspeedL, machL;

      rhoL = INITVALUE;
      epsL = INITVALUE;
      velxL = INITVALUE;
      velyL = INITVALUE;
      velzL = INITVALUE;
    
      gxxL = INITVALUE;
      gxyL = INITVALUE;
      gxzL = INITVALUE;
      gyyL = INITVALUE;
      gyzL = INITVALUE;
      gzzL = INITVALUE;
    
      dpdrho = INITVALUE;
      dpdeps = INITVALUE;
    
      // Load Gridfunctions
      CCTK_INT index  =  CCTK_GFINDEX3D(cctkGH,i,j,k);
    
    
      rhoL = rho[index];
      epsL = eps[index];
      pressL = press[index];
    
      velxL = velx[index];
      velyL = vely[index];
      velzL = velz[index];
    
      gxxL = gxx[index];
      gxyL = gxy[index];
      gxzL = gxz[index];
      gyyL = gyy[index];
      gyzL = gyz[index];
      gzzL = gzz[index];
    
      // Main calculation
      dpdrho = EOS_DPressByDRho(*eos_handle,rhoL,epsL);
      dpdeps = EOS_DPressByDEps(*eos_handle,rhoL,epsL);
      enthalpy = 1. + epsL + pressL/rhoL;
    
      cs2 = (dpdrho + pressL*dpdeps/(rhoL*rhoL)) / enthalpy;
    
      v2 = gxxL*SQR(velxL) + gyyL*SQR(velyL) + gzzL*SQR(velzL)
    	+ 2.*gxyL*velxL*velyL + 2.*gxzL*velxL*velzL + 2.*gyzL*velyL*velzL;
    
      // Checks on the physicality
      if (cs2 < 0 || cs2 > 1.) {
    	CCTK_VInfo(CCTK_THORNSTRING,"Speed of sound is unphysical at (%g,%g,%g)!!",
                              x[index],y[index],z[index]);
      }
      soundspeedL = sqrt(cs2);
      
      // Fill in the gridfunction 
      SoundSpeed[index] = soundspeedL;
      MachNum[index] = ( cs2 <0 || cs2 > 1 ? -1 : sqrt(v2)/soundspeedL );

  } CCTK_ENDLOOP3_ALL( soundspeed_setGFs );


} 
