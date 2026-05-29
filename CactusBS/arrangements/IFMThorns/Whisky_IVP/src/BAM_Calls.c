 /*@@
   @file      BAM_Calls.c
   @date      Sun Jul  7 14:26:46 2002
   @author    Ian Hawke
   @desc 
   Interfaces to BAM.
   @enddesc 
   @version  $Header$
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "NaNChecker.h"

#include "CactusPUGH/PUGH/src/include/pugh.h"
#include "AEIThorns/BAM_Elliptic/src/bbmg.h"

static const char *rcsid = "$Header: /EUNetwork/Whisky/Whisky_IVP/"
                           "src/BAM_Calls.c,v 1.1.1.1 2002/07/07 "
                           "11:58:15 hawke Exp $";

CCTK_FILEVERSION(Whisky_IVP_BAM_Calls_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

void BAM_Elliptic_VecL(cGH *cctkGH, int *nu, int *nN, CCTK_REAL tol);
void BAM_Elliptic_Powers(cGH *cctkGH, pGH *GH, int *u, double *p,
                         int n, CCTK_REAL tol);

void Call_BAM_VecLap(cGH *cctkGH, CCTK_REAL Elliptic_Tolerance,
                                  CCTK_INT IVP_Check_NaNs);
void Call_BAM_Nonlin(cGH *cctkGH, CCTK_REAL Elliptic_Tolerance,
                                  CCTK_INT IVP_Check_NaNs);

void CCTK_FCALL CCTK_FNAME(Call_BAM_VecLap)
                          (cGH **cctkGH, CCTK_REAL *Elliptic_Tolerance,
                                         CCTK_INT *IVP_Check_NaNs);
void CCTK_FCALL CCTK_FNAME(Call_BAM_Nonlin)
                          (cGH **cctkGH, CCTK_REAL *Elliptic_Tolerance,
                                         CCTK_INT *IVP_Check_NaNs);

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

void IVP_NaNChecker_before_BAM(cGH *cctkGH)
{
  NaNChecker_CheckVarsForNaN(cctkGH,-1,"Whisky_IVP::IVP_Vx","both","abort");
  NaNChecker_CheckVarsForNaN(cctkGH,-1,"Whisky_IVP::IVP_Vy","both","abort");
  NaNChecker_CheckVarsForNaN(cctkGH,-1,"Whisky_IVP::IVP_Vz","both","abort");
  NaNChecker_CheckVarsForNaN(cctkGH,-1,"Whisky_IVP::IVP_Vx_Source",
                                       "both","abort");
  NaNChecker_CheckVarsForNaN(cctkGH,-1,"Whisky_IVP::IVP_Vy_Source",
                                       "both","abort");
  NaNChecker_CheckVarsForNaN(cctkGH,-1,"Whisky_IVP::IVP_Vz_Source",
                                       "both","abort");

  NaNChecker_CheckVarsForNaN(cctkGH,-1,"Whisky_IVP::IVP_uxx",
                                       "both","abort");
  NaNChecker_CheckVarsForNaN(cctkGH,-1,"Whisky_IVP::IVP_uyy",
                                       "both","abort");
  NaNChecker_CheckVarsForNaN(cctkGH,-1,"Whisky_IVP::IVP_uzz",
                                       "both","abort");
  NaNChecker_CheckVarsForNaN(cctkGH,-1,"Whisky_IVP::IVP_uxy",
                                       "both","abort");
  NaNChecker_CheckVarsForNaN(cctkGH,-1,"Whisky_IVP::IVP_uyz",
                                       "both","abort");
  NaNChecker_CheckVarsForNaN(cctkGH,-1,"Whisky_IVP::IVP_uxz",
                                       "both","abort");
  
  NaNChecker_CheckVarsForNaN(cctkGH,-1,"ADMBase::gxx",
                                       "both","abort");
  NaNChecker_CheckVarsForNaN(cctkGH,-1,"ADMBase::gyy",
                                       "both","abort");
  NaNChecker_CheckVarsForNaN(cctkGH,-1,"ADMBase::gzz",
                                       "both","abort");
  NaNChecker_CheckVarsForNaN(cctkGH,-1,"ADMBase::gxy",
                                       "both","abort");
  NaNChecker_CheckVarsForNaN(cctkGH,-1,"ADMBase::gyz",
                                       "both","abort");
  NaNChecker_CheckVarsForNaN(cctkGH,-1,"ADMBase::gxz",
                                       "both","abort");
  
/*  NaNChecker_CheckVarsForNaN(cctkGH,-1,"StaticConformal::Psi",
                                       "both","abort");
*/
}

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    Call_BAM_VecLap
   @date       Sun Jul  7 14:40:28 2002
   @author     Ian Hawke
   @desc 
   Interface to BAM. Calls the vector laplace solver to solve
   the momentum constraints.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void Call_BAM_VecLap(cGH *cctkGH, CCTK_REAL Elliptic_Tolerance,
                                  CCTK_INT IVP_Check_NaNs)
{

  CCTK_INT vars[3], sources[3];
  
  if (IVP_Check_NaNs)
  {
     IVP_NaNChecker_before_BAM(cctkGH);
  }

  vars[0] = CCTK_VarIndex("Whisky_IVP::IVP_Vx");
  vars[1] = CCTK_VarIndex("Whisky_IVP::IVP_Vy");
  vars[2] = CCTK_VarIndex("Whisky_IVP::IVP_Vz");
  
  sources[0] = CCTK_VarIndex("Whisky_IVP::IVP_Vx_Source");
  sources[1] = CCTK_VarIndex("Whisky_IVP::IVP_Vy_Source");
  sources[2] = CCTK_VarIndex("Whisky_IVP::IVP_Vz_Source");
  
  BAM_Elliptic_VecL(cctkGH, vars, sources, Elliptic_Tolerance);

  return;

}

 /*@@
   @routine    Call_BAM_Nonlin
   @date       Sun Jul  7 14:41:10 2002
   @author     Ian Hawke
   @desc 
   Interface to BAM. Calls the nonlinear power solver to solve
   the energy constraint.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void Call_BAM_Nonlin(cGH *cctkGH, CCTK_REAL Elliptic_Tolerance,
                                  CCTK_INT IVP_Check_NaNs)
{
  CCTK_INT vars[6];
  CCTK_REAL powers[6];

  pGH *GH = (pGH *) cctkGH->extensions [CCTK_GHExtensionHandle ("PUGH")];

  if (IVP_Check_NaNs)
  {
     IVP_NaNChecker_before_BAM(cctkGH);
  }
  
  vars[0] = CCTK_VarIndex("Whisky_IVP::IVP_Psi");
  vars[1] = CCTK_VarIndex("Whisky_IVP::IVP_Mlinear");
  vars[2] = CCTK_VarIndex("Whisky_IVP::IVP_Nsource");
  vars[3] = CCTK_VarIndex("Whisky_IVP::IVP_Psi5Source");
  vars[4] = CCTK_VarIndex("Whisky_IVP::IVP_Psim7Source");
  vars[5] = CCTK_VarIndex("Whisky_IVP::IVP_Psim3Source");
  
  powers[0] = 0.0;
  powers[1] = 0.0;
  powers[2] = 0.0;
  powers[3] = 5.0;
  powers[4] = -7.0;
  powers[5] = -3.0;
  
  BAM_Elliptic_Powers(cctkGH, GH, vars, powers, 6, Elliptic_Tolerance);

  return;

}

 /*@@
   @routine    Call_BAM_*
   @date       Sun Jul  7 14:42:16 2002
   @author     Ian Hawke
   @desc 
   Fortran wrappers for the above routines.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void CCTK_FCALL CCTK_FNAME(Call_BAM_VecLap)
                          (cGH **cctkGH, CCTK_REAL *Elliptic_Tolerance,
                                         CCTK_INT *IVP_Check_NaNs)
{
  Call_BAM_VecLap(*cctkGH, *Elliptic_Tolerance, *IVP_Check_NaNs);
  return;
}

void CCTK_FCALL CCTK_FNAME(Call_BAM_Nonlin)
                          (cGH **cctkGH, CCTK_REAL *Elliptic_Tolerance,
                                         CCTK_INT *IVP_Check_NaNs)
{
  Call_BAM_Nonlin(*cctkGH, *Elliptic_Tolerance, *IVP_Check_NaNs);
  return;
}


/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

