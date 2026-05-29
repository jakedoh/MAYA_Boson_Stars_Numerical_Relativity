 /*@@
   @file      Symmetry.c
   @date      Tue Jul 30 11:15:23 CEST 2002
   @author    David Rideout
   @desc 
   Symmetry registration stuff for ADMAnalysis
   @enddesc
   @version $Header$
 @@*/



#include "cctk.h"

#include "cctk_Arguments.h"
#include "Symmetry.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusEinstein_ADMAnalysis_Symmetry_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/
void ADMAnalysis_RegisterSymmetry(CCTK_ARGUMENTS);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************

 ********************************************************************/

 /*@@
   @routine    ADMAnalysis_RegisterSymmetry
   @date       Tue Jul 30 11:31:50 CEST 2002
   @author     David Rideout
   @desc 
   Scheduled routine to register symmetry of Ricci.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
 @@*/
void ADMAnalysis_RegisterSymmetry(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_ADMAnalysis_RegisterSymmetry;
  int sym[3];
  int ierr;

  /*CCTK_INFO("Registering symmetries for ricci_scalar and ricci_tensor");*/
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  ierr=SetCartSymVN(cctkGH, sym, "ADMAnalysis::trK");
  if(ierr)
    CCTK_VWarn (0, __LINE__, __FILE__,"ADMAnalysis", "Error returned from function SetCartSymVN");
  ierr=SetCartSymVN(cctkGH, sym, "ADMAnalysis::detg");
  if(ierr)
    CCTK_VWarn (0, __LINE__, __FILE__,"ADMAnalysis", "Error returned from function SetCartSymVN");
  ierr=SetCartSymVN(cctkGH, sym, "ADMAnalysis::Ricci");
  if(ierr)
    CCTK_VWarn (0, __LINE__, __FILE__,"ADMAnalysis", "Error returned from function SetCartSymVN");
  ierr=SetCartSymVN(cctkGH, sym, "ADMAnalysis::Ricci11");
  if(ierr)
    CCTK_VWarn (0, __LINE__, __FILE__,"ADMAnalysis", "Error returned from function SetCartSymVN");
  ierr=SetCartSymVN(cctkGH, sym, "ADMAnalysis::Ricci22");
  if(ierr)
    CCTK_VWarn (0, __LINE__, __FILE__,"ADMAnalysis", "Error returned from function SetCartSymVN");
  ierr=SetCartSymVN(cctkGH, sym, "ADMAnalysis::Ricci33");
  if(ierr)
    CCTK_VWarn (0, __LINE__, __FILE__,"ADMAnalysis", "Error returned from function SetCartSymVN");
  sym[0] = -1; sym[1] = -1; sym[2] = 1;
  ierr=SetCartSymVN(cctkGH, sym, "ADMAnalysis::Ricci12");
  if(ierr)
    CCTK_VWarn (0, __LINE__, __FILE__,"ADMAnalysis", "Error returned from function SetCartSymVN");
  sym[0] = -1; sym[1] = 1; sym[2] = -1;
  ierr=SetCartSymVN(cctkGH, sym, "ADMAnalysis::Ricci13");
  if(ierr)
    CCTK_VWarn (0, __LINE__, __FILE__,"ADMAnalysis", "Error returned from function SetCartSymVN");
  sym[0] = 1; sym[1] = -1; sym[2] = -1;
  ierr=SetCartSymVN(cctkGH, sym, "ADMAnalysis::Ricci23");
  if(ierr)
    CCTK_VWarn (0, __LINE__, __FILE__,"ADMAnalysis", "Error returned from function SetCartSymVN");
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
