 /*@@
   @file      Ricci.c
   @date      Tue May 14 04:38:49 2002
   @author    Ian Hawke
   @desc 
   Routine to calculate the Ricci tensor and scalar. These are taken 
   straight from ADMMacros. 
   @enddesc 
   @version $Header$
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"

#include "ADMAnalysis.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(EinsteinBase_ADMAnalysis_Ricci_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void ADMAnalysis_Ricci(CCTK_ARGUMENTS);

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
   @routine    ADMAnalysis_Ricci
   @date       Tue May 14 21:17:26 2002
   @author     Ian Hawke
   @desc 
   Uses the macros to calculate the Ricci tensor and scalar.
   Unfortunately at the moment it calculates the scalar even
   when you don't want it to.
   @enddesc 
   @calls     ADMAnalysis_Trace
              ADMMacros/src/macro/RICCI_*.h
   @calledby   
   @history 
 
   @endhistory 

@@*/

void ADMAnalysis_Ricci(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS_ADMAnalysis_Ricci;

  CCTK_INT lsh[3], ash[3];
  CCTK_INT i,j,k, d, ijk, di, dj, dk;

#include "EinsteinBase/ADMMacros/src/macro/RICCI_declare.h"

  for (k = 1; k < cctk_lsh[2]-1; k++) 
  {
    for (j = 1; j < cctk_lsh[1]-1; j++) 
    {
      for (i = 1; i < cctk_lsh[0]-1; i++) 
      {
    
        ijk = CCTK_GFINDEX3D(cctkGH,i,j,k);
        di = ijk - CCTK_GFINDEX3D(cctkGH,i-1,j,k);
        dj = ijk - CCTK_GFINDEX3D(cctkGH,i,j-1,k);
        dk = ijk - CCTK_GFINDEX3D(cctkGH,i,j,k-1);

#include "EinsteinBase/ADMMacros/src/macro/RICCI_guts.h"
  
        Ricci11[ijk] = RICCI_RXX;
        Ricci12[ijk] = RICCI_RXY;
        Ricci13[ijk] = RICCI_RXZ;
        Ricci22[ijk] = RICCI_RYY;
        Ricci23[ijk] = RICCI_RYZ;
        Ricci33[ijk] = RICCI_RZZ;
        
      }
    }
  }

#include "EinsteinBase/ADMMacros/src/macro/RICCI_undefine.h"

  for (d=0; d<3; ++d) ash[d] = cctk_ash[d];
  ADMAnalysis_Trace(ash, gxx, gxy, gxz, gyy, gyz, gzz, 
                    Ricci11, Ricci12, Ricci13, Ricci22, Ricci23, Ricci33,
                    Ricci, detg);

  return;
  
}

void ADMAnalysis_Ricci_Boundaries(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_ADMAnalysis_Ricci_Boundaries;

  CCTK_INT err;

  /* Apply Flat Boundary Condition */
  err = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                                  "ADMAnalysis::ricci_scalar", "Flat");
  err += Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                                  "ADMAnalysis::ricci_tensor", "Flat");
  if (err < 0)
  {
    CCTK_WARN(2,"Error in applying flat boundary condition to Ricci tensor");
  }

  /* WARNING:  Only flat boundary conditions are used here.  If the
   *  value of the Ricci tensor is going to be used for something other
   *  than output, then the boundaries should be handled more properly.
   */
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
