 /*@@
   @header    Startup.c
   @date      
   @author    
   @desc 
   Register known elliptic interfaces
   @enddesc
   @version $Header$
 @@*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "CactusElliptic/EllBase/src/EllBase.h"
#include "CactusElliptic/EllBase/src/Ell_DBstructure.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusElliptic_EllSOR_Startup_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

void EllSOR_Register(CCTK_ARGUMENTS);
void SORConfMetric(cGH *GH, 
                   int *MetricPsiI, 
                   int FieldI, 
                   int MI,
                   int NI, 
                   CCTK_REAL *AbsTol, 
                   CCTK_REAL *RelTol);
void SORMetric(cGH *GH, 
               int *MetricI, 
               int FieldI, 
               int MI,
               int NI, 
               CCTK_REAL *AbsTol, 
               CCTK_REAL *RelTol);
void SORFlat(cGH *GH, 
             int FieldI,
             int MI, 
             int NI, 
             CCTK_REAL *AbsTol, 
             CCTK_REAL *RelTol);


/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void EllSOR_Register(CCTK_ARGUMENTS);
 
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
  @routine    EllSOR_Register
  @date       
  @author     
  @desc 
  Scheduled routine to register the SOR solvers SORConfMetric and 
  SORFlatunder with the name "sor" for the Elliptic Classes 
  LinConfMetric and LinFlat
  @enddesc 
  @calls     
  @calledby   
  @history 
  
  @endhistory 
  
 @@*/
void EllSOR_Register(CCTK_ARGUMENTS) 
{

  DECLARE_CCTK_ARGUMENTS 
  DECLARE_CCTK_PARAMETERS 

  int err;

  if (Ell_RegisterSolver(SORConfMetric,"sor","Ell_LinConfMetric")
        != ELL_SUCCESS)
  {
    CCTK_WARN(0,
              "EllSOR_Register: Failed to register sor for Ell_LinConfMetric");
  }

  if (Ell_RegisterSolver(SORMetric,"sor","Ell_LinMetric")
        != ELL_SUCCESS)
  {
    CCTK_WARN(0,
              "EllSOR_Register: Failed to register sor for Ell_LinMetric");
  }

  if (Ell_RegisterSolver(SORFlat,"sor","Ell_LinFlat") != ELL_SUCCESS)
  {
    CCTK_WARN(0,"EllSOR_Register: Failed to register sor for Ell_LinFlat");
  }

  /* These "keys" have to be same in other elliptic solvers and in 
     the routines that sets them ! */

  /* Register boundaries which SOR can handle */
  err = Ell_CreateKey(CCTK_VARIABLE_STRING,"EllLinFlat::Bnd::Robin");
  if (err == ELLCREATE_FAILED)
  {
    CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
               "EllSOR_Register: Failed to create key EllLinFlat::Bnd::Robin (Error %d)",err);
  }
  err = Ell_CreateKey(CCTK_VARIABLE_STRING,"EllLinFlat::Bnd::Const");
  if (err == ELLCREATE_FAILED)
  {
    CCTK_WARN(1,
              "EllSOR_Register: Failed to create key EllLinFlat::Bnd::Const");
  }

  /* Create a key for the maximum number of iterations allowed. */
  err = Ell_CreateKey(CCTK_VARIABLE_INT, "Ell::SORmaxit");
  if (err == ELLCREATE_FAILED)
  {
    CCTK_WARN(0, "EllSOR_Register: Failed to create key Ell::SORmaxit");
  } 

  /* Create a key for the type of acceleration to be used. */
  err = Ell_CreateKey(CCTK_VARIABLE_STRING, "Ell::SORaccel");
  if (err == ELLCREATE_FAILED)
  {
    CCTK_WARN(0, "EllSOR_Register: Failed to create key Ell::SORaccel");
  } 

}
