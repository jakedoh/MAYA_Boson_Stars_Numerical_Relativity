 /*@@
   @file      Register.c
   @date      6 May 2003
   @author    David Rideout
   @desc 
              Register implemented boundary conditions.
   @enddesc 
   @version   $Header$
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "SampleBnd.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusExamples_SampleBoundary_Register_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Aliased Routine Prototypes ***********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void SampleBoundary_RegisterBCs(CCTK_ARGUMENTS);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     Aliased Routines   ***********************
 ********************************************************************/

/********************************************************************
 *********************     Scheduled Routines   *********************
 ********************************************************************/

 /*@@
   @routine    SampleBoundary_RegisterBCs
   @date       6 May 2003
   @author     David Rideout
   @desc 
               Register all boundary conditions implemented by this thorn.
   @enddesc 
   @calls      
   @history 
   @endhistory
   @var        CCTK_ARGUMENTS
   @vdesc      Cactus argument list
   @vtype      CCTK_*
   @vio        in
   @endvar
   @returntype void
@@*/

void SampleBoundary_RegisterBCs(CCTK_ARGUMENTS)
{
  int err;  
  
  err = Boundary_RegisterPhysicalBC((CCTK_POINTER) cctkGH, 
                                    (phys_bc_fn_ptr) &BndLinExtrap,
                                    "LinExtrap");
  if (err)
  {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Error %d when registering routine to handle \"LinExtrap\" "
                 "boundary condition", err);
  }
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
