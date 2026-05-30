 /*@@
   @file      MaskBounds.c
   @date      Tue Aug 25 19:19:48 EDT 2009
   @author    Roland Haas
   @desc 
   Register boundary conditions with thorn Boundary, based on Ian Hinder's and
   Sascha Husa's Kranc's code 
   @enddesc 
 @@*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "cctk.h"

#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "Symmetry.h"
#include "Boundary.h"
#include "SpaceMask.h"


/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void MaskBounds(CCTK_ARGUMENTS);

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

/*@@
   @routine    MaskBounds
   @date       Tue Aug 25 19:19:48 EDT 2009
   @author     Roland Haas
   @desc 
   Register boundary conditions with thorn Boundary, based on Ian Hinder's and 
   Sascha Husa's Kranc's code (for symmetry boundaries)
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
@@*/

void MaskBounds(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  
  CCTK_INT ierr = 0;
  
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                    "SpaceMask::mask", "none");
  if (ierr < 0)
     CCTK_WARN(0, "Failed to register 'none' BC for SpaceMask::mask!");
  
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                    "SpaceMask::space_mask_group", "none");
  if (ierr < 0)
     CCTK_WARN(0, "Failed to register 'none' BC for SpaceMask::space_mask_group!");

  return;
}
