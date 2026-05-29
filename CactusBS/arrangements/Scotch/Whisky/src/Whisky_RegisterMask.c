 /*@@
   @file      Whisky_RegisterMask.c
   @date      Sun Jan 26 01:55:25 2003
   @author    Ian Hawke
   @desc 
   Routines to register states with SpaceMask.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "SpaceMask.h"

#include <stdio.h>
#include <stdlib.h>

 /*@@
   @routine    Whisky_RegisterMask
   @date       Sun Jan 26 01:56:06 2003
   @author     
   @desc 
   Register the different mask states with the SpaceMask thorn.
   
   At the moment, the recognized states and values are

   Hydro_Atmosphere (in_atmosphere, not_in_atmosphere)
   Hydro_RiemannProblem (trivial, not_trivial)

   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

int Whisky_RegisterMask(void)
{
  
  DECLARE_CCTK_PARAMETERS;
  
  int ierr;

  const char *atmosphere_list[2] = {"in_atmosphere","not_in_atmosphere"};
  const char *rp_list[2] = {"trivial","not_trivial"};
  const char *alive_list[2] = {"alive","dead"}; /* NB: alive==0 */
  const char *kill_list[2] = {"will_live", "will_die"};
  
  ierr = SpaceMask_RegisterType("Hydro_Atmosphere", 2, atmosphere_list);
  if (ierr)
  {
    CCTK_WARN(0, "Failed to register the atmosphere with the mask!");
  }

  ierr = SpaceMask_RegisterType("Hydro_RiemannProblemX", 2, rp_list);
  if (ierr)
  {
    CCTK_WARN(0, "Failed to register the x Riemann Problem with the mask!");
  }

  ierr = SpaceMask_RegisterType("Hydro_RiemannProblemY", 2, rp_list);
  if (ierr)
  {
    CCTK_WARN(0, "Failed to register the y Riemann Problem with the mask!");
  }

  ierr = SpaceMask_RegisterType("Hydro_RiemannProblemZ", 2, rp_list);
  if (ierr)
  {
    CCTK_WARN(0, "Failed to register the z Riemann Problem with the mask!");
  }

  ierr = SpaceMask_RegisterType("Hydro_AliveCells", 2, alive_list);
  if (ierr)
  {
    CCTK_WARN(0, "Failed to register the dead/alive indicators with the mask!");
  }

  ierr = SpaceMask_RegisterType("Hydro_KillCells", 2, kill_list);
  if (ierr)
  {
    CCTK_WARN(0, "Failed to register the kill/letlive indicators with the mask!");
  }

  return 0;
}

void Whisky_SetupDescriptors(CCTK_ARGUMENTS)
{
  
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  *atmosphere_field_descriptor=SpaceMask_GetTypeBits("Hydro_Atmosphere");
  *atmosphere_atmosp_descriptor=SpaceMask_GetStateBits("Hydro_Atmosphere",
                                                    "in_atmosphere");
  *atmosphere_normal_descriptor  = SpaceMask_GetStateBits("Hydro_Atmosphere",
                                                    "not_in_atmosphere");
}

