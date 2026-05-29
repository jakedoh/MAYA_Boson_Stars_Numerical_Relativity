 /*@@
   @file      RegisterMask.c
   @desc 
   Routines to register states with SpaceMask.
   Copied from Whisky.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "SpaceMask.h"

namespace AHMask {

extern "C" void AHMask_RegisterMask(CCTK_ARGUMENTS)
{
  
  DECLARE_CCTK_PARAMETERS;
  
  int ierr;
  // the ordering below should make sure that a zeroed out mask is in an all outside mask
  const char * const ahmask_list[3] = {space_mask_outside_value,space_mask_inside_value,space_mask_buffer_value};
  
  // SpaceMask is not smart enough to detect if we want to register a type that
  // already exists, so we are checking for this situation here
  if (SpaceMask_GetTypeBits(mask_bitfield_name) != 0)
  {
    CCTK_WARN(2, "AHFinder type mask already registered, not registering it again.");
    return;
  }    

  ierr = SpaceMask_RegisterType(mask_bitfield_name, sizeof(ahmask_list)/sizeof(ahmask_list[0]), ahmask_list);
  if (ierr)
  {
    CCTK_WARN(0, "Failed to register the AHFinder type mask with the mask!");
  }

}

}
