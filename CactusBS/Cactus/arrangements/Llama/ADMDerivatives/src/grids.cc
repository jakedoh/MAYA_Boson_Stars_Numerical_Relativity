#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "util_String.h"
#include "util_Table.h"  

#include "Symmetry.h"
#include "GlobalDerivative.h"

//#include "util.h"
#include "util_String.h"
#include "util_Table.h"

#include "grids.hh"

// Copied from CTGBase/src/util.c

void
ADMDerivatives_GetPhysicalBoundaries(cGH *cctkGH, int* do_bc)
{
//  const cGH* const cctkGH = (const cGH*) GH;  
//  const cGH* cctkGH = (const cGH*) GH;  
  // const cGH *cctkGH = (const cGH *) GH;


  DECLARE_CCTK_ARGUMENTS;

//  int patch = MultiPatch_GetMap(cctkGH);

  if (CCTK_IsFunctionAliased("MultiPatch_GetBbox"))
    {
      CCTK_INT bbox[6];
      MultiPatch_GetBbox(cctkGH, 6, bbox);
      for (int i=0; i<6; ++i)
	do_bc[i] = bbox[i];
    }
  else
    {
      CCTK_INT symbnd[6];

      CCTK_INT symtable = SymmetryTableHandleForGrid(cctkGH);
      if (symtable < 0)
	CCTK_WARN(0, "Could not get symmetry table for grid.");

      int ierr = Util_TableGetIntArray(symtable, 6, symbnd, "symmetry_handle");
      if (ierr != 6)
	CCTK_WARN(0, "Could not get symmetry table handle for grid.");
      for (int i=0; i<6; ++i)
	do_bc[i] = (cctk_bbox[i] && symbnd[i]<0);
    }

  return;
}


CCTK_INT ADMDerivatives_GetGridRanges(cGH *cctkGH, int* istart,
                                      int* iend)
{
  // const cGH* const cctkGH = GH;  
  DECLARE_CCTK_ARGUMENTS;
  int i;
  int do_bc[6];
  
  if (!CCTK_IsFunctionAliased("MultiPatch_GetBbox"))
  {
      ADMDerivatives_GetPhysicalBoundaries (cctkGH, do_bc);
      
      for (i=0; i<3; ++i)
         {
            if (do_bc[2*i])
              istart[i] = 0;
            else
              istart[i] = cctk_nghostzones[i];
            if (do_bc[2*i+1])
              iend[i] = cctk_lsh[i];
            else
              iend[i] = cctk_lsh[i] - cctk_nghostzones[i];
         }
      
      return 0;
  }
  
  CCTK_INT bbox[6];
  CCTK_INT nboundaryzones[6];
  CCTK_INT is_internal[6];
  CCTK_INT is_staggered[6];
  CCTK_INT shiftout[6];
  CCTK_INT symbnd[6];
  CCTK_INT is_ipbnd[6];
  
  CCTK_INT dir = 0;
  CCTK_INT symtable = 0;
  CCTK_INT face = 0;
  CCTK_INT npoints = 0;
  CCTK_INT iret = 0;
  CCTK_INT ierr = 0;
  
  
  if (CCTK_IsFunctionAliased ("MultiPatch_GetBbox")) {
    ierr = MultiPatch_GetBbox (cctkGH, 6, bbox);
    if (ierr != 0)
      CCTK_WARN(0, "Could not obtain bbox specification");
  } else {
    for (dir = 0; dir < 6; dir++)
    {
      bbox[dir] = 0;
    }
  }

  if (CCTK_IsFunctionAliased ("MultiPatch_GetBoundarySpecification")) {
    CCTK_INT const map = MultiPatch_GetMap (cctkGH);
    if (map < 0)
      CCTK_WARN(0, "Could not obtain boundary specification");
    ierr = MultiPatch_GetBoundarySpecification
      (map, 6, nboundaryzones, is_internal, is_staggered, shiftout);
    if (ierr != 0)
      CCTK_WARN(0, "Could not obtain boundary specification");
  } else if (CCTK_IsFunctionAliased ("GetBoundarySpecification")) {
    ierr = GetBoundarySpecification
      (6, nboundaryzones, is_internal, is_staggered, shiftout);
    if (ierr != 0)
      CCTK_WARN(0, "Could not obtain boundary specification");
  } else {
    CCTK_WARN(0, "Could not obtain boundary specification");
  }

  symtable = SymmetryTableHandleForGrid(cctkGH);
  if (symtable < 0)
  {
    CCTK_WARN(0, "Could not obtain symmetry table");
  }

  iret = Util_TableGetIntArray(symtable, 6, symbnd, "symmetry_handle");
  if (iret != 6) CCTK_WARN (0, "Could not obtain symmetry information");

  for (dir = 0; dir < 6; dir++)
  {
    is_ipbnd[dir] = (!cctk_bbox[dir]);
  }

  for (dir = 0; dir < 3; dir++)
  {
    for (face = 0; face < 2; face++)
    {
      CCTK_INT index = dir*2 + face;
      if (is_ipbnd[index])
      {
        /* Inter-processor boundary */
        npoints = cctk_nghostzones[dir];
      }
      else
      {
        /* Symmetry or physical boundary */
        npoints = nboundaryzones[index];

      }

      switch(face)
      {
      case 0: /* Lower boundary */
        istart[dir] = npoints;
        break;
      case 1: /* Upper boundary */
        iend[dir] = cctk_lsh[dir] - npoints;
        break;
      default:
        CCTK_WARN(0, "internal error");
      }
    }
  }


  return 0;
}
