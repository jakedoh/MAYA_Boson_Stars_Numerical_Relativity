/*  file produced by user hinder, 15/11/2006 */
/*  Produced with Mathematica Version 5.2 for Linux (June 20, 2005) */

/*  Mathematica script written by Ian Hinder and Sascha Husa */

/*  $Id$ */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Faces.h"
#include "util_Table.h"
#include "Symmetry.h"


/* the boundary treatment is split into 3 steps:    */
/* 1. excision                                      */
/* 2. symmetries                                    */
/* 3. "other" boundary conditions, e.g. radiative */

/* to simplify scheduling and testing, the 3 steps  */
/* are currently applied in separate functions      */


void WeylScal4_CheckBoundaries(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  return;
}

void WeylScal4_ApplyBoundConds(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_INT ierr = 0;

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                                   "WeylScal4::Psi4r_group", "none");
  if (ierr < 0)
    CCTK_WARN(-1, "Failed to register Psi4r_group_bound BC for WeylScal4::Psi4r_group!");
  
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1, 
                                   "WeylScal4::Psi4i_group", "none");
  if (ierr < 0)
    CCTK_WARN(-1, "Failed to register Psi4i_group_bound BC for WeylScal4::Psi4i_group!");
  
  return;
}



/* template for entries in parameter file:
*/

