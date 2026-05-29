/*  file produced by user hinder, 15/11/2006 */
/*  Produced with Mathematica Version 5.2 for Linux (June 20, 2005) */

/*  Mathematica script written by Ian Hinder and Sascha Husa */

/*  $Id$ */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

void WeylScal4_RegisterSymmetries(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  
  /* array holding symmetry definitions */
  CCTK_INT sym[3];
  
  
  /* Register symmetries of grid functions */
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "WeylScal4::Psi4i");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "WeylScal4::Psi4r");
  
}
