/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

extern "C" void EvolveScalarFields_RegisterSymmetries(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_EvolveScalarFields_RegisterSymmetries
  DECLARE_CCTK_ARGUMENTS_CHECKED(EvolveScalarFields_RegisterSymmetries);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  /* array holding symmetry definitions */
  int sym[3];
  
  /* Register symmetries of grid functions */
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "EvolveScalarFields::phia");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "EvolveScalarFields::pia");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "EvolveScalarFields::phib");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "EvolveScalarFields::pib");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "EvolveScalarFields::gLocal11");
  
  sym[0] = -1;
  sym[1] = -1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "EvolveScalarFields::gLocal21");
  
  sym[0] = -1;
  sym[1] = 1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "EvolveScalarFields::gLocal31");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "EvolveScalarFields::gLocal22");
  
  sym[0] = 1;
  sym[1] = -1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "EvolveScalarFields::gLocal32");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "EvolveScalarFields::gLocal33");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "EvolveScalarFields::ttrace");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "EvolveScalarFields::SFTtt");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "EvolveScalarFields::rhoE");
  
}
