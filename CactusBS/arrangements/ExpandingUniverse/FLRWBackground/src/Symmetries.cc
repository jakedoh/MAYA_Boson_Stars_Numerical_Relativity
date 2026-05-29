
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"
#include "FLRWBackground.h"


extern "C" void FLRWBackground_RegisterSymmetries(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int sym[3];
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;

  SetCartSymVN(cctkGH , sym , "HydroBase::rho");
  SetCartSymVN(cctkGH , sym , "HydroBase::eps");
  SetCartSymVN(cctkGH , sym , "HydroBase::press");
  SetCartSymVN(cctkGH , sym , "HydroBase::vel[0]");
  SetCartSymVN(cctkGH , sym , "HydroBase::vel[1]");
  SetCartSymVN(cctkGH , sym , "HydroBase::vel[2]");
}
