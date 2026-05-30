#include <iostream>

#include "ShiftTracker.hh"

// ShiftTracker_Interp:
//   This evaluates what the shift etc. should be at the 
//   coordinate location, x.  This requires interpolation.
//   The InterpGridArrays API is used for this.
extern "C" void ShiftTracker_Interp(CCTK_ARGUMENTS);
