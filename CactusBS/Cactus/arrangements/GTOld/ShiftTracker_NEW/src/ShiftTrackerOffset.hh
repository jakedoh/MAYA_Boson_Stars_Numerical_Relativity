#include <iostream>

#include "ShiftTrackerInterpolate.hh"

// ShiftTracker_Offset:
//   This offsets the shift by a constant value so that it is zero at the puncture. 
extern "C" void ShiftTracker_Offset(CCTK_ARGUMENTS);
