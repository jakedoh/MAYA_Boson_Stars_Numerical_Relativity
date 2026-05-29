#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "loopcontrol.h"

// ShiftTracker_Init:
//   This simply initializes the tracker locations to those
//   specified in the parameter file by the user.
extern "C" void ShiftTracker_Offset(CCTK_ARGUMENTS);

extern "C" void ShiftTracker_Offset(CCTK_ARGUMENTS) {
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS

   //const CCTK_REAL dx CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_SPACE(0);   
   //if (dx < 0.1) return;
  
   //static int done = 0;
   //if (done) return;
   //done = 1;
 
   int tracker_index = 0;
   double vx = st_vx[tracker_index];
   double vy = st_vy[tracker_index];
   double vz = st_vz[tracker_index];
  
   CCTK_VINFO("ShiftTracker_Offset: running!");
 
   for (int i = 0; i < cctk_ash[0] * cctk_ash[1] * cctk_ash[2]; i++) {
     	    beta1[i] -= vx;
      	    beta2[i] -= vy;
	    beta3[i] -= vz;
	    //CCTK_VInfo(CCTK_THORNSTRING,
            // "betax %g", beta1);
	 }

   //st_vx[tracker_index] = 0;
   //st_vy[tracker_index] = 0;
   //st_vz[tracker_index] = 0;

 //   CCTK_LOOP3(StopPunctureOffset, i, j, k, 0, 0, 0,
 //          cctk_lsh[0], cctk_lsh[1], cctk_lsh[2],
 //          cctk_ash[0], cctk_ash[1], cctk_ash[2])
//    {
//    const int ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
//    betax[ind] = 1.0;
//    betay[ind] = 1.0;
//    betaz[ind] = 1.0;
//    }
//    CCTK_ENDLOOP3(StopPunctureOffset);
//    CCTK_SyncGroup(cctkGH, "admbase::shift");  
   
}   
