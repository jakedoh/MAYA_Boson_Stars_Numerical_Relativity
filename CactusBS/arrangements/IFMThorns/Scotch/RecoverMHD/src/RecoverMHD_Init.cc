/* Startup and Initialization */
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

extern "C" void RecoverMHD_Startup(CCTK_ARGUMENTS);
extern "C" void RecoverMHD_PctInfo(CCTK_ARGUMENTS);

extern "C" void RecoverMHD_Startup(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS;
   DECLARE_CCTK_PARAMETERS;

   CCTK_RegisterBanner("RecoverMHD: Recovery a la Duez et al");

}


/* Set the centers for Faber Limit */
extern "C" void RecoverMHD_PctInfo(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS;
   DECLARE_CCTK_PARAMETERS;

   int s_ind, sr_ind;

   for (int pct=0; pct<number_of_punctures; pct++) {

       radii[pct] = -1.;
       center_x[pct] = 0;
       center_y[pct] = 0;
       center_z[pct] = 0;

       if ( CCTK_EQUALS(puncture_center_from[pct],"parameter") && CCTK_EQUALS(puncture_radius_from[pct],"parameter") )
       {

	  radii[pct] = puncture_radius[pct];
	  center_x[pct] = puncture_center_x[pct];
	  center_y[pct] = puncture_center_y[pct];
	  center_z[pct] = puncture_center_z[pct];

       } else {

	  s_ind = puncture_surface[pct];
          sr_ind = ( puncture_radius_surface[pct] >= 0 ? puncture_radius_surface[pct] : s_ind );

          if ( sf_active[s_ind] != 0 ) {

		if ( sf_valid[s_ind] < 0 ){
			CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
			  "Surface index %d is not a valid surface this timestep. Not applying limit for puncture %d.",s_ind,pct);
		} else {

		    center_x[pct] = sf_centroid_x[s_ind];
		    center_y[pct] = sf_centroid_y[s_ind];
		    center_z[pct] = sf_centroid_z[s_ind];

		    if ( CCTK_EQUALS(puncture_radius_from[pct],"spherical surface") ) {
                         if ( sr_ind == s_ind ) {
                            radii[pct] = sf_min_radius[s_ind];
                         } else {
                            radii[pct] = ( sf_valid[sr_ind]>=0 && sf_active[sr_ind]!=0 ) ? sf_min_radius[sr_ind] : puncture_radius[pct];
                         }
		    } else {
			 radii[pct] = puncture_radius[pct];
		    }
                    radii[pct] *= puncture_radius_factor[pct];
		
		}
	  } else {

		CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
			   "Surface index %d not active! Not applying limit for puncture %d.",s_ind,pct);

	  }
       }
   }
}

/* Reset the fail mask everywhere */
extern "C" void RecoverMHD_ResetFailMask(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS;

   CCTK_LOOP3_ALL( init_con2prim_failmask, cctkGH, i,j,k ) {
//   #pragma omp parallel for
//   for( int k=0; k<cctk_lsh[2]; k++ )
//      for( int j=0; j<cctk_lsh[1]; j++ )
//         for( int i=0; i<cctk_lsh[0]; i++ ) {
            CCTK_INT idx=CCTK_GFINDEX3D(cctkGH,i,j,k);
            failed_con2prim_mask[idx] = 0;
   } CCTK_ENDLOOP3_ALL( init_con2prim_failmask ); 

}
