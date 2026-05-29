#include <cctk.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include "GRHydro_Analysis_initialize.h"

/* -----------------------------------------------

   S Y M M E T R Y   I N I T I A L I Z A T I O N  

----------------------------------------------  */

struct GlobalGRHydroAnalysis GlobalGRHydroAnalysisDetect;

void GRHydroAnalysis_Init(CCTK_ARGUMENTS);
void GRHydroAnalysis_Init(CCTK_ARGUMENTS){
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    char *xsym,*ysym,*zsym;

    GlobalGRHydroAnalysisDetect.isPUGH    = 0;
    GlobalGRHydroAnalysisDetect.symXref   = 0;
    GlobalGRHydroAnalysisDetect.symYref   = 0;
    GlobalGRHydroAnalysisDetect.symZref   = 0; 
    GlobalGRHydroAnalysisDetect.symRot90  = 0;
    GlobalGRHydroAnalysisDetect.symRot180 = 0;
    GlobalGRHydroAnalysisDetect.UnitVolume = cctk_delta_space[0] * cctk_delta_space[1] * cctk_delta_space[2] ; 

    if ( CCTK_IsThornActive("PUGH") !=0 ) {
       CCTK_INFO("PUGH is active: results may be wrong !");
    }
    if ( CCTK_IsThornActive("RotatingSymmetry90") !=0 ) {
       CCTK_INFO("RotatingSymmetry90 is active.");
       // --------------------------------------------------------------
       // This correspond to the sym-tranformation: (x,y,z) -> (-y,x,z) 
       // need to be applyed 4 times to have the identity
       // --------------------------------------------------------------
       GlobalGRHydroAnalysisDetect.symRot90=1;
       GlobalGRHydroAnalysisDetect.UnitVolume *= 4.0;
    }
    if ( CCTK_IsThornActive("RotatingSymmetry180") !=0 ) {
       CCTK_INFO("RotatingSymmetry180 is active");
       // --------------------------------------------------------------
       // This correspond to the sym-tranformation: (x,y,z) -> (-x,-y.z) 
       // need to be applyed 2 times to have the identity
       // ---------------------------------------------------------------
       GlobalGRHydroAnalysisDetect.symRot180=1;
       GlobalGRHydroAnalysisDetect.UnitVolume *= 2.0;
    }
    if ( CCTK_IsThornActive("ReflectionSymmetry")  !=0 ) {
        xsym=CCTK_ParameterValString("reflection_x","ReflectionSymmetry");
        ysym=CCTK_ParameterValString("reflection_y","ReflectionSymmetry");
        zsym=CCTK_ParameterValString("reflection_z","ReflectionSymmetry");
        if (CCTK_Equals(xsym,"yes") == 1 ) {
            CCTK_INFO("ReflectionSymmetry::reflection_x==yes is active");
            // --------------------------------------------------------------
            // This correspond to the sym-tranformation: (x,y,z) -> (-x,y.z) 
            // need to be applyed 2 times to have the identity
            // ---------------------------------------------------------------
            GlobalGRHydroAnalysisDetect.symXref=1;
            GlobalGRHydroAnalysisDetect.UnitVolume *= 2.0;
        }
        if (CCTK_Equals(ysym,"yes") == 1 ) {
            CCTK_INFO("ReflectionSymmetry::reflection_y==yes is active");
            // --------------------------------------------------------------
            // This correspond to the sym-tranformation: (x,y,z) -> (x,-y.z) 
            // need to be applyed 2 times to have the identity
            // ---------------------------------------------------------------
            GlobalGRHydroAnalysisDetect.symYref=1;
            GlobalGRHydroAnalysisDetect.UnitVolume *= 2.0;
        }
        if (CCTK_Equals(zsym,"yes") == 1 ) {
            CCTK_INFO("ReflectionSymmetry::reflection_z==yes is active");
            // --------------------------------------------------------------
            // This correspond to the sym-tranformation: (x,y,z) -> (x,y.-z) 
            // need to be applyed 2 times to have the identity
            // ---------------------------------------------------------------
            GlobalGRHydroAnalysisDetect.symZref=1;
            GlobalGRHydroAnalysisDetect.UnitVolume *= 2.0;
        }
        free(xsym);
        free(ysym);
        free(zsym);
    }

    /* --   MASK info -------------------  */
    if ( masktype == 0) {
        CCTK_INFO("Computation will not be masked for atmosphere");
        GlobalGRHydroAnalysisDetect.masktype  =0;
        GlobalGRHydroAnalysisDetect.maskvalue =0.0;
    } else {
        CCTK_INFO("Masked not yet implemented");
        GlobalGRHydroAnalysisDetect.masktype  =masktype; 
        GlobalGRHydroAnalysisDetect.maskvalue =maskvalue;
    }

    /* --   -----------------------------   */

    /* --   Output info on the computation  */
    if ( Verbosity > 0 ) { 
      char strbuffer[256];
      sprintf(strbuffer,"Coarser frid dx=%g and UNIT volume correction factor is %g (CU)",cctk_delta_space[0],GlobalGRHydroAnalysisDetect.UnitVolume);
      CCTK_INFO(strbuffer);
    } 


}

