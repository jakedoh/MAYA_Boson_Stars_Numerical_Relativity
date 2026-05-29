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

   C Y L I N D R I C A L   D I A G N O S T I C 
   C O M P U T A T I O N 

----------------------------------------------  */

void GRHydroAnalysisComputeCylDiagnostic_local(CCTK_ARGUMENTS);
void GRHydroAnalysisComputeCylDiagnostic_global(CCTK_ARGUMENTS);


void GRHydroAnalysisComputeCylDiagnostic_local(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if ( Verbosity > 0 ) {
       char strbuffer[256];
       sprintf(strbuffer,"Compute local CylDiagnostic at level: dx is %g (CU)",CCTK_DELTA_SPACE(0));
       CCTK_INFO(strbuffer);
    } 

    CCTK_REAL * restrict M0   = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0, 0) ] ); 
    CCTK_REAL * restrict M    = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0, 1) ] ); 
    CCTK_REAL * restrict J    = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0, 2) ] ); 
    CCTK_REAL * restrict T    = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0, 3) ] ); 
    CCTK_REAL * restrict En   = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0, 4) ] ); 
    CCTK_REAL * restrict M1re = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0, 5) ] ); 
    CCTK_REAL * restrict M1im = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0, 6) ] ); 
    CCTK_REAL * restrict M2re = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0, 7) ] ); 
    CCTK_REAL * restrict M2im = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0, 8) ] ); 
    CCTK_REAL * restrict M3re = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0, 9) ] ); 
    CCTK_REAL * restrict M3im = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0,10) ] );
    CCTK_REAL * restrict M4re = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0,11) ] );
    CCTK_REAL * restrict M4im = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0,12) ] );
    CCTK_REAL * restrict MA   = &(TempPoolGF1[ CCTK_VECTGFINDEX3D(cctkGH, 0, 0, 0,13) ] ); 

    CCTK_REAL * restrict velx = &(vel[CCTK_VECTGFINDEX3D(cctkGH,0,0,0,0)]);
    CCTK_REAL * restrict vely = &(vel[CCTK_VECTGFINDEX3D(cctkGH,0,0,0,1)]);
    CCTK_REAL * restrict velz = &(vel[CCTK_VECTGFINDEX3D(cctkGH,0,0,0,2)]);

    int i,j,k;
#pragma omp for 
    for(i=0; i<cctk_lsh[0]; i++)
      for(j=0; j<cctk_lsh[1]; j++)
        for(k=0; k<cctk_lsh[2]; k++)
          {
            CCTK_INT i3D  = CCTK_GFINDEX3D(cctkGH, i, j, k);
            // CCTK_INT i3Dx = CCTK_VECTGFINDEX3D(cctkGH, i,j,k,0);
            // CCTK_INT i3Dy = CCTK_VECTGFINDEX3D(cctkGH, i,j,k,1);
            // CCTK_INT i3Dz = CCTK_VECTGFINDEX3D(cctkGH, i,j,k,2);

            CCTK_REAL LOCgxx,LOCgxy,LOCgxz,LOCgyy,LOCgyz,LOCgzz;
            CCTK_REAL LOCalp,LOCbetax,LOCbetay,LOCbetaz;
            CCTK_REAL LOCrho,LOCeps,LOCpress,LOCvelx,LOCvely,LOCvelz;
            CCTK_REAL LOCx,LOCy,LOCz,LOCRcyl;

            CCTK_REAL LOCsqrtg, LOCW;
            CCTK_REAL LOCu0,LOCux,LOCuy,LOCuz,LOCuD0,LOCuDx,LOCuDy,LOCuDz;
            CCTK_REAL LOCbetaDx,LOCbetaDy,LOCbetaDz;
            CCTK_REAL LOCrhoStar,LOCe0;
            CCTK_REAL TMPRcyl,TMPsin1,TMPcos1,TMPsin2,TMPcos2,TMPsin3,TMPcos3,TMPsin4,TMPcos4;
            CCTK_REAL TMPT0Dx,TMPT0Dy,TMPT0D0,TMPTmuDmu,TMPT0DphiOverRcyl,TMPRcylOMEGA,TMPT0Dphi;

            LOCx = x[i3D];
            LOCy = y[i3D];
            LOCz = z[i3D];

            LOCalp = alp[i3D];
            LOCbetax = betax[i3D];
            LOCbetay = betay[i3D];
            LOCbetaz = betaz[i3D];
            LOCgxx = gxx[i3D];
            LOCgxy = gxy[i3D];
            LOCgxz = gxz[i3D];
            LOCgyy = gyy[i3D];
            LOCgyz = gyz[i3D];
            LOCgzz = gzz[i3D];

            LOCrho  = rho[i3D];
            LOCeps  = eps[i3D];
            LOCpress= press[i3D];
            LOCvelx = velx[i3D]; // vel[i3Dx];
            LOCvely = vely[i3D]; // vel[i3Dy];
            LOCvelz = velz[i3D]; // vel[i3Dz];
 
            // ---------------------------------------------
            // 
            // ---------------------------------------------

            LOCsqrtg = sqrt(LOCgxx * LOCgyy * LOCgzz
               + 2.0 * LOCgxy * LOCgxz * LOCgyz
               - LOCgxy * LOCgxy * LOCgzz
               - LOCgxz * LOCgxz * LOCgyy 
               - LOCgyz * LOCgyz * LOCgxx);

            LOCW    = 1.0 /sqrt( 1 
                 - LOCgxx * LOCvelx * LOCvelx - LOCgyy * LOCvely * LOCvely - LOCgzz * LOCvelz * LOCvelz  
                 - 2.0 *( LOCgxy * LOCvelx * LOCvely + LOCgxz * LOCvelx * LOCvelz + LOCgyz * LOCvely * LOCvelz )
                           );

            //  -------------------------------------------
            //  The contra-variant components of the 4-velocity
            //  :>> u^\mu
            //  ------------------------------------------
            LOCu0 = LOCW / LOCalp;
            LOCux = LOCW * (LOCvelx - LOCbetax / LOCalp);
            LOCuy = LOCW * (LOCvely - LOCbetay / LOCalp);
            LOCuz = LOCW * (LOCvelz - LOCbetaz / LOCalp);
            //  -------------------------------------------
            //  The covariant components of the 4-velocity
            //  :>> u_\mu
            //  ------------------------------------------
            LOCuD0 = LOCW * ( -LOCalp 
               + LOCgxx * LOCvelx * LOCbetax + LOCgxy * LOCvelx * LOCbetay + LOCgxz * LOCvelx * LOCbetaz 
               + LOCgxy * LOCvely * LOCbetax + LOCgyy * LOCvely * LOCbetay + LOCgyz * LOCvely * LOCbetaz 
               + LOCgxz * LOCvelz * LOCbetax + LOCgyz * LOCvelz * LOCbetay + LOCgzz * LOCvelz * LOCbetaz 
               );
            LOCuDx = LOCW * (LOCgxx * LOCvelx + LOCgxy * LOCvely + LOCgxz * LOCvelz);
            LOCuDy = LOCW * (LOCgxy * LOCvelx + LOCgyy * LOCvely + LOCgyz * LOCvelz);
            LOCuDz = LOCW * (LOCgxz * LOCvelx + LOCgyz * LOCvely + LOCgzz * LOCvelz);
            //  -------------------------------------------
            //  The covariant components of the spatial
            //  shift :>> \beta_i
            //  ------------------------------------------
            LOCbetaDx = LOCgxx * LOCbetax + LOCgxy * LOCbetay + LOCgxz * LOCbetaz;
            LOCbetaDy = LOCgxy * LOCbetax + LOCgyy * LOCbetay + LOCgyz * LOCbetaz;
            LOCbetaDz = LOCgxz * LOCbetax + LOCgyz * LOCbetay + LOCgzz * LOCbetaz;
            //  -------------------------------------
            //  Now we defined the conserved density
            //  and some of the component of mater 
            //  energy momentum tensor
            //  -------------------------------------
            
            LOCrhoStar = LOCsqrtg * LOCW * LOCrho;
            LOCe0      = LOCrho * (1 + LOCeps) + LOCpress ;

            //  Plese not that e0 is \rho*h (h=hentalpy)
            // -----------------------------------------------------
            //  Real computation of the variables to be integrated
            // -----------------------------------------------------
            // Now we are ready to define quantity useful
            // to define the modes and circular components
            // -------------------------------------------
            // // sin# = sin(#*\phi)
            // // cos# = cos(#*\phi)
            // // -------------------------------------------

            TMPRcyl = sqrt(LOCx * LOCx + LOCy * LOCy);
            TMPsin1 = ( (TMPRcyl > GRHydro_Analysis_FuzzyTiny) ? LOCy / TMPRcyl : 0.0 ) ;
            TMPcos1 = ( (TMPRcyl > GRHydro_Analysis_FuzzyTiny) ? LOCx / TMPRcyl : 0.0 ) ;
            TMPsin2 = 2 * TMPcos1 * TMPsin1;                    
            TMPcos2 = TMPcos1 * TMPcos1 - TMPsin1 * TMPsin1;  
            TMPsin3 = TMPcos1 * TMPsin2 + TMPsin1 * TMPcos2;
            TMPcos3 = TMPcos1 * TMPcos2 - TMPsin1 * TMPsin2;
            TMPsin4 = 2.0 * TMPcos2 * TMPsin2;
            TMPcos4 = TMPcos2 * TMPcos2 - TMPsin2 * TMPsin2;

            TMPT0Dx   = LOCe0 * LOCu0 * LOCuDx;
            TMPT0Dy   = LOCe0 * LOCu0 * LOCuDy;
            TMPT0D0   = LOCe0 * LOCu0 * LOCuD0  + LOCpress;
            TMPTmuDmu = 3.0 * LOCpress - LOCrho * (1.0+LOCeps);

            TMPT0DphiOverRcyl  = (TMPcos1 * TMPT0Dy - TMPsin1 * TMPT0Dx);
            TMPRcylOMEGA       = (TMPcos1 * LOCuy   - TMPsin1 * LOCux) / LOCu0;
            TMPT0Dphi          = TMPRcyl * TMPT0DphiOverRcyl;

            // -----------------------------------------------------------
            // The following quantities are local scalar and, indeed, 
            // they integated  values is effected by the impossed
            // symmetry only through the corrected UnitVolume
            // -----------------------------------------------------------

            M0 [i3D] = LOCrhoStar;
            M  [i3D] = (-2.0 * TMPT0D0+ TMPTmuDmu) * LOCsqrtg * LOCalp;
            if (rho[i3D] > Atmosphere_Cutoff)
              MA [i3D] = (-2.0 * TMPT0D0+ TMPTmuDmu) * LOCsqrtg * LOCalp;
            else
              MA [i3D] = 0.0;
            J  [i3D] = TMPT0Dphi * LOCsqrtg * LOCalp;
            T  [i3D] = 0.5 * TMPRcylOMEGA *  TMPT0DphiOverRcyl * LOCsqrtg * LOCalp;
            En [i3D] = LOCrhoStar * LOCeps ;

            // -----------------------------------------------------------
            // The following quantities are quite involved and the
            // Symmetries correction is better performed at this level
            // -----------------------------------------------------------
            if ( GlobalGRHydroAnalysisDetect.symRot180 == 1 ) {
                M1re [i3D] = 0.0;
                M1im [i3D] = 0.0;
                M2re [i3D] = LOCrhoStar * TMPcos2 ;
                M2im [i3D] = LOCrhoStar * TMPsin2 ;
                M3re [i3D] = 0.0;
                M3im [i3D] = 0.0;
                M4re [i3D] = LOCrhoStar * TMPcos4 ;
                M4im [i3D] = LOCrhoStar * TMPsin4 ;
            } else if ( GlobalGRHydroAnalysisDetect.symRot90  == 1 ) {
                M1re [i3D] = 0.0;
                M1im [i3D] = 0.0; 
                M2re [i3D] = 0.0;
                M2im [i3D] = 0.0;
                M3re [i3D] = 0.0;
                M3im [i3D] = 0.0;
                M4re [i3D] = LOCrhoStar * TMPcos4 ;
                M4im [i3D] = LOCrhoStar * TMPsin4 ;
            } else if ( GlobalGRHydroAnalysisDetect.symXref   == 1 ) {
                M1re [i3D] = 0.0 ;
                M1im [i3D] = LOCrhoStar * TMPsin1 ;
                M2re [i3D] = LOCrhoStar * TMPcos2 ;
                M2im [i3D] = 0.0 ;
                M3re [i3D] = 0.0 ;
                M3im [i3D] = LOCrhoStar * TMPsin3 ;
                M4re [i3D] = LOCrhoStar * TMPcos4 ;
                M4im [i3D] = 0.0 ;
            } else if ( GlobalGRHydroAnalysisDetect.symYref   == 1 ) {
                M1re [i3D] = LOCrhoStar * TMPcos1 ;
                M1im [i3D] = 0.0 ;
                M2re [i3D] = LOCrhoStar * TMPcos2 ;
                M2im [i3D] = 0.0 ;
                M3re [i3D] = LOCrhoStar * TMPcos3 ;
                M3im [i3D] = 0.0 ; 
                M4re [i3D] = LOCrhoStar * TMPcos4 ;
                M4im [i3D] = 0.0 ; 
            } else {
                M1re [i3D] = LOCrhoStar * TMPcos1 ;
                M1im [i3D] = LOCrhoStar * TMPsin1 ;
                M2re [i3D] = LOCrhoStar * TMPcos2 ;
                M2im [i3D] = LOCrhoStar * TMPsin2 ;
                M3re [i3D] = LOCrhoStar * TMPcos3 ;
                M3im [i3D] = LOCrhoStar * TMPsin3 ;
                M4re [i3D] = LOCrhoStar * TMPcos4 ;
                M4im [i3D] = LOCrhoStar * TMPsin4 ;
            }
          }
}

void GRHydroAnalysisComputeCylDiagnostic_global(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    CCTK_INT reduction_handle;

    if ( Verbosity > 0 ) {
       char strbuffer[256];
       sprintf(strbuffer,"Reduce CylDiagnostic_global: Corser Volume factor is %g (CU)",GlobalGRHydroAnalysisDetect.UnitVolume);
       CCTK_INFO(strbuffer);
    }  

    reduction_handle = CCTK_ReductionHandle("sum");
    if (reduction_handle < 0)
      CCTK_WARN(0, "Unable to get reduction handle.");

    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,CCTK_VARIABLE_REAL,hydro_Cyl_M0   , 1,CCTK_VarIndex("TempPool::TempPoolGF1[0]")))  CCTK_WARN(0, "Error while reducing 0");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,CCTK_VARIABLE_REAL,hydro_Cyl_M    , 1,CCTK_VarIndex("TempPool::TempPoolGF1[1]")))  CCTK_WARN(0, "Error while reducing 1");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,CCTK_VARIABLE_REAL,hydro_Cyl_J    , 1,CCTK_VarIndex("TempPool::TempPoolGF1[2]")))  CCTK_WARN(0, "Error while reducing 2");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,CCTK_VARIABLE_REAL,hydro_Cyl_T    , 1,CCTK_VarIndex("TempPool::TempPoolGF1[3]")))  CCTK_WARN(0, "Error while reducing 3");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,CCTK_VARIABLE_REAL,hydro_Cyl_En   , 1,CCTK_VarIndex("TempPool::TempPoolGF1[4]")))  CCTK_WARN(0, "Error while reducing 4");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,CCTK_VARIABLE_REAL,hydro_Cyl_M1_re, 1,CCTK_VarIndex("TempPool::TempPoolGF1[5]")))  CCTK_WARN(0, "Error while reducing 5");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,CCTK_VARIABLE_REAL,hydro_Cyl_M1_im, 1,CCTK_VarIndex("TempPool::TempPoolGF1[6]")))  CCTK_WARN(0, "Error while reducing 6");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,CCTK_VARIABLE_REAL,hydro_Cyl_M2_re, 1,CCTK_VarIndex("TempPool::TempPoolGF1[7]")))  CCTK_WARN(0, "Error while reducing 7");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,CCTK_VARIABLE_REAL,hydro_Cyl_M2_im, 1,CCTK_VarIndex("TempPool::TempPoolGF1[8]")))  CCTK_WARN(0, "Error while reducing 8");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,CCTK_VARIABLE_REAL,hydro_Cyl_M3_re, 1,CCTK_VarIndex("TempPool::TempPoolGF1[9]")))  CCTK_WARN(0, "Error while reducing 9");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,CCTK_VARIABLE_REAL,hydro_Cyl_M3_im, 1,CCTK_VarIndex("TempPool::TempPoolGF1[10]")))  CCTK_WARN(0, "Error while reducing 10");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,CCTK_VARIABLE_REAL,hydro_Cyl_M4_re, 1,CCTK_VarIndex("TempPool::TempPoolGF1[11]")))  CCTK_WARN(0, "Error while reducing 11");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,CCTK_VARIABLE_REAL,hydro_Cyl_M4_im, 1,CCTK_VarIndex("TempPool::TempPoolGF1[12]")))  CCTK_WARN(0, "Error while reducing 12");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,CCTK_VARIABLE_REAL,hydro_Cyl_MA   , 1,CCTK_VarIndex("TempPool::TempPoolGF1[13]")))  CCTK_WARN(0, "Error while reducing 1");

    *hydro_Cyl_M0    *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *hydro_Cyl_M     *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *hydro_Cyl_J     *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *hydro_Cyl_T     *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *hydro_Cyl_En    *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *hydro_Cyl_M1_re *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *hydro_Cyl_M1_im *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *hydro_Cyl_M2_re *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *hydro_Cyl_M2_im *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *hydro_Cyl_M3_re *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *hydro_Cyl_M3_im *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *hydro_Cyl_M4_re *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *hydro_Cyl_M4_im *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *hydro_Cyl_MA    *= GlobalGRHydroAnalysisDetect.UnitVolume;

    if ( GlobalGRHydroAnalysisDetect.symXref   == 1 ) {
      *hydro_Cyl_M1_re *= 0.0;
      *hydro_Cyl_M2_re *= 0.0;
      *hydro_Cyl_M3_re *= 0.0;
      *hydro_Cyl_M4_re *= 0.0;
    }
    if ( GlobalGRHydroAnalysisDetect.symYref   == 1 ) {
      *hydro_Cyl_M1_im *= 0.0;
      *hydro_Cyl_M2_im *= 0.0;
      *hydro_Cyl_M3_im *= 0.0;
      *hydro_Cyl_M4_im *= 0.0;
    }
    if ( GlobalGRHydroAnalysisDetect.symZref   == 1 ) {
    } 
    if ( GlobalGRHydroAnalysisDetect.symRot90  == 1 ) {
      *hydro_Cyl_M1_re *= 0.0;
      *hydro_Cyl_M1_im *= 0.0;
      *hydro_Cyl_M2_re *= 0.0;
      *hydro_Cyl_M2_im *= 0.0;
      *hydro_Cyl_M3_re *= 0.0;
      *hydro_Cyl_M3_im *= 0.0;
    }
    if ( GlobalGRHydroAnalysisDetect.symRot180 == 1 ) {
      *hydro_Cyl_M1_re *= 0.0;
      *hydro_Cyl_M1_im *= 0.0;
      *hydro_Cyl_M3_re *= 0.0;
      *hydro_Cyl_M3_im *= 0.0;
    }

    *hydro_Cyl_W    = *hydro_Cyl_En + *hydro_Cyl_T + *hydro_Cyl_M0 - *hydro_Cyl_M;
    *hydro_Cyl_beta = *hydro_Cyl_T / *hydro_Cyl_W;

}

