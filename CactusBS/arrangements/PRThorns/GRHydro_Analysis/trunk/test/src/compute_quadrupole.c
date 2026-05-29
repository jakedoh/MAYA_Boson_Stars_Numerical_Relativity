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

   Q U D R U P O L E   C O M P U T A T I O N 

----------------------------------------------  */

void GRHydroAnalysisComputeQuadrupole_local(CCTK_ARGUMENTS);
void GRHydroAnalysisComputeQuadrupole_global(CCTK_ARGUMENTS);

void GRHydroAnalysisComputeQuadrupole_local(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if ( Verbosity > 0 ) {
       char strbuffer[256];
       sprintf(strbuffer,"Compute local Quadrupole at level: dx is %g (CU)",CCTK_DELTA_SPACE(0));
       CCTK_INFO(strbuffer);
    } 

    CCTK_REAL * restrict Qxx = & (TempPoolGF1[CCTK_VECTGFINDEX3D(cctkGH,0,0,0, 0)] );
    CCTK_REAL * restrict Qxy = & (TempPoolGF1[CCTK_VECTGFINDEX3D(cctkGH,0,0,0, 1)] );
    CCTK_REAL * restrict Qxz = & (TempPoolGF1[CCTK_VECTGFINDEX3D(cctkGH,0,0,0, 2)] );
    CCTK_REAL * restrict Qyy = & (TempPoolGF1[CCTK_VECTGFINDEX3D(cctkGH,0,0,0, 3)] );
    CCTK_REAL * restrict Qyz = & (TempPoolGF1[CCTK_VECTGFINDEX3D(cctkGH,0,0,0, 4)] );
    CCTK_REAL * restrict Qzz = & (TempPoolGF1[CCTK_VECTGFINDEX3D(cctkGH,0,0,0, 5)] );
    CCTK_REAL * restrict Mx  = & (TempPoolGF1[CCTK_VECTGFINDEX3D(cctkGH,0,0,0, 6)] );
    CCTK_REAL * restrict My  = & (TempPoolGF1[CCTK_VECTGFINDEX3D(cctkGH,0,0,0, 7)] );
    CCTK_REAL * restrict Mz  = & (TempPoolGF1[CCTK_VECTGFINDEX3D(cctkGH,0,0,0, 8)] );
    CCTK_REAL * restrict M0  = & (TempPoolGF1[CCTK_VECTGFINDEX3D(cctkGH,0,0,0, 9)] );

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
            
            CCTK_REAL LOC_Vx,LOC_Vy,LOC_Vz;
            CCTK_REAL LOCdens,LOCx,LOCy,LOCz;

            // ----------------------------------------------------
            // We can advoid using Whisky or GRHydro by explicitely
            // compute dens from the variable define in:
            // * ADMBase and HydroBase
            // ----------------------------------------------------

            CCTK_REAL LOCgxx,LOCgxy,LOCgxz,LOCgyy,LOCgyz,LOCgzz;
            CCTK_REAL LOCsqrtg, LOCW;

            LOCx =  x[i3D];
            LOCy =  y[i3D];
            LOCz =  z[i3D];

            // ----------------------------------------------------
            // LOCdens =  dens[i3D]; (for just using HydroBase)
            // ----------------------------------------------------
            LOCgxx = gxx[i3D];
            LOCgxy = gxy[i3D];
            LOCgxz = gxz[i3D];
            LOCgyy = gyy[i3D];
            LOCgyz = gyz[i3D];
            LOCgzz = gzz[i3D];
            LOCsqrtg = sqrt(LOCgxx * LOCgyy * LOCgzz
               + 2.0 * LOCgxy * LOCgxz * LOCgyz
               - LOCgxy * LOCgxy * LOCgzz
               - LOCgxz * LOCgxz * LOCgyy 
               - LOCgyz * LOCgyz * LOCgxx);
            // LOCW = 1.0 / sqrt(1.0
            //     - LOCgxx * vel[i3Dx] * vel[i3Dx] 
            //     - LOCgyy * vel[i3Dy] * vel[i3Dy] 
            //     - LOCgzz * vel[i3Dz] * vel[i3Dz] 
            //     - 2.0 *( LOCgxy * vel[i3Dx] * vel[i3Dy]
            //            + LOCgxz * vel[i3Dx] * vel[i3Dz]
            //            + LOCgyz * vel[i3Dy] * vel[i3Dz] )
	    //         );
            LOCW = 1.0 / sqrt(1.0
                - LOCgxx * velx[i3D] * velx[i3D] 
                - LOCgyy * vely[i3D] * vely[i3D] 
                - LOCgzz * velz[i3D] * velz[i3D] 
                - 2.0 *( LOCgxy * velx[i3D] * vely[i3D]
                       + LOCgxz * velx[i3D] * velz[i3D]
                       + LOCgyz * vely[i3D] * velz[i3D] )
	            );
            LOCdens = LOCsqrtg * LOCW * rho[i3D]; 
            // ----------------------------------------------------
            // LOCdens =  dens[i3D];  (for just using HydroBase)
            // ----------------------------------------------------

            // LOC_Vx = alp[i3D] * vel[i3Dx] - betax[i3D];
            // LOC_Vy = alp[i3D] * vel[i3Dy] - betay[i3D];
            // LOC_Vz = alp[i3D] * vel[i3Dz] - betaz[i3D];

            // LOC_Vx = alp[i3D] * velx[i3D] - betax[i3D];
            // LOC_Vy = alp[i3D] * vely[i3D] - betay[i3D];
            // LOC_Vz = alp[i3D] * velz[i3D] - betaz[i3D];

            // -----------------------------------------------------------
            // The following quantities are the component of locals 
            // 2-tensor. The correction can be easely erformed after 
            // reduction.
            // -----------------------------------------------------------

	    M0 [i3D] = LOCdens;

	    Mx [i3D] = LOCdens * LOCx;
	    My [i3D] = LOCdens * LOCy;
	    Mz [i3D] = LOCdens * LOCz;

	    Qxx [i3D] = LOCdens * LOCx * LOCx;
	    Qxy [i3D] = LOCdens * LOCx * LOCy;
	    Qxz [i3D] = LOCdens * LOCx * LOCz;
	    Qyy [i3D] = LOCdens * LOCy * LOCy;
	    Qyz [i3D] = LOCdens * LOCy * LOCz;
	    Qzz [i3D] = LOCdens * LOCz * LOCz;

	  }

}

void GRHydroAnalysisComputeQuadrupole_global(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    CCTK_INT reduction_handle;

    reduction_handle = CCTK_ReductionHandle("sum");
    if (reduction_handle < 0)
      CCTK_WARN(0, "Unable to get reduction handle.");

    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,CCTK_VARIABLE_REAL,hydro_Qxx , 1,CCTK_VarIndex("TempPool::TempPoolGF1[0]")))  CCTK_WARN(0, "Error while reducing");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,CCTK_VARIABLE_REAL,hydro_Qxy , 1,CCTK_VarIndex("TempPool::TempPoolGF1[1]")))  CCTK_WARN(0, "Error while reducing");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,CCTK_VARIABLE_REAL,hydro_Qxz , 1,CCTK_VarIndex("TempPool::TempPoolGF1[2]")))  CCTK_WARN(0, "Error while reducing");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,CCTK_VARIABLE_REAL,hydro_Qyy , 1,CCTK_VarIndex("TempPool::TempPoolGF1[3]")))  CCTK_WARN(0, "Error while reducing");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,CCTK_VARIABLE_REAL,hydro_Qyz , 1,CCTK_VarIndex("TempPool::TempPoolGF1[4]")))  CCTK_WARN(0, "Error while reducing");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,CCTK_VARIABLE_REAL,hydro_Qzz , 1,CCTK_VarIndex("TempPool::TempPoolGF1[5]")))  CCTK_WARN(0, "Error while reducing");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,CCTK_VARIABLE_REAL,hydro_Mx  , 1,CCTK_VarIndex("TempPool::TempPoolGF1[6]")))  CCTK_WARN(0, "Error while reducing");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,CCTK_VARIABLE_REAL,hydro_My  , 1,CCTK_VarIndex("TempPool::TempPoolGF1[7]")))  CCTK_WARN(0, "Error while reducing");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,CCTK_VARIABLE_REAL,hydro_Mz  , 1,CCTK_VarIndex("TempPool::TempPoolGF1[8]")))  CCTK_WARN(0, "Error while reducing");
    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,CCTK_VARIABLE_REAL,hydro_M0  , 1,CCTK_VarIndex("TempPool::TempPoolGF1[9]")))  CCTK_WARN(0, "Error while reducing");

    *hydro_Qxx  *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *hydro_Qyy  *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *hydro_Qzz  *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *hydro_Qxy  *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *hydro_Qxz  *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *hydro_Qyz  *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *hydro_Mx   *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *hydro_My   *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *hydro_Mz   *= GlobalGRHydroAnalysisDetect.UnitVolume;
    *hydro_M0   *= GlobalGRHydroAnalysisDetect.UnitVolume;

    if ( GlobalGRHydroAnalysisDetect.symXref   == 1 ) {
      *hydro_Qxy   *= 0.0;
      *hydro_Qxz   *= 0.0;
      *hydro_Mx    *= 0.0;
    }
    if ( GlobalGRHydroAnalysisDetect.symYref   == 1 ) {
      *hydro_Qxy   *= 0.0;
      *hydro_Qyz   *= 0.0;
      *hydro_My    *= 0.0;
    }
    if ( GlobalGRHydroAnalysisDetect.symZref   == 1 ) {
      *hydro_Qxz   *= 0.0;
      *hydro_Qyz   *= 0.0;
      *hydro_Mz    *= 0.0;
    } 
    if ( GlobalGRHydroAnalysisDetect.symRot90  == 1 ) {
      CCTK_REAL tmp;
      tmp  = 0.5*(*hydro_Qxx + *hydro_Qyy);
      *hydro_Qxx   *= tmp;
      *hydro_Qyy   *= tmp;
      *hydro_Qxy   *= 0.0;
      *hydro_Qxz   *= 0.0;
      *hydro_Qyz   *= 0.0;
      *hydro_Mx    *= 0.0;
      *hydro_My    *= 0.0;
    }
    if ( GlobalGRHydroAnalysisDetect.symRot180 == 1 ) {
      *hydro_Qxz   *= 0.0;
      *hydro_Qyz   *= 0.0;
      *hydro_Mx    *= 0.0;
      *hydro_My    *= 0.0;
    }

    *hydro_Xcm = *hydro_Mx / *hydro_M0; 
    *hydro_Ycm = *hydro_My / *hydro_M0; 
    *hydro_Zcm = *hydro_Mz / *hydro_M0; 

    if (GRHydro_Analysis_SphericalSurfaceIndex >= 0) {
      sf_active  [GRHydro_Analysis_SphericalSurfaceIndex] = 1;
      sf_valid   [GRHydro_Analysis_SphericalSurfaceIndex] = 1;
      sf_origin_x[GRHydro_Analysis_SphericalSurfaceIndex] = *hydro_Xcm;
      sf_origin_y[GRHydro_Analysis_SphericalSurfaceIndex] = *hydro_Ycm;
      sf_origin_z[GRHydro_Analysis_SphericalSurfaceIndex] = *hydro_Zcm;
      sf_centroid_x[GRHydro_Analysis_SphericalSurfaceIndex] = *hydro_Xcm;
      sf_centroid_y[GRHydro_Analysis_SphericalSurfaceIndex] = *hydro_Ycm;
      sf_centroid_z[GRHydro_Analysis_SphericalSurfaceIndex] = *hydro_Zcm;
    }

    if ( Verbosity > 0 ) {
       char strbuffer[256];
       sprintf(strbuffer,"Reduce Quadrupole: Corser Volume factor is %g (CU)",GlobalGRHydroAnalysisDetect.UnitVolume);
       CCTK_INFO(strbuffer);
    }  

}



