 /*@@
   @file      volume_integral.c
   @date      Wed Oct 5 10:18:34 2005
   @author    Frank Loeffler and Luca Baiotti
   @desc 
          Computes the ADM mass as a volume integral.
   @enddesc 
 @@*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "SpaceMask.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

void ADMMass_Volume(CCTK_ARGUMENTS);

void ADMMass_Volume_Global(CCTK_ARGUMENTS);


void ADMMass_Volume(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    CCTK_INT i,j,k, ijk, ghost, ti, tj, tk;
    CCTK_REAL u[3][3], dg[3][3][3];
    CCTK_REAL radius;

    CCTK_INT mask_descriptor = -1, state_descriptor_outside = -1;

    int avoid_origin_parameter;

#include "EinsteinBase/ADMMacros/src/macro/UPPERMET_declare.h"

    /* grid-function strides */
    const CCTK_INT di = 1;
    const CCTK_INT dj = cctk_lsh[0];
    const CCTK_INT dk = cctk_lsh[0]*cctk_lsh[1];

    /* deonminators for derivatives */
    const CCTK_REAL OneOverTwoDX = 1.0 / (2.0 * CCTK_DELTA_SPACE(0));
    const CCTK_REAL OneOverTwoDY = 1.0 / (2.0 * CCTK_DELTA_SPACE(1));
    const CCTK_REAL OneOverTwoDZ = 1.0 / (2.0 * CCTK_DELTA_SPACE(2));

    if (ADMMass_Excise_Horizons)
    {
      mask_descriptor = SpaceMask_GetTypeBits("OutsideMask");
      if (mask_descriptor < 0)
          CCTK_WARN(0, "Thorn OutsideMask not activated, but "
                       "ADMMass_Excise_Horizons requires it.");
      state_descriptor_outside = SpaceMask_GetStateBits("OutsideMask",
                                                        "outside");
      if (state_descriptor_outside < 0)
          CCTK_WARN(0, "Error in obtaining OutsideMask state descriptors");
    }

    if (ADMMass_use_surface_distance_as_volume_radius &&
        (ADMMass_volume_radius[*ADMMass_LoopCounter] < 0.0))
        radius = ADMMass_surface_distance[*ADMMass_LoopCounter];
    else
        radius = ADMMass_volume_radius[*ADMMass_LoopCounter];

    if ((radius <= 0.0)&&(!ADMMass_use_all_volume_as_volume_radius))
    {
        CCTK_WARN(2, "radius < 0 / not set, not calculating "
                     "the volume integral to get the ADM mass.");
        return;
    }

    if (CCTK_IsThornActive("PUGH"))
      {
        const void *ghost_ptr = CCTK_ParameterGet("ghost_size","PUGH",NULL);
        assert( ghost_ptr != NULL );
        ghost = *(const int *)ghost_ptr;
      }
    else /* carpet */ 
      {
        const void *ghost_ptr = CCTK_ParameterGet("ghost_size","Carpet",NULL);
        assert( ghost_ptr != NULL );
        ghost = *(const int *)ghost_ptr;
      }

    for(i=0; i<cctk_lsh[0]; i++)
     for(j=0; j<cctk_lsh[1]; j++)
      for(k=0; k<cctk_lsh[2]; k++)   
      {
          ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
          ADMMass_VolumeMass_GF[ijk] = 0.0;
          ADMMass_VolumeMass_pot_x[ijk] = 0.0;
          ADMMass_VolumeMass_pot_y[ijk] = 0.0;
          ADMMass_VolumeMass_pot_z[ijk] = 0.0;
      }

    for(i=1; i<cctk_lsh[0]-1; i++)
     for(j=1; j<cctk_lsh[1]-1; j++)
      for(k=1; k<cctk_lsh[2]-1; k++)   
      {
          ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

#include "EinsteinBase/ADMMacros/src/macro/UPPERMET_guts.h"
          u[0][0] = UPPERMET_UXX;
          u[0][1] = UPPERMET_UXY;
          u[0][2] = UPPERMET_UXZ;
          u[1][0] = UPPERMET_UXY;
          u[1][1] = UPPERMET_UYY;
          u[1][2] = UPPERMET_UYZ;
          u[2][0] = UPPERMET_UXZ;
          u[2][1] = UPPERMET_UYZ;
          u[2][2] = UPPERMET_UZZ;

          dg[0][0][0] = ( gxx[di+ijk] - gxx[-di+ijk] ) * OneOverTwoDX;
          dg[0][0][1] = ( gxx[dj+ijk] - gxx[-dj+ijk] ) * OneOverTwoDY;
          dg[0][0][2] = ( gxx[dk+ijk] - gxx[-dk+ijk] ) * OneOverTwoDZ;

          dg[0][1][0] = ( gxy[di+ijk] - gxy[-di+ijk] ) * OneOverTwoDX;
          dg[0][1][1] = ( gxy[dj+ijk] - gxy[-dj+ijk] ) * OneOverTwoDY;
          dg[0][1][2] = ( gxy[dk+ijk] - gxy[-dk+ijk] ) * OneOverTwoDZ;

          dg[0][2][0] = ( gxz[di+ijk] - gxz[-di+ijk] ) * OneOverTwoDX;
          dg[0][2][1] = ( gxz[dj+ijk] - gxz[-dj+ijk] ) * OneOverTwoDY;
          dg[0][2][2] = ( gxz[dk+ijk] - gxz[-dk+ijk] ) * OneOverTwoDZ;

          dg[1][0][0] = dg[0][1][0];
          dg[1][0][1] = dg[0][1][1];
          dg[1][0][2] = dg[0][1][2];

          dg[1][1][0] = ( gyy[di+ijk] - gyy[-di+ijk] ) * OneOverTwoDX;
          dg[1][1][1] = ( gyy[dj+ijk] - gyy[-dj+ijk] ) * OneOverTwoDY;
          dg[1][1][2] = ( gyy[dk+ijk] - gyy[-dk+ijk] ) * OneOverTwoDZ;

          dg[1][2][0] = ( gyz[di+ijk] - gyz[-di+ijk] ) * OneOverTwoDX;
          dg[1][2][1] = ( gyz[dj+ijk] - gyz[-dj+ijk] ) * OneOverTwoDY;
          dg[1][2][2] = ( gyz[dk+ijk] - gyz[-dk+ijk] ) * OneOverTwoDZ;

          dg[2][0][0] = dg[0][2][0];
          dg[2][0][1] = dg[0][2][1];
          dg[2][0][2] = dg[0][2][2];

          dg[2][1][0] = dg[1][2][0];
          dg[2][1][1] = dg[1][2][1];
          dg[2][1][2] = dg[1][2][2];
          
          dg[2][2][0] = ( gzz[di+ijk] - gzz[-di+ijk] ) * OneOverTwoDX;
          dg[2][2][1] = ( gzz[dj+ijk] - gzz[-dj+ijk] ) * OneOverTwoDY;
          dg[2][2][2] = ( gzz[dk+ijk] - gzz[-dk+ijk] ) * OneOverTwoDZ;

          for (ti = 0; ti < 3; ti++)
           for (tj = 0; tj < 3; tj++)
            for (tk = 0; tk < 3; tk++)
            {
                ADMMass_VolumeMass_pot_x[ijk] +=
                    u[ti][tj] * u[tk][0] *
                  ( dg[ti][tk][tj] - dg[ti][tj][tk] );
                ADMMass_VolumeMass_pot_y[ijk] +=
                    u[ti][tj] * u[tk][1] *
                  ( dg[ti][tk][tj] - dg[ti][tj][tk] );
                ADMMass_VolumeMass_pot_z[ijk] +=
                    u[ti][tj] * u[tk][2] *
                  ( dg[ti][tk][tj] - dg[ti][tj][tk] );
            }
          ADMMass_VolumeMass_pot_x[ijk] *= alp[ijk] * sqrt(DETG_DETG);
          ADMMass_VolumeMass_pot_y[ijk] *= alp[ijk] * sqrt(DETG_DETG);
          ADMMass_VolumeMass_pot_z[ijk] *= alp[ijk] * sqrt(DETG_DETG);
      }

    /* Do not compute in ghost zones */
    for(i=ghost; i<cctk_lsh[0]-ghost; i++)
     for(j=ghost; j<cctk_lsh[1]-ghost; j++)
      for(k=ghost; k<cctk_lsh[2]-ghost; k++)
      {
          ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

          if ((ADMMass_use_all_volume_as_volume_radius)||
              ((!ADMMass_Excise_Horizons ||
               SpaceMask_CheckStateBits(space_mask, ijk,
                                        mask_descriptor,
                                        state_descriptor_outside)) &&
              ((x[ijk]-ADMMass_x_pos[*ADMMass_LoopCounter])*
               (x[ijk]-ADMMass_x_pos[*ADMMass_LoopCounter]) +
               (y[ijk]-ADMMass_y_pos[*ADMMass_LoopCounter])*
               (y[ijk]-ADMMass_y_pos[*ADMMass_LoopCounter]) +
               (z[ijk]-ADMMass_z_pos[*ADMMass_LoopCounter])*
               (z[ijk]-ADMMass_z_pos[*ADMMass_LoopCounter]) <=
               radius * radius)))
          {
              ADMMass_VolumeMass_GF[ijk] =
               ((ADMMass_VolumeMass_pot_x[CCTK_GFINDEX3D(cctkGH,i+1,j,k)]-
                 ADMMass_VolumeMass_pot_x[CCTK_GFINDEX3D(cctkGH,i-1,j,k)])*
                OneOverTwoDX+
                (ADMMass_VolumeMass_pot_y[CCTK_GFINDEX3D(cctkGH,i,j+1,k)]-
                 ADMMass_VolumeMass_pot_y[CCTK_GFINDEX3D(cctkGH,i,j-1,k)])*
                OneOverTwoDY+
                (ADMMass_VolumeMass_pot_z[CCTK_GFINDEX3D(cctkGH,i,j,k+1)]-
                 ADMMass_VolumeMass_pot_z[CCTK_GFINDEX3D(cctkGH,i,j,k-1)])*
                OneOverTwoDZ);
          }
      }


    /* Carpet does the following itself, but pugh does not. */

    if (CCTK_IsThornActive("PUGH"))     
    {
      /* for short hand, right and left modified stencil values */ 
        const int lst  = ghost;
        const int rst0 = cctk_lsh[0] - ghost;
        const int rst1 = cctk_lsh[1] - ghost;
        const int rst2 = cctk_lsh[2] - ghost;
                
        /*find the value of the "avoid_origin"  parameter */
        const void *avoid_origin_parameter_ptr =
            CCTK_ParameterGet("avoid_origin","CartGrid3D",NULL);
        assert( avoid_origin_parameter_ptr != NULL );
        avoid_origin_parameter = *(const CCTK_INT *)avoid_origin_parameter_ptr;
    
        /* if  avoid_origin = yes (default), then do not devide by two. 
           This gives the correct result for the trapezoidal integration rule
           on the symmetry boundaries, while (even if it is wrong not to devide
           by two the points on the faces that are 
           physical outer boundaries) it should not spoil the result. */
        
        if (! avoid_origin_parameter)
        {        
	  if (cctk_lbnd[2] == 0)
	    for(j = lst; j <= rst1; j++)
	      for(i = lst; i <= rst0; i++)
		{
		  ijk = CCTK_GFINDEX3D(cctkGH, i, j, lst);
		  ADMMass_VolumeMass_GF[ijk] *= 0.5;
		}
	  /* Recompute the integrand on the symmetry and physical boundaries:
	     devide by 2 on the boundary faces cctk_lbnd values start from zero,
	     so must add one to get the total number of points.
	     If cctk_ubnd[2]=n, it means that it is the (n+1)-st point
	     (C notation) */
	  
	  if (cctk_ubnd[2]+1 == cctk_gsh[2])
	    for(j = lst; j <= rst1; j++)
	      for(i = lst; i <= rst0; i++)
		{
		  ijk = CCTK_GFINDEX3D(cctkGH, i, j, rst2);
		  ADMMass_VolumeMass_GF[ijk] *= 0.5;
		}       
	  if (cctk_lbnd[1] == 0)
	    for(k = lst; k<= rst2; k++)
	      for(i = lst; i<= rst0; i++)
		{
		  ijk = CCTK_GFINDEX3D(cctkGH, i, lst, k);
		  ADMMass_VolumeMass_GF[ijk] *= 0.5;
		}
	  if (cctk_ubnd[1]+1 == cctk_gsh[1])
	    for(k = lst; k<= rst2; k++)
	      for(i = lst; i<= rst0; i++)
		{
		  ijk = CCTK_GFINDEX3D(cctkGH, i, rst1, k);
		  ADMMass_VolumeMass_GF[ijk] *= 0.5;
		}
	  if (cctk_lbnd[0] == 0)
	    for(k = lst; k <= rst2; k++)
	      for(j = lst; j <= rst1; j++)
		{
		  ijk = CCTK_GFINDEX3D(cctkGH, lst, j, k);
		  ADMMass_VolumeMass_GF[ijk] *= 0.5;
		}
	  if (cctk_ubnd[0]+1 == cctk_gsh[0])
	    for(k = lst; k <= rst2; k++)
	      for(j = lst; j <= rst1; j++)
		{
		  ijk = CCTK_GFINDEX3D(cctkGH, rst0, j, k);
		  ADMMass_VolumeMass_GF[ijk] *= 0.5;
		}
	} /* end if avoid_origin */
    } /* end if PUGH */
    
    *grid_spacing_product = cctk_delta_space[0]*
                            cctk_delta_space[1]*
                            cctk_delta_space[2];
     
#include "EinsteinBase/ADMMacros/src/macro/UPPERMET_undefine.h"
}      

void ADMMass_Volume_Global(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    CCTK_INT reduction_handle;

    reduction_handle = CCTK_ReductionHandle("sum");
    if (reduction_handle < 0)
        CCTK_WARN(0, "Unable to get reduction handle.");

    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                    CCTK_VARIABLE_REAL,
                    &ADMMass_VolumeMass[*ADMMass_LoopCounter], 1,
                    CCTK_VarIndex("ADMMass::ADMMass_VolumeMass_GF")))
        CCTK_WARN(0, "Error while reducing ADMMass_VolumeMass_GF");

    ADMMass_VolumeMass[*ADMMass_LoopCounter] *=
      *grid_spacing_product / (16.0*PI);

    if (ADMMass_use_surface_distance_as_volume_radius)
    {
        if (ADMMass_Debug)
        CCTK_VInfo(CCTK_THORNSTRING,
                   "detector %d: ADM mass  %g, volume (sphere of radius) %g",
                   (int)*ADMMass_LoopCounter,
                   ADMMass_VolumeMass[*ADMMass_LoopCounter],
                   ADMMass_surface_distance[*ADMMass_LoopCounter]);
    }
    else if (ADMMass_use_all_volume_as_volume_radius)
    {
        if (ADMMass_Debug)
        CCTK_VInfo(CCTK_THORNSTRING," detector %d: ADM mass %g, volume: the whole grid",
               (int)*ADMMass_LoopCounter, ADMMass_VolumeMass[*ADMMass_LoopCounter]);
    }
    else
    {
        if (ADMMass_Debug)
        CCTK_VInfo(CCTK_THORNSTRING," detector %d: ADM mass %g, volume (sphere of radius) %g",
               (int)*ADMMass_LoopCounter, ADMMass_VolumeMass[*ADMMass_LoopCounter],
               ADMMass_volume_radius[*ADMMass_LoopCounter]);
    }
}

