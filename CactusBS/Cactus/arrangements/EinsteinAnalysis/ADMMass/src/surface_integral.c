 /*@@
   @file      surface_integral.c
   @date      Wed Oct 5 10:02:12 2005
   @author    Frank Loeffler and Luca Baiotti
   @desc 
          Computes the ADM mass as a surface integral.
   @enddesc 
 @@*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

int find_closest(const cGH *cctkGH, const int *cctk_lsh,
                 const CCTK_REAL *cctk_delta_space, int ghost,
                 CCTK_REAL *coord, CCTK_REAL coord_min, int dir);

void ADMMass_Surface(CCTK_ARGUMENTS);

void ADMMass_Surface_Global(CCTK_ARGUMENTS);

void ADMMass_Surface_Lapse(CCTK_ARGUMENTS);

void ADMMass_Surface_Lapse_Global(CCTK_ARGUMENTS);




int find_closest(const cGH *cctkGH, const int *cctk_lsh,
                 const CCTK_REAL *cctk_delta_space, int ghost,
                 CCTK_REAL *coord, CCTK_REAL coord_min, int dir)
{
    int i, ijk, min_i = -1;
    CCTK_REAL min = 1.e100;
    
    for(i=ghost; i<cctk_lsh[dir]-ghost; i++)
    {
        ijk = CCTK_GFINDEX3D(cctkGH, (dir==0)?i:0, (dir==1)?i:0, (dir==2)?i:0);

        if (fabs(coord[ijk] - coord_min) < min)
        {
            min = fabs(coord[ijk] - coord_min);
            min_i = i;
        }
    }
    return min_i;
}

void ADMMass_Surface(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    CCTK_INT i,j,k, ijk, ierr, ghost, ti, tj, tk, tl;
    CCTK_INT  x_min_i, x_max_i, y_min_j, y_max_j, z_min_k, z_max_k;
    CCTK_REAL ds[3];
    CCTK_REAL u[3][3], dg[3][3][3];

    CCTK_REAL physical_min[3];
    CCTK_REAL physical_max[3];
    CCTK_REAL interior_min[3];
    CCTK_REAL interior_max[3];
    CCTK_REAL exterior_min[3];
    CCTK_REAL exterior_max[3];
    CCTK_REAL spacing[3];

#include "EinsteinBase/ADMMacros/src/macro/UPPERMET_declare.h"

    /* grid-function strides for ADMMacros */
    const CCTK_INT di = 1;
    const CCTK_INT dj = cctk_lsh[0];
    const CCTK_INT dk = cctk_lsh[0]*cctk_lsh[1];

    /* deonminators for derivatives */
    const CCTK_REAL OneOverTwoDX = 1.0 / (2.0 * CCTK_DELTA_SPACE(0));
    const CCTK_REAL OneOverTwoDY = 1.0 / (2.0 * CCTK_DELTA_SPACE(1));
    const CCTK_REAL OneOverTwoDZ = 1.0 / (2.0 * CCTK_DELTA_SPACE(2));

    const CCTK_REAL oneDX =  CCTK_DELTA_SPACE(0);
    const CCTK_REAL oneDY =  CCTK_DELTA_SPACE(1);
    const CCTK_REAL oneDZ =  CCTK_DELTA_SPACE(2);


    x_min_i = -INT_MAX;
    x_max_i =  INT_MAX; 
    y_min_j = -INT_MAX; 
    y_max_j =  INT_MAX;
    z_min_k = -INT_MAX;
    z_max_k =  INT_MAX;

    if (ADMMass_distance_from_grid_boundary[*ADMMass_LoopCounter] > 0.0)
    { 
        if (!CCTK_EQUALS(CCTK_ParameterValString("type", "cartgrid3d"),
                         "coordbase"))
            CCTK_WARN(0,"This thorn used with the ADMMass_distance_from_grid_boundary parameter requires to set coordinates through coordbase.");

        /* Find the physical coordinates of the boundaries */
        ierr = GetDomainSpecification( 3, physical_min, physical_max,
                                          interior_min, interior_max,
                                          exterior_min, exterior_max, spacing);
        if (ierr)
          CCTK_VWarn ( 0, __LINE__, __FILE__, "CartGrid3D", "error returned from function GetDomainSpecification");
        *ADMMass_box_x_min =
            physical_min[0] +
            ADMMass_distance_from_grid_boundary[*ADMMass_LoopCounter];
        *ADMMass_box_x_max =
            physical_max[0] -
            ADMMass_distance_from_grid_boundary[*ADMMass_LoopCounter];
        *ADMMass_box_y_min =
            physical_min[1] +
            ADMMass_distance_from_grid_boundary[*ADMMass_LoopCounter]; 
        *ADMMass_box_y_max =
            physical_max[1] -
            ADMMass_distance_from_grid_boundary[*ADMMass_LoopCounter];
        *ADMMass_box_z_min =
            physical_min[2] +
            ADMMass_distance_from_grid_boundary[*ADMMass_LoopCounter];
        *ADMMass_box_z_max =
            physical_max[2] -
            ADMMass_distance_from_grid_boundary[*ADMMass_LoopCounter];
    }
    else if (ADMMass_surface_distance[*ADMMass_LoopCounter] > 0.0)
    {
        *ADMMass_box_x_min = ADMMass_x_pos[*ADMMass_LoopCounter] -
                             ADMMass_surface_distance[*ADMMass_LoopCounter];
        *ADMMass_box_x_max = ADMMass_x_pos[*ADMMass_LoopCounter] +
                             ADMMass_surface_distance[*ADMMass_LoopCounter];
        *ADMMass_box_y_min = ADMMass_y_pos[*ADMMass_LoopCounter] -
                             ADMMass_surface_distance[*ADMMass_LoopCounter];
        *ADMMass_box_y_max = ADMMass_y_pos[*ADMMass_LoopCounter] +
                             ADMMass_surface_distance[*ADMMass_LoopCounter];
        *ADMMass_box_z_min = ADMMass_z_pos[*ADMMass_LoopCounter] -
                             ADMMass_surface_distance[*ADMMass_LoopCounter];
        *ADMMass_box_z_max = ADMMass_z_pos[*ADMMass_LoopCounter] +
                             ADMMass_surface_distance[*ADMMass_LoopCounter];
    }
    else
    {
        *ADMMass_box_x_min = ADMMass_x_min[*ADMMass_LoopCounter];
        *ADMMass_box_x_max = ADMMass_x_max[*ADMMass_LoopCounter];
        *ADMMass_box_y_min = ADMMass_y_min[*ADMMass_LoopCounter];
        *ADMMass_box_y_max = ADMMass_y_max[*ADMMass_LoopCounter];
        *ADMMass_box_z_min = ADMMass_z_min[*ADMMass_LoopCounter];
        *ADMMass_box_z_max = ADMMass_z_max[*ADMMass_LoopCounter];
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

    ijk = CCTK_GFINDEX3D(cctkGH, cctk_ubnd[0]-cctk_lbnd[0],
                                 cctk_ubnd[1]-cctk_lbnd[1],
                                 cctk_ubnd[2]-cctk_lbnd[2] );

    if ( (x[0] < *ADMMass_box_x_min) && (*ADMMass_box_x_min < x[ijk]) )
        x_min_i = find_closest(cctkGH, cctk_lsh, cctk_delta_space, ghost,
                               x, *ADMMass_box_x_min, 0);
    if ( (x[0] < *ADMMass_box_x_max) && (*ADMMass_box_x_max < x[ijk]) )
        x_max_i = find_closest(cctkGH, cctk_lsh, cctk_delta_space, ghost,
                               x, *ADMMass_box_x_max, 0);
    if ( (y[0] < *ADMMass_box_y_min) && (*ADMMass_box_y_min < y[ijk]) )
        y_min_j = find_closest(cctkGH, cctk_lsh, cctk_delta_space, ghost,
                               y, *ADMMass_box_y_min, 1);
    if ( (y[0] < *ADMMass_box_y_max) && (*ADMMass_box_y_max < y[ijk]) )
        y_max_j = find_closest(cctkGH, cctk_lsh, cctk_delta_space, ghost,
                               y, *ADMMass_box_y_max, 1);
    if ( (z[0] < *ADMMass_box_z_min) && (*ADMMass_box_z_min < z[ijk]) )
        z_min_k = find_closest(cctkGH, cctk_lsh, cctk_delta_space, ghost,
                               z, *ADMMass_box_z_min, 2);
    if ( (z[0] < *ADMMass_box_z_max) && (*ADMMass_box_z_max < z[ijk]) )
        z_max_k = find_closest(cctkGH, cctk_lsh, cctk_delta_space, ghost,
                               z, *ADMMass_box_z_max, 2);

    for(i=0; i<cctk_lsh[0]; i++)
     for(j=0; j<cctk_lsh[1]; j++)
      for(k=0; k<cctk_lsh[2]; k++)
      {
          ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
          ADMMass_SurfaceMass_GF[ijk] = 0.0;
      }
    for(i=ghost; i<cctk_lsh[0]-ghost; i++)
     for(j=ghost; j<cctk_lsh[1]-ghost; j++)
      for(k=ghost; k<cctk_lsh[2]-ghost; k++)
      {
          ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
          
          /* delimit the cube on whose surface we want to integrate */
          if ((i >= x_min_i) &&
              (i <= x_max_i) &&
              (j >= y_min_j) &&
              (j <= y_max_j) &&
              (k >= z_min_k) &&
              (k <= z_max_k))
          {
              ds[0] = 0.0;
              ds[1] = 0.0;
              ds[2] = 0.0;

              /* Select the points on the surfaces of the requested cube */
              CCTK_INT n_bound = 0;
              if (i == x_min_i)
              {
                  n_bound++;
                  ds[0] = -oneDY*oneDZ;
              }
              if (i == x_max_i)
              {
                  n_bound++;
                  ds[0] = oneDY*oneDZ;
              }

              if (j == y_min_j)
              {
                  n_bound++;
                  ds[1] = -oneDX*oneDZ;
              }
              if (j == y_max_j)
              {
                  n_bound++;
                  ds[1] = oneDX*oneDZ;
              }

              if (k == z_min_k)
              {
                  n_bound++;
                  ds[2] = -oneDX*oneDY;
              }
              if (k == z_max_k)
              {
                  n_bound++;
                  ds[2] = oneDX*oneDY;
              }

              /* Take care of corners and edges */
              if (n_bound == 2)
                  for (ti=0; ti<3; ti++)
                      ds[ti] /= 2;
              if (n_bound == 3)
                  for (ti=0; ti<3; ti++)
                      ds[ti] /= 4;

              if ((ds[0] != 0.0) || (ds[1] != 0.0) || (ds[2] != 0.0))
              {
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
                     for (tl = 0; tl < 3; tl++)
                         ADMMass_SurfaceMass_GF[ijk] +=
                           u[ti][tj] * u[tk][tl] *
                             ( dg[ti][tk][tj] - dg[ti][tj][tk] ) * ds[tk];
                  ADMMass_SurfaceMass_GF[ijk] *= sqrt(DETG_DETG);
              }
          }
      }

#include "EinsteinBase/ADMMacros/src/macro/UPPERMET_undefine.h"
}

void ADMMass_Surface_Global(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    CCTK_INT reduction_handle;

    reduction_handle = CCTK_ReductionHandle("sum");
    if (reduction_handle < 0)
        CCTK_WARN(0, "Unable to get reduction handle.");

    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                    CCTK_VARIABLE_REAL,
                    &ADMMass_SurfaceMass[*ADMMass_LoopCounter], 1,
                    CCTK_VarIndex("ADMMass::ADMMass_SurfaceMass_GF")))
        CCTK_WARN(0, "Error while reducing ADMMass_SurfaceMass_GF");

    ADMMass_SurfaceMass[*ADMMass_LoopCounter] /=  16.0*PI;

    if (ADMMass_Debug)
    CCTK_VInfo(CCTK_THORNSTRING,
               "detector %d: ADM mass  %g, surface (cube): xmin %g, xmax %g, ymin %g, ymax %g, zmin %g, zmax %g",
               (int)*ADMMass_LoopCounter,
               ADMMass_SurfaceMass[*ADMMass_LoopCounter],
               *ADMMass_box_x_min,*ADMMass_box_x_max,
               *ADMMass_box_y_min,*ADMMass_box_y_max,
               *ADMMass_box_z_min,*ADMMass_box_z_max);
}

void ADMMass_Surface_Lapse(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    CCTK_INT i,j,k, ijk, ghost;

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

    for(i=ghost; i<cctk_lsh[0]-ghost; i++)
     for(j=ghost; j<cctk_lsh[1]-ghost; j++)
      for(k=ghost; k<cctk_lsh[2]-ghost; k++)
      {
          ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
          ADMMass_SurfaceMass_GF[ijk] *= alp[ijk];
      }
}

void ADMMass_Surface_Lapse_Global(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    CCTK_INT reduction_handle;

    reduction_handle = CCTK_ReductionHandle("sum");
    if (reduction_handle < 0)
        CCTK_WARN(0, "Unable to get reduction handle.");

    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                    CCTK_VARIABLE_REAL,
                    &ADMMass_SurfaceMass_Lapse[*ADMMass_LoopCounter], 1,
                    CCTK_VarIndex("ADMMass::ADMMass_SurfaceMass_GF")))
        CCTK_WARN(0, "Error while reducing ADMMass_SurfaceMass_GF");

    ADMMass_SurfaceMass_Lapse[*ADMMass_LoopCounter] /=  16.0*PI;    

    if (ADMMass_Debug)
    CCTK_VInfo(CCTK_THORNSTRING,
               "detector %d: ADM mass with lapse  %g, surface (cube): xmin %g, xmax %g, ymin %g, ymax %g, zmin %g, zmax %g",
               (int)*ADMMass_LoopCounter,
               ADMMass_SurfaceMass_Lapse[*ADMMass_LoopCounter],
               *ADMMass_box_x_min,*ADMMass_box_x_max,
               *ADMMass_box_y_min,*ADMMass_box_y_max,
               *ADMMass_box_z_min,*ADMMass_box_z_max);
}

