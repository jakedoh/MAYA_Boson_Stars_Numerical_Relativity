
/* Copyright 2013 Peter Diener, Nils Dorband, Roland Haas, Ian Hinder,
Christian Ott, Denis Pollney, Thomas Radke, Christian Reisswig, Erik
Schnetter, Barry Wardell and Burkhard Zink

This file is part of Llama.

Llama is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 2 of the License, or (at your
option) any later version.

Llama is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with Llama.  If not, see <http://www.gnu.org/licenses/>. */

/*
 * 2-patch system with a cylinderical grid (patch 0) and an additional Cartesian
 * column (patch 1) filling the central axis.
 *
 * On the cylindrical patch:
 *
 * (a,b,c) == (r,phi,z)
 * X = r*cos(phi)
 * Y = r**sin(phi)
 * Z = z
 *
 * On the column patch: 
 * X = x
 * Y = y
 * Z = z

 * The local coordinate ranges on the cylindrical patch are:
 * r:     [cylinder_inner_radius,cylinder_outer_radius]
 * phi:   [-Pi,Pi]
 * z:     [cylinder_zmin,cylinder_zmax]

 * The local coordinate ranges on the column patch are:
 * x: [-cylinder_inner_radius,cylinder_inner_radius]
 * y: [-cylinder_inner_radius,cylinder_inner_radius]
 * z: [cylinder_zmin,cylinder_zmax]
 *
 */
#include <cassert>
#include <cmath>
#include <iostream>

#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctk_Arguments.h>

#include "coordinates.hh"
#include "patchsystem.hh"

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
#define DECLARE_CCTK_ARGUMENTS_CHECKED(func) DECLARE_CCTK_ARGUMENTS
#endif


namespace Coordinates
{
  
  using namespace std;
  
  // Determine the extent of the intermediate coordinate system.  This
  // is the local coordinate system of the coarsest grid.
  static
  void
  get_intermediate_coordinates (CCTK_INT const map,
                                CCTK_INT const size,
                                CCTK_REAL * const exterior_min,
                                CCTK_REAL * const exterior_max,
                                CCTK_REAL * const h)
  {
    /*CCTK_REAL physical_min[3];
    CCTK_REAL physical_max[3];
    CCTK_INT ncells[3];
    for (int d=0; d<3; ++d) {
      physical_min[d] = patch_descriptions.at(map).xmin[d];
      physical_max[d] = patch_descriptions.at(map).xmax[d];
      ncells[d] = patch_descriptions.at(map).ncells[d];
      h[d] = (physical_max[d] - physical_min[d]) / ncells[d];
    }
    
    // Calculate the exterior boundaries
    CCTK_REAL interior_min[3];
    CCTK_REAL interior_max[3];
    CCTK_INT const ierr = MultiPatch_ConvertFromPhysicalBoundary
      (map, 3,
       physical_min, physical_max,
       interior_min, interior_max,
       exterior_min, exterior_max,
       h);
    assert (not ierr);*/
    
    CCTK_REAL physical_min[3];
    CCTK_REAL physical_max[3];
    CCTK_REAL interior_min[3];
    CCTK_REAL interior_max[3];
    CCTK_INT const ierr = Coordinates_GetDomainSpecification (map,
                                               size,
                                               physical_min,
                                               physical_max,
                                               interior_min,
                                               interior_max,
                                               exterior_min,
                                               exterior_max,
                                               h);

    assert (not ierr);

  }

  static
  void
  get_local_coordinates (cGH const * const cctkGH,
                         CCTK_REAL const * const x0,
                         CCTK_REAL const * const dx,
                         int       const * const ipos,
                         CCTK_REAL       * const a)
  {
    CCTK_REAL a0[3], da[3];
    for (int d=0; d<3; ++d) {
      // grid spacing on this refinement level
      da[d] = dx[d] / cctkGH->cctk_levfac[d];
      // offset of this refinement level
      a0[d] =
        x0[d] + da[d] * cctkGH->cctk_levoff[d] / cctkGH->cctk_levoffdenom[d];
      // coordinate of this grid point
      a[d] = a0[d] + da[d] * (cctkGH->cctk_lbnd[d] + ipos[d]);
    }
  }



  extern "C"
  void
  Coordinates_SetGlobalCoords_CylinderColumn(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS_CHECKED(Coordinates_SetGlobalCoords_CylinderColumn);
    DECLARE_CCTK_PARAMETERS;
    
    CCTK_INT const map = MultiPatch_GetMap (cctkGH);
    assert (map >= 0);
    
    CCTK_INFO("Setting up global coordinates for cylinder+column.");

    // Determine the extent of the intermediate (see comment above)
    // coordinate system
    CCTK_REAL exterior_min[3];
    CCTK_REAL exterior_max[3];
    CCTK_REAL h[3];
    get_intermediate_coordinates
      (map, 3,
       exterior_min, exterior_max,
       h);

    switch (map) {
      
    case CC_CYLINDER: {
#pragma omp parallel for collapse(3)
      for (int k=0; k<cctk_lsh[2]; ++k)
        for (int j=0; j<cctk_lsh[1]; ++j)
          for (int i=0; i<cctk_lsh[0]; ++i)
            {
              const int ipos[3] = { i, j, k };
              const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
              CCTK_REAL a[3];
              get_local_coordinates (cctkGH, exterior_min, h, ipos, a);

              CCTK_REAL const sinphi = sin(a[1]);
              CCTK_REAL const cosphi = cos(a[1]);

              CCTK_REAL const R = a[0];
	      CCTK_REAL const zloc = a[2];

              x[ijk] = R*cosphi;
              y[ijk] = R*sinphi;
              z[ijk] = zloc;
              r[ijk] = sqrt(x[ijk]*x[ijk] + y[ijk]*y[ijk] + z[ijk]*z[ijk]);
//              printf("(i,j,k) = %d,%d,%d\n",i,j,k);
//              printf("(x,y,z,r) = %f,%f,%f,%f\n",x[ijk],y[ijk],z[ijk],r[ijk]);
            }
      break;
    }
      
    case CC_COLUMN: {
#pragma omp parallel for collapse(3)
      for (int k=0; k<cctk_lsh[2]; ++k)
        for (int j=0; j<cctk_lsh[1]; ++j)
          for (int i=0; i<cctk_lsh[0]; ++i)
            {
              const int ipos[3] = { i, j, k };
              const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
              CCTK_REAL a[3];
              get_local_coordinates (cctkGH, exterior_min, h, ipos, a);
            
              x[ijk] = a[0];
              y[ijk] = a[1];
              z[ijk] = a[2];
              r[ijk] = sqrt(x[ijk]*x[ijk] + y[ijk]*y[ijk] + z[ijk]*z[ijk]);
//              printf("(i,j,k) = %d,%d,%d\n",i,j,k);
//              printf("(x,y,z,r) = %f,%f,%f,%f\n",x[ijk],y[ijk],z[ijk],r[ijk]);
            }
      break;
    }
      
    default:
      CCTK_WARN(0, "Patch number must be 0..1");
      
    }
  }

  extern "C"
  void
  Coordinates_SetJacobian_CylinderColumn(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS_CHECKED(Coordinates_SetJacobian_CylinderColumn);
    DECLARE_CCTK_PARAMETERS;

    using namespace std;

    int i, j, k;

    CCTK_INFO("Setting jacobian: cylinder+column");

    *general_coordinates = 1;
    *interpolate_boundary_points = 0; /* not sure about this - safer to say no */
    *jacobian_state = 1;
    *jacobian_derivative_state = 1;

    const CCTK_INT patch = Coordinates_GetMap(cctkGH);

    switch (patch)
      {
      case CC_CYLINDER:
        {
          for (k=0; k<cctk_lsh[2]; ++k)
            for (j=0; j<cctk_lsh[1]; ++j)
              for (i=0; i<cctk_lsh[0]; ++i)
                {
                  const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

		  const CCTK_REAL R2 = x[ijk]*x[ijk] + y[ijk]*y[ijk];
		  const CCTK_REAL R = sqrt(R2);
		  const CCTK_REAL R3 = R * R2;
		  const CCTK_REAL R4 = R2 * R2;

                  J11[ijk] = x[ijk] / R;
                  J12[ijk] = y[ijk] / R;
                  J13[ijk] = 0;
                  J21[ijk] = -y[ijk] / R2;
                  J22[ijk] = x[ijk] / R2;
                  J23[ijk] = 0;
                  J31[ijk] = 0;
                  J32[ijk] = 0;
                  J33[ijk] = 1.0;

                  dJ111[ijk] = 1.0 / R - x[ijk] / R3;
                  dJ112[ijk] = -x[ijk] * y[ijk] / R3;
                  dJ113[ijk] = 0;
                  dJ122[ijk] = 1.0 / R - y[ijk] / R3;
                  dJ123[ijk] = 0;
                  dJ133[ijk] = 0;

                  dJ211[ijk] = 2 * x[ijk] * y[ijk] / R4;
                  dJ212[ijk] = (y[ijk] * y[ijk] - x[ijk] * x[ijk]) / R4;
                  dJ213[ijk] = 0;
                  dJ222[ijk] = -2 * x[ijk] * y[ijk] / R4;
                  dJ223[ijk] = 0;
                  dJ233[ijk] = 0;

                  dJ311[ijk] = 0;
                  dJ312[ijk] = 0;
                  dJ313[ijk] = 0;
                  dJ322[ijk] = 0;
                  dJ323[ijk] = 0;
                  dJ333[ijk] = 0;
                }
        }
        break;

      case CC_COLUMN:
        {
          for (k=0; k<cctk_lsh[2]; ++k)
            for (j=0; j<cctk_lsh[1]; ++j)
              for (i=0; i<cctk_lsh[0]; ++i)
                {
                  const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

                  J11[ijk] = 1.0;
                  J12[ijk] = 0;
                  J13[ijk] = 0;
                  J21[ijk] = 0;
                  J22[ijk] = 1.0;
                  J23[ijk] = 0;
                  J31[ijk] = 0;
                  J32[ijk] = 0;
                  J33[ijk] = 1.0;
		  
                  dJ111[ijk] = 0;
                  dJ112[ijk] = 0;
                  dJ113[ijk] = 0;
                  dJ122[ijk] = 0;
                  dJ123[ijk] = 0;
                  dJ133[ijk] = 0;

                  dJ211[ijk] = 0;
                  dJ212[ijk] = 0;
                  dJ213[ijk] = 0;
                  dJ222[ijk] = 0;
                  dJ223[ijk] = 0;
                  dJ233[ijk] = 0;

                  dJ311[ijk] = 0;
                  dJ312[ijk] = 0;
                  dJ313[ijk] = 0;
                  dJ322[ijk] = 0;
                  dJ323[ijk] = 0;
                  dJ333[ijk] = 0;
                }
        }
        break;

      default:
        CCTK_WARN(0, "Patch number must be 0..1");

      }
    
    return;
  }
  
  
  extern "C"
  CCTK_INT
  global_to_local_CylinderColumn
  (CCTK_REAL const x, CCTK_REAL const y, CCTK_REAL const z,
   CCTK_REAL localcoord[ndirs])
  {
    DECLARE_CCTK_PARAMETERS;
    
    CCTK_REAL R = sqrt(x * x + y * y);

    if (R < cylinder_inner_radius) 
      {
	localcoord[0] = x;
	localcoord[1] = y;
	localcoord[2] = z;
	return CC_COLUMN;
      }
    else
      {
	localcoord[0] = R;
	localcoord[1] = atan2(y,x);
	localcoord[2] = z;
	return CC_CYLINDER;
      }

    CCTK_WARN(0, "CylinderColumn: Provided coordinates outside of patches");
    return -1;
  }

  extern "C"
  void
  dadx_CylinderColumn
  (CCTK_INT const patch,
   CCTK_REAL const xp, CCTK_REAL const yp, CCTK_REAL const zp,
   CCTK_REAL J[ndirs][ndirs])
  {
    DECLARE_CCTK_PARAMETERS;

    switch (patch)
      {
      case CC_CYLINDER:
        {
	  const CCTK_REAL R2 = xp * xp + yp * yp;
	  const CCTK_REAL R = sqrt(R2);

	  J[0][0] = xp / R;
	  J[0][1] = yp / R;
	  J[0][2] = 0;
	  J[1][0] = -yp / R2;
	  J[1][1] = xp / R2;
	  J[1][2] = 0;
	  J[2][0] = 0;
	  J[2][1] = 0;
	  J[2][2] = 1.0;	  
        }
        break;
      
      case CC_COLUMN:
        {
	  J[0][0] = 1.0;
	  J[0][1] = 0;
	  J[0][2] = 0;
	  J[1][0] = 0;
	  J[1][1] = 1.0;
	  J[1][2] = 0;
	  J[2][0] = 0;
	  J[2][1] = 0;
	  J[2][2] = 1.0;	  
        }
        break;
      
      default:
        {
          CCTK_WARN(0,"Patch number must be 0..1");
        }
      }
  } 
  
  extern "C"
  void
  ddadxdx_CylinderColumn
  (CCTK_INT const patch,
   CCTK_REAL const xp, CCTK_REAL const yp, CCTK_REAL const zp,
   CCTK_REAL dJ[ndirs][ndirs][ndirs])
  {
    DECLARE_CCTK_PARAMETERS;

    switch (patch)
      {
      
      case CC_CYLINDER:
        {
	  const CCTK_REAL R2 = xp * xp + yp * yp;
	  const CCTK_REAL R = sqrt(R2);
	  const CCTK_REAL R3 = R * R2;
	  const CCTK_REAL R4 = R2 * R2;

	  dJ[0][0][0] = 1.0 / R - xp / R3;
	  dJ[0][0][1] = -xp * yp / R3;
	  dJ[0][0][2] = 0;
	  dJ[0][1][0] = dJ[0][0][1];
	  dJ[0][1][1] = 1.0 / R - yp / R3;
	  dJ[0][1][2] = 0;
	  dJ[0][2][0] = dJ[0][0][2];
	  dJ[0][2][1] = dJ[0][1][2];
	  dJ[0][2][2] = 0;

	  dJ[1][0][0] = 2 * xp * yp / R4;
	  dJ[1][0][1] = (yp * yp - xp * xp) / R4;
          dJ[1][0][2] = 0;
	  dJ[1][1][0] = dJ[1][0][1];
	  dJ[1][1][1] =-2 * xp * yp / R4;
	  dJ[1][1][2] = 0;
	  dJ[1][2][0] = dJ[1][0][2];
	  dJ[1][2][1] = dJ[1][1][2];	  
	  dJ[1][2][2] = 0;

	  dJ[2][0][0] = 0;
	  dJ[2][0][1] = 0;
	  dJ[2][0][2] = 0;
	  dJ[2][1][0] = dJ[2][0][1];
	  dJ[2][1][1] = 0;
	  dJ[2][1][2] = 0;
	  dJ[2][2][0] = dJ[2][0][2];
	  dJ[2][2][1] = dJ[2][1][2];	  
	  dJ[2][2][2] = 0;	  
        }
        break;

      case CC_COLUMN:
	{
	  dJ[0][0][0] = 0;
	  dJ[0][0][1] = 0;
	  dJ[0][0][2] = 0;
	  dJ[0][1][0] = dJ[0][0][1];
	  dJ[0][1][1] = 0;
	  dJ[0][1][2] = 0;
	  dJ[0][2][0] = dJ[0][0][2];
	  dJ[0][2][1] = dJ[0][1][2];	  
	  dJ[0][2][2] = 0;	  
	  
	  dJ[1][0][0] = 0;
	  dJ[1][0][1] = 0;
	  dJ[1][0][2] = 0;
	  dJ[1][1][0] = dJ[1][0][1];
	  dJ[1][1][1] = 0;
	  dJ[1][1][2] = 0;
	  dJ[1][2][0] = dJ[1][0][2];
	  dJ[1][2][1] = dJ[1][1][2];	  
	  dJ[1][2][2] = 0;	  

	  dJ[2][0][0] = 0;
	  dJ[2][0][1] = 0;
	  dJ[2][0][2] = 0;
	  dJ[2][1][0] = dJ[2][0][1];
	  dJ[2][1][1] = 0;
	  dJ[2][1][2] = 0;
	  dJ[2][2][0] = dJ[2][0][2];
	  dJ[2][2][1] = dJ[2][1][2];	  
	  dJ[2][2][2] = 0;	  
	}
	break;
      
      default:
        {
          CCTK_WARN(0,"Patch number must be 0..1");
        }
      }
  }

}
