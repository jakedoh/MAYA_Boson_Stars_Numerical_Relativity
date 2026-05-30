
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
  
  // Determine the extent of the global coordinate system.  This
  // is the global coordinate system of the coarsest grid.
  static
  void
  get_global_coordinates (CCTK_INT const map,
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
  Coordinates_SetGlobalCoords_CylinderInBox(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS_CHECKED(Coordinates_SetGlobalCoords_CylinderInBox);
    DECLARE_CCTK_PARAMETERS;
    
    CCTK_INT const map = MultiPatch_GetMap (cctkGH);
    assert (map >= 0);
    
    CCTK_INFO("Setting up global coordinates for a cylinder in a box.");

    // Determine the extent of the global (see comment above)
    // coordinate system
    CCTK_REAL exterior_min[3];
    CCTK_REAL exterior_max[3];
    CCTK_REAL h[3];
    get_global_coordinates
      (map, 3,
       exterior_min, exterior_max,
       h);

    switch (map) {
      
    case CIB_BOX: {
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
              r[ijk] = sqrt (pow(x[ijk],2) + pow(y[ijk],2) + pow(z[ijk],2));
            }
      break;
    }
      
    case CIB_CYLINDER: {
#pragma omp parallel for collapse(3)
      for (int k=0; k<cctk_lsh[2]; ++k)
        for (int j=0; j<cctk_lsh[1]; ++j)
          for (int i=0; i<cctk_lsh[0]; ++i)
            {
              const int ipos[3] = { i, j, k };
              const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
              CCTK_REAL a[3];
              get_local_coordinates (cctkGH, exterior_min, h, ipos, a);
              
              x[ijk] = a[0] * cos(a[1]);
              y[ijk] = a[0] * sin(a[1]);
              z[ijk] = a[2];
              r[ijk] = sqrt (pow(x[ijk],2) + pow(y[ijk],2) + pow(z[ijk],2));
            }
      break;
    }
      
    default:
      CCTK_WARN(0, "Patch number must be 0..1");
      
    }
  }

  extern "C"
  void
  Coordinates_SetJacobian_CylinderInBox(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS_CHECKED(Coordinates_SetJacobian_CylinderInBox);
    DECLARE_CCTK_PARAMETERS;
    
    CCTK_INFO("Setting local coordinates: cylinder in box");
    
    *general_coordinates = 1;
    *interpolate_boundary_points = 0; /* not sure about this - safer to say no */
    *jacobian_state = 1;   
    *jacobian_derivative_state = 1;
 
    CCTK_INT const map = MultiPatch_GetMap (cctkGH);
    assert (map >= 0);
    
    // Determine the extent of the global (see comment above)
    // coordinate system
    CCTK_REAL exterior_min[3];
    CCTK_REAL exterior_max[3];
    CCTK_REAL h[3];
    get_global_coordinates
      (map, 3,
       exterior_min, exterior_max,
       h);
    
    switch (map) {
      
    case CIB_BOX: {
#pragma omp parallel for collapse(3)
      for (int k=0; k<cctk_lsh[2]; ++k)
        for (int j=0; j<cctk_lsh[1]; ++j)
          for (int i=0; i<cctk_lsh[0]; ++i)
            {
              const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
              
              J11[ijk] = 1.0;
              J12[ijk] = 0.0;
              J13[ijk] = 0.0;
              J21[ijk] = 0.0;
              J22[ijk] = 1.0;
              J23[ijk] = 0.0;
              J31[ijk] = 0.0;
              J32[ijk] = 0.0;
              J33[ijk] = 1.0;
              
              
              dJ111[ijk] = 0.0;
              dJ112[ijk] = 0.0;
              dJ113[ijk] = 0.0;
              dJ122[ijk] = 0.0;
              dJ123[ijk] = 0.0;
              dJ133[ijk] = 0.0;
              
              dJ211[ijk] = 0.0;
              dJ212[ijk] = 0.0;
              dJ213[ijk] = 0.0;
              dJ222[ijk] = 0.0;
              dJ223[ijk] = 0.0;
              dJ233[ijk] = 0.0;
              
              dJ311[ijk] = 0.0;
              dJ312[ijk] = 0.0;
              dJ313[ijk] = 0.0;
              dJ322[ijk] = 0.0;
              dJ323[ijk] = 0.0;
              dJ333[ijk] = 0.0;
            }
      break;
    }
      
    case CIB_CYLINDER: {
#pragma omp parallel for collapse(3)
      for (int k=0; k<cctk_lsh[2]; ++k)
        for (int j=0; j<cctk_lsh[1]; ++j)
          for (int i=0; i<cctk_lsh[0]; ++i)
          {
            const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
            
            const CCTK_REAL xp = x[ijk];
            const CCTK_REAL yp = y[ijk];
            
            const CCTK_REAL rhop = sqrt (pow(xp,2) + pow(yp,2));
            
            // x = a cos b
            // y = a sin b
            // z = c
            
            // a = sqrt (x^2 + y^2)
            // b = atan2 (y, x)
            // c = z
            
            // J^a_i = dx^a/dx^i
            // Jax = x / a
            // Jay = y / a
            // Jaz = 0
            // Jbx = +y / a^2
            // Jby = -x / a^2
            // Jbz = 0
            // Jcx = 0
            // Jcy = 0
            // Jcz = 1
            
            J11[ijk] = xp / rhop;
            J12[ijk] = yp / rhop;
            J13[ijk] = 0;
            J21[ijk] = + yp / pow(rhop,2);
            J22[ijk] = - xp / pow(rhop,2);
            J23[ijk] = 0;
            J31[ijk] = 0;
            J32[ijk] = 0;
            J33[ijk] = 1;
            
            // J^a_i,j = ddx^a / dx^i dx^j
            
            dJ111[ijk] = pow(yp,2) / pow(rhop,3);
            dJ112[ijk] = - xp * yp / pow(rhop,3);
            dJ113[ijk] = 0;
            dJ122[ijk] = pow(xp,2) / pow(rhop,3);
            dJ123[ijk] = 0;
            dJ133[ijk] = 0;
            
            dJ211[ijk] = - 2 * xp * yp / pow(rhop,4);
            dJ212[ijk] = (pow(xp,2) - pow(yp,2)) / pow(rhop,4);
            dJ213[ijk] = 0;
            dJ222[ijk] = + 2 * xp * yp / pow(rhop,4);
            dJ223[ijk] = 0;
            dJ233[ijk] = 0;
            
            dJ311[ijk] = 0;
            dJ312[ijk] = 0;
            dJ313[ijk] = 0;
            dJ322[ijk] = 0;
            dJ323[ijk] = 0;
            dJ333[ijk] = 0;
          }
      break;
    }
      
    default:
      CCTK_WARN(0, "Patch number must be 0..1");
      
    }
  }
  
  
  
  extern "C"
  CCTK_INT
  global_to_local_CylinderInBox
  (CCTK_REAL const xp, CCTK_REAL const yp, CCTK_REAL const zp,
   CCTK_REAL a[ndirs])
  {
    DECLARE_CCTK_PARAMETERS;
    
    const CCTK_REAL rhop2 = pow(xp,2) + pow(yp,2);
    
    if (rhop2 > pow (transition_radius, 2))
      {
        a[0] = xp;
        a[1] = yp;
        a[2] = zp;
	
        return CIB_BOX;
      }
    else
      {
        a[0] = sqrt (rhop2);
        a[1] = atan2 (yp, xp);
        a[2] = zp;
        
        return CIB_CYLINDER;
      }
    
    CCTK_WARN(0, "CylinderInBox: Provided coordinates outside of patches");
    return -1;
  }
    
  extern "C"
  void
  dadx_CylinderInBox
  (CCTK_INT const patch,
   CCTK_REAL const xp, CCTK_REAL const yp, CCTK_REAL const zp,
   CCTK_REAL J[ndirs][ndirs])
  {
    DECLARE_CCTK_PARAMETERS;
    
    switch (patch) {
      
    case CIB_BOX: {
      J[0][0] = 1.0;
      J[0][1] = 0.0;
      J[0][2] = 0.0;
      J[1][0] = 0.0;
      J[1][1] = 1.0;
      J[1][2] = 0.0;
      J[2][0] = 0.0;
      J[2][1] = 0.0;
      J[2][2] = 1.0;
      break;
    }
      
    case CIB_CYLINDER: {
      const CCTK_REAL rhop = sqrt (pow(xp,2) + pow(yp,2));
      
      J[0][0] = xp / rhop;
      J[0][1] = yp / rhop;
      J[0][2] = 0;
      J[1][0] = + yp / pow(rhop,2);
      J[1][1] = - xp / pow(rhop,2);
      J[1][2] = 0;
      J[2][0] = 0;
      J[2][1] = 0;
      J[2][2] = 1;
      break;
    }
      
    default:
      CCTK_WARN(0,"Patch number must be 0..1");
    }
  } 
  
  extern "C"
  void
  ddadxdx_CylinderInBox
  (CCTK_INT const patch,
   CCTK_REAL const xp, CCTK_REAL const yp, CCTK_REAL const zp,
   CCTK_REAL dJ[ndirs][ndirs][ndirs])
  {
    DECLARE_CCTK_PARAMETERS;
    
    switch (patch) {
      
    case CIB_BOX: {
      dJ[0][0][0] = 0.0;
      dJ[0][0][1] = 0.0;
      dJ[0][0][2] = 0.0;
      dJ[0][1][1] = 0.0;
      dJ[0][1][2] = 0.0;
      dJ[0][2][2] = 0.0;
      
      dJ[1][0][0] = 0.0;
      dJ[1][0][1] = 0.0;
      dJ[1][0][2] = 0.0;
      dJ[1][1][1] = 0.0;
      dJ[1][1][2] = 0.0;
      dJ[1][2][2] = 0.0;
      
      dJ[2][0][0] = 0.0;
      dJ[2][0][1] = 0.0;
      dJ[2][0][2] = 0.0;
      dJ[2][1][1] = 0.0;
      dJ[2][1][2] = 0.0;
      dJ[2][2][2] = 0.0;
      
      dJ[0][1][0] = dJ[0][0][1];
      dJ[0][2][0] = dJ[0][0][2];
      dJ[0][2][1] = dJ[0][1][2];
      dJ[1][1][0] = dJ[1][0][1];
      dJ[1][2][0] = dJ[1][0][2];
      dJ[1][2][1] = dJ[1][1][2];
      dJ[2][1][0] = dJ[2][0][1];
      dJ[2][2][0] = dJ[2][0][2];
      dJ[2][2][1] = dJ[2][1][2];
      break;
    }
      
    case CIB_CYLINDER: {
      const CCTK_REAL rhop = sqrt (pow(xp,2) + pow(yp,2));
      
      dJ[0][0][0] = pow(yp,2) / pow(rhop,3);
      dJ[0][0][1] = - xp * yp / pow(rhop,3);
      dJ[0][0][2] = 0;
      dJ[0][1][1] = pow(xp,2) / pow(rhop,3);
      dJ[0][1][2] = 0;
      dJ[0][2][2] = 0;
      
      dJ[1][0][0] = - 2 * xp * yp / pow(rhop,4);
      dJ[1][0][1] = (pow(xp,2) - pow(yp,2)) / pow(rhop,4);
      dJ[1][0][2] = 0;
      dJ[1][1][1] = + 2 * xp * yp / pow(rhop,4);
      dJ[1][1][2] = 0;
      dJ[1][2][2] = 0;
      
      dJ[2][0][0] = 0;
      dJ[2][0][1] = 0;
      dJ[2][0][2] = 0;
      dJ[2][1][1] = 0;
      dJ[2][1][2] = 0;
      dJ[2][2][2] = 0;
      
      dJ[0][1][0] = dJ[0][0][1];
      dJ[0][2][0] = dJ[0][0][2];
      dJ[0][2][1] = dJ[0][1][2];
      dJ[1][1][0] = dJ[1][0][1];
      dJ[1][2][0] = dJ[1][0][2];
      dJ[1][2][1] = dJ[1][1][2];
      dJ[2][1][0] = dJ[2][0][1];
      dJ[2][2][0] = dJ[2][0][2];
      dJ[2][2][1] = dJ[2][1][2];
      break;
    }
      
    default:
      CCTK_WARN(0,"Patch number must be 0..1");
    }
  }
}
