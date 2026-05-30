
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

#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctk_Arguments.h>

#include "patchsystem.hh"

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
#define DECLARE_CCTK_ARGUMENTS_CHECKED(func) DECLARE_CCTK_ARGUMENTS
#endif

namespace Coordinates
{
  
  using namespace std;
  
  // Determine the extent of the global coordinate system
  static
  void
  get_global_coordinates (CCTK_INT const map,
                          CCTK_INT const size,
                          CCTK_INT * const ncells,
                          CCTK_REAL * const physical_min,
                          CCTK_REAL * const physical_max,
                          CCTK_REAL * const interior_min,
                          CCTK_REAL * const interior_max,
                          CCTK_REAL * const exterior_min,
                          CCTK_REAL * const exterior_max,
                          CCTK_REAL * const spacing)
  {
    /*for (int d=0; d<3; ++d) {
      physical_min[d] = patch_descriptions.at(map).xmin[d];
      physical_max[d] = patch_descriptions.at(map).xmax[d];
      ncells[d] = patch_descriptions.at(map).ncells[d];
      spacing[d] = (physical_max[d] - physical_min[d]) / ncells[d];
    }
    
    // Calculate the exterior boundaries
    CCTK_INT const ierr = MultiPatch_ConvertFromPhysicalBoundary
      (map, 3,
       physical_min, physical_max,
       interior_min, interior_max,
       exterior_min, exterior_max,
       spacing);
    assert (not ierr);*/

    CCTK_INT const ierr = Coordinates_GetDomainSpecification (map,
                                               size,
                                               physical_min,
                                               physical_max,
                                               interior_min,
                                               interior_max,
                                               exterior_min,
                                               exterior_max,
                                               spacing);

    assert (not ierr);

  }
  
  
  
  extern "C"
  void
  Coordinates_SetGlobalCoords_TwoPatchCartesian(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS_CHECKED(Coordinates_SetGlobalCoords_TwoPatchCartesian);
    DECLARE_CCTK_PARAMETERS;
    
    CCTK_INT const map = MultiPatch_GetMap (cctkGH);
    assert (map >= 0);
    
    // Determine the extent of the global coordinate system
    CCTK_INT ncells[3];
    CCTK_REAL physical_min[3];
    CCTK_REAL physical_max[3];
    CCTK_REAL interior_min[3];
    CCTK_REAL interior_max[3];
    CCTK_REAL exterior_min[3];
    CCTK_REAL exterior_max[3];
    CCTK_REAL spacing[3];
    get_global_coordinates
      (map, 3,
       ncells,
       physical_min, physical_max,
       interior_min, interior_max,
       exterior_min, exterior_max,
       spacing);
    
    // Set the global coordinates
    for (int k=0; k<cctk_lsh[2]; ++k)
      for (int j=0; j<cctk_lsh[1]; ++j)
        for (int i=0; i<cctk_lsh[0]; ++i) 
        {
          int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
          
          x[ijk] = exterior_min[0] + (cctk_lbnd[0] + i) * spacing[0];
          y[ijk] = exterior_min[1] + (cctk_lbnd[1] + j) * spacing[1];
          z[ijk] = exterior_min[2] + (cctk_lbnd[2] + k) * spacing[2];
          r[ijk] = sqrt (pow(x[ijk],2) + pow(y[ijk],2) + pow(z[ijk],2));
        }
    
  }
  
  
  
  extern "C"
  void
  Coordinates_SetJacobian_TwoPatchCartesian(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS_CHECKED(Coordinates_SetJacobian_TwoPatchCartesian);
    DECLARE_CCTK_PARAMETERS;
    
    CCTK_INFO("Calculating 'two patch cartesian' jacobian components.");

    *general_coordinates = 1;
    *interpolate_boundary_points = 1;
    *jacobian_state = 1;   
    *jacobian_derivative_state = 1;
 
    // Set Jacobian
    for (int k=0; k<cctk_lsh[2]; ++k)
      for (int j=0; j<cctk_lsh[1]; ++j)
        for (int i=0; i<cctk_lsh[0]; ++i)
        {
          int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
          
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
  }
  
  
  
  extern "C"
  CCTK_INT
  global_to_local_TwoPatchCartesian
  (CCTK_REAL const x, CCTK_REAL const y, CCTK_REAL const z,
   CCTK_REAL localcoord[ndirs])
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (x < patch_one_xmax)
    {
      localcoord[0] = x;
      localcoord[1] = y;
      localcoord[2] = z;
      return 0;
    } 
    else
    {
      localcoord[0] = x;
      localcoord[1] = y;
      localcoord[2] = z;
      return 1;
    }
    
    CCTK_WARN(0, "Twopatch: Provided coordinates outside of patches");
    return -1;
  }
  
  
  
  extern "C"
  void
  dadx_TwoPatchCartesian
  (CCTK_INT const patch,
   CCTK_REAL const a, CCTK_REAL const b, CCTK_REAL const c,
   CCTK_REAL J[ndirs][ndirs])
  {
    DECLARE_CCTK_PARAMETERS;
    
    J[0][0] = 1.0;
    J[0][1] = 0.0;
    J[0][2] = 0.0;
    J[1][0] = 0.0;
    J[1][1] = 1.0;
    J[1][2] = 0.0;
    J[2][0] = 0.0;
    J[2][1] = 0.0;
    J[2][2] = 1.0;
  }
  
  
  
  extern "C"
  void
  ddadxdx_TwoPatchCartesian
  (CCTK_INT const patch,
   CCTK_REAL const a, CCTK_REAL const b, CCTK_REAL const c,
   CCTK_REAL dJ[ndirs][ndirs][ndirs])
  {
    for (int i=0; i < 3; ++i)
      for (int j=0; j < 3; ++j)
        for (int k=0; k < 3; ++k)
          dJ[k][j][i] = 0;
  }
  
} // namespace coordinates 
