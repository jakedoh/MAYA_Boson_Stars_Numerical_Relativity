
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

#include <carpet.hh>

#include "patchsystem.hh"


namespace Coordinates {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  // Return the current map, or -1 if in level or global mode
  extern "C"
  CCTK_INT
  Coordinates_GetMap (CCTK_POINTER_TO_CONST const cctkGH_)
  {
    cGH const * const cctkGH = static_cast<cGH const *> (cctkGH_);
    struct CarpetGH const * const carpetGH = GetCarpetGH (cctkGH);
    return carpetGH->map;
  }
  
  
  
  // Return the number of maps
  extern "C"
  CCTK_INT
  Coordinates_GetMaps (CCTK_POINTER_TO_CONST const cctkGH_)
  {
    cGH const * const cctkGH = static_cast<cGH const *> (cctkGH_);
    struct CarpetGH const * const carpetGH = GetCarpetGH (cctkGH);
    return carpetGH->maps;
  }
  
  
  
  // Return which boundaries of the current patch are outer
  // boundaries.  This counts both Cactus outer boundaries and
  // multipatch outer boundaries.
  extern "C"
  CCTK_INT
  Coordinates_GetBbox (CCTK_POINTER_TO_CONST const cctkGH_,
                       CCTK_INT const size,
                       CCTK_INT * const bbox)
  {
    cGH const * restrict const cctkGH = static_cast <cGH const *> (cctkGH_);
    struct CarpetGH const * restrict const carpetGH = GetCarpetGH (cctkGH);
    
    assert (size == ndirs * nfaces);
    assert (bbox);
    
    int const patch = carpetGH->map;
    
    if (patch < 0) {
      // We are not in singlemap mode; indicate an error
      assert (0);
      return -1;
    }
    
    for (int dir = 0; dir < ndirs; ++ dir) {
      for (int face = 0; face < nfaces; ++ face) {
        assert (patch >= 0 and patch < npatches);
        bbox[2*dir+face] =
          patch_descriptions[patch].is_outer_boundary[dir][face];
      }
    }
    
    return 0;
  }
  
  
  
  extern "C"
  CCTK_INT
  Coordinates_GetSystemSpecification (CCTK_INT * const maps)
  {
    assert (npatches > 0);
    * maps = npatches;
    return 0;
  }
  
  
  
  extern "C"
  CCTK_INT
  Coordinates_GetBoundarySpecification (CCTK_INT const map,
                                        CCTK_INT const size,
                                        CCTK_INT * const nboundaryzones,
                                        CCTK_INT * const is_internal,
                                        CCTK_INT * const is_staggered,
                                        CCTK_INT * const shiftout)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (map >= 0 and map < npatches);
    assert (size == ndirs * nfaces);
    assert (nboundaryzones);
    assert (is_internal);
    assert (is_staggered);
    assert (shiftout);
    
    for (int dir = 0; dir < ndirs; ++ dir) {
      for (int face = 0; face < nfaces; ++ face) {
        if (patch_descriptions[map].is_outer_boundary[dir][face]) {
          
          // Outer boundary
          nboundaryzones[2*dir+face] = patch_descriptions[map].bndry_type[dir][face] == 0 ?
                                       outer_boundary_size :
                                       patch_boundary_size;  /// if this is a symmetry boundary, use patch_boundary_size (since it is an internal boundary and not a true outer bndry)
          is_internal   [2*dir+face] = internal_outer_boundaries;
          is_staggered  [2*dir+face] = stagger_outer_boundaries;
          shiftout      [2*dir+face] = (shiftout_outer_boundaries +
                                           patch_descriptions[map].bndry_type[dir][face] != 0 ?
                                           additional_symmetry_size : 0);
          
        } else {
          
          // Inter-patch boundary
          nboundaryzones[2*dir+face] = patch_boundary_size;
          is_internal   [2*dir+face] = 0;
          is_staggered  [2*dir+face] = stagger_patch_boundaries;
          shiftout      [2*dir+face] = (additional_overlap_size +
                                        !stagger_patch_boundaries);
          
        }
      }
    }
    
    return 0;
  }
  
  
  
  extern "C"
  CCTK_INT
  Coordinates_GetDomainSpecification (CCTK_INT const map,
                                      CCTK_INT const size,
                                      CCTK_REAL * const physical_min,
                                      CCTK_REAL * const physical_max,
                                      CCTK_REAL * const interior_min,
                                      CCTK_REAL * const interior_max,
                                      CCTK_REAL * const exterior_min,
                                      CCTK_REAL * const exterior_max,
                                      CCTK_REAL * const spacing)
  {
    assert (map >= 0 and map < npatches);
    assert (size == ndirs);
    assert (physical_min);
    assert (physical_max);
    assert (interior_min);
    assert (interior_max);
    assert (exterior_min);
    assert (exterior_max);
    assert (spacing);
    
    CCTK_INT nboundaryzones[ndirs * nfaces];
    CCTK_INT is_internal[ndirs * nfaces];
    CCTK_INT is_staggered[ndirs * nfaces];
    CCTK_INT shiftout[ndirs * nfaces];
    
    CCTK_INT const ierr = MultiPatch_GetBoundarySpecification
      (map, ndirs * nfaces,
       nboundaryzones, is_internal, is_staggered, shiftout);
    assert (not ierr);
    
    for (int dir = 0; dir < ndirs; ++ dir) {
      physical_min[dir] = patch_descriptions[map].xmin[dir];
      physical_max[dir] = patch_descriptions[map].xmax[dir];
      CCTK_REAL const stagger_lo = is_staggered[2*dir  ] ? 0.5 : 0.0;
      CCTK_REAL const stagger_hi = is_staggered[2*dir+1] ? 0.5 : 0.0;
      CCTK_REAL const ncells =
        patch_descriptions[map].ncells[dir] + stagger_lo + stagger_hi;
      if (fabs(ncells) > 1.0e-12) {
        spacing[dir] = (physical_max[dir] - physical_min[dir]) / ncells;
      } else {
        spacing[dir] = 1.0e-3;  // arbitrary (but better not too large)
      }
    }
    
    return MultiPatch_ConvertFromPhysicalBoundary (map,
                                                   size,
                                                   physical_min, physical_max,
                                                   interior_min, interior_max,
                                                   exterior_min, exterior_max,
                                                   spacing);
  }
  
  
  
  extern "C"
  CCTK_INT
  Coordinates_ConvertFromPhysicalBoundary (CCTK_INT const map,
                                           CCTK_INT const size,
                                           CCTK_REAL const * const physical_min,
                                           CCTK_REAL const * const physical_max,
                                           CCTK_REAL * const interior_min,
                                           CCTK_REAL * const interior_max,
                                           CCTK_REAL * const exterior_min,
                                           CCTK_REAL * const exterior_max,
                                           CCTK_REAL const * const spacing)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (map >= 0 and map < npatches);
    assert (size == ndirs);
    assert (physical_min);
    assert (physical_max);
    assert (interior_min);
    assert (interior_max);
    assert (exterior_min);
    assert (exterior_max);
    assert (spacing);
    
    CCTK_INT nboundaryzones[ndirs * nfaces];
    CCTK_INT is_internal[ndirs * nfaces];
    CCTK_INT is_staggered[ndirs * nfaces];
    CCTK_INT shiftout[ndirs * nfaces];
    
    CCTK_INT const ierr = MultiPatch_GetBoundarySpecification
      (map, ndirs * nfaces,
       nboundaryzones, is_internal, is_staggered, shiftout);
    assert (not ierr);

     for (int d=0; d<ndirs; ++d) {
      exterior_min[d] = physical_min[d] - spacing[d] *
        (+ (is_internal[2*d] ? 0 : nboundaryzones[2*d] - 1)
         + (is_staggered[2*d] ? 0.5 : 0.0)
         + shiftout[2*d]);
      exterior_max[d] = physical_max[d] + spacing[d] *
        (+ (is_internal[2*d+1] ? 0 : nboundaryzones[2*d+1] - 1)
         + (is_staggered[2*d+1] ? 0.5 : 0.0)
         + shiftout[2*d+1]);
      
      interior_min[d] = exterior_min[d] + spacing[d] * nboundaryzones[2*d];
      interior_max[d] = exterior_max[d] - spacing[d] * nboundaryzones[2*d+1];
    }
    
    return 0;
  }
  
  
  
  extern "C"
  CCTK_INT Coordinates_MapIsCartesian (CCTK_INT const map)
  {
    DECLARE_CCTK_PARAMETERS
  
    if (CCTK_Equals(coordinate_system, "cartesian"))
      return 1;
    
    if (CCTK_Equals(coordinate_system, "twopatchcartesian"))
      return 0;

    if (CCTK_Equals(coordinate_system, "TwoPatchDistorted"))
      return 0;
    
    if (CCTK_Equals(coordinate_system, "Thornburg04"))
    {
      if (map == 0)
        return 1;
      else
        return 0;
    }

    if (CCTK_Equals(coordinate_system, "Thornburg13"))
    {
      if (map == 0)
        return 1;
      else
        return 0;
    }

    if (CCTK_Equals(coordinate_system, "Thornburg04nc"))
    {
      return 0;
    }
    
    if (CCTK_Equals(coordinate_system, "CylinderInBox"))
    {
      if (map == 0)
        return 1;
      else
        return 0;
    }
    
    if (CCTK_Equals(coordinate_system, "Sphere+Column"))
    {
      return 0;
    }

    if (CCTK_Equals(coordinate_system, "Cylinder+Column"))
    {
      return 0;
    }
    
    return 1;
  }
  
  
  extern "C"
  CCTK_INT
  Coordinates_ProvidesThornburg04 ()
  {
    DECLARE_CCTK_PARAMETERS
  
    if (CCTK_Equals(coordinate_system, "Thornburg04"))
	return 1;
    
    return 0;
  }
  
   extern "C"
  CCTK_INT
  Coordinates_ProvidesThornburg13 ()
  {
    DECLARE_CCTK_PARAMETERS
  
    if (CCTK_Equals(coordinate_system, "Thornburg13"))
	return 1;
    
    return 0;
  }
  
  
  extern "C"
  CCTK_INT
  Coordinates_GetInnerRadius ()
  {
    DECLARE_CCTK_PARAMETERS
  
    if (CCTK_Equals(coordinate_system, "Thornburg04") || CCTK_Equals(coordinate_system, "Thornburg04nc") || CCTK_Equals(coordinate_system, "Thornburg13"))
	return sphere_inner_radius;
    
    return 0;
  }
  
  
  extern "C"
  CCTK_INT
  Coordinates_GetSymmetrySpecification (CCTK_INT const map,
                                        CCTK_INT const size,
                                        CCTK_INT * const type,
                                        CCTK_INT * const symmetry_overlap,
                                        CCTK_INT * const symmetry_dir)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (map >= 0 and map < npatches);
    assert (size == ndirs * nfaces);
    assert (type);
    assert (symmetry_overlap);
    assert (symmetry_dir);
    
    *symmetry_overlap = additional_symmetry_size;
    
    for (int dir = 0; dir < ndirs; ++ dir) {
      for (int face = 0; face < nfaces; ++ face) {
        
          type[2*dir+face] = patch_descriptions[map].bndry_type[dir][face];
          symmetry_dir[2*dir+face] = patch_descriptions[map].sym_dir[dir][face];
          
      }
    }
    
    return 0;
  }
  
} // namespace Coordinates
