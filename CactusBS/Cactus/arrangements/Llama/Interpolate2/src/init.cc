
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
#include <vector>
#include <iostream>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "mask.h"

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
#define DECLARE_CCTK_ARGUMENTS_CHECKED(func) DECLARE_CCTK_ARGUMENTS
#endif



namespace Interpolate2 {
  
  using namespace std;
  
  
  
  extern "C"
  void
  Interpolate2Init (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS_CHECKED(Interpolate2Init );
    DECLARE_CCTK_PARAMETERS;
    
    
    
    // Find the interpolation sources for each of my grid points
    
    assert (cctk_dim == 3);
    CCTK_INT const npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
    
    vector<CCTK_POINTER_TO_CONST> globalcoords (cctk_dim);
    globalcoords.at(0) = x;
    globalcoords.at(1) = y;
    globalcoords.at(2) = z;
    
    {
      CCTK_INT const ierr =
        MultiPatch_GlobalToLocal (cctkGH, cctk_dim, npoints,
                                  & globalcoords.front(),
                                  Sn,
                                  NULL,
                                  NULL,
                                  NULL);
      assert (not ierr);
    }
    
    
    
    // Get some information about this patch
    
    CCTK_INT const map = MultiPatch_GetMap (cctkGH);
    assert (map >= 0);
    
    CCTK_INT nboundaryzones[6];
    CCTK_INT is_internal[6];
    CCTK_INT is_staggered[6];
    CCTK_INT shiftout[6];
    {
      CCTK_INT const ierr = MultiPatch_GetBoundarySpecification
        (map, 6,
         nboundaryzones, is_internal, is_staggered, shiftout);
      assert (not ierr);
    }
    
    if (CCTK_MyProc(cctkGH) == 0) {
      for (int d=0; d<3; ++d) {
        for (int f=0; f<2; ++f) {
          if (nboundaryzones[2*d+f] < cctk_nghostzones[d]) {
            CCTK_VWarn (CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                        "The %s %c boundary on patch %d has fewer boundary zones than ghost zones.  This may trip the local interpolator later on.  If so, and if this is not an error, then please report this, this can be corrected.",
                        (f==0 ? "lower" : "upper"), "xyz"[d],
                        (int) map);
          }
        }
      }
    }
    
    CCTK_INT bbox[6];
    {
      CCTK_INT const ierr = MultiPatch_GetBbox (cctkGH, 6, bbox);
      assert (not ierr);
    }
    
    // Ignore all points which are not part of the inter-patch
    // boundary, and test the points which are on the inter-patch
    // boundary
    
    enum boundary_t { bnd_none, bnd_ip, bnd_mp, bnd_outer };
    
    CCTK_REAL num_interior       = 0;
    CCTK_REAL num_ip_boundary    = 0;
    CCTK_REAL num_mp_boundary    = 0;
    CCTK_REAL num_outer_boundary = 0;
    bool thereisaproblem = false;
    
#pragma omp parallel for reduction(+: num_mp_boundary, num_outer_boundary, num_interior, num_ip_boundary) reduction(|: thereisaproblem)
    for (int k=0; k<cctk_lsh[2]; ++k) {
      for (int j=0; j<cctk_lsh[1]; ++j) {
        for (int i=0; i<cctk_lsh[0]; ++i) {
          int const ind3d = CCTK_GFINDEX3D (cctkGH, i, j, k);
          int const ipos[3] = { i, j, k };
          
          // Find out whether this point is on an inter-processor
          // boundary, a multi-patch boundary, or on an outer boundary
          // of one of the faces.
          boundary_t boundary[3][2];
          for (int d=0; d<3; ++d) {
            for (int f=0; f<2; ++f) {
              const bool is_interprocessor_boundary = not cctk_bbox[2*d+f];
              const bool is_outer_boundary = bbox[2*d+f];
              // default value
              boundary[d][f] = bnd_none;
              if (is_interprocessor_boundary) {
                // inter-processor boundary
                if (f==0) {
                  // lower boundary
                  if (ipos[d] < cctk_nghostzones[d]) {
                    boundary[d][f] = bnd_ip;
                  }
                } else {
                  // upper boundary
                  if (ipos[d] >= cctk_lsh[d] - cctk_nghostzones[d]) {
                    boundary[d][f] = bnd_ip;
                  }
                }
              } else if (is_outer_boundary) {
                // outer (physical) boundary
                if (f==0) {
                  // lower boundary
                  if (ipos[d] < nboundaryzones[2*d+f]) {
                    boundary[d][f] = bnd_outer;
                  }
                } else {
                  // upper boundary
                  if (ipos[d] >= cctk_lsh[d] - nboundaryzones[2*d+f]) {
                    boundary[d][f] = bnd_outer;
                  }
                }
              } else {
                // multi-patch boundary
                if (f==0) {
                  // lower boundary
                  if (ipos[d] < nboundaryzones[2*d+f]) {
                    boundary[d][f] = bnd_mp;
                  }
                } else {
                  // upper boundary
                  if (ipos[d] >= cctk_lsh[d] - nboundaryzones[2*d+f]) {
                    boundary[d][f] = bnd_mp;
                  }
                }
              }
            }
          }
          
          // Default value
          boundary_t this_boundary = bnd_none;
          
          // Apply an outer boundary condition to this point if it is
          // on an outer boundary in at least one direction.
          for (int d=0; d<3; ++d) {
            for (int f=0; f<2; ++f) {
              if (boundary[d][f] == bnd_outer) {
                this_boundary = bnd_outer;
              }
            }
          }
          
          // If this point is not on an outer boundary, then
          // interpolate it if it is on a multi-patch boundary.
          // If it is on an outer boundary, interpolate it anyway if
          // it is safe to do so and it is on a multi-patch boundary.
          // 
          // Note that synchronisation happens before interpolation,
          // so that interpolated values cannot be synchronised.
          // (This may be something that can be optimised later on.)
          //
          // Note that this implicitly requires that either the outer
          // boundary condition is applied before the multi-patch
          // boundary condition (unlikely), or that the outer boundary
          // does not require a boundary condition (e.g. because is
          // evolved in time).
          if (this_boundary == bnd_none or *interpolate_boundary_points) {
            for (int d=0; d<3; ++d) {
              for (int f=0; f<2; ++f) {
                if (boundary[d][f] == bnd_mp) {
                  this_boundary = bnd_mp;
                }
              }
            }
          }
          
          // Otherwise, synchronise a point if it is a ghost point in
          // at least one direction.
          if (this_boundary == bnd_none) {
            for (int d=0; d<3; ++d) {
              for (int f=0; f<2; ++f) {
                if (boundary[d][f] == bnd_ip) {
                  this_boundary = bnd_ip;
                }
              }
            }
          }
          
          // If this is an interior point, then maybe we want to
          // interpolate anyway
          CCTK_REAL const xyr2 = pow(x[ind3d],2) + pow(y[ind3d],2);
          if (map == 0 and
              this_boundary == bnd_none and
              ((r[ind3d] >= fill_patch0_radius_min and
                r[ind3d] <= fill_patch0_radius_max) or
               (xyr2 >= pow(fill_patch0_xyradius_min,2) and
                xyr2 <= pow(fill_patch0_xyradius_max,2))))
          {
            this_boundary = bnd_mp;
          }
          
          switch (this_boundary) {
          case bnd_none:
            // Set the source patch to -1 to indicate that this point
            // should not be interpolated because it is in the
            // interior of the domain.
            Sn[ind3d] = SOURCE_PATCH_INTERIOR;
            ++ num_interior;
            break;
          case bnd_outer:
            // Set the source patch to -2 to indicate that this point
            // should not be interpolated, but an outer boundary
            // conditions should be applied.
            Sn[ind3d] = SOURCE_PATCH_OUTER;
            ++ num_outer_boundary;
            break;
          case bnd_mp:
            // This point needs to be interpolated.  Check whether
            // this point can indeed be interpolated.
            thereisaproblem |= Sn[ind3d] < 0;
            ++ num_mp_boundary;
            break;
          case bnd_ip:
            // Set the source patch to -3 to indicate that this point
            // should not be interpolated, but should be synchronised
            // instead.
            Sn[ind3d] = SOURCE_PATCH_GHOST;
            ++ num_ip_boundary;
            break;
          default:
            assert (0);
          }
          
        }
      }
    }
    
    if (thereisaproblem) {
      CCTK_WARN (CCTK_WARN_ABORT,
                 "Cannot determine interpolation source patch for all grid points.  Maybe the patches do not overlap enough?");
    }
    
    // This cannot be done in local mode.  If statistics should be
    // output, they need to be reduced in level or global mode.
#if 0
    {
      CCTK_REAL local[4], global[4];
      local[0] = num_interior;
      local[1] = num_ip_boundary;
      local[2] = num_outer_boundary;
      local[3] = num_mp_boundary;
      int const sum = CCTK_ReductionArrayHandle ("sum");
      assert (sum >= 0);
      int const ierr =
        CCTK_ReduceLocArrayToArray1D
        (cctkGH, 0, sum, local, global, 4, CCTK_VARIABLE_REAL);
      assert (not ierr);
      if (CCTK_MyProc (cctkGH) == 0) {
        num_interior       = global[0];
        num_ip_boundary    = global[1];
        num_outer_boundary = global[2];
        num_mp_boundary    = global[3];
        CCTK_REAL const num_total = cctk_gsh[0] * cctk_gsh[1] * cctk_gsh[2];
        CCTK_VInfo (CCTK_THORNSTRING,
                    "Llama global grid point statistics:");
        CCTK_VInfo (CCTK_THORNSTRING,
                    "   interior points:             %8.1f (%4.1f%%)",
                    double (num_interior),
                    double (num_interior/num_total*100));
        CCTK_VInfo (CCTK_THORNSTRING,
                    "   inter-processor ghosts:      %8.1f (%4.1f%%)",
                    double (num_ip_boundary),
                    double (num_ip_boundary/num_total*100));
        CCTK_VInfo (CCTK_THORNSTRING,
                    "   multi-patch boundary points: %8.1f (%4.1f%%)",
                    double (num_mp_boundary),
                    double (num_mp_boundary/num_total*100));
        CCTK_VInfo (CCTK_THORNSTRING,
                    "   outer boundary points:       %8.1f (%4.1f%%)",
                    double (num_outer_boundary),
                    double (num_outer_boundary/num_total*100));
        CCTK_VInfo (CCTK_THORNSTRING,
                    "   total:                       %8.1f (%4.1f%%)",
                    double (num_total),
                    double (100));
      }
    }
#endif
    
  }
  
  
  
} // namespace Interpolate2
