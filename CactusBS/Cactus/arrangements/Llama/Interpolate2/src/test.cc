
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
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loopcontrol.h>

#include "mask.h"

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
#define DECLARE_CCTK_ARGUMENTS_CHECKED(func) DECLARE_CCTK_ARGUMENTS
#endif



namespace Interpolate2 {
  
  using namespace std;
  
  
  
  extern "C"
  void
  Interpolate2TestInit (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS_CHECKED(Interpolate2TestInit);
    DECLARE_CCTK_PARAMETERS;
    
    int const reflevel = GetRefinementLevel (cctkGH);
    // int const map = MultiPatch_GetMap (cctkGH);
    
    bool have_interpatch_boundary = false;
#pragma omp parallel reduction(||: have_interpatch_boundary)
    LC_LOOP3 (Interpolate2TestInit,
              i, j, k,
              0, 0, 0,
              cctk_lsh[0], cctk_lsh[1], cctk_lsh[2],
              cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
    {
      int const ind3d = CCTK_GFINDEX3D (cctkGH, i, j, k);
      // Initialise all interior grid points with zero, and all
      // inter-patch boundary points with a non-zero value.
      if (Sn[ind3d] >= 0) {
        have_interpatch_boundary = true;
        test[ind3d] = poison;
        // test[ind3d] = (((1.0 * 100 + map) * 100 + i) * 100 + j) * 100 + k;
      } else {
        // Interior points, outer boundary points, and ghost points
        // are all legal sources for interpolation
        test[ind3d] = 0.0;
      }
    } LC_ENDLOOP3 (Interpolate2TestInit);
    
    if (reflevel > 0 and have_interpatch_boundary) {
      CCTK_VWarn (CCTK_WARN_ABORT,
                  __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Refinement level %d contains inter-patch boundaries; this is currently not supported",
                  reflevel);
    }
  }
  
  
  
  extern "C"
  void
  Interpolate2TestSelectBCs (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS_CHECKED(Interpolate2TestSelectBCs);
    DECLARE_CCTK_PARAMETERS;
    
    CCTK_INT const faces          = CCTK_ALL_FACES;
    CCTK_INT const boundary_width = 0;
    CCTK_INT const table_handle   = -1;
    
    CCTK_INT const ierr =
      Boundary_SelectVarForBC (cctkGH, faces, boundary_width, table_handle,
                               "Interpolate2::Test", "none");
    assert (not ierr);
  }
  
  
  
  extern "C"
  void
  Interpolate2TestCheck (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS_CHECKED(Interpolate2TestCheck);
    DECLARE_CCTK_PARAMETERS;
    
    int const reflevel = GetRefinementLevel (cctkGH);
    int const map = MultiPatch_GetMap(cctkGH);
    
    int error_points = 0;

    LC_LOOP3 (Interpolate2TestCheck,
              i, j, k,
              0, 0, 0,
              cctk_lsh[0], cctk_lsh[1], cctk_lsh[2],
              cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
    {
      int const ind3d = CCTK_GFINDEX3D (cctkGH, i, j, k);
      // if (fabs(test[ind3d]) > 1e-6) {
      if (test[ind3d] != 0.0) {
        ++ error_points;
        CCTK_VInfo (CCTK_THORNSTRING,
                    "Grid point [rl=%d,m=%d:%d,%d,%d] at (%g,%g,%g) is interpolated from an inter-patch boundary point (source_patch=%d, test value=%.17g)",
                    reflevel, map, (int) i, (int) j, (int) k,
                    (double) x[ind3d], (double) y[ind3d], (double) z[ind3d],
                    (int) Sn[ind3d], (double) test[ind3d]);
      }
    } LC_ENDLOOP3 (Interpolate2TestCheck);
    
    if (error_points > 0) {
      int const warnlevel =
        continue_if_selftest_fails ? CCTK_WARN_ALERT : CCTK_WARN_ABORT;
      CCTK_VWarn (warnlevel,
                  __LINE__, __FILE__, CCTK_THORNSTRING,
                  "There were %d grid points which are interpolated from boundary points.  You may need to increase Coordinates::additional_overlap_size.",
                  error_points
                  );
    }
  }
  
} // namespace Interpolate2
