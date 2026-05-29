
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

#include "patchsystem.hh"

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
#define DECLARE_CCTK_ARGUMENTS_CHECKED(func) DECLARE_CCTK_ARGUMENTS
#endif


namespace Coordinates {
  
  using namespace std;
  
  extern "C"
  void
  Coordinates_RegisterSymmetry (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS_CHECKED(Coordinates_RegisterSymmetry );
    DECLARE_CCTK_PARAMETERS;
    
    // Register all faces
    CCTK_INT faces[ndirs*nfaces];
    CCTK_INT width[ndirs*nfaces];
    for (int f=0; f<ndirs*nfaces; ++f) {
      faces[f] = 1;
      width[f] = patch_boundary_size;
    }
    
    CCTK_INT const handle = SymmetryRegister ("multipatch");
    if (handle < 0) {
      CCTK_WARN (CCTK_WARN_ABORT,
                 "Could not register symmetry boundary condition");
    }
    
    CCTK_INT const ierr1 =
      SymmetryRegisterGrid (cctkGH, handle, faces, width);
    if (ierr1 < 0) {
      CCTK_WARN (CCTK_WARN_ABORT,
                 "Could not register the symmetry boundaries -- probably some other thorn has already registered the same boundary faces for a different symmetry");
    }
    
    CCTK_INT ierr2 = 0;
    if (CCTK_EQUALS(symmetry, "full")) {
      // This is the full domain:
      // use standard Coordinates symmetry interpolator as defined within this thorn
      ierr2 =
        SymmetryRegisterGridInterpolator
          (cctkGH, handle, Coordinates_SymmetryInterpolate);
    } else {
      // This is only part of the full domain:
      // use external symmetry interpolator that can handle symmetry conditions (i.e. reflection, rotation)
      ierr2 =
        SymmetryRegisterGridInterpolator
          (cctkGH, handle, CoordinatesSymmetry_Interpolate);
    }
    if (ierr2 < 0) {
      CCTK_WARN (CCTK_WARN_ABORT,
                   "Could not register the symmetry interpolator");
    }
  }
  
} // namespace Coordinates
