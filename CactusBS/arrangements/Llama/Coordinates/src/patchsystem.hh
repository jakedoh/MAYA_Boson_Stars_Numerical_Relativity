
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

#ifndef PATCHSYSTEM_HH
#define PATCHSYSTEM_HH

#include <vector>

#include <cctk.h>

namespace Coordinates {
  
  using namespace std;
  
  
  
  // For ndirs: 0 = x, 1 = y, 2 = z
  int const ndirs  = 3;
  
  // For nfaces: 0 = lower or -, 1 = upper or +
  int const nfaces = 2;
  
  
  
  // This describes the "physical" domain, which is the domain before
  // discretisation.  Depending on the boundary specification, there
  // may be grid points outside of this domain, e.g. in the overlap
  // and boundary zones
  struct patch_description_t {
    CCTK_REAL xmin[ndirs];
    CCTK_REAL xmax[ndirs];
    int ncells[ndirs];
    bool is_outer_boundary[ndirs][nfaces];
    int bndry_type[ndirs][nfaces];
    int sym_dir[ndirs][nfaces];
  };
  
  
  
  extern int npatches;
  
  extern vector <patch_description_t> patch_descriptions;
  
  CCTK_INT
  SymmetryInterpolate (CCTK_POINTER_TO_CONST const cctkGH_,
                       CCTK_INT const N_dims,
                       CCTK_INT const local_interp_handle,
                       CCTK_INT const param_table_handle,
                       CCTK_INT const coord_system_handle,
                       CCTK_INT const N_interp_points,
                       CCTK_INT const interp_coords_type,
                       CCTK_POINTER_TO_CONST const interp_coords[],
                       CCTK_INT const N_input_arrays,
                       CCTK_INT const input_array_indices[],
                       CCTK_INT const N_output_arrays,
                       CCTK_INT const output_array_types[],
                       CCTK_POINTER const output_arrays[],
                       CCTK_INT const faces);
  
} // namespace Coordinates

#endif  // #ifndef PATCHSYSTEM_HH
