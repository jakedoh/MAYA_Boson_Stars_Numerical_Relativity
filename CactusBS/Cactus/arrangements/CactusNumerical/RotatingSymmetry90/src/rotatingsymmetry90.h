#ifndef ROTATINGSYMMETRY90_H
#define ROTATINGSYMMETRY90_H

#include "cctk.h"
#include "cctk_Arguments.h"

void Rot90_RegisterSymmetry (CCTK_ARGUMENTS);

int BndRot90VI (cGH const * restrict const cctkGH,
                int const nvars,
                CCTK_INT const * restrict const vis);

void Rot90_ApplyBC (CCTK_ARGUMENTS);

CCTK_INT
Rot90_SymmetryInterpolate (CCTK_POINTER_TO_CONST const cctkGH_,
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

#endif /* ! defined ROTATINGSYMMETRY90_H */
