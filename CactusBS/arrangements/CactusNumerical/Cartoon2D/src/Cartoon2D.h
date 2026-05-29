/* $Header$ */

#include "cctk.h"

/* Function prototypes */
#ifdef __cplusplus 
extern "C" {
#endif

int BndCartoon2DVI(const cGH *GH, int tensortype, int prolongtype, int vi);

CCTK_INT
Cartoon2D_SymmetryInterpolate (CCTK_POINTER_TO_CONST const cctkGH_,
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

#ifdef __cplusplus 
}
#endif
