#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "reflection.h"

void
ReflectionSymmetry_Register (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_ReflectionSymmetry_Register;
  DECLARE_CCTK_PARAMETERS;
  
  int do_reflection[6];
  
  CCTK_INT nboundaryzones[6];
  CCTK_INT is_internal[6];
  CCTK_INT is_staggered[6];
  CCTK_INT shiftout[6];
  
  CCTK_INT handle;
  CCTK_INT faces[6];
  CCTK_INT width[6];
  int f;
  CCTK_INT ierr;
  
  do_reflection[0] = reflection_x;
  do_reflection[1] = reflection_upper_x;
  do_reflection[2] = reflection_y;
  do_reflection[3] = reflection_upper_y;
  do_reflection[4] = reflection_z;
  do_reflection[5] = reflection_upper_z;
  
  /* Get the boundary specification */
  ierr = GetBoundarySpecification
    (6, nboundaryzones, is_internal, is_staggered, shiftout);
  if (ierr < 0)
  {
    CCTK_WARN (0, "Could not get the boundary specification");
  }
  
  for (f=0; f<6; ++f)
  {
    if (do_reflection[f])
    {
      faces[f] = 1;
      width[f] = nboundaryzones[f/2];
    }
    else
    {
      faces[f] = 0;
      width[f] = 0;
    }
  }
  
  handle = SymmetryRegister ("reflection_symmetry");
  if (handle < 0)
  {
    CCTK_WARN (0, "Could not register symmetry boundary condition");
  }
  
  ierr = SymmetryRegisterGrid (cctkGH, handle, faces, width);
  if (ierr < 0)
  {
    CCTK_WARN (0, "Could not register the symmetry boundaries -- probably some other thorn has already registered the same boundary faces for a different symmetry");
  }
  
  ierr = SymmetryRegisterGridInterpolator
    (cctkGH, handle, ReflectionSymmetry_Interpolate);
  if (ierr < 0)
  {
    CCTK_WARN (0, "Could not register the symmetry interpolator");
  }
}
