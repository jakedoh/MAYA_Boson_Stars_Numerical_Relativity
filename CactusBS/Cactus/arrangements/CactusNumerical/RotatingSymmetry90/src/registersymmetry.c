#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "rotatingsymmetry90.h"


void Rot90_RegisterSymmetry (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT nboundaryzones[6];
  CCTK_INT is_internal[6];
  CCTK_INT is_staggered[6];
  CCTK_INT shiftout[6];
  
  int f;
  CCTK_INT handle;
  CCTK_INT faces[6];
  CCTK_INT width[6];
  CCTK_INT ierr;
  
  /* Get the boundary specification */
  ierr = GetBoundarySpecification
    (6, nboundaryzones, is_internal, is_staggered, shiftout);
  if (ierr < 0)
  {
    CCTK_WARN (0, "Could not get the boundary specification");
  }
  
  for (f=0; f<6; ++f) {
    faces[f] = 0;
    width[f] = 0;
  }
  
  faces[0] = 1;
  width[0] = nboundaryzones[0];
  faces[2] = 1;
  width[2] = nboundaryzones[2];
  
  handle = SymmetryRegister ("rotating_symmetry_90");
  if (handle < 0) {
    CCTK_WARN (0, "Could not register symmetry boundary condition");
  }
  
  ierr = SymmetryRegisterGrid (cctkGH, handle, faces, width);
  if (ierr < 0) {
    CCTK_WARN (0, "Could not register the symmetry boundaries -- probably some other thorn has already registered the same boundary faces for a different symmetry");
  }
  
  ierr = SymmetryRegisterGridInterpolator
    (cctkGH, handle, Rot90_SymmetryInterpolate);
  if (ierr < 0) {
    CCTK_WARN (0, "Could not register the symmetry interpolator");
  }
}
