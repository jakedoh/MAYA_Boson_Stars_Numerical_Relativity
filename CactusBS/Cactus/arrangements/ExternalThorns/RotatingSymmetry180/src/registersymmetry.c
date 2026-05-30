/* $Header$ */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "rotatingsymmetry180.h"



void Rot180_RegisterSymmetry (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  int f;
  CCTK_INT handle;
  CCTK_INT faces[6];
  CCTK_INT width[6];
  CCTK_INT ierr;
  
  for (f=0; f<6; ++f) {
    faces[f] = 0;
    width[f] = 0;
  }
  
  faces[0] = 1;
  width[0] = cctk_nghostzones[0];
  
  handle = SymmetryRegister ("rotating_symmetry_180");
  if (handle < 0) {
    CCTK_WARN (0, "Could not register symmetry boundary condition");
  }
  
  ierr = SymmetryRegisterGrid (cctkGH, handle, faces, width);
  if (ierr < 0) {
    CCTK_WARN (0, "Could not register the symmetry boundaries -- probably some other thorn has already registered the same boundary faces for a different symmetry");
  }
  
  ierr = SymmetryRegisterGridInterpolator
    (cctkGH, handle, Rot180_SymmetryInterpolate);
  if (ierr < 0) {
    CCTK_WARN (0, "Could not register the symmetry interpolator");
  }
}
