
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Faces.h"
#include "util_Table.h"
#include "Symmetry.h"

// Set `none` boundary conditions on BSSN RHSs, as these are set via NewRad.
void Baikal_BoundaryConditions_evolved_gfs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::aDD00GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::aDD00GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::aDD01GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::aDD01GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::aDD02GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::aDD02GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::aDD11GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::aDD11GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::aDD12GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::aDD12GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::aDD22GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::aDD22GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::alphaGF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::alphaGF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::betU0GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::betU0GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::betU1GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::betU1GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::betU2GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::betU2GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::cfGF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::cfGF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::hDD00GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::hDD00GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::hDD01GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::hDD01GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::hDD02GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::hDD02GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::hDD11GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::hDD11GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::hDD12GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::hDD12GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::hDD22GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::hDD22GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::lambdaU0GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::lambdaU0GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::lambdaU1GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::lambdaU1GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::lambdaU2GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::lambdaU2GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::trKGF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::trKGF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::vetU0GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::vetU0GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::vetU1GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::vetU1GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "Baikal::vetU2GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::vetU2GF!");

}

// Set `flat` boundary conditions on BSSN constraints, similar to what Lean does.
void Baikal_BoundaryConditions_aux_gfs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;

  CCTK_INT bndsize;  // <- bndsize = number of ghostzones
  if(FD_order == 2) {
    bndsize = 2;
  } else if(FD_order == 4) {
    bndsize = 3;
  } else {
    CCTK_ERROR("Error: Inside Baikal_BoundaryConditions_aux_gfs(): Chose unsupported FD_order.");
  }
  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, bndsize, -1, "Baikal::HGF", "flat");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::HGF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, bndsize, -1, "Baikal::MU0GF", "flat");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::MU0GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, bndsize, -1, "Baikal::MU1GF", "flat");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::MU1GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, bndsize, -1, "Baikal::MU2GF", "flat");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for Baikal::MU2GF!");
}
