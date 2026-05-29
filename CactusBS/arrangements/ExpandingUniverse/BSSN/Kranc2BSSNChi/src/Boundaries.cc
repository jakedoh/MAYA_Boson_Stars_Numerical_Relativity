/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Faces.h"
#include "util_Table.h"
#include "Symmetry.h"


/* the boundary treatment is split into 3 steps:    */
/* 1. excision                                      */
/* 2. symmetries                                    */
/* 3. "other" boundary conditions, e.g. radiative */

/* to simplify scheduling and testing, the 3 steps  */
/* are currently applied in separate functions      */


extern "C" void Kranc2BSSNChi_CheckBoundaries(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  return;
}

extern "C" void Kranc2BSSNChi_SelectBoundConds(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  
  if (CCTK_EQUALS(h_group_bound, "none"  ) ||
      CCTK_EQUALS(h_group_bound, "static") ||
      CCTK_EQUALS(h_group_bound, "flat"  ) ||
      CCTK_EQUALS(h_group_bound, "zero"  ) )
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::h_group", h_group_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register h_group_bound BC for Kranc2BSSNChi::h_group!");
  }
  
  if (CCTK_EQUALS(A_group_bound, "none"  ) ||
      CCTK_EQUALS(A_group_bound, "static") ||
      CCTK_EQUALS(A_group_bound, "flat"  ) ||
      CCTK_EQUALS(A_group_bound, "zero"  ) )
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::A_group", A_group_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register A_group_bound BC for Kranc2BSSNChi::A_group!");
  }
  
  if (CCTK_EQUALS(chi_group_bound, "none"  ) ||
      CCTK_EQUALS(chi_group_bound, "static") ||
      CCTK_EQUALS(chi_group_bound, "flat"  ) ||
      CCTK_EQUALS(chi_group_bound, "zero"  ) )
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::chi_group", chi_group_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register chi_group_bound BC for Kranc2BSSNChi::chi_group!");
  }
  
  if (CCTK_EQUALS(K_group_bound, "none"  ) ||
      CCTK_EQUALS(K_group_bound, "static") ||
      CCTK_EQUALS(K_group_bound, "flat"  ) ||
      CCTK_EQUALS(K_group_bound, "zero"  ) )
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::K_group", K_group_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register K_group_bound BC for Kranc2BSSNChi::K_group!");
  }
  
  if (CCTK_EQUALS(Gam_group_bound, "none"  ) ||
      CCTK_EQUALS(Gam_group_bound, "static") ||
      CCTK_EQUALS(Gam_group_bound, "flat"  ) ||
      CCTK_EQUALS(Gam_group_bound, "zero"  ) )
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::Gam_group", Gam_group_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register Gam_group_bound BC for Kranc2BSSNChi::Gam_group!");
  }
  
  if (CCTK_EQUALS(alpha_group_bound, "none"  ) ||
      CCTK_EQUALS(alpha_group_bound, "static") ||
      CCTK_EQUALS(alpha_group_bound, "flat"  ) ||
      CCTK_EQUALS(alpha_group_bound, "zero"  ) )
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::alpha_group", alpha_group_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register alpha_group_bound BC for Kranc2BSSNChi::alpha_group!");
  }
  
  if (CCTK_EQUALS(beta_group_bound, "none"  ) ||
      CCTK_EQUALS(beta_group_bound, "static") ||
      CCTK_EQUALS(beta_group_bound, "flat"  ) ||
      CCTK_EQUALS(beta_group_bound, "zero"  ) )
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::beta_group", beta_group_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register beta_group_bound BC for Kranc2BSSNChi::beta_group!");
  }
  
  if (CCTK_EQUALS(betat_group_bound, "none"  ) ||
      CCTK_EQUALS(betat_group_bound, "static") ||
      CCTK_EQUALS(betat_group_bound, "flat"  ) ||
      CCTK_EQUALS(betat_group_bound, "zero"  ) )
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::betat_group", betat_group_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register betat_group_bound BC for Kranc2BSSNChi::betat_group!");
  }
  
  if (CCTK_EQUALS(h11_bound, "none"  ) ||
      CCTK_EQUALS(h11_bound, "static") ||
      CCTK_EQUALS(h11_bound, "flat"  ) ||
      CCTK_EQUALS(h11_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::h11", h11_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register h11_bound BC for Kranc2BSSNChi::h11!");
  }
  
  if (CCTK_EQUALS(h21_bound, "none"  ) ||
      CCTK_EQUALS(h21_bound, "static") ||
      CCTK_EQUALS(h21_bound, "flat"  ) ||
      CCTK_EQUALS(h21_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::h21", h21_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register h21_bound BC for Kranc2BSSNChi::h21!");
  }
  
  if (CCTK_EQUALS(h31_bound, "none"  ) ||
      CCTK_EQUALS(h31_bound, "static") ||
      CCTK_EQUALS(h31_bound, "flat"  ) ||
      CCTK_EQUALS(h31_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::h31", h31_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register h31_bound BC for Kranc2BSSNChi::h31!");
  }
  
  if (CCTK_EQUALS(h22_bound, "none"  ) ||
      CCTK_EQUALS(h22_bound, "static") ||
      CCTK_EQUALS(h22_bound, "flat"  ) ||
      CCTK_EQUALS(h22_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::h22", h22_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register h22_bound BC for Kranc2BSSNChi::h22!");
  }
  
  if (CCTK_EQUALS(h32_bound, "none"  ) ||
      CCTK_EQUALS(h32_bound, "static") ||
      CCTK_EQUALS(h32_bound, "flat"  ) ||
      CCTK_EQUALS(h32_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::h32", h32_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register h32_bound BC for Kranc2BSSNChi::h32!");
  }
  
  if (CCTK_EQUALS(h33_bound, "none"  ) ||
      CCTK_EQUALS(h33_bound, "static") ||
      CCTK_EQUALS(h33_bound, "flat"  ) ||
      CCTK_EQUALS(h33_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::h33", h33_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register h33_bound BC for Kranc2BSSNChi::h33!");
  }
  
  if (CCTK_EQUALS(A11_bound, "none"  ) ||
      CCTK_EQUALS(A11_bound, "static") ||
      CCTK_EQUALS(A11_bound, "flat"  ) ||
      CCTK_EQUALS(A11_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::A11", A11_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register A11_bound BC for Kranc2BSSNChi::A11!");
  }
  
  if (CCTK_EQUALS(A21_bound, "none"  ) ||
      CCTK_EQUALS(A21_bound, "static") ||
      CCTK_EQUALS(A21_bound, "flat"  ) ||
      CCTK_EQUALS(A21_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::A21", A21_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register A21_bound BC for Kranc2BSSNChi::A21!");
  }
  
  if (CCTK_EQUALS(A31_bound, "none"  ) ||
      CCTK_EQUALS(A31_bound, "static") ||
      CCTK_EQUALS(A31_bound, "flat"  ) ||
      CCTK_EQUALS(A31_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::A31", A31_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register A31_bound BC for Kranc2BSSNChi::A31!");
  }
  
  if (CCTK_EQUALS(A22_bound, "none"  ) ||
      CCTK_EQUALS(A22_bound, "static") ||
      CCTK_EQUALS(A22_bound, "flat"  ) ||
      CCTK_EQUALS(A22_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::A22", A22_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register A22_bound BC for Kranc2BSSNChi::A22!");
  }
  
  if (CCTK_EQUALS(A32_bound, "none"  ) ||
      CCTK_EQUALS(A32_bound, "static") ||
      CCTK_EQUALS(A32_bound, "flat"  ) ||
      CCTK_EQUALS(A32_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::A32", A32_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register A32_bound BC for Kranc2BSSNChi::A32!");
  }
  
  if (CCTK_EQUALS(A33_bound, "none"  ) ||
      CCTK_EQUALS(A33_bound, "static") ||
      CCTK_EQUALS(A33_bound, "flat"  ) ||
      CCTK_EQUALS(A33_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::A33", A33_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register A33_bound BC for Kranc2BSSNChi::A33!");
  }
  
  if (CCTK_EQUALS(chi_bound, "none"  ) ||
      CCTK_EQUALS(chi_bound, "static") ||
      CCTK_EQUALS(chi_bound, "flat"  ) ||
      CCTK_EQUALS(chi_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::chi", chi_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register chi_bound BC for Kranc2BSSNChi::chi!");
  }
  
  if (CCTK_EQUALS(K_bound, "none"  ) ||
      CCTK_EQUALS(K_bound, "static") ||
      CCTK_EQUALS(K_bound, "flat"  ) ||
      CCTK_EQUALS(K_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::K", K_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register K_bound BC for Kranc2BSSNChi::K!");
  }
  
  if (CCTK_EQUALS(Gam1_bound, "none"  ) ||
      CCTK_EQUALS(Gam1_bound, "static") ||
      CCTK_EQUALS(Gam1_bound, "flat"  ) ||
      CCTK_EQUALS(Gam1_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::Gam1", Gam1_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register Gam1_bound BC for Kranc2BSSNChi::Gam1!");
  }
  
  if (CCTK_EQUALS(Gam2_bound, "none"  ) ||
      CCTK_EQUALS(Gam2_bound, "static") ||
      CCTK_EQUALS(Gam2_bound, "flat"  ) ||
      CCTK_EQUALS(Gam2_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::Gam2", Gam2_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register Gam2_bound BC for Kranc2BSSNChi::Gam2!");
  }
  
  if (CCTK_EQUALS(Gam3_bound, "none"  ) ||
      CCTK_EQUALS(Gam3_bound, "static") ||
      CCTK_EQUALS(Gam3_bound, "flat"  ) ||
      CCTK_EQUALS(Gam3_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::Gam3", Gam3_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register Gam3_bound BC for Kranc2BSSNChi::Gam3!");
  }
  
  if (CCTK_EQUALS(alpha_bound, "none"  ) ||
      CCTK_EQUALS(alpha_bound, "static") ||
      CCTK_EQUALS(alpha_bound, "flat"  ) ||
      CCTK_EQUALS(alpha_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::alpha", alpha_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register alpha_bound BC for Kranc2BSSNChi::alpha!");
  }
  
  if (CCTK_EQUALS(beta1_bound, "none"  ) ||
      CCTK_EQUALS(beta1_bound, "static") ||
      CCTK_EQUALS(beta1_bound, "flat"  ) ||
      CCTK_EQUALS(beta1_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::beta1", beta1_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register beta1_bound BC for Kranc2BSSNChi::beta1!");
  }
  
  if (CCTK_EQUALS(beta2_bound, "none"  ) ||
      CCTK_EQUALS(beta2_bound, "static") ||
      CCTK_EQUALS(beta2_bound, "flat"  ) ||
      CCTK_EQUALS(beta2_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::beta2", beta2_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register beta2_bound BC for Kranc2BSSNChi::beta2!");
  }
  
  if (CCTK_EQUALS(beta3_bound, "none"  ) ||
      CCTK_EQUALS(beta3_bound, "static") ||
      CCTK_EQUALS(beta3_bound, "flat"  ) ||
      CCTK_EQUALS(beta3_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::beta3", beta3_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register beta3_bound BC for Kranc2BSSNChi::beta3!");
  }
  
  if (CCTK_EQUALS(betat1_bound, "none"  ) ||
      CCTK_EQUALS(betat1_bound, "static") ||
      CCTK_EQUALS(betat1_bound, "flat"  ) ||
      CCTK_EQUALS(betat1_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::betat1", betat1_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register betat1_bound BC for Kranc2BSSNChi::betat1!");
  }
  
  if (CCTK_EQUALS(betat2_bound, "none"  ) ||
      CCTK_EQUALS(betat2_bound, "static") ||
      CCTK_EQUALS(betat2_bound, "flat"  ) ||
      CCTK_EQUALS(betat2_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::betat2", betat2_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register betat2_bound BC for Kranc2BSSNChi::betat2!");
  }
  
  if (CCTK_EQUALS(betat3_bound, "none"  ) ||
      CCTK_EQUALS(betat3_bound, "static") ||
      CCTK_EQUALS(betat3_bound, "flat"  ) ||
      CCTK_EQUALS(betat3_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "Kranc2BSSNChi::betat3", betat3_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register betat3_bound BC for Kranc2BSSNChi::betat3!");
  }
  
  if (CCTK_EQUALS(h_group_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_h_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_h_group_bound < 0) handle_h_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_h_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_h_group_bound , h_group_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_h_group_bound ,h_group_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_h_group_bound, 
                      "Kranc2BSSNChi::h_group", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::h_group!");
  
  }
  
  if (CCTK_EQUALS(A_group_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_A_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_A_group_bound < 0) handle_A_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_A_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_A_group_bound , A_group_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_A_group_bound ,A_group_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_A_group_bound, 
                      "Kranc2BSSNChi::A_group", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::A_group!");
  
  }
  
  if (CCTK_EQUALS(chi_group_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_chi_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_chi_group_bound < 0) handle_chi_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_chi_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_chi_group_bound , chi_group_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_chi_group_bound ,chi_group_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_chi_group_bound, 
                      "Kranc2BSSNChi::chi_group", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::chi_group!");
  
  }
  
  if (CCTK_EQUALS(K_group_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_K_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_K_group_bound < 0) handle_K_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_K_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_K_group_bound , K_group_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_K_group_bound ,K_group_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_K_group_bound, 
                      "Kranc2BSSNChi::K_group", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::K_group!");
  
  }
  
  if (CCTK_EQUALS(Gam_group_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_Gam_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_Gam_group_bound < 0) handle_Gam_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Gam_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_Gam_group_bound , Gam_group_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_Gam_group_bound ,Gam_group_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Gam_group_bound, 
                      "Kranc2BSSNChi::Gam_group", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::Gam_group!");
  
  }
  
  if (CCTK_EQUALS(alpha_group_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_alpha_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_alpha_group_bound < 0) handle_alpha_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_alpha_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_alpha_group_bound , alpha_group_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_alpha_group_bound ,alpha_group_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_alpha_group_bound, 
                      "Kranc2BSSNChi::alpha_group", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::alpha_group!");
  
  }
  
  if (CCTK_EQUALS(beta_group_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_beta_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_beta_group_bound < 0) handle_beta_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_beta_group_bound , beta_group_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_beta_group_bound ,beta_group_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta_group_bound, 
                      "Kranc2BSSNChi::beta_group", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::beta_group!");
  
  }
  
  if (CCTK_EQUALS(betat_group_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_betat_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_betat_group_bound < 0) handle_betat_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_betat_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_betat_group_bound , betat_group_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_betat_group_bound ,betat_group_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_betat_group_bound, 
                      "Kranc2BSSNChi::betat_group", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::betat_group!");
  
  }
  
  if (CCTK_EQUALS(h11_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_h11_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_h11_bound < 0) handle_h11_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_h11_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_h11_bound , h11_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_h11_bound ,h11_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_h11_bound, 
                      "Kranc2BSSNChi::h11", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::h11!");
  
  }
  
  if (CCTK_EQUALS(h21_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_h21_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_h21_bound < 0) handle_h21_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_h21_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_h21_bound , h21_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_h21_bound ,h21_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_h21_bound, 
                      "Kranc2BSSNChi::h21", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::h21!");
  
  }
  
  if (CCTK_EQUALS(h31_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_h31_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_h31_bound < 0) handle_h31_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_h31_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_h31_bound , h31_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_h31_bound ,h31_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_h31_bound, 
                      "Kranc2BSSNChi::h31", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::h31!");
  
  }
  
  if (CCTK_EQUALS(h22_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_h22_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_h22_bound < 0) handle_h22_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_h22_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_h22_bound , h22_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_h22_bound ,h22_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_h22_bound, 
                      "Kranc2BSSNChi::h22", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::h22!");
  
  }
  
  if (CCTK_EQUALS(h32_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_h32_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_h32_bound < 0) handle_h32_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_h32_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_h32_bound , h32_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_h32_bound ,h32_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_h32_bound, 
                      "Kranc2BSSNChi::h32", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::h32!");
  
  }
  
  if (CCTK_EQUALS(h33_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_h33_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_h33_bound < 0) handle_h33_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_h33_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_h33_bound , h33_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_h33_bound ,h33_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_h33_bound, 
                      "Kranc2BSSNChi::h33", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::h33!");
  
  }
  
  if (CCTK_EQUALS(A11_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_A11_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_A11_bound < 0) handle_A11_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_A11_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_A11_bound , A11_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_A11_bound ,A11_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_A11_bound, 
                      "Kranc2BSSNChi::A11", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::A11!");
  
  }
  
  if (CCTK_EQUALS(A21_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_A21_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_A21_bound < 0) handle_A21_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_A21_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_A21_bound , A21_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_A21_bound ,A21_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_A21_bound, 
                      "Kranc2BSSNChi::A21", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::A21!");
  
  }
  
  if (CCTK_EQUALS(A31_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_A31_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_A31_bound < 0) handle_A31_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_A31_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_A31_bound , A31_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_A31_bound ,A31_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_A31_bound, 
                      "Kranc2BSSNChi::A31", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::A31!");
  
  }
  
  if (CCTK_EQUALS(A22_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_A22_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_A22_bound < 0) handle_A22_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_A22_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_A22_bound , A22_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_A22_bound ,A22_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_A22_bound, 
                      "Kranc2BSSNChi::A22", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::A22!");
  
  }
  
  if (CCTK_EQUALS(A32_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_A32_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_A32_bound < 0) handle_A32_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_A32_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_A32_bound , A32_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_A32_bound ,A32_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_A32_bound, 
                      "Kranc2BSSNChi::A32", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::A32!");
  
  }
  
  if (CCTK_EQUALS(A33_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_A33_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_A33_bound < 0) handle_A33_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_A33_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_A33_bound , A33_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_A33_bound ,A33_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_A33_bound, 
                      "Kranc2BSSNChi::A33", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::A33!");
  
  }
  
  if (CCTK_EQUALS(chi_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_chi_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_chi_bound < 0) handle_chi_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_chi_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_chi_bound , chi_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_chi_bound ,chi_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_chi_bound, 
                      "Kranc2BSSNChi::chi", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::chi!");
  
  }
  
  if (CCTK_EQUALS(K_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_K_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_K_bound < 0) handle_K_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_K_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_K_bound , K_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_K_bound ,K_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_K_bound, 
                      "Kranc2BSSNChi::K", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::K!");
  
  }
  
  if (CCTK_EQUALS(Gam1_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_Gam1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_Gam1_bound < 0) handle_Gam1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Gam1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_Gam1_bound , Gam1_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_Gam1_bound ,Gam1_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Gam1_bound, 
                      "Kranc2BSSNChi::Gam1", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::Gam1!");
  
  }
  
  if (CCTK_EQUALS(Gam2_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_Gam2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_Gam2_bound < 0) handle_Gam2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Gam2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_Gam2_bound , Gam2_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_Gam2_bound ,Gam2_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Gam2_bound, 
                      "Kranc2BSSNChi::Gam2", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::Gam2!");
  
  }
  
  if (CCTK_EQUALS(Gam3_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_Gam3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_Gam3_bound < 0) handle_Gam3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Gam3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_Gam3_bound , Gam3_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_Gam3_bound ,Gam3_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Gam3_bound, 
                      "Kranc2BSSNChi::Gam3", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::Gam3!");
  
  }
  
  if (CCTK_EQUALS(alpha_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_alpha_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_alpha_bound < 0) handle_alpha_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_alpha_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_alpha_bound , alpha_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_alpha_bound ,alpha_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_alpha_bound, 
                      "Kranc2BSSNChi::alpha", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::alpha!");
  
  }
  
  if (CCTK_EQUALS(beta1_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_beta1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_beta1_bound < 0) handle_beta1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_beta1_bound , beta1_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_beta1_bound ,beta1_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta1_bound, 
                      "Kranc2BSSNChi::beta1", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::beta1!");
  
  }
  
  if (CCTK_EQUALS(beta2_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_beta2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_beta2_bound < 0) handle_beta2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_beta2_bound , beta2_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_beta2_bound ,beta2_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta2_bound, 
                      "Kranc2BSSNChi::beta2", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::beta2!");
  
  }
  
  if (CCTK_EQUALS(beta3_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_beta3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_beta3_bound < 0) handle_beta3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_beta3_bound , beta3_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_beta3_bound ,beta3_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta3_bound, 
                      "Kranc2BSSNChi::beta3", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::beta3!");
  
  }
  
  if (CCTK_EQUALS(betat1_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_betat1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_betat1_bound < 0) handle_betat1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_betat1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_betat1_bound , betat1_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_betat1_bound ,betat1_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_betat1_bound, 
                      "Kranc2BSSNChi::betat1", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::betat1!");
  
  }
  
  if (CCTK_EQUALS(betat2_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_betat2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_betat2_bound < 0) handle_betat2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_betat2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_betat2_bound , betat2_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_betat2_bound ,betat2_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_betat2_bound, 
                      "Kranc2BSSNChi::betat2", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::betat2!");
  
  }
  
  if (CCTK_EQUALS(betat3_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_betat3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_betat3_bound < 0) handle_betat3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_betat3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_betat3_bound , betat3_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_betat3_bound ,betat3_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_betat3_bound, 
                      "Kranc2BSSNChi::betat3", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for Kranc2BSSNChi::betat3!");
  
  }
  
  if (CCTK_EQUALS(h_group_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_h_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_h_group_bound < 0) handle_h_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_h_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_h_group_bound ,h_group_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_h_group_bound, 
                      "Kranc2BSSNChi::h_group", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for Kranc2BSSNChi::h_group!");
  
  }
  
  if (CCTK_EQUALS(A_group_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_A_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_A_group_bound < 0) handle_A_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_A_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_A_group_bound ,A_group_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_A_group_bound, 
                      "Kranc2BSSNChi::A_group", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for Kranc2BSSNChi::A_group!");
  
  }
  
  if (CCTK_EQUALS(chi_group_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_chi_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_chi_group_bound < 0) handle_chi_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_chi_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_chi_group_bound ,chi_group_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_chi_group_bound, 
                      "Kranc2BSSNChi::chi_group", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for Kranc2BSSNChi::chi_group!");
  
  }
  
  if (CCTK_EQUALS(K_group_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_K_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_K_group_bound < 0) handle_K_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_K_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_K_group_bound ,K_group_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_K_group_bound, 
                      "Kranc2BSSNChi::K_group", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for Kranc2BSSNChi::K_group!");
  
  }
  
  if (CCTK_EQUALS(Gam_group_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_Gam_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_Gam_group_bound < 0) handle_Gam_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Gam_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_Gam_group_bound ,Gam_group_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Gam_group_bound, 
                      "Kranc2BSSNChi::Gam_group", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for Kranc2BSSNChi::Gam_group!");
  
  }
  
  if (CCTK_EQUALS(alpha_group_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_alpha_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_alpha_group_bound < 0) handle_alpha_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_alpha_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_alpha_group_bound ,alpha_group_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_alpha_group_bound, 
                      "Kranc2BSSNChi::alpha_group", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for Kranc2BSSNChi::alpha_group!");
  
  }
  
  if (CCTK_EQUALS(beta_group_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_beta_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_beta_group_bound < 0) handle_beta_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_beta_group_bound ,beta_group_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta_group_bound, 
                      "Kranc2BSSNChi::beta_group", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for Kranc2BSSNChi::beta_group!");
  
  }
  
  if (CCTK_EQUALS(betat_group_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_betat_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_betat_group_bound < 0) handle_betat_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_betat_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_betat_group_bound ,betat_group_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_betat_group_bound, 
                      "Kranc2BSSNChi::betat_group", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for Kranc2BSSNChi::betat_group!");
  
  }
  
  if (CCTK_EQUALS(h11_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_h11_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_h11_bound < 0) handle_h11_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_h11_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_h11_bound ,h11_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_h11_bound, 
                      "Kranc2BSSNChi::h11", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for Kranc2BSSNChi::h11!");
  
  }
  
  if (CCTK_EQUALS(h21_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_h21_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_h21_bound < 0) handle_h21_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_h21_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_h21_bound ,h21_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_h21_bound, 
                      "Kranc2BSSNChi::h21", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for Kranc2BSSNChi::h21!");
  
  }
  
  if (CCTK_EQUALS(h31_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_h31_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_h31_bound < 0) handle_h31_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_h31_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_h31_bound ,h31_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_h31_bound, 
                      "Kranc2BSSNChi::h31", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for Kranc2BSSNChi::h31!");
  
  }
  
  if (CCTK_EQUALS(h22_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_h22_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_h22_bound < 0) handle_h22_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_h22_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_h22_bound ,h22_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_h22_bound, 
                      "Kranc2BSSNChi::h22", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for Kranc2BSSNChi::h22!");
  
  }
  
  if (CCTK_EQUALS(h32_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_h32_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_h32_bound < 0) handle_h32_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_h32_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_h32_bound ,h32_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_h32_bound, 
                      "Kranc2BSSNChi::h32", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for Kranc2BSSNChi::h32!");
  
  }
  
  if (CCTK_EQUALS(h33_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_h33_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_h33_bound < 0) handle_h33_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_h33_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_h33_bound ,h33_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_h33_bound, 
                      "Kranc2BSSNChi::h33", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for Kranc2BSSNChi::h33!");
  
  }
  
  if (CCTK_EQUALS(A11_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_A11_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_A11_bound < 0) handle_A11_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_A11_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_A11_bound ,A11_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_A11_bound, 
                      "Kranc2BSSNChi::A11", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for Kranc2BSSNChi::A11!");
  
  }
  
  if (CCTK_EQUALS(A21_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_A21_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_A21_bound < 0) handle_A21_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_A21_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_A21_bound ,A21_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_A21_bound, 
                      "Kranc2BSSNChi::A21", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for Kranc2BSSNChi::A21!");
  
  }
  
  if (CCTK_EQUALS(A31_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_A31_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_A31_bound < 0) handle_A31_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_A31_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_A31_bound ,A31_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_A31_bound, 
                      "Kranc2BSSNChi::A31", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for Kranc2BSSNChi::A31!");
  
  }
  
  if (CCTK_EQUALS(A22_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_A22_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_A22_bound < 0) handle_A22_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_A22_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_A22_bound ,A22_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_A22_bound, 
                      "Kranc2BSSNChi::A22", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for Kranc2BSSNChi::A22!");
  
  }
  
  if (CCTK_EQUALS(A32_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_A32_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_A32_bound < 0) handle_A32_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_A32_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_A32_bound ,A32_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_A32_bound, 
                      "Kranc2BSSNChi::A32", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for Kranc2BSSNChi::A32!");
  
  }
  
  if (CCTK_EQUALS(A33_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_A33_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_A33_bound < 0) handle_A33_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_A33_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_A33_bound ,A33_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_A33_bound, 
                      "Kranc2BSSNChi::A33", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for Kranc2BSSNChi::A33!");
  
  }
  
  if (CCTK_EQUALS(chi_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_chi_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_chi_bound < 0) handle_chi_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_chi_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_chi_bound ,chi_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_chi_bound, 
                      "Kranc2BSSNChi::chi", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for Kranc2BSSNChi::chi!");
  
  }
  
  if (CCTK_EQUALS(K_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_K_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_K_bound < 0) handle_K_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_K_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_K_bound ,K_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_K_bound, 
                      "Kranc2BSSNChi::K", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for Kranc2BSSNChi::K!");
  
  }
  
  if (CCTK_EQUALS(Gam1_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_Gam1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_Gam1_bound < 0) handle_Gam1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Gam1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_Gam1_bound ,Gam1_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Gam1_bound, 
                      "Kranc2BSSNChi::Gam1", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for Kranc2BSSNChi::Gam1!");
  
  }
  
  if (CCTK_EQUALS(Gam2_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_Gam2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_Gam2_bound < 0) handle_Gam2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Gam2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_Gam2_bound ,Gam2_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Gam2_bound, 
                      "Kranc2BSSNChi::Gam2", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for Kranc2BSSNChi::Gam2!");
  
  }
  
  if (CCTK_EQUALS(Gam3_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_Gam3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_Gam3_bound < 0) handle_Gam3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Gam3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_Gam3_bound ,Gam3_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Gam3_bound, 
                      "Kranc2BSSNChi::Gam3", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for Kranc2BSSNChi::Gam3!");
  
  }
  
  if (CCTK_EQUALS(alpha_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_alpha_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_alpha_bound < 0) handle_alpha_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_alpha_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_alpha_bound ,alpha_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_alpha_bound, 
                      "Kranc2BSSNChi::alpha", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for Kranc2BSSNChi::alpha!");
  
  }
  
  if (CCTK_EQUALS(beta1_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_beta1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_beta1_bound < 0) handle_beta1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_beta1_bound ,beta1_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta1_bound, 
                      "Kranc2BSSNChi::beta1", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for Kranc2BSSNChi::beta1!");
  
  }
  
  if (CCTK_EQUALS(beta2_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_beta2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_beta2_bound < 0) handle_beta2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_beta2_bound ,beta2_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta2_bound, 
                      "Kranc2BSSNChi::beta2", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for Kranc2BSSNChi::beta2!");
  
  }
  
  if (CCTK_EQUALS(beta3_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_beta3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_beta3_bound < 0) handle_beta3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_beta3_bound ,beta3_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta3_bound, 
                      "Kranc2BSSNChi::beta3", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for Kranc2BSSNChi::beta3!");
  
  }
  
  if (CCTK_EQUALS(betat1_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_betat1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_betat1_bound < 0) handle_betat1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_betat1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_betat1_bound ,betat1_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_betat1_bound, 
                      "Kranc2BSSNChi::betat1", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for Kranc2BSSNChi::betat1!");
  
  }
  
  if (CCTK_EQUALS(betat2_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_betat2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_betat2_bound < 0) handle_betat2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_betat2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_betat2_bound ,betat2_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_betat2_bound, 
                      "Kranc2BSSNChi::betat2", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for Kranc2BSSNChi::betat2!");
  
  }
  
  if (CCTK_EQUALS(betat3_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_betat3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_betat3_bound < 0) handle_betat3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_betat3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_betat3_bound ,betat3_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_betat3_bound, 
                      "Kranc2BSSNChi::betat3", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for Kranc2BSSNChi::betat3!");
  
  }
  return;
}



/* template for entries in parameter file:
#$bound$#Kranc2BSSNChi::h_group_bound       = "skip"
#$bound$#Kranc2BSSNChi::h_group_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::h_group_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::h_group_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::A_group_bound       = "skip"
#$bound$#Kranc2BSSNChi::A_group_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::A_group_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::A_group_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::chi_group_bound       = "skip"
#$bound$#Kranc2BSSNChi::chi_group_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::chi_group_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::chi_group_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::K_group_bound       = "skip"
#$bound$#Kranc2BSSNChi::K_group_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::K_group_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::K_group_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::Gam_group_bound       = "skip"
#$bound$#Kranc2BSSNChi::Gam_group_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::Gam_group_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::Gam_group_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::alpha_group_bound       = "skip"
#$bound$#Kranc2BSSNChi::alpha_group_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::alpha_group_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::alpha_group_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::beta_group_bound       = "skip"
#$bound$#Kranc2BSSNChi::beta_group_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::beta_group_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::beta_group_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::betat_group_bound       = "skip"
#$bound$#Kranc2BSSNChi::betat_group_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::betat_group_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::betat_group_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::h11_bound       = "skip"
#$bound$#Kranc2BSSNChi::h11_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::h11_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::h11_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::h21_bound       = "skip"
#$bound$#Kranc2BSSNChi::h21_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::h21_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::h21_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::h31_bound       = "skip"
#$bound$#Kranc2BSSNChi::h31_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::h31_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::h31_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::h22_bound       = "skip"
#$bound$#Kranc2BSSNChi::h22_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::h22_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::h22_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::h32_bound       = "skip"
#$bound$#Kranc2BSSNChi::h32_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::h32_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::h32_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::h33_bound       = "skip"
#$bound$#Kranc2BSSNChi::h33_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::h33_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::h33_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::A11_bound       = "skip"
#$bound$#Kranc2BSSNChi::A11_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::A11_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::A11_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::A21_bound       = "skip"
#$bound$#Kranc2BSSNChi::A21_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::A21_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::A21_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::A31_bound       = "skip"
#$bound$#Kranc2BSSNChi::A31_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::A31_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::A31_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::A22_bound       = "skip"
#$bound$#Kranc2BSSNChi::A22_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::A22_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::A22_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::A32_bound       = "skip"
#$bound$#Kranc2BSSNChi::A32_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::A32_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::A32_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::A33_bound       = "skip"
#$bound$#Kranc2BSSNChi::A33_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::A33_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::A33_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::chi_bound       = "skip"
#$bound$#Kranc2BSSNChi::chi_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::chi_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::chi_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::K_bound       = "skip"
#$bound$#Kranc2BSSNChi::K_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::K_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::K_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::Gam1_bound       = "skip"
#$bound$#Kranc2BSSNChi::Gam1_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::Gam1_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::Gam1_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::Gam2_bound       = "skip"
#$bound$#Kranc2BSSNChi::Gam2_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::Gam2_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::Gam2_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::Gam3_bound       = "skip"
#$bound$#Kranc2BSSNChi::Gam3_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::Gam3_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::Gam3_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::alpha_bound       = "skip"
#$bound$#Kranc2BSSNChi::alpha_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::alpha_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::alpha_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::beta1_bound       = "skip"
#$bound$#Kranc2BSSNChi::beta1_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::beta1_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::beta1_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::beta2_bound       = "skip"
#$bound$#Kranc2BSSNChi::beta2_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::beta2_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::beta2_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::beta3_bound       = "skip"
#$bound$#Kranc2BSSNChi::beta3_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::beta3_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::beta3_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::betat1_bound       = "skip"
#$bound$#Kranc2BSSNChi::betat1_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::betat1_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::betat1_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::betat2_bound       = "skip"
#$bound$#Kranc2BSSNChi::betat2_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::betat2_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::betat2_bound_scalar = 0.0

#$bound$#Kranc2BSSNChi::betat3_bound       = "skip"
#$bound$#Kranc2BSSNChi::betat3_bound_speed = 1.0
#$bound$#Kranc2BSSNChi::betat3_bound_limit = 0.0
#$bound$#Kranc2BSSNChi::betat3_bound_scalar = 0.0

*/

