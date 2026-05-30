/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Faces.h"
#include "util_Table.h"
#include "Symmetry.h"
#include "Kranc.hh"


/* the boundary treatment is split into 3 steps:    */
/* 1. excision                                      */
/* 2. symmetries                                    */
/* 3. "other" boundary conditions, e.g. radiative */

/* to simplify scheduling and testing, the 3 steps  */
/* are currently applied in separate functions      */


extern "C" void ML_WaveToy_CheckBoundaries(CCTK_ARGUMENTS)
{
#ifdef DECLARE_CCTK_ARGUMENTS_ML_WaveToy_CheckBoundaries
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_WaveToy_CheckBoundaries);
#else
  DECLARE_CCTK_ARGUMENTS;
#endif
  DECLARE_CCTK_PARAMETERS;
  
  return;
}

extern "C" void ML_WaveToy_SelectBoundConds(CCTK_ARGUMENTS)
{
#ifdef DECLARE_CCTK_ARGUMENTS_ML_WaveToy_SelectBoundConds
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_WaveToy_SelectBoundConds);
#else
  DECLARE_CCTK_ARGUMENTS;
#endif
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  
  if (CCTK_EQUALS(WT_u_bound, "none"  ) ||
      CCTK_EQUALS(WT_u_bound, "static") ||
      CCTK_EQUALS(WT_u_bound, "flat"  ) ||
      CCTK_EQUALS(WT_u_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_WaveToy::WT_u", WT_u_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register WT_u_bound BC for ML_WaveToy::WT_u!");
  }
  
  if (CCTK_EQUALS(WT_rho_bound, "none"  ) ||
      CCTK_EQUALS(WT_rho_bound, "static") ||
      CCTK_EQUALS(WT_rho_bound, "flat"  ) ||
      CCTK_EQUALS(WT_rho_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_WaveToy::WT_rho", WT_rho_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register WT_rho_bound BC for ML_WaveToy::WT_rho!");
  }
  
  if (CCTK_EQUALS(u_bound, "none"  ) ||
      CCTK_EQUALS(u_bound, "static") ||
      CCTK_EQUALS(u_bound, "flat"  ) ||
      CCTK_EQUALS(u_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_WaveToy::u", u_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register u_bound BC for ML_WaveToy::u!");
  }
  
  if (CCTK_EQUALS(rho_bound, "none"  ) ||
      CCTK_EQUALS(rho_bound, "static") ||
      CCTK_EQUALS(rho_bound, "flat"  ) ||
      CCTK_EQUALS(rho_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_WaveToy::rho", rho_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register rho_bound BC for ML_WaveToy::rho!");
  }
  
  if (CCTK_EQUALS(WT_u_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_WT_u_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_WT_u_bound < 0) handle_WT_u_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_WT_u_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_WT_u_bound , WT_u_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_WT_u_bound ,WT_u_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_WT_u_bound, 
                      "ML_WaveToy::WT_u", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_WaveToy::WT_u!");
  
  }
  
  if (CCTK_EQUALS(WT_rho_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_WT_rho_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_WT_rho_bound < 0) handle_WT_rho_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_WT_rho_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_WT_rho_bound , WT_rho_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_WT_rho_bound ,WT_rho_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_WT_rho_bound, 
                      "ML_WaveToy::WT_rho", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_WaveToy::WT_rho!");
  
  }
  
  if (CCTK_EQUALS(u_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_u_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_u_bound < 0) handle_u_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_u_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_u_bound , u_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_u_bound ,u_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_u_bound, 
                      "ML_WaveToy::u", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_WaveToy::u!");
  
  }
  
  if (CCTK_EQUALS(rho_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_rho_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_rho_bound < 0) handle_rho_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_rho_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_rho_bound , rho_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_rho_bound ,rho_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_rho_bound, 
                      "ML_WaveToy::rho", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_WaveToy::rho!");
  
  }
  
  if (CCTK_EQUALS(WT_u_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_WT_u_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_WT_u_bound < 0) handle_WT_u_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_WT_u_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_WT_u_bound ,WT_u_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_WT_u_bound, 
                      "ML_WaveToy::WT_u", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for ML_WaveToy::WT_u!");
  
  }
  
  if (CCTK_EQUALS(WT_rho_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_WT_rho_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_WT_rho_bound < 0) handle_WT_rho_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_WT_rho_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_WT_rho_bound ,WT_rho_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_WT_rho_bound, 
                      "ML_WaveToy::WT_rho", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for ML_WaveToy::WT_rho!");
  
  }
  
  if (CCTK_EQUALS(u_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_u_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_u_bound < 0) handle_u_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_u_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_u_bound ,u_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_u_bound, 
                      "ML_WaveToy::u", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_WaveToy::u!");
  
  }
  
  if (CCTK_EQUALS(rho_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_rho_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_rho_bound < 0) handle_rho_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_rho_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_rho_bound ,rho_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_rho_bound, 
                      "ML_WaveToy::rho", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_WaveToy::rho!");
  
  }
  return;
}



/* template for entries in parameter file:
#$bound$#ML_WaveToy::WT_u_bound       = "skip"
#$bound$#ML_WaveToy::WT_u_bound_speed = 1.0
#$bound$#ML_WaveToy::WT_u_bound_limit = 0.0
#$bound$#ML_WaveToy::WT_u_bound_scalar = 0.0

#$bound$#ML_WaveToy::WT_rho_bound       = "skip"
#$bound$#ML_WaveToy::WT_rho_bound_speed = 1.0
#$bound$#ML_WaveToy::WT_rho_bound_limit = 0.0
#$bound$#ML_WaveToy::WT_rho_bound_scalar = 0.0

#$bound$#ML_WaveToy::u_bound       = "skip"
#$bound$#ML_WaveToy::u_bound_speed = 1.0
#$bound$#ML_WaveToy::u_bound_limit = 0.0
#$bound$#ML_WaveToy::u_bound_scalar = 0.0

#$bound$#ML_WaveToy::rho_bound       = "skip"
#$bound$#ML_WaveToy::rho_bound_speed = 1.0
#$bound$#ML_WaveToy::rho_bound_limit = 0.0
#$bound$#ML_WaveToy::rho_bound_scalar = 0.0

*/

