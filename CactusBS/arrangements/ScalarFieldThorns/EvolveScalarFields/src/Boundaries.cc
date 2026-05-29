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


extern "C" void EvolveScalarFields_CheckBoundaries(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  return;
}

extern "C" void EvolveScalarFields_SelectBoundConds(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  
  if (CCTK_EQUALS(phia_group_bound, "none"  ) ||
      CCTK_EQUALS(phia_group_bound, "static") ||
      CCTK_EQUALS(phia_group_bound, "flat"  ) ||
      CCTK_EQUALS(phia_group_bound, "zero"  ) )
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "EvolveScalarFields::phia_group", phia_group_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register phia_group_bound BC for EvolveScalarFields::phia_group!");
  }
  
  if (CCTK_EQUALS(pia_group_bound, "none"  ) ||
      CCTK_EQUALS(pia_group_bound, "static") ||
      CCTK_EQUALS(pia_group_bound, "flat"  ) ||
      CCTK_EQUALS(pia_group_bound, "zero"  ) )
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "EvolveScalarFields::pia_group", pia_group_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register pia_group_bound BC for EvolveScalarFields::pia_group!");
  }
  
  if (CCTK_EQUALS(phib_group_bound, "none"  ) ||
      CCTK_EQUALS(phib_group_bound, "static") ||
      CCTK_EQUALS(phib_group_bound, "flat"  ) ||
      CCTK_EQUALS(phib_group_bound, "zero"  ) )
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "EvolveScalarFields::phib_group", phib_group_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register phib_group_bound BC for EvolveScalarFields::phib_group!");
  }
  
  if (CCTK_EQUALS(pib_group_bound, "none"  ) ||
      CCTK_EQUALS(pib_group_bound, "static") ||
      CCTK_EQUALS(pib_group_bound, "flat"  ) ||
      CCTK_EQUALS(pib_group_bound, "zero"  ) )
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "EvolveScalarFields::pib_group", pib_group_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register pib_group_bound BC for EvolveScalarFields::pib_group!");
  }
  
  if (CCTK_EQUALS(phia_bound, "none"  ) ||
      CCTK_EQUALS(phia_bound, "static") ||
      CCTK_EQUALS(phia_bound, "flat"  ) ||
      CCTK_EQUALS(phia_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "EvolveScalarFields::phia", phia_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register phia_bound BC for EvolveScalarFields::phia!");
  }
  
  if (CCTK_EQUALS(pia_bound, "none"  ) ||
      CCTK_EQUALS(pia_bound, "static") ||
      CCTK_EQUALS(pia_bound, "flat"  ) ||
      CCTK_EQUALS(pia_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "EvolveScalarFields::pia", pia_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register pia_bound BC for EvolveScalarFields::pia!");
  }
  
  if (CCTK_EQUALS(phib_bound, "none"  ) ||
      CCTK_EQUALS(phib_bound, "static") ||
      CCTK_EQUALS(phib_bound, "flat"  ) ||
      CCTK_EQUALS(phib_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "EvolveScalarFields::phib", phib_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register phib_bound BC for EvolveScalarFields::phib!");
  }
  
  if (CCTK_EQUALS(pib_bound, "none"  ) ||
      CCTK_EQUALS(pib_bound, "static") ||
      CCTK_EQUALS(pib_bound, "flat"  ) ||
      CCTK_EQUALS(pib_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "EvolveScalarFields::pib", pib_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register pib_bound BC for EvolveScalarFields::pib!");
  }
  
  if (CCTK_EQUALS(phia_group_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_phia_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_phia_group_bound < 0) handle_phia_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_phia_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_phia_group_bound , phia_group_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_phia_group_bound ,phia_group_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_phia_group_bound, 
                      "EvolveScalarFields::phia_group", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for EvolveScalarFields::phia_group!");
  
  }
  
  if (CCTK_EQUALS(pia_group_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_pia_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_pia_group_bound < 0) handle_pia_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_pia_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_pia_group_bound , pia_group_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_pia_group_bound ,pia_group_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_pia_group_bound, 
                      "EvolveScalarFields::pia_group", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for EvolveScalarFields::pia_group!");
  
  }
  
  if (CCTK_EQUALS(phib_group_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_phib_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_phib_group_bound < 0) handle_phib_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_phib_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_phib_group_bound , phib_group_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_phib_group_bound ,phib_group_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_phib_group_bound, 
                      "EvolveScalarFields::phib_group", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for EvolveScalarFields::phib_group!");
  
  }
  
  if (CCTK_EQUALS(pib_group_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_pib_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_pib_group_bound < 0) handle_pib_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_pib_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_pib_group_bound , pib_group_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_pib_group_bound ,pib_group_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_pib_group_bound, 
                      "EvolveScalarFields::pib_group", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for EvolveScalarFields::pib_group!");
  
  }
  
  if (CCTK_EQUALS(phia_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_phia_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_phia_bound < 0) handle_phia_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_phia_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_phia_bound , phia_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_phia_bound ,phia_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_phia_bound, 
                      "EvolveScalarFields::phia", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for EvolveScalarFields::phia!");
  
  }
  
  if (CCTK_EQUALS(pia_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_pia_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_pia_bound < 0) handle_pia_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_pia_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_pia_bound , pia_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_pia_bound ,pia_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_pia_bound, 
                      "EvolveScalarFields::pia", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for EvolveScalarFields::pia!");
  
  }
  
  if (CCTK_EQUALS(phib_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_phib_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_phib_bound < 0) handle_phib_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_phib_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_phib_bound , phib_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_phib_bound ,phib_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_phib_bound, 
                      "EvolveScalarFields::phib", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for EvolveScalarFields::phib!");
  
  }
  
  if (CCTK_EQUALS(pib_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_pib_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_pib_bound < 0) handle_pib_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_pib_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_pib_bound , pib_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_pib_bound ,pib_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_pib_bound, 
                      "EvolveScalarFields::pib", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for EvolveScalarFields::pib!");
  
  }
  
  if (CCTK_EQUALS(phia_group_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_phia_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_phia_group_bound < 0) handle_phia_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_phia_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_phia_group_bound ,phia_group_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_phia_group_bound, 
                      "EvolveScalarFields::phia_group", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for EvolveScalarFields::phia_group!");
  
  }
  
  if (CCTK_EQUALS(pia_group_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_pia_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_pia_group_bound < 0) handle_pia_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_pia_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_pia_group_bound ,pia_group_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_pia_group_bound, 
                      "EvolveScalarFields::pia_group", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for EvolveScalarFields::pia_group!");
  
  }
  
  if (CCTK_EQUALS(phib_group_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_phib_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_phib_group_bound < 0) handle_phib_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_phib_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_phib_group_bound ,phib_group_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_phib_group_bound, 
                      "EvolveScalarFields::phib_group", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for EvolveScalarFields::phib_group!");
  
  }
  
  if (CCTK_EQUALS(pib_group_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_pib_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_pib_group_bound < 0) handle_pib_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_pib_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_pib_group_bound ,pib_group_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_pib_group_bound, 
                      "EvolveScalarFields::pib_group", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for EvolveScalarFields::pib_group!");
  
  }
  
  if (CCTK_EQUALS(phia_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_phia_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_phia_bound < 0) handle_phia_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_phia_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_phia_bound ,phia_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_phia_bound, 
                      "EvolveScalarFields::phia", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for EvolveScalarFields::phia!");
  
  }
  
  if (CCTK_EQUALS(pia_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_pia_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_pia_bound < 0) handle_pia_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_pia_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_pia_bound ,pia_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_pia_bound, 
                      "EvolveScalarFields::pia", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for EvolveScalarFields::pia!");
  
  }
  
  if (CCTK_EQUALS(phib_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_phib_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_phib_bound < 0) handle_phib_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_phib_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_phib_bound ,phib_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_phib_bound, 
                      "EvolveScalarFields::phib", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for EvolveScalarFields::phib!");
  
  }
  
  if (CCTK_EQUALS(pib_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_pib_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_pib_bound < 0) handle_pib_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_pib_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_pib_bound ,pib_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_pib_bound, 
                      "EvolveScalarFields::pib", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for EvolveScalarFields::pib!");
  
  }
  return;
}



/* template for entries in parameter file:
#$bound$#EvolveScalarFields::phia_group_bound       = "skip"
#$bound$#EvolveScalarFields::phia_group_bound_speed = 1.0
#$bound$#EvolveScalarFields::phia_group_bound_limit = 0.0
#$bound$#EvolveScalarFields::phia_group_bound_scalar = 0.0

#$bound$#EvolveScalarFields::pia_group_bound       = "skip"
#$bound$#EvolveScalarFields::pia_group_bound_speed = 1.0
#$bound$#EvolveScalarFields::pia_group_bound_limit = 0.0
#$bound$#EvolveScalarFields::pia_group_bound_scalar = 0.0

#$bound$#EvolveScalarFields::phib_group_bound       = "skip"
#$bound$#EvolveScalarFields::phib_group_bound_speed = 1.0
#$bound$#EvolveScalarFields::phib_group_bound_limit = 0.0
#$bound$#EvolveScalarFields::phib_group_bound_scalar = 0.0

#$bound$#EvolveScalarFields::pib_group_bound       = "skip"
#$bound$#EvolveScalarFields::pib_group_bound_speed = 1.0
#$bound$#EvolveScalarFields::pib_group_bound_limit = 0.0
#$bound$#EvolveScalarFields::pib_group_bound_scalar = 0.0

#$bound$#EvolveScalarFields::phia_bound       = "skip"
#$bound$#EvolveScalarFields::phia_bound_speed = 1.0
#$bound$#EvolveScalarFields::phia_bound_limit = 0.0
#$bound$#EvolveScalarFields::phia_bound_scalar = 0.0

#$bound$#EvolveScalarFields::pia_bound       = "skip"
#$bound$#EvolveScalarFields::pia_bound_speed = 1.0
#$bound$#EvolveScalarFields::pia_bound_limit = 0.0
#$bound$#EvolveScalarFields::pia_bound_scalar = 0.0

#$bound$#EvolveScalarFields::phib_bound       = "skip"
#$bound$#EvolveScalarFields::phib_bound_speed = 1.0
#$bound$#EvolveScalarFields::phib_bound_limit = 0.0
#$bound$#EvolveScalarFields::phib_bound_scalar = 0.0

#$bound$#EvolveScalarFields::pib_bound       = "skip"
#$bound$#EvolveScalarFields::pib_bound_speed = 1.0
#$bound$#EvolveScalarFields::pib_bound_limit = 0.0
#$bound$#EvolveScalarFields::pib_bound_scalar = 0.0

*/

