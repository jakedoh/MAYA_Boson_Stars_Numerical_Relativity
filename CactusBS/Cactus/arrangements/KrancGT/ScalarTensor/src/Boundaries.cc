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


extern "C" void ScalarTensor_CheckBoundaries(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  return;
}

extern "C" void ScalarTensor_SelectBoundConds(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  
  if (CCTK_EQUALS(STphi_group_bound, "none"  ) ||
      CCTK_EQUALS(STphi_group_bound, "static") ||
      CCTK_EQUALS(STphi_group_bound, "flat"  ) ||
      CCTK_EQUALS(STphi_group_bound, "zero"  ) )
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ScalarTensor::STphi_group", STphi_group_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register STphi_group_bound BC for ScalarTensor::STphi_group!");
  }
  
  if (CCTK_EQUALS(STpi_group_bound, "none"  ) ||
      CCTK_EQUALS(STpi_group_bound, "static") ||
      CCTK_EQUALS(STpi_group_bound, "flat"  ) ||
      CCTK_EQUALS(STpi_group_bound, "zero"  ) )
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ScalarTensor::STpi_group", STpi_group_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register STpi_group_bound BC for ScalarTensor::STpi_group!");
  }
  
  if (CCTK_EQUALS(STphi_bound, "none"  ) ||
      CCTK_EQUALS(STphi_bound, "static") ||
      CCTK_EQUALS(STphi_bound, "flat"  ) ||
      CCTK_EQUALS(STphi_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ScalarTensor::STphi", STphi_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register STphi_bound BC for ScalarTensor::STphi!");
  }
  
  if (CCTK_EQUALS(STpi_bound, "none"  ) ||
      CCTK_EQUALS(STpi_bound, "static") ||
      CCTK_EQUALS(STpi_bound, "flat"  ) ||
      CCTK_EQUALS(STpi_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ScalarTensor::STpi", STpi_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register STpi_bound BC for ScalarTensor::STpi!");
  }
  
  if (CCTK_EQUALS(STphi_group_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_STphi_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_STphi_group_bound < 0) handle_STphi_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_STphi_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_STphi_group_bound , STphi_group_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_STphi_group_bound ,STphi_group_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_STphi_group_bound, 
                      "ScalarTensor::STphi_group", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ScalarTensor::STphi_group!");
  
  }
  
  if (CCTK_EQUALS(STpi_group_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_STpi_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_STpi_group_bound < 0) handle_STpi_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_STpi_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_STpi_group_bound , STpi_group_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_STpi_group_bound ,STpi_group_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_STpi_group_bound, 
                      "ScalarTensor::STpi_group", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ScalarTensor::STpi_group!");
  
  }
  
  if (CCTK_EQUALS(STphi_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_STphi_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_STphi_bound < 0) handle_STphi_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_STphi_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_STphi_bound , STphi_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_STphi_bound ,STphi_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_STphi_bound, 
                      "ScalarTensor::STphi", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ScalarTensor::STphi!");
  
  }
  
  if (CCTK_EQUALS(STpi_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_STpi_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_STpi_bound < 0) handle_STpi_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_STpi_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_STpi_bound , STpi_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_STpi_bound ,STpi_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_STpi_bound, 
                      "ScalarTensor::STpi", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ScalarTensor::STpi!");
  
  }
  
  if (CCTK_EQUALS(STphi_group_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_STphi_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_STphi_group_bound < 0) handle_STphi_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_STphi_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_STphi_group_bound ,STphi_group_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_STphi_group_bound, 
                      "ScalarTensor::STphi_group", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for ScalarTensor::STphi_group!");
  
  }
  
  if (CCTK_EQUALS(STpi_group_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_STpi_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_STpi_group_bound < 0) handle_STpi_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_STpi_group_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_STpi_group_bound ,STpi_group_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_STpi_group_bound, 
                      "ScalarTensor::STpi_group", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for ScalarTensor::STpi_group!");
  
  }
  
  if (CCTK_EQUALS(STphi_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_STphi_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_STphi_bound < 0) handle_STphi_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_STphi_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_STphi_bound ,STphi_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_STphi_bound, 
                      "ScalarTensor::STphi", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ScalarTensor::STphi!");
  
  }
  
  if (CCTK_EQUALS(STpi_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_STpi_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_STpi_bound < 0) handle_STpi_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_STpi_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_STpi_bound ,STpi_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_STpi_bound, 
                      "ScalarTensor::STpi", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ScalarTensor::STpi!");
  
  }
  return;
}



/* template for entries in parameter file:
#$bound$#ScalarTensor::STphi_group_bound       = "skip"
#$bound$#ScalarTensor::STphi_group_bound_speed = 1.0
#$bound$#ScalarTensor::STphi_group_bound_limit = 0.0
#$bound$#ScalarTensor::STphi_group_bound_scalar = 0.0

#$bound$#ScalarTensor::STpi_group_bound       = "skip"
#$bound$#ScalarTensor::STpi_group_bound_speed = 1.0
#$bound$#ScalarTensor::STpi_group_bound_limit = 0.0
#$bound$#ScalarTensor::STpi_group_bound_scalar = 0.0

#$bound$#ScalarTensor::STphi_bound       = "skip"
#$bound$#ScalarTensor::STphi_bound_speed = 1.0
#$bound$#ScalarTensor::STphi_bound_limit = 0.0
#$bound$#ScalarTensor::STphi_bound_scalar = 0.0

#$bound$#ScalarTensor::STpi_bound       = "skip"
#$bound$#ScalarTensor::STpi_bound_speed = 1.0
#$bound$#ScalarTensor::STpi_bound_limit = 0.0
#$bound$#ScalarTensor::STpi_bound_scalar = 0.0

*/

