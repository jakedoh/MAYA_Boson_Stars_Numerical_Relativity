 /*@@
   @file      TestGlobalReduce.c
   @date      Tues Sep  7 16:00:29 2004
   @author    Yaakoub El Khamra
   @desc
              Test the global reduction implementation
   @enddesc
   @version   $Header$
 @@*/

#include <stdio.h>

#include "cctk.h"
#include "util_Table.h"

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusTest_TestGlobalReduce_TestGlobalReduce_c);


/********************************************************************
 *********************  Scheduled Routine Prototypes  ***************
 ********************************************************************/


void TestGlobalReduceC_Maximum(CCTK_ARGUMENTS);

void TestGlobalReduceC_Minimum(CCTK_ARGUMENTS);

void TestGlobalReduceC_Sum(CCTK_ARGUMENTS);

void TestGlobalReduceC_Avg(CCTK_ARGUMENTS);

void TestGlobalReduceC_GF(CCTK_ARGUMENTS);

void TestGlobalReduceC_Weighted_Maximum(CCTK_ARGUMENTS);

void TestGlobalReduceC_Weighted_Minimum(CCTK_ARGUMENTS);

void TestGlobalReduceC_Weighted_Sum(CCTK_ARGUMENTS);

void TestGlobalReduceC_Weighted_Avg(CCTK_ARGUMENTS);


/********************************************************************
 *********************  Scheduled Routines  *************************
 ********************************************************************/

/*@@
   @routine TestGlobalReduceC_Maximum
   @author     Yaakoub Y El Khamra
   @date       Tues Sep  7 16:00:29 2004
   @desc
               Tests global maximum reduction
   @enddesc
   @calls      CCTK_LocalArrayReductionHandle
               CCTK_ReduceGridArrays
               CCTK_VarIndex

   @returntype int
   @returndesc
                  return error code from CCTK_ReduceGridArrays
                  else return 0
   @endreturndesc
 @@*/
void TestGlobalReduceC_Maximum(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int handle, ierr;

  /* input array void * pointers */
  CCTK_INT input_array[2];

  /* array containing types of data */
  CCTK_INT input_array_type_codes[2];

  /* output data */
  void*  reduction_3D_value[2];
  
  /* print what we are testing for */
  CCTK_Info("TestGlobalReduce", "Testing for maximum value");

  /* make sure all pointers point to valid memory locations */  
  reduction_3D_value[0] = real_max_3D;
  reduction_3D_value[1] = int_max_3D; 

  input_array_type_codes[0]  = CCTK_VARIABLE_REAL;
  input_array_type_codes[1]  = CCTK_VARIABLE_INT;

  /* 3D Array */
  handle = CCTK_LocalArrayReductionHandle("maximum");
  input_array[0] = CCTK_VarIndex("TestGlobalReduce::grid_real");
  input_array[1] = CCTK_VarIndex("TestGlobalReduce::grid_int");

if (handle >= 0)
  {
    ierr =  CCTK_ReduceGridArrays (cctkGH, dest_proc, handle, -1, 2, input_array , 2, input_array_type_codes, reduction_3D_value);
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestGlobalReduce: Invalid reduction operator '%s' with handle '%d'", 
                "maximum", handle);
  }
}

/*@@
   @routine TestGlobalReduceC_Minimum
   @author     Yaakoub Y El Khamra
   @date       Tues Sep  7 16:00:29 2004
   @desc
               Tests global minimum reduction
   @enddesc
   @calls      CCTK_LocalArrayReductionHandle
               CCTK_ReduceGridArrays
               CCTK_VarIndex

   @returntype int
   @returndesc
                  return error code from CCTK_ReduceGridArrays
                  else return 0
   @endreturndesc
 @@*/
void TestGlobalReduceC_Minimum(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int handle, ierr;

  /* input array void * pointers */
  CCTK_INT input_array[2];

  /* array containing types of data */
  CCTK_INT input_array_type_codes[2];

  /* output data */
  void*  reduction_3D_value[2];
  
  /* print what we are testing for */
  CCTK_Info("TestGlobalReduce", "Testing for minimum value");
  
  /* make sure all pointers point to valid memory locations */  
  reduction_3D_value[0] = real_min_3D;
  reduction_3D_value[1] = int_min_3D; 

  input_array_type_codes[0]  = CCTK_VARIABLE_REAL;
  input_array_type_codes[1]  = CCTK_VARIABLE_INT;

  /* 3D Array */
  handle = CCTK_LocalArrayReductionHandle("minimum");
  input_array[0] = CCTK_VarIndex("TestGlobalReduce::grid_real");
  input_array[1] = CCTK_VarIndex("TestGlobalReduce::grid_int");

if (handle >= 0)
  {
    ierr =  CCTK_ReduceGridArrays (cctkGH, dest_proc, handle, -1, 2, input_array , 2, input_array_type_codes, reduction_3D_value);
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestGlobalReduce: Invalid reduction operator '%s' with handle '%d'", 
                "minimum", handle);
  }
}

/*@@
   @routine TestGlobalReduceC_Sum
   @author     Yaakoub Y El Khamra
   @date       Tues Sep  7 16:00:29 2004
   @desc
               Tests global sum reduction
   @enddesc
   @calls      CCTK_LocalArrayReductionHandle
               CCTK_ReduceGridArrays
               CCTK_VarIndex

   @returntype int
   @returndesc
                  return error code from CCTK_ReduceGridArrays
                  else return 0
   @endreturndesc
 @@*/
void TestGlobalReduceC_Sum(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int handle, ierr;

  /* input array void * pointers */
  CCTK_INT input_array[2];

  /* array containing types of data */
  CCTK_INT input_array_type_codes[2];

  /* output data */
  void*  reduction_3D_value[2];
  
  /* print what we are testing for */
  CCTK_Info("TestGlobalReduce", "Testing for sum");
  
  /* make sure all pointers point to valid memory locations */  
  reduction_3D_value[0] = real_sum_3D;
  reduction_3D_value[1] = int_sum_3D; 

  input_array_type_codes[0]  = CCTK_VARIABLE_REAL;
  input_array_type_codes[1]  = CCTK_VARIABLE_INT;

  /* 3D Array */
  handle = CCTK_LocalArrayReductionHandle("sum");
  input_array[0] = CCTK_VarIndex("TestGlobalReduce::grid_real");
  input_array[1] = CCTK_VarIndex("TestGlobalReduce::grid_int");

  if (handle >= 0)
  {
    ierr =  CCTK_ReduceGridArrays (cctkGH, dest_proc, handle, -1, 2, input_array , 2, input_array_type_codes, reduction_3D_value);
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestLocalReduce: Invalid reduction operator '%s' with handle '%d'", 
                "sum", handle);
  }
}

/*@@
   @routine TestGlobalReduceC_Avg
   @author     Yaakoub Y El Khamra
   @date       Tues Sep  7 16:00:29 2004
   @desc
               Tests global average (mean) reduction
   @enddesc
   @calls      CCTK_LocalArrayReductionHandle
               CCTK_ReduceGridArrays
               CCTK_VarIndex

   @returntype int
   @returndesc
                  return error code from CCTK_ReduceGridArrays
                  else return 0
   @endreturndesc
 @@*/
void TestGlobalReduceC_Avg(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int handle, ierr;

  /* input array void * pointers */
  CCTK_INT input_array[2];

  /* array containing types of data */
  CCTK_INT input_array_type_codes[2];

  /* output data */
  void*  reduction_3D_value[2];
  
  /* print what we are testing for */
  CCTK_Info("TestGlobalReduce", "Testing for Average ");

  /* make sure all pointers point to valid memory locations */  
  reduction_3D_value[0] = real_avg_3D;
  reduction_3D_value[1] = int_avg_3D; 

  input_array_type_codes[0]  = CCTK_VARIABLE_REAL;
  input_array_type_codes[1]  = CCTK_VARIABLE_INT;

  /* 3D Array */
  handle = CCTK_LocalArrayReductionHandle("mean");
  input_array[0] = CCTK_VarIndex("TestGlobalReduce::grid_real");
  input_array[1] = CCTK_VarIndex("TestGlobalReduce::grid_int");

  if (handle >= 0)
  {
    ierr =  CCTK_ReduceGridArrays (cctkGH, dest_proc, handle, -1, 2, input_array , 2, input_array_type_codes, reduction_3D_value);
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestLocalReduce: Invalid reduction operator '%s' with handle '%d'", 
                "sum", handle);
  }
}


/*@@
   @routine TestGlobalReduceC_GF
   @author     Yaakoub Y El Khamra
   @date       Tues Sep  7 16:00:29 2004
   @desc
               Tests grid function reduction (sum and max)
   @enddesc
   @calls      CCTK_LocalArrayReductionHandle
               CCTK_ReduceGridArrays
               CCTK_VarIndex

   @returntype int
   @returndesc
                  return error code from CCTK_ReduceGridArrays
                  else return 0
   @endreturndesc
 @@*/
void TestGlobalReduceC_GF(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int handle, ierr;

  /* input array void * pointers */
  CCTK_INT input_array[1];

  /* array containing types of data */
  CCTK_INT input_array_type_codes[1];

  /* output data */
  void*  reduction_3D_value[1];

  input_array_type_codes[0]  = CCTK_VARIABLE_REAL;
  input_array[0] = CCTK_VarIndex("TestGlobalReduce::real_test_gf");

  CCTK_Info("TestGlobalReduce", "Testing for maximum value of GF");
  /* make sure all pointers point to valid memory locations */  
  reduction_3D_value[0] = gf_max;

  /* Max reduction */
  handle = CCTK_LocalArrayReductionHandle("maximum");

  if (handle >= 0)
  {
    ierr =  CCTK_ReduceGridArrays (cctkGH, dest_proc, handle, -1, 1, input_array , 1, input_array_type_codes, reduction_3D_value);
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestLocalReduce: Invalid reduction operator '%s' with handle '%d'", 
                "max", handle);
  }

  CCTK_Info("TestGlobalReduce", "Testing for sum value of GF");
  /* make sure all pointers point to valid memory locations */  
  reduction_3D_value[0] = gf_sum;

  /* Sum reduction */

  input_array[0] = CCTK_VarIndex("TestGlobalReduce::real_sum_gf");
  handle = CCTK_LocalArrayReductionHandle("sum");

  if (handle >= 0)
  {
    ierr =  CCTK_ReduceGridArrays (cctkGH, dest_proc, handle, -1, 1, input_array , 1, input_array_type_codes, reduction_3D_value);
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestLocalReduce: Invalid reduction operator '%s' with handle '%d'", 
                "sum", handle);
  }
}

/*@@
   @routine TestGlobalReduceC_Norms
   @author     Yaakoub Y El Khamra
   @date       
   @desc
               Tests norm 1-4 reductions
   @enddesc
   @calls      CCTK_LocalArrayReductionHandle
               CCTK_ReduceGridArrays
               CCTK_VarIndex

   @returntype int
   @returndesc
                  return error code from CCTK_ReduceGridArrays
                  else return 0
   @endreturndesc
 @@*/
void TestGlobalReduceC_Norms(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int handle, ierr;

  /* input array void * pointers */
  CCTK_INT input_array[1];

  /* array containing types of data */
  CCTK_INT input_array_type_codes[1];

  /* output data */
  void*  reduction_3D_value[1];

  input_array_type_codes[0]  = CCTK_VARIABLE_REAL;
  input_array[0] = CCTK_VarIndex("TestGlobalReduce::grid_real");
  
  CCTK_Info("TestGlobalReduce", "Testing norm1, norm2, norm3, norm4 and norminf value of GF");

  /* Nrom1 reduction */
  handle = CCTK_LocalArrayReductionHandle("norm1");
  reduction_3D_value[0] = gf_norm1;
  
  if (handle >= 0)
  {
    ierr =  CCTK_ReduceGridArrays (cctkGH, dest_proc, handle, -1, 1, input_array , 1, input_array_type_codes, reduction_3D_value);
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestLocalReduce: Invalid reduction operator '%s' with handle '%d'", 
                "norm1", handle);
  }

  /* Nrom2 reduction */
  handle = CCTK_LocalArrayReductionHandle("norm2");
  reduction_3D_value[0] = gf_norm2;
  
  if (handle >= 0)
  {
    ierr =  CCTK_ReduceGridArrays (cctkGH, dest_proc, handle, -1, 1, input_array , 1, input_array_type_codes, reduction_3D_value);
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestLocalReduce: Invalid reduction operator '%s' with handle '%d'", 
                "norm2", handle);
  }
  
  /* Nrom3 reduction */
  handle = CCTK_LocalArrayReductionHandle("norm3");
  reduction_3D_value[0] = gf_norm3;
  
  if (handle >= 0)
  {
    ierr =  CCTK_ReduceGridArrays (cctkGH, dest_proc, handle, -1, 1, input_array , 1, input_array_type_codes, reduction_3D_value);
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestLocalReduce: Invalid reduction operator '%s' with handle '%d'", 
                "norm3", handle);
  }

  /* Nrom4 reduction */
  handle = CCTK_LocalArrayReductionHandle("norm4");
  reduction_3D_value[0] = gf_norm4;
  
  if (handle >= 0)
  {
    ierr =  CCTK_ReduceGridArrays (cctkGH, dest_proc, handle, -1, 1, input_array , 1, input_array_type_codes, reduction_3D_value);
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestLocalReduce: Invalid reduction operator '%s' with handle '%d'", 
                "norm4", handle);
  }
}


/*@@
   @routine TestGlobalReduceC_Scalar
   @author     Yaakoub Y El Khamra
   @date       
   @desc
               Tests scalar reduction
   @enddesc
   @calls      CCTK_LocalArrayReductionHandle
               CCTK_ReduceGridArrays
               CCTK_VarIndex

   @returntype int
   @returndesc
                  return error code from CCTK_ReduceGridArrays
                  else return 0
   @endreturndesc
 @@*/
void TestGlobalReduceC_Scalar(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int handle, ierr;
  int nprocs = 1;
  
  /* input array void * pointers */
  CCTK_INT input_array[1];

  /* array containing types of data */
  CCTK_INT input_array_type_codes[1];

  /* output data */
  void*  value[1];

  input_array_type_codes[0]  = CCTK_VARIABLE_REAL;
  input_array[0] = CCTK_VarIndex("TestGlobalReduce::num_to_reduce");

  CCTK_Info("TestGlobalReduce", "Testing for maximum scalar reduction");
  /* make sure all pointers point to valid memory locations */  
  value[0] = max_value;

  /* Max reduction */
  handle = CCTK_LocalArrayReductionHandle("maximum");

  if (handle >= 0)
  {
    ierr =  CCTK_ReduceGridArrays (cctkGH, dest_proc, handle, -1, 1, input_array , 1, input_array_type_codes, value);
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestLocalReduce: Invalid reduction operator '%s' with handle '%d'", 
                "max", handle);
  }
  
  CCTK_Info("TestGlobalReduce", "Testing for sum scalar reduction");
  /* make sure all pointers point to valid memory locations */  
  value[0] = sum_value;
  nprocs = CCTK_nProcs(cctkGH);
  
  /* Max reduction */
  handle = CCTK_LocalArrayReductionHandle("sum");

  if (handle >= 0)
  {
    ierr =  CCTK_ReduceGridArrays (cctkGH, dest_proc, handle, -1, 1, input_array , 1, input_array_type_codes, value);
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestLocalReduce: Invalid reduction operator '%s' with handle '%d'", 
                "sum", handle);
  }
  
  *sum_value = *sum_value / nprocs;
}

/*@@
   @routine TestGlobalReduceC_Weighted_Maximum
   @author     Yaakoub Y El Khamra
   @date       Tues Sep  7 16:00:29 2004
   @desc
               Tests global maximum reduction
   @enddesc
   @calls      CCTK_LocalArrayReductionHandle
               CCTK_ReduceGridArrays
               CCTK_VarIndex

   @returntype int
   @returndesc
                  return error code from CCTK_ReduceGridArrays
                  else return 0
   @endreturndesc
 @@*/
void TestGlobalReduceC_Weighted_Maximum(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int handle, ierr, weight_index, par_handle;

  /* input array void * pointers */
  CCTK_INT input_array[2];

  /* array containing types of data */
  CCTK_INT input_array_type_codes[2];

  /* output data */
  void*  reduction_3D_value[2];
  
  /* print what we are testing for */
  CCTK_Info("TestGlobalReduce", "Testing for maximum value");

  /* Create parameter table and store weight */
  par_handle = Util_TableCreate (UTIL_TABLE_FLAGS_DEFAULT);
  weight_index = CCTK_VarIndex("TestGlobalReduce::weight");
  ierr = Util_TableSetInt(par_handle, 1, "weight_on");
  ierr = Util_TableSetInt(par_handle, weight_index, "weight_variable_index");

  /* make sure all pointers point to valid memory locations */  
  reduction_3D_value[0] = weighted_real_max;
  reduction_3D_value[1] = weighted_int_max; 

  input_array_type_codes[0]  = CCTK_VARIABLE_REAL;
  input_array_type_codes[1]  = CCTK_VARIABLE_INT;

  /* 3D Array */
  handle = CCTK_LocalArrayReductionHandle("maximum");
  input_array[0] = CCTK_VarIndex("TestGlobalReduce::grid_real");
  input_array[1] = CCTK_VarIndex("TestGlobalReduce::grid_int");

if (handle >= 0)
  {
    ierr =  CCTK_ReduceGridArrays (cctkGH, dest_proc, handle, weight_index, 2, input_array , 2, input_array_type_codes, reduction_3D_value);
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestGlobalReduce: Invalid reduction operator '%s' with handle '%d'", 
                "maximum", handle);
  }
}

/*@@
   @routine TestGlobalReduceC_Weighted_Minimum
   @author     Yaakoub Y El Khamra
   @date       Tues Sep  7 16:00:29 2004
   @desc
               Tests global minimum reduction
   @enddesc
   @calls      CCTK_LocalArrayReductionHandle
               CCTK_ReduceGridArrays
               CCTK_VarIndex

   @returntype int
   @returndesc
                  return error code from CCTK_ReduceGridArrays
                  else return 0
   @endreturndesc
 @@*/
void TestGlobalReduceC_Weighted_Minimum(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int handle, ierr, weight_index, par_handle;

  /* input array void * pointers */
  CCTK_INT input_array[2];

  /* array containing types of data */
  CCTK_INT input_array_type_codes[2];

  /* output data */
  void*  reduction_3D_value[2];
  
  /* print what we are testing for */
  CCTK_Info("TestGlobalReduce", "Testing for minimum value");

  /* Create parameter table and store weight */
  par_handle = Util_TableCreate (UTIL_TABLE_FLAGS_DEFAULT);
  weight_index = CCTK_VarIndex("TestGlobalReduce::weight");
  ierr = Util_TableSetInt(par_handle, 1, "weight_on");
  ierr = Util_TableSetInt(par_handle, weight_index, "weight_variable_index");

  /* make sure all pointers point to valid memory locations */  
  reduction_3D_value[0] = weighted_real_min;
  reduction_3D_value[1] = weighted_int_min; 

  input_array_type_codes[0]  = CCTK_VARIABLE_REAL;
  input_array_type_codes[1]  = CCTK_VARIABLE_INT;

  /* 3D Array */
  handle = CCTK_LocalArrayReductionHandle("minimum");
  input_array[0] = CCTK_VarIndex("TestGlobalReduce::grid_real");
  input_array[1] = CCTK_VarIndex("TestGlobalReduce::grid_int");

if (handle >= 0)
  {
    ierr =  CCTK_ReduceGridArrays (cctkGH, dest_proc, handle, par_handle, 2, input_array , 2, input_array_type_codes, reduction_3D_value);
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestGlobalReduce: Invalid reduction operator '%s' with handle '%d'", 
                "minimum", handle);
  }
}

/*@@
   @routine TestGlobalReduceC_Weighted_Sum
   @author     Yaakoub Y El Khamra
   @date       Tues Sep  7 16:00:29 2004
   @desc
               Tests global sum reduction
   @enddesc
   @calls      CCTK_LocalArrayReductionHandle
               CCTK_ReduceGridArrays
               CCTK_VarIndex

   @returntype int
   @returndesc
                  return error code from CCTK_ReduceGridArrays
                  else return 0
   @endreturndesc
 @@*/
void TestGlobalReduceC_Weighted_Sum(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int handle, ierr, weight_index, par_handle;

  /* input array void * pointers */
  CCTK_INT input_array[2];

  /* array containing types of data */
  CCTK_INT input_array_type_codes[2];

  /* output data */
  void*  reduction_3D_value[2];
  
  /* print what we are testing for */
  CCTK_Info("TestGlobalReduce", "Testing for sum");
 
  /* Create parameter table and store weight */
  par_handle = Util_TableCreate (UTIL_TABLE_FLAGS_DEFAULT);
  weight_index = CCTK_VarIndex("TestGlobalReduce::weight");
  ierr = Util_TableSetInt(par_handle, 1, "weight_on");
  ierr = Util_TableSetInt(par_handle, weight_index, "weight_variable_index");

  /* make sure all pointers point to valid memory locations */  
  reduction_3D_value[0] = weighted_real_sum;
  reduction_3D_value[1] = weighted_int_sum; 

  input_array_type_codes[0]  = CCTK_VARIABLE_REAL;
  input_array_type_codes[1]  = CCTK_VARIABLE_INT;

  /* 3D Array */
  handle = CCTK_LocalArrayReductionHandle("sum");
  input_array[0] = CCTK_VarIndex("TestGlobalReduce::grid_real");
  input_array[1] = CCTK_VarIndex("TestGlobalReduce::grid_int");

  if (handle >= 0)
  {
    ierr =  CCTK_ReduceGridArrays (cctkGH, dest_proc, handle, par_handle, 2, input_array , 2, input_array_type_codes, reduction_3D_value);
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestLocalReduce: Invalid reduction operator '%s' with handle '%d'", 
                "sum", handle);
  }
}

/*@@
   @routine TestGlobalReduceC_Weighted_Avg
   @author     Yaakoub Y El Khamra
   @date       Tues Sep  7 16:00:29 2004
   @desc
               Tests global average (mean) reduction
   @enddesc
   @calls      CCTK_LocalArrayReductionHandle
               CCTK_ReduceGridArrays
               CCTK_VarIndex

   @returntype int
   @returndesc
                  return error code from CCTK_ReduceGridArrays
                  else return 0
   @endreturndesc
 @@*/
void TestGlobalReduceC_Weighted_Avg(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int handle, ierr, weight_index, par_handle;

  /* input array void * pointers */
  CCTK_INT input_array[2];

  /* array containing types of data */
  CCTK_INT input_array_type_codes[2];

  /* output data */
  void*  reduction_3D_value[2];
  
  /* print what we are testing for */
  CCTK_Info("TestGlobalReduce", "Testing for Average ");

  /* Create parameter table and store weight */
  par_handle = Util_TableCreate (UTIL_TABLE_FLAGS_DEFAULT);
  weight_index = CCTK_VarIndex("TestGlobalReduce::weight");
  ierr = Util_TableSetInt(par_handle, 1, "weight_on");
  ierr = Util_TableSetInt(par_handle, weight_index, "weight_variable_index");

  /* make sure all pointers point to valid memory locations */  
  reduction_3D_value[0] = weighted_real_avg;
  reduction_3D_value[1] = weighted_int_avg;

  input_array_type_codes[0]  = CCTK_VARIABLE_REAL;
  input_array_type_codes[1]  = CCTK_VARIABLE_INT;

  /* 3D Array */
  handle = CCTK_LocalArrayReductionHandle("mean");
  input_array[0] = CCTK_VarIndex("TestGlobalReduce::grid_real");
  input_array[1] = CCTK_VarIndex("TestGlobalReduce::grid_int");

  if (handle >= 0)
  {
    ierr =  CCTK_ReduceGridArrays (cctkGH, dest_proc, handle, par_handle, 2, input_array , 2, input_array_type_codes, reduction_3D_value);
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestLocalReduce: Invalid reduction operator '%s' with handle '%d'", 
                "sum", handle);
  }
}
