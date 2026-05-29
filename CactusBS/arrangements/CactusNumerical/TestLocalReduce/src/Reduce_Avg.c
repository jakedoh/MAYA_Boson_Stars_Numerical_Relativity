 /*@@
   @file      Reduce_Avg.c
   @date      Tue Aug 31 12:48:07 2004
   @author    Yaakoub Y El-Khamra
   @desc 
   Test the local reduction interface and operators.
   @enddesc 
   @version $Header$
 @@*/

#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusTest_TestLocalReduce_Reduce_avg_c);

/********************************************************************
 *********************  Scheduled Routine Prototypes  ***************
 ********************************************************************/

void TestLocalReduceC_Avg(CCTK_ARGUMENTS);

/********************************************************************
 *********************  Scheduled Routines  *************************
 ********************************************************************/

/*@@
   @routine TestLocalReduceC_Avg
   @author     Yaakoub Y El Khamra
   @date       Tues Sep  7 16:00:29 2004
   @desc
               Tests local average (mean) reduction
   @enddesc
   @calls      CCTK_LocalArrayReductionHandle
               CCTK_ReduceLocalArrays

   @returntype void
 @@*/
void TestLocalReduceC_Avg(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int handle,ierr;
  int i,j,k, index1, index2;

  /* input array void * pointers */
  const void* input_array[2];

  /* arrays containing data */
  CCTK_REAL * real_3D_array;
  CCTK_INT  * int_3D_array;

  CCTK_REAL * real_2D_array;
  CCTK_INT  * int_2D_array;
  
  CCTK_REAL * real_1D_array;
  CCTK_INT  * int_1D_array;

  /* arrays containing dimensions of data */
  CCTK_INT input_array_dims3D[3];
  CCTK_INT input_array_dims2D[2];
  CCTK_INT input_array_dims1D[1];

  /* array containing types of data */
  CCTK_INT input_array_type_codes[2];

  /* output data */
  void*  reduction_3D_value[2];
  void*  reduction_2D_value[2];
  void*  reduction_1D_value[2];
  
  /* allocate memory to arrays */
  real_3D_array = (CCTK_REAL *) malloc (array_nx*array_ny*array_nz* sizeof(CCTK_REAL));
  int_3D_array = (CCTK_INT *) malloc (array_nx*array_ny*array_nz* sizeof(CCTK_REAL));

  real_2D_array = (CCTK_REAL *) malloc (array_nx*array_ny* sizeof(CCTK_REAL));
  int_2D_array = (CCTK_INT *) malloc (array_nx*array_ny* sizeof(CCTK_REAL));

  real_1D_array = (CCTK_REAL *) malloc (array_nx* sizeof(CCTK_REAL));
  int_1D_array = (CCTK_INT *) malloc (array_nx* sizeof(CCTK_REAL));

  /* setting dimensions */
  input_array_dims3D[0] = array_nx;
  input_array_dims3D[1] = array_ny;
  input_array_dims3D[2] = array_nz;

  input_array_dims2D[0] = array_nx;
  input_array_dims2D[1] = array_ny;

  input_array_dims1D[0] = array_nx;

  /* filling in values in 3,2 and 1D arrays */
  for (i=0;i<array_nx;i++)
  {
    real_1D_array[i] = (CCTK_REAL)(i+1);
    int_1D_array[i]  = (CCTK_INT)(i+2);
    for(j=0;j<array_ny;j++)
    {
      index1 = i + array_nx*j;
      real_2D_array[index1] = (CCTK_REAL)((i+1)*(j+1));
      int_2D_array[index1]  = (CCTK_INT)((i+2)*(j+2));
      for(k=0;k<array_nz;k++)
      {
        index2 = i + array_nx*(j+array_ny*k);
        real_3D_array[index2] = (CCTK_REAL)((i+1)*(j+1)*(k+1));
        int_3D_array[index2]  = (CCTK_INT)((i+2)*(j+2)*(k+2));
      }
    }
  }
      
  /* make sure all pointers point to valid memory locations */  
  reduction_3D_value[0] = real_avg_3D;
  reduction_3D_value[1] = int_avg_3D;
  reduction_2D_value[0] = real_avg_2D;
  reduction_2D_value[1] = int_avg_2D;
  reduction_1D_value[0] = real_avg_1D;
  reduction_1D_value[1] = int_avg_1D;
  

  input_array_type_codes[0]  = CCTK_VARIABLE_REAL;
  input_array_type_codes[1]  = CCTK_VARIABLE_INT;

  /* 3D Array */
  handle = CCTK_LocalArrayReductionHandle("mean");
  input_array[0] = real_3D_array;
  input_array[1] = int_3D_array;

  if (handle >= 0)
  {
    ierr =  CCTK_ReduceLocalArrays (3, handle, 0, 2, input_array_dims3D, input_array_type_codes,input_array , 2, input_array_type_codes, reduction_3D_value);
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestLocalReduce: Invalid reduction operator '%s' with handle '%d'", 
                "sum", handle);
  }

  /* 2D Array */
  handle = CCTK_LocalArrayReductionHandle("mean");
  input_array[0] = real_2D_array;
  input_array[1] = int_2D_array;

  if (handle >= 0)
  {
    ierr =  CCTK_ReduceLocalArrays (2, handle, 0, 2, input_array_dims3D, input_array_type_codes,input_array , 2, input_array_type_codes, reduction_2D_value);

  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestLocalReduce: Invalid reduction operator '%s' with handle '%d'", 
                "sum", handle);
  }

  /* 1D Array */
  handle = CCTK_LocalArrayReductionHandle("mean");
  input_array[0] = real_1D_array;
  input_array[1] = int_1D_array;

  if (handle >= 0)
  {
    ierr =  CCTK_ReduceLocalArrays (1, handle, 0, 2, input_array_dims3D, input_array_type_codes,input_array , 2, input_array_type_codes, reduction_1D_value);

  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestLocalReduce: Invalid reduction operator '%s' with handle '%d'", 
                "sum", handle);
  }
  
  /* freeing memory */
  free (real_3D_array );
  free (int_3D_array );
  free (real_2D_array );
  free (int_2D_array );
  free (real_1D_array );
  free (int_1D_array );

  return;

}
