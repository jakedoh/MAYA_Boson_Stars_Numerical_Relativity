#include <stdio.h>

#include "cctk.h"
#include "cctk_Arguments.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusTest_TestReduce_Reduce_Maximum_c)

void TestReduceC_Maximum(CCTK_ARGUMENTS);

void TestReduceC_Maximum(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS

  int handle;
  int index;
  int ierr;
  CCTK_REAL reduction_value;
  CCTK_INT ireduction_value;

  CCTK_INFO("TESTING MAXIMUM REDUCTIONS FROM C");

  /* 3D GF */

  index = CCTK_VarIndex("TestReduce::phi_gf3");
  handle = CCTK_ReductionHandle("maximum");

  if (handle > 0)
  {
    ierr =  CCTK_Reduce (cctkGH, 
                         0, 
                         handle, 
                         1, 
                         CCTK_VARIABLE_REAL,
                         &reduction_value, 
                         1, 
                         index);

    if (ierr == 0 && CCTK_MyProc(cctkGH)==0)
    {
      printf("Maximum value of phi_gf3 is %f\n",reduction_value);
    }
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestReduce: Invalid reduction operator '%s'", 
                "maximum");
  }

  index = CCTK_VarIndex("TestReduce::phi_igf3");
  handle = CCTK_ReductionHandle("maximum");

  if (handle > 0)
  {
    ierr =  CCTK_Reduce (cctkGH, 
                         0, 
                         handle, 
                         1, 
                         CCTK_VARIABLE_INT,
                         &ireduction_value, 
                         1, 
                         index);

    if (ierr == 0 && CCTK_MyProc(cctkGH)==0)
    {
      printf("Maximum value of phi_igf3 is %d\n",(int)ireduction_value);
    }
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestReduce: Invalid reduction operator '%s'", 
                "maximum");
  }

  /* 2D GF */

  index = CCTK_VarIndex("TestReduce::phi_gf2");
  handle = CCTK_ReductionHandle("maximum");

  if (handle > 0)
  {
    ierr =  CCTK_Reduce (cctkGH, 
                         0, 
                         handle, 
                         1, 
                         CCTK_VARIABLE_REAL,
                         &reduction_value, 
                         1, 
                         index);

    if (ierr == 0 && CCTK_MyProc(cctkGH)==0)
    {
      printf("Maximum value of phi_gf2 is %f\n",reduction_value);
    }
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestReduce: Invalid reduction operator '%s'", 
                "maximum");
  }

  index = CCTK_VarIndex("TestReduce::phi_igf2");
  handle = CCTK_ReductionHandle("maximum");

  if (handle > 0)
  {
    ierr =  CCTK_Reduce (cctkGH, 
                         0, 
                         handle, 
                         1, 
                         CCTK_VARIABLE_INT,
                         &ireduction_value, 
                         1, 
                         index);

    if (ierr == 0 && CCTK_MyProc(cctkGH)==0)
    {
      printf("Maximum value of phi_igf2 is %d\n",(int)ireduction_value);
    }
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestReduce: Invalid reduction operator '%s'", 
                "maximum");
  }


  /* 1D GF */

  index = CCTK_VarIndex("TestReduce::phi_gf1");  
  handle = CCTK_ReductionHandle("maximum");

  if (handle > 0)
  {
    ierr =  CCTK_Reduce (cctkGH, 
                         0, 
                         handle, 
                         1, 
                         CCTK_VARIABLE_REAL,
                         &reduction_value, 
                         1, 
                         index);

    if (ierr == 0 && CCTK_MyProc(cctkGH)==0)
    {
      printf("Maximum value of phi_gf1 is %f\n",reduction_value);
    }
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestReduce: Invalid reduction operator '%s'", 
                "maximum");
  }

  index = CCTK_VarIndex("TestReduce::phi_igf1");  
  handle = CCTK_ReductionHandle("maximum");

  if (handle > 0)
  {
    ierr =  CCTK_Reduce (cctkGH, 
                         0, 
                         handle, 
                         1, 
                         CCTK_VARIABLE_INT,
                         &ireduction_value, 
                         1, 
                         index);

    if (ierr == 0 && CCTK_MyProc(cctkGH)==0)
    {
      printf("Maximum value of phi_igf1 is %d\n",(int)ireduction_value);
    }
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestReduce: Invalid reduction operator '%s'", 
                "maximum");
  }


  /* Grid Scalar */

  index = CCTK_VarIndex("TestReduce::phi_scalar");
  handle = CCTK_ReductionHandle("maximum");

  if (handle > 0)
  {
    ierr =  CCTK_Reduce (cctkGH, 
                         0, 
                         handle, 
                         1, 
                         CCTK_VARIABLE_REAL,
                         &reduction_value, 
                         1, 
                         index);

    if (ierr == 0 && CCTK_MyProc(cctkGH)==0)
    {
      printf("Maximum value of phi_scalar is %f\n",reduction_value);
    }
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestReduce: Invalid reduction operator '%s'", 
                "maximum");
  }

  index = CCTK_VarIndex("TestReduce::phi_iscalar");
  handle = CCTK_ReductionHandle("maximum");

  if (handle > 0)
  {
    ierr =  CCTK_Reduce (cctkGH, 
                         0, 
                         handle, 
                         1, 
                         CCTK_VARIABLE_INT,
                         &ireduction_value, 
                         1, 
                         index);

    if (ierr == 0 && CCTK_MyProc(cctkGH)==0)
    {
      printf("Maximum value of phi_iscalar is %d\n",(int)ireduction_value);
    }
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "TestReduce: Invalid reduction operator '%s'", 
                "maximum");
  }

  return;

}
