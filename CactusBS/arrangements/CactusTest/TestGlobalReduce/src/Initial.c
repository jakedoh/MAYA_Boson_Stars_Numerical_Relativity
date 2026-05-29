 /*@@
   @file      Initial.c
   @date      Tues Sep  7 16:00:29 2004
   @author    Yaakoub El Khamra
   @desc
              Setup the grid arrays we want to test reduction with
   @enddesc
   @version   $Id$
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <math.h>

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusTest_TestGlobalReduce_Initial_c)

#define SQR(x) ((x) * (x))

/********************************************************************
 *********************  Scheduled Routine Prototypes  ***************
 ********************************************************************/
void TestGlobalReduce_Initial(CCTK_ARGUMENTS);

/********************************************************************
 *********************  Scheduled Routines  *************************
 ********************************************************************/

/*@@
   @routine TestGlobalReduce_Initial
   @author     Yaakoub Y El Khamra
   @date       Tues Sep  7 16:00:29 2004
   @desc
               Sets up the grid arrays and fills them with values
   @enddesc
   @calls      CCTK_GroupubndVI
               CCTK_GrouplbndVI

   @returntype void
 @@*/
void TestGlobalReduce_Initial(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int sum_index;
  int i,j,k;
  int ierr;
  int lsh[3];
  int min[3]; 

  CCTK_REAL dx,dy,dz;
  CCTK_REAL X,Y,Z;
  CCTK_REAL R;

  /* Set up shorthands */
  dx = CCTK_DELTA_SPACE(0);
  dy = CCTK_DELTA_SPACE(1);
  dz = CCTK_DELTA_SPACE(2);
  
  ierr = CCTK_GrouplshVN(cctkGH, 3, lsh, "TestGlobalReduce::grid_real");
  if (ierr < 0)
    CCTK_WARN(0, "Error reading group lsh");
  ierr = CCTK_GrouplbndVN(cctkGH, 3, min, "TestGlobalReduce::grid_real");
  if (ierr < 0)
    CCTK_WARN(0, "Error reading group lbnd");

  for(k=0; k<lsh[2]; k++)
  {
    for(j=0; j<lsh[1]; j++)
    {
      for(i=0; i<lsh[0]; i++)
      {      
        sum_index = i+ lsh[0]*(j + lsh[1]*k);
        grid_real[sum_index] = ((i+min[0]+1)*(j+min[1]+1)*(k+min[2]+1));
        grid_int [sum_index] = ((i+min[0]+2)*(j+min[1]+2)*(k+min[2]+2));
        weight   [sum_index] = uniform_weight_value;
      }
    }
  }

  for(k=0; k<cctk_lsh[2]; k++)
  {
    for(j=0; j<cctk_lsh[1]; j++)
    {
      for(i=0; i<cctk_lsh[0]; i++)
      {      
        sum_index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        X = x[sum_index];
	Y = y[sum_index];
	Z = z[sum_index];
        R = sqrt( X*X + Y*Y + Z*Z );
        real_test_gf[sum_index] = (amplitude*exp( -SQR( (R - radius) / sigma ) ));
        real_sum_gf[sum_index] = 0.0;
      }
    }
  }

  for(k=0; k<cctk_lsh[2]; k++)
  {
    for(j=0; j<cctk_lsh[1]; j++)
    {
      for(i=0; i<cctk_lsh[0]; i++)
      {      
        sum_index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        X = x[sum_index];
	Y = y[sum_index];
	Z = z[sum_index];

        if (fabs(Z-1.0)<1e-14)
          continue;
        if (fabs(Y-1.0)<1e-14)
          continue;
        if (fabs(X-1.0)<1e-14)
          continue;
        
        real_sum_gf[sum_index] = (dx*dy*dz*(X+Y+Z+dx/2+dy/2+dz/2));
      }
    }
  }
  
  /* set the value of the scalar to reduce to an arbitrary value */
  num_to_reduce[0]=3.14159;
  
  return;
}
