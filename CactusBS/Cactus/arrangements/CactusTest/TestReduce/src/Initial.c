#include "cctk.h"
#include "cctk_Arguments.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusTest_TestReduce_Initial_c)

void TestReduce_Initial(CCTK_ARGUMENTS);

void TestReduce_Initial(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS

  int index;
  int i,j,k;

  for(k=0; k<cctk_lsh[2]; k++)
  {
    for(j=0; j<cctk_lsh[1]; j++)
    {
      for(i=0; i<cctk_lsh[0]; i++)
      {
        index =  CCTK_GFINDEX3D(cctkGH,i,j,k);
        
        phi_gf3[index] = (float)((i+1)*(j+1)*(k+1));
        phi_igf3[index] = ((i+1)*(j+1)*(k+1));
      }
    }
  }

  for(j=0; j<cctk_lsh[1]; j++)
  {
    for(i=0; i<cctk_lsh[0]; i++)
    {
      index =  CCTK_GFINDEX2D(cctkGH,i,j);
      
      phi_gf2[index] = (float)((i+1)*(j+1));
      phi_igf2[index] = ((i+1)*(j+1));
    }
  }

  for(i=0; i<cctk_lsh[0]; i++)
  {
    index =  CCTK_GFINDEX1D(cctkGH,i);
    
    phi_gf1[index] = (float)((i+1));
    phi_igf1[index] = ((i+1));
  }

  *phi_scalar = 1.0;
  *phi_iscalar = 1;

  return;

}
