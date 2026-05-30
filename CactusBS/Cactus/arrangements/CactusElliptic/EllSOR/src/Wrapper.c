 /*@@
   @file      Wrapper.c
   @date      Tue Aug 24 12:50:07 1999
   @author    Gerd Lanfermann
   @desc 
   The C wrapper, which calles the core Fortran routine, which 
   performs the actual solve.
   We cannot derive the pointers to the GF data from the indices in 
   Fortran. So we do this here in C and then pass the everything 
   over to the Fortran routine.

   This wrapper is registers with the Elliptic solver registry 
   (not the Fortran file) , as coded up in ./CactusElliptic/EllBase
   @enddesc 
 @@*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_FortranString.h"

#include "EllBase.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusElliptic_EllSOR_Wrapper_c)

int SORFlat3D(cGH *GH, int FieldIndex, int MIndex, int NIndex,
              CCTK_REAL *AbsTol, CCTK_REAL *RelTol);
int SORConfMetric3D(cGH *GH, int *MetricPsiI, int conformal,
                    int FieldIndex, int MIndex, int NIndex,
                    CCTK_REAL *AbsTol, CCTK_REAL *RelTol);
int SORConfMetric(cGH *GH, 
                  int *MetricPsiI, 
                  int FieldIndex, 
                  int MIndex, 
                  int NIndex, 
                  CCTK_REAL *AbsTol,
                  CCTK_REAL *RelTol);
int SORFlat(cGH *GH, 
            int FieldIndex, 
            int MIndex, 
            int NIndex, 
            CCTK_REAL *AbsTol, 
            CCTK_REAL *RelTol);

/*@@
   @routine    SORConfMetric
   @date       Tue Sep 26 11:31:42 2000
   @author     Gerd Lanfermann
   @desc 
     elliptic solver wrapper which provides a interface to the 
     different n-dimensional SOR solvers for the conformal metric, 
     of which only 3D is coded currently.

     This wrapper is registered and if it is being called with 
     a n-dim. grid function, it goes of and picks the correct solver.

     We pass in the arguments that are neccessary for this class of elliptic eq. 
     this solver is intended to solve. See ./CactusElliptic/EllBase/src/ for the
     classes of elliptic eq. 
   @enddesc 
   @calls      
   @calledby   
   @history 
 
   @endhistory 

@@*/

int SORConfMetric(cGH *GH, 
                  int *MetricPsiI, 
                   int FieldIndex, 
                   int MIndex, 
                   int NIndex, 
                   CCTK_REAL *AbsTol,
                   CCTK_REAL *RelTol) 
{
  int retval = ELL_NOSOLVER;

  switch (CCTK_GroupDimFromVarI(FieldIndex))
  {
    case 1:  
      CCTK_WARN(0,"SORConfMetric: No 1D SOR solver implemented");
      break;
    case 2:
      CCTK_WARN(0,"SORConfMetric: No 2D SOR solver implemented");
      break;
    case 3:
      retval = SORConfMetric3D(GH, MetricPsiI, 1, 
                               FieldIndex, MIndex, NIndex,
                               AbsTol, RelTol);
      break;
    default:
      CCTK_WARN(1,"SORConfMetric: Solver only implemented for 3D");
      break;
  }

  return retval;
}
 
int SORMetric(cGH *GH, 
              int *MetricI, 
              int FieldIndex, 
              int MIndex, 
              int NIndex, 
              CCTK_REAL *AbsTol,
              CCTK_REAL *RelTol) 
{
  int retval = ELL_NOSOLVER;

  switch (CCTK_GroupDimFromVarI(FieldIndex))
  {
    case 1:  
      CCTK_WARN(0,"SORMetric: No 1D SOR solver implemented");
      break;
    case 2:
      CCTK_WARN(0,"SORMetric: No 2D SOR solver implemented");
      break;
    case 3:
      retval = SORConfMetric3D(GH, MetricI, 0, 
                           FieldIndex, MIndex, NIndex,
                           AbsTol, RelTol);
      break;
    default:
      CCTK_WARN(1,"SORMetric: Solver only implemented for 3D");
      break;
  }

  return retval;
}
 
/*@@
   @routine    SORFlat
   @date       Tue Sep 26 11:31:42 2000
   @author     Gerd Lanfermann
   @desc 
     Elliptic solver wrapper which provides a interface to the 
     different n-dimensional SOR solvers (flat case), 
     of which only 3D is coded currently.

     This wrapper is registered and if it is being called with 
     a n-dimensional grid function, it then picks the correct solver.

     We pass in the arguments that are necessary for this class of 
     elliptic equations this solver is intended to solve. 
     See ./CactusElliptic/EllBase/src/ for the classes of elliptic equations.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

int SORFlat(cGH *GH, 
            int FieldIndex, 
            int MIndex, 
            int NIndex, 
            CCTK_REAL *AbsTol, 
            CCTK_REAL *RelTol) 
{
  int retval;

  retval = ELL_NOSOLVER;

  switch (CCTK_GroupDimFromVarI(FieldIndex))
  {
    case 1:  
      CCTK_WARN(1,"SORFlat: No 1D SOR solver implemented");
      break;
    case 2:
      CCTK_WARN(1,"SORFlat: No 2D SOR solver implemented");
      break;
    case 3:
      retval = SORFlat3D(GH, FieldIndex, MIndex, NIndex, AbsTol, RelTol);
      break;
  default:
      CCTK_WARN(1,"SORFlat: Solver only implemented for 3D");
      break;
  }

  return retval;

}
