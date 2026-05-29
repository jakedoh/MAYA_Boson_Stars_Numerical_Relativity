/*@@
  @file       Caller.c
   @date      Thu 13 June 2002
   @author    Thomas Radke
   @desc
              C interface of the TestFpointerNULL thorn.
   @enddesc
   @version   $Id$
 @@*/

#include <stdio.h>
#include <string.h>

#include "cctk.h"
#include "cctk_FortranString.h"

/* the rcsid and the macro to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusTest_TestFpointerNULL_Caller_c)


/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/
void TestFpointerNULL_Caller (cGH *GH);
void CCTK_FCALL CCTK_FNAME (TestFpointerNULL_Callee)
                           (CCTK_INT *null_scalar, CCTK_INT *scalar,
                            CCTK_REAL *null_array, CCTK_REAL *array,
                            const int *lsh);


 /*@@
   @routine    TestFpointerNULL_Caller
   @author     Thomas Radke
   @date       Thu 13 June 2002
   @desc
               C function which initializes a scalar and a 2D array to zero
               and passes them on to a Fortran routine for modification.
               In addition, a dummy scalar and array are passed as NULL
               pointers. These should not be touched by the Fortran routine.
   @enddesc
   @calls      TestFpointerNULL_Callee

   @var        GH
   @vdesc      Pointer to CCTK GH
   @vtype      const cGH *
   @vio        in
   @endvar
@@*/
void TestFpointerNULL_Caller (cGH *GH)
{
  CCTK_INT *scalar;
  CCTK_REAL *array;
  int lsh[2];


  /* get the grid variables' pointers */
  scalar = (CCTK_INT *) CCTK_VarDataPtr (GH, 0, "TestFpointerNULL::scalar");
  array = (CCTK_REAL *) CCTK_VarDataPtr (GH, 0, "TestFpointerNULL::array");
  CCTK_GrouplshVI (GH, 2, lsh, CCTK_VarIndex ("TestFpointerNULL::array"));

  /* initialize the grid variables */
  *scalar = 0;
  memset (array, 0, lsh[0] * lsh[1] * sizeof (CCTK_REAL));

  /* call the Fortran routine */
  CCTK_FNAME (TestFpointerNULL_Callee) (NULL, scalar, NULL, array, lsh);
}
