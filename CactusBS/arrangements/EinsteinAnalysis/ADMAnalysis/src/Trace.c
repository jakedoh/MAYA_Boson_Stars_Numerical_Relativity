 /*@@
   @file      Trace.c
   @date      Thu Apr 25 16:35:04 2002
   @author    Tom Goodale
   @desc 
   Calculates the trace of things
   @enddesc 
   @version $Header$
 @@*/

#include "cctk.h"

#include "ADMAnalysis.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusEinstein_ADMAnalysis_Trace_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

#define SQR(a) ((a)*(a))

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    ADMAnalysis_Trace
   @date       Thu Apr 25 16:49:57 2002
   @author     Tom Goodale
   @desc 
   Calculates the trace of a symmetric 2 index 3-tensor.
   Optionally returns the trace of the metric at each point.
   @enddesc 
   @calls     
   @calledby   
   @history 
   @hdate Thu Apr 25 16:50:38 2002 @hauthor Tom Goodale
   @hdesc Took old evlatrK routine and made general 
   @endhistory 
    @var     ash
   @vdesc   grid size
   @vtype   const int *
   @vio     in
   @var     g11
   @vdesc   the 11 component of the metric tensor
   @vtype   const CCTK_REAL *
   @vio     in
   @var     g12
   @vdesc   the 12 component of the metric tensor
   @vtype   const CCTK_REAL *
   @vio     in
   @var     g13
   @vdesc   the 13 component of the metric tensor
   @vtype   const CCTK_REAL *
   @vio     in
   @var     g22
   @vdesc   the 22 component of the metric tensor
   @vtype   const CCTK_REAL *
   @vio     in
   @var     g23
   @vdesc   the 23 component of the metric tensor
   @vtype   const CCTK_REAL *
   @vio     in
   @var     g33
   @vdesc   the 33 component of the metric tensor
   @vtype   const CCTK_REAL *
   @vio     in
   @var     tensor11
   @vdesc   the 11 component of the tensor
   @vtype   CCTK_REAL *
   @vio     out
   @var     tensor12
   @vdesc   the 12 component of the tensor
   @vtype   CCTK_REAL *
   @vio     out
   @var     tensor13
   @vdesc   the 13 component of the tensor
   @vtype   CCTK_REAL *
   @vio     out
   @var     tensor22
   @vdesc   the 22 component of the tensor
   @vtype   CCTK_REAL *
   @vio     out
   @var     tensor23
   @vdesc   the 23 component of the tensor
   @vtype   CCTK_REAL *
   @vio     out
   @var     tensor33
   @vdesc   the 33 component of the tensor
   @vtype   CCTK_REAL *
   @vio     out
   @var     detg
   @vdesc   the determinant of the metric tensor
   @vtype   CCTK_REAL *
   @vio     out
   @vcomment 
   Will be ignored if NULL.
   @endvar 

 @@*/
void ADMAnalysis_Trace(const CCTK_INT *ash,
                       const CCTK_REAL *g11,
                       const CCTK_REAL *g12,
                       const CCTK_REAL *g13,
                       const CCTK_REAL *g22,
                       const CCTK_REAL *g23,
                       const CCTK_REAL *g33,
                       CCTK_REAL *tensor11,
                       CCTK_REAL *tensor12,
                       CCTK_REAL *tensor13,
                       CCTK_REAL *tensor22,
                       CCTK_REAL *tensor23,
                       CCTK_REAL *tensor33,
                       CCTK_REAL *trace,
                       CCTK_REAL *detg)
{
  
  int  i;
  CCTK_REAL det, u11, u12, u22, u13, u23, u33, two;
  CCTK_REAL g11p, g12p, g13p, g22p, g23p, g33p;
  
  two = 2.0;

  /* loop over all the gridpoints */
  for(i = 0; i< ash[0]*ash[1]*ash[2];i++)  
  {
       
    /* get the metric */
    g11p = g11[i];
    g12p = g12[i];
    g13p = g13[i];
    g22p = g22[i];
    g23p = g23[i];
    g33p = g33[i];
      
    /* compute determinant */
    det=-(g13p*g13p*g22p)+2.*g12p*g13p*g23p-g11p*g23p*
      g23p-g12p*g12p*g33p+g11p*g22p*g33p;
    
    /* invert metric. This is the conformal upper metric */
    u11=(-SQR(g23p) + g22p*g33p)/det;
    u12=(g13p*g23p - g12p*g33p)/det;
    u22=(-SQR(g13p) + g11p*g33p)/det;
    u13=(-g13p*g22p + g12p*g23p)/det;
    u23=(g12p*g13p - g11p*g23p)/det;
    u33=(-SQR(g12p) + g11p*g22p)/det;
    
    /* Calculate trK */
    trace[i] = (u11*tensor11[i] + u22*tensor22[i] +
                u33*tensor33[i]+ two*u12*tensor12[i] +
                two*u13*tensor13[i] + two*u23*tensor23[i]);
    if(detg)
    {
      detg[i]= det;
    }
  }
}
/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

