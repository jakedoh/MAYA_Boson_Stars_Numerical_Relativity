 /*@@
   @file      WaveBinary.c
   @date      
   @author    Cactus Maintainers
   @desc 
   Add a binary source term to the 3D scalar wave equation
   @enddesc 
   @version $Header$
 @@*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusWave_WaveBinarySource_WaveBinary_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/


/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/


/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/


/********************************************************************
 *****************   External Routines Prototype ********************
 ********************************************************************/

void WaveBinaryC(CCTK_ARGUMENTS);


/********************************************************************
 **********************   External Routines  ************************
 ********************************************************************/

void WaveBinaryC(CCTK_ARGUMENTS)
{
  
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static int firstcall=1;

  const CCTK_REAL binary_pi=3.141592653589793;
  CCTK_REAL charge_factor;
  CCTK_REAL xsm,ysm,zsm,rad2m;
  CCTK_REAL xsp,ysp,zsp,rad2p;

  int i,j,k,vindex;

  charge_factor = 3.0*binary_charge/(4.0*binary_pi*pow(binary_size,3));

  /* Calculate the centers of the binary sources  */
  xsm = - binary_radius * cos(binary_omega*cctk_time);
  ysm = - binary_radius * sin(binary_omega*cctk_time);
  zsm = 0;
  xsp = + binary_radius * cos(binary_omega*cctk_time);
  ysp = + binary_radius * sin(binary_omega*cctk_time);
  zsp = 0;

  if ((CCTK_EQUALS(binary_verbose, "yes") && firstcall) ||
      CCTK_EQUALS(binary_verbose, "debug"))
  {
    CCTK_VInfo(CCTK_THORNSTRING,
               "Charge factor: %g\n",(double)charge_factor);
    CCTK_VInfo(CCTK_THORNSTRING,
               "Charge center (-): %g %g %g\n",
               (double)xsm,(double)ysm,(double)zsm);
    CCTK_VInfo(CCTK_THORNSTRING,
               "Charge center (+): %g %g %g\n",
               (double)xsp,(double)ysp,(double)zsp);
    CCTK_VInfo(CCTK_THORNSTRING,
               "Charge extent: %g\n",(double)binary_size);
  }
  firstcall=0;

  for (k=0;k<cctk_lsh[2];k++)
  {
    for (j=0;j<cctk_lsh[1];j++)
    {
      for (i=0;i<cctk_lsh[0];i++)
      {
        vindex = CCTK_GFINDEX3D(cctkGH,i,j,k);

        rad2m =
          pow(x[vindex]-xsm,2) + pow(y[vindex]-ysm,2) + pow(z[vindex]-zsm,2);
        rad2p =
          pow(x[vindex]-xsp,2) + pow(y[vindex]-ysp,2) + pow(z[vindex]-zsp,2);

        /* Note that both sources are positive, leading to a net
           monopole moment  */
        if (rad2m<pow(binary_size,2))
        {
          phi[vindex] += 0.5*pow(CCTK_DELTA_TIME,2)*charge_factor;
        }
        if (rad2p<pow(binary_size,2))
        {
          phi[vindex] += 0.5*pow(CCTK_DELTA_TIME,2)*charge_factor;
        }
      }
    }
  }

  /* Note that we do not need to sync anything, since each grid patch
     has filled out its ghostzones  */
}
