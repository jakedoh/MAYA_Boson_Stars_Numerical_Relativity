 /*@@
   @file    WaveToy.c
   @date
   @author  Tom Goodale
   @desc
            Evolution routines for the wave equation solver
   @enddesc
   @version $Id$
 @@*/


#include <stddef.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
  DECLARE_CCTK_ARGUMENTS_CHECKED(Func) DECLARE_CCTK_ARGUMENTS
#endif

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusWave_WaveToyC_WaveToy_c);

void WaveToyC_Boundaries(CCTK_ARGUMENTS);
void  WaveToyC_Evolution(CCTK_ARGUMENTS);

 /*@@
   @routine WaveToyC_Evolution
   @date
   @author  Tom Goodale
   @desc
            Evolution for the wave equation
   @enddesc
   @calls   CCTK_SyncGroup, WaveToyC_Boundaries
@@*/

void  WaveToyC_Evolution(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_CHECKED(WaveToyC_Evolution);

  int i,j,k;
  int vindex;
  int istart, jstart, kstart, iend, jend, kend;
  CCTK_REAL dx,dy,dz,dt,dx2,dy2,dz2,dt2;
  CCTK_REAL dx2i,dy2i,dz2i;

  CCTK_REAL factor;

  /* Set up shorthands */
  dx = CCTK_DELTA_SPACE(0);
  dy = CCTK_DELTA_SPACE(1);
  dz = CCTK_DELTA_SPACE(2);
  dt = CCTK_DELTA_TIME;

  dx2 = dx*dx;
  dy2 = dy*dy;
  dz2 = dz*dz;
  dt2 = dt*dt;

  dx2i = 1.0/dx2;
  dy2i = 1.0/dy2;
  dz2i = 1.0/dz2;

  istart = 1;
  jstart = 1;
  kstart = 1;

  iend = cctk_lsh[0]-1;
  jend = cctk_lsh[1]-1;
  kend = cctk_lsh[2]-1;

  /* Do the evolution */
  factor = 2*(1 - (dt2)*(dx2i + dy2i + dz2i));

  for (k=kstart; k<kend; k++)
  {
    for (j=jstart; j<jend; j++)
    {
      for (i=istart; i<iend; i++)
      {
        vindex = CCTK_GFINDEX3D(cctkGH,i,j,k);

        phi[vindex] = factor*
                     phi_p[vindex] - phi_p_p[vindex]
                     + (dt2) *
                   ( ( phi_p[CCTK_GFINDEX3D(cctkGH,i+1,j  ,k  )]
                      +phi_p[CCTK_GFINDEX3D(cctkGH,i-1,j  ,k  )] )*dx2i
                    +( phi_p[CCTK_GFINDEX3D(cctkGH,i  ,j+1,k  )]
                      +phi_p[CCTK_GFINDEX3D(cctkGH,i  ,j-1,k  )] )*dy2i
                    +( phi_p[CCTK_GFINDEX3D(cctkGH,i  ,j  ,k+1)]
                      +phi_p[CCTK_GFINDEX3D(cctkGH,i  ,j,  k-1)] )*dz2i);
      }
    }
  }
}


 /*@@
   @routine WaveToyC_Boundaries
   @date
   @author  Tom Goodale
   @desc
            Boundary conditions for the wave equation
   @enddesc
@@*/

void WaveToyC_Boundaries(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;

  const char *bctype;

  bctype = NULL;
  if (CCTK_EQUALS(bound,"flat") || CCTK_EQUALS(bound,"static") ||
      CCTK_EQUALS(bound,"radiation") || CCTK_EQUALS(bound,"robin") ||
      CCTK_EQUALS(bound,"none"))
  {
    bctype = bound;
  }
  else if (CCTK_EQUALS(bound,"zero"))
  {
    bctype = "scalar";
  }

  /* Uses all default arguments, so invalid table handle -1 can be passed */
  if (bctype && Boundary_SelectVarForBC (cctkGH, CCTK_ALL_FACES, 1, -1,
                                         "wavetoy::phi", bctype) < 0)
  {
    CCTK_ERROR("Failed to register bound BC for wavemol::scalarevolvemol_scalar!");
  }
}
