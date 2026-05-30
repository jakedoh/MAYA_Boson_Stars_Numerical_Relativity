 /*@@
   @file    WaveToy.cc
   @date
   @author  Tom Goodale
   @desc
            Evolution routines for the wave equation solver
   @enddesc
   @version $Id$
 @@*/

#include <stddef.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#define val(gridfunc,i,j,k)  gridfunc[CCTK_GFINDEX3D(cctkGH,i,j,k)]

 /*@@
   @routine    WaveToyC_Evolution
   @date
   @author     Tom Goodale
   @desc
               Evolution for the wave equation
   @enddesc
@@*/

extern "C" void  WaveToyCXX_Evolution(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  // Set up shorthands

  CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  CCTK_REAL dz = CCTK_DELTA_SPACE(2);
  CCTK_REAL dt = CCTK_DELTA_TIME;

  CCTK_REAL dx2 = dx*dx;
  CCTK_REAL dy2 = dy*dy;
  CCTK_REAL dz2 = dz*dz;
  CCTK_REAL dt2 = dt*dt;

  CCTK_REAL dx2i = 1.0/dx2;
  CCTK_REAL dy2i = 1.0/dy2;
  CCTK_REAL dz2i = 1.0/dz2;

  CCTK_REAL factor = 2*(1 - (dt2)*(dx2i + dy2i + dz2i));

  int istart = 1;
  int jstart = 1;
  int kstart = 1;

  int iend = cctk_lsh[0]-1;
  int jend = cctk_lsh[1]-1;
  int kend = cctk_lsh[2]-1;

  //
  // Do the evolution
  //

  for (int k=kstart; k<kend; k++)
  {
    for (int j=jstart; j<jend; j++)
    {
      for (int i=istart; i<iend; i++)
      {
        int vindex = CCTK_GFINDEX3D(cctkGH,i,j,k);

        phi[vindex] =
	  factor*phi_p[vindex] - phi_p_p[vindex]
	  + dt2 *
	  ( (  val( phi_p, i+1,j  ,k  ) + val( phi_p, i-1,j  ,k)   )*dx2i
	    +( val( phi_p, i  ,j+1,k  ) + val( phi_p, i  ,j-1,k)   )*dy2i
	    +( val( phi_p, i  ,j  ,k+1) + val( phi_p, i  ,j,  k-1) )*dz2i
	    );
      }
    }
  }
}

 /*@@
   @routine    WaveToyC_Boundaries
   @date
   @author     Tom Goodale
   @desc
               Boundary conditions for the wave equation
   @enddesc
@@*/

extern "C" void WaveToyCXX_Boundaries(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
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

  // Uses all default arguments, so invalid table handle -1 can be passed
  if (bctype && Boundary_SelectVarForBC (cctkGH, CCTK_ALL_FACES, 1, -1,
                                         "wavetoy::phi", bctype) < 0)
  {
    CCTK_WARN (0,"WaveToyCXX_Boundaries: Error selecting boundary condition");
  }
}
