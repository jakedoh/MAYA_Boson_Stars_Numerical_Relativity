/* $Header$ */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



#include "noise.h"

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_FortranString.h"

#include "Symmetry.h"



/* #define DEBUG_BOUNDARY 1 */

static int ApplyBndNoise (const cGH *GH,
			  int stencil_dir,
			  const int *stencil,
			  int dir,
			  int first_var,
			  int num_vars);



int
BndNoiseVI (const cGH *GH, const int *stencil, int vi)
{
  int retval;
  retval = ApplyBndNoise (GH, -1, stencil, 0, vi, 1);
  return retval;
}

void
CCTK_FCALL CCTK_FNAME (BndNoiseVI) (int *ierr, const cGH **GH,
				    const int *stencil, const int *vi)
{
  *ierr = BndNoiseVI (*GH, stencil, *vi);
}

int
BndNoiseVN (const cGH *GH, const int *stencil, const char *vn)
{
  int vi, retval;
  vi = CCTK_VarIndex(vn);
  retval = BndNoiseVI (GH, stencil, vi);
  return retval;
}

void
CCTK_FCALL CCTK_FNAME (BndNoiseVN) (int *ierr, const cGH **GH,
				    const int *stencil, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (vn);
  *ierr = BndNoiseVN (*GH, stencil, vn);
  free (vn);
}



int
BndNoiseGI (const cGH *GH, const int *stencil, int gi)
{
  int first_vi, retval;

  first_vi = CCTK_FirstVarIndexI (gi);
  if (first_vi >= 0)
    {
      retval = ApplyBndNoise (GH, -1, stencil, 0, first_vi,
			      CCTK_NumVarsInGroupI (gi));
    }
  else
    {
      CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Invalid group index %d in BndNoiseGI", gi);
      retval = -1;
    }

  return (retval);
}



void
CCTK_FCALL CCTK_FNAME (BndNoiseGI) (int *ierr, const cGH **GH,
				    const int *stencil, const int *gi)
{
  *ierr = BndNoiseGI (*GH, stencil, *gi);
}


int
BndNoiseGN (const cGH *GH, const int *stencil, const char *gn)
{
  int gi, retval;

  gi = CCTK_GroupIndex (gn);
  if (gi >= 0)
    {
      retval = BndNoiseGI (GH, stencil, gi);
    }
  else
    {
      CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Invalid group name '%s' in BndNoiseGN", gn);
      retval = -1;
    }

  return (retval);
}

void CCTK_FCALL CCTK_FNAME (BndNoiseGN)
     (int *ierr,
      const cGH **GH,
      const int *stencil,
      ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (gn)
    *ierr = BndNoiseGN (*GH, stencil, gn);
  free (gn);
}





static int ApplyBndNoise (const cGH *GH,
			  int stencil_dir,
			  const int *stencil,
			  int dir,
			  int first_var,
			  int num_vars)
{
  DECLARE_CCTK_PARAMETERS;
  int i, j, k;
  int var, gindex, gdim, timelvl;
  int doBC[2*MAXDIM], lsh[MAXDIM];
  SymmetryGHex *sGHex;

  /* This argument is unused an undocumented; better make sure people
     don't try to use it for something.  */
  assert (stencil_dir == -1);

  /* get the group index of the variables */
  gindex = CCTK_GroupIndexFromVarI (first_var);

  /* get the number of dimensions */
  gdim      = CCTK_GroupDimI (gindex);

  /* make sure we can deal with this number of dimensions */
  if (gdim > MAXDIM)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "ApplyBndNoise: Variable dimension of %d not supported", gdim);
      return (-1);
    }

  /* check the direction parameter */
  if (abs (dir) > gdim)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "ApplyBndNoise: direction %d greater than dimension %d",
		  dir, gdim);
      return (-2);
    }

  /* initialize arrays for variables with less dimensions than MAXDIM
     so that we can use the INDEX_3D macro later on */
  for (i = gdim; i < MAXDIM; i++)
    {
      lsh[i] = 1;
    }

  /* get the current timelevel */
  timelvl = 0;

  /* see if we have a symmetry array */
  sGHex = (SymmetryGHex *) CCTK_GHExtension (GH, "Symmetry");


  /* now loop over all variables */
  for (var = first_var; var < first_var + num_vars; var++)
    {
      /* Apply condition if:
	 + boundary is not a symmetry boundary (no symmetry or unset(=unsed))
	 + boundary is a physical boundary
	 + have enough grid points
      */
      memset (doBC, 1, sizeof (doBC));
      if (sGHex)
	{
	  for (i = 0; i < 2 * gdim; i++)
	    {
	      doBC[i] = sGHex->GFSym[var][i] == GFSYM_NOSYM ||
		sGHex->GFSym[var][i] == GFSYM_UNSET;
	    }
	}
      for (i = 0; i < gdim; i++)
	{
	  lsh[i]       = GH->cctk_lsh[i];
	  doBC[i*2]   &= GH->cctk_lsh[i] > 1 && GH->cctk_bbox[i*2];
	  doBC[i*2+1] &= GH->cctk_lsh[i] > 1 && GH->cctk_bbox[i*2+1];
	  if (dir != 0)
	    {
	      doBC[i*2]   &= i+1 == -dir;
	      doBC[i*2+1] &= i+1 ==  dir;
	    }
	}

      /* now apply the boundaries face by face */
      if (gdim > 0)
	{
#ifdef DEBUG_BOUNDARY
	  if (doBC[0])
	    {
	      printf("Boundary: Applying lower x noise to boundary\n");
	    }
	  if (doBC[1])
	    {
	      printf("Boundary: Applying upper x noise to boundary\n");
	    }
#endif /* DEBUG_BOUNDARY */
	  /* lower x */
	  BOUNDARY_NOISE (doBC[0], stencil[0], lsh[1], lsh[2],
			  i, j, k);
	  /* upper x */
	  BOUNDARY_NOISE (doBC[1], stencil[0], lsh[1], lsh[2],
			  lsh[0]-i-1, j, k);

	}
      if (gdim > 1)

	{
#ifdef DEBUG_BOUNDARY
	  if (doBC[2])
	    {
	      printf("Boundary: Applying lower y noise to boundary\n");
	    }
	  if (doBC[3])
	    {
	      printf("Boundary: Applying upper y noise to boundary\n");
	    }
#endif /* DEBUG_BOUNDARY */
	  /* lower y */
	  BOUNDARY_NOISE (doBC[2], lsh[0], stencil[1], lsh[2],
			  i, j, k);
	  /* upper y */
	  BOUNDARY_NOISE (doBC[3], lsh[0], stencil[1], lsh[2],
			  i, lsh[1]-j-1, k);
	}
      if (gdim > 2)
	{
#ifdef DEBUG_BOUNDARY
	  if (doBC[4])
	    {
	      printf("Boundary: Applying lower z noise to boundary\n");
	    }
	  if (doBC[5])
	    {
	      printf("Boundary: Applying upper z noise to boundary\n");
	    }
#endif /* DEBUG_BOUNDARY */
	  /* lower z */
	  BOUNDARY_NOISE (doBC[4], lsh[0], lsh[1], stencil[2],
			  i, j, k);
	  /* upper z */
	  BOUNDARY_NOISE (doBC[5], lsh[0], lsh[1], stencil[2],
			  i, j, lsh[2]-k-1);
	}
    }

  return(0);
}


static void
add_bc_noise_to_var (int idx, const char* optstring, void* cctkGH)
{
  DECLARE_CCTK_PARAMETERS;
  cGH* GH = cctkGH;
  int sw[3];

  /* Change type from CCTK_INT to int */
  sw[0] = noise_stencil[0];
  sw[1] = noise_stencil[1];
  sw[2] = noise_stencil[2];

  if (noise_boundaries[0]) ApplyBndNoise (GH, -1, sw, -1, idx, 1);
  if (noise_boundaries[1]) ApplyBndNoise (GH, -1, sw, +1, idx, 1);
  if (noise_boundaries[2]) ApplyBndNoise (GH, -1, sw, -2, idx, 1);
  if (noise_boundaries[3]) ApplyBndNoise (GH, -1, sw, +2, idx, 1);
  if (noise_boundaries[4]) ApplyBndNoise (GH, -1, sw, -3, idx, 1);
  if (noise_boundaries[5]) ApplyBndNoise (GH, -1, sw, +3, idx, 1);
}


void
bc_noise(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

/*   Boundary_MakeSureThatTheSelectionIsEmpty(); */
  if (CCTK_TraverseString(bc_vars, add_bc_noise_to_var, cctkGH,
      		    CCTK_GROUP_OR_VAR) < 0)
    {
      CCTK_WARN (1, "Failed to parse 'Noise::bc_vars' parameter");
    }
}
