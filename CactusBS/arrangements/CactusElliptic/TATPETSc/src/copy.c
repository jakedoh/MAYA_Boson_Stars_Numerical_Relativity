/* (C) 2001-04-18 Erik Schnetter <schnetter@uni-tuebingen.de> */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include <petscdmda.h>

#include <cctk.h>
#include <cctk_Parameters.h>

#include "TATPETSc.h"



int TATPETSc_copy (Vec x, void *userptr, TATdir dir, TATvarset varset)
{
  DECLARE_CCTK_PARAMETERS;
  
  userdata *user = (userdata*)userptr;
  const cGH *cctkGH = user->cctkGH;
  
  const PetscInt nvars = user->nvars;

  PetscInt ierr;
  
  assert (user->magic==MAGIC);
  
  if (veryverbose) CCTK_INFO ("*** TATPETSc_copy");
  
  assert (dir==TATcopyin || dir==TATcopyout);
  assert (varset==TATcopyvars || varset==TATcopyvals);
  
  if (veryverbose) {
    CCTK_VInfo (CCTK_THORNSTRING, "copy %s %s PETSc",
		varset==TATcopyvars ? "variables" : "function values",
		dir==TATcopyin ? "into" : "out of");
  }
  
  
  
  assert (user);
  assert (user->nvars >= 0);
  assert (user->var);
  assert (user->val);
  
  
  
  /* Get grid boundaries */
  
  /* Cactus grid extent */
  PetscInt NI[DIM];          /* global extent */
  PetscInt ni[DIM];          /* local extent including ghost zones */
  PetscInt gi[DIM];          /* number of ghost zones */
  PetscInt i0[DIM];          /* local-to-global offset */
  PetscInt bi[DIM], ei[DIM]; /* lower and upper number of ghost zones */
  PetscInt mi[DIM];          /* local extent without ghost zones */
  
  for (PetscInt d=0; d<user->dyndata.dim; ++d) {
    
    /* Local Cactus boundaries */
    ni[d] = user->dyndata.lsh[d];
    i0[d] = user->dyndata.lbnd[d];
    gi[d] = user->dyndata.nghostzones[d];
    bi[d] = user->dyndata.bbox[2*d  ] ? user->nboundaryzones[2*d  ] : gi[d];
    ei[d] = user->dyndata.bbox[2*d+1] ? user->nboundaryzones[2*d+1] : gi[d];
    mi[d] = ni[d]-bi[d]-ei[d];
    
    /* Global Cactus boundaries */
    NI[d] = (user->dyndata.gsh[d]
	     - (user->nboundaryzones[2*d] + user->nboundaryzones[2*d+1]));
  }
  
  for (PetscInt d=user->dyndata.dim; d<DIM; ++d) {
    
    /* Local Cactus boundaries */
    ni[d] = 1;
    i0[d] = 0;
    gi[d] = 0;
    bi[d] = 0;
    ei[d] = 0;
    mi[d] = ni[d]-bi[d]-ei[d];
    
    /* Global Cactus boundaries */
    NI[d] = 1;
  }
  
  /* Local PETSc boundaries */
  
  /* PETSc grid extent */
  PetscInt XM[DIM];             /* global extent */
  PetscInt xs[DIM], xm[DIM];    /* offset, local extent */
  
  ierr = DMDAGetCorners (user->da, &xs[0],&xs[1],&xs[2], &xm[0],&xm[1],&xm[2]);
  CHKERRQ(ierr);
  
  /* Global PETSc boundaries */
  ierr = DMDAGetInfo (user->da,
                      PETSC_NULL,
                      &XM[0],&XM[1],&XM[2],
                      PETSC_NULL,PETSC_NULL,PETSC_NULL,
                      PETSC_NULL,PETSC_NULL,
                      PETSC_NULL,PETSC_NULL,PETSC_NULL,
                      PETSC_NULL);
  CHKERRQ(ierr);
  
  if (veryverbose) {
    CCTK_VInfo (CCTK_THORNSTRING,
		"Dimension: dim=%d", user->dyndata.dim);
    CCTK_VInfo (CCTK_THORNSTRING,
		"Global Cactus shape: extent=[%td,%td,%td]",
		(ptrdiff_t)NI[0],(ptrdiff_t)NI[1],(ptrdiff_t)NI[2]);
    CCTK_VInfo (CCTK_THORNSTRING,
		"Global PETSc shape: extent=[%td,%td,%td]",
		(ptrdiff_t)XM[0],(ptrdiff_t)XM[1],(ptrdiff_t)XM[2]);
    CCTK_VInfo (CCTK_THORNSTRING,
		"Local Cactus shape: start=[%td,%td,%td], extent=[%td,%td,%td]",
		(ptrdiff_t)i0[0],(ptrdiff_t)i0[1],(ptrdiff_t)i0[2],
                (ptrdiff_t)ni[0],(ptrdiff_t)ni[1],(ptrdiff_t)ni[2]);
    CCTK_VInfo (CCTK_THORNSTRING,
		"Local PETSc shape: start=[%td,%td,%td], extent=[%td,%td,%td]",
		(ptrdiff_t)xs[0],(ptrdiff_t)xs[1],(ptrdiff_t)xs[2],
                (ptrdiff_t)xm[0],(ptrdiff_t)xm[1],(ptrdiff_t)xm[2]);
  }
  
  /* Global shapes must match */
  for (PetscInt d=0; d<user->dyndata.dim; ++d) {
    assert (XM[d] == NI[d]);
  }
  
  /* Local shapes must match */
  for (PetscInt d=0; d<user->dyndata.dim; ++d) {
    assert (xs[d] == i0[d]);
    assert (xm[d] == mi[d]);
  }
  
  
  
  /* Get Cactus variable pointers */
  CCTK_REAL * restrict * restrict var = malloc(sizeof(*var) * nvars);
  assert (var);
  for (PetscInt n=0; n<nvars; ++n) {
    switch (varset) {
    case TATcopyvars:
      /* variables */
      var[n] = CCTK_VarDataPtrI(user->cctkGH, 0, user->var[n]);
      assert (var[n]);
      break;
    case TATcopyvals:
      /* values */
      var[n] = CCTK_VarDataPtrI(user->cctkGH, 0, user->val[n]);
      assert (var[n]);
      break;
    default:
      assert (0);
    }
  }
  
  
  
  switch (dir) {
  case TATcopyin: {
    /* copy in */
    /* Address of PETSc vector */
    double *xx;

    /* Get local PETSc vector */
    if (veryverbose) CCTK_INFO ("VecGetArray");
    ierr = VecGetArray (x, &xx);
    CHKERRQ(ierr);

    /* Copy Cactus variable into PETSc variable */
    for (PetscInt k=0; k<mi[2]; ++k) {
      for (PetscInt j=0; j<mi[1]; ++j) {
	for (PetscInt i=0; i<mi[0]; ++i) {
	  for (PetscInt n=0; n<nvars; ++n) {
	    PetscInt ind = bi[0]+i + ni[0]*(bi[1]+j + ni[1]*(bi[2]+k));
	    xx[n+nvars*(i+mi[0]*(j+mi[1]*k))] = var[n][ind];
	  }
	}
      }
    }
    
    if (veryverbose) CCTK_INFO ("Destroy");
    ierr = VecRestoreArray (x, &xx);
    CHKERRQ(ierr);
    } // TATcopyin
    break;
    
  case TATcopyout: {
    /* copy out */
    /* Address of PETSc vector */
    const double *xx;
    
    /* Get local PETSc vector */
    if (veryverbose) CCTK_INFO ("VecGetArray");
    ierr = VecGetArrayRead (x, &xx);
    CHKERRQ(ierr);

    /* Copy Cactus variable */
    for (PetscInt k=0; k<mi[2]; ++k) {
      for (PetscInt j=0; j<mi[1]; ++j) {
	for (PetscInt i=0; i<mi[0]; ++i) {
	  for (PetscInt n=0; n<nvars; ++n) {
	    PetscInt ind = bi[0]+i + ni[0]*(bi[1]+j + ni[1]*(bi[2]+k));
	    var[n][ind] = xx[n+nvars*(i+mi[0]*(j+mi[1]*k))];
	  }
	}
      }
    }
    
    if (veryverbose) CCTK_INFO ("Destroy");
    ierr = VecRestoreArrayRead (x, &xx);
    CHKERRQ(ierr);

    /* This loop possibly does too much work: it might synchronise
       groups several times */
    for (PetscInt n=0; n<user->nvars; ++n) {
      CCTK_SyncGroupWithVarI (cctkGH, user->var[n]);
    }
    } // TATcopyout
    break;
    
  default:
    CCTK_VERROR ("Unknown copy direction: %d", (int)dir);
  }
  
  
  
  /* Clean up */
  free ((void*)var);
  var = NULL;
  
  
  
  if (veryverbose) CCTK_INFO ("*** TATPETSc_copy done.");
  
  return 0;
}
