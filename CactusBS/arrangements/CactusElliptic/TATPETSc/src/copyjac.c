/* (C) 2003-11-10 Erik Schnetter <schnetter@aei.mpg.de> */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include <petscdmda.h>

#include <cctk.h>
#include <cctk_Parameters.h>

#include "TATPETSc.h"



int TATPETSc_copyjac (Mat J, void *userptr)
{
  DECLARE_CCTK_PARAMETERS;
  
  userdata *user = (userdata*)userptr;
  const cGH *cctkGH = user->cctkGH;
  
  const PetscInt nvars = user->nvars;
  
  PetscInt ierr;
  
  
  
  assert (user->magic==MAGIC);
  
  if (veryverbose) CCTK_INFO ("*** TATPETSc_copyjac");
  
  
  
  assert (user);
  assert (user->nvars >= 0);
  assert (user->jac);
  
  
  
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
  PetscInt xs[DIM], xm[DIM];    /* offset, local extent */
  ierr = DMDAGetCorners (user->da, &xs[0],&xs[1],&xs[2], &xm[0],&xm[1],&xm[2]);
  CHKERRQ(ierr);
  
  /* Global PETSc boundaries */
  PetscInt XM[DIM];             /* global extent */
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
    var[n] = CCTK_VarDataPtrI(user->cctkGH, 0, user->jac0[n]);
    assert (var[n]);
  }
  
  
  
  /* Copy Cactus variable into PETSc variable */
  PetscInt * restrict rows;
  ierr = PetscMalloc(nvars*sizeof(int),&rows);CHKERRQ(ierr);
  PetscInt * restrict cols;
  ierr = PetscMalloc(27*nvars*sizeof(int),&cols);CHKERRQ(ierr);
  PetscScalar * restrict vals;
  assert (sizeof(PetscScalar) == sizeof(CCTK_REAL));
  ierr = PetscMalloc(27*nvars*sizeof(PetscScalar),&vals);CHKERRQ(ierr);
  AO ao;
  ierr = DMDAGetAO(user->da,&ao);CHKERRQ(ierr);
  
#if PETSC_VERSION_MAJOR < 3
  ierr = MatSetOption(J,MAT_ROWS_SORTED);CHKERRQ(ierr);
  ierr = MatSetOption(J,MAT_COLUMNS_SORTED);CHKERRQ(ierr);
/*   ierr = MatSetOption(J,MAT_IGNORE_ZERO_ENTRIES);CHKERRQ(ierr); */
#endif
  
  assert (nvars==1);
  
  for (PetscInt k=0; k<mi[2]; ++k) {
    for (PetscInt j=0; j<mi[1]; ++j) {
      for (PetscInt i=0; i<mi[0]; ++i) {
        
        for (PetscInt n=0; n<nvars; ++n) {
          
          const PetscInt rowind =
            n + nvars * (i0[0]+i + NI[0] * (i0[1]+j + NI[1] * (i0[2]+k)));
          rows[n] = rowind;
          
          for (PetscInt dk=-1; dk<=1; ++dk) {
            for (PetscInt dj=-1; dj<=1; ++dj) {
              for (PetscInt di=-1; di<=1; ++di) {
                const PetscInt m = di+1 + 3*(dj+1 + 3*(dk+1));
                const PetscInt ind =
                  bi[0]+i + ni[0] * (bi[1]+j + ni[1] * (bi[2]+k + ni[2] * m));
                const PetscInt colind =
                  n + nvars * (i0[0]+i+di + NI[0] * (i0[1]+j+dj +
                                                     NI[1] * (i0[2]+k+dk)));
                if (i+di>=0 && i+di<mi[0] && j+dj>=0 && j+dj<mi[1] && k+dk>=0 &&
                    k+dk<mi[2])
                {
                  cols[n+nvars*m] = colind;
                  vals[n+nvars*m] = var[n][ind];
                } else {
                  assert (var[n][ind] == 0);
                  cols[n+nvars*m] = -1;
                  vals[n+nvars*m] = 0;
                }
              }
            }
          }
        }
        
        ierr = AOApplicationToPetsc(ao,nvars,rows);CHKERRQ(ierr);
        ierr = AOApplicationToPetsc(ao,27*nvars,cols);CHKERRQ(ierr);
        ierr = MatSetValues (J,nvars,rows,27*nvars,cols,vals,INSERT_VALUES);
        CHKERRQ(ierr);
        
      }
    }
  }
  ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);  
  ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);  
  ierr = MatSetOption(J,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE);CHKERRQ(ierr);
  
  ierr = PetscFree(vals);CHKERRQ(ierr);
  ierr = PetscFree(rows);CHKERRQ(ierr);
  ierr = PetscFree(cols);CHKERRQ(ierr);
  
  
  
  /* Clean up */
  free (var);
  var = NULL;
  
  
  
  if (veryverbose) CCTK_INFO ("*** TATPETSc_copyjac done.");
  
  return 0;
}
