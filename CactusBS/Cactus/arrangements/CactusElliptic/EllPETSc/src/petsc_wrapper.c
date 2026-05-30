
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Flesh.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusElliptic_EllPETSc_petsc_wrapper_c)

void petsc_confmetric(cGH *GH, int *MetricPsiI, int FieldIndex, 
                      int MIndex, int NIndex, 
                      CCTK_REAL *AbsTol, CCTK_REAL *RelTol);

/* The wrapper functions for the core PETSc solver, that performs the 
   actual solve. One wrapper function for each class, because these 
   functions are registered with a fixed set up arguments. Still, they 
   all call petsc_confmetric_solver as the core solver. */

/* Wrapper function for the class of elliptic equations that needs a metric */
void petsc_metric(cGH *GH, int *MetricI, int FieldIndex, 
                  int MIndex, int NIndex, 
                  CCTK_REAL *AbsTol, CCTK_REAL *RelTol) {

/* Arrays for the M/N size info */
  
  int petsc_confmetric_solver(cGH *GH, int *MetricI, int MetricISize, 
                               int FieldIndex, int MIndex, int NIndex, 
                               CCTK_REAL *AbsTol, CCTK_REAL *RelTol);
  
  
  if (GH->cctk_dim>3) 
  {
   CCTK_WARN(0,"This elliptic solver implementation does not do dimension>3!");
  }

  petsc_confmetric_solver(GH, MetricI, 6, FieldIndex, 
                          MIndex, NIndex, 
                          AbsTol, RelTol);
}
  
/* wrapper function for the class of elliptic equations, that needs a conf.
   factor */
void petsc_confmetric(cGH *GH, int *MetricPsiI, int FieldIndex, 
                      int MIndex, int NIndex, 
                      CCTK_REAL *AbsTol, CCTK_REAL *RelTol) {
  
  int petsc_confmetric_solver(cGH *, int *, int, int, int, int, 
                               CCTK_REAL *, CCTK_REAL *);

 /* petsc_confmetric_solver expects an metric array, in the case of this 
    equation class,  needs conf.factor as last entry -> size 7 */
  petsc_confmetric_solver(GH, MetricPsiI, 7,
                          FieldIndex, MIndex, NIndex, 
                          AbsTol, RelTol);
}
