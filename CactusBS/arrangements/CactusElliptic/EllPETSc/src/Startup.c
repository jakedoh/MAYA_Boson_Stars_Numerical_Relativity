#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cctk.h>
#include <cctk_Parameters.h>

#include "EllBase.h"

void petsc_confmetric(cGH *GH, int *MetricPsiI, int *FieldI, int *MI,
                      int *NI, int *AbsTol, int *RelTol);

void petsc_metric(cGH *GH, int *MetricI, int *FieldI, int *MI,
                  int *NI, int *AbsTol, int *RelTol);

void petsc_flat(cGH *GH, int *FieldIndex, int *MIndex, int *NIndex, 
                int *AbsTol, int *RelTol);

/* Registration of the petsc solvers with the Elliptic solver registry.
   This routine registers petsc_confmetric under the name "petsc" for 
   the class of elliptic equations "LinConfMetric" */
   
void EllPETSc_Register(cGH *GH) 
{
  DECLARE_CCTK_PARAMETERS 
  
  Ell_RegisterSolver(petsc_confmetric,"petsc","Ell_LinConfMetric");
  Ell_RegisterSolver(petsc_metric,"petsc","Ell_LinMetric");
  Ell_RegisterSolver(petsc_flat,"petsc","Ell_LinFlat");
}
