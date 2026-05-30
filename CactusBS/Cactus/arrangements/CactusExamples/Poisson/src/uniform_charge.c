#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <util_ErrorCodes.h>
#include <util_Table.h>

#include <carpet.h>

#include <TATelliptic.h>



static int calc_residual (const cGH * const cctkGH,
                          int const options_table,
                          void * const userdata);
static int apply_bounds (const cGH * const cctkGH,
                         int const options_table,
                         void * const userdata);
static void apply_bounds_level (CCTK_ARGUMENTS);



void Poisson_prepare (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_Poisson_prepare;
  DECLARE_CCTK_PARAMETERS;
  
  /* Initial data for the solver */
#pragma omp parallel for collapse(3)
  for (int k=0; k<cctk_lsh[2]; ++k) {
    for (int j=0; j<cctk_lsh[1]; ++j) {
      for (int i=0; i<cctk_lsh[0]; ++i) {
	int ipos = CCTK_GFINDEX3D(cctkGH,i,j,k);
	phi[ipos] = 0;
      }
    }
  }
}



void Poisson_solve (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_Poisson_solve;
  DECLARE_CCTK_PARAMETERS;
  
  int ierr;
  
  /* Grid variables for the solver */
#define NVAR 1                  /* number of equations to solve */
  int var_ind[NVAR];		/* index of variable */
  int res_ind[NVAR];		/* index of residual */
  var_ind[0] = CCTK_VarIndex ("Poisson::phi");
  res_ind[0] = CCTK_VarIndex ("Poisson::res");
  
  /* Options for the solver */
  int options_table = Util_TableCreateFromString (options);
  assert (options_table>=0);
  
  /* Call solver */
  ierr = TATelliptic_CallSolver (cctkGH, var_ind, res_ind, NVAR,
				 options_table,
				 calc_residual, apply_bounds, 0,
				 solver);
  if (ierr!=0) {
    CCTK_WARN (CCTK_WARN_ALERT, "Failed to solve elliptic equation");
  }
  
  ierr = Util_TableDestroy (options_table);
  assert (!ierr);
}



/* Caculate the residual */
int calc_residual (const cGH * const cctkGH,
                   int const options_table,
                   void * const userdata)
{
  DECLARE_CCTK_ARGUMENTS_Poisson_solve;
  DECLARE_CCTK_PARAMETERS;
  
  /* charge density */
  CCTK_REAL rho = charge / (4.0/3.0 * M_PI * pow(radius,3));
  
  /* offsets for 3D array layout */
  int di = CCTK_GFINDEX3D(cctkGH,1,0,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  int dj = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  int dk = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  
  CCTK_REAL idx2[3];		/* inverse squared grid spacing */
  for (int d=0; d<3; ++d) {
    idx2[d] = 1.0 / pow(CCTK_DELTA_SPACE(d), 2);
  }
  
  for (int d=0; d<3; ++d) {
    assert (cctk_nghostzones[d] >= 1);
  }

  /* Initialize residual to zero:
     We only do this to ensure the boundaries of the residual are
     defined, and are not nan. That in turn is only needed to make
     output look nicer. */
#pragma omp parallel for collapse(3)
  for (int k=0; k<cctk_lsh[2]; ++k) {
    for (int j=0; j<cctk_lsh[1]; ++j) {
      for (int i=0; i<cctk_lsh[0]; ++i) {
	int ipos = CCTK_GFINDEX3D(cctkGH,i,j,k);
        res[ipos] = 0.0;
      }
    }
  }
  
  /* Calculate residual:
     res = Laplace phi + 4 pi rho */
#pragma omp parallel for collapse(3)
  for (int k=cctk_nghostzones[2]; k<cctk_lsh[2]-cctk_nghostzones[2]; ++k) {
    for (int j=cctk_nghostzones[1]; j<cctk_lsh[1]-cctk_nghostzones[1]; ++j) {
      for (int i=cctk_nghostzones[0]; i<cctk_lsh[0]-cctk_nghostzones[0]; ++i) {
	int ipos = CCTK_GFINDEX3D(cctkGH,i,j,k);
	res[ipos] = 
	  + (phi[ipos-di] - 2*phi[ipos] + phi[ipos+di]) * idx2[0]
	  + (phi[ipos-dj] - 2*phi[ipos] + phi[ipos+dj]) * idx2[1]
	  + (phi[ipos-dk] - 2*phi[ipos] + phi[ipos+dk]) * idx2[2]
	  + (r[ipos]<radius ? 4 * M_PI * rho : 0);
      }
    }
  }
  
  /* Success */
  return 0;
}



/* Apply the boundary conditions to the solution */
int apply_bounds (const cGH * const cctkGH,
                  int const options_table,
                  void * const userdata)
{
  DECLARE_CCTK_ARGUMENTS_Poisson_solve;
  DECLARE_CCTK_PARAMETERS;
  
  int ierr;
  
  /* Switch to level mode */
  ierr = CallLevelFunction ((cGH *)cctkGH, apply_bounds_level);
  assert (!ierr);
  
  /* Success */
  return 0;
}

void apply_bounds_level (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_Poisson_solve;
  DECLARE_CCTK_PARAMETERS;
  
  int ierr;
  
  /* Call boundary condition schedule group */
  ierr = CallScheduleGroup (cctkGH, "Poisson_boundaries");
  assert (!ierr);
}

void Poisson_boundaries_select (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_Poisson_boundaries_select;
  DECLARE_CCTK_PARAMETERS;
  
  int ierr;
  
  /* Select Direchlet boundary conditions */
  ierr = Boundary_SelectGroupForBC
    (cctkGH, CCTK_ALL_FACES, 1, -1, "Poisson::potential", "scalar");
  assert (!ierr);
}
