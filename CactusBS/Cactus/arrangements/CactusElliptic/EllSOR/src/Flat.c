 /*@@
   @file      Flat.c
   @date      Tue Aug 24 12:50:07 1999
   @author    Gerd Lanfermann
   @desc 
   SOR solver for 3D flat equation
   @enddesc 
 @@*/

/*#define DEBUG_ELLIPTIC*/

#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "Boundary.h"
#include "Symmetry.h"
#include "Ell_DBstructure.h"
#include "EllBase.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusElliptic_EllSOR_Flat_c)
 
int SORFlat3D(cGH *GH, int FieldIndex, int MIndex, int NIndex,
              CCTK_REAL *AbsTol, CCTK_REAL *RelTol);

int SORFlat3D(cGH *GH, int FieldIndex, int MIndex, int NIndex,
              CCTK_REAL *AbsTol, CCTK_REAL *RelTol)
{
  DECLARE_CCTK_PARAMETERS

  int retval = ELL_SUCCESS;

  /* The pointer to the data fields */
  CCTK_REAL *Mlin=NULL;
  CCTK_REAL *Nlin=NULL;   
  CCTK_REAL *var =NULL; 

  /* shortcuts for deltas,etc. */
  CCTK_REAL dx,dy,dz; 

  /* Some physical variables */
  CCTK_REAL omega, resnorm, residual, glob_residual, rjacobian; 
  CCTK_REAL finf;
  CCTK_INT npow;
  CCTK_REAL tol;

  /* Iteration and stepping  variables */
  int sorit;
  CCTK_INT maxit;
  int i, j, k;
  int origin_sign;

  /* stencil index */
  int ijk;
  int ipjk, ijpk, ijkp, imjk, ijmk, ijkm;

  /* Coefficients for the stencil...  */
  CCTK_REAL ac,ac_orig,aw,ae,an,as,at,ab;
 
  /* Miscellaneous */
  int sum_handle;
  int sw[3];
  int Mstorage=0, Nstorage=0;
  static int firstcall = 1;
  CCTK_REAL  dx2rec, dy2rec, dz2rec; 
  char *robin = NULL;
  char *sor_accel = NULL;
  int ierr = 0;
  const void* input_array[1];
  void*  reduction_value[1];
  CCTK_INT input_array_dim[1];
  CCTK_INT input_array_type_codes[1];

  /* Boundary conditions */
  int use_robin = 0;
  
  input_array[0]             = (CCTK_REAL *)&resnorm;
  input_array_dim[0]         = 0;
  reduction_value[0]         = (CCTK_REAL *)&glob_residual;
  input_array_type_codes[0]  = CCTK_VARIABLE_REAL;
      
  /* Avoid compiler warnings */
  RelTol = RelTol;

  /* Get the reduction handle */
  sum_handle = CCTK_LocalArrayReductionHandle("sum");
  if (sum_handle<0) 
  {
    CCTK_WARN(1,"SORFlat3D: Cannot get reduction handle for sum operation");
  }

  /* IF Robin BCs are set, prepare for a boundary call:
     setup stencil width and get Robin constants (set by the routine
     which is calling the solver interface) */

  if (Ell_GetStrKey(&robin,"EllLinFlat::Bnd::Robin") < 0)
  {
    CCTK_WARN(2,"SORFlat3D: Robin keys not set");
  }
  else
  {
    if (CCTK_Equals(robin,"yes"))
    {
      use_robin = 1;

      sw[0]=1; 
      sw[1]=1; 
      sw[2]=1;

      if (Ell_GetRealKey(&finf, "EllLinFlat::Bnd::Robin::inf") == 
          ELLGET_NOTSET)
      {
        CCTK_WARN(1,"SORFlat3D: EllLinFlat::Bnd::Robin::inf not set, "
                  "setting to 1");
        finf = 1;
      }
      if (Ell_GetIntKey(&npow, "EllLinFlat::Bnd::Robin::falloff") 
          == ELLGET_NOTSET)
      {
        CCTK_WARN(1,"SORFlat3D: EllLinFlat::Bnd::Robin::falloff not set, "
                  "setting to 1");
        npow = 1;
      }
    }
  }

  /* Get the maximum number of iterations allowed. */
  if (Ell_GetIntKey(&maxit, "Ell::SORmaxit") == ELLGET_NOTSET)
  {
    CCTK_WARN(1,"SORFlat3D: Ell::SORmaxit not set. "
              "Maximum allowed iterations being set to 100.");
    maxit = 100;
  }
  
  /* Only supports absolute tolerance */
  tol   = AbsTol[0];

  /* Get the type of acceleration to use. */
  if (Ell_GetStrKey(&sor_accel, "Ell::SORaccel") == ELLGET_NOTSET)
  {
    const char tmpstr[6] = "const";
    CCTK_WARN(3, "SORFlat3D: Ell::SORaccel not set. "
              "Omega being set to a constant value of 1.8.");
    sor_accel = strdup(tmpstr);
  }

  /* Things to do only once! */
  /* TODO: Need to handle this differently, since it may be called
           for different reasons by the same code. */
  if (firstcall==1) 
  {
    if (CCTK_Equals(elliptic_verbose, "yes"))
    {
      if (CCTK_Equals(sor_accel, "cheb"))
      {
        CCTK_INFO("SOR with Chebyshev acceleration");
      }
      else if (CCTK_Equals(sor_accel, "const"))
      {
        CCTK_INFO("SOR with hardcoded omega = 1.8");
      }
      else if (CCTK_Equals(sor_accel, "none"))
      {
        CCTK_INFO("SOR with unaccelerated relaxation (omega = 1)");
      }
      else
      {
         CCTK_WARN(0, "SORFlat3D: sor_accel not set");
      }
    }
    firstcall = 0;
  }

  /* Get the data ptr of these GFs, They all have to be
     on the same timelevel; if we have a negative index for M/N, 
     this GF is not set,  there for don't even look for it and flag it  */

  var = (CCTK_REAL *)CCTK_VarDataPtrI(GH, 0, FieldIndex);
  if (var == NULL)
  {
    CCTK_WARN(0, "SORFlat3D: Unable to get pointer to GF variable");
  }

  if (MIndex>=0)  
  { 
    Mlin = (CCTK_REAL *)CCTK_VarDataPtrI(GH, 0, MIndex);
    if (Mlin)
    {
      Mstorage = 1;
    }
    else
    {
      CCTK_WARN(0, "SORFlat3D: Unable to get pointer to M");
    }
  }

  if (NIndex>=0) 
  {
    Nlin = (CCTK_REAL *)CCTK_VarDataPtrI(GH, 0, NIndex);
    if (Nlin)
    {
      Nstorage = 1;
    }
    else
    {
      CCTK_WARN(0, "SORFlat3D: Unable to get pointer to N");
    }
  }

  /* Shortcuts */
  dx   = GH->cctk_delta_space[0];
  dy   = GH->cctk_delta_space[1];
  dz   = GH->cctk_delta_space[2];
  
  dx2rec = 1.0/(dx*dx);
  dy2rec = 1.0/(dy*dy);
  dz2rec = 1.0/(dz*dz);

  ae = dx2rec;
  aw = dx2rec;
  an = dy2rec;
  as = dy2rec;
  at = dz2rec;
  ab = dz2rec;

  ac_orig = -2.0*dx2rec - 2.0*dy2rec - 2.0*dz2rec;

  /* Initialize omega. */
  /* TODO: Make it so that the value of the constant omega can be set. */
  if (CCTK_Equals(sor_accel, "const"))
  {
    omega = 1.8;
  }
  else
  {
    omega = 1.;
  }

  /* Set the spectral radius of the Jacobi iteration. */
  /* TODO: I think this can be easily computed for flat metrics? */
  rjacobian =  0.999;

  /* Compute whether the origin of this processor's sub grid is "red" or
     "black". */
  origin_sign = (GH->cctk_lbnd[0] + GH->cctk_lbnd[1] + GH->cctk_lbnd[2])%2;

  /* start at 1 for historic (Fortran) reasons */
  for (sorit=1; sorit<=maxit; sorit++) 
  {
    resnorm = 0.;
    
    for (k=1; k<GH->cctk_lsh[2]-1; k++) 
    {
      for (j=1; j<GH->cctk_lsh[1]-1; j++)     
      {
        for (i=1; i<GH->cctk_lsh[0]-1; i++)    
        {
          if ((origin_sign + i + j + k)%2 == sorit%2)
          {
            ac = ac_orig;

            ijk   = CCTK_GFINDEX3D(GH,i  ,j  ,k  );
            ipjk  = CCTK_GFINDEX3D(GH,i+1,j  ,k  );
            imjk  = CCTK_GFINDEX3D(GH,i-1,j  ,k  );
            ijpk  = CCTK_GFINDEX3D(GH,i  ,j+1,k  );
            ijmk  = CCTK_GFINDEX3D(GH,i  ,j-1,k  );
            ijkp  = CCTK_GFINDEX3D(GH,i  ,j  ,k+1);
            ijkm  = CCTK_GFINDEX3D(GH,i  ,j  ,k-1);
          
            if (Mstorage)
            {
              ac += Mlin[ijk];
            }

            residual = ac * var[ijk]
                + ae*var[ipjk] + aw*var[imjk]
                + an*var[ijpk] + as*var[ijmk]
                + at*var[ijkp] + ab*var[ijkm];

            if (Nstorage)
            {
              residual += Nlin[ijk];
            }

            resnorm += fabs(residual);

            var[ijk] -= omega*residual/ac; 
          }
        }
      }
    }
    
    /* reduction operation on processor-local residual values */
    ierr = CCTK_ReduceArraysGlobally(GH, -1,sum_handle, -1, 1,input_array,0,
                          input_array_dim,
                          input_array_type_codes,
                          1,
                          input_array_type_codes,
                          reduction_value);

    if (ierr<0) 
    {
      CCTK_WARN(1,"SORFlat3D: Reduction of Norm failed");
    }

#ifdef DEBUG_ELLIPTIC
    printf("it: %d SOR solve residual %f (tol %f)\n",sorit,glob_residual,tol);
#endif

    glob_residual /= (GH->cctk_gsh[0]*GH->cctk_gsh[1]*GH->cctk_gsh[2]);


#ifdef DEBUG_ELLIPTIC
    printf("it: %d SOR solve residual %f (tol %f)\n",sorit,glob_residual,tol);
#endif

    /* apply symmetry boundary conditions within loop */    
    if (CartSymVI(GH,FieldIndex)<0)
    { 
      CCTK_WARN(1,"SORFlat3D: CartSymVI failed in EllSOR loop");
      break;
    }

    /* apply Robin boundary conditions within loop */
    if (use_robin)
    {
      if (BndRobinVI(GH, sw, finf, npow,  FieldIndex)<0)
      { 
        CCTK_WARN(1,"SORFlat3D: BndRobinVI failed in EllSOR loop");
        break;
      }
    }

    /* synchronization of grid variable */
    CCTK_SyncGroupWithVarI(GH, FieldIndex);

    /* Leave iteration loop if tolerance criterium is met */
    if (glob_residual<tol)
    {
      break;
    }

    /* Update omega if using Chebyshev acceleration. */
    if (CCTK_Equals(sor_accel, "cheb"))
    {
      if (sorit == 1)
      {
        omega = 1. / (1. - 0.5 * rjacobian * rjacobian);
      }
      else
      {
         omega = 1. / (1. - 0.25 * rjacobian * rjacobian * omega);
      }
    }
    
  }

  if (elliptic_verbose)
  {
    CCTK_VInfo("EllSOR","Required SOR tolerence was set at %f",tol);
    CCTK_VInfo("EllSOR","Achieved SOR residual was %f",glob_residual);
    CCTK_VInfo("EllSOR","Number of iterations was %d (max: %d)",
               (int)sorit,(int)maxit);
  }
  
  if (glob_residual>tol) 
  {
    CCTK_WARN(1,"SORFlat3D: SOR Solver did not converge");
    retval = ELL_NOCONVERGENCE;
  }

  if (use_robin) 
  {
    free(robin);
  }
  if (sor_accel)
  { 
    free(sor_accel);
  }

  return retval;
}
