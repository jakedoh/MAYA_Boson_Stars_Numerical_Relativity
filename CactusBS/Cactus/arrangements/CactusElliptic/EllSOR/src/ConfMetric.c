 /*@@
   @file      ConfMetric.c
   @date      Tue Sep 26 11:29:18 2000
   @author    Gerd Lanfermann
   @desc 
     Provides sor solver engine routines
   @enddesc 
 @@*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "Boundary.h"
#include "Symmetry.h"
#include "Ell_DBstructure.h"
#include "EllBase.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusElliptic_EllSOR_ConfMetric_c)

#define SQR(a) ((a)*(a))

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

int SORConfMetric3D(cGH *GH, int *MetricPsiI, int conformal,
                    int FieldIndex, int MIndex, int NIndex,
                    CCTK_REAL *AbsTol, CCTK_REAL *RelTol);

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/* static int firstcall = 1; */

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/


 /*@@
   @routine    SORConfMetric3D
   @date       Tue Sep 26 11:28:08 2000
   @author     Joan Masso, Paul Walker, Gerd Lanfermann
   @desc 
   This is a standalone sor solver engine, 
   called by the wrapper functions
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

int SORConfMetric3D(cGH *GH, int *MetricPsiI, int conformal,
                    int FieldIndex, int MIndex, int NIndex,
                    CCTK_REAL *AbsTol, CCTK_REAL *RelTol)
{
  DECLARE_CCTK_PARAMETERS  

  int retval = ELL_SUCCESS;
  int timelevel;
  int origin_sign;

  /* The pointer to the metric fields */
  CCTK_REAL *gxx =NULL, *gxy =NULL; 
  CCTK_REAL *gxz =NULL, *gyy =NULL; 
  CCTK_REAL *gyz =NULL, *gzz =NULL; 
  CCTK_REAL *psi =NULL;
  CCTK_REAL *Mlin=NULL, *Nlin=NULL;   
  CCTK_REAL *var =NULL; 

  /* The inverse metric, allocated here */
  CCTK_REAL *uxx=NULL, *uyy=NULL, 
            *uzz=NULL, *uxz=NULL,
            *uxy=NULL, *uyz=NULL;

  /* shortcuts for metric, psi, deltas, etc. */
  CCTK_REAL pm4, p12, det;

  CCTK_REAL dx,dy,dz;
  CCTK_REAL dxdx, dydy, dzdz, 
            dxdy, dxdz, dydz;

  /* Some physical variables */
  CCTK_REAL omega, resnorm=0.0, residual=0.0;
  CCTK_REAL glob_residual=0.0, rjacobian=0.0; 
  CCTK_REAL finf;
  CCTK_INT npow;
  CCTK_REAL tol;

  /* Iteration / stepping  variables */
  int sorit;
  CCTK_INT maxit; 
  int i,j,k;
  int nxyz;
  
  /* 19 point stencil index */
  int ijk;
  int ipjk, ijpk, ijkp, imjk, ijmk, ijkm;
  int ipjpk, ipjmk, imjpk, imjmk;
  int ipjkp, ipjkm, imjkp, imjkm;
  int ijpkp, ijpkm, ijmkp, ijmkm;
 
  /* 19 point stencil coefficients */
  CCTK_REAL ac;
  CCTK_REAL ae,aw,an,as,at,ab;
  CCTK_REAL ane, anw, ase, asw, ate, atw, abe, abw;
  CCTK_REAL atn, ats, abn, absol;

  /* Miscellaneous */
  int sum_handle=-1;
  int sw[3], ierr;
  int Mstorage=0, Nstorage=0;
  size_t varsize;
  CCTK_REAL detrec_pm4, sqrtdet;
  char *robin=NULL;
  char *sor_accel=NULL;
  const void* input_array[1];
  void*  reduction_value[1];
  CCTK_INT input_array_dim[1];
  CCTK_INT input_array_type_codes[1];

  input_array[0]             = &resnorm;
  reduction_value[0]         = &glob_residual;
  input_array_type_codes[0]  = CCTK_VARIABLE_REAL;
  
  /* Get the reduction handle */
  sum_handle = CCTK_LocalArrayReductionHandle("sum");
  if (sum_handle<0) 
  {
    CCTK_WARN(1,"SORConfMetric3D: "
              "Cannot get handle for sum reduction");
    retval = ELL_FAILURE;
  }

  /* If Robin BCs are set, prepare for a boundary call:
     setup stencil width and get Robin constants (set by the routine
     which is calling the solver interface) */
  if (conformal)
  {
    if (Ell_GetStrKey(&robin, "EllLinConfMetric::Bnd::Robin")< 0)
    {
      CCTK_WARN(2,"SORConfMetric3D: Robin keys not set for LinConfMetric solver");
    }
  }
  else
  {
    if (Ell_GetStrKey(&robin, "EllLinMetric::Bnd::Robin")< 0)
    {
      CCTK_WARN(2,"SORConfMetric3D: Robin keys not set for LinMetric solver");
    }
  }
  
  if (robin && CCTK_EQUALS(robin,"yes")) 
  { 
    sw[0]=1; 
    sw[1]=1; 
    sw[2]=1;
    
    if (conformal)
    {
      ierr = Ell_GetRealKey(&finf, "EllLinConfMetric::Bnd::Robin::inf");
      if (ierr < 0)
      {
        CCTK_WARN(1,"EllLinConfMetric::Bnd::Robin::inf is not set");
      }
      ierr = Ell_GetIntKey (&npow, "EllLinConfMetric::Bnd::Robin::falloff");
      if (ierr < 0)
      {
        CCTK_WARN(1,"EllLinConfMetric::Bnd::Robin::falloff is not set");
      }
    }
    else
    {
      ierr = Ell_GetRealKey(&finf, "EllLinMetric::Bnd::Robin::inf");
      if (ierr < 0)
      {
        CCTK_WARN(1,"EllLinMetric::Bnd::Robin::inf is not set");
      }
      ierr = Ell_GetIntKey (&npow, "EllLinMetric::Bnd::Robin::falloff");
      if (ierr < 0)
      {
        CCTK_WARN(1,"EllLinMetric::Bnd::Robin::falloff is not set");
      }
    }
  }

  /* Get the maximum number of iterations allowed. */
  if (Ell_GetIntKey(&maxit, "Ell::SORmaxit") == ELLGET_NOTSET)
  {
    CCTK_WARN(1,"SORConfMetric3D: Ell::SORmaxit not set. "
              "Maximum allowed iterations being set to 100.");
    maxit = 100;
  }

  /* Only supports absolute tolerance */
  tol   = AbsTol[0];

  /* Get the type of acceleration to use. */
  if (Ell_GetStrKey(&sor_accel, "Ell::SORaccel") == ELLGET_NOTSET)
  {
    const char tmpstr[6] = "const";
    CCTK_WARN(3,"SORConfMetric3D: Ell::SORaccel not set. "
              "Omega being set to a constant value of 1.8.");
    sor_accel = strdup(tmpstr);
  }

  if (CCTK_Equals(elliptic_verbose, "yes"))
  {
    if (conformal)
    {
      CCTK_VInfo(CCTK_THORNSTRING,"SOR solver for linear conformal metric (to tolerance %f)",tol);
    }
    else
    {
      CCTK_VInfo(CCTK_THORNSTRING,"SOR solver for linear metric (to tolerance %f)",tol);
    }
  }

  /* Things to do only once! 
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
        CCTK_WARN(0, "SORConfMetric3D: sor_accel not set");
      }
    }
    firstcall = 0;
    }*/

  /* 
     Get the data ptr of these GFs.
     They all have to be on the same timelevel. 
     Derive the metric data pointer from the index array. 
     Note the ordering in the metric   
  */

  timelevel = 0;
  var = (CCTK_REAL*) CCTK_VarDataPtrI(GH, timelevel, FieldIndex);
  if (!var)
  {
    CCTK_WARN(0,"Invalid data for Field");
  }

  timelevel = 0;
  gxx = (CCTK_REAL*) CCTK_VarDataPtrI(GH, timelevel, MetricPsiI[0]);
  gxy = (CCTK_REAL*) CCTK_VarDataPtrI(GH, timelevel, MetricPsiI[1]);
  gxz = (CCTK_REAL*) CCTK_VarDataPtrI(GH, timelevel, MetricPsiI[2]);
  gyy = (CCTK_REAL*) CCTK_VarDataPtrI(GH, timelevel, MetricPsiI[3]);
  gyz = (CCTK_REAL*) CCTK_VarDataPtrI(GH, timelevel, MetricPsiI[4]);
  gzz = (CCTK_REAL*) CCTK_VarDataPtrI(GH, timelevel, MetricPsiI[5]);
  if (!(gxx&&gyy&&gzz&&gxy&&gxz&&gyz))
  {
    CCTK_WARN(0,"Invalid data for Metric");
  }

  if (conformal)
  {
    timelevel = 0;
    psi = NULL;
    psi = (CCTK_REAL*) CCTK_VarDataPtrI(GH, timelevel, MetricPsiI[6]);
    if (!psi)
    {
      CCTK_WARN(0,"Invalid data for conformal factor");
    }

  }

  /* if we have a negative index for M/N, this GF is not needed, 
     there for don't even look for it. when index positive,
     set flag Mstorage=1, dito for N */
  if (MIndex>=0)  
  { 
    timelevel = 0;
    Mlin = (CCTK_REAL*) CCTK_VarDataPtrI(GH,timelevel,MIndex);
    if (Mlin)
    {
      Mstorage = 1;
    }
    else
    {
      CCTK_WARN(0, "SORConfMetric3D: Unable to get pointer to M.");
    }
  }
  if (NIndex>=0) 
  {
    timelevel = 0;
    Nlin = (CCTK_REAL*) CCTK_VarDataPtrI(GH,timelevel,NIndex);
    if (Nlin)
    {
      Nstorage = 1;
    }
    else
    {
      CCTK_WARN(0, "SORConfMetric3D: Unable to get pointer to N.");
    }
  }

  /* Shortcuts */
  dx   = GH->cctk_delta_space[0];
  dy   = GH->cctk_delta_space[1];
  dz   = GH->cctk_delta_space[2];
  nxyz = GH->cctk_lsh[0]*GH->cctk_lsh[1]*GH->cctk_lsh[2];

  /* Allocate the inverse metric */
  varsize = (size_t)CCTK_VarTypeSize(CCTK_VarTypeI(FieldIndex))*nxyz;
  uxx = (CCTK_REAL*) malloc(varsize);
  uxy = (CCTK_REAL*) malloc(varsize);
  uxz = (CCTK_REAL*) malloc(varsize);
  uyy = (CCTK_REAL*) malloc(varsize);
  uyz = (CCTK_REAL*) malloc(varsize);
  uzz = (CCTK_REAL*) malloc(varsize);

  if (!uxx || !uxy || !uxz || !uyy || !uyz || !uzz)
  {
    CCTK_WARN(0,"SORConfMetric3D: Memory allocation failed for upper metric");
  }
  
  /* calculate the differential coefficient */
  dxdx = 0.5/(dx*dx);
  dydy = 0.5/(dy*dy);
  dzdz = 0.5/(dz*dz);
  dxdy = 0.25/(dx*dy);
  dxdz = 0.25/(dx*dz);
  dydz = 0.25/(dy*dz);

  /* Calculate the inverse metric */
  for (ijk=0; ijk<nxyz; ijk++) 
  {
    det = -(SQR(gxz[ijk])*gyy[ijk]) + 
      2*gxy[ijk]*gxz[ijk]*gyz[ijk] - 
      gxx[ijk]*SQR(gyz[ijk])  -
      SQR(gxy[ijk])*gzz[ijk] + 
      gxx[ijk]*gyy[ijk]*gzz[ijk];
    
    if (conformal) 
    {
      pm4 = 1.0/pow(psi[ijk],4.0);
      p12 = pow(psi[ijk],12.0);
    } 
    else 
    {
      pm4 = 1.0;
      p12 = 1.0;
    }

    /* try to avoid constly division */
    detrec_pm4 = 1.0/det*pm4;
    sqrtdet    = sqrt(det);

    uxx[ijk]=(-SQR(gyz[ijk])     + gyy[ijk]*gzz[ijk])*detrec_pm4;
    uxy[ijk]=( gxz[ijk]*gyz[ijk] - gxy[ijk]*gzz[ijk])*detrec_pm4;
    uxz[ijk]=(-gxz[ijk]*gyy[ijk] + gxy[ijk]*gyz[ijk])*detrec_pm4;
    uyy[ijk]=(-SQR(gxz[ijk])     + gxx[ijk]*gzz[ijk])*detrec_pm4;
    uyz[ijk]=( gxy[ijk]*gxz[ijk] - gxx[ijk]*gyz[ijk])*detrec_pm4;
    uzz[ijk]=(-SQR(gxy[ijk])     + gxx[ijk]*gyy[ijk])*detrec_pm4;
        
    det    = det*p12;
    sqrtdet= sqrt(det);
   
    /* Rescaling for the uppermetric solver */
    if (Mstorage) 
    {
      Mlin[ijk] = Mlin[ijk]*sqrt(det);
    }
    if (Nstorage)
    {
      Nlin[ijk] = Nlin[ijk]*sqrt(det);
    }
    
    uxx[ijk]=uxx[ijk]*dxdx*sqrtdet;
    uyy[ijk]=uyy[ijk]*dydy*sqrtdet;
    uzz[ijk]=uzz[ijk]*dzdz*sqrtdet;
    uxy[ijk]=uxy[ijk]*dxdy*sqrtdet;
    uxz[ijk]=uxz[ijk]*dxdz*sqrtdet;
    uyz[ijk]=uyz[ijk]*dydz*sqrtdet;

  }     

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
  rjacobian =  0.999;

  /* Compute whether the origin of this processor's 
     sub grid is "red" or "black". */

  origin_sign = (GH->cctk_lbnd[0] + GH->cctk_lbnd[1] + GH->cctk_lbnd[2])%2;

  /* start at 1 for historic (Fortran) reasons */
  for (sorit=1; sorit<=maxit; sorit++)
  {
    resnorm = 0.0;

    for (k=1; k<GH->cctk_lsh[2]-1; k++) 
    {
      for (j=1; j<GH->cctk_lsh[1]-1; j++)     
      {
        for (i=1; i<GH->cctk_lsh[0]-1; i++)    
        {
          if ((origin_sign+i+j+k)%2 == sorit%2)
          {
            ijk   = CCTK_GFINDEX3D(GH,i  ,j  ,k  );

            ipjk  = CCTK_GFINDEX3D(GH,i+1,j  ,k  );
            imjk  = CCTK_GFINDEX3D(GH,i-1,j  ,k  );
            ijpk  = CCTK_GFINDEX3D(GH,i  ,j+1,k  );
            ijmk  = CCTK_GFINDEX3D(GH,i  ,j-1,k  );
            ijkp  = CCTK_GFINDEX3D(GH,i  ,j  ,k+1);
            ijkm  = CCTK_GFINDEX3D(GH,i  ,j  ,k-1);

            ipjpk = CCTK_GFINDEX3D(GH,i+1,j+1,k  );
            ipjmk = CCTK_GFINDEX3D(GH,i+1,j-1,k  );
            imjpk = CCTK_GFINDEX3D(GH,i-1,j+1,k  );
            imjmk = CCTK_GFINDEX3D(GH,i-1,j-1,k  );
          
            ijpkp = CCTK_GFINDEX3D(GH,i  ,j+1,k+1);
            ijpkm = CCTK_GFINDEX3D(GH,i  ,j+1,k-1);
            ijmkp = CCTK_GFINDEX3D(GH,i  ,j-1,k+1);
            ijmkm = CCTK_GFINDEX3D(GH,i  ,j-1,k-1);
          
            ipjkp = CCTK_GFINDEX3D(GH,i+1,j  ,k+1);
            ipjkm = CCTK_GFINDEX3D(GH,i+1,j  ,k-1);
            imjkp = CCTK_GFINDEX3D(GH,i-1,j  ,k+1);
            imjkm = CCTK_GFINDEX3D(GH,i-1,j  ,k-1);

            ac = -1.0*uxx[ipjk] - 2.0*uxx[ijk] - 1.0*uxx[imjk]
                 -1.0*uyy[ijpk] - 2.0*uyy[ijk] - 1.0*uyy[ijmk]
                 -1.0*uzz[ijkp] - 2.0*uzz[ijk] - 1.0*uzz[ijkm];
          
            if (Mstorage) 
            {
              ac += Mlin[ijk];
            }

            ae  = uxx[ipjk]+uxx[ijk];
            aw  = uxx[imjk]+uxx[ijk];
            an  = uyy[ijpk]+uyy[ijk];
            as  = uyy[ijmk]+uyy[ijk];
            at  = uzz[ijkp]+uzz[ijk];
            ab  = uzz[ijkm]+uzz[ijk];
          
            ane = uxy[ijpk]+uxy[ipjk];
            anw =-uxy[imjk]-uxy[ijpk];
            ase =-uxy[ipjk]-uxy[ijmk];
            asw = uxy[imjk]+uxy[ijmk];
          
            ate = uxz[ijkp]+uxz[ipjk];
            atw =-uxz[imjk]-uxz[ijkp];
            abe =-uxz[ipjk]-uxz[ijkm];
            abw = uxz[imjk]+uxz[ijkm];
          
            atn = uyz[ijpk]+uyz[ijkp];               
            ats =-uyz[ijkp]-uyz[ijmk];
            abn =-uyz[ijkm]-uyz[ijpk];
            absol = uyz[ijkm]+uyz[ijmk];
          
            residual = ac * var[ijk]
                + ae *var[ipjk]  +  aw*var[imjk]
                + an *var[ijpk]  +  as*var[ijmk]
                + at *var[ijkp]  +  ab*var[ijkm];
          
            residual = residual 
                + ane*var[ipjpk] + anw*var[imjpk]; 
          
            residual = residual 
                + ase*var[ipjmk] + asw*var[imjmk];

            residual = residual    
                + ate*var[ipjkp] + atw*var[imjkp] 
                + abe*var[ipjkm] + abw*var[imjkm]
                + atn*var[ijpkp] + ats*var[ijmkp]
                + abn*var[ijpkm] + absol*var[ijmkm];

            if (Nstorage)
            {
              residual +=Nlin[ijk];
            }

            resnorm  = resnorm + fabs(residual);

            var[ijk] = var[ijk] - omega*residual/ac; 

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
      CCTK_WARN(1,"SORConfMetric3D: Reduction of residual  failed");
    }

    glob_residual = glob_residual /  
      (GH->cctk_gsh[0]*GH->cctk_gsh[1]*GH->cctk_gsh[2]);

    /* apply symmetry boundary conditions within loop */    
    if (CartSymVI(GH,FieldIndex)<0) 
    {
      CCTK_WARN(1,"SORConfMetric3D: CartSymVI failed in EllSOR loop");
    }

    /* apply Robin boundary conditions within loop */
    if (robin && CCTK_EQUALS(robin,"yes")) 
    {
      if (BndRobinVI(GH, sw, finf, npow,  FieldIndex)<0)
      {
        CCTK_WARN(1, "SORConfMetric3D: BndRobinVI failed in EllSOR loop");
        break;
      }
    }

    /* synchronization of grid variable */
    ierr = CCTK_SyncGroupWithVarI(GH, FieldIndex);
    if (ierr < 0)
    {
      CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
                 "ConfMetric3D: Synchronization of %s failed (Error: %d)",
                 CCTK_VarName(FieldIndex),ierr);
    }

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

    if (CCTK_Equals(elliptic_verbose,"debug"))
    {
      CCTK_VInfo(CCTK_THORNSTRING,"  .. residual at %f",glob_residual);
    }
  }

  if (CCTK_Equals(elliptic_verbose,"yes"))
  {
    CCTK_VInfo(CCTK_THORNSTRING,"  ... achieved SOR residual %f (at %d iterations)",glob_residual,sorit);
  }

  if (glob_residual>tol) 
  {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING, 
               "SOR Solver did not converge (reached tol:%f at %d iterations, requested tol:%f",glob_residual,sorit,tol);
    retval = ELL_NOCONVERGENCE;
  }

  if (uxx) 
  {
    free(uxx); 
  }
  if (uyy) 
  {
    free(uyy); 
  }
  if (uzz)
  { 
    free(uzz);
  }
  if (uxy) 
  {
    free(uxy); 
  }
  if (uxz) 
  {
    free(uxz); 
  }
  if (uyz) 
  {
    free(uyz);
  }
  if (robin) 
  {
    free(robin);
  }
  if (sor_accel) 
  {
    free(sor_accel);
  }

  return retval;
}
