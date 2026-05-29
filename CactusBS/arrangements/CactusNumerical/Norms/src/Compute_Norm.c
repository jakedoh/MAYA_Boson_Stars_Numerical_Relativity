#include <math.h>
#include <stdio.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"


static const char * rcsid = "$Header$";

CCTK_FILEVERSION(Norms_Compute_Norm_c);


void Norms_Compute_Norms (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int vind; /* variable index */
  CCTK_REAL *data;
  int index,indexp, i,j,k;
  int ierr;
  int istart, iend, jstart, jend, kstart, kend;
  CCTK_REAL dx,dy,dz;
  int target_proc,reduction_handle;
  CCTK_REAL result;

  if (verbose>1)
    CCTK_INFO("computing norms");

  *norm=0.;
  result=0.;

  dx = CCTK_DELTA_SPACE(0);
  dy = CCTK_DELTA_SPACE(1);
  dz = CCTK_DELTA_SPACE(2);

  istart = cctk_nghostzones[0];
  jstart = cctk_nghostzones[1];
  kstart = cctk_nghostzones[2];

  iend = cctk_lsh[0] - cctk_nghostzones[0];
  jend = cctk_lsh[1] - cctk_nghostzones[1];
  kend = cctk_lsh[2] - cctk_nghostzones[2];

  /* initialise boundaries */
  for (k=0;k<cctk_lsh[2];k++) {
    for (j=0;j<cctk_lsh[1];j++) {
      for (i=0;i<cctk_lsh[0];i++) {
        index=CCTK_GFINDEX3D(cctkGH,i,j,k);
        diff_term[index]=0.;
      }
    }
  }

  target_proc = 0;
  reduction_handle = CCTK_ReductionArrayHandle("norm2");
  if (reduction_handle < 0)
       CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                        "can't get sum-reduction handle!");

  if (verbose>2)
    fprintf(stderr,"   nr1stvars %d nr2ndvars %d\n",
            (int)*nr1stvars,(int)*nr2ndvars);

  /* L2-norm of 1st order variables */
  *norm_u=0.;
  for (vind=0; vind< *nr1stvars; vind++) {
    ierr=CCTK_Reduce(cctkGH,target_proc,reduction_handle,1,
                  CCTK_VARIABLE_REAL,&result,1,varindices_1st[vind]);
    if (ierr<0)
        CCTK_WARN(1,"reduction failed");
    if (verbose>4) {
      fprintf(stderr,"     vind=%d vi=%d var=%s result=%1.8g\n",
              vind,(int)varindices_1st[vind],
              CCTK_VarName(varindices_1st[vind]),(double)result);
    }
    result=result-bgvals_1st[vind];
    *norm_u=*norm_u+result*result;
    *norm=*norm+result*result;
  }

  if (verbose>2)
    fprintf(stderr,"  after 1st order vars norm=%16.8g\n",(double)*norm);

  /* L2-norm of 2nd order variables */
  *norm_v=0.;
  for (vind=0; vind< *nr2ndvars; vind++) {
    ierr=CCTK_Reduce(cctkGH,target_proc,reduction_handle,1,
                  CCTK_VARIABLE_REAL,&result,1,varindices_2nd[vind]);
    if (ierr<0)
      CCTK_WARN(1,"reduction failed");
    if (verbose>4) {
      fprintf(stderr,"     vind=%d vi=%d var=%s result=%1.8g\n",
              vind,(int)varindices_2nd[vind],
              CCTK_VarName(varindices_2nd[vind]),(double)result);
    }
    result=result-bgvals_2nd[vind];
    *norm_v=*norm_v+result*result;
    *norm=*norm+result*result;
  }

  if (verbose>2)
    fprintf(stderr,"  after 2nd order vars norm=%16.8g\n",(double)*norm);

  /* L2-norm of D_{+i} terms */
  *norm_grad=0.;
  /* d_x */
  for (vind=0; vind< *nr2ndvars; vind++) {
    data=CCTK_VarDataPtrI(cctkGH,0,varindices_2nd[vind]);
    for (k=kstart;k<kend;k++) {
      for (j=jstart;j<jend;j++) {
        for (i=istart;i<iend;i++) {
          index =CCTK_GFINDEX3D(cctkGH,i,j,k);
          indexp=CCTK_GFINDEX3D(cctkGH,i+1,j,k);
          diff_term[index]=(data[indexp]-data[index])/dx;
        }
      }
    }
    ierr=CCTK_Reduce(cctkGH,target_proc,reduction_handle,1,
                  CCTK_VARIABLE_REAL,&result,1,
                  CCTK_VarIndex("norms::diff_term"));
    if (ierr<0)
      CCTK_WARN(1,"reduction failed");
    if (verbose>4) {
      fprintf(stderr,"     vind=%d vi=%d var=%s result=%1.8g\n",
              vind,(int)varindices_2nd[vind],
              CCTK_VarName(varindices_2nd[vind]),(double)result);
    }
    *norm_grad=*norm_grad+result*result;
    *norm=*norm+result*result;
  }

  if (verbose>2)
    fprintf(stderr,"  after dx 2nd order terms norm=%16.8g\n",(double)*norm);

  /* d_y */
  for (vind=0; vind< *nr2ndvars; vind++) {
    data=CCTK_VarDataPtrI(cctkGH,0,varindices_2nd[vind]);
    for (k=kstart;k<kend;k++) {
      for (j=jstart;j<jend;j++) {
        for (i=istart;i<iend;i++) {
          index =CCTK_GFINDEX3D(cctkGH,i,j,k);
          indexp=CCTK_GFINDEX3D(cctkGH,i,j+1,k);
          diff_term[index]=(data[indexp]-data[index])/dy;
        }
      }
    }
    ierr=CCTK_Reduce(cctkGH,target_proc,reduction_handle,1,
                  CCTK_VARIABLE_REAL,&result,1,
                  CCTK_VarIndex("norms::diff_term"));
    if (ierr<0)
      CCTK_WARN(1,"reduction failed");
    if (verbose>4) {
      fprintf(stderr,"     vind=%d vi=%d var=%s result=%1.8g\n",
              vind,(int)varindices_2nd[vind],
              CCTK_VarName(varindices_2nd[vind]),(double)result);
    }
    *norm_grad=*norm_grad+result*result;
    *norm=*norm+result*result;
  }

  if (verbose>2)
    fprintf(stderr,"  after dy 2nd order terms norm=%16.8g\n",(double)*norm);

  /* d_z */
  for (vind=0; vind< *nr2ndvars; vind++) {
    data=CCTK_VarDataPtrI(cctkGH,0,varindices_2nd[vind]);
    for (k=kstart;k<kend;k++) {
      for (j=jstart;j<jend;j++) {
        for (i=istart;i<iend;i++) {
          index =CCTK_GFINDEX3D(cctkGH,i,j,k);
          indexp=CCTK_GFINDEX3D(cctkGH,i,j,k+1);
          diff_term[index]=(data[indexp]-data[index])/dz;
        }
      }
    }
    ierr=CCTK_Reduce(cctkGH,target_proc,reduction_handle,1,
                  CCTK_VARIABLE_REAL,&result,1,
                  CCTK_VarIndex("norms::diff_term"));
    if (ierr<0)
      CCTK_WARN(1,"reduction failed");
    if (verbose>4) {
      fprintf(stderr,"     vind=%d vi=%d var=%s result=%1.8g\n",
              vind,(int)varindices_2nd[vind],
              CCTK_VarName(varindices_2nd[vind]),(double)result);
    }
    *norm_grad=*norm_grad+result*result;
    *norm=*norm+result*result;
  }

  if (verbose>2)
    fprintf(stderr,"  after dz 2nd order terms norm=%16.8g\n",(double)*norm);

  if (verbose>2)
    fprintf(stderr,"  after 2nd order vars (incl. D_{+i}): norm=%16.8g\n",
            (double)*norm);
  if (verbose>2)
    fprintf(stderr,"final norm^2=%16.8g on this CPU\n",(double)*norm);

  *norm_u=sqrt(*norm_u);
  *norm_v=sqrt(*norm_v);
  *norm_grad=sqrt(*norm_grad);
  *norm=sqrt(*norm);

  if (verbose>1) {
    CCTK_VInfo(CCTK_THORNSTRING,
               "norm= %16.8g (norm_u=%f, norm_v=%f, norm_grad=%f)",
               (double)*norm,(double)*norm_u,(double)*norm_v,(double)*norm_grad);
  }
}
