 /*@@
   @file      Analysis.c
   @date      Thu Apr 25 18:46:27 2002
   @author    Tom Goodale
   @desc 
   Routines to do the various analysis tasks.
   @enddesc
   @version $Header$
 @@*/

#include "cctk.h"

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "ADMAnalysis.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusEinstein_ADMAnalysis_Analysis_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

#define SQR(a) ((a)*(a))

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void ADMAnalysis_EvaltrK(CCTK_ARGUMENTS);

void ADMAnalysis_MetricCartToSphere(CCTK_ARGUMENTS);

void ADMAnalysis_CurvCartToSphere(CCTK_ARGUMENTS);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    ADMAnalysis_EvaltrK
   @date       Thu Apr 25 19:11:32 2002
   @author     Tom Goodale
   @desc 
   Scheduled routine to evaluate the trace of the extrinsic curvature
   and the determinant of the metric.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
@@*/
void ADMAnalysis_EvaltrK(CCTK_ARGUMENTS)
{
  int i, d;
  CCTK_INT ash[3];
  DECLARE_CCTK_ARGUMENTS_ADMAnalysis_EvaltrK;
  DECLARE_CCTK_PARAMETERS;

  for (d=0; d<3; ++d) ash[d] = cctk_ash[d];
  ADMAnalysis_Trace(ash,
                    gxx,
                    gxy,
                    gxz,
                    gyy,
                    gyz,
                    gzz,
                    kxx,
                    kxy,
                    kxz,
                    kyy,
                    kyz,
                    kzz,
                    trK,
                    detg);

  /* If we have a conformal metric, must convert trK to non-conformal form */
  if(CCTK_EQUALS(metric_type, "static conformal") && *conformal_state > 0)
  {
    for(i = 0; i< cctk_ash[0]*cctk_ash[1]*cctk_ash[2]; i++)  
    {
      trK[i] = trK[i] /(SQR(psi[i])*SQR(psi[i]));
    }
  }
}

/*@@
  @routine    ADMAnalysis_MetricCartToSphere
  @date       Thu Apr 25 19:11:32 2002
  @author     Tom Goodale
  @desc 
  Scheduled routine to evaluate the components of the metric
  in spherical coordinates
  @enddesc 
  @calls     
  @calledby   
  @history 
  
  @endhistory 
  
@@*/
void ADMAnalysis_MetricCartToSphere(CCTK_ARGUMENTS)
{
  int d;
  CCTK_INT ash[3];
  DECLARE_CCTK_ARGUMENTS_ADMAnalysis_MetricCartToSphere;
  DECLARE_CCTK_PARAMETERS;

  for (d=0; d<3; ++d) ash[d] = cctk_ash[d];
  ADMAnalysis_CartToSphere(ash,
                           normalize_dtheta_dphi,
                           x,
                           y,
                           z,
                           r,
                           gxx,
                           gxy,
                           gxz,
                           gyy,
                           gyz,
                           gzz,
                           grr,
                           grq,
                           grp,
                           gqq,
                           gqp,
                           gpp);
}

/*@@
  @routine    ADMAnalysis_MetricCartToSphere
  @date       Thu Apr 25 19:11:32 2002
  @author     Tom Goodale
  @desc 
  Scheduled routine to evaluate the components of the extrinsic curvature
  in spherical coordinates
  @enddesc 
  @calls     
  @calledby   
  @history 
  
  @endhistory 
  
@@*/
void ADMAnalysis_CurvCartToSphere(CCTK_ARGUMENTS)
{
  int d;
  CCTK_INT ash[3];
  DECLARE_CCTK_ARGUMENTS_ADMAnalysis_CurvCartToSphere;
  DECLARE_CCTK_PARAMETERS;

  for (d=0; d<3; ++d) ash[d] = cctk_ash[d];
  ADMAnalysis_CartToSphere(ash,
                           normalize_dtheta_dphi,
                           x,
                           y,
                           z,
                           r,
                           kxx,
                           kxy,
                           kxz,
                           kyy,
                           kyz,
                           kzz,
                           krr,
                           krq,
                           krp,
                           kqq,
                           kqp,
                           kpp);
}
/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

