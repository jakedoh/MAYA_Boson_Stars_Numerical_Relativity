 /*@@
   @file      CartToSphere.c
   @date      Thu Apr 25 16:04:53 2002
   @author    Tom Goodale
   @desc 
   
   @enddesc
   @version $Header
 @@*/

#include "cctk.h"

#include <math.h>

#include "ADMAnalysis.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusEinstein_ADMAnalysis_CartToSphere_c)

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
   @routine    ADMAnalysis_CartToSphere
   @date       Thu Apr 25 16:17:58 2002
   @author     Tom Goodale
   @desc 
   Calculates the spherical components (r,q,p) corresponding to 
   the cartesian components (xyz) of a symmetric spatial tensor.   
   q corresponds to theta, and p to phi.
   @enddesc 
   @calls     
   @calledby   
   @history 
   @hdate Thu Apr 25 16:20:07 2002 @hauthor Tom Goodale
   @hdesc Took old metric_ and curv_cartosphere routines and
          made this generic one.
   @endhistory 
   @var     ash
   @vdesc   the local shape of the grid
   @vtype   const int *
   @vio     in
   @var     project_rdtheta_rdphi
   @vdesc   If set project angular components onto r*dtheta and 
   r*sin(theta)*dphi instead of dtheta and dphi
   @vtype   int
   @vio     in
   @var     x
   @vdesc   the x coordinate
   @vtype   const CCTK_REAL *
   @vio     in
   @var     y
   @vdesc   the y coordinate
   @vtype   const CCTK_REAL *
   @vio     in
   @var     z
   @vdesc   the z coordinate
   @vtype   const CCTK_REAL *
   @vio     in
   @var     r
   @vdesc   the r coordinate
   @vtype   const CCTK_REAL *
   @vio     in
   @var     cart_xx
   @vdesc   the xx component of the tensor
   @vtype   const CCTK_REAL *
   @vio     in
   @var     cart_xy
   @vdesc   the xy component of the tensor
   @vtype   const CCTK_REAL *
   @vio     in
   @var     cart_xz
   @vdesc   the xz component of the tensor
   @vtype   const CCTK_REAL *
   @vio     in
   @var     cart_yy
   @vdesc   the yy component of the tensor
   @vtype   const CCTK_REAL *
   @vio     in
   @var     cart_yz
   @vdesc   the yz component of the tensor
   @vtype   const CCTK_REAL *
   @vio     in
   @var     cart_zz
   @vdesc   the zz component of the tensor
   @vtype   const CCTK_REAL *
   @vio     in
   @var     sphere_rr
   @vdesc   the rr component of the tensor
   @vtype   const CCTK_REAL *
   @vio     in
   @var     sphere_rq
   @vdesc   the rq component of the tensor
   @vtype   const CCTK_REAL *
   @vio     in
   @var     sphere_rp
   @vdesc   the rp component of the tensor
   @vtype   const CCTK_REAL *
   @vio     in
   @var     sphere_qq
   @vdesc   the qq component of the tensor
   @vtype   const CCTK_REAL *
   @vio     in
   @var     sphere_qp
   @vdesc   the qp component of the tensor
   @vtype   const CCTK_REAL *
   @vio     in
   @var     sphere_pp
   @vdesc   the pp component of the tensor
   @vtype   const CCTK_REAL *
   @vio     in

 @@*/
void ADMAnalysis_CartToSphere(const CCTK_INT *ash,
                              CCTK_INT project_rdtheta_rdphi,
                              const CCTK_REAL *x,
                              const CCTK_REAL *y,
                              const CCTK_REAL *z,
                              const CCTK_REAL *r,
                              const CCTK_REAL *cart_xx,
                              const CCTK_REAL *cart_xy,
                              const CCTK_REAL *cart_xz,
                              const CCTK_REAL *cart_yy,
                              const CCTK_REAL *cart_yz,
                              const CCTK_REAL *cart_zz,
                              CCTK_REAL *sphere_rr,
                              CCTK_REAL *sphere_rq,
                              CCTK_REAL *sphere_rp,
                              CCTK_REAL *sphere_qq,
                              CCTK_REAL *sphere_qp,
                              CCTK_REAL *sphere_pp)
{
  int i;
  CCTK_REAL cost;
  CCTK_REAL sint;
  CCTK_REAL sinp;
  CCTK_REAL cosp;
  CCTK_REAL rvalue;
  CCTK_REAL sxy;

  CCTK_REAL    txx,txy,txz,tyy,tyz,tzz;
    
  /* loop over all the grid */
  for(i = 0; i < ash[0]*ash[1]*ash[2]; i++)  
  {
    
     txx = cart_xx[i];
     txy = cart_xy[i];
     txz = cart_xz[i];
     tyy = cart_yy[i];
     tyz = cart_yz[i];
     tzz = cart_zz[i];
     rvalue  = r[i];
     sxy = sqrt( SQR(x[i]) + SQR(y[i]));
       
  /* be careful with r=0 and xy plane */
     if (rvalue==0.0) 
     {
       cost = 1.0;
       sint = 0.0; 
       sinp = 0.0;
       cosp = 1.0;
     } 
     else if (sxy==0) 
     {
       cost = 1.0; 
       sint = 0.0;
       sinp = 0.0;
       cosp = 1.0;
     } 
     else 
     {
       cost = z[i]/rvalue;
       sint = sxy/rvalue;
       sinp = y[i]/sxy;
       cosp = x[i]/sxy;
     }
     
     sphere_rr[i]= 
         tyy*SQR(sinp)*SQR(sint)+
         2*cosp*txy*sinp*SQR(sint)+
         SQR(cosp)*txx*SQR(sint)+
         2*cost*tyz*sinp*sint+
         2*cosp*cost*txz*sint+
         SQR(cost)*tzz;
               
     sphere_qq[i] = 
         (tzz*SQR(sint)+
          (-2*cost*tyz*sinp-
           2*cosp*cost*txz)*sint+
          SQR(cost)*tyy*SQR(sinp)+
          2*cosp*SQR(cost)*txy*sinp
          +SQR(cosp)*SQR(cost)*txx);
     
     if (!project_rdtheta_rdphi) 
     {
       sphere_qq[i] *= SQR(r[i]);
     }
     
     sphere_pp[i] = 
         (txx*SQR(sinp)-
          2*cosp*txy*sinp+
          SQR(cosp)*tyy);
     
     if (!project_rdtheta_rdphi)  
     {
       sphere_pp[i] *= SQR(r[i]) * SQR(sint);
     }
               
     sphere_rq[i] = 
         (cost*tyy*SQR(sinp)*sint+
          2*cosp*cost*txy*sinp*sint-
          cost*tzz*sint+
          SQR(cosp)*cost*txx*sint+
          2*SQR(cost)*tyz*sinp-
          tyz*sinp+
          2*cosp*SQR(cost)*txz-
          cosp*txz);
     
     if (!project_rdtheta_rdphi) 
     {
       sphere_rq[i] *= r[i];
     }

     sphere_rp[i] = 
         ((-txy*SQR(sinp)+
           (cosp*tyy-cosp*txx)*sinp+
           SQR(cosp)*txy)*sint-
          cost*txz*sinp+cosp*cost*tyz);
     
     if (!project_rdtheta_rdphi)  
     {
       sphere_rp[i] *= r[i] * sint;
     }
     
     sphere_qp[i] =
         ((txz*sinp-cosp*tyz)*sint+
          cost*(-txy*SQR(sinp)+
                cosp*(tyy-txx)*sinp+SQR(cosp)*txy));

     if (!project_rdtheta_rdphi) 
     {
       sphere_qp[i] *= SQR(r[i]) * sint;
     }
  }
}
/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

