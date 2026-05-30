 /*@@
   @header    ADMAnalysis.h
   @date      Thu Apr 25 16:36:20 2002
   @author    Tom Goodale
   @desc 
   
   @enddesc
   @version $Header$
 @@*/

#ifndef _ADMAnalysis_H_
#define _ADMAnalysis_H_ 1

#ifdef __cplusplus
extern "C" 
{
#endif

void ADMAnalysis_CartToSphere(const CCTK_INT *ash,
                              CCTK_INT r2norm,
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
                              CCTK_REAL *sphere_pp);

void ADMAnalysis_Trace(const CCTK_INT *ash,
                       const CCTK_REAL *g11,
                       const CCTK_REAL *g12,
                       const CCTK_REAL *g13,
                       const CCTK_REAL *g22,
                       const CCTK_REAL *g23,
                       const CCTK_REAL *g33,
                       CCTK_REAL *tensor11,
                       CCTK_REAL *tensor12,
                       CCTK_REAL *tensor13,
                       CCTK_REAL *tensor22,
                       CCTK_REAL *tensor23,
                       CCTK_REAL *tensor33,
                       CCTK_REAL *trace,
                       CCTK_REAL *detg);

#ifdef __cplusplus
}
#endif

#endif /* _ADMAnalysis_H_ */
