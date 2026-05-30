 /*@@
   @header    EllBase.h
   @date      
   @author    Gerd Lanferman
   @desc 
   Basic Elliptic solver functions.   
   @enddesc 
   @version $Header$
 @@*/

#ifndef _ELLBASE_H_
#define _ELLBASE_H_

#define ELL_SUCCESS       0
#define ELL_NOSOLVER      -1
#define ELL_NOCONVERGENCE -2
#define ELL_NOCLASS       -3
#define ELL_SOLVEREXISTS  -4
#define ELL_CLASSEXISTS   -5
#define ELL_FAILURE       -6

#ifdef CCODE

/* Argumennt structure for the four different types of elliptic solvers
   provided at this point. Difference is MetricI, MetricPsiI and StencilGFI, 
   which are arrays holding the grid function indices of the metric, 
   metric+psi or the 27 stencil grid functions, respectively */

#define LINELL_FLAT3D_ARGS \
          cGH *GH, \
          CCTK_REAL tolerance, \
          int FieldIndex, \
          int MIndex, \
          int NIndex \

#ifdef __cplusplus
extern "C"
{
#endif

int Ell_RegisterSolver(void (*function),  
                       const char *sname, 
                       const char *eqname);

#ifdef __cplusplus
}
#endif

#endif

#endif /* _ELLBASE_H_ */
