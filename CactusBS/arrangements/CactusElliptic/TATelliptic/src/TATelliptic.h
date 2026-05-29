#ifndef TATELLIPTIC_H
#define TATELLIPTIC_H

#ifdef __cplusplus
extern "C" {
#endif

#include "cctk.h"

/* A calcfunc is a user supplied function that does some calculation.
   When a calcfunc returns non-zero, the solving is aborted. */
typedef int (* calcfunc) (const cGH * cctkGH,
			  int options_table,
			  void * userdata);

/* A solvefunc is a solver that the user can call.  The user has to
   provide calcres, a function that calculates the residual, and
   applybnds, a function that aplies the boundary conditions to the
   solution.  calcres has to calculate the residual everywhere in the
   domain, including the boundaries.  A solvefunc returns the solution
   in var with the boundaries applied, and also returns the current
   residual in res.  A non-zero return value indicates that an error
   occurred. */
typedef int (* solvefunc) (const cGH * cctkGH,
			   const int * var,
			   const int * res,
			   int nvars,
			   int options_table,
			   calcfunc calcres,
			   calcfunc applybnds,
			   void * userdata);

int TATelliptic_RegisterSolver (solvefunc solver,
				const char * solvername);

solvefunc TATelliptic_GetSolver (const char * solvername);

int TATelliptic_CallSolver (const cGH * cctkGH,
			    const int * var,
			    const int * res,
			    int nvars,
			    int options_table,
			    calcfunc calcres,
			    calcfunc applybnds,
			    void * userdata,
			    const char * solvername);

#ifdef __cplusplus
}
#endif

#endif /* TATELLIPTIC_H */
