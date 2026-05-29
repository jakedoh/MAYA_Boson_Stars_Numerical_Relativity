#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cctk.h>
#include <cctk_FortranString.h>

#include <TATelliptic.h>



/* Linked list with all registered solvers */
static struct solverinfo {
  struct solverinfo * next;	/* linked list cdr pointer */
  solvefunc solver;		/* the solver */
  const char * solvername;	/* (unique) name of the solver */
} * solvers = 0;



/* Register a solver */
int TATelliptic_RegisterSolver (const solvefunc solver,
				const char * const solvername)
{
  struct solverinfo * p;
  char * s;
  if (!solver) {
    CCTK_WARN (CCTK_WARN_ALERT,
               "Cannot register solver: Solver function pointer is null");
    return -1;
  }
  assert (solver);
  if (!solvername) {
    CCTK_WARN (CCTK_WARN_ALERT,
               "Cannot register solver: Solver name pointer is null");
    return -1;
  }
  assert (solvername);
  for (p = solvers; p; p = p->next) {
    if (strcmp(solvername, p->solvername) == 0) {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Cannot register solver: A solver with the same name \"%s\" has already been registered", solvername);
      return -2;
    }
  }
  p = malloc(sizeof(struct solverinfo));
  assert (p);
  p->solver = solver;
  s = malloc(strlen(solvername)+1);
  assert (s);
  strcpy (s, solvername);
  p->solvername = s;
  p->next = solvers;
  solvers = p;
  return 0;
}



/* Get a solver that has been registered */
solvefunc TATelliptic_GetSolver (const char * const solvername)
{
  struct solverinfo * p;
  if (!solvername) {
    CCTK_WARN (CCTK_WARN_ALERT,
               "Cannot return solver: Pointer to requested solver name is null");
    return 0;
  }
  assert (solvername);
  for (p = solvers; p; p = p->next) {
    if (strcmp(solvername, p->solvername) == 0) {
      return p->solver;
    }
  }
  return 0;
}



/* Call a solver that has been registered (the solver must exist) */
int TATelliptic_CallSolver (const cGH * const cctkGH,
			    const int * const var,
			    const int * const res,
			    const int nvars,
			    const int options_table,
			    const calcfunc calcres,
			    const calcfunc applybnds,
			    void * const userdata,
			    const char * const solvername)
{
  solvefunc solver;
  if (!solvername) {
    CCTK_WARN (CCTK_WARN_ALERT,
               "Cannot call solver: Pointer to solver name is null");
    return -1;
  }
  solver = TATelliptic_GetSolver (solvername);
  if (!solver) {
    CCTK_VWarn (CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
		"Cannot call solver: There is no solver with the requested name \"%s\".  (Maybe the solver thorn has not been activated.)", solvername);
    return -2;
  }
  assert (solver);
  return solver (cctkGH, var, res, nvars, options_table,
		 calcres, applybnds, userdata);
}



void CCTK_FCALL CCTK_FNAME(TATelliptic_CallSolver)
     (int * const ierr,
      const cGH * const * const cctkGH,
      const int * const var,
      const int * const res,
      const int nvars,
      const int options_table,
      const calcfunc calcres,
      const calcfunc applybnds,
      void * const userdata,
      ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(solvername)
  *ierr = TATelliptic_CallSolver (*cctkGH, var, res, nvars, options_table,
				  calcres, applybnds, userdata, solvername);
  free (solvername);
}
