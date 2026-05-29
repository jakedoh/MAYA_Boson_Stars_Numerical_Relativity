 /*@@
   @header    Ell_Register.h
   @date      
   @author    Gerd Lanferman
   @desc 
   Registration routines for elliptic classes and solvers   
   @enddesc 
   @version $Header$
 @@*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "EllBase.h"

#include "StoreNamedData.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusElliptic_EllBase_Ell_Register_c)

static pNamedData *EqNameDB;

int Ell_RegisterSolver(void (*function),  
                       const char *sname, 
                       const char *eqname); 

int Ell_RegisterEq(void *(function)(const char *, void*), const char *eqname);


/* Ell_RegisterEq takes a routine ("function") and registers that routine 
   under the name "eqname" in the EqNameDB"
   Application: Call Ell_Register with the routine that registers a solver 
   for a elliptic equation class. */


int Ell_RegisterEq(void *(function)(const char *, void*), const char *eqname) 
{

  DECLARE_CCTK_PARAMETERS

  int retval = ELL_FAILURE;

  /* Register if function not already there with this name */

  if (!GetNamedData(EqNameDB, eqname))
  {
    if (StoreNamedData(&EqNameDB, eqname, (void*)function)==0)
    {
      if CCTK_EQUALS(elliptic_verbose,"yes") 
      {
        CCTK_VInfo(CCTK_THORNSTRING,"Registered elliptic class: %s",eqname);
      }
      retval = ELL_SUCCESS;
    }
    else
    {
      CCTK_VInfo(CCTK_THORNSTRING,"Failed to register elliptic class: %s",
                 eqname);
    }
  }
  else
  {
    CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
               "Elliptic class %s already registered",eqname);
    retval = ELL_CLASSEXISTS;
  }

  return retval;

  

}

/* Ell_RegisterSolver takes a routine ("function") and registers that 
   routine under the name  "sname" with the database specified by "eqname".
   So, how do we get to the database, after all it needs to be hardcoded
   as pNamedData ? Well, we know its name and in Ell_RegisterEq we have 
   registered a function under the equation class name. 
   We now get that function and use that function to register the solver.
   Sounds confusing, well it is. The advantage is, that somebody can come 
   up with a new equation class and can keep the database (the pNamedData
   declaration) in his own routine and does not have to put it in a central
   place. Amen*/

int Ell_RegisterSolver(void (*function),  
                       const char *sname, 
                       const char *eqname) 
{

  DECLARE_CCTK_PARAMETERS

  int retval=ELL_FAILURE;
  int ierr;
  int (*fn)(void *, const char *);

  fn = (int(*)(void (*func), const char *solvename))
    GetNamedData(EqNameDB, eqname);

  if (fn)
  {

    ierr = fn(function,sname);

    if (ierr==ELL_SUCCESS)
    {

      if CCTK_EQUALS(elliptic_verbose,"yes") 
      {
        CCTK_VInfo(CCTK_THORNSTRING,
                   "Registered elliptic solver %s for %s",sname,eqname);
      }
      
      retval = ELL_SUCCESS;
    }
    else if (ierr==ELL_SOLVEREXISTS)
    {
      CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                 "Ell_RegisterSolver: "
                 "Registered second solver %s for %s",sname,eqname);
      retval = ELL_SOLVEREXISTS;
    }
    else 
    {
      CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
                 "Ell_RegisterSolver: "
                 "Failed to register solver %s for %s",sname,eqname);
    }

  }
  else
  {
    CCTK_WARN(0,"Ell_RegisterSolver: Cannot get function in EqName");
    retval = ELL_NOCLASS;
  }
    
  return retval; 
}
