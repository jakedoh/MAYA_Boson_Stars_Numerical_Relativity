 /*@@
   @file      Output.c
   @date      Tue Jan 9 1999
   @author    Gabrielle Allen
   @desc
              Functions to deal with output of variables in HDF5 file format
   @enddesc
   @version   $Id$
 @@*/


#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "util_String.h"
#include "ioHDF5GH.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusPUGHIO_IOHDF5_Output_c)


/********************************************************************
 ************************    Typedefs   *****************************
 ********************************************************************/
/* typedef for a structure to be passed to the function callback called by
   IOHDF5_OutputVarAs() via CCTK_TraverseString() */
typedef struct
{
  const cGH *GH;
  const char *fullname;
  const char *alias;
  int retval;
} callback_info_t;


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
static void OutputVarAs (int vindex, const char *optstring, void *arg);


/*@@
   @routine    IOHDF5_OutputGH
   @date       Sat March 6 1999
   @author     Gabrielle Allen
   @desc
               Loops over all variables and outputs them if necessary
   @enddesc

   @calls      IOHDF5_TimeFor
               IOHDF5_Write

   @var        GH
   @vdesc      pointer to CCTK GH
   @vtype      const cGH *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               the number of variables which were output at this iteration
               (or 0 if it wasn't time to output yet)
   @endreturndesc
@@*/
int IOHDF5_OutputGH (const cGH *GH)
{
  int vindex, retval;
  const ioHDF5GH *myGH;


  retval = 0;
  myGH = CCTK_GHExtension (GH, "IOHDF5");

  /* loop over all variables */
  for (vindex = CCTK_NumVars () - 1; vindex >= 0; vindex--)
  {
    if (IOHDF5_TimeFor (GH, vindex) &&
        IOHDF5_Write (GH, vindex, CCTK_VarName (vindex)) == 0)
    {
      /* register variable as having output this iteration */
      myGH->out_last[vindex] = GH->cctk_iteration;
      retval++;
    }
  }

  return (retval);
}


/*@@
   @routine    IOHDF5_OutputVarAs
   @date       Sat March 6 1999
   @author     Gabrielle Allen
   @desc
               Unconditional output of a variable using the IOHDF5 I/O method.
               Since the 'fullname' argument may have an option string appended
               to the actual variable's fullname, we just use
               CCTK_TraverseString() to split both, and let the callback do
               the actual work.
   @enddesc

   @calls      CCTK_TraverseString

   @var        GH
   @vdesc      pointer to CCTK GH
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        fullname
   @vdesc      complete name of variable to output, plus an optional options
               string (in curly braces)
   @vtype      const char *
   @vio        in
   @endvar
   @var        alias
   @vdesc      alias name of variable to output
               (used to generate output filename)
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine IOHDF5_Write
   @endreturndesc
@@*/
int IOHDF5_OutputVarAs (const cGH *GH, const char *fullname, const char *alias)
{
  callback_info_t info;


  info.GH = GH;
  info.fullname = fullname;
  info.alias = alias;
  CCTK_TraverseString (fullname, OutputVarAs, &info, CCTK_VAR);
  return (info.retval);
}


/*@@
   @routine    IOHDF5_TimeFor
   @date       Sat March 6 1999
   @author     Gabrielle Allen
   @desc
               Decides if it is time to output a variable
               using the IOHDF5 I/O method.
   @enddesc

   @calls      IOHDF5_CheckSteerableParameters

   @var        GH
   @vdesc      pointer to CCTK GH
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        vindex
   @vdesc      index of variable
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               1 if output should take place at this iteration, or<BR>
               0 if not
   @endreturndesc
@@*/
int IOHDF5_TimeFor (const cGH *GH, int vindex)
{
  int retval;
  char *fullname;
  ioHDF5GH *myGH;


  myGH = CCTK_GHExtension (GH, "IOHDF5");
  IOHDF5_CheckSteerableParameters (GH, myGH);

  /* check if this variable should be output */
  retval = myGH->requests[vindex] && myGH->requests[vindex]->out_every > 0 &&
           GH->cctk_iteration % myGH->requests[vindex]->out_every == 0;
  if (retval)
  {
    /* check if variable was not already output this iteration */
    if (myGH->out_last[vindex] == GH->cctk_iteration)
    {
      fullname = CCTK_FullName (vindex);
      CCTK_VWarn (6, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Already done IOHDF5 output for variable '%s' in current "
                  "iteration (probably via triggers)", fullname);
      free (fullname);
      retval = 0;
    }
  }

  return (retval);
}


/*@@
   @routine    IOHDF5_TriggerOutput
   @date       Sat March 6 1999
   @author     Gabrielle Allen
   @desc
               Triggers the output of a variable using the IOHDF5 I/O method.
   @enddesc

   @calls      IOHDF5_Write

   @var        GH
   @vdesc      pointer to CCTK GH
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        vindex
   @vdesc      index of variable
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine IOHDF5_Write
   @endreturndesc
@@*/
int IOHDF5_TriggerOutput (const cGH *GH, int vindex)
{
  int retval;
  ioHDF5GH *myGH;
  char *fullname;
  const char *varname;
  DECLARE_CCTK_PARAMETERS


  varname = CCTK_VarName (vindex);

  myGH = CCTK_GHExtension (GH, "IOHDF5");

  if (CCTK_Equals (verbose, "full"))
  {
    fullname = CCTK_FullName (vindex);
    CCTK_VInfo (CCTK_THORNSTRING, "IOHDF5_TriggerOutput: "
                                  "output of (varname, fullname) "
                                  "= ('%s', '%s')", varname, fullname);
    free (fullname);
  }

  /* do the output */
  retval = IOHDF5_Write (GH, vindex, varname);

  if (retval == 0)
  {
    /* register variable as having output this iteration */
    myGH->out_last[vindex] = GH->cctk_iteration;
  }

  return (retval);
}


/*@@
   @routine    IOHDF5_CheckSteerableParameters
   @date       Mon Oct 10 2000
   @author     Thomas Radke
   @desc
               Checks if IOHDF5 steerable parameters were changed
               and does appropriate re-evaluation.
   @enddesc

   @calls      IOHDF5Util_ParseVarsForOutput

   @var        GH
   @vdesc      pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        myGH
   @vdesc      pointer to IOHDF5 grid hierarchy
   @vtype      ioHDF5GH *
   @vio        inout
   @endvar
@@*/
void IOHDF5_CheckSteerableParameters (const cGH *GH, ioHDF5GH *myGH)
{
  int i;
  char *fullname, *msg;
  DECLARE_CCTK_PARAMETERS


  /* how often to output */
  i = myGH->out_every_default;
  myGH->out_every_default = out_every >= 0 ? out_every : io_out_every;

  /* report if frequency changed */
  if (myGH->out_every_default != i && ! CCTK_Equals (verbose, "none"))
  {
    if (myGH->out_every_default > 0)
    {
      CCTK_VInfo (CCTK_THORNSTRING, "Periodic HDF5 output every %d "
                  "iterations", myGH->out_every_default);
    }
    else
    {
      CCTK_INFO ("Periodic HDF5 output turned off");
    }
  }

  /* re-parse the 'IOHDF5::out_vars' parameter if it was changed */
  if (strcmp (out_vars, myGH->out_vars) || myGH->out_every_default != i)
  {
    IOUtil_ParseVarsForOutput (GH, CCTK_THORNSTRING, "IOHDF5::out_vars",
                               myGH->stop_on_parse_errors, out_vars,
                               myGH->out_every_default, -1.0, myGH->requests);

    if (myGH->out_every_default == i || ! CCTK_Equals (verbose, "none"))
    {
      msg = NULL;
      for (i = CCTK_NumVars () - 1; i >= 0; i--)
      {
        if (myGH->requests[i])
        {
          fullname = CCTK_FullName (i);
          if (! msg)
          {
            Util_asprintf (&msg, "Periodic HDF5 output requested for '%s'",
                           fullname);
          }
          else
          {
            char *tmp = msg;
            Util_asprintf (&msg, "%s, '%s'", msg, fullname);
            free (tmp);
          }
          free (fullname);
        }
      }
      if (msg)
      {
        CCTK_INFO (msg);
        free (msg);
      }
    }

    /* save the last setting of 'IOHDF5::out_vars' parameter */
    free (myGH->out_vars);
    myGH->out_vars = strdup (out_vars);
  }
}


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
/*@@
   @routine    OutputVarAs
   @date       Fri 6 June 2003
   @author     Thomas Radke
   @desc
               Unconditional output of a variable using the IOHDF5 I/O method.
   @enddesc

   @calls      IOUtil_ParseVarsForOutput
               IOHDF5_Write
               IOUtil_FreeIORequest

   @var        index
   @vdesc      variable index of a variable found in the 'fullname' argument
               string given to IOHDF5_OutputVarAs()
   @vtype      int
   @vio        in
   @endvar
   @var        optstring
   @vdesc      option string attached to the variable in the fullname
   @vtype      const char *
   @vio        in
   @endvar
   @var        arg
   @vdesc      user-supplied argument to the callback function
   @vtype      void *
   @vio        in
   @endvar
@@*/
static void OutputVarAs (int vindex, const char *optstring, void *arg)
{
  int oneshot;
  const ioHDF5GH *myGH;
  callback_info_t *info = arg;
  DECLARE_CCTK_PARAMETERS


  /* don't warn about unused function arguments */
  (void) (optstring + 0);

  if (CCTK_Equals (verbose, "full"))
  {
    CCTK_VInfo (CCTK_THORNSTRING, "IOHDF5_OutputVarAs: fullname, alias, "
                "index = (%s, %s, %d)", info->fullname, info->alias, vindex);
  }

  /* check whether the variable already has an I/O request entry */
  myGH = CCTK_GHExtension (info->GH, "IOHDF5");
  oneshot = myGH->requests[vindex] == NULL;
  if (oneshot)
  {
    IOUtil_ParseVarsForOutput (info->GH, CCTK_THORNSTRING, "IOHDF5::out_vars",
                               0, info->fullname, 1, -1.0, myGH->requests);
  }

  /* do the output */
  info->retval = IOHDF5_Write (info->GH, vindex, info->alias);

  if (oneshot && myGH->requests[vindex])
  {
    IOUtil_FreeIORequest (&myGH->requests[vindex]);
  }
}
