 /*@@
   @file      Output.c
   @date      Thu May 11 2000
   @author    Thomas Radke
   @desc
              Functions to deal with Jpeg output of variables.
   @enddesc
   @version   $Id$
 @@*/

#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "util_String.h"
#include "CactusBase/IOUtil/src/ioutil_Utils.h"
#include "ioJpegGH.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusIO_IOJpeg_Output_c)

/* define this if you want debug output */
/* #define IOJPEG_DEBUG 1 */


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
static int CheckOutputVar (int vindex);
static void SetOutputFlag (int vindex, const char *optstring, void *arg);


/*@@
   @routine    IOJpeg_OutputGH
   @date       Thu 18 April 2002
   @author     Thomas Radke
   @desc
               Loops over all variables and outputs them if necessary
   @enddesc
   @calls      IOJpeg_TimeFor
               IOJpeg_Write

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
int IOJpeg_OutputGH (const cGH *GH)
{
  int vindex, retval;
  const ioJpegGH *myGH;


  retval = 0;
  myGH = CCTK_GHExtension (GH, "IOJpeg");

  /* loop over all variables */
  for (vindex = CCTK_NumVars () - 1; vindex >= 0; vindex--)
  {
    if (IOJpeg_TimeFor (GH, vindex) &&
        IOJpeg_Write (GH, vindex, CCTK_VarName (vindex)) == 0)
    {
      /* register variable as having output this iteration */
      myGH->out_last[vindex] = GH->cctk_iteration;
      retval++;
    }
  }

  return (retval);
}


/*@@
   @routine    IOJpeg_OutputVarAs
   @date       Thu 18 April 2002
   @author     Thomas Radke
   @desc
               Unconditional output of a variable using the IOJpeg I/O method.
   @enddesc
   @calls      IOJpeg_Write

   @var        GH
   @vdesc      pointer to CCTK GH
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        fullname
   @vdesc      complete name of variable to output
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
               return code of @seeroutine IOJpeg_Write, or<BR>
               -1 if variable cannot be output by IOJpeg
   @endreturndesc
@@*/
int IOJpeg_OutputVarAs (const cGH *GH, const char *fullname, const char *alias)
{
  int vindex, retval;


  retval = -1;
  vindex = CCTK_VarIndex (fullname);

  if (CheckOutputVar (vindex) == 0)
  {
    retval = IOJpeg_Write (GH, vindex, alias);
  }

  return (retval);
}


 /*@@
   @routine    IOJpeg_TimeFor
   @date       Thu 18 April 2002
   @author     Thomas Radke
   @desc
               Decides if it is time to output a variable
               using the IOJpeg I/O method.
   @enddesc

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
int IOJpeg_TimeFor (const cGH *GH, int vindex)
{
  int result;
  char *fullname;
  ioJpegGH *myGH;


  myGH = CCTK_GHExtension (GH, "IOJpeg");
  IOJpeg_CheckSteerableParameters (myGH);

  /* check if this variable should be output */
  result = myGH->out_every[vindex] > 0 &&
           GH->cctk_iteration % myGH->out_every[vindex] == 0;
  if (result)
  {
    /* check if variable wasn't already output this iteration */
    if (myGH->out_last[vindex] == GH->cctk_iteration)
    {
      fullname = CCTK_FullName (vindex);
      CCTK_VWarn (5, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Already done IOJpeg output for '%s' in current "
                  "iteration (probably via triggers)", fullname);
      free (fullname);
      result = 0;
    }
  }

  return (result);
}


/*@@
   @routine    IOJpeg_TriggerOutput
   @date       Thu 18 April 2002
   @author     Thomas Radke
   @desc
               Triggers the output of a variable using the IOJpeg I/O method.
   @enddesc
   @calls      IOJpeg_Write

   @var        GH
   @vdesc      pointer to CCTK GH
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        vindex
   @vdesc      index of variable to output
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine IOJpeg_Write
   @endreturndesc
@@*/
int IOJpeg_TriggerOutput (const cGH *GH, int vindex)
{
  int retval;
  const ioJpegGH *myGH;


  /* do the  output */
  retval = IOJpeg_Write (GH, vindex, CCTK_VarName (vindex));
  if (retval == 0)
  {
    /* register variable as having output this iteration */
    myGH = CCTK_GHExtension (GH, "IOJpeg");
    myGH->out_last[vindex] = GH->cctk_iteration;
  }

  return (retval);
}


/*@@
   @routine    IOJpeg_CheckSteerableParameters
   @date       Mon Oct 10 2000
   @author     Thomas Radke
   @desc
               Checks if IOJpeg steerable parameters were changed
               and does appropriate re-evaluation.
   @enddesc

   @calls      CCTK_TraverseString

   @var        myGH
   @vdesc      pointer to IOJpeg grid hierarchy
   @vtype      ioJpegG *
   @vio        inout
   @endvar
@@*/
void IOJpeg_CheckSteerableParameters (ioJpegGH *myGH)
{
  int i, num_vars;
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
      CCTK_VInfo (CCTK_THORNSTRING, "Periodic IOJpeg output every %d "
                  "iterations", myGH->out_every_default);
    }
    else
    {
      CCTK_INFO ("Periodic IOJpeg output turned off");
    }
  }

  /* re-parse the 'IOJpeg::out_vars' parameter if it was changed */
  if (strcmp (out_vars, myGH->out_vars) || myGH->out_every_default != i)
  {
    num_vars = CCTK_NumVars ();
    memset (myGH->out_every, 0, num_vars * sizeof (CCTK_INT));
    if (CCTK_TraverseString (out_vars, SetOutputFlag, myGH,
                             CCTK_GROUP_OR_VAR) < 0)
    {
      CCTK_WARN (myGH->stop_on_parse_errors ? 0 : 1,
                 "error while parsing parameter 'IOJpeg::out_vars'");
    }

    if (myGH->out_every_default == i || ! CCTK_Equals (verbose, "none"))
    {
      msg = NULL;
      for (i = 0; i < num_vars; i++)
      {
        if (myGH->out_every[i] > 0)
        {
          fullname = CCTK_FullName (i);
          if (! msg)
          {
            Util_asprintf (&msg, "Periodic IOJpeg output requested for '%s'",
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

    /* save the last setting of 'IOJpeg::out_vars' parameter */
    free (myGH->out_vars);
    myGH->out_vars = strdup (out_vars);
  }
}


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
/* check if this variable can be output (static conditions) */
static int CheckOutputVar (int vindex)
{
  int groupindex;
  cGroup groupinfo;
  char *fullname;
  const char *errormsg;


  /* get the variable group information */
  groupindex = CCTK_GroupIndexFromVarI (vindex);
  CCTK_GroupData (groupindex, &groupinfo);

  errormsg = NULL;
  if (groupinfo.dim < 2 || groupinfo.dim > 3)
  {
    errormsg = "dim != [2,3]";
  }
  else if (groupinfo.grouptype != CCTK_GF && groupinfo.grouptype != CCTK_ARRAY)
  {
    errormsg = "not a grid function or array";
  }
  else if (! strncmp (CCTK_VarTypeName (groupinfo.vartype), "CCTK_COMPLEX", 12))
  {
    errormsg = "CCTK_COMPLEX variables cannot be dealt with";
  }

  if (errormsg)
  {
    fullname = CCTK_FullName (vindex);
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "No IOJpeg output for '%s': %s", fullname, errormsg);
    free (fullname);
  }

  return (errormsg != NULL);
}


/* callback for CCTK_TraverseString() to set the output flag
   for the given variable */
static void SetOutputFlag (int vindex, const char *optstring, void *arg)
{
  const ioJpegGH *myGH = (const ioJpegGH *) arg;


  if (CheckOutputVar (vindex) == 0)
  {
    myGH->out_every[vindex] = myGH->out_every_default;

    if (optstring)
    {
      IOUtil_ParseOutputFrequency (CCTK_THORNSTRING, "IOJpeg::out_vars",
                                   myGH->stop_on_parse_errors,
                                   vindex, optstring,
                                   &myGH->out_every[vindex], NULL);
    }
  }
}
