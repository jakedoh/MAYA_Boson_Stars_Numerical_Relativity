 /*@@
   @file      Output.c
   @date      Tue 14 May 2002
   @author    Thomas Radke
   @desc
              Functions registered as I/O method "SampleIO"
   @enddesc
   @version   $Id$
 @@*/

#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "util_String.h"
#include "CactusBase/IOUtil/src/ioutil_Utils.h"
#include "ioSampleIOMethodGH.h"

/* the rcs ID and its dummy function to use it (this is so that you can grep
   for a version number of this source file in a Cactus executable) */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusExamples_SampleIO_Output_c)


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
static int CheckOutputVar (int vindex, const cGH *GH);
static void SetOutputFrequency (int vindex, const char *optstring, void *arg);


/*@@
   @routine    SampleIO_OutputGH
   @date       Tue 14 May 2002
   @author     Thomas Radke
   @desc
               Loops over all variables and outputs them if necessary
   @enddesc
   @calls      SampleIO_TimeFor
               SampleIO_Write

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
int SampleIO_OutputGH (const cGH *GH)
{
  int vindex, retval;
  const ioSampleIOMethodGH *myGH;


  retval = 0;
  myGH = CCTK_GHExtension (GH, "SampleIO");

  /* loop over all variables */
  for (vindex = CCTK_NumVars () - 1; vindex >= 0; vindex--)
  {
    if (SampleIO_TimeFor (GH, vindex) &&
        SampleIO_Write (GH, vindex) == 0)
    {
      /* register variable as having output this iteration */
      myGH->out_last[vindex] = GH->cctk_iteration;
      retval++;
    }
  }

  return (retval);
}


/*@@
   @routine    SampleIO_OutputVarAs
   @date       Tue 14 May 2002
   @author     Thomas Radke
   @desc
               Unconditional output of a variable using the "SampleIO"
               I/O method
   @enddesc
   @calls      SampleIO_Write

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
   @vtype      const char *
   @vio        unused
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine SampleIO_Write, or<BR>
               -1 if variable cannot be output by "SampleIO"
   @endreturndesc
@@*/
int SampleIO_OutputVarAs (const cGH *GH, const char *fullname,
                          const char *alias)
{
  int vindex, retval;


  /* avoid compiler warning about unused function argument */
  (void) (alias + 0);

  retval = -1;
  vindex = CCTK_VarIndex (fullname);

  /* make sure this variable can be output by "SampleIO" */
  if (CheckOutputVar (vindex, GH))
  {
    retval = SampleIO_Write (GH, vindex);
  }

  return (retval);
}


 /*@@
   @routine    SampleIO_TimeFor
   @date       Tue 14 May 2002
   @author     Thomas Radke
   @desc
               Decides if it is time to output a variable using the "SampleIO"
               I/O method.
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
int SampleIO_TimeFor (const cGH *GH, int vindex)
{
  int result;
  char *fullname;
  ioSampleIOMethodGH *myGH;


  /* check whether steerable parameters have changed */
  myGH = CCTK_GHExtension (GH, "SampleIO");
  SampleIO_CheckSteerableParameters (GH, myGH);

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
                  "Already done SampleIO output for '%s' in current "
                  "iteration (probably via triggers)", fullname);
      free (fullname);
      result = 0;
    }
  }

  return (result);
}


/*@@
   @routine    SampleIO_TriggerOutput
   @date       Tue 14 May 2002
   @author     Thomas Radke
   @desc
               Triggers the output of a variable using the "SampleIO" I/O method
   @enddesc
   @calls      SampleIO_Write

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
               return code of @seeroutine SampleIO_Write
   @endreturndesc
@@*/
int SampleIO_TriggerOutput (const cGH *GH, int vindex)
{
  int retval;
  const ioSampleIOMethodGH *myGH;


  /* do the output */
  retval = SampleIO_Write (GH, vindex);
  if (retval == 0)
  {
    /* register variable as having output this iteration */
    myGH = CCTK_GHExtension (GH, "SampleIO");
    myGH->out_last[vindex] = GH->cctk_iteration;
  }

  return (retval);
}


/*@@
   @routine    SampleIO_CheckSteerableParameters
   @date       Mon Oct 10 2000
   @author     Thomas Radke
   @desc
               Checks if SampleIO steerable parameters were changed
               and does appropriate re-evaluation.
   @enddesc

   @calls      CCTK_TraverseString

   @var        GH
   @vdesc      pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        myGH
   @vdesc      pointer to SampleIO grid hierarchy
   @vtype      ioSampleIOMethodGH *
   @vio        inout
   @endvar
@@*/
void SampleIO_CheckSteerableParameters (const cGH *GH, ioSampleIOMethodGH *myGH)
{
  int i, num_vars;
  char *fullname, *msg;
  DECLARE_CCTK_PARAMETERS


  /* check the 'out_every' parameter from both implementations
     "SampleIO" and "IO" (the first overriding the latter) */
  i = myGH->out_every_default;
  myGH->out_every_default = out_every >= 0 ? out_every : io_out_every;

  /* report if periodic output frequency changed */
  if (myGH->out_every_default != i && ! CCTK_Equals (verbose, "none"))
  {
    if (myGH->out_every_default > 0)
    {
      CCTK_VInfo (CCTK_THORNSTRING, "Periodic SampleIO output every "
                  "%d iterations", myGH->out_every_default);
    }
    else
    {
      CCTK_INFO ("Periodic SampleIO output turned off");
    }
  }

  /* re-parse the 'SampleIO::out_vars' parameter if it was changed */
  if (strcmp (out_vars, myGH->out_vars) || myGH->out_every_default != i)
  {
    num_vars = CCTK_NumVars ();

    /* reset all variables' output frequencies to zero (disable output) */
    memset (myGH->out_every, 0, num_vars * sizeof (CCTK_INT));

    /* set output frequency for selected variables by traversing the 'out_vars'
       string */
    if (CCTK_TraverseString (out_vars, SetOutputFrequency, &GH,
                             CCTK_GROUP_OR_VAR) < 0)
    {
      CCTK_WARN (myGH->stop_on_parse_errors ? 0 : 1,
                 "error while parsing parameter 'SampleIO::out_vars'");
    }

    /* report the variables to do periodic output for */
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
            Util_asprintf (&msg, "Periodic SampleIO output requested "
                           "for '%s'", fullname);
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

    /* save the last setting of the 'SampleIO::out_vars' parameter */
    free (myGH->out_vars);
    myGH->out_vars = strdup (out_vars);
  }
}


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
/*
 * check if this variable can be output by I/O method "SampleIO"
 * A variable can be output if
 *   - it is a grid function or array
 *   - it is not of any of the CCTK_COMPLEX data types
 *   - it has 3 dimensions
 *   - its extents are larger than the selected point position to output
 */
static int CheckOutputVar (int vindex, const cGH *GH)
{
  int groupindex;
  cGroup groupinfo;
  char *fullname;
  int extent[3];
  const char *errormsg;
  DECLARE_CCTK_PARAMETERS


  /* get the variable group information */
  groupindex = CCTK_GroupIndexFromVarI (vindex);
  CCTK_GroupData (groupindex, &groupinfo);

  errormsg = NULL;
  if (groupinfo.grouptype != CCTK_GF && groupinfo.grouptype != CCTK_ARRAY)
  {
    errormsg = "is not a grid function or array";
  }
  else if (! strncmp (CCTK_VarTypeName (groupinfo.vartype), "CCTK_COMPLEX", 12))
  {
    errormsg = "CCTK_COMPLEX variables cannot be dealt with";
  }
  else if (groupinfo.dim != 3)
  {
    errormsg = "grid function/array dimension is other than 3";
  }
  else
  {
    /* get the extents of the variable */
    CCTK_GroupgshGI (GH, 3, extent, groupindex);
    if (extent[0] <= point_x || extent[1] <= point_y || extent[2] <= point_z)
    {
      errormsg = "grid point location exceeds grid function/array extents";
    }
  }

  if (errormsg)
  {
    fullname = CCTK_FullName (vindex);
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "No SampleIO output for '%s': %s", fullname, errormsg);
    free (fullname);
  }

  return (errormsg == NULL);
}


/*
 * callback for CCTK_TraverseString() to set the output flag
 * for the given variable
 */
static void SetOutputFrequency (int vindex, const char *optstring, void *arg)
{
  const cGH *GH = *(const cGH **) arg;
  const ioSampleIOMethodGH *myGH;


  /* only set frequency if this variable can be output */
  if (CheckOutputVar (vindex, GH))
  {
    /* take the default output frequency, or set the value from the option
       string if supplied */
    myGH = CCTK_GHExtension (GH, "SampleIO");
    myGH->out_every[vindex] = myGH->out_every_default;

    if (optstring)
    {
      IOUtil_ParseOutputFrequency (CCTK_THORNSTRING, "SampleIO::out_vars",
                                   myGH->stop_on_parse_errors, vindex,
                                   optstring, &myGH->out_every[vindex], NULL);
    }
  }
}
