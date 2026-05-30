 /*@@
   @file      Startup.c
   @date      Thu 18 April 2002
   @author    Thomas Radke
   @desc
              Startup routines for IOJpeg.
   @enddesc
   @version   $Id$
 @@*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_IOMethods.h"
#include "cctk_Parameters.h"
#include "CactusBase/IOUtil/src/ioutil_Utils.h"
#include "ioJpegGH.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusIO_IOJpeg_Startup_c)


/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/

/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
static void *SetupGH (tFleshConfig *config, int conv_level, cGH *GH);


 /*@@
   @routine   IOJpeg_Startup
   @date      Thu 18 April 2002
   @author    Thomas Radke
   @desc
              The startup registration routine for IOJpeg.
              Registers the GH extensions needed for IOJpeg,
              along with its setup routine, if there are grid variables
              of at least dimensionality of 2.
   @enddesc
   @calls     CCTK_RegisterGHExtension
              CCTK_RegisterGHExtensionSetupGH
@@*/
int IOJpeg_Startup (void)
{
  if (CCTK_MaxDim () >= 2)
  {
    CCTK_RegisterGHExtensionSetupGH (CCTK_RegisterGHExtension ("IOJpeg"),
                                     SetupGH);
  }

  return 0;
}


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
 /*@@
   @routine   SetupGH
   @date      Thu 18 April 2002
   @author    Thomas Radke
   @desc
              Allocates and sets up IOJpeg's GH extension structure
   @enddesc

   @calls     CCTK_RegisterIOMethod
              CCTK_RegisterIOMethodOutputGH
              CCTK_RegisterIOMethodOutputVarAs
              CCTK_RegisterIOMethodTimeToOutput
              CCTK_RegisterIOMethodTriggerOutput

   @var       config
   @vdesc     the CCTK configuration as provided by the flesh
   @vtype     tFleshConfig *
   @vio       unused
   @endvar
   @var       conv_level
   @vdesc     the convergence level
   @vtype     int
   @vio       unused
   @endvar
   @var       GH
   @vdesc     pointer to CCTK grid hierarchy
   @vtype     cGH *
   @vio       in
   @endvar

   @returntype void *
   @returndesc
               pointer to the allocated GH extension structure
   @endreturndesc
@@*/
static void *SetupGH (tFleshConfig *config, int conv_level, cGH *GH)
{
  int i, numvars;
  ioJpegGH *myGH;
  const char *my_out_dir;
  DECLARE_CCTK_PARAMETERS


  /* suppress compiler warnings about unused variables */
  (void) (config + 0);
  (void) (conv_level + 0);

  /* allocate the GH extension and its components */
  myGH = (ioJpegGH *) malloc (sizeof (ioJpegGH));
  if (! myGH)
  {
    CCTK_WARN (0, "Unable to allocate memory for GH");
  }

  /* register the IOJpeg routines as I/O method "IOJpeg"  */
  i = CCTK_RegisterIOMethod ("IOJpeg");
  CCTK_RegisterIOMethodOutputGH (i, IOJpeg_OutputGH);
  CCTK_RegisterIOMethodOutputVarAs (i, IOJpeg_OutputVarAs);
  CCTK_RegisterIOMethodTimeToOutput (i, IOJpeg_TimeFor);
  CCTK_RegisterIOMethodTriggerOutput (i, IOJpeg_TriggerOutput);

  numvars = CCTK_NumVars ();
  myGH->out_every = malloc (numvars * sizeof (CCTK_INT));
  myGH->out_last = malloc (numvars * sizeof (int));

  for (i = 0; i < numvars; i++)
  {
    myGH->out_last[i] = -1;
  }

  myGH->out_vars = strdup ("");
  myGH->out_every_default = out_every - 1;

  /* get the name for IOJpeg's output directory */
  my_out_dir = *out_dir ? out_dir : io_out_dir;

  /* omit the directory if it's the current working dir */
  if (strcmp (my_out_dir, ".") == 0)
  {
    myGH->out_dir = strdup ("");
  }
  else
  {
    myGH->out_dir = (char *) malloc (strlen (my_out_dir) + 2);
    sprintf (myGH->out_dir, "%s/", my_out_dir);
  }

  /* create the output dir */
  i = IOUtil_CreateDirectory (GH, myGH->out_dir, 0, 0);
  if (i < 0)
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Problem creating IOJpeg output directory '%s'",
                myGH->out_dir);
  }
  else if (i >= 0 && CCTK_Equals (verbose, "full"))
  {
    CCTK_VInfo (CCTK_THORNSTRING, "IOJpeg: Output to directory '%s'",
                myGH->out_dir);
  }

  myGH->stop_on_parse_errors = strict_io_parameter_check;
  if (! CCTK_Equals (verbose, "none"))
  {
    CCTK_INFO ("I/O Method 'IOJpeg' registered: output of 2D jpeg images of "
               "grid functions/arrays");
  }
  IOJpeg_CheckSteerableParameters (myGH);
  myGH->stop_on_parse_errors = 0;

  return (myGH);
}
