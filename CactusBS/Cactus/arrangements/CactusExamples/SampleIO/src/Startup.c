 /*@@
   @file      Startup.c
   @date      Tue 14 May 2002
   @author    Thomas Radke
   @desc
              Startup routines for thorn SampleIO
   @enddesc
   @version   $Id$
 @@*/

#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "ioSampleIOMethodGH.h"

/* the rcs ID and its dummy function to use it (this is so that you can grep
   for a version number of this source file in a Cactus executable) */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusExamples_SampleIO_Startup_c)


/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/
int SampleIO_Startup (void);


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
static void *SetupGH (tFleshConfig *config, int conv_level, cGH *GH);


 /*@@
   @routine   SampleIO_Startup
   @date      Tue 14 May 2002
   @author    Thomas Radke
   @desc
              The startup registration routine for thorn SampleIO.
              Registers the GH extension for this thorn along with its
              setup routine.
   @enddesc
   @calls     CCTK_RegisterGHExtension
              CCTK_RegisterGHExtensionSetupGH
@@*/
int SampleIO_Startup (void)
{
  int extension;


  /* register a new GH extension with a unique name */
  extension = CCTK_RegisterGHExtension ("SampleIO");
  if (extension >= 0)
  {
    CCTK_RegisterGHExtensionSetupGH (extension, SetupGH);
  }
  else
  {
    CCTK_WARN (1, "Couldn't register GH extension 'SampleIO'. "
                  "This I/O method will not be registered");
  }
  return 0;
}


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
 /*@@
   @routine   SetupGH
   @date      Tue 14 May 2002
   @author    Thomas Radke
   @desc
              Allocates and sets up SampleIO's GH extension structure
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
   @vio       unused
   @endvar

   @returntype void *
   @returndesc
               pointer to the allocated GH extension structure
   @endreturndesc
@@*/
static void *SetupGH (tFleshConfig *config, int conv_level, cGH *GH)
{
  int i, numvars;
  ioSampleIOMethodGH *myGH;
  DECLARE_CCTK_PARAMETERS


  /* suppress compiler warnings about unused variables */
  (void) (config + 0);
  (void) (conv_level + 0);
  (void) (GH + 0);

  /* allocate the GH extension and its components */
  myGH = (ioSampleIOMethodGH *) malloc (sizeof (ioSampleIOMethodGH));
  if (! myGH)
  {
    CCTK_WARN (0, "Unable to allocate memory for GH");
  }

  /* register routines as a new I/O method "SampleIO" */
  i = CCTK_RegisterIOMethod ("SampleIO");
  CCTK_RegisterIOMethodOutputGH (i, SampleIO_OutputGH);
  CCTK_RegisterIOMethodOutputVarAs (i, SampleIO_OutputVarAs);
  CCTK_RegisterIOMethodTimeToOutput (i, SampleIO_TimeFor);
  CCTK_RegisterIOMethodTriggerOutput (i, SampleIO_TriggerOutput);

  /* announcing this new I/O method to the user */
  if (! CCTK_Equals (verbose, "none"))
  {
    CCTK_INFO ("I/O Method 'SampleIO' registered");
    CCTK_INFO ("SampleIO: Print the value of 3D grid functions/arrays "
               "at a chosen location to screen");
  }

  /* now allocate and initialize the GH extension structure's components */
  numvars = CCTK_NumVars ();
  myGH->out_every = (CCTK_INT *) malloc (numvars * sizeof (CCTK_INT));
  myGH->out_last = (int *) malloc (numvars * sizeof (int));

  for (i = 0; i < numvars; i++)
  {
    myGH->out_last[i] = -1;
  }

  myGH->out_vars = strdup ("");
  myGH->out_every_default = out_every - 1;

  myGH->stop_on_parse_errors = strict_io_parameter_check;
  if (! CCTK_Equals (verbose, "none"))
  {
    CCTK_INFO ("I/O Method 'SampleIO' registered: output of individual grid "
               "points from grid variables to screen");
  }
  SampleIO_CheckSteerableParameters (GH, myGH);
  myGH->stop_on_parse_errors = 0;

  return (myGH);
}
