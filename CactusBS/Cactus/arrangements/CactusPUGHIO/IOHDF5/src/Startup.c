 /*@@
   @file      Startup.c
   @date      Fri May 21 1999
   @author    Thomas Radke
   @desc
              Startup and termination routines for IOHDF5.
   @enddesc
   @version   $Id$
 @@*/


#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "CactusBase/IOUtil/src/ioGH.h"
#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"
#include "ioHDF5GH.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusPUGHIO_IOHDF5_Startup_c)


/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/
int IOHDF5_Startup (void);


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
static void *SetupGH (tFleshConfig *config, int conv_level, cGH *GH);


 /*@@
   @routine   IOHDF5_Startup
   @date      Fri May 21 1999
   @author    Thomas Radke
   @desc
              The startup registration routine for IOHDF5.
              Registers the GH extensions needed for IOHDF5
              along with its setup routine.
   @enddesc

   @calls     CCTK_RegisterGHExtension
              CCTK_RegisterGHExtensionSetupGH
@@*/
int IOHDF5_Startup (void)
{
  CCTK_RegisterGHExtensionSetupGH (CCTK_RegisterGHExtension ("IOHDF5"),SetupGH);
  return 0;
}


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
 /*@@
   @routine    SetupGH
   @date       Mon Jun 19 2000
   @author     Thomas Radke
   @desc
               Allocates and sets up IOHDF5's GH extension structure.
   @enddesc

   @calls      CCTK_RegisterIOMethod
               CCTK_RegisterIOMethodOutputGH
               CCTK_RegisterIOMethodOutputVarAs
               CCTK_RegisterIOMethodTimeToOutput
               CCTK_RegisterIOMethodTriggerOutput
               CCTK_TimerCreate
               CCTK_TimerDestroyI
               CCTK_TimerResetI
               IOUtil_CreateDirectory

   @var        config
   @vdesc      the CCTK configuration as provided by the flesh
   @vtype      tFleshConfig *
   @vio        usused
   @endvar
   @var        conv_level
   @vdesc      the convergence level
   @vtype      int
   @vio        unused
   @endvar
   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      cGH *
   @vio        in
   @endvar

   @returntype void *
   @returndesc
               pointer to the new GH extension structure
   @endreturndesc
@@*/
static void *SetupGH (tFleshConfig *config, int conv_level, cGH *GH)
{
  int i, numvars;
  ioHDF5GH *myGH;
  const ioGH *ioUtilGH;
  const char *my_out_dir;
  const char *timer_names[4] = {"IOHDF5 time to dump parameters",
                                "IOHDF5 time to dump variables",
                                "IOHDF5 total time to checkpoint",
                                "IOHDF5 time to recover"};
  DECLARE_CCTK_PARAMETERS


  /* suppress compiler warnings about unused variables */
  (void) (config + 0);
  (void) (conv_level + 0);

  /* register the IOHDF5 routines as a new I/O method */
  i = CCTK_RegisterIOMethod ("IOHDF5");
  CCTK_RegisterIOMethodOutputGH (i, IOHDF5_OutputGH);
  CCTK_RegisterIOMethodOutputVarAs (i, IOHDF5_OutputVarAs);
  CCTK_RegisterIOMethodTimeToOutput (i, IOHDF5_TimeFor);
  CCTK_RegisterIOMethodTriggerOutput (i, IOHDF5_TriggerOutput);

  /* register the IOHDF5 recovery routine with thorn IOUtil */
  if (IOUtil_RegisterRecover ("IOHDF5 recovery", IOHDF5_Recover) < 0)
  {
    CCTK_WARN (1, "Failed to register IOHDF5 recovery routine");
  }

  /* allocate a new GH extension structure */
  numvars = CCTK_NumVars ();
  myGH            = malloc (sizeof (ioHDF5GH));
  myGH->out_last  = malloc (numvars * sizeof (int));
  myGH->requests  = calloc (numvars, sizeof (ioRequest *));
  myGH->cp_filename_index = 0;
  myGH->checkpoint_keep = abs (checkpoint_keep);
  myGH->cp_filename_list = calloc (myGH->checkpoint_keep, sizeof (char *));
  myGH->out_vars = strdup ("");
  myGH->out_every_default = out_every - 1;
  myGH->last_checkpoint_iteration = -1;

  myGH->stop_on_parse_errors = strict_io_parameter_check;
  if (! CCTK_Equals (verbose, "none"))
  {
    CCTK_INFO ("I/O Method 'IOHDF5' registered: output of grid variables and "
               "hyperslabs thereof to HDF5 files");
  }
  IOHDF5_CheckSteerableParameters (GH, myGH);
  myGH->stop_on_parse_errors = 0;

  /* get the name of IOHDF5's output directory */
  my_out_dir = out_dir;
  if (*my_out_dir == 0)
  {
    my_out_dir = io_out_dir;
  }

  /* skip the directory pathname if output goes into current directory */
  if (strcmp (my_out_dir, "."))
  {
    i = strlen (my_out_dir);
    if (CCTK_Equals (out_mode, "onefile") || ! strstr (my_out_dir, "%u"))
    {
      myGH->out_dir = malloc (i + 2);
      strcpy (myGH->out_dir, my_out_dir);
      myGH->out_dir[i] = '/';
      myGH->out_dir[i+1] = 0;
    }
    else
    {
      myGH->out_dir = malloc (i + 20);
      sprintf (myGH->out_dir, my_out_dir, CCTK_MyProc (GH));
      strcat (myGH->out_dir, "/");
    }
  }
  else
  {
    myGH->out_dir = strdup ("");
  }

  /* create the output directory */
  ioUtilGH = CCTK_GHExtension (GH, "IO");
  i = IOUtil_CreateDirectory (GH, myGH->out_dir,
                              ! CCTK_Equals (out_mode, "onefile"),
                              ioUtilGH->ioproc);
  if (i < 0)
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Problem creating HDF5 output directory '%s'", myGH->out_dir);
  }
  else if (i > 0 && CCTK_Equals (verbose, "full"))
  {
    CCTK_VInfo (CCTK_THORNSTRING,
                "HDF5 output directory '%s' already exists", myGH->out_dir);
  }

  for (i = 0; i < numvars; i++)
  {
    myGH->out_last [i] = -1;
  }

  myGH->open_output_files = NULL;

  /* create timers if timing info was requested */
  myGH->print_timing_info = print_timing_info;
  if (myGH->print_timing_info)
  {
    for (i = 0; i < IOHDF5_NUM_TIMERS; i++)
    {
      if ((myGH->timers[i] = CCTK_TimerCreate (timer_names[i])) < 0)
      {
        break;
      }
    }
    if (i != IOHDF5_NUM_TIMERS)
    {
      CCTK_WARN (1, "Could not create timers for checkpoint/recovery ! "
                    "No timer information will be available.");
      while (--i >= 0)
      {
        CCTK_TimerDestroyI (myGH->timers[i]);
      }
      myGH->print_timing_info = 0;
    }
    else
    {
      CCTK_TimerResetI (myGH->timers[CP_TOTAL_TIMER]);
      CCTK_TimerResetI (myGH->timers[RECOVERY_TIMER]);
    }
  }

  return (myGH);
}
