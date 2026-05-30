 /*@@
   @file      DumpGH.c
   @date      Wed Jun 10 14:13:35 1998
   @author    Paul Walker
   @desc
              Checkpoint routines scheduled at CCTK_CPINITIAL, CCTK_CHECKPOINT,
              and CCTK_TERMINATE.
              They check the IO checkpointing parameters and - if it's time
              to do so - call the routine which finally creates a checkpoint.
   @enddesc
   @version   $Id$
 @@*/


#include "cctk.h"
#include "cctk_Parameters.h"
#include "CactusBase/IOUtil/src/ioGH.h"
#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"
#include "CactusPUGHIO/IOHDF5Util/src/ioHDF5UtilGH.h"
#ifdef CCTK_MPI
#include "CactusPUGH/PUGH/src/include/pugh.h"
#endif
#include "ioHDF5GH.h"

#include <stdlib.h>
#include <string.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusPUGHIO_IOHDF5_DumpGH_c)


/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/
void IOHDF5_InitialDataCheckpoint (cGH *GH);
void IOHDF5_EvolutionCheckpoint (cGH *GH);
void IOHDF5_TerminationCheckpoint (cGH *GH);


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
static int Checkpoint (const cGH *GH, int called_from);


 /*@@
   @routine    IOHDF5_InitialDataCheckpoint
   @date       Fri Aug 21 14:46:28 1998
   @author     Gerd Lanfermann
   @desc
               This routine is registered at CCTK_CPINITIAL.
               It checks if initial data should be checkpointed.
   @enddesc

   @calls      Checkpoint

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
@@*/
void IOHDF5_InitialDataCheckpoint (cGH *GH)
{
  DECLARE_CCTK_PARAMETERS


  if (checkpoint && checkpoint_ID)
  {
    if (! CCTK_Equals (verbose, "none"))
    {
      CCTK_INFO ("---------------------------------------------------------");
      CCTK_INFO ("Dumping initial data checkpoint");
      CCTK_INFO ("---------------------------------------------------------");
    }
    Checkpoint (GH, CP_INITIAL_DATA);
  }
}


 /*@@
   @routine    IOHDF5_EvolutionCheckpoint
   @date       Fri Aug 21 14:38:25 1998
   @author     Gabrielle Allen
   @desc
               This routine is registered at CCTK_CHECKPOINT.
               It periodically checks if it's time to checkpoint evolution data.
   @enddesc

   @calls      Checkpoint

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
@@*/
void IOHDF5_EvolutionCheckpoint (cGH *GH)
{
  DECLARE_CCTK_PARAMETERS


  if (checkpoint &&
      ((checkpoint_every > 0 && GH->cctk_iteration % checkpoint_every == 0) ||
       checkpoint_next))
  {
    if (! CCTK_Equals (verbose, "none"))
    {
      CCTK_INFO ("---------------------------------------------------------");
      CCTK_VInfo (CCTK_THORNSTRING, "Dumping periodic checkpoint at "
                  "iteration %d", GH->cctk_iteration);
      CCTK_INFO ("---------------------------------------------------------");
    }
    Checkpoint (GH, CP_EVOLUTION_DATA);

    /* reset the 'checkpoint_next' parameter */
    if (checkpoint_next)
    {
      CCTK_ParameterSet ("checkpoint_next", CCTK_THORNSTRING, "no");
    }
  }
}


 /*@@
   @routine    IOHDF5_TerminationCheckpoint
   @date       Fri Aug 21 14:40:21 1998
   @author     Gabrielle Allen
   @desc
               This routine is registered at CCTK_TERMINATE.
               It checks if the last iteration should be checkpointed.
   @enddesc

   @calls      Checkpoint

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
@@*/
void IOHDF5_TerminationCheckpoint (cGH *GH)
{
  const ioHDF5GH *myGH;
  DECLARE_CCTK_PARAMETERS


  myGH = CCTK_GHExtension (GH, "IOHDF5");
  if (checkpoint && checkpoint_on_terminate && myGH)
  {
    if (myGH->last_checkpoint_iteration < GH->cctk_iteration)
    {
      if (! CCTK_Equals (verbose, "none"))
      {
        CCTK_INFO ("---------------------------------------------------------");
        CCTK_VInfo (CCTK_THORNSTRING, "Dumping termination checkpoint at "
                    "iteration %d", GH->cctk_iteration);
        CCTK_INFO ("---------------------------------------------------------");
      }
      Checkpoint (GH, CP_EVOLUTION_DATA);
    }
    else if (! CCTK_Equals (verbose, "none"))
    {
      CCTK_INFO ("---------------------------------------------------------");
      CCTK_VInfo (CCTK_THORNSTRING, "Termination checkpoint already dumped "
                  "as last evolution checkpoint at iteration %d",
                  myGH->last_checkpoint_iteration);
      CCTK_INFO ("---------------------------------------------------------");
    }
  }
}


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
 /*@@
   @routine    Checkpoint
   @date       Fri Aug 21 15:13:13 1998
   @author     Paul Walker
   @desc
               The heart of checkpointing.
               Called by the different wrappers, this routine creates
               a new checkpoint file and then dumps away using the
               dump routines from IOHDF5Util.
   @enddesc

   @calls      IOHDF5Util_DumpGH

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        called_from
   @vdesc      flag indicating where this routine was called from
               (either CHECKPOINT_ID, CHECKPOINT, or filereader)
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               -1 if checkpoint file could not be created
               or returncode of @seeroutine IOUtilHDF5_DumpGH
   @endreturndesc
@@*/
static int Checkpoint (const cGH *GH, int called_from)
{
  hid_t file;
  int myproc;
  CCTK_INT retval;
  ioGH *ioUtilGH;
  ioHDF5GH *myGH;
#ifdef CCTK_MPI
  CCTK_INT tmp;
#endif
  char *filename, *tempname;
  const char *timer_descriptions[] = {"Time to dump parameters: ",
                                      "Time to dump datasets:   ",
                                      "Total time to checkpoint:"};
  DECLARE_CCTK_PARAMETERS


  retval = 0;
  file = -1;
  myproc = CCTK_MyProc (GH);

  ioUtilGH = CCTK_GHExtension (GH, "IO");
  myGH = CCTK_GHExtension (GH, "IOHDF5");

  /* check if IOHDF5 was registered as I/O method */
  if (myGH == NULL)
  {
    CCTK_WARN (2, "No IOHDF5 I/O methods registered");
    return (-1);
  }

  /* start the CP_TOTAL_TIMER timer */
  if (myGH->print_timing_info)
  {
    CCTK_TimerStartI (myGH->timers[CP_TOTAL_TIMER]);
  }

  /* get the filenames for both the temporary and real checkpoint file */
  filename = IOUtil_AssembleFilename (GH, NULL, "", ".h5", called_from,
                                      myproc / ioUtilGH->ioproc_every,
                                      ioUtilGH->unchunked);
  tempname = IOUtil_AssembleFilename (GH, NULL, ".tmp", ".h5", called_from,
                                      myproc / ioUtilGH->ioproc_every,
                                      ioUtilGH->unchunked);

  /* Now open the file */
  if (myproc == ioUtilGH->ioproc)
  {
    if (CCTK_Equals (verbose, "full"))
    {
      CCTK_VInfo (CCTK_THORNSTRING, "Creating checkpoint file '%s'", tempname);
    }

    file = H5Fcreate (tempname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Can't open checkpoint file '%s'. Checkpointing is skipped.",
                  tempname);
      retval = -1;
    }
  }

  /* dump the GH */
  if (retval == 0)
  {
    retval = IOHDF5Util_DumpGH (GH,
                                myGH->print_timing_info ? myGH->timers : NULL,
                                file);
  }

#ifdef CCTK_MPI
  /* find out whether all IO processors succeeded in writing the checkpoint */
  if (ioUtilGH->nioprocs > 1)
  {
    tmp = retval;
    CACTUS_MPI_ERROR (MPI_Allreduce (&tmp, &retval, 1, PUGH_MPI_INT, MPI_SUM,
                                     PUGH_pGH (GH)->PUGH_COMM_WORLD));
  }
#endif

  /* close the temporary checkpoint file and rename it to the real file */
  if (myproc == ioUtilGH->ioproc)
  {
    if (CCTK_Equals (verbose, "full"))
    {
      CCTK_VInfo (CCTK_THORNSTRING, "Closing temporary checkpoint file '%s' "
                                    "and renaming it to '%s'",
                                    tempname, filename);
    }

    if (file >= 0)
    {
      HDF5_ERROR (H5Fclose (file));
    }

    /* delete the oldest checkpoint file if requested
       and put the new filename into the ring buffer */
    if (retval == 0)
    {
#ifdef _WIN32
      /* Windows' rename(2) routine isn't POSIX compatible
         in that it would unlink an existing file :-( */
      unlink (filename);
#endif
      if (rename (tempname, filename))
      {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Could not rename temporary checkpoint file '%s' to '%s'",
                    tempname, filename);
        retval = -1;
      }
      else if (called_from == CP_EVOLUTION_DATA)
      {
        if (checkpoint_keep > 0)
        {
          if (myGH->cp_filename_list[myGH->cp_filename_index])
          {
            remove (myGH->cp_filename_list[myGH->cp_filename_index]);
            free (myGH->cp_filename_list[myGH->cp_filename_index]);
            myGH->cp_filename_list[myGH->cp_filename_index] = NULL;
          }
          myGH->cp_filename_list[myGH->cp_filename_index] = strdup (filename);
          myGH->cp_filename_index = (myGH->cp_filename_index+1) % checkpoint_keep;
          if (myGH->checkpoint_keep != checkpoint_keep)
          {
            char **cp_filename_list = calloc (checkpoint_keep, sizeof (char *));
            int min = myGH->checkpoint_keep < checkpoint_keep ?
                      myGH->checkpoint_keep : checkpoint_keep;
            while (min-- > 0)
            {
              cp_filename_list[min] = myGH->cp_filename_list[min];
            }
            free (myGH->cp_filename_list);
            myGH->cp_filename_list = cp_filename_list;
            myGH->checkpoint_keep = checkpoint_keep;
          }
        }
      }
    }
  }

  /* stop the CP_TOTAL_TIMER timer and print timing information */
  if (myGH->print_timing_info)
  {
    CCTK_TimerStopI (myGH->timers[CP_TOTAL_TIMER]);
    IOUtil_PrintTimings ("Timing information for checkpointing with IOHDF5:",
                         3, myGH->timers, timer_descriptions);
  }

  /* save the iteration number of this checkpoint */
  myGH->last_checkpoint_iteration = GH->cctk_iteration;

  /* free allocated resources */
  free (tempname);
  free (filename);

  return (retval);
}
