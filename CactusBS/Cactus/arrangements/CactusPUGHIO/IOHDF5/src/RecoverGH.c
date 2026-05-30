 /*@@
   @file      RecoverGH.c
   @date      Fri Jun 19 09:14:22 1998
   @author    Tom Goodale
   @desc
              Contains the routines to recover from a HDF5 checkpoint file.

              Currently can recover from:
                (1) One file containing recombined data
                (2) Multiple unrecombined files, where the current
                    number of processors and IO processors
                    match those used to write the files.
   @enddesc
   @version   $Id$
 @@*/


#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "CactusPUGH/PUGH/src/include/pugh.h"
#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"
#include "ioHDF5GH.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusPUGHIO_IOHDF5_RecoverGH_c)


/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/
int IOHDF5_RecoverParameters (void);


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
static int OpenFile (cGH *GH, const char *basefilename, int called_from,
                     fileinfo_t *fileinfo);


 /*@@
   @routine    IOHDF5_Recover
   @date       Fri Jun 19 09:22:52 1998
   @author     Tom Goodale
   @desc
               Recovers a GH from an HDF5 file.
               This routine is registered with IOUtil
               as IOHDF5's recovery routine.
   @enddesc
   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      cGH
   @vio        in
   @endvar
   @var        basefilename
   @vdesc      the basefilename of the file to recover from
               The file suffix is appended by the routine.
   @vtype      const char *
   @vio        in
   @endvar
   @var        called_from
   @vdesc      flag indicating where this routine was called from
               (either RECOVER or filereader)
   @vtype      int
   @vio        in
   @endvar

   @calls      OpenFile
               IOHDF5Util_RecoverParameters
               IOHDF5Util_RecoverVariables
               IOHDF5Util_RecoverGHextensions
               CCTK_TimerStartI
               CCTK_TimerStopI
               IOUtil_PrintTimings

   @returntype int
   @returndesc
               -1 if there is no valid HDF5 file,
               or the returncode of
                 @seeroutine IOHDF5Util_RecoverParameters or
                 @seeroutine IOHDF5Util_RecoverVariables or
                 @seeroutine IOHDF5Util_RecoverGHextensions
   @endreturndesc
@@*/
int IOHDF5_Recover (cGH *GH, const char *basefilename, int called_from)
{
  int result;
  ioHDF5GH *myGH;
  static fileinfo_t fileinfo;  /* this is static because info is passed from
                                  CP_RECOVERY_PARAMETERS to CP_RECOVERY_DATA */
  const char *timer_description = "Time to recover:";
  DECLARE_CCTK_PARAMETERS


  result = 0;

  /* start the recovery timer if we were called at CCTK_RECOVER */
  myGH = (ioHDF5GH *) CCTK_GHExtension (GH, "IOHDF5");
  if (myGH && myGH->print_timing_info)
  {
    CCTK_TimerStartI (myGH->timers[RECOVERY_TIMER]);
  }

  /* open the file if it wasn't already opened at CCTK_RECOVER_PARAMETERS */
  /* FIXME Gab ... asymmetric levfac */
  if (called_from == CP_RECOVER_PARAMETERS ||
      called_from == FILEREADER_DATA ||
      (GH && (GH->cctk_levfac[0] > 1 || GH->cctk_convlevel > 0)))
  {
    if (OpenFile (GH, basefilename, called_from, &fileinfo) < 0)
    {
      return (-1);
    }
  }
  else
  {
    /* This is the case for CP_RECOVER_DATA.
       CCTK_RECOVER_PARAMETERS must have been called before
       and set up the file info structure. */
    if (! fileinfo.is_HDF5_file)
    {
      return (-1);
    }
  }

  /* if called at CCTK_RECOVER_PARAMETERS
     just do this and return (keeping the file open) */
  if (called_from == CP_RECOVER_PARAMETERS)
  {
    return (IOHDF5Util_RecoverParameters (&fileinfo));
  }

  /* Recover GH extensions */
  if (called_from == CP_RECOVER_DATA)
  {
    if (CCTK_Equals (verbose, "full"))
    {
      CCTK_INFO ("Recovering GH extensions");
    }
    result = IOHDF5Util_RecoverGHextensions (GH, &fileinfo);
  }

  if (! result)
  {
    /* Recover variables */
    if (CCTK_Equals (verbose, "full"))
    {
      CCTK_VInfo (CCTK_THORNSTRING, "Recovering %schunked data with ioproc %d, "
                  "ioproc_every %d", fileinfo.unchunked ? "un" : "",
                  fileinfo.ioproc, fileinfo.ioproc_every);
    }
    result = IOHDF5Util_RecoverVariables (GH, &fileinfo);
  }

  /* Close the file and remove it if requested by the user */
  if (CCTK_MyProc (GH) == fileinfo.ioproc)
  {
    if (CCTK_Equals (verbose, "full"))
    {
      if (called_from == CP_RECOVER_DATA)
      {
        CCTK_VInfo (CCTK_THORNSTRING, "Closing checkpoint file '%s' after "
                                      "successful recovery", fileinfo.filename);
      }
      else
      {
        CCTK_VInfo (CCTK_THORNSTRING, "Closing data file '%s'",
                    fileinfo.filename);
      }
    }
    HDF5_ERROR (H5Fclose (fileinfo.file));

    if (called_from == CP_RECOVER_DATA && recover_and_remove)
    {
      if (CCTK_Equals (verbose, "full"))
      {
        CCTK_VInfo (CCTK_THORNSTRING, "Old checkpoint file '%s' will be removed"
                                      " after next IO::checkpoint_keep "
                                      "successful checkpoints",
                                      fileinfo.filename);
      }
      myGH->cp_filename_list[myGH->cp_filename_index] =
        strdup (fileinfo.filename);
      myGH->cp_filename_index = (myGH->cp_filename_index+1) % checkpoint_keep;
    }
  }

  /* free the allocated filename */
  free (fileinfo.filename);
  fileinfo.filename = NULL;

  /* stop total recovery timer and print timing info */
  if (called_from == CP_RECOVER_DATA && myGH->print_timing_info)
  {
    CCTK_TimerStopI (myGH->timers[RECOVERY_TIMER]);
    IOUtil_PrintTimings ("Timing information for recovery in IOHDF5:",
                         1, &myGH->timers[RECOVERY_TIMER], &timer_description);
  }

  /* print an info message */
  if (called_from == CP_RECOVER_DATA)
  {
    CCTK_VInfo (CCTK_THORNSTRING,
                "Restarting simulation at iteration %d (physical time %g)",
                GH->cctk_iteration, (double) GH->cctk_time);
  }

  return (result);
}


 /*@@
   @routine    IOHDF5_RecoverParameters
   @date       Thu Apr 13 2000
   @author     Thomas Radke
   @desc
               Recovers the parameters from an HDF5 checkpoint file.
               This routine is scheduled at CCTK_RECOVER_PARAMETERS.

               Note that it cannot be registered with IOUtil to be scheduled
               from there (as done with the IOHDF5_Recover routine) because
               the registration mechanism isn't activated yet
               at CCTK_RECOVER_PARAMETERS.
               Instead we call the generic parameter recovery routine
               from IOUtil here, and just pass the necessary callback function
               and its arguments.

               Note also that this routine doesn't get passed any parameters,
               not even a GH, because this doesn't exist yet at the time it is
               being called.
   @enddesc

   @calls      IOUtil_RecoverParameters

   @returntype int
   @returndesc
               returncode of @seeroutine IOUtil_RecoverParameters
   @endreturndesc
@@*/
int IOHDF5_RecoverParameters (void)
{
  return (IOUtil_RecoverParameters (IOHDF5_Recover, ".h5", "HDF5"));
}


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
 /*@@
   @routine    OpenFile
   @date       Tue Oct 10 2000
   @author     Thomas Radke
   @desc
               Open a HDF5 file given by its basefilename.
               The basefilename is expanded into a full filename by calling
               IOUtil_AssembleFilename(). Both chunked and unchunked
               filenames are tested.
               If a file of that name could be opened it checks whether
               we can recover from that file.
               The file information is then broadcasted by the IO processor(s)
               to all other processors so that everyone knows how to proceed.
   @enddesc

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      cGH *
   @vio        in
   @endvar
   @var        basefilename
   @vdesc      basefilename of the HDF5 file to recover from
               For streamed files this should be of format 'host:port'.
   @vtype      int
   @vio        in
   @endvar
   @var        called_from
   @vdesc      flag indicating where this routine was called from
   @vtype      int
   @vio        in
   @endvar
   @var        fileinfo
   @vdesc      pointer to structure describing the file
   @vtype      fileinfo_t *
   @vio        out
   @endvar

   @returntype int
   @returndesc
               0 for success, -1 if file could not be opened
   @endreturndesc
@@*/
static int OpenFile (cGH *GH, const char *basefilename, int called_from,
                     fileinfo_t *fileinfo)
{
  int nprocs, myproc;
  hid_t group, version_attr;
  char *filename;
#ifdef CCTK_MPI
  MPI_Comm comm;
  CCTK_INT4 info[4];
#endif
  DECLARE_CCTK_PARAMETERS


#ifdef CCTK_MPI
  /* Get the communicator for broadcasting the info structure */
  /* NOTE: When recovering parameters thorn PUGH is not yet initialized
           so that we have to use MPI_COMM_WORLD in this case */
  comm = CCTK_GHExtensionHandle ("PUGH") < 0 ?
           MPI_COMM_WORLD : PUGH_pGH (GH)->PUGH_COMM_WORLD;
#endif

  /* identify myself */
  nprocs = CCTK_nProcs (GH);
  myproc = CCTK_MyProc (GH);

  /* Examine basefile to find out whether we are recovering from
   * one or multiple files and whether the data are chunked or not.
   *
   * This is done by processor 0 only since this is always an IO processor
   * and a corresponding file must exist in all cases.
   */

  /* Determine name of base file
     NOTE: As we don't know whether the file is chunked or not
           we need to try both file names. */
  /* at first try with unchunked mode */
  fileinfo->unchunked = 1;
  filename = IOUtil_AssembleFilename (GH, basefilename, "", ".h5", called_from,
                                      0, fileinfo->unchunked);

  if (myproc == 0)
  {
    if (CCTK_Equals (verbose, "full"))
    {
      CCTK_VInfo (CCTK_THORNSTRING, "Opening file '%s'", filename);
    }

    /* turn automatic error printing off again during check */
    H5E_BEGIN_TRY
    {
      /* Check the filetype */
      if (H5Fis_hdf5 (filename) > 0)
      {
        fileinfo->is_HDF5_file = 1;
      }
      else
      {
        CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "No valid HDF5 file '%s' found", filename);

        /* now try with chunked mode */
        fileinfo->unchunked = 0;
        free (filename);
        filename = IOUtil_AssembleFilename (GH, basefilename, "", ".h5",
                                            called_from, 0,fileinfo->unchunked);

        if (CCTK_Equals (verbose, "full"))
        {
          CCTK_VInfo (CCTK_THORNSTRING, "Trying now file '%s'...", filename);
        }

        fileinfo->is_HDF5_file = H5Fis_hdf5 (filename) > 0;
      }
    } H5E_END_TRY;
  }

  /* Okay, we have the complete filename. Let's read the file now. */
  if (myproc == 0 && fileinfo->is_HDF5_file)
  {
    fileinfo->is_HDF5_file = 0;

    fileinfo->file = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fileinfo->file < 0)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Cannot open file '%s' !", filename);
    }
    else if ((group = H5Gopen (fileinfo->file, GLOBAL_ATTRIBUTES_GROUP)) < 0)
    {
      CCTK_WARN (1, "Can't find global attributes group. "
                    "Is this really a Cactus HDF5 datafile ?");
    }
    else
    {
       /* Determine how the data was written */
      READ_ATTRIBUTE (group, "nprocs", H5T_NATIVE_INT, &fileinfo->nprocs);
      READ_ATTRIBUTE (group, "unchunked", H5T_NATIVE_INT, &fileinfo->unchunked);
      READ_ATTRIBUTE (group, "ioproc_every", H5T_NATIVE_INT,
                      &fileinfo->ioproc_every);

      /* check if there exists a version attribute
         For this we temporarily turn off automatic error printing. */
      H5E_BEGIN_TRY
      {
        version_attr = H5Aopen_name (group, "Cactus version");
      } H5E_END_TRY;

      fileinfo->has_version = version_attr >= 0;
      if (version_attr >= 0)
      {
        HDF5_ERROR (H5Aclose (version_attr));
      }

      HDF5_ERROR (H5Gclose (group));

      /* If we recover from multiple files the number of
       * writing processors must match the number of reading
       * processors, and the total number of processors must match.
       */
      if ((fileinfo->ioproc_every == nprocs && nprocs > 1) ||
          fileinfo->unchunked)
      {
        if (! CCTK_Equals (verbose, "none"))
        {
          CCTK_VInfo (CCTK_THORNSTRING, "Recovering from one %schunked file",
                                        fileinfo->unchunked ? "un" : "");
        }
        fileinfo->ioproc_every = nprocs;
        fileinfo->is_HDF5_file = 1;
      }
      else
      {
        if (fileinfo->nprocs != nprocs)
        {
          CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Must restart on %d processors with chunked files "
                      "or recombine them", fileinfo->nprocs);
        }
        else
        {
          if (! CCTK_Equals (verbose, "none"))
          {
            CCTK_VInfo (CCTK_THORNSTRING, "Recovering from %d chunked files",
                        nprocs / fileinfo->ioproc_every +
                        (nprocs % fileinfo->ioproc_every ? 1 : 0));
          }
          fileinfo->is_HDF5_file = 1;
        }
      }
    }
  }

  if (myproc == 0 && ! fileinfo->is_HDF5_file)
  {
    CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                "No valid HDF5 file '%s' found", filename);
  }

#ifdef CCTK_MPI
  /* Broadcast the file information to all processors
     Need to convert everything into CCTK_INTs which can be communicated. */
  info[0] = fileinfo->is_HDF5_file;
  info[1] = fileinfo->unchunked;
  info[2] = fileinfo->ioproc_every;
  info[3] = fileinfo->has_version;
  CACTUS_MPI_ERROR (MPI_Bcast (info, 4, PUGH_MPI_INT4, 0, comm));
  fileinfo->is_HDF5_file = info[0];
  fileinfo->unchunked = info[1];
  fileinfo->ioproc_every = info[2];
  fileinfo->has_version = info[3];
#endif

  if (fileinfo->is_HDF5_file)
  {

    /* Determine the IO processors for each node and the corresponding
       checkpoint file */
    fileinfo->ioproc = myproc - (myproc % fileinfo->ioproc_every);
    filename = IOUtil_AssembleFilename (GH, basefilename, "", ".h5",called_from,
                                        fileinfo->ioproc/fileinfo->ioproc_every,
                                        fileinfo->unchunked);

    /* Open chunked files on other IO processors */
    if (myproc == fileinfo->ioproc && myproc != 0)
    {
      if (CCTK_Equals (verbose, "full"))
      {
        CCTK_VInfo (CCTK_THORNSTRING, "Opening chunked file '%s' on "
                                      "processor %d", filename, myproc);
      }

      /* Open file, make sure the file is valid */
      fileinfo->file = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
      if (fileinfo->file < 0)
      {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Cannot open file '%s' on processor %d", filename, myproc);
        fileinfo->is_HDF5_file = 0;
      }
    }
#ifdef CCTK_MPI
    /* Finally check whether all files have valid recovery files */
    info[0] = fileinfo->is_HDF5_file;
    CACTUS_MPI_ERROR (MPI_Allreduce (&info[0], &info[1], 1,
                                     PUGH_MPI_INT4, MPI_LAND, comm));
    fileinfo->is_HDF5_file = info[1];
#endif
  }

  /* set the filename in the info structure */
  fileinfo->filename = filename;

  /* return 0 for success otherwise negative */
  return (fileinfo->is_HDF5_file ? 0 : -1);
}
