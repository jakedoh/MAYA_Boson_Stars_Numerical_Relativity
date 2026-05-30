/*@@
   @file      Write.c
   @date      Fri May 21 1999
   @author    Thomas Radke, Gerd Lanfermann
   @desc
              File handling routines for IOHDF5.
   @enddesc
   @version   $Id$
 @@*/


#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "StoreNamedData.h"
#include "CactusBase/IOUtil/src/ioGH.h"
#include "CactusBase/IOUtil/src/ioutil_AdvertisedFiles.h"

#include "ioHDF5GH.h"
#if defined(CCTK_MPI) && defined(H5_HAVE_PARALLEL)
#include "CactusPUGH/PUGH/src/include/pugh.h"   /* PUGH_COMM_WORLD */
#endif


/* the rcs ID and its dummy funtion to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusPUGHIO_IOHDF5_Write_c)


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
static char *IOHDF5_GetFilename (const cGH *GH,
                                 int vindex,
                                 const char *name,
                                 int *is_new_file);


/*@@
   @routine    IOHDF5_Write
   @date       Tue Oct 10 2000
   @author     Thomas Radke
   @desc
               Opens the HDF5 output file, calls the dump routine,
               and closes it again.
   @enddesc

   @calls      IOHDF5Util_DumpVar

   @var        GH
   @vdesc      Pointer to CCTK GH
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        vindex
   @vdesc      index of variable
   @vtype      int
   @vio        in
   @endvar
   @var        alias
   @vdesc      alias name of variable to output
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine IOHDF5Util_DumpVar, or<BR>
               -1 if no storage was assigned to variable<BR>
               -2 if filename couldn't be built from alias
   @endreturndesc
@@*/
int IOHDF5_Write (const cGH *GH, int vindex, const char *alias)
{
  const ioGH *ioUtilGH;
  ioHDF5GH *myGH;
  hid_t file, plist;
  int myproc, is_new_file, retval;
  char *filename, *fullname;
  DECLARE_CCTK_PARAMETERS


  ioUtilGH = CCTK_GHExtension (GH, "IO");
  myGH     = CCTK_GHExtension (GH, "IOHDF5");

  /* check if variable has storage assigned */
  if (! CCTK_QueryGroupStorageI (GH, CCTK_GroupIndexFromVarI (vindex)))
  {
    fullname = CCTK_FullName (vindex);
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "No IOHDF5 output for '%s' (no storage)", fullname);
    free (fullname);
    return (-1);
  }

  /* get the filename for output */
  filename = IOHDF5_GetFilename (GH, vindex, alias, &is_new_file);
  if (! filename)
  {
    fullname = CCTK_FullName (vindex);
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Unable to construct file name for '%s'", fullname);
    free (fullname);
    return (-2);
  }

  /* open the output file on all I/O processors or - in case of parallel HDF5 -
     on every processor */
  file = -1;
  myproc = CCTK_MyProc (GH);
#ifdef H5_HAVE_PARALLEL
  if ((myproc == ioUtilGH->ioproc || ioUtilGH->unchunked) &&
      (CCTK_GroupTypeFromVarI (vindex) != CCTK_SCALAR || myproc == 0))
#else
  if (myproc == ioUtilGH->ioproc &&
      (CCTK_GroupTypeFromVarI (vindex) != CCTK_SCALAR || myproc == 0))
#endif
  {
    if (CCTK_Equals (verbose, "full"))
    {
      CCTK_VInfo (CCTK_THORNSTRING, "%s HDF5 output file '%s'",
                  is_new_file ? "Opening" : "Appending", filename);
    }

    HDF5_ERROR (plist = H5Pcreate (H5P_FILE_ACCESS));
#if defined(CCTK_MPI) && defined(H5_HAVE_PARALLEL)
    if (ioUtilGH->unchunked)
    {
      HDF5_ERROR (H5Pset_fapl_mpio (plist, PUGH_pGH (GH)->PUGH_COMM_WORLD,
                                    MPI_INFO_NULL));
    }
#endif

    if (is_new_file)
    {
      HDF5_ERROR (file = H5Fcreate (filename, H5F_ACC_TRUNC, H5P_DEFAULT,
                                      plist));
    }
    else
    {
      HDF5_ERROR (file = H5Fopen (filename, H5F_ACC_RDWR, plist));
    }
    if (file < 0)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Couldn't %s IOHDF5 output file '%s'",
                  is_new_file ? "create" : "open", filename);
    }

    HDF5_ERROR (H5Pclose (plist));
  }

  /* output global attributes */
  if (is_new_file && file >= 0)
  {
    /* all GH extension variables and parameter stuff which is not
       specific to any dataset goes into dedicated groups */
    if (! CCTK_Equals (out_save_parameters, "no"))
    {
      IOHDF5Util_DumpParameters (GH, CCTK_Equals (out_save_parameters, "all"),
                                 file);
    }

    IOHDF5Util_DumpGHExtensions (GH, file);
  }

  /* output the data */
  retval = IOHDF5Util_DumpVar (GH, myGH->requests[vindex], file);

  /* close the file */
  if (file >= 0)
  {
    if (CCTK_Equals (verbose, "full"))
    {
      CCTK_INFO ("Closing HDF5 output file from this iteration");
    }
    HDF5_ERROR (H5Fclose (file));
  }

  /* free the filename if it was not registered for re-opening */
  if (out_timesteps_per_file > 0)
  {
    free (filename);
  }

  return (retval);
}


/*@@
   @routine    IOHDF5_GetFilename
   @author     Paul Walker
   @date       Feb 1997
   @desc
               Builds the filename for output.
               The names of all output files are stored in the filename
               database.
               The routine first searches in this database if there is already
               registered a filename for alias. If so, this name is returned.
               Otherwise it builds the filename from alias and stores it in the
               database.
               Generalized for Hyperslab extraction (GL)
   @enddesc

   @var        GH
   @vdesc      Pointer to CCTK GH
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        vindex
   @vdesc      index of variable to output
   @vtype      int
   @vio        in
   @endvar
   @var        varname
   @vdesc      alias name of variable to output
   @vtype      const char *
   @vio        in
   @endvar
   @var        is_new_file
   @vdesc      flag indicating if the file is to be created
   @vtype      int
   @vio        out
   @endvar

   @returntype char *
   @returndesc
               the allocated string containing the filename
   @endreturndesc
@@*/
static char *IOHDF5_GetFilename (const cGH *GH,
                                 int vindex,
                                 const char *varname,
                                 int *is_new_file)
{
  const ioGH *ioUtilGH;
  ioHDF5GH *myGH;
  int myproc, vdim, result;
  /* FIXME: make these strings dynamic */
  char extra[256], extradir[256];
  char *filename, *outputdir, *fullname, *tmp;
  ioAdvertisedFileDesc advertised_file;
  DECLARE_CCTK_PARAMETERS


  /* get GH extensions for IOUtil and IOHDF5 */
  ioUtilGH = CCTK_GHExtension (GH, "IO");
  myGH     = CCTK_GHExtension (GH, "IOHDF5");

  filename = GetNamedData (myGH->open_output_files, varname);
  if (filename)
  {
    /* set flag to indicate that file should be opened in append mode */
    *is_new_file = 0;
    return (filename);
  }

  extra[0] = extradir[0]= '\0';

  myproc = CCTK_MyProc (GH);
  vdim   = CCTK_GroupDimFromVarI (vindex);

  if (out_timesteps_per_file > 0)
  {
    tmp = extra;
    sprintf (extra, "%s.time_%7.3f", extra, GH->cctk_time);

    /* And be sure to replace any spaces in the filename with an _ */
    do
    {
      if (*tmp == ' ')
      {
        *tmp = '_';
      }
    } while (*++tmp);
  }


  /*
   *  OUTPUT ONE FILE FOR EACH N PROCESSORS
   *  -------------------------------------
   *
   *  If only one output file, the single file for each GV is written
   *  in the default output directory (ie. extradir = ""). Otherwise
   *  a directory is created for each output GV to hold the multiple files.
   */
  if (! ioUtilGH->unchunked && ioUtilGH->ioproc_every < CCTK_nProcs (GH) &&
      CCTK_GroupTypeFromVarI (vindex) != CCTK_SCALAR)
  {
    /* Add the output processor number to the extra string */
    sprintf (extra, "%s.file_%d", extra, myproc / ioUtilGH->ioproc_every);

    /* If necessary create the output directory */
    outputdir = malloc (strlen (myGH->out_dir) + strlen (varname) + 5);
    sprintf (outputdir, "%s%s_%dd", myGH->out_dir, varname, vdim);

    result = IOUtil_CreateDirectory (GH, outputdir,
                                     ! CCTK_Equals (out_mode, "onefile"),
                                     ioUtilGH->ioproc);
    if (result < 0)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Problem creating output directory '%s'", outputdir);
    }
    else if (result > 0 && CCTK_Equals (verbose, "full"))
    {
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Output directory '%s' already exists", outputdir);
    }

    free (outputdir);

#ifdef CCTK_MPI
    /* Wait for all processors to catch up */
    CCTK_Barrier (GH);
#endif

    /* extradir is the relative output directory */
    sprintf (extradir, "%s_%dd/", varname, vdim);
  }

  /* CREATE THE COMPLETE OUTPUT FILENAME
     ----------------------------------- */
  filename = malloc (8 + strlen (myGH->out_dir) + strlen (extradir) +
                     strlen (varname) + strlen (extra) + strlen (GH->identity));
  sprintf (filename, "%s%s%s%s%s.h5",
           myGH->out_dir, extradir, varname, extra, GH->identity);

  /* no need to store filenames if used only once */
  *is_new_file = 1;
  if (out_timesteps_per_file < 0)
  {
    if (myproc == ioUtilGH->ioproc)
    {
      if (! IO_TruncateOutputFiles (GH))
      {
        H5E_BEGIN_TRY
        {
          *is_new_file = H5Fis_hdf5 (filename) <= 0;
        } H5E_END_TRY;
      }
    }
    StoreNamedData (&myGH->open_output_files, varname, filename);
  }

  /* advertise the file for downloading */
  fullname = CCTK_FullName (vindex);
  advertised_file.slice = "";
  advertised_file.thorn = CCTK_THORNSTRING;
  advertised_file.varname = fullname;
  advertised_file.description = "Full-dimensional variable contents";
  advertised_file.mimetype = "data/x-hdf5";

  IOUtil_AdvertiseFile (GH, filename, &advertised_file);

  free (fullname);

  return (filename);
}
