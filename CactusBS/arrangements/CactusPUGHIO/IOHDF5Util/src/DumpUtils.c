 /*@@
   @file      DumpUtils.c
   @date      Fri Oct 6 2000
   @author    Thomas Radke
   @desc
              Utility routines for dumping data into HDF5 files.
   @enddesc
   @version   $Id$
 @@*/


#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Version.h"
#include "cctk_Parameters.h"
#include "util_Table.h"
#include "CactusBase/IOUtil/src/ioGH.h"
#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"
#include "ioHDF5UtilGH.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusPUGHIO_IOHDF5Util_DumpUtils_c)


/*@@
   @routine    IOHDF5Util_DumpGH
   @date       Tue Oct 10 2000
   @author     Thomas Radke
   @desc
               Dump all variables of the given GH along with the global
               attributes and parameters to the given HDF5 checkpoint file.
   @enddesc
   @var        GH
   @vdesc      pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        timers
   @vdesc      array of timers for timing the checkpointing
   @vtype      const int *
   @vio        in
   @endvar
   @var        file
   @vdesc      HDF5 file handle
   @vtype      hid_t
   @vio        in
   @endvar

   @returntype int
   @returndesc
               the accumulated return code from @seeroutine IOHDF5Util_DumpVar
               (should be 0 for success, or negative if some error occured)
   @endreturndesc
@@*/
int IOHDF5Util_DumpGH (const cGH *GH, const int *timers, hid_t file)
{
  int len, dim, groupvarsize, first_vindex, gindex, retval;
  cGroup gdata;
  cGroupDynamicData gdynamicdata;
  char *fullname;
  ioRequest *request;
  DECLARE_CCTK_PARAMETERS


  retval = 0;

  /* start CP_PARAMETERS_TIMER timer */
  if (timers)
  {
    CCTK_TimerResetI (timers[CP_PARAMETERS_TIMER]);
    CCTK_TimerStartI (timers[CP_PARAMETERS_TIMER]);
  }

  /* now start dumping away */
  if (file >= 0)
  {
    /* all GH extension variables and parameter stuff which is not
       specific to any dataset goes into dedicated groups */
    IOHDF5Util_DumpParameters (GH, 1, file);
    IOHDF5Util_DumpGHExtensions (GH, file);
  }

  if (CCTK_Equals (verbose, "full"))
  {
    CCTK_INFO ("Dumping Grid Variables ...");
  }

  /* stop CP_PARAMETERS_TIMER timer and start CP_VARIABLES_TIMER */
  if (timers)
  {
    CCTK_TimerStopI  (timers[CP_PARAMETERS_TIMER]);
    CCTK_TimerResetI (timers[CP_VARIABLES_TIMER]);
    CCTK_TimerStartI (timers[CP_VARIABLES_TIMER]);
  }

  /* ... now the variables, sorted by groups */
  for (gindex = CCTK_NumGroups () - 1; gindex >= 0; gindex--)
  {
    /* skip empty groups */
    if (CCTK_NumVarsInGroupI (gindex) <= 0)
    {
      continue;
    }

    /* only dump groups which have storage assigned */
    if (CCTK_QueryGroupStorageI (GH, gindex) <= 0)
    {
      continue;
    }

    /* get the number of allocated timelevels */
    CCTK_GroupData (gindex, &gdata);
    gdata.numtimelevels = 0;
    gdata.numtimelevels = CCTK_GroupStorageIncrease (GH, 1, &gindex,
                                                     &gdata.numtimelevels,NULL);

    /* dump all timelevels except the oldest (for multi-level groups) */
    if (gdata.numtimelevels > 1)
    {
      gdata.numtimelevels--;
    }

    /* skip zero-sized groups */
    CCTK_GroupDynamicData (GH, gindex, &gdynamicdata);
    groupvarsize = 1;
    for (dim = 0; dim < gdynamicdata.dim; dim++)
    {
      groupvarsize *= gdynamicdata.gsh[dim];
    }
    if (groupvarsize <= 0)
    {
      continue;
    }

    len = Util_TableGetString (gdata.tagstable, 0, NULL, "checkpoint");
    if (len > 0)
    {
      char* value = malloc (len + 1);
      Util_TableGetString (gdata.tagstable, len + 1, value, "checkpoint");
      if (len == sizeof ("no") - 1 && CCTK_Equals (value, "no"))
      {
        continue;
      }
      else if (! CCTK_Equals (value, "yes"))
      {
        char* groupname = CCTK_GroupName (gindex);
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Ignoring unknown checkpoint tag '%s' for group '%s'",
                    value, groupname);
        free (groupname);
      }
      free (value);
    }

    first_vindex = CCTK_FirstVarIndexI (gindex);

    /* get the default I/O request for this group */
    request = IOUtil_DefaultIORequest (GH, first_vindex, 1, -1.0);

    /* disable checking for old data objects,
       disable datatype conversion and downsampling,
       request hyperslab output including ghostzones */
    request->check_exist = 0;
    request->hdatatype = gdata.vartype;
    for (request->hdim = 0; request->hdim < request->vdim; request->hdim++)
    {
      request->downsample[request->hdim] = 1;
    }
    request->with_ghostzones = 1;

    /* loop over all variables in this group */
    for (request->vindex = first_vindex;
         request->vindex < first_vindex + gdata.numvars;
         request->vindex++)
    {
      /* loop over all allocated timelevels of this variable */
      for (request->timelevel = 0;
           request->timelevel < gdata.numtimelevels;
           request->timelevel++)
      {
        if (CCTK_Equals (verbose, "full") && file >= 0)
        {
          fullname = CCTK_FullName (request->vindex);
          CCTK_VInfo (CCTK_THORNSTRING, "  %s (timelevel %d)",
                      fullname, request->timelevel);
          free (fullname);
        }

        retval += IOHDF5Util_DumpVar (GH, request, file);
      }

    } /* end of loop over all variables */

    /* free I/O request */
    IOUtil_FreeIORequest (&request);

  } /* end of loop over all groups */

  /* stop CP_VARIABLES_TIMER timer */
  if (timers)
  {
    CCTK_TimerStopI (timers[CP_VARIABLES_TIMER]);
  }

  return (retval);
}


/*@@
   @routine    IOHDF5Util_DumpCommonAttributes
   @date       May 21 1999
   @author     Thomas Radke
   @desc
               Add "Common" attributes to the given dataset, these are:
               <ul>
                 <li> the variable's groupname
                 <li> the grouptype
                 <li> number of timelevels
                 <li> global grid size
                 <li> simulation time
                 <li> bounding box info (if the variable has a coordinate
                      system associated with it)
               </ul>
   @enddesc
   @var        GH
   @vdesc      pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        request
   @vdesc      reference to the I/O request description
   @vtype      const ioRequest *
   @vio        in
   @endvar
   @var        dataset
   @vdesc      the object handle where the attributes should be attached to
   @vtype      hid_t
   @vio        in
   @endvar
@@*/
void IOHDF5Util_DumpCommonAttributes (const cGH *GH, const ioRequest *request,
                                      hid_t object)
{
  int hdim, vdim, coord_system_handle;
  char *groupname;
  CCTK_INT iattr;
  CCTK_INT *coord_handles;
  CCTK_REAL *rattr;
  const ioHDF5UtilGH *myGH;
  DECLARE_CCTK_PARAMETERS


  myGH = CCTK_GHExtension (GH, "IOHDF5Util");

  /* attributes describing the variable */
  groupname = CCTK_GroupNameFromVarI (request->vindex);
  WRITE_ATTRIBUTE ("groupname", groupname, object, myGH, 0, myGH->HDF5_STRING);
  coord_system_handle = -1;
  if (CCTK_IsFunctionAliased ("Coord_GroupSystem"))
  {
    coord_system_handle = Coord_GroupSystem (GH, groupname);
  }
  free (groupname);
  iattr = CCTK_GroupTypeFromVarI (request->vindex);
  WRITE_ATTRIBUTE ("grouptype", &iattr, object, myGH, 0, HDF5_INT);
  iattr = CCTK_MaxActiveTimeLevelsVI (GH, request->vindex);
  WRITE_ATTRIBUTE ("ntimelevels", &iattr, object, myGH, 0, HDF5_INT);
  WRITE_ATTRIBUTE ("global_size", request->hsize, object, myGH, request->hdim,
                   HDF5_INT);
  WRITE_ATTRIBUTE ("time", &GH->cctk_time, object, myGH, 0, HDF5_REAL);

  /* write bbox attributes if we have coordinate system info */
  coord_handles = malloc (request->vdim * sizeof (CCTK_INT));
  if (coord_system_handle >= 0 &&
      Util_TableGetIntArray (coord_system_handle, request->vdim,
                             coord_handles, "COORDINATES") >= 0)
  {
    rattr = calloc (3 * request->vdim, sizeof (CCTK_REAL));

    for (vdim = 0; vdim < request->vdim; vdim++)
    {
      for (hdim = 0; hdim < request->hdim; hdim++)
      {
        if (request->direction[hdim*request->hdim + vdim])
        {
          Util_TableGetReal (coord_handles[vdim], &rattr[hdim+0*request->vdim],
                             "COMPMIN");
          if (Util_TableGetReal (coord_handles[vdim],
                                 &rattr[hdim+2*request->vdim], "DELTA") >= 0)
          {
            rattr[hdim+2*request->vdim] *= request->downsample[hdim];
            rattr[hdim+1*request->vdim]  = rattr[hdim+0*request->vdim];
            rattr[hdim+1*request->vdim] +=
              ((request->extent[hdim] + request->downsample[hdim]-1)
               / request->downsample[hdim] - 1)
              * rattr[hdim+2*request->vdim];
          }
        }
      }
    }

    WRITE_ATTRIBUTE ("origin", rattr + 0*vdim, object, myGH, vdim, HDF5_REAL);
    WRITE_ATTRIBUTE ("min_ext", rattr + 0*vdim, object, myGH, vdim, HDF5_REAL);
    WRITE_ATTRIBUTE ("max_ext", rattr + 1*vdim, object, myGH, vdim, HDF5_REAL);
    WRITE_ATTRIBUTE ("delta", rattr + 2*vdim, object, myGH, vdim, HDF5_REAL);

    free (rattr);
  }

  free (coord_handles);
}


 /*@@
   @routine    IOHDF5Util_DumpParameters
   @date       Thu Jan 20 2000
   @author     Thomas Radke
   @desc
               Collects the parameters of all active implementations
               into a single string and writes it as a single dataset
               into the CACTUS_PARAMETERS_GROUP group in the HDF5 file.

               Note that we used to write the parameters string as a single
               attribute attached to the CACTUS_PARAMETERS_GROUP group.
               This caused problems with very long strings (HDF5 has a 16k-limit
               on the total size of attributes).
   @enddesc
   @var        GH
   @vdesc      pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        all
   @vdesc      flag indicating whether to save all parameters or just the ones
               which have been set before
   @vtype      int
   @vio        in
   @endvar
   @var        file
   @vdesc      the HDF5 file handle
   @vtype      hid_t
   @vio        in
   @endvar
@@*/
void IOHDF5Util_DumpParameters (const cGH *GH, int all, hid_t file)
{
  char *parameters;
  hid_t group, dataspace, dataset;
  hsize_t size;
  DECLARE_CCTK_PARAMETERS


  if (CCTK_Equals (verbose, "full"))
  {
    CCTK_INFO ("Dumping Parameters ...");
  }

  /* get the parameter string buffer */
  parameters = IOUtil_GetAllParameters (GH, all);

  if (parameters)
  {
    size = strlen (parameters) + 1;
    HDF5_ERROR (group = H5Gcreate (file, CACTUS_PARAMETERS_GROUP, 0));
    HDF5_ERROR (dataspace = H5Screate_simple (1, &size, NULL));
    HDF5_ERROR (dataset = H5Dcreate (group, ALL_PARAMETERS, H5T_NATIVE_UCHAR,
                                     dataspace, H5P_DEFAULT));
    HDF5_ERROR (H5Dwrite (dataset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, parameters));
    HDF5_ERROR (H5Dclose (dataset));
    HDF5_ERROR (H5Sclose (dataspace));
    HDF5_ERROR (H5Gclose (group));

    free (parameters);
  }
}


 /*@@
   @routine    IOHDF5Util_DumpGHExtensions
   @date       Thu Jan 20 2000
   @author     Thomas Radke
   @desc
               Attaches the current values of important flesh and IO variables
               as attributes to the GLOBAL_ATTRIBUTES_GROUP group
               in the HDF5 file.
   @enddesc
   @var        GH
   @vdesc      pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        file
   @vdesc      the HDF5 file handle
   @vtype      hid_t
   @vio        in
   @endvar
@@*/
void IOHDF5Util_DumpGHExtensions (const cGH *GH, hid_t file)
{
  int value;
  hid_t group;
  char buffer[128];
  const char *version;
  const ioGH *ioUtilGH;
  const ioHDF5UtilGH *myGH;
  DECLARE_CCTK_PARAMETERS


  if (CCTK_Equals (verbose, "full"))
  {
    CCTK_INFO ("Dumping GH extensions ...");
  }

  ioUtilGH = CCTK_GHExtension (GH, "IO");
  myGH = CCTK_GHExtension (GH, "IOHDF5Util");

  HDF5_ERROR (group = H5Gcreate (file, GLOBAL_ATTRIBUTES_GROUP, 0));

  value = CCTK_MainLoopIndex ();
  WRITE_ATTRIBUTE ("main_loop_index", &value, group, myGH, 0, H5T_NATIVE_INT);
  value = CCTK_nProcs (GH);
  WRITE_ATTRIBUTE ("nprocs", &value, group, myGH, 0, H5T_NATIVE_INT);
  WRITE_ATTRIBUTE ("ioproc_every", &ioUtilGH->ioproc_every, group, myGH, 0,
                   H5T_NATIVE_INT);
  WRITE_ATTRIBUTE ("unchunked", &ioUtilGH->unchunked, group, myGH, 0,
                   H5T_NATIVE_INT);
  WRITE_ATTRIBUTE ("cctk_time", &GH->cctk_time, group, myGH, 0, HDF5_REAL);
  WRITE_ATTRIBUTE ("cctk_delta_time", &GH->cctk_delta_time, group, myGH, 0,
                   HDF5_REAL);
  WRITE_ATTRIBUTE ("cctk_iteration", &GH->cctk_iteration, group, myGH, 0,
                   H5T_NATIVE_INT);
  version = CCTK_FullVersion ();
  WRITE_ATTRIBUTE ("Cactus version", version, group, myGH, 0,myGH->HDF5_STRING);

  /* add the parameter filename and the creation date
     as file identification attributes */
  if (CCTK_Equals (out_fileinfo, "parameter filename") ||
      CCTK_Equals (out_fileinfo, "all"))
  {
    buffer[0] = 0;
    CCTK_ParameterFilename (sizeof (buffer), buffer);
    WRITE_ATTRIBUTE ("parameter file", buffer, group, myGH,0,myGH->HDF5_STRING);
  }
  if (CCTK_Equals (out_fileinfo, "creation date") ||
      CCTK_Equals (out_fileinfo, "all"))
  {
    buffer[0] = 0;
    Util_CurrentDate (sizeof (buffer), buffer);
    value = strlen (buffer) + 1;
    buffer[value-1] = ' ';
    Util_CurrentTime (sizeof (buffer) - value, buffer + value);
    WRITE_ATTRIBUTE ("creation date", buffer, group, myGH, 0,myGH->HDF5_STRING);
  }

  HDF5_ERROR (H5Gclose (group));
}


/*@@
   @routine    IOHDF5Util_DataType
   @author     Thomas Radke
   @date       Mon 11 June 2001
   @desc
               Returns the HDF5 datatype for a given CCTK datatype
   @enddesc

   @var        myGH
   @vdesc      pointer to IOHDF5Util's GH extensions
   @vtype      ioHDF5UtilGH *
   @vio        in
   @endvar
   @var        cctk_type
   @vdesc      CCTK datatype
   @vtype      int
   @vio        in
   @endvar

   @returntype hid_t
   @returndesc
               the appropriate HDF5 datatype for success, or -1 otherwise
   @endreturndesc
@@*/
hid_t IOHDF5Util_DataType (const ioHDF5UtilGH *myGH, int cctk_type)
{
  hid_t retval;


  switch (cctk_type)
  {
    case CCTK_VARIABLE_BYTE:      retval = HDF5_BYTE; break;
    case CCTK_VARIABLE_INT:       retval = HDF5_INT; break;
    case CCTK_VARIABLE_REAL:      retval = HDF5_REAL; break;
    case CCTK_VARIABLE_COMPLEX:   retval = myGH->HDF5_COMPLEX; break;
#ifdef HAVE_CCTK_INT1
    case CCTK_VARIABLE_INT1:      retval = HDF5_INT1; break;
#endif
#ifdef HAVE_CCTK_INT2
    case CCTK_VARIABLE_INT2:      retval = HDF5_INT2; break;
#endif
#ifdef HAVE_CCTK_INT4
    case CCTK_VARIABLE_INT4:      retval = HDF5_INT4; break;
#endif
#ifdef HAVE_CCTK_INT8
    case CCTK_VARIABLE_INT8:      retval = HDF5_INT8; break;
#endif
#ifdef HAVE_CCTK_REAL4
    case CCTK_VARIABLE_REAL4:     retval = HDF5_REAL4; break;
    case CCTK_VARIABLE_COMPLEX8:  retval = myGH->HDF5_COMPLEX8; break;
#endif
#ifdef HAVE_CCTK_REAL8
    case CCTK_VARIABLE_REAL8:     retval = HDF5_REAL8; break;
    case CCTK_VARIABLE_COMPLEX16: retval = myGH->HDF5_COMPLEX16; break;
#endif
#ifdef HAVE_CCTK_REAL16
    case CCTK_VARIABLE_REAL16:    retval = HDF5_REAL16; break;
    case CCTK_VARIABLE_COMPLEX32: retval = myGH->HDF5_COMPLEX32; break;
#endif

    default: CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                         "Unsupported CCTK variable datatype %d", cctk_type);
             retval = -1;
  }

  return (retval);
}
