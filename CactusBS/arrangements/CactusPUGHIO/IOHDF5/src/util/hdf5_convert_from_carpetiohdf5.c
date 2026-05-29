 /*@@
   @file      hdf5_convert_from_carpetiohdf5.c
   @date      Tue 11 October 2005
   @author    Thomas Radke
   @desc
              Utility program to convert CarpetIOHDF5 datafiles
              into IOHDF5 datafiles.
   @enddesc
   @version   $Id$
 @@*/

#define  MAXDIM       3
#define  MAXNAMESIZE  100

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/* Cactus includes (defines CCTK_FILEVERSION) */
#include "cctk.h"
#include "cctk_Version.h"

/* IOHDF5 includes */
#include "../../../IOHDF5Util/src/ioHDF5UtilGH.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusPUGHIO_IOHDF5_util_hdf5_convert_from_carpetiohdf5_c)


/*****************************************************************************/
/*                           macro definitions                               */
/*****************************************************************************/
/* uncomment the following to get some debugging output */
/* #define IOHDF5_DEBUG  1 */

/* metadata group in CarpetIOHDF5 files */
/* FIXME: should rather include the proper CarpetIOHDF5 header file */
#define METADATA_GROUP "Parameters and Global Attributes"

/* macro to do an HDF5 call, check its return code, and print a warning
   in case of an error */
#define CHECK_ERROR(hdf5_call)                                                \
          do                                                                  \
          {                                                                   \
            hid_t _error_code = hdf5_call;                                    \
                                                                              \
                                                                              \
            if (_error_code < 0)                                              \
            {                                                                 \
              fprintf (stderr, "WARNING: line %d: HDF5 call '%s' returned "   \
                               "error code %d\n",                             \
                                __LINE__, #hdf5_call, (int)_error_code);      \
              exit (-1);                                                      \
            }                                                                 \
          } while (0)

/*****************************************************************************/
/*                           global variables                                */
/*****************************************************************************/
static int   nprocs = 0;           /* number of processors for IOHDF5 chunked data */
static int   unchunked = 0;        /* flag indicating chunked/unchunked data */
static char *pathname = NULL;      /* pathname of the current object */
static int ntimelevels = 1;        /* value of the timelevel attribute */
static const char *groupname = NULL; /* value of the groupname attribute */

/*****************************************************************************/
/*                           local function prototypes                       */
/*****************************************************************************/
static herr_t CopyObject (hid_t copy_from, const char *objectname, void *arg);
static herr_t CopyAttribute (hid_t src, const char *attr_name, void *arg);



 /*@@
   @routine    main
   @date       Tue 11 October 2005
   @author     Thomas Radke
   @desc
               Main routine of the CarpetIOHDF5-to-IOHDF5 converter
   @enddesc

   @var        argc
   @vdesc      number of command line arguments
   @vtype      int
   @vio        in
   @endvar
   @var        argv
   @vdesc      command line arguments
   @vtype      char *[]
   @vio        in
   @endvar

   @returntype int
   @returndesc
               0 for success, negative return values indicate an error
   @endreturndesc
@@*/
int main (int argc, char **argv)
{
  hid_t infile, outfile, attr, attrDataspace;
  hid_t group, hdf5String;
  const char *cactus_version = "Cactus Beta16"; /* CCTK_FullVersion (); */


  if (argc != 6)
  {
    fprintf (stderr, "Usage: %s <groupname> <ntimelevels> <nprocs> <inputfile> <outputfile>\n", argv[0]);
    fprintf (stderr, "   eg. %s grid::coordinates 1 4 gxx-CarpetIOHDF5.h5 gxx-IOHDF5.h5\n", argv[0]);
    return (0);
  }

  groupname = argv[1];
  ntimelevels = atoi (argv[2]);
  nprocs = atoi (argv[3]);

  infile = H5Fopen (argv[4], H5F_ACC_RDONLY, H5P_DEFAULT);
  if (infile < 0)
  {
    fprintf (stderr, "Could not open CarpetIOHDF5 input file '%s'\n", argv[4]);
    return (-1);
  }

  outfile = H5Fcreate (argv[5], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (outfile < 0)
  {
    fprintf (stderr, "Could not create HDF5 output file '%s'\n", argv[5]);
    return (-1);
  }

  printf ("\n  ---------------------------------------\n"
            "  Cactus CarpetIOHDF5-to-IOHDF5 Converter\n"
            "  ---------------------------------------\n");

  CHECK_ERROR (hdf5String = H5Tcopy (H5T_C_S1));

#if 0
  /* read the metadata group from the CarpetIOHDF5 file */
  CHECK_ERROR (group = H5Gopen (infile, METADATA_GROUP));

  CHECK_ERROR (attr = H5Aopen_name (group, "nioprocs"));
  CHECK_ERROR (H5Aread (attr, H5T_NATIVE_INT, &nprocs));
  CHECK_ERROR (H5Aclose (attr));

  CHECK_ERROR (H5Gclose (group));
#endif

  /* create the metadata group in the IOHDF5 file */
  CHECK_ERROR (group = H5Gcreate (outfile, GLOBAL_ATTRIBUTES_GROUP, 0));
  CHECK_ERROR (attrDataspace = H5Screate (H5S_SCALAR));

  CHECK_ERROR (attr = H5Acreate (group, "nprocs", H5T_NATIVE_INT,
                                 attrDataspace, H5P_DEFAULT));
  CHECK_ERROR (H5Awrite (attr, H5T_NATIVE_INT, &nprocs));
  CHECK_ERROR (H5Aclose (attr));

  CHECK_ERROR (attr = H5Acreate (group, "ioproc_every", H5T_NATIVE_INT,
                                 attrDataspace, H5P_DEFAULT));
  CHECK_ERROR (H5Awrite (attr, H5T_NATIVE_INT, &nprocs));
  CHECK_ERROR (H5Aclose (attr));

  CHECK_ERROR (attr = H5Acreate (group, "unchunked", H5T_NATIVE_INT,
                                 attrDataspace, H5P_DEFAULT));
  unchunked = nprocs == 1;
  CHECK_ERROR (H5Awrite (attr, H5T_NATIVE_INT, &unchunked));
  CHECK_ERROR (H5Aclose (attr));

  CHECK_ERROR (H5Tset_size (hdf5String, strlen (cactus_version)));
  CHECK_ERROR (attr = H5Acreate (group, "Cactus version", hdf5String,
                                 attrDataspace, H5P_DEFAULT));
  CHECK_ERROR (H5Awrite (attr, hdf5String, cactus_version));
  CHECK_ERROR (H5Aclose (attr));

  CHECK_ERROR (H5Gclose (group));

  /* copy datasets as we traverse through the input file */
  pathname = "";
  CHECK_ERROR (H5Giterate (infile, "/", NULL, CopyObject, &outfile));

  CHECK_ERROR (H5Sclose (attrDataspace));
  CHECK_ERROR (H5Tclose (hdf5String));
  CHECK_ERROR (H5Fclose (outfile));
  CHECK_ERROR (H5Fclose (infile));

  return (0);
}


/*****************************************************************************/
/*                           local routines                                  */
/*****************************************************************************/
 /*@@
   @routine    CopyObject
   @date       Thu 10 Jan 2002
   @author     Thomas Radke
   @desc
               Iterator recursively called by H5Giterate() for every object
               in the input file
   @enddesc

   @calls      CopyAttribute

   @var        from
   @vdesc      identifier for the group the current object belongs to
   @vtype      hid_t
   @vio        in
   @endvar
   @var        objectname
   @vdesc      name of the current object
   @vtype      const char *
   @vio        in
   @endvar
   @var        _to
   @vdesc      user-supplied argument indicating the output object identifier
   @vtype      hid_t
   @vio        in
   @endvar

   @returntype int
   @returndesc
               0 - continue the iteration for following group objects
               1 - short-curcuit, no further iteration of this group
   @endreturndesc
@@*/
static herr_t CopyObject (hid_t from,
                          const char *objectname,
                          void *_to)
{  
  hid_t to, attr, datatype, dataspace, dataset;
  H5G_stat_t objectinfo;
  char *current_pathname, *name;
  size_t objectsize;
  int i, iteration, timelevel;
  char *data;
   
 
  /* skip the CarpetIOHDF5 metadata group */
  if (! strcmp (objectname, METADATA_GROUP))
  {
    return (0);
  }

  /* check the type of the current object */
  CHECK_ERROR (H5Gget_objinfo (from, objectname, 0, &objectinfo));
  if (objectinfo.type != H5G_DATASET)
  {
    fprintf (stderr, "WARNING: Found object '%s/%s' which is not a dataset ! "
                     "Object will not be copied.\n", pathname, objectname);
    return (0);
  }

  /* build the full pathname for the current to object to process */
  current_pathname = pathname;
  pathname = (char *) malloc (strlen (current_pathname) +
                              strlen (objectname) + 2);
  sprintf (pathname, "%s/%s", current_pathname, objectname);
  printf ("  copying dataset '%s'\n", pathname);

  /* get the output object identifier */
  to = *(hid_t *) _to;

  CHECK_ERROR (from = H5Dopen (from, objectname));
  CHECK_ERROR (attr = H5Aopen_name (from, "name"));
  CHECK_ERROR (datatype = H5Aget_type (attr));
  CHECK_ERROR (objectsize = H5Tget_size (datatype));
  data = malloc (objectsize + 1);
  CHECK_ERROR (H5Aread (attr, datatype, data));
  data[objectsize] = 0;
  CHECK_ERROR (H5Tclose (datatype));
  CHECK_ERROR (H5Aclose (attr));

  CHECK_ERROR (attr = H5Aopen_name (from, "timestep"));
  CHECK_ERROR (H5Aread (attr, H5T_NATIVE_INT, &iteration));
  CHECK_ERROR (H5Aclose (attr));

  CHECK_ERROR (attr = H5Aopen_name (from, "group_timelevel"));
  CHECK_ERROR (H5Aread (attr, H5T_NATIVE_INT, &timelevel));
  CHECK_ERROR (H5Aclose (attr));

  /* build the new name for the PUGH dataset */
  name = malloc (objectsize + 100);
  sprintf (name, "%s timelevel %d at iteration %d", data, timelevel, iteration);
  free (data);

  CHECK_ERROR (dataspace = H5Dget_space (from));
  CHECK_ERROR (datatype = H5Dget_type (from));
  objectsize = H5Tget_size (datatype); 
  if (H5Sis_simple (dataspace)) {
    objectsize *= H5Sget_simple_extent_npoints (dataspace);
  }

  data = malloc (objectsize);
  CHECK_ERROR (H5Dread (from, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data));

  if (! unchunked) {
    CHECK_ERROR (to = H5Gcreate (to, name, 0));
  }
  dataset = -1;
  for (i = 0; i < nprocs; i++) {
    char chunkname[20];
    sprintf (chunkname, "chunk%d", i);
    CHECK_ERROR (dataset = H5Dcreate (to, unchunked ? name : chunkname,
                                      datatype, dataspace, H5P_DEFAULT));
    CHECK_ERROR (H5Dwrite (dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           data));
    if (i + 1 < nprocs) {
      CHECK_ERROR (H5Dclose (dataset));
    }
  }
  free (name);
  free (data);
  CHECK_ERROR (H5Sclose (dataspace));
  CHECK_ERROR (H5Tclose (datatype));

  to = unchunked ? dataset : to;

  /* copy Carpet attributes */
  CHECK_ERROR (H5Aiterate (from, NULL, CopyAttribute, &to));
  CHECK_ERROR (H5Dclose (from));

  /* add attributes which aren't stored by Carpet but needed by PUGH */
  CHECK_ERROR (dataspace = H5Screate (H5S_SCALAR));

  const int grouptype = 403; /* CCTK_ARRAY */
  CHECK_ERROR (attr = H5Acreate (to, "grouptype", H5T_NATIVE_INT,
                                 dataspace, H5P_DEFAULT));
  CHECK_ERROR (H5Awrite (attr, H5T_NATIVE_INT, &grouptype));
  CHECK_ERROR (H5Aclose (attr));

  CHECK_ERROR (attr = H5Acreate (to, "ntimelevels", H5T_NATIVE_INT,
                                 dataspace, H5P_DEFAULT));
  CHECK_ERROR (H5Awrite (attr, H5T_NATIVE_INT, &ntimelevels));
  CHECK_ERROR (H5Aclose (attr));

  CHECK_ERROR (datatype = H5Tcopy (H5T_C_S1));
  CHECK_ERROR (H5Tset_size (datatype, strlen (groupname)));
  CHECK_ERROR (attr = H5Acreate (to, "groupname", datatype,
                                 dataspace, H5P_DEFAULT));
  CHECK_ERROR (H5Awrite (attr, datatype, groupname));
  CHECK_ERROR (H5Aclose (attr));
  CHECK_ERROR (H5Tclose (datatype));

  CHECK_ERROR (H5Sclose (dataspace));

  if (unchunked) {
    CHECK_ERROR (H5Dclose (to));
  } else {
    CHECK_ERROR (H5Gclose (to));
  }

  /* reset the pathname */
  free (pathname);
  pathname = current_pathname;

  return (0);
}


 /*@@
   @routine    CopyAttribute
   @date       Thu 10 Jan 2002
   @author     Thomas Radke
   @desc
               Iterator recursively called by H5Aiterate() for every attribute
               of an object (dataset or group)
   @enddesc

   @var        from
   @vdesc      identifier for the group or dataset to read the attribute from
   @vtype      hid_t
   @vio        in
   @endvar
   @var        attrname
   @vdesc      name of the current attribute
   @vtype      const char *
   @vio        in
   @endvar
   @var        _to
   @vdesc      user-supplied argument indicating the group or dataset
               to copy the attribute to
   @vtype      hid_t
   @vio        in
   @endvar

   @returntype int
   @returndesc
               0 - continue the iteration for following attributes
   @endreturndesc
@@*/
static herr_t CopyAttribute (hid_t from,
                             const char *attrname,
                             void *_to)
{ 
  hid_t attr, datatype, dataspace, to;
  size_t attrsize;
  void *value;


  /* get the target group/dataset identifier */
  to = *(hid_t *) _to;
  
  /* open the attribute given by its name, get type, dataspace, and value
     and just copy it */
  CHECK_ERROR (attr = H5Aopen_name (from, attrname));
  CHECK_ERROR (datatype = H5Aget_type (attr));
  CHECK_ERROR (dataspace = H5Aget_space (attr));
  attrsize = H5Tget_size (datatype);
  if (H5Sis_simple (dataspace) > 0)
  { 
    attrsize *= H5Sget_simple_extent_npoints (dataspace);
  }
  if (attrsize > 0)
  { 
    value = malloc (attrsize);
    CHECK_ERROR (H5Aread (attr, datatype, value));
    CHECK_ERROR (H5Aclose (attr));
    CHECK_ERROR (attr = H5Acreate (to, attrname, datatype, dataspace,
                                   H5P_DEFAULT));
    CHECK_ERROR (H5Awrite (attr, datatype, value));
    free (value);
  }
  CHECK_ERROR (H5Aclose (attr));
  CHECK_ERROR (H5Sclose (dataspace));
  CHECK_ERROR (H5Tclose (datatype));
  
  return (0);
}
