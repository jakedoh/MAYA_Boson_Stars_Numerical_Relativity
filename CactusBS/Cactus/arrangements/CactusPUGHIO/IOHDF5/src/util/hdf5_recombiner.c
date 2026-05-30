 /*@@
   @file      hdf5_recombiner.c
   @date      Sat 24 Feb 2001
   @author    Thomas Radke
   @desc
              This utility program recombines chunked Cactus output file(s)
              in HDF5 file format into a single unchunked HDF5 file.
   @enddesc
   @version   $Id$
 @@*/

#include "cctk.h"       /* HAVE_UNISTD_H */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>     /* sysconf(3) */
#endif

#include "CactusPUGHIO/IOHDF5Util/src/ioHDF5UtilGH.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusPUGHIO_IOHDF5_util_hdf5_recombiner_c)


/*****************************************************************************/
/*                           macro definitions                               */
/*****************************************************************************/
/* number of file descriptors that should be reserved for system usage
   (these are at least stdin, stdout, stderr; some MPI implementations
    may need additional descriptors for internal use) */
#define RESERVED_FILE_DESCRIPTORS  5

/* the size (in bytes) of the data sieving buffer for input files
   (should be larger than HDF5's default of 64 kB) */
#define SIEVE_BUFFERSIZE  (8 * 1024 * 1024)

/* macro to do an HDF5 call, check its return code, and print a warning
   in case of an error */
#define CHECK_ERROR(hdf5_call)                                                \
          {                                                                   \
            hid_t _error_code = hdf5_call;                                    \
                                                                              \
                                                                              \
            if (_error_code < 0)                                              \
            {                                                                 \
              fprintf (stderr, "WARNING: line %d: HDF5 call '%s' returned "   \
                               "error code %d\n",                             \
                                __LINE__, #hdf5_call, (int)_error_code);      \
              nerrors++;                                                      \
            }                                                                 \
          }


/*****************************************************************************/
/*                           global variables                                */
/*****************************************************************************/
/* NOTE: although it isn't good programming practice
         we make these variables global for convenience
         since they are accessed from recursively or
         indirectly called routines which only get passed
         a single user-supplied argument */
static char *pathname = NULL;      /* pathname of the current object */
static hid_t *infiles = NULL;      /* list of input file handles */
static hid_t outfile = -1 ;        /* output file handle */
static hid_t sieve_plist = -1 ;    /* input file data access property list */
static int max_filedescriptors = 0;/* maximum number of open files */
static char **infilenames = NULL;  /* list of input filenames */
static int nprocs = 0;             /* total number of processors */
static int ioproc_every = 0;       /* I/O was done on every N'th processor */
static int ninfiles = 0;           /* number of input files */
static unsigned int nerrors = 0;   /* global error counter */
static int single_precision = 0;   /* recombine fp data in single precision */

/*****************************************************************************/
/*                           local function prototypes                       */
/*****************************************************************************/
static herr_t CopyObject (hid_t copy_from,
                          const char *objectname,
                          void *arg);
static herr_t CopyAttribute (hid_t src,
                             const char *attr_name,
                             void *arg);
static int RecombineGroupData (const char *groupname);


 /*@@
   @routine    main
   @date       Sat 24 Feb 2001
   @author     Thomas Radke
   @desc
               Main routine of the HDF5 recombiner
   @enddesc

   @calls      CopyObject

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
int main (int argc, const char *argv[])
{
  int unchunked;
  char *tmp, *template;
  hid_t infile, group, attr, dataspace;


  /* say hello to the user */
  printf ("\n   -----------------------------"
          "\n   Cactus 4 HDF5 File Recombiner"
          "\n   -----------------------------\n\n");

  if (argc > 1 && strcmp (argv[1], "-single_precision") == 0)
  {
    single_precision = 1;
    if (argc-- == 4)
    {
      argv[1] = argv[2];
      argv[2] = argv[3];
    }
  }

  /* give some help if called with incorrect number of parameters */
  if (argc != 3)
  {
    fprintf (stderr, "Usage: %s [-single_precision] <chunked_infile0> <unchunked_outfile>\n"
                     "   eg, %s alp.file_0.h5 alp.h5\n\n", argv[0], argv[0]);
    return (0);
  }

  /* set the data sieving buffer to a much higher value than its default (64kB)
     in order to minimize the number of low-level reads when reading individual
     chunks as hyperslabs into the slice buffer in memory */
  CHECK_ERROR (sieve_plist = H5Pcreate (H5P_FILE_ACCESS));
  CHECK_ERROR (H5Pset_sieve_buf_size (sieve_plist, SIEVE_BUFFERSIZE));

  /* open (first) input file */
  H5E_BEGIN_TRY
  {
    infile = H5Fopen (argv[1], H5F_ACC_RDONLY, sieve_plist);
  } H5E_END_TRY
  if (infile < 0)
  {
    fprintf (stderr, "ERROR: Cannot open HDF5 input file '%s' !\n\n", argv[1]);
    return (-1);
  }
  /* read the file set info from the GLOBAL_ATTRIBUTES_GROUP group */
  CHECK_ERROR (group = H5Gopen (infile, GLOBAL_ATTRIBUTES_GROUP));
  if (group < 0)
  {
    fprintf (stderr, "ERROR: Cannot open attribute group '"
                     GLOBAL_ATTRIBUTES_GROUP "' in input file '%s' ! "
                     "Is this really a Cactus 4 HDF5 data file ?\n\n",
                     argv[1]);
    return (-1);
  }
  tmp = NULL;
  if (tmp == NULL)
  {
    attr = H5Aopen_name (group, "nprocs");
    if (attr < 0)
    {
      tmp = "Cannot find attribute 'nprocs'";
    }
    else
    {
      if (H5Aread (attr, H5T_NATIVE_INT, &nprocs) < 0)
      {
        tmp = "Cannot read attribute 'nprocs'";
      }
      CHECK_ERROR (H5Aclose (attr));
    }
  }
  if (tmp == NULL)
  {
    attr = H5Aopen_name (group, "ioproc_every");
    if (attr < 0)
    {
      tmp = "Cannot find attribute 'ioproc_every'";
    }
    else
    {
      if (H5Aread (attr, H5T_NATIVE_INT, &ioproc_every) < 0)
      {
        tmp = "Cannot read attribute 'ioproc_every'";
      }
      CHECK_ERROR (H5Aclose (attr));
    }
  }
  if (tmp == NULL)
  {
    attr = H5Aopen_name (group, "unchunked");
    if (attr < 0)
    {
      tmp = "Cannot find attribute 'unchunked'";
    }
    else
    {
      if (H5Aread (attr, H5T_NATIVE_INT, &unchunked) < 0)
      {
        tmp = "Cannot read attribute 'unchunked'";
      }
      CHECK_ERROR (H5Aclose (attr));
    }
  }
  CHECK_ERROR (H5Gclose (group));
  if (tmp)
  {
    fprintf (stderr, "ERROR: %s in attribute group '"
                     GLOBAL_ATTRIBUTES_GROUP "' of input file '%s' ! "
                     "Is this really a Cactus 4 HDF5 data file ?\n\n",
                     tmp, argv[1]);
    return (-1);
  }

  /* check if the input file already contains unchunked data */
  if (unchunked)
  {
    fprintf (stderr, "WARNING: Nothing to do ! Input file '%s' already "
                     "contains unchunked data.\n\n", argv[1]);
    return (0);
  }

  /* create output file */
  outfile = H5Fcreate (argv[2], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (outfile < 0)
  {
    fprintf (stderr, "ERROR: Cannot create output file '%s' !\n\n", argv[2]);
    return (-1);
  }

  /* get the number of chunked input files to process */
  ninfiles = nprocs / ioproc_every;
  if (nprocs % ioproc_every)
  {
    ninfiles++;
  }

  /* determine maximum number of files opened simultaneously */
  max_filedescriptors = sysconf (_SC_OPEN_MAX);
  if (max_filedescriptors < 0)
  {
    fprintf (stderr, "WARNING: Cannot determine filehandle limit ! "
                     "Assuming no limit...\n");
  }
  /* subtract by number of reserved file descriptors
     and the one for the output file */
  max_filedescriptors -= RESERVED_FILE_DESCRIPTORS + 1;
  if (max_filedescriptors < 0)
  {
    max_filedescriptors = ninfiles;
  }

  printf ("Recombining HDF5 data from %d processors, "
          "output was written to %d file(s)\n\n", nprocs, ninfiles);

  /* allocate arrays for input file descriptors and filenames */
  infiles = malloc (ninfiles * sizeof (hid_t));
  infiles[0] = infile;
  infilenames = malloc (ninfiles * sizeof (char *));

  /* now get all chunked input filenames and check that they can be opened */
  if (ninfiles == 1)
  {
    /* not much to be done here */
    infilenames[0] = strdup (argv[1]);
  }
  else
  {
    /* close the first file (it might not be the one written by processor 0) */
    CHECK_ERROR (H5Fclose (infiles[0]));

    /* get the basename of input file(s) */
    tmp = strstr (argv[1], ".file_");
    if (tmp == NULL)
    {
      fprintf (stderr, "ERROR: Cannot parse input file name '%s' ! "
                       "Is this really a chunked Cactus HDF5 input file ?\n\n",
                       argv[1]);
      return (-1);
    }

    /* build the input filename template */
    template = malloc (strlen (argv[1]) + 2);
    strcpy (template, argv[1]);
    template[tmp - argv[1] + 6] = 0;
    strcat (template, "%d.h5");

    /* now loop through all the files */
    for (infile = 0; infile < ninfiles; infile++)
    {
      /* build the input filename */
      infilenames[infile] = malloc (strlen (template) + 10);
      sprintf (infilenames[infile], template, infile);
      infiles[infile] = H5Fopen (infilenames[infile], H5F_ACC_RDONLY,
                                 sieve_plist);
      if (infiles[infile] < 0)
      {
        fprintf (stderr, "ERROR: Cannot open chunked HDF5 input file '%s' !\n",
                 infilenames[infile]);
        return (-1);
      }

      /* close file if filehandle limit would be exceeded */
      if (infile > max_filedescriptors)
      {
        CHECK_ERROR (H5Fclose (infiles[infile]));
      }
    }

    free (template);
  }

  /* do the recombination by iterating over all objects
     in the (first) input file written by processor 0 */
  pathname = "";
  CHECK_ERROR (H5Giterate (infiles[0], "/", NULL, CopyObject, &outfile));

  /* now reset the file info attributes in the GLOBAL_ATTRIBUTES_GROUP group
     to indicate unchunked file data */
  CHECK_ERROR (group = H5Gopen (outfile, GLOBAL_ATTRIBUTES_GROUP));
  CHECK_ERROR (H5Adelete (group, "unchunked"));
  CHECK_ERROR (H5Adelete (group, "nprocs"));
  CHECK_ERROR (H5Adelete (group, "ioproc_every"));

  unchunked = 1;
  nprocs = ioproc_every = 1;
  CHECK_ERROR (dataspace = H5Screate (H5S_SCALAR));
  CHECK_ERROR (attr = H5Acreate (group, "unchunked", H5T_NATIVE_INT,
                                 dataspace, H5P_DEFAULT));
  CHECK_ERROR (H5Awrite (attr, H5T_NATIVE_INT, &unchunked));
  CHECK_ERROR (H5Aclose (attr));
  CHECK_ERROR (attr = H5Acreate (group, "nprocs", H5T_NATIVE_INT,
                                 dataspace, H5P_DEFAULT));
  CHECK_ERROR (H5Awrite (attr, H5T_NATIVE_INT, &nprocs));
  CHECK_ERROR (H5Aclose (attr));
  CHECK_ERROR (attr = H5Acreate (group, "ioproc_every", H5T_NATIVE_INT,
                                 dataspace, H5P_DEFAULT));
  CHECK_ERROR (H5Awrite (attr, H5T_NATIVE_INT, &ioproc_every));
  CHECK_ERROR (H5Aclose (attr));
  CHECK_ERROR (H5Sclose (dataspace));
  CHECK_ERROR (H5Gclose (group));

  /* finally, close all open files */
  for (infile = 0; infile < ninfiles; infile++)
  {
    if (infile <= max_filedescriptors)
    {
      CHECK_ERROR (H5Fclose (infiles[infile]));
    }
    free (infilenames[infile]);
  }
  free (infiles);
  CHECK_ERROR (H5Fclose (outfile));

  /* report status */
  if (nerrors == 0)
  {
    printf ("\n\n   *** All Cactus data successfully recombined. ***\n\n");
  }
  else
  {
    fprintf (stderr, "\n\n   *** WARNING: %d errors occured during "
                     "recombination. ***\n\n", nerrors);
  }

  return (0);
}


/*****************************************************************************/
/*                           local routines                                  */
/*****************************************************************************/
 /*@@
   @routine    CopyObject
   @date       Sat 24 Feb 2001
   @author     Thomas Radke
   @desc
               Iterator recursively called by H5Giterate() for every object
               in the input file
               It copies the current object or invokes the recombiner on it.
   @enddesc

   @calls      CopyAttribute
               RecombineGroupData

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
  hid_t to, datatype, dataspace;
  H5G_stat_t objectinfo;
  char *current_pathname;
  int retval;
  size_t objectsize;
  void *data;


  /* build the full pathname for the current to object to process */
  current_pathname = pathname;
  pathname = malloc (strlen (current_pathname) + strlen (objectname) + 2);
  sprintf (pathname, "%s/%s", current_pathname, objectname);

  /* get the output object identifier */
  retval = 0;
  to = *(hid_t *) _to;

  /* check the type of the current object */
  CHECK_ERROR (H5Gget_objinfo (from, objectname, 0, &objectinfo));
  if (objectinfo.type == H5G_GROUP)
  {
    /* try to recombine data within this group */
    if (RecombineGroupData (pathname) <= 0)
    {
      /* this group didn't contain chunked data
         so just copy it as is */
      printf ("  copying group '%s'\n", pathname);

      CHECK_ERROR (from = H5Gopen (from, objectname));
      CHECK_ERROR (to = H5Gcreate (to, objectname, 0));
      /* iterate over all objects in the (first) input file */
      CHECK_ERROR (H5Giterate (from, ".", NULL, CopyObject, &to));
      CHECK_ERROR (H5Aiterate (from, NULL, CopyAttribute, &to));
      CHECK_ERROR (H5Gclose (to));
      CHECK_ERROR (H5Gclose (from));
    }
  }
  else if (objectinfo.type == H5G_DATASET)
  {
    /* current object is an unchunked dataset - copy it as is */
    printf ("  copying dataset '%s'\n", pathname);

    CHECK_ERROR (from = H5Dopen (from, objectname));
    CHECK_ERROR (datatype = H5Dget_type (from));

    /* recombine floating-point data in single precision if requested */
    if (single_precision && H5Tget_class (datatype) == H5T_FLOAT)
    {
      CHECK_ERROR (H5Tclose (datatype));
      CHECK_ERROR (datatype = H5Tcopy (H5T_NATIVE_FLOAT));
    }

    CHECK_ERROR (dataspace = H5Dget_space (from));
    CHECK_ERROR (to = H5Dcreate (to, objectname, datatype, dataspace,
                                 H5P_DEFAULT));
    objectsize = H5Sget_select_npoints (dataspace) * H5Tget_size (datatype);
    if (objectsize > 0)
    {
      data = malloc (objectsize);
      CHECK_ERROR (H5Dread (from, datatype, H5S_ALL, H5S_ALL,H5P_DEFAULT,data));
      CHECK_ERROR (H5Dwrite (to, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,data));
      free (data);
    }
    CHECK_ERROR (H5Aiterate (from, NULL, CopyAttribute, &to));
    CHECK_ERROR (H5Dclose (to));
    CHECK_ERROR (H5Dclose (from));
    CHECK_ERROR (H5Sclose (dataspace));
    CHECK_ERROR (H5Tclose (datatype));
  }
  else
  {
    fprintf (stderr, "WARNING: Found object '%s' which is neither a "
                     "group nor a dataset ! Object will not be copied.\n",
                     pathname);
    nerrors++;
  }

  /* reset the pathname */
  free (pathname);
  pathname = current_pathname;

  return (retval);
}


 /*@@
   @routine    CopyAttribute
   @date       Sat 24 Feb 2001
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
   @var        name
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
                             const char *name,
                             void *_to)
{
  hid_t attr, datatype, dataspace, to;
  size_t attrsize;
  void *value;


  /* get the target group/dataset identifier */
  to = *(hid_t *) _to;

  /* open the attribute given by its name, get type, dataspace, and value
     and just copy it */
  CHECK_ERROR (attr = H5Aopen_name (from, name));
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
    CHECK_ERROR (attr = H5Acreate (to, name, datatype, dataspace, H5P_DEFAULT));
    CHECK_ERROR (H5Awrite (attr, datatype, value));
    free (value);
  }
  CHECK_ERROR (H5Aclose (attr));
  CHECK_ERROR (H5Sclose (dataspace));
  CHECK_ERROR (H5Tclose (datatype));

  return (0);
}


 /*@@
   @routine    RecombineGroupData
   @date       Sat 24 Feb 2001
   @author     Thomas Radke
   @desc
               Checks whether the passed group indicates a group with chunked
               datasets. If so it will recombine these from all chunked input
               files.
               In order to minimize the number of write I/O calls, the
               recombination is done in z-slices (the last changing slice
               in a 3D dataset). Each z-slice is recombined in memory from the
               individual chunks, and then written out by a single H5Dwrite().
   @enddesc

   @calls      CopyAttribute

   @var        groupname
   @vdesc      name of the current group to process
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               0 - group does not contain chunked data
               1 - chunked group data successfully recombined
   @endreturndesc
@@*/
static int RecombineGroupData (const char *groupname)
{
  int infile, chunk, nchunks;
  size_t size;
  hid_t group, attr, datatype, dataspace;
  hid_t chunked_datatype, chunked_dataspace, chunked_dataset;
  hid_t unchunked_datatype, unchunked_dataspace, unchunked_dataset;
  hid_t slice_dataspace;
  hsize_t tmp1, ndims, chunk_ndims, *dims, *chunk_dims, *slice_dims;
#if (H5_VERS_MAJOR == 1 && \
     (H5_VERS_MINOR < 6 || (H5_VERS_MINOR == 6 && H5_VERS_RELEASE < 4)))
  hssize_t tmp2, *chunk_origin, *slice_origin;
#else
  hsize_t tmp2, *chunk_origin, *slice_origin;
#endif
  void *data;
  char *chunkname;
  H5T_class_t class;


  /* open the group given by its name */
  CHECK_ERROR (group = H5Gopen (infiles[0], groupname));

  /* a group with chunked data must have a 'global_size' attribute */
  H5E_BEGIN_TRY
  {
    attr = H5Aopen_name (group, "global_size");
  } H5E_END_TRY;
  if (attr < 0)
  {
    CHECK_ERROR (H5Gclose (group));
    return (0);
  }

  /* read the global size of the unchunked dataset */
  datatype = H5Aget_type (attr);
  dataspace = H5Aget_space (attr);
  ndims = H5Sget_simple_extent_npoints (dataspace);
  if (H5Tget_class (datatype) != H5T_INTEGER)
  {
    fprintf (stderr, "WARNING: 'global_size' attribute of group '%s' is not "
                     "of type integer !\n", groupname);
    ndims = 0;
  }
  dims = NULL;
  if (ndims > 0)
  {
    dims = calloc (ndims, sizeof (hsize_t));
    CHECK_ERROR (H5Aread (attr, H5T_NATIVE_HSIZE, dims));
  }
  CHECK_ERROR (H5Sclose (dataspace));
  CHECK_ERROR (H5Tclose (datatype));
  CHECK_ERROR (H5Aclose (attr));

  /* return if the 'global_size' attribute couldn't be read */
  if (ndims <= 0)
  {
    CHECK_ERROR (H5Gclose (group));
    return (0);
  }

  /* allocate string buffer holding the possible dataset chunk names */
  chunkname = malloc (strlen (groupname) + sizeof ("chunk") + 20);

  /* the unchunked dataset's dimensions are taken from the 'global_size'
     attribute, with the least changing dimension being the first element */
  for (chunk_ndims = 0; chunk_ndims < ndims/2; chunk_ndims++)
  {
    tmp1 = dims[chunk_ndims];
    dims[chunk_ndims] = dims[ndims - chunk_ndims - 1];
    dims[ndims - chunk_ndims - 1] = tmp1;
  }
  CHECK_ERROR (unchunked_dataspace = H5Screate_simple (ndims, dims,
                                                       NULL));

  /* allocate buffers to read a 'chunk_origin' attribute and the chunk dims */
  chunk_origin = calloc (2 * ndims, sizeof (*chunk_origin));
  chunk_dims   = calloc (2 * ndims, sizeof (hsize_t));
  slice_origin = chunk_origin + ndims;
  slice_dims   = chunk_dims + ndims;

  unchunked_dataset = unchunked_datatype = slice_dataspace = -1;
  data = NULL;

  /* now read all the chunks from all input files and write them into the
     unchunked output dataset as a hyperslab */
  for (infile = 0; infile < ninfiles; infile++)
  {
    /* get the number of chunks in this file */
    nchunks = ioproc_every;
    if (infile == ninfiles - 1 && nprocs % ioproc_every)
    {
      nchunks = nprocs % ninfiles;
    }

    /* re-open the file if it was closed before */
    if (infile > max_filedescriptors)
    {
#ifdef DEBUG
      printf ("reopening input file '%s'\n", infilenames[infile]);
#endif
      CHECK_ERROR (infiles[infile] = H5Fopen (infilenames[infile],
                                              H5F_ACC_RDONLY, sieve_plist));
    }

    /* loop over all chunks of this input file */
    for (chunk = 0; chunk < nchunks; chunk++)
    {
      /* build the object name of this chunk */
      sprintf (chunkname, "%s/chunk%d", pathname, chunk);

      /* try to open the dataset chunk */
      H5E_BEGIN_TRY
      {
        chunked_dataset = H5Dopen (infiles[infile], chunkname);
      } H5E_END_TRY;
      if (chunked_dataset < 0)
      {
        continue;
      }

      /* read the 'chunk_origin' attribute of this chunk */
      CHECK_ERROR (attr = H5Aopen_name (chunked_dataset, "chunk_origin"));
      CHECK_ERROR (datatype = H5Aget_type (attr));
      CHECK_ERROR (dataspace = H5Aget_space (attr));
      class = H5Tget_class (datatype);
      chunk_ndims = H5Sget_simple_extent_npoints (dataspace);
      if (chunk_ndims == ndims)
      {
#if (H5_VERS_MAJOR == 1 && \
     (H5_VERS_MINOR < 6 || (H5_VERS_MINOR == 6 && H5_VERS_RELEASE < 4)))
        CHECK_ERROR (H5Aread (attr, H5T_NATIVE_HSSIZE, chunk_origin));
#else
        CHECK_ERROR (H5Aread (attr, H5T_NATIVE_HSIZE, chunk_origin));
#endif
      }
      CHECK_ERROR (H5Sclose (dataspace));
      CHECK_ERROR (H5Tclose (datatype));
      CHECK_ERROR (H5Aclose (attr));

      /* check consistency */
      if (class != H5T_INTEGER)
      {
        fprintf (stderr, "WARNING: 'chunk_origin' attribute of dataset '%s' "
                         "is not of type integer !\n", chunkname);
        CHECK_ERROR (H5Dclose (chunked_dataset));
        nerrors++;
        continue;
      }
      if (chunk_ndims != ndims)
      {
        fprintf (stderr, "WARNING: dimensions of 'chunk_origin' attribute of "
                         "dataset '%s' don't match with dataset dimensions !\n",
                         chunkname);
        CHECK_ERROR (H5Dclose (chunked_dataset));
        nerrors++;
        continue;
      }

      /* the first time through, create the unchunked dataset */
      if (unchunked_dataset < 0)
      {
        printf ("  recombining dataset '%s'\n", pathname);

        /* the unchunked dataset gets the same datatype as the first chunk */
        CHECK_ERROR (unchunked_datatype = H5Dget_type (chunked_dataset));

        /* recombine floating-point data in single precision if requested */
        if (single_precision && H5Tget_class (unchunked_datatype) == H5T_FLOAT)
        {
          CHECK_ERROR (H5Tclose (unchunked_datatype));
          CHECK_ERROR (unchunked_datatype = H5Tcopy (H5T_NATIVE_FLOAT));
        }

        CHECK_ERROR (unchunked_dataset = H5Dcreate (outfile, pathname,
                                                    unchunked_datatype,
                                                    unchunked_dataspace,
                                                    H5P_DEFAULT));
      }

      /* check the datatype of this chunk to be consistent with previous ones */
      CHECK_ERROR (chunked_datatype = H5Dget_type (chunked_dataset));
      class = H5Tget_class (chunked_datatype);
      CHECK_ERROR (H5Tclose (chunked_datatype));
      if (class != H5Tget_class (unchunked_datatype))
      {
        fprintf (stderr, "WARNING: datatype class of chunk '%s' differs from "
                         "first chunk's datatype ! Chunk will be omitted.\n",
                         chunkname);
        CHECK_ERROR (H5Dclose (chunked_dataset));
        nerrors++;
        continue;
      }

      /* get the chunk dims */
      CHECK_ERROR (chunked_dataspace = H5Dget_space (chunked_dataset));
      chunk_ndims = H5Sget_simple_extent_ndims (chunked_dataspace);
      if (chunk_ndims != ndims)
      {
        fprintf (stderr, "WARNING: dimensions of 'chunk_origin' attribute of "
                         "dataset '%s' don't match with dataset dimensions !\n",
                         chunkname);
        CHECK_ERROR (H5Sclose (chunked_dataspace));
        CHECK_ERROR (H5Dclose (chunked_dataset));
        nerrors++;
        continue;
      }
      CHECK_ERROR (H5Sget_simple_extent_dims (chunked_dataspace, chunk_dims,
                                              NULL));

      /* check the chunk's dataspace to be a simple one */
      if (H5Sis_simple (chunked_dataspace) <= 0)
      {
        fprintf (stderr, "WARNING: dataset '%s' is not a simple "
                 "multidimensional dataset (dataset will be ignored) !\n",
                 chunkname);
        CHECK_ERROR (H5Sclose (chunked_dataspace));
        CHECK_ERROR (H5Dclose (chunked_dataset));
        nerrors++;
        continue;
      }

      /* HDF5 needs the least changing dimension first */
      for (chunk_ndims = 0; chunk_ndims < ndims/2; chunk_ndims++)
      {
        tmp2 = chunk_origin[chunk_ndims];
        chunk_origin[chunk_ndims] = chunk_origin[ndims - chunk_ndims - 1];
        chunk_origin[ndims - chunk_ndims - 1] = tmp2;
      }

      /* give some info output */
      printf ("    - file %d chunk %d\n", infile, chunk + ioproc_every*infile);
      printf ("      chunk dimensions: [%d", (int) chunk_dims[ndims - 1]);
      for (chunk_ndims = 1; chunk_ndims < ndims; chunk_ndims++)
      {
        printf (", %d", (int) chunk_dims[ndims - chunk_ndims - 1]);
      }
      printf ("]     chunk origin: [%d", (int) chunk_origin[ndims - 1]);
      for (chunk_ndims = 1; chunk_ndims < ndims; chunk_ndims++)
      {
        printf (", %d", (int) chunk_origin[ndims - chunk_ndims - 1]);
      }
      printf ("]\n");

      if (slice_dataspace < 0 &&
          chunk_origin[ndims - 1] == 0 &&
          (ndims == 1 || chunk_origin[ndims - 2] == 0))
      {
        memset (slice_origin, 0, ndims * sizeof (*slice_origin));
        memcpy (slice_dims, chunk_dims, ndims * sizeof (hsize_t));
        slice_dims[ndims - 1] = dims[ndims - 1];
        if (ndims > 1)
        {
          slice_dims[ndims - 2] = dims[ndims - 2];
        }
        CHECK_ERROR (slice_dataspace = H5Screate_simple (ndims, slice_dims,
                                                         NULL));
        size = H5Tget_size (unchunked_datatype) *
               (size_t) H5Sget_simple_extent_npoints (slice_dataspace);
        data = malloc (size);
        if (! data)
        {
          fprintf (stderr, "WARNING: couldn't allocate %d bytes to recombine "
                   "slice of dataset '%s' !\n", (int) size, chunkname);
          CHECK_ERROR (H5Sclose (slice_dataspace));
          CHECK_ERROR (H5Sclose (chunked_dataspace));
          CHECK_ERROR (H5Dclose (chunked_dataset));
          nerrors++;
          continue;
        }
      }
      slice_origin[ndims - 1] = chunk_origin[ndims - 1];
      if (ndims > 1)
      {
        slice_origin[ndims - 2] = chunk_origin[ndims - 2];
      }

      if (slice_dataspace >= 0)
      {
        CHECK_ERROR (H5Sselect_hyperslab (slice_dataspace, H5S_SELECT_SET,
                                          slice_origin, NULL, chunk_dims,NULL));
        /* set the buffer size to read the full hyperslab */
        CHECK_ERROR (H5Dread (chunked_dataset, unchunked_datatype,
                              slice_dataspace, H5S_ALL, H5P_DEFAULT, data));
      }
      CHECK_ERROR (H5Dclose (chunked_dataset));
      CHECK_ERROR (H5Sclose (chunked_dataspace));

      /* if this was the last chunk for the current z-slice
         then write it back now */
      if (slice_dataspace >= 0 &&
          chunk_dims[ndims - 1] + chunk_origin[ndims - 1] == dims[ndims - 1] &&
          (ndims == 1 ||
           chunk_dims[ndims - 2] + chunk_origin[ndims - 2] == dims[ndims - 2]))
      {
        /* reset the x/y-offsets to zero */
        chunk_origin[ndims - 1] = 0;
        if (ndims > 1)
        {
          chunk_origin[ndims - 2] = 0;
        }
        CHECK_ERROR (H5Sselect_hyperslab (unchunked_dataspace, H5S_SELECT_SET,
                                          chunk_origin, NULL, slice_dims,NULL));
        CHECK_ERROR (H5Sselect_all (slice_dataspace));
        CHECK_ERROR (H5Dwrite (unchunked_dataset, unchunked_datatype,
                               slice_dataspace, unchunked_dataspace,
                               H5P_DEFAULT, data));
        free (data);
        CHECK_ERROR (H5Sclose (slice_dataspace));
        slice_dataspace = -1;
      }

    } /* end of loop over all chunks in this input file */

    /* close input file if filehandle limit would be exceeded */
    if (infile > max_filedescriptors)
    {
#ifdef DEBUG
      printf ("temporarily closing input file '%s'\n", infilenames[infile]);
#endif
      CHECK_ERROR (H5Fclose (infiles[infile]));
    }
  } /* end of loop over all input files */

  /* close objects and free allocated resources */
  if (unchunked_dataset >= 0)
  {
    /* copy all group attributes to the new dataset */
    CHECK_ERROR (H5Aiterate (group, NULL, CopyAttribute, &unchunked_dataset));
    CHECK_ERROR (H5Dclose (unchunked_dataset));
    CHECK_ERROR (H5Tclose (unchunked_datatype));
  }
  CHECK_ERROR (H5Sclose (unchunked_dataspace));
  CHECK_ERROR (H5Gclose (group));
  free (chunk_dims);
  free (chunk_origin);
  free (chunkname);
  free (dims);

  /* indicate no further processing of this group in H5Giterate() */
  return (1);
}
