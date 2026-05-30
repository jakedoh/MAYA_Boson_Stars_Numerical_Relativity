 /*@@
   @file      hdf5_convert_from_sdf.c
   @date      Tue 24 August 2004
   @author    Thomas Radke
   @desc
              Utility program to convert SDF datafiles into IOHDF5
              datafiles.
   @enddesc
   @version   $Id$
 @@*/

#define  MAXDIM       3
#define  MAXNAMESIZE  100

#include <stdio.h>
#include <stdlib.h>


/* Cactus includes (defines CCTK_FILEVERSION) */
#include "cctk.h"

/* SDF includes */
#include <bbhutil.h>
#include <sdf_priv.h>

/* HDF5 include */
#define H5_USE_16_API 1
#include <hdf5.h>

#define GLOBAL_ATTRIBUTES_GROUP "Global Attributes"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusPUGHIO_IOHDF5_util_hdf5_convert_from_sdf_c)


/*****************************************************************************/
/*                           macro definitions                               */
/*****************************************************************************/
/* uncomment the following to get some debugging output */
/* #define IOHDF5_DEBUG  1 */

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


 /*@@
   @routine    main
   @date       Tue 24 August 2004
   @author     Thomas Radke
   @desc
               Main routine of the SDF-to-HDF5 converter
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
  int i, j;
  int version, rank, dsize, csize;
  const int grouptype = CCTK_GF;
  int *dims;
  double timestep;
  char *varname, *coordname, *tag;
  const char *groupname;
  gft_sdf_file_data *infile;
  double *coords, *data, *bbox, *origin, *delta;
  hid_t outfile, dataset, dataspace, attr, attrDataspace;
  hid_t group, hdf5String;
  hsize_t hdf5Dims[MAXDIM];
  int iteration, timelevel, ntimelevels;
  char hdf5DatasetName[2 * MAXNAMESIZE];
  const char *cactus_version = "Cactus Beta14";


  if (argc != 7)
  {
    fprintf (stderr, "Usage: %s <groupname> <iteration> <timelevel> "
             "<ntimelevels> <inputfile> <outputfile>\n", argv[0]);
    fprintf (stderr, "   eg. %s ADMBase::metric 0 0 3 foobar.sdf gxx.h5\n",
             argv[0]);
    return (0);
  }

  groupname = argv[1];
  iteration = atoi (argv[2]);
  timelevel = atoi (argv[3]);
  ntimelevels = atoi (argv[4]);

  infile = gft_open_sdf_file (argv[5]);
  if (! infile)
  {
    printf ("Could not open SDF input file '%s'\n", argv[5]);
    return (-1);
  }

  outfile = H5Fcreate (argv[6], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (outfile < 0)
  {
    fprintf (stderr, "Could not create HDF5 output file '%s'\n", argv[6]);
    return (-1);
  }

  printf ("\n  ----------------------------\n"
            "  Cactus SDF-to-HDF5 Converter\n"
            "  ----------------------------\n");

  CHECK_ERROR (hdf5String = H5Tcopy (H5T_C_S1));

  /* add a dummy GLOBAL_ATTRIBUTES_GROUP so that the HDF5 file is recognized as
     unchunked Cactus data */
  CHECK_ERROR (group = H5Gcreate (outfile, GLOBAL_ATTRIBUTES_GROUP, 0));
  CHECK_ERROR (attrDataspace = H5Screate (H5S_SCALAR));

  CHECK_ERROR (attr = H5Acreate (group, "nprocs", H5T_NATIVE_INT,
                                      attrDataspace, H5P_DEFAULT));
  i = 1;
  CHECK_ERROR (H5Awrite (attr, H5T_NATIVE_INT, &i));
  CHECK_ERROR (H5Aclose (attr));

  CHECK_ERROR (attr = H5Acreate (group, "ioproc_every", H5T_NATIVE_INT,
                                      attrDataspace, H5P_DEFAULT));
  i = 1;
  CHECK_ERROR (H5Awrite (attr, H5T_NATIVE_INT, &i));
  CHECK_ERROR (H5Aclose (attr));

  CHECK_ERROR (attr = H5Acreate (group, "unchunked", H5T_NATIVE_INT,
                                      attrDataspace, H5P_DEFAULT));
  i = 1;
  CHECK_ERROR (H5Awrite (attr, H5T_NATIVE_INT, &i));
  CHECK_ERROR (H5Aclose (attr));

  CHECK_ERROR (H5Tset_size (hdf5String, strlen (cactus_version)));
  CHECK_ERROR (attr = H5Acreate (group, "Cactus version", hdf5String,
                                 attrDataspace, H5P_DEFAULT));
  CHECK_ERROR (H5Awrite (attr, hdf5String, cactus_version));
  CHECK_ERROR (H5Aclose (attr));

  CHECK_ERROR (H5Gclose (group));

  while (low_read_sdf_stream (1, infile->fp, &timestep, &version, &rank, &dsize,
                              &csize, &varname, &coordname, &tag, &dims, &bbox,
                              &coords, &data))
  {
    printf ("Processing dataset '%s' (timestep %f)\n", varname, timestep);

    /* convert from int to hsize_t */
    for (j = 0; j < rank; j++)
    {
      hdf5Dims[rank-j-1] = dims[j];
    }
    CHECK_ERROR (dataspace = H5Screate_simple (rank, hdf5Dims, NULL));

    sprintf (hdf5DatasetName, "/%s timelevel %d at iteration %d",
             varname, timelevel, iteration);
    CHECK_ERROR (dataset = H5Dcreate (outfile, hdf5DatasetName,
                                      H5T_NATIVE_DOUBLE,dataspace,H5P_DEFAULT));
    CHECK_ERROR (H5Dwrite (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                           H5P_DEFAULT, data));
    CHECK_ERROR (H5Sclose (dataspace));

    /* attach necessary attributes */
    CHECK_ERROR (H5Tset_size (hdf5String, strlen (groupname)));
    CHECK_ERROR (attr = H5Acreate (dataset, "groupname", hdf5String,
                                   attrDataspace, H5P_DEFAULT));
    CHECK_ERROR (H5Awrite (attr, hdf5String, groupname));
    CHECK_ERROR (H5Aclose (attr));

    CHECK_ERROR (attr = H5Acreate (dataset, "grouptype", H5T_NATIVE_INT,
                                   attrDataspace, H5P_DEFAULT));
    CHECK_ERROR (H5Awrite (attr, H5T_NATIVE_INT, &grouptype));
    CHECK_ERROR (H5Aclose (attr));

    CHECK_ERROR (attr = H5Acreate (dataset, "ntimelevels", H5T_NATIVE_INT,
                                   attrDataspace, H5P_DEFAULT));
    CHECK_ERROR (H5Awrite (attr, H5T_NATIVE_INT, &ntimelevels));
    CHECK_ERROR (H5Aclose (attr));

    CHECK_ERROR (attr = H5Acreate (dataset, "time", H5T_NATIVE_DOUBLE,
                                   attrDataspace, H5P_DEFAULT));
    CHECK_ERROR (H5Awrite (attr, H5T_NATIVE_DOUBLE, &timestep));
    CHECK_ERROR (H5Aclose (attr));

    /* set origin/delta from bbox information */
    origin = malloc (2*rank * sizeof (double));
    delta  = origin + rank;
    for (j = 0; j < rank; j++)
    {
      origin[j] = bbox[2*j+0];
      delta[j]  = (bbox[2*j+1] - bbox[2*j+0]) / (dims[j] - 1);
    }
    hdf5Dims[0] = rank;
    CHECK_ERROR (dataspace = H5Screate_simple (1, hdf5Dims, NULL));
    CHECK_ERROR (attr = H5Acreate (dataset, "origin", H5T_NATIVE_DOUBLE,
                                   dataspace, H5P_DEFAULT));
    CHECK_ERROR (H5Awrite (attr, H5T_NATIVE_DOUBLE, origin));
    CHECK_ERROR (H5Aclose (attr));
    CHECK_ERROR (attr = H5Acreate (dataset, "delta", H5T_NATIVE_DOUBLE,
                                   dataspace, H5P_DEFAULT));
    CHECK_ERROR (H5Awrite (attr, H5T_NATIVE_DOUBLE, delta));
    CHECK_ERROR (H5Aclose (attr));
    CHECK_ERROR (H5Sclose (dataspace));
    CHECK_ERROR (H5Dclose (dataset));

    free (varname);
    free (coordname);
    free (tag);
    free (dims);
    free (bbox);
    free (coords);
    free (data);
    free (origin);

    /*** FIXME: increment iteration number to distinguish datasets ***/
    iteration++;
  }

  CHECK_ERROR (H5Sclose (attrDataspace));
  CHECK_ERROR (H5Tclose (hdf5String));
  CHECK_ERROR (H5Fclose (outfile));
  gsfd_close (infile);

  return (0);
}
