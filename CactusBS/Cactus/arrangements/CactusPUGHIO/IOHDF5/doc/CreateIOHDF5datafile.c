 /*@@
   @file      CreateIOHDF5datafile.c
   @date      Mon 12 Mar 2001
   @author    Thomas Radke
   @desc
              Example program to create an unchunked HDF5 datafile
              with a single dataset which can be read as input data
              into Cactus.
   @enddesc
   @version   $Id$
 @@*/


#include <stdio.h>
#include <stdlib.h>

#define H5_USE_16_API 1
#include <hdf5.h>


/* the name of our sample data file */
#define DATAFILENAME  "x.h5"

/* the number of dimensions of our sample data array and its size */
#define  NDIM         3
#define  NSIZE        20

/* maximum size of the dataset name */
#define  MAXNAMESIZE  100

/* the name of the attributes' group describing how the data file was created
   and what type of Cactus variable the sample dataset is (a grid function)
   These definitions were taken from Cactus header files - for ease of use
   we didn't include these here. */
#define GLOBAL_ATTRIBUTES_GROUP  "Global Attributes"
#define CCTK_GF                  2

/* the Cactus version that this datafile is (upwards) compatible with 
   (versions below beta10 use a different timelevel scheme) */
#define CACTUS_VERSION  "4.0.b10"

/* a simple macro to do an HDF5 call with return code checking
   in case of an error it will issue an error message and exit */
#define CHECK_ERROR(hdf5_call)                                                \
          {                                                                   \
            hid_t _error_code = hdf5_call;                                    \
                                                                              \
                                                                              \
            if (_error_code < 0)                                              \
            {                                                                 \
              fprintf (stderr, "ERROR: line %d: HDF5 call '%s' returned "     \
                               "error code %d\n",                             \
                                __LINE__, #hdf5_call, (int)_error_code);      \
              return (-1);                                                    \
            }                                                                 \
          }


 /*@@
   @routine    main
   @date       Mon 12 Mar 2001
   @author     Thomas Radke
   @desc
               Main routine creating a sample HDF5 datafile
   @enddesc

   @returntype int
   @returndesc
               0 for success, negative return values indicate an error
   @endreturndesc
@@*/
int main (void)
{
  void *data;
  int i, elements;
  int iteration, grouptype;
  int timelevel, ntimelevels;
  int nprocs, ioproc_every, unchunked;
  hsize_t dims[NDIM];
  char datasetname[MAXNAMESIZE], *varname, *groupname;
  hid_t datafile, group, attr, dataset, dataspace, stringtype;


  /* create a datafile (truncate if already exists) */
  datafile = H5Fcreate (DATAFILENAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (datafile < 0)
  {
    fprintf (stderr, "Could not create output file '%s'\n", DATAFILENAME);
    return (-1);
  }

  /* set the dimensions of our sample data array
     count the number of elements */
  elements = 1;
  for (i = 0; i < NDIM; i++)
  {
    dims[i] = NSIZE;
    elements *= NSIZE;
  }
  /* allocate the data array
     we are lazy here and only initialize it to zero */
  data = calloc (elements, sizeof (double));

  CHECK_ERROR (stringtype = H5Tcopy (H5T_C_S1));


  /**************************************************************************/
  /* create a group with attributes describing how the data in this file    */
  /* was written                                                            */
  /**************************************************************************/

  CHECK_ERROR (dataspace = H5Screate (H5S_SCALAR));
  CHECK_ERROR (group = H5Gcreate (datafile, GLOBAL_ATTRIBUTES_GROUP, 0));

  /* we are writing unchunked data */
  unchunked = 1;
  CHECK_ERROR (attr = H5Acreate (group, "unchunked", H5T_NATIVE_INT,
                                 dataspace, H5P_DEFAULT));
  CHECK_ERROR (H5Awrite (attr, H5T_NATIVE_INT, &unchunked));
  CHECK_ERROR (H5Aclose (attr));

  /* the number of processors isn't really needed here
     (only for chunked data) */
  nprocs = 1;
  CHECK_ERROR (attr = H5Acreate (group, "nprocs", H5T_NATIVE_INT,
                                 dataspace, H5P_DEFAULT));
  CHECK_ERROR (H5Awrite (attr, H5T_NATIVE_INT, &nprocs));
  CHECK_ERROR (H5Aclose (attr));

  /* the number of I/O processors isn't really needed here
     (only for chunked data) */
  ioproc_every = 1;
  CHECK_ERROR (attr = H5Acreate (group, "ioproc_every", H5T_NATIVE_INT,
                                 dataspace, H5P_DEFAULT));
  CHECK_ERROR (H5Awrite (attr, H5T_NATIVE_INT, &ioproc_every));
  CHECK_ERROR (H5Aclose (attr));

  /* the version of Cactus that this datafile is (upwards) compatible with */
  CHECK_ERROR (H5Tset_size (stringtype, strlen (CACTUS_VERSION)));
  CHECK_ERROR (attr = H5Acreate (group, "Cactus version", stringtype,
                                 dataspace, H5P_DEFAULT));
  CHECK_ERROR (H5Awrite (attr, stringtype, CACTUS_VERSION));
  CHECK_ERROR (H5Aclose (attr));

  CHECK_ERROR (H5Gclose (group));
  CHECK_ERROR (H5Sclose (dataspace));


  /**************************************************************************/
  /* write your data as a dataset into the file                             */
  /* the dataset name template has the format                               */
  /*   "<full_varname> timelevel <timelevel> at iteration <iteration>"      */
  /**************************************************************************/

  /* the name of the grid variable (as specified in the interface.ccl file) */
  varname = "grid::x";
  /* the timelevel of the variable (set to 0 if there's just one timelevel) */
  timelevel = 0;
  /* iteration isn't really used here */
  iteration = 0;
  sprintf (datasetname, "%s timelevel %d at iteration %d",
           varname, timelevel, iteration);

  /* dataspace is given by dims[] */
  CHECK_ERROR (dataspace = H5Screate_simple (NDIM, dims, NULL));
  CHECK_ERROR (dataset = H5Dcreate (datafile, datasetname, H5T_NATIVE_DOUBLE,
                                    dataspace, H5P_DEFAULT));

  /* write the data */
  CHECK_ERROR (H5Dwrite (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                         H5P_DEFAULT, data));
  CHECK_ERROR (H5Sclose (dataspace));


  /**************************************************************************/
  /* add the necessary attributes describing the dataset                    */
  /* as a Cactus grid variable                                              */
  /**************************************************************************/
  /* the variable's group name (as specified in the interface.ccl file) */
  groupname = "grid::coordinates";
  CHECK_ERROR (dataspace = H5Screate (H5S_SCALAR));
  CHECK_ERROR (H5Tset_size (stringtype, strlen (groupname)));
  CHECK_ERROR (attr = H5Acreate (dataset, "groupname", stringtype,
                                 dataspace, H5P_DEFAULT));
  CHECK_ERROR (H5Awrite (attr, stringtype, groupname));
  CHECK_ERROR (H5Aclose (attr));

  /* the variable's group type (as specified in the interface.ccl file) */
  grouptype = CCTK_GF;
  CHECK_ERROR (attr = H5Acreate (dataset, "grouptype", H5T_NATIVE_INT,
                                 dataspace, H5P_DEFAULT));
  CHECK_ERROR (H5Awrite (attr, H5T_NATIVE_INT, &grouptype));
  CHECK_ERROR (H5Aclose (attr));

  /* the number of timelevels of the variable (as specified in the
     interface.ccl file) */
  ntimelevels = 1;
  CHECK_ERROR (attr = H5Acreate (dataset, "ntimelevels", H5T_NATIVE_INT,
                                 dataspace, H5P_DEFAULT));
  CHECK_ERROR (H5Awrite (attr, H5T_NATIVE_INT, &ntimelevels));
  CHECK_ERROR (H5Aclose (attr));
  CHECK_ERROR (H5Sclose (dataspace));


  /* close remaining HDF5 objects and free allocated resources */
  CHECK_ERROR (H5Dclose (dataset));
  CHECK_ERROR (H5Tclose (stringtype));
  CHECK_ERROR (H5Fclose (datafile));
  free (data);

  return (0);
}
