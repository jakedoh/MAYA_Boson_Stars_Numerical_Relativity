 /*@@
   @file      hdf5_convert_from_ieeeio.c
   @date      Fri 14 Dec 2001
   @author    Thomas Radke
   @desc
              Utility program to convert IOFlexIO datafiles into IOHDF5
              datafiles.
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

/* FlexIO includes */
#include "IOProtos.h"
#include "IEEEIO.h"

/* HDF5 include */
#define H5_USE_16_API 1
#include <hdf5.h>

#define GLOBAL_ATTRIBUTES_GROUP "Global Attributes"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusPUGHIO_IOHDF5_util_hdf5_convert_from_ieeeio_c)


int main (int argc, char **argv)
{
  int i, j;
  IOFile infile;
  int nDatasets, nAttributes;
  int ieeeDatatype, rank, dims[MAXDIM];
  int attrNumberType;
  Long attrLen;
  char attrName[MAXNAMESIZE];
  void *data, *attrData;
  hid_t outfile, dataset, dataspace, attribute, attrDataspace, hdf5Datatype;
  hid_t group, hdf5String;
  hsize_t hdf5Dims[MAXDIM];
  int iteration, timelevel;
  char ieeeDatasetName[MAXNAMESIZE], hdf5DatasetName[2 * MAXNAMESIZE];


  if (argc <= 2) {
    printf ("Usage: %s <inputfile> <outputfile>\n", argv[0]);
    printf ("   eg. %s foo.ieee bar.h5\n", argv[0]);
    return (0);
  }

  infile = IEEEopen (argv[1], "r");
  if (! IOisValid (infile)) {
    printf ("Could not open input file '%s'\n", argv[1]);
    return (-1);
  }

  outfile = H5Fcreate (argv[2], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (outfile < 0)
  {
    fprintf (stderr, "Could not create output file '%s'\n", argv[2]);
    return (-1);
  }

  /* add a dummy GLOBAL_ATTRIBUTES_GROUP so that the HDF5 file is recognized as
     unchunked Cactus data */
  group = H5Gcreate (outfile, GLOBAL_ATTRIBUTES_GROUP, 0);
  if (group < 0)
  {
    fprintf ( stderr, "IOHDF5:hdf5_convert_from_ieeeio.c: group is out of bounds");
    return (-1);
  }
  attrDataspace = H5Screate (H5S_SCALAR);

  attribute = H5Acreate (group, "nprocs", H5T_NATIVE_INT, attrDataspace,
                         H5P_DEFAULT);
  if (attribute < 0)
  {
    fprintf ( stderr, "IOHDF5:hdf5_convert_from_ieeeio.c: attribute is out of bounds");
    return (-1);
  }
  i = 1;
  if(! (H5Awrite (attribute, H5T_NATIVE_INT, &i) >= 0))
  {
    fprintf ( stderr, "IOHDF5:hdf5_convert_from_ieeeio.c: error returned from H5Awrite");
    return (-1);
  }
  H5Aclose (attribute);

  attribute = H5Acreate (group, "ioproc_every", H5T_NATIVE_INT, attrDataspace,
                         H5P_DEFAULT);
  if (attribute < 0)
  {
    fprintf ( stderr, "IOHDF5:hdf5_convert_from_ieeeio.c: attribute is out of bounds");
    return (-1);
  }
  i = 1;
  if(! (H5Awrite (attribute, H5T_NATIVE_INT, &i) >= 0))
  {
    fprintf ( stderr, "IOHDF5:hdf5_convert_from_ieeeio.c: error returned from function H5Awrite");
    return (-1);
  }
  H5Aclose (attribute);

  attribute = H5Acreate (group, "unchunked", H5T_NATIVE_INT, attrDataspace,
                         H5P_DEFAULT);
  if (attribute < 0)
  {
    fprintf ( stderr, "IOHDF5:hdf5_convert_from_ieeeio.c: attribute is out of bounds");
    return (-1);
  }
  i = 1;
  if(! (H5Awrite (attribute, H5T_NATIVE_INT, &i) >= 0))
  {
    fprintf ( stderr, "IOHDF5:hdf5_convert_from_ieeeio.c: error returned from function H5Awrite");
    return (-1);
  }
  H5Aclose (attribute);

  H5Sclose (attrDataspace);
  H5Gclose (group);

  hdf5String = H5Tcopy (H5T_C_S1);

  nDatasets = IOnDatasets (infile);
  printf ("Input file contains %d datasets\n", nDatasets);

  for (i = 0; i < nDatasets; i++)
  {
    if (IOreadInfo (infile, &ieeeDatatype, &rank, dims, MAXDIM) <= 0)
    {
      fprintf (stderr, "Cannot read info of datatset %d, skipping ...\n", i);
      continue;
    }

    /* retrieve name and iteration attribute of dataset */
    strcpy (attrName, "name");
    j = IOreadAttributeInfo (infile, attrName, &attrNumberType, &attrLen);
    if (j < 0)
    {
      sprintf (ieeeDatasetName, "Dataset %d", i);
      fprintf (stderr,
               "Cannot find name attribute of dataset %d, using default name '%s' ...\n", i, ieeeDatasetName);
    }
    else if (IOreadAttribute (infile, j, ieeeDatasetName) < 0)
    {
      fprintf (stderr,
               "Cannot read name attribute of dataset %d, skipping ...\n", i);
      continue;
    }
    strcpy (attrName, "iteration");
    j = IOreadAttributeInfo (infile, attrName, &attrNumberType, &attrLen);
    if (j < 0)
    {
      iteration = 0;
      fprintf (stderr, "Cannot find iteration attribute of dataset %d, "
               "using default iteration number 0 ...\n", i);
    }
    else if (IOreadAttribute (infile, j, &iteration) < 0)
    {
      fprintf (stderr, "Cannot read iteration attribute of dataset %d, "
               "skipping ...\n", i);
      continue;
    }
    strcpy (attrName, "timelevel");
    j = IOreadAttributeInfo (infile, attrName, &attrNumberType, &attrLen);
    if (j < 0)
    {
      timelevel = 0;
      fprintf (stderr, "Cannot find timelevel attribute of dataset %d,\n"
               "assuming Cactus 3.x dataset with implicite timelevel 0 "
               "...\n", i);
    }
    else if (IOreadAttribute (infile, j, &timelevel) < 0)
    {
      fprintf (stderr, "Cannot read timelevel attribute of dataset %d, "
               "skipping ...\n", i);
      continue;
    }

    switch (ieeeDatatype)
    {
      case BYTE:
      case CHAR:
        hdf5Datatype = H5T_NATIVE_CHAR; break;

      case FLOAT32:
        hdf5Datatype = H5T_NATIVE_FLOAT; break;

      case FLOAT64:
        hdf5Datatype = H5T_NATIVE_DOUBLE; break;

      case INT32:
        hdf5Datatype = H5T_NATIVE_INT; break;

      default:
        fprintf (stderr, "Unknown datatype %d for dataset %d, skipping ...\n",
                 ieeeDatatype, i);
        continue;
    }

    /* convert ordering (FlexIO uses fortran order, HDF5 uses C order) */
    for (j = 0; j < rank; j++)
    {
      hdf5Dims[j] = dims[rank - j - 1];
    }
    dataspace = H5Screate_simple (rank, hdf5Dims, NULL);
    if (dataspace < 0)
    {
      fprintf ( stderr, "IOHDF5:hdf5_convert_from_ieeeio.c: dataspace is out of bounds");
      return (-1);
    }

    sprintf (hdf5DatasetName, "/%s timelevel %d at iteration %d",
             ieeeDatasetName, timelevel, iteration);
    dataset = H5Dcreate (outfile, hdf5DatasetName, hdf5Datatype, dataspace,
                         H5P_DEFAULT);
    if (dataset < 0)
    {
      fprintf ( stderr, "IOHDF5:hdf5_convert_from_ieeeio.c: dataset is out of bounds");
      return (-1);
    }

    data = malloc (IOnBytes (ieeeDatatype, rank, dims));
    if (data == NULL)
    {
      fprintf (stderr, "Could not allocate %d bytes to read in dataset\n",
               IOnBytes (ieeeDatatype, rank, dims));
      return (-1);
    }
    IOread (infile, data);
    if(! (H5Dwrite (dataset, hdf5Datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      data) >= 0))
    {
      fprintf ( stderr, "IOHDF5:hdf5_convert_from_ieeeio.c: error returned from function H5Dwrite");
      return (-1);
    }

    nAttributes = IOnAttributes (infile);
    printf ("Processing dataset '%s' with %d attributes\n",
            ieeeDatasetName, nAttributes);

    for (j = 0; j < nAttributes; j++)
    {
      int len;

      if (IOreadIndexedAttributeInfo (infile, j, attrName, &attrNumberType,
                                      &attrLen, MAXNAMESIZE) < 0)
      {
#if  0
        /* it's only a bug if IOreadIndexedAttributeInfo() returns -1 */
        printf ("Cannot read info of attribute %d, skipping ...\n", j);
        continue;
#endif
      }

      len = (int) attrLen;
      /* I don't know why but I have to add one byte - otherwise
         it crashes sometimes when freeing this memory */
      attrData = malloc (1 + IOnBytes (attrNumberType, 1, &len));
      if (IOreadAttribute (infile, j, attrData) < 0)
      {
        fprintf (stderr, "Cannot read value of attribute %d, skipping ...\n",j);
        continue;
      }

      switch (attrNumberType)
      {
        case BYTE:
        case CHAR:
          hdf5Datatype = hdf5String;
          if (len > 1)
          {
            len--;              /* cut trailing '\0' */
          }
          H5Tset_size (hdf5Datatype, len);
          attrDataspace = H5Screate (H5S_SCALAR);
          break;

        case FLOAT64:
          hdf5Datatype = H5T_NATIVE_DOUBLE;
          hdf5Dims[0] = len;
          attrDataspace = H5Screate_simple (1, hdf5Dims, NULL);
          break;

        case INT32:
          hdf5Datatype = H5T_NATIVE_INT;
          hdf5Dims[0] = len;
          attrDataspace = H5Screate_simple (1, hdf5Dims, NULL);
          break;

        default:
          fprintf (stderr, "Unknown datatype %d for attribute %s, skipping "
                   "...\n", attrNumberType, attrName);
          continue;
      }

      if (attrDataspace < 0)
      {
        fprintf ( stderr, "IOHDF5:hdf5_convert_from_ieeeio.c: attrDataspace is out of bounds");
        return (-1);
      }

      /* turn off error messages about an already existing attribute */
      H5E_BEGIN_TRY
      {
        H5Adelete (dataset, attrName);
      } H5E_END_TRY;

      attribute = H5Acreate (dataset, attrName, hdf5Datatype,
                             attrDataspace, H5P_DEFAULT);
      if (attribute < 0)
      {
        fprintf ( stderr, "IOHDF5:hdf5_convert_from_ieeeio.c: attribute is out of bounds");
        return (-1);
      }
      if(! (H5Awrite (attribute, hdf5Datatype, attrData) >= 0))
      {
        fprintf ( stderr, "IOHDF5:hdf5_convert_from_ieeeio.c: error returned from function H5Awrite");
        return (-1);
      }
      H5Sclose (attrDataspace);
      H5Aclose (attribute);

      free (attrData);
    }

    H5Dclose (dataset);
    H5Sclose (dataspace);
    free (data);
  }

  H5Tclose (hdf5String);
  H5Fclose (outfile);
  IOclose (infile);

  return (0);
}
