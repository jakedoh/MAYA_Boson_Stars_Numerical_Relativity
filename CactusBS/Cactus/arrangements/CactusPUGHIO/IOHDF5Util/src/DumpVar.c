/*@@
   @file      DumpVar.c
   @date      Fri May 21 1999
   @author    Thomas Radke
   @desc
              Routines for writing variables into HDF5 data or checkpoint files.
              These routines are used by other HDF5 I/O methods.
   @enddesc
   @version   $Id$
 @@*/


#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "util_Table.h"
#include "cctk_Parameters.h"
#include "CactusPUGH/PUGH/src/include/pugh.h"
#include "CactusBase/IOUtil/src/ioGH.h"
#include "ioHDF5UtilGH.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusPUGHIO_IOHDF5Util_DumpVar_c)


/********************************************************************
 ********************    Macro Definitions   ************************
 ********************************************************************/
/* uncomment the following line to enable some debugging output */
/* #define IOHDF5UTIL_DEBUG 1 */

/* tag base for MPI messages (this may break on more than 2000 processors) */
#define MPITAGBASE 20000


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
static int WriteGS (const cGH *GH, const ioRequest *request, const char *name,
                    hid_t file);
static int WriteGA (const cGH *GH, const ioRequest *request, const char *name,
                    hid_t file);
static void WriteData (const cGH *GH, const ioRequest *request,const char *name,
                       const void *data, int proc, hid_t file);
#if defined(CCTK_MPI) && defined(H5_HAVE_PARALLEL)
static void WriteDataCollective (const cGH *GH, const ioRequest *request,
                                 const char *name, const void *data,hid_t file);
#endif


/*@@
   @routine    IOHDF5Util_DumpVar
   @date       16 Apr 1999
   @author     Thomas Radke
   @desc
               Generic dump routine, just calls the type-specific dump routine.
   @enddesc

   @calls      WriteGS
               WriteGA

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
   @var        file
   @vdesc      the HDF5 file handle
   @vtype      hid_t
   @vio        in
   @endvar

   @returntype int
   @returndesc
               -1 if variable has unknown group type
               or return code of either @seeroutine WriteGS or
                 @seeroutine WriteGA
   @endreturndesc
@@*/
int IOHDF5Util_DumpVar (const cGH *GH, const ioRequest *request, hid_t file)
{
  int gtype, retval;
  char *fullname, *objectname;


  /* build the unique name for the file object to write */
  fullname = CCTK_FullName (request->vindex);
  objectname = malloc (strlen (fullname) + 80);
  sprintf (objectname, "%s timelevel %d at iteration %d",
           fullname, request->timelevel, GH->cctk_iteration);

  /* check whether the object already exists */
  if (request->check_exist && file >= 0)
  {
    H5E_BEGIN_TRY
    {
      H5Gunlink (file, objectname);
    } H5E_END_TRY;
  }

  /* branch to the appropriate dump routine */
  gtype = CCTK_GroupTypeFromVarI (request->vindex);
  if (gtype == CCTK_SCALAR)
  {
    retval = WriteGS (GH, request, objectname, file);
  }
  else if (gtype == CCTK_ARRAY || gtype == CCTK_GF)
  {
    retval = WriteGA (GH, request, objectname, file);
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Invalid group type %d for variable '%s'", gtype, fullname);
    retval = -1;
  }

  /* clean up */
  free (objectname);
  free (fullname);

  return (retval);
}


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
/*@@
   @routine    WriteGS
   @date       May 21 1999
   @author     Thomas Radke
   @desc
               Writes a CCTK_SCALAR type variable into a HDF5 file.
               Only the value from processor 0 is written here.
               If other processors have different values for this scalar
               variable a warning will be issued.
   @enddesc

   @calls      IOHDF5Util_DumpCommonAttributes

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
   @var        name
   @vdesc      name of the dataset to write
   @vtype      const char *
   @vio        in
   @endvar
   @var        file
   @vdesc      the HDF5 file to dump to
   @vtype      hid_t
   @vio        in
   @endvar

   @returntype int
   @returndesc
               0 for success, or -1 if file handle is invalid
   @endreturndesc
@@*/
static int WriteGS (const cGH *GH, const ioRequest *request, const char *name,
                    hid_t file)
{
  int i, myproc, nprocs, hdatatypesize, retval;
  char *buffer, *fullname;
  const ioGH *ioUtilGH;
  const ioHDF5UtilGH *myGH;
  hid_t dataset, hdf5type;


  ioUtilGH = CCTK_GHExtension (GH, "IO");
  myproc = CCTK_MyProc (GH);
  nprocs = CCTK_nProcs (GH);
  fullname = CCTK_FullName (request->vindex);

  hdatatypesize = CCTK_VarTypeSize (request->hdatatype);
  buffer = calloc (nprocs, hdatatypesize);
  memcpy (buffer + myproc*hdatatypesize,
          CCTK_VarDataPtrI (GH, request->timelevel, request->vindex),
          hdatatypesize);

  if (nprocs > 1)
  {
    i = CCTK_ReductionHandle ("sum");
    if (i >= 0)
    {
      i = CCTK_ReduceArray (GH, -1, i, nprocs, request->hdatatype,
                            buffer, 1, 1, request->hdatatype, nprocs, buffer);
    }
    if (i < 0)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "WriteGS: Cannot check whether values on different "
                  "processors are the same for grid scalar '%s'", fullname);

      /* copy this processor's value to the start of buffer */
      memcpy (buffer, buffer + myproc*hdatatypesize, hdatatypesize);
    }
    else
    {
      retval = 0;
      for (i = 1; i < nprocs; i++)
      {
        retval |= memcmp (buffer, buffer + i*hdatatypesize, hdatatypesize);
      }
      if (retval)
      {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "WriteGS: value of grid scalar variable '%s' (timelevel %d)"
                    " differs between processors, only value from processor 0 "
                    "will be written", fullname, request->timelevel);
      }
    }
  }

  /* only I/O processors write data */
  if (myproc != ioUtilGH->ioproc)
  {
    retval = 0;
  }
  else if (file < 0)
  {
    retval = -1;
  }
  else
  {
    myGH = CCTK_GHExtension (GH, "IOHDF5Util");
    hdf5type = IOHDF5Util_DataType (myGH, request->hdatatype);
    HDF5_ERROR (dataset = H5Dcreate (file, name, hdf5type,
                                     myGH->scalar_dataspace, H5P_DEFAULT));
    HDF5_ERROR (H5Dwrite (dataset, hdf5type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          buffer));

    /* scalars have size 0 */
    request->hsize[0] = 0;
    IOHDF5Util_DumpCommonAttributes (GH, request, dataset);

    HDF5_ERROR (H5Dclose (dataset));
    retval = 0;
  }

  free (buffer);
  free (fullname);

  return (retval);
}


/*@@
   @routine    WriteGA
   @date       May 21 1999
   @author     Paul Walker
   @desc
               Writes a grid array into a HDF5 file.
   @enddesc

   @calls      Hyperslab_LocalMappingByIndex
               Hyperslab_FreeMapping
               Hyperslab_Get
               WriteData
               WriteDataCollective

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
   @var        name
   @vdesc      name of the dataset to write
   @vtype      const char *
   @vio        in
   @endvar
   @var        file
   @vdesc      the HDF5 file to dump to
   @vtype      hid_t
   @vio        in
   @endvar

   @returntype int
   @returndesc
                0 for success, or<BR>
               -1 if hyperslab mapping couldn't be defined, or<BR>
               -2 if hyperslab couldn't be extracted
   @endreturndesc
@@*/
static int WriteGA (const cGH *GH, const ioRequest *request, const char *name,
                    hid_t file)
{
  const ioGH *ioUtilGH;
  int i, myproc, mapping, hdatasize, table, retval;
  void *hdata;
  char *fullname;
#ifdef CCTK_MPI
  int nprocs;
  void *tmpd;
  MPI_Comm comm;
  MPI_Status ms;
  MPI_Datatype mpitype;
#endif
  DECLARE_CCTK_PARAMETERS


  /* define the hyperslab mapping */
  table = -1;
  if (request->with_ghostzones)
  {
    table = Util_TableCreateFromString ("with_ghostzones = 1");
    if (table < 0)
    {
      CCTK_WARN (1, "Failed to hyperslab parameter create table from string");
    }
  }

  mapping = Hyperslab_LocalMappingByIndex (GH, request->vindex,
                                           request->hdim,
                                           request->direction,
                                           request->origin,
                                           request->extent,
                                           request->downsample,
                                           table, NULL,
                                           request->hsize_chunk,
                                           request->hsize,
                                           request->hoffset);
  if (table >= 0)
  {
    Util_TableDestroy (table);
  }

  if (mapping < 0)
  {
    fullname = CCTK_FullName (request->vindex);
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Failed to define hyperslab mapping for variable '%s'",
                fullname);
    free (fullname);
    return (-1);
  }

  /* calculate the size of the hyperslab */
  request->hsize_chunk[request->hdim] = 1;
  for (i = 0; i < request->hdim; i++)
  {
    request->hsize_chunk[request->hdim] *= request->hsize_chunk[i];
  }

  /* get the hyperslab */
  hdatasize = CCTK_VarTypeSize (request->hdatatype);
  hdata = request->hsize_chunk[request->hdim] > 0 ?
          malloc (request->hsize_chunk[request->hdim] * hdatasize) : NULL;
  retval = Hyperslab_Get (GH, mapping, -1, request->vindex, request->timelevel,
                          request->hdatatype, hdata);

  /* release the mapping structure */
  Hyperslab_FreeMapping (mapping);

  if (retval)
  {
    fullname = CCTK_FullName (request->vindex);
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Failed to extract hyperslab for variable '%s'", fullname);
    free (fullname);
    if (hdata)
    {
      free (hdata);
    }
    return (-2);
  }

  ioUtilGH = CCTK_GHExtension (GH, "IO");

#ifdef CCTK_MPI
#ifdef H5_HAVE_PARALLEL
  if (ioUtilGH->unchunked)
  {
    WriteDataCollective (GH, request, name, hdata, file);
    if (hdata)
    {
      free (hdata);
    }
    return (0);
  }
#endif
#endif

  /* dump data held on I/O processor */
  myproc = CCTK_MyProc (GH);
  if (myproc == ioUtilGH->ioproc)
  {
    WriteData (GH, request, name, hdata, myproc, file);
  }
#ifdef CCTK_MPI
  nprocs = CCTK_nProcs (GH);
  comm = PUGH_pGH (GH)->PUGH_COMM_WORLD;
  mpitype = PUGH_MPIDataType (PUGH_pGH (GH), request->hdatatype);

  if (myproc == ioUtilGH->ioproc)
  {
    /* dump data from all other processors */
    for (i = myproc+1; i < myproc+ioUtilGH->ioproc_every && i < nprocs; i++)
    {
      /* receive geometry (this assumes the geometry arrays
         to be contiguous starting at hoffset */
      CACTUS_MPI_ERROR (MPI_Recv (request->hoffset, 3*request->hdim + 1,
                                  PUGH_MPI_INT, i, 2*i + MPITAGBASE + 1,
                                  comm, &ms));

      /* receive data */
      tmpd = NULL;
      if (request->hsize_chunk[request->hdim] > 0)
      {
        tmpd = malloc (request->hsize_chunk[request->hdim] * hdatasize);
        CACTUS_MPI_ERROR (MPI_Recv (tmpd, request->hsize_chunk[request->hdim],
                                    mpitype, i, 2*i + MPITAGBASE, comm, &ms));
      }

      /* write data */
      WriteData (GH, request, name, tmpd, i, file);

      if (tmpd)
      {
        free (tmpd);
      }

    } /* end loop over processors */

  }
  else
  {
    /* send geometry (this assumes the geometry arrays to be contiguous
       starting at hoffset) */
    CACTUS_MPI_ERROR (MPI_Send (request->hoffset, 3*request->hdim + 1,
                                PUGH_MPI_INT, ioUtilGH->ioproc,
                                2*myproc + MPITAGBASE + 1, comm));
    /* send data */
    if (request->hsize_chunk[request->hdim] > 0)
    {
      CACTUS_MPI_ERROR (MPI_Send (hdata, request->hsize_chunk[request->hdim],
                                  mpitype, ioUtilGH->ioproc,
                                  2*myproc + MPITAGBASE, comm));
#ifdef IOHDF5UTIL_DEBUG
      printf ("Processor %d sent %d data points\n",
              myproc, request->hsize_chunk[request->hdim]);
#endif
    }
  }

  /* wait for every processor to catch up */
  CCTK_Barrier (GH);
#endif

  /* free allocated resources */
  if (hdata)
  {
    free (hdata);
  }

  return (retval);
}


/*@@
   @routine    WriteData
   @author     Thomas Radke
   @date       May 21 1999
   @desc
               All I/O processors dump the data of their group into the file,
               either as one chunk per processor or unchunked
               (depending on parameter "unchunked").
   @enddesc

   @calls      IOHDF5Util_DumpCommonAttributes

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
   @var        name
   @vdesc      name of the dataset to write
   @vtype      const char *
   @vio        in
   @endvar
   @var        data
   @vdesc      pointer to the chunk to dump
   @vtype      const void *
   @vio        in
   @endvar
   @var        proc
   @vdesc      the processor which's chunk is to be dumped
   @vtype      int
   @vio        in
   @endvar
   @var        file
   @vdesc      the HDF5 file handle
   @vtype      hid_t
   @vio        in
   @endvar
@@*/
static void WriteData (const cGH *GH, const ioRequest *request,const char *name,
                       const void *data, int proc, hid_t file)
{
  DECLARE_CCTK_PARAMETERS
  int i, myproc;
  ioGH *ioUtilGH;
  ioHDF5UtilGH *myGH;
  hid_t hdf5type, group, dataset, memspace, filespace, plist;
  char *chunkname;
#if (H5_VERS_MAJOR == 1 && \
     (H5_VERS_MINOR < 6 || (H5_VERS_MINOR == 6 && H5_VERS_RELEASE < 4)))
  hssize_t *chunk_origin;
#else
  hsize_t *chunk_origin;
#endif
  hsize_t *chunk_dims, *file_dims;
  hsize_t buffersize;
  const int compression_lvl = request->compression_level >= 0 ?
                              request->compression_level : compression_level;


  /* immediately return if file handle is invalid */
  if (file < 0)
  {
    return;
  }

  ioUtilGH = CCTK_GHExtension (GH, "IO");
  myGH = CCTK_GHExtension (GH, "IOHDF5Util");

  /* copy the size arrays from CCTK_INT to appropriate types
     note that HDF5 wants elements in reverse order */
  chunk_origin = malloc (request->hdim * sizeof (*chunk_origin));
  chunk_dims   = malloc (2*request->hdim * sizeof (*chunk_dims));
  file_dims    = chunk_dims + request->hdim;
  for (i = 0; i < request->hdim; i++)
  {
    chunk_origin[i] = request->hoffset[request->hdim - 1 - i];
    file_dims   [i] = request->hsize[request->hdim - 1 - i];
    chunk_dims  [i] = request->hsize_chunk[request->hdim - 1 - i];
  }

  myproc = CCTK_MyProc (GH);

  memspace = -1;
  if (data)
  {
    /* create the memspace according to chunk dims */
    HDF5_ERROR (memspace = H5Screate_simple (request->hdim, chunk_dims, NULL));
  }

  hdf5type = IOHDF5Util_DataType (myGH, request->hdatatype);
  if (ioUtilGH->unchunked)
  {
    /* create the (global) filespace and set the hyperslab for the chunk */
    HDF5_ERROR (filespace = H5Screate_simple (request->hdim, file_dims, NULL));
    HDF5_ERROR (H5Sselect_hyperslab (filespace, H5S_SELECT_SET, chunk_origin,
                                     NULL, chunk_dims, NULL));

    /* the I/O processor creates the dataset and adds the common attributes
       when writing its own data, otherwise the dataset is reopened */
    if (proc == myproc)
    {
      /* enable compression for chunked dataset if compression was requested */
      HDF5_ERROR (plist = H5Pcreate (H5P_DATASET_CREATE));
      if (compression_lvl)
      {
        HDF5_ERROR (H5Pset_chunk (plist, request->hdim, chunk_dims));
        HDF5_ERROR (H5Pset_shuffle (plist));
        HDF5_ERROR (H5Pset_deflate (plist, compression_lvl));
      }
      HDF5_ERROR (dataset = H5Dcreate (file, name, hdf5type, filespace, plist));
      HDF5_ERROR (H5Pclose (plist));
      IOHDF5Util_DumpCommonAttributes (GH, request, dataset);
    }
    else
    {
      HDF5_ERROR (dataset = H5Dopen (file, name));
    }

    if (memspace >= 0)
    {
      /* increase the buffer size if the default isn't sufficient */
      HDF5_ERROR (plist = H5Pcreate (H5P_DATASET_XFER));
      buffersize = H5Dget_storage_size (dataset);
      if (buffersize > H5Pget_buffer (plist, NULL, NULL))
      {
        HDF5_ERROR (H5Pset_buffer (plist, buffersize, NULL, NULL));
      }

      /* write the data */
      HDF5_ERROR (H5Dwrite (dataset, hdf5type, memspace, filespace, plist,
                            data));

      /* close the transfer property list */
      HDF5_ERROR (H5Pclose (plist));
    }
    /* close the file dataspace */
    HDF5_ERROR (H5Sclose (filespace));
  }
  else
  {
    /* the I/O processor creates the chunk group and adds common attributes */
    if (proc == myproc)
    {
      HDF5_ERROR (group = H5Gcreate (file, name, 0));
      IOHDF5Util_DumpCommonAttributes (GH, request, group);
      HDF5_ERROR (H5Gclose (group));
    }

    dataset = -1;
    if (memspace >= 0)
    {
      /* now the chunk dataset for each processor is created within the group */
      chunkname = malloc (strlen (name) + 20);
      sprintf (chunkname, "%s/chunk%d", name, proc - myproc);
      HDF5_ERROR (plist = H5Pcreate (H5P_DATASET_CREATE));
      /* enable compression for chunked dataset if compression was requested */
      if (compression_level)
      {
        HDF5_ERROR (H5Pset_chunk (plist, request->hdim, chunk_dims));
        HDF5_ERROR (H5Pset_shuffle (plist));
        HDF5_ERROR (H5Pset_deflate (plist, compression_level));
      }
      /* create the chunk dataset and dump the chunk data */
      HDF5_ERROR (dataset = H5Dcreate (file, chunkname, hdf5type, memspace,
                                       plist));
      HDF5_ERROR (H5Pclose (plist));
      HDF5_ERROR (H5Dwrite (dataset, hdf5type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            data));

      /* add the "origin" attribute for the chunk */
      WRITE_ATTRIBUTE ("chunk_origin", request->hoffset, dataset, myGH,
                       request->hdim, HDF5_INT);

      free (chunkname);
    }
  }

  /* close the dataset and the memspace */
  if (dataset >= 0)
  {
    HDF5_ERROR (H5Dclose (dataset));
  }
  if (memspace >= 0)
  {
    HDF5_ERROR (H5Sclose (memspace));
  }

  /* free allocated resources */
  free (chunk_origin);
  free (chunk_dims);
}


#if defined(CCTK_MPI) && defined(H5_HAVE_PARALLEL)
/*@@
   @routine    WriteDataCollective
   @author     Thomas Radke
   @date       May 21 1999
   @desc
               All processors dump their data into an unchunked file.
   @enddesc

   @calls      IOHDF5Util_DumpCommonAttributes

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
   @var        name
   @vdesc      name of the dataset to write
   @vtype      const char *
   @vio        in
   @endvar
   @var        data
   @vdesc      pointer to the chunk to dump
   @vtype      const void *
   @vio        in
   @endvar
   @var        file
   @vdesc      the HDF5 file handle
   @vtype      hid_t
   @vio        in
   @endvar
@@*/
static void WriteDataCollective (const cGH *GH, const ioRequest *request,
                                 const char *name, const void *data, hid_t file)
{
  DECLARE_CCTK_PARAMETERS
  int i;
  hid_t hdf5type, dataset, memspace, filespace, plist;
#if (H5_VERS_MAJOR == 1 && \
     (H5_VERS_MINOR < 6 || (H5_VERS_MINOR == 6 && H5_VERS_RELEASE < 4)))
  hssize_t *chunk_origin;
#else
  hsize_t *chunk_origin;
#endif
  hsize_t *chunk_dims, *file_dims;
  hsize_t buffersize;
  const ioHDF5UtilGH *myGH;
  const int compression_lvl = request->compression_level >= 0 ?
                              request->compression_level : compression_level;


  /* immediately return if file handle is invalid */
  if (file < 0)
  {
    return;
  }

  /* copy the size arrays from CCTK_INT to appropriate types
     note that HDF5 wants elements in reverse order */
  chunk_origin = malloc (request->hdim * sizeof (*chunk_origin));
  chunk_dims   = malloc (2*request->hdim * sizeof (*chunk_dims));
  file_dims    = chunk_dims + request->hdim;
  for (i = 0; i < request->hdim; i++)
  {
    file_dims   [i] = request->hsize[request->hdim - 1 - i];
    chunk_origin[i] = request->hoffset[request->hdim - 1 - i];
    chunk_dims  [i] = request->hsize_chunk[request->hdim - 1 - i];
  }

  /* create the memspace according to chunk dims */
  HDF5_ERROR (memspace = H5Screate_simple (request->hdim, chunk_dims, NULL));

  /* create the (global) filespace and set the hyperslab for the chunk */
  HDF5_ERROR (filespace = H5Screate_simple (request->hdim, file_dims, NULL));
  HDF5_ERROR (H5Sselect_hyperslab (filespace, H5S_SELECT_SET, chunk_origin,
                                   NULL, chunk_dims, NULL));

  /* the I/O processor creates the dataset and adds the common attributes
     when writing its own data, otherwise the dataset is reopened */
  myGH = CCTK_GHExtension (GH, "IOHDF5Util");
  hdf5type = IOHDF5Util_DataType (myGH, request->hdatatype);

  /* enable compression for chunked dataset if compression was requested */
  HDF5_ERROR (plist = H5Pcreate (H5P_DATASET_CREATE));
  if (compression_lvl)
  {
    HDF5_ERROR (H5Pset_chunk (plist, request->hdim, chunk_dims));
    HDF5_ERROR (H5Pset_shuffle (plist));
    HDF5_ERROR (H5Pset_deflate (plist, compression_lvl));
  }
  HDF5_ERROR (dataset = H5Dcreate (file, name, hdf5type, filespace, plist));
  HDF5_ERROR (H5Pclose (plist));
  IOHDF5Util_DumpCommonAttributes (GH, request, dataset);

  /* increase the buffer size if the default isn't sufficient */
  HDF5_ERROR (plist = H5Pcreate (H5P_DATASET_XFER));
  HDF5_ERROR (H5Pset_dxpl_mpio (plist, H5FD_MPIO_COLLECTIVE));
  buffersize = H5Dget_storage_size (dataset);
  if (buffersize > H5Pget_buffer (plist, NULL, NULL))
  {
    HDF5_ERROR (H5Pset_buffer (plist, buffersize, NULL, NULL));
  }

  /* write the data */
  HDF5_ERROR (H5Dwrite (dataset, hdf5type, memspace, filespace, plist, data));

  /* close resources */
  HDF5_ERROR (H5Pclose (plist));
  HDF5_ERROR (H5Sclose (filespace));
  HDF5_ERROR (H5Dclose (dataset));
  HDF5_ERROR (H5Sclose (memspace));

  /* free allocated resources */
  free (chunk_origin);
  free (chunk_dims);
}
#endif /* CCTK_MPI && H5_HAVE_PARALLEL */
