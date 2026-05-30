 /*@@
   @file      RecoverVar.c
   @date      Thu Jun 18 16:34:59 1998
   @author    Tom Goodale
   @desc 
              Routines to recover variables from a given HDF5 data or
              checkpoint file.
              These routines are used by other HDF5 IO methods.
   @enddesc 
   @version   $Id$
 @@*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "util_String.h"
#include "CactusBase/IOUtil/src/ioGH.h"
#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"
#include "CactusPUGH/PUGH/src/include/pugh.h"
#include "ioHDF5UtilGH.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusPUGHIO_IOHDF5Util_RecoverVar_c)


/********************************************************************
 ********************    Macro Definitions   ************************
 ********************************************************************/
/* tag base for MPI messages */
#define MPITAGBASE 1001        


/********************************************************************
 ********************    Internal Typedefs   ************************
 ********************************************************************/
typedef struct
{
  const cGH *GH;
  int ioproc;
  int ioproc_every;
  int unchunked;
  int has_version;
  int num_recovered;
} iterate_info_t;

typedef struct
{
  iterate_info_t *it_info;
  int element_size;
  hid_t hdf5type;
#ifdef CCTK_MPI
  MPI_Datatype mpi_type;
#endif
} recover_info_t;


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
static herr_t processDataset (hid_t group, const char *datasetname, void *arg);


 /*@@
   @routine    IOHDF5Util_RecoverVariables
   @date       Fri Jun 19 09:19:48 1998
   @author     Tom Goodale
   @desc 
               Reads in data from an open HDF5 file.
   @enddesc 
   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      cGH *
   @vio        in
   @endvar
   @var        fileinfo
   @vdesc      pointer to info structure describing the HDF5 file
   @vtype      const fileinfo_t *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               the total number of recovered variables
   @endreturndesc
@@*/
int IOHDF5Util_RecoverVariables (cGH *GH, const fileinfo_t *fileinfo)
{
  iterate_info_t info;
#ifdef CCTK_MPI
  pGH *pughGH;
  CCTK_INT var_info[3];
  MPI_Status ms;
  MPI_Datatype mpi_type;
  int vindex, timelevel, proc, npoints;
#endif
  DECLARE_CCTK_PARAMETERS


  info.GH = GH;
  info.ioproc = fileinfo->ioproc;
  info.unchunked = fileinfo->unchunked;
  info.ioproc_every = fileinfo->ioproc_every;
  info.has_version = fileinfo->has_version;
  info.num_recovered = 0;

#ifdef CCTK_MPI
  pughGH = PUGH_pGH (GH);

  /* now process the datasets:
     All IO processors read the datasets from their checkpoint file,
     verify their contents and communicate them to the non-I/O processors. */

  /* At first the code for the IO processors.
     This holds also for the single processor case. */
  if (CCTK_MyProc (GH) == fileinfo->ioproc)
  {
#endif /* CCTK_MPI */

    /* iterate over all datasets starting from "/" in the HDF5 file */
    HDF5_ERROR (H5Giterate (fileinfo->file, "/", NULL, processDataset, &info));

#ifdef CCTK_MPI
    /* To signal completion to the non-IO processors
       an invalid variable index is broadcast. */
    var_info[0] = -1;
    var_info[1] = info.num_recovered;
    for (proc = 1; proc < fileinfo->ioproc_every; proc++)
    for (proc = fileinfo->ioproc + 1;
         proc < fileinfo->ioproc + fileinfo->ioproc_every &&
         proc < CCTK_nProcs (GH);
         proc++)
    {
      CACTUS_MPI_ERROR (MPI_Send (var_info, 3, PUGH_MPI_INT, proc,
                        MPITAGBASE, pughGH->PUGH_COMM_WORLD));
    }
  }
  else
  {
    /* And here the code for non-I/O processors: */
    /* They don't know how many datasets there are, because the I/O processors
       could skip some on the fly during their consistency checks.
       The I/O Processor sends the index of the variable to be processed next.
       So, all non-I/O processors execute a loop where the termination condition
       is when an invalid index was received.
    */
    while (1)
    {
      /* receive the next variable index from my IO processor */
      CACTUS_MPI_ERROR (MPI_Recv (var_info, 3, PUGH_MPI_INT, fileinfo->ioproc,
                        MPITAGBASE, pughGH->PUGH_COMM_WORLD, &ms));
      vindex = var_info[0]; timelevel = var_info[1]; npoints = var_info[2];

      /* check for termination condition */
      if (vindex < 0)
      {
        info.num_recovered = var_info[1];
        break;
      }

      mpi_type = PUGH_MPIDataType (pughGH, CCTK_VarTypeI (vindex));
      if (! mpi_type)
      {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Unsupported datatype %d", CCTK_VarTypeI (vindex));
        continue;
      }

      /* receive following data from my IO processor */
      if (npoints > 0)
      {
        CACTUS_MPI_ERROR (MPI_Recv (GH->data[vindex][timelevel], npoints,
                                    mpi_type, fileinfo->ioproc, MPITAGBASE,
                                    pughGH->PUGH_COMM_WORLD, &ms));
      }
    }
  }
#endif /* CCTK_MPI */

  return (info.num_recovered);
}


/* NOTE: Although we could read the GH extensions from multiple recovery files
         in parallel, this is done only on by processor 0 here.
         Broadcasting the GH extensions is found faster than
         sending it in a loop from each I/O processor to all the non I/Os
         (don't have subcommunicators yet) */
int IOHDF5Util_RecoverGHextensions (cGH *GH, const fileinfo_t *fileinfo)
{
  hid_t group;
  CCTK_REAL realBuffer[2];
  CCTK_INT4 int4Buffer[3];


  if (CCTK_MyProc (GH) == 0)
  {
    /* all the important global attributes and GH extensions
       are stored in the GLOBAL_ATTRIBUTES_GROUP group */
    group = H5Gopen (fileinfo->file, GLOBAL_ATTRIBUTES_GROUP);
    int4Buffer[0] = group >= 0;

    if (int4Buffer[0])
    {
      READ_ATTRIBUTE (group, "cctk_iteration", HDF5_INT4, &int4Buffer[1]);
      READ_ATTRIBUTE (group, "main_loop_index", HDF5_INT4, &int4Buffer[2]);
      READ_ATTRIBUTE (group, "cctk_time", HDF5_REAL, &realBuffer[0]);
      // pre rev 173 checkpoints do not contain cctk_delta_time. When reading
      // such a checkpoint, we output a warning and continue to use the old
      // cctk_delta_time. For this we need to initialize realBuffer[1] since we
      // cannot tell if READ_ATTRIBUTE succeeded or not.
      realBuffer[1] = GH->cctk_delta_time;
      READ_ATTRIBUTE (group, "cctk_delta_time", HDF5_REAL, &realBuffer[1]);

      HDF5_ERROR (H5Gclose (group));
    }
    else
    {
      CCTK_WARN (1, "Can't find global attributes group. "
                    "Is this really a Cactus HDF5 datafile ?");
    }
  }

#ifdef CCTK_MPI
  /* Broadcast the GH extensions to all processors */
  /* NOTE: We have to use MPI_COMM_WORLD here 
     because PUGH_COMM_WORLD is not yet set up at parameter recovery time.
     We also assume that PUGH_MPI_INT4 is a compile-time defined datatype. */
  CACTUS_MPI_ERROR (MPI_Bcast (int4Buffer, 3, PUGH_MPI_INT4, 0,MPI_COMM_WORLD));
  if (int4Buffer[0])
  {
    CACTUS_MPI_ERROR (MPI_Bcast (realBuffer, 2, PUGH_MPI_REAL, 0,
                                 MPI_COMM_WORLD));
  }
#endif

  if (int4Buffer[0])
  {
    GH->cctk_time = realBuffer[0];
    GH->cctk_delta_time = realBuffer[1];
    GH->cctk_iteration = int4Buffer[1];
    CCTK_SetMainLoopIndex ((int) int4Buffer[2]);
  }

  /* return 0 for success otherwise negative */
  return (int4Buffer[0] ? 0 : -1);
}


/* NOTE: Although we could read the parameters from multiple recovery files
         in parallel, this is done only on by processor 0 here.
         Broadcasting the complete parameter string is found faster than
         sending it in a loop from each IO processor to all the non IOs
         (don't have subcommunicators yet) */
int IOHDF5Util_RecoverParameters (const fileinfo_t *fileinfo)
{
  hid_t group, attr, dataset, atype;
  char *parameters;
  CCTK_INT4 parameterSize;
  DECLARE_CCTK_PARAMETERS


  /* To make the compiler happy */
  group = attr = dataset = -1;
  parameters = NULL;
  parameterSize = 0;

  if (CCTK_MyProc (NULL) == 0)
  {
    if (CCTK_Equals (verbose, "full"))
    {
      CCTK_VInfo (CCTK_THORNSTRING, "Recovering parameters from checkpoint "
                                    "file '%s'", fileinfo->filename);
    }

    /* the parameters are stored in the CACTUS_PARAMETERS_GROUP group
       in the attribute (old style) or dataset (new style) ALL_PARAMETERS */
    group = H5Gopen (fileinfo->file, CACTUS_PARAMETERS_GROUP);
    if (group > 0)
    {
      H5E_BEGIN_TRY
      {
        attr = H5Aopen_name (group, ALL_PARAMETERS);
        dataset = H5Dopen (group, ALL_PARAMETERS);
      } H5E_END_TRY
    }

    if (group >= 0 && (attr >= 0 || dataset >= 0))
    {
      if (attr >= 0)
      {
        HDF5_ERROR (atype = H5Aget_type (attr));
        parameterSize = H5Tget_size (atype);
        parameters = malloc (parameterSize + 1);
        HDF5_ERROR (H5Aread (attr, atype, parameters));
        parameters[parameterSize] = 0;
        HDF5_ERROR (H5Tclose (atype));
      }
      else
      {
        parameterSize = H5Dget_storage_size (dataset);
        parameters = malloc (parameterSize);
        HDF5_ERROR (H5Dread (dataset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL,
                             H5P_DEFAULT, parameters));
      }
    }
    else
    {
      CCTK_WARN (1, "Can't find global parameters. "
                    "Is this really a Cactus HDF5 checkpoint file ?");
    }

    if (attr >= 0)
    {
      HDF5_ERROR (H5Aclose (attr));
    }
    if (dataset >= 0)
    {
      HDF5_ERROR (H5Dclose (dataset));
    }
    if (group >= 0)
    {
      HDF5_ERROR (H5Gclose (group));
    }
  }

#ifdef CCTK_MPI
  /* Broadcast the parameter buffer size to all processors */
  /* NOTE: We have to use MPI_COMM_WORLD here 
     because PUGH_COMM_WORLD is not yet set up at parameter recovery time.
     We also assume that PUGH_MPI_INT4 is a compile-time defined datatype. */
  CACTUS_MPI_ERROR (MPI_Bcast (&parameterSize, 1, PUGH_MPI_INT4, 0,
                    MPI_COMM_WORLD));
#endif

  if (parameterSize > 0)
  {
#ifdef CCTK_MPI
    if (CCTK_MyProc (NULL) != 0)
    {
      parameters = malloc (parameterSize + 1);
    }

    CACTUS_MPI_ERROR (MPI_Bcast (parameters, parameterSize + 1, PUGH_MPI_CHAR,
                      0, MPI_COMM_WORLD));
#endif

    IOUtil_SetAllParameters (parameters);

    free (parameters);
  }

  /* return positive value for success otherwise negative */
  return (parameterSize > 0 ? 1 : -1);
}


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
/* local routine GetCommonAttributes() reads in the next dataset's attributes
   and verifies them:

  * checks if variable <vindex> belongs to the same group
  * checks the group data info:
    - group type
    - variable type
    - ntimelevels
    - sizes (rank, dimensions) according to chunking mode

 If there is a mismatch a warning (warning level 2) is printed and
 a negative value is returned to indicate that this dataset should be ignored.
 If successful, the group type and the possibly corrected timelevel
 to restore are stored in {*gtype, *timelevel}, and 0 is returned.
*/
static int GetCommonAttributes (const cGH *GH,
                                hid_t dataset,
                                int vindex,
                                const char *fullname,
                                int unchunked,
                                int *grouptype,
                                int *timelevel,
                                int is_group,
                                int has_version)
{
  cGroup group_static_data;
  cGroupDynamicData group_dynamic_data;
  int i, flag;
  const int *dims;
  int groupindex;
  hid_t datatype, dataspace;
  hsize_t rank_stored, *dims_stored;
  int grouptype_stored, numtimelevels_stored;
  char *groupname, *msg, *oldmsg;
  char groupname_stored[128];


  /* read and verify the group name */
  READ_ATTRIBUTE (dataset, "groupname", H5T_C_S1, groupname_stored);
  groupname = CCTK_GroupNameFromVarI (vindex);
  i = CCTK_Equals (groupname_stored, groupname);
  free (groupname);
  if (! i)
  {
    CCTK_WARN (2, "Groupnames don't match");
    return (-1);
  }

  /* verify group type, variable type, dims, sizes and ntimelevels */
  READ_ATTRIBUTE (dataset, "grouptype", H5T_NATIVE_INT, &grouptype_stored);
  /* be backwards compatible */
  switch (grouptype_stored)
  {
  case 1: grouptype_stored = CCTK_SCALAR; break;
  case 2: grouptype_stored = CCTK_GF; break;
  case 3: grouptype_stored = CCTK_ARRAY; break;
  }
  READ_ATTRIBUTE (dataset, "ntimelevels", H5T_NATIVE_INT,&numtimelevels_stored);

  /* get the group data */
  groupindex = CCTK_GroupIndex (groupname_stored);
  if (CCTK_GroupData (groupindex, &group_static_data) != 0)
  {
    CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot get static group data for '%s'", fullname);
    return (-1);
  }

  /* now check the group data against the information in the checkpoint file */
  if (group_static_data.grouptype != grouptype_stored)
  {
    CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Group types don't match for '%s'", fullname);
    return (-1);
  }
  if (group_static_data.numtimelevels != numtimelevels_stored)
  {
    CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Number of timelevels don't match for '%s'", fullname);
    return (-1);
  }
  /* increment the timelevel for data from old checkpoint files */
  if (! has_version && group_static_data.numtimelevels > 1)
  {
    (*timelevel)++;
  }
  /* open the first chunk to determine datatype, dims and sizes
     if the dataset is a chunk group */
  if (is_group)
  {
    HDF5_ERROR (dataset = H5Dopen (dataset, "chunk0"));
  }
  HDF5_ERROR (datatype = H5Dget_type (dataset));

  /* The CCTK variable type defines do not correlate with the HDF5 defines
     so compare them explicitely here. */
  flag = (H5Tget_class (datatype) == H5T_FLOAT &&
          strncmp (CCTK_VarTypeName (group_static_data.vartype),
                   "CCTK_VARIABLE_REAL", 18) == 0) ||
         (H5Tget_class (datatype) == H5T_INTEGER &&
          (strncmp (CCTK_VarTypeName (group_static_data.vartype),
                    "CCTK_VARIABLE_INT", 17) == 0 ||
           strcmp (CCTK_VarTypeName (group_static_data.vartype),
                   "CCTK_VARIABLE_BYTE") == 0)) ||
         (H5Tget_class (datatype) == H5T_COMPOUND &&
          strncmp (CCTK_VarTypeName (group_static_data.vartype),
                   "CCTK_VARIABLE_COMPLEX", 21) == 0);

  HDF5_ERROR (H5Tclose (datatype));
  if (! flag)
  {
    CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Variable types don't match for '%s'", fullname);
    return (-1);
  }

  /* verify the dims and sizes */
  HDF5_ERROR (dataspace = H5Dget_space (dataset));
  HDF5_ERROR (rank_stored = H5Sget_simple_extent_ndims (dataspace));
  dims_stored = NULL;
  if (rank_stored > 0)
  {
    dims_stored = malloc (rank_stored * sizeof (hsize_t));
    HDF5_ERROR (H5Sget_simple_extent_dims (dataspace, dims_stored, NULL));
  }
  HDF5_ERROR (H5Sclose (dataspace));

  flag = group_static_data.dim != (int) rank_stored;
  if (group_static_data.grouptype == CCTK_ARRAY ||
      group_static_data.grouptype == CCTK_GF)
  {
    if (CCTK_GroupDynamicData (GH, groupindex, &group_dynamic_data) != 0)
    {
      CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Cannot get dynamic group data for '%s'", fullname);
      return (-1);
    }
    dims = unchunked ? group_dynamic_data.gsh : group_dynamic_data.lsh;
    for (i = 0; i < group_static_data.dim; i++)
    {
      if (dims[group_static_data.dim - i - 1] != (int) dims_stored[i])
      {
        flag = 1;
      }
    }
  }

  if (flag)
  {
    if (group_static_data.dim != (int) rank_stored)
    {
      CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Variable dimensions don't match for '%s', got %d, expected "
                  "%d", fullname, (int) rank_stored, group_static_data.dim);
    }
    else
    {
      msg = NULL;
      Util_asprintf (&msg, "Variable sizes don't match for '%s', got (%d",
                     fullname, (int) dims_stored[0]);
      for (i = 1; i < group_static_data.dim; i++)
      {
        oldmsg = msg;
        Util_asprintf (&msg, "%s, %d", msg, (int) dims_stored[i]);
        free (oldmsg);
      }
      dims = unchunked ? group_dynamic_data.gsh : group_dynamic_data.lsh;
      oldmsg = msg;
      Util_asprintf (&msg, "%s), expected (%d", msg,
                     dims[group_static_data.dim - 1]);
      free (oldmsg);
      for (i = 1; i < group_static_data.dim; i++)
      {
        oldmsg = msg;
        Util_asprintf (&msg, "%s, %d", msg, dims[group_static_data.dim-i-1]);
        free (oldmsg);
      }
      oldmsg = msg;
      Util_asprintf (&msg, "%s)", msg);
      free (oldmsg);
      CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING, "%s", msg);
      free (msg);
    }
    return (-1);
  }
  if (dims_stored)
  {
    free (dims_stored);
  }

  if (! CCTK_QueryGroupStorageI (GH, groupindex))
  {
    CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Can't read into variable '%s': no storage", fullname);
    return (-1);
  }

  /* close the first chunk if the dataset is a chunk group */
  if (is_group)
  {
    HDF5_ERROR (H5Dclose (dataset));
  }

  *grouptype = group_static_data.grouptype;

  return (0);
}


static int IOHDF5Util_RestoreGS (hid_t dataset, int vindex, int timelevel,
                                 recover_info_t *rec_info)
{
  void *data;
#ifdef CCTK_MPI
  int proc;
  CCTK_INT var_info[3];
#endif


  data = CCTK_VarDataPtrI (rec_info->it_info->GH, timelevel, vindex);

  /* read the data into the local variable ... */
  HDF5_ERROR (H5Dread (dataset, rec_info->hdf5type, H5S_ALL, H5S_ALL,
                       H5P_DEFAULT, data));

#ifdef CCTK_MPI
  /* ... and communicate it for the MPI case */

  /* set the variable's index and the timelevel */
  var_info[0] = vindex; var_info[1] = timelevel; var_info[2] = 1;

  /* send info and data to the non-IO processors */
  for (proc = rec_info->it_info->ioproc + 1;
       proc < rec_info->it_info->ioproc + rec_info->it_info->ioproc_every &&
       proc < CCTK_nProcs (rec_info->it_info->GH);
       proc++)
  {
    CACTUS_MPI_ERROR (MPI_Send (var_info, 3, PUGH_MPI_INT, proc, MPITAGBASE,
                            PUGH_pGH (rec_info->it_info->GH)->PUGH_COMM_WORLD));
    CACTUS_MPI_ERROR (MPI_Send (data, 1, rec_info->mpi_type, proc, MPITAGBASE,
                            PUGH_pGH (rec_info->it_info->GH)->PUGH_COMM_WORLD));
  }
#endif /* CCTK_MPI */

  return (0);
}


static int IOHDF5Util_RestoreGA (hid_t dataset, int vindex, int timelevel,
                                 recover_info_t *rec_info)
{
#ifdef CCTK_MPI
  int i, dim, proc, npoints;
  CCTK_INT var_info[3];
  pGH *pughGH;
  void *buffer, *data;
  hid_t filespace, memspace, chunk;
  pGExtras *extras;
  char chunkname[32];
#if (H5_VERS_MAJOR == 1 && \
     (H5_VERS_MINOR < 6 || (H5_VERS_MINOR == 6 && H5_VERS_RELEASE < 4)))
  hssize_t *chunk_origin;
#else
  hsize_t *chunk_origin;
#endif
  hsize_t *chunk_dims;
#endif


  /* single processor case is easy: just read the whole dataset */
  if (CCTK_nProcs (rec_info->it_info->GH) == 1)
  {
    HDF5_ERROR (H5Dread (dataset, rec_info->hdf5type, H5S_ALL, H5S_ALL,
                         H5P_DEFAULT,
                         rec_info->it_info->GH->data[vindex][timelevel]));
    return (0);
  }

#ifdef CCTK_MPI

  /* Get the handle for PUGH extensions */
  pughGH = PUGH_pGH (rec_info->it_info->GH);

  /* Get the pGExtras pointer as a shortcut */
  extras = ((pGA ***) pughGH->variables)[vindex][timelevel]->extras;

  /* get the dimension of the variable */
  dim = CCTK_GroupDimFromVarI (vindex);
  chunk_origin = malloc (dim * sizeof (*chunk_origin));
  chunk_dims = malloc (dim * sizeof (*chunk_dims));

  /* allocate memory for the biggest chunk */
  npoints = 1;
  for (proc = rec_info->it_info->ioproc + 1;
       proc < rec_info->it_info->ioproc + rec_info->it_info->ioproc_every &&
       proc < CCTK_nProcs (rec_info->it_info->GH);
       proc++)
  {
    if (npoints < extras->rnpoints[proc])
    {
      npoints = extras->rnpoints[proc];
    }
  }
  buffer = npoints > 0 ? malloc (npoints * rec_info->element_size) : NULL;

  /* set the variable's index and timelevel to restore */
  var_info[0] = vindex; var_info[1] = timelevel;

  /* now loop over the group of processors associated to each IO processor */
  for (proc = rec_info->it_info->ioproc;
       proc < rec_info->it_info->ioproc + rec_info->it_info->ioproc_every &&
       proc < CCTK_nProcs (rec_info->it_info->GH);
       proc++)
  {
    /* read own data directly into variable */
    if (proc == rec_info->it_info->ioproc)
    {
      data = rec_info->it_info->GH->data[vindex][timelevel];
    }
    else
    {
      data = buffer;
    }

    if (extras->rnpoints[proc] <= 0)
    {
      /* chunk contains no data - fall through */
    }
    else if (! rec_info->it_info->unchunked)
    {
      /* Chunked data is stored as separate chunk datasets within a group.
         So open, read and close the separate chunks here. */
      sprintf (chunkname, "chunk%d", proc - rec_info->it_info->ioproc);
      HDF5_ERROR (chunk = H5Dopen (dataset, chunkname));
      HDF5_ERROR (H5Dread (chunk, rec_info->hdf5type, H5S_ALL, H5S_ALL,
                           H5P_DEFAULT, data));
      HDF5_ERROR (H5Dclose (chunk));

    }
    else
    {
      /* Unchunked data is read as one hyperslab per processor.
         So prepare the memspace and the filespace and read the hyperslab. */
      for (i = 0; i < dim; i++)
      {
        chunk_dims[dim - 1 - i] = extras->rnsize[proc][i];
        chunk_origin[dim - 1 - i] = extras->lb[proc][i];
      }

      HDF5_ERROR (filespace = H5Dget_space (dataset));
      HDF5_ERROR (memspace = H5Screate_simple (dim, chunk_dims, NULL));
      HDF5_ERROR (H5Sselect_hyperslab (filespace, H5S_SELECT_SET,
                                       chunk_origin, NULL, chunk_dims, NULL));

      HDF5_ERROR (H5Dread (dataset, rec_info->hdf5type, memspace, filespace,
                           H5P_DEFAULT, data));

      HDF5_ERROR (H5Sclose (memspace));
      HDF5_ERROR (H5Sclose (filespace));
    }

    /* send the index and the data to the non-IO processors */
    if (proc != rec_info->it_info->ioproc)
    {
      var_info[2] = extras->rnpoints[proc];
      CACTUS_MPI_ERROR (MPI_Send (var_info, 3, PUGH_MPI_INT, proc,
                                  MPITAGBASE, pughGH->PUGH_COMM_WORLD));
      if (extras->rnpoints[proc] > 0)
      {
        CACTUS_MPI_ERROR (MPI_Send (data, extras->rnpoints[proc],
                                    rec_info->mpi_type, proc, MPITAGBASE,
                                    pughGH->PUGH_COMM_WORLD));
      }
    }
  }

  if (buffer)
  {
    free (buffer);
  }
  free (chunk_dims);
  free (chunk_origin);
#endif /* CCTK_MPI */

  return (0);
}


static herr_t processDataset (hid_t group, const char *datasetname, void *arg)
{
  const ioGH *ioUtilGH;
  const ioHDF5UtilGH *myGH;
  int vindex, vtype, gtype, timelevel, iteration, is_group, retval;
  iterate_info_t *it_info = arg;
  recover_info_t rec_info;
  hid_t dataset;
  H5G_stat_t object_info;
  char *fullname;
  DECLARE_CCTK_PARAMETERS


  /* skip the global attributes and GH extensions groups */
  if (! strcmp (datasetname, CACTUS_PARAMETERS_GROUP) ||
      ! strcmp (datasetname, GLOBAL_ATTRIBUTES_GROUP))
  {
    return (0);
  }

  retval = 0;

  /* decompose the datasetname, ignore the iteration number */
  fullname = malloc (strlen (datasetname));
  if (sscanf (datasetname, "%[^ ] timelevel %d at iteration %d",
              fullname, &timelevel, &iteration) != 3)
  {
    CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot parse datasetname '%s'", datasetname);
    retval = -1;
  }

  /* check if there is a matching variable */
  vindex = CCTK_VarIndex (fullname);
  if (vindex < 0)
  {
    CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                "No matching variable found for '%s'", fullname);
    retval = -1;
  }

  if (retval < 0)
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Ignoring dataset '%s'", datasetname);
    free (fullname);
    return (0);
  }

  ioUtilGH = CCTK_GHExtension (it_info->GH, "IO");
  myGH = CCTK_GHExtension (it_info->GH, "IOHDF5Util");

  /* if we read in initial data via the file reader interface
     check whether the user wants to have this variable with a specific
     iteration number restored */
  if (ioUtilGH->do_inVars && ioUtilGH->do_inVars[vindex] >= 0 &&
      ioUtilGH->do_inVars[vindex] != iteration + 1)
  {
    if (CCTK_Equals (verbose, "full"))
    {
      CCTK_VInfo (CCTK_THORNSTRING, "Ignoring dataset '%s' for file reader "
                  "recovery", datasetname);
    }
    free (fullname);
    return (0);
  }

  HDF5_ERROR (H5Gget_objinfo (group, datasetname, 0, &object_info));
  is_group = object_info.type == H5G_GROUP;
  if (is_group)
  {
    HDF5_ERROR (dataset = H5Gopen (group, datasetname));
  }
  else
  {
    HDF5_ERROR (dataset = H5Dopen (group, datasetname));
  }

  /* read in the dataset's attributes and verify them */
  if (GetCommonAttributes (it_info->GH, dataset, vindex, fullname,
                           it_info->unchunked, &gtype, &timelevel,
                           is_group, it_info->has_version) < 0)
  { 
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Ignoring dataset '%s'", datasetname);
  }
  else
  {
    if (CCTK_Equals (verbose, "full") || 
        (ioUtilGH->do_inVars && ! CCTK_Equals (verbose, "none")))
    {
      CCTK_VInfo (CCTK_THORNSTRING, "Restoring variable '%s' (timelevel %d, "
                  "cctk_iteration %d)", fullname, timelevel, iteration);
    }

    vtype = CCTK_VarTypeI (vindex);
    rec_info.it_info = it_info;
    rec_info.element_size = CCTK_VarTypeSize (vtype);
#ifdef CCTK_MPI
    rec_info.mpi_type = PUGH_MPIDataType (PUGH_pGH (it_info->GH), vtype);
#endif
    rec_info.hdf5type = IOHDF5Util_DataType (myGH, vtype);
    if (rec_info.element_size <= 0 ||
#ifdef CCTK_MPI
        rec_info.mpi_type == MPI_DATATYPE_NULL ||
#endif
        rec_info.hdf5type < 0)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Unsupported variable datatype %d", vtype);
    }
    else if (gtype != CCTK_SCALAR && gtype != CCTK_GF && gtype != CCTK_ARRAY)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Unsupported group type %d", gtype);
    }
    else
    {
      /* Read in the data */
      if (gtype == CCTK_SCALAR)
      {
        IOHDF5Util_RestoreGS (dataset, vindex, timelevel, &rec_info);
      }
      else
      {
        IOHDF5Util_RestoreGA (dataset, vindex, timelevel, &rec_info);
      }

      /* increment counter for total number of recovered variables */
      it_info->num_recovered++;
    }
  }

  if (is_group)
  {
    HDF5_ERROR (H5Gclose (dataset));
  }
  else
  {
    HDF5_ERROR (H5Dclose (dataset));
  }

  free (fullname);

  return (0);
}
