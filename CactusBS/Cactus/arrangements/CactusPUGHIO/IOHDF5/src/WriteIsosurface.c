/*@@
   @file      WriteIsosurface.c
   @date      Sat Mar 25 2000
   @author    Thomas Radke
   @desc 
              Output routine for isosurfaces in HDF5 fiber bundle format.
   @enddesc 
   @version   $Id$
 @@*/


#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "ioHDF5GH.h"

/* the rcs ID and its dummy funtion to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusPUGHIO_IOHDF5_WriteIsosurface_c)


#define SOCKET_URL_PREFIX      "sock://"
#define SOCKET_URL_PREFIX_LEN  (sizeof (SOCKET_URL_PREFIX) - 1)

/* this struct stores the current slice info for further queries */
typedef struct {
  CCTK_INT iteration;
  CCTK_INT timelevel;
} sliceInfo_t;


/*@@
  @routine  IOHDF5_WriteIsosurface
  @author   Thomas Radke
  @date     Sat Mar 25 2000
  @calledby thorn CactusPUGHIO/IsoSurfacer
  @var     
  @vdesc  
  @vtype 
  @vio  
  @vcomment
  @endvar
@@*/
int IOHDF5_WriteIsosurface (cGH *GH,
                            const char *filename,
                            const char *GVname,
                            CCTK_INT iteration,
                            CCTK_INT timelevel,
                            CCTK_REAL isoval,
                            CCTK_REAL minval,
                            CCTK_REAL maxval,
                            int nTriangles,
                            CCTK_INT4 *triangles,
                            int nVertices,
                            CCTK_REAL4 *vertices)
{
  hid_t h5file, plist;
  hid_t aspace, dataspace, dataset;
  hid_t group0, group1, group2, group3, group4;
  hsize_t dim [2];
  H5FD_stream_fapl_t fapl;
  char sliceName [40], *isoVarname;
  int isSocketURL;
  sliceInfo_t *sliceInfo;
  hobj_ref_t reference;
  static pNamedData *isoFilesList = NULL;
  static int first_called = 1;
  static hid_t hdf5_string;

  /* FIXME: should use the predefined type from iohdf5GH
            but we cannot be sure that thorn IOHDF5 was activated */
  if (first_called) {
    first_called = 0;
    /* predefine a C string datatype */
    CACTUS_IOHDF5_ERROR (hdf5_string = H5Tcopy (H5T_C_S1));
  }

  CACTUS_IOHDF5_ERROR (plist = H5Pcreate (H5P_FILE_ACCESS));

  isSocketURL = strncasecmp (filename, SOCKET_URL_PREFIX, SOCKET_URL_PREFIX_LEN) == 0;
  if (isSocketURL)
  {
#ifndef H5_HAVE_STREAM
    CCTK_WARN (1, "The configured HDF5 installation doesn't include the Stream "
                  "VFD. No HDF5 streaming output available !");
    CACTUS_IOHDF5_ERROR (H5Pclose (plist));
    return (-1);
#else
    /* set VFD to stream and open the file */
    fapl.increment = 0;
    fapl.socket = -1;
    fapl.do_socket_io = 1;
    fapl.backlog = 5;
    fapl.broadcast_fn = NULL;
    fapl.broadcast_arg = NULL;

    CACTUS_IOHDF5_ERROR (H5Pset_fapl_stream (plist, &fapl));
    sliceInfo = NULL;
#endif
  }
  else
  {
    sliceInfo = (sliceInfo_t *) GetNamedData (isoFilesList, filename);
  }

  /* create/open the file */
  if (! sliceInfo)
    CACTUS_IOHDF5_ERROR (h5file = H5Fcreate (filename, H5F_ACC_TRUNC,
                                             H5P_DEFAULT, plist));
  else
    CACTUS_IOHDF5_ERROR (h5file = H5Fopen (filename, H5F_ACC_RDWR, plist));
  if (h5file < 0) {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Could not create/open output file '%s'", filename);
    return (-1);
  }

  /* create a dataspace for storing scalar attributes */
  CACTUS_IOHDF5_ERROR (aspace = H5Screate (H5S_SCALAR));

  /* open the highest-level hierarchy which is group "Bundle" */
  if (! sliceInfo)
    CACTUS_IOHDF5_ERROR (group0 = H5Gcreate (h5file, "Bundle", 0));
  else
    CACTUS_IOHDF5_ERROR (group0 = H5Gopen (h5file, "Bundle"));

  /* next level in hierarchy is the "Slice [<iteration>] [<timelevel>]" group */
  sprintf (sliceName, "Slice [%d] [%d]", (int) iteration, (int) timelevel);
  if (! sliceInfo ||
      sliceInfo->iteration != iteration || sliceInfo->timelevel != timelevel) {
    if (sliceInfo) {
      sliceInfo->iteration = iteration;
      sliceInfo->timelevel = timelevel;
    }
    CACTUS_IOHDF5_ERROR (group1 = H5Gcreate (group0, sliceName, 0));
    WRITE_ATTRIBUTE ("time", &GH->cctk_time, group1, aspace, 0,IOHDF5_REAL);
    WRITE_ATTRIBUTE ("timestep", &iteration, group1, aspace, 0, IOHDF5_INT);
    WRITE_ATTRIBUTE ("timelevel", &timelevel, group1, aspace, 0,IOHDF5_INT);
  } else
    CACTUS_IOHDF5_ERROR (group1 = H5Gopen (group0, sliceName));

  /* third level is a group for each grid variable with specific isovalue */
  isoVarname = (char *) malloc (strlen (GVname) + 40);
  sprintf (isoVarname, "Isolevel [%s] [%.10f]", GVname, (double) isoval);
  CACTUS_IOHDF5_ERROR (group2 = H5Gcreate (group1, isoVarname, 0));
  free (isoVarname);

  /* fourth level is group "PointTopology" */
  CACTUS_IOHDF5_ERROR (group3 = H5Gcreate (group2, "PointTopology", 0));

  /* fifth level is group "Representation" */
  CACTUS_IOHDF5_ERROR (group4 = H5Gcreate (group3, "Representation [CartesianChart]", 0));

#if 0
  if (IsoL.Ncolorinfo > 0) {
    dim [0] = (hsize_t) IsoL.Ncolorinfo;
    CACTUS_IOHDF5_ERROR (dataspace = H5Screate_simple (1, dim, NULL));
    CACTUS_IOHDF5_ERROR (dataset = H5Dcreate (h5file, "ColorMap",
                         IOHDF5_REAL4, dataspace, H5P_DEFAULT));
    CACTUS_IOHDF5_ERROR (H5Dwrite (dataset, IOHDF5_REAL4, H5S_ALL,
                         H5S_ALL, H5P_DEFAULT, &IsoL.colorinfo));
    CACTUS_IOHDF5_ERROR (H5Dclose (dataset));
    CACTUS_IOHDF5_ERROR (H5Sclose (dataspace));
  }
#endif

  /* write vertices */
  dim [0] = (hsize_t) nVertices;
  dim [1] = 3;
  CACTUS_IOHDF5_ERROR (dataspace = H5Screate_simple (2, dim, NULL));
  CACTUS_IOHDF5_ERROR (dataset = H5Dcreate (group4, "Positions",
                       IOHDF5_REAL4, dataspace, H5P_DEFAULT));
  CACTUS_IOHDF5_ERROR (H5Dwrite (dataset, IOHDF5_REAL4, H5S_ALL,
                       H5S_ALL, H5P_DEFAULT, vertices));
  WRITE_ATTRIBUTE ("data_min", &minval, dataset, aspace, 0, IOHDF5_REAL);
  WRITE_ATTRIBUTE ("data_max", &maxval, dataset, aspace, 0, IOHDF5_REAL);
#if 0
  WRITE_ATTRIBUTE ("instance", &isoindex, isogroup, aspace, 0, IOHDF5_INT);
#endif
  CACTUS_IOHDF5_ERROR (H5Rcreate (&reference, group4, "Positions", H5R_OBJECT, -1));
  CACTUS_IOHDF5_ERROR (H5Dclose (dataset));
  CACTUS_IOHDF5_ERROR (H5Sclose (dataspace));
  WRITE_ATTRIBUTE ("Domain", "CartesianChart", group4, aspace, 0, hdf5_string);
  CACTUS_IOHDF5_ERROR (H5Gclose (group4));
  CACTUS_IOHDF5_ERROR (H5Gclose (group3));

  /* fourth level is group "Connectivity" */
  CACTUS_IOHDF5_ERROR (group3 = H5Gcreate (group2, "Connectivity", 0));

  /* fifth level is group "Representation" */
  CACTUS_IOHDF5_ERROR (group4 = H5Gcreate (group3, "Representation [PointTopology]", 0));

  /* write polygons */
  dim [0] = (hsize_t) nTriangles;
  dim [1] = 3;
  CACTUS_IOHDF5_ERROR (dataspace = H5Screate_simple (2, dim, NULL));
  CACTUS_IOHDF5_ERROR (dataset = H5Dcreate (group4, "Positions",
                       IOHDF5_INT4, dataspace, H5P_DEFAULT));
  CACTUS_IOHDF5_ERROR (H5Dwrite (dataset, IOHDF5_INT4, H5S_ALL,
                       H5S_ALL, H5P_DEFAULT, triangles));
  CACTUS_IOHDF5_ERROR (H5Dclose (dataset));
  CACTUS_IOHDF5_ERROR (H5Sclose (dataspace));
  WRITE_ATTRIBUTE ("Domain", &reference, group4, aspace, 0, H5T_STD_REF_OBJ);
  CACTUS_IOHDF5_ERROR (H5Glink(group4, H5G_LINK_SOFT,
     "../../PointTopology/Representation [CartesianChart]/Positions", "Domain"));
  CACTUS_IOHDF5_ERROR (H5Sclose (aspace));
  CACTUS_IOHDF5_ERROR (H5Gclose (group4));
  CACTUS_IOHDF5_ERROR (H5Gclose (group3));
  CACTUS_IOHDF5_ERROR (H5Gclose (group2));
  CACTUS_IOHDF5_ERROR (H5Gclose (group1));
  CACTUS_IOHDF5_ERROR (H5Gclose (group0));

  CACTUS_IOHDF5_ERROR (H5Fclose (h5file));

  /* if this was a newly created file add it to the list of open files */
  if (! isSocketURL && ! sliceInfo) {
    sliceInfo_t *newInfo = (sliceInfo_t *) malloc (sizeof (sliceInfo_t));

    newInfo->iteration = iteration;
    newInfo->timelevel = timelevel;
    StoreNamedData (&isoFilesList, filename, newInfo);
  }

  return (0);
}
