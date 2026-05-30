 /*@@
   @file      Startup.c
   @date      Fri Oct 6 2000
   @author    Thomas Radke
   @desc
              Startup and termination routines for IOHDF5Util.
   @enddesc
   @version   $Id$
@@*/


#include <stdlib.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "ioHDF5UtilGH.h"


/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusPUGHIO_IOHDF5Util_Startup_c)


/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/
int IOHDF5Util_Startup (void);
void IOHDF5Util_Terminate (cGH *GH);


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
static void *IOHDF5Util_SetupGH (tFleshConfig *config, int conv_level, cGH *GH);


 /*@@
   @routine    IOHDF5Util_Startup
   @date       Fri Oct 6 2000
   @author     Thomas Radke
   @desc
               The startup registration routine for IOHDF5Util.
               Registers the GH extensions needed for IOHDF5Util
               and the routine to set it up.
   @enddesc

   @calls      CCTK_RegisterGHExtension
               CCTK_RegisterGHExtensionSetupGH
@@*/
int IOHDF5Util_Startup (void)
{
  /* check dynamic thorn dependencies */
  if (CCTK_GHExtensionHandle ("IO") < 0)
  {
    CCTK_WARN (1, "Thorn IOUtil was not activated. "
                  "No HDF5 IO methods will be registered.");
  }
  else if (CCTK_GHExtensionHandle ("PUGH") < 0)
  {
    CCTK_WARN (1, "Thorn PUGH was not activated. "
                  "No HDF5 IO methods will be registered.");
  }
  else
  {
    CCTK_RegisterGHExtensionSetupGH (CCTK_RegisterGHExtension ("IOHDF5Util"),
                                     IOHDF5Util_SetupGH);
  }
  return 0;
}


 /*@@
   @routine    IOHDF5Util_Terminate
   @date       Fri Oct 6 2000
   @author     Thomas Radke
   @desc
               The termination registration routine for IOHDF5Util.
               It frees any resources kept on the GH extensions.
   @enddesc

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
@@*/
void IOHDF5Util_Terminate (cGH *GH)
{
  ioHDF5UtilGH *myGH;


  myGH = (ioHDF5UtilGH *) CCTK_GHExtension (GH, "IOHDF5Util");
  if (myGH)
  {
    /* close the dataspaces used to write scalar and array attributes */
    if (myGH->scalar_dataspace >= 0)
    {
      HDF5_ERROR (H5Sclose (myGH->scalar_dataspace));
    }
    if (myGH->array_dataspace >= 0)
    {
      HDF5_ERROR (H5Sclose (myGH->array_dataspace));
    }

    /* close the predefined complex and string datatypes */
    if (myGH->HDF5_COMPLEX >= 0)
    {
      HDF5_ERROR (H5Tclose (myGH->HDF5_COMPLEX));
    }
#ifdef HAVE_CCTK_REAL4
    if (myGH->HDF5_COMPLEX8 >= 0)
    {
      HDF5_ERROR (H5Tclose (myGH->HDF5_COMPLEX8));
    }
#endif
#ifdef HAVE_CCTK_REAL8
    if (myGH->HDF5_COMPLEX16 >= 0)
    {
      HDF5_ERROR (H5Tclose (myGH->HDF5_COMPLEX16));
    }
#endif
#ifdef HAVE_CCTK_REAL16
    if (myGH->HDF5_COMPLEX32 >= 0)
    {
      HDF5_ERROR (H5Tclose (myGH->HDF5_COMPLEX32));
    }
#endif
    if (myGH->HDF5_STRING >= 0)
    {
      HDF5_ERROR (H5Tclose (myGH->HDF5_STRING));
    }
  }
}


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
 /*@@
   @routine    IOHDF5Util_SetupGH
   @date       Fri Oct 6 2000
   @author     Thomas Radke
   @desc
               Allocates and sets up IOHDF5Util's GH extension structure.
   @enddesc

   @calls      CCTK_RegisterIOMethod
               CCTK_RegisterIOMethodOutputGH
               CCTK_RegisterIOMethodOutputVarAs
               CCTK_RegisterIOMethodTimeToOutput
               CCTK_RegisterIOMethodTriggerOutput

   @var        config
   @vdesc      the CCTK configuration as provided by the flesh
   @vtype      tFleshConfig *
   @vio        unused
   @endvar
   @var        conv_level
   @vdesc      the convergence level
   @vtype      int
   @vio        unused
   @endvar
   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      cGH *
   @vio        unused
   @endvar

   @returntype void *
   @returndesc
               pointer to the allocated GH extension structure
   @endreturndesc
@@*/
static void *IOHDF5Util_SetupGH (tFleshConfig *config, int conv_level, cGH *GH)
{
  ioHDF5UtilGH *myGH;
  DECLARE_CCTK_PARAMETERS


  /* suppress compiler warnings about unused variables */
  (void) (config + 0);
  (void) (conv_level + 0);
  (void) (GH + 0);

  /* allocate the GH extension */
  myGH = (ioHDF5UtilGH *) malloc (sizeof (ioHDF5UtilGH));

  /* predefine dataspaces for writing scalar and array attributes */
  /* the dimension of the array dataspace is set when used */
  HDF5_ERROR (myGH->scalar_dataspace = H5Screate (H5S_SCALAR));
  HDF5_ERROR (myGH->array_dataspace  = H5Screate (H5S_SIMPLE));

  /* predefine HDF5_COMPLEX datatypes */
  HDF5_ERROR (myGH->HDF5_COMPLEX =
              H5Tcreate (H5T_COMPOUND, sizeof (CCTK_COMPLEX)));
  HDF5_ERROR (H5Tinsert (myGH->HDF5_COMPLEX, "real",
                         0, HDF5_REAL));
  HDF5_ERROR (H5Tinsert (myGH->HDF5_COMPLEX, "imag",
                         sizeof (CCTK_REAL), HDF5_REAL));
#ifdef HAVE_CCTK_REAL4
  HDF5_ERROR (myGH->HDF5_COMPLEX8 =
              H5Tcreate (H5T_COMPOUND, sizeof (CCTK_COMPLEX8)));
  HDF5_ERROR (H5Tinsert (myGH->HDF5_COMPLEX8, "real",
                         0, HDF5_REAL4));
  HDF5_ERROR (H5Tinsert (myGH->HDF5_COMPLEX8, "imag",
                         sizeof (CCTK_REAL4), HDF5_REAL4));
#endif
#ifdef HAVE_CCTK_REAL8
  HDF5_ERROR (myGH->HDF5_COMPLEX16 =
              H5Tcreate (H5T_COMPOUND, sizeof (CCTK_COMPLEX16)));
  HDF5_ERROR (H5Tinsert (myGH->HDF5_COMPLEX16, "real",
                         0, HDF5_REAL8));
  HDF5_ERROR (H5Tinsert (myGH->HDF5_COMPLEX16, "imag",
                         sizeof (CCTK_REAL8), HDF5_REAL8));
#endif
#ifdef HAVE_CCTK_REAL16
  HDF5_ERROR (myGH->HDF5_COMPLEX32 =
              H5Tcreate (H5T_COMPOUND, sizeof (CCTK_COMPLEX32)));
  HDF5_ERROR (H5Tinsert (myGH->HDF5_COMPLEX32, "real",
                         0, HDF5_REAL16));
  HDF5_ERROR (H5Tinsert (myGH->HDF5_COMPLEX32, "imag",
                         sizeof (CCTK_REAL16), HDF5_REAL16));
#endif

  /* predefine a C string datatype */
  HDF5_ERROR (myGH->HDF5_STRING = H5Tcopy (H5T_C_S1));

  return (myGH);
}
