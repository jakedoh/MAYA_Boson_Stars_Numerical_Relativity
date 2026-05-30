/*@@
  @file      ChooseOutput.c
  @date      Thu 18 April 2002
  @author    Thomas Radke
  @desc
             Choose what 2D planes to output by IOJpeg.
  @enddesc

  @version   $Id$
 @@*/

#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "CactusBase/IOUtil/src/ioutil_Utils.h"
#include "ioJpegGH.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusIO_IOJpeg_ChooseOutput_c)


/********************************************************************
 ********************    Macro Definitions   ************************
 ********************************************************************/
/* macro to choose origin according actual parameter settings:
     1. Indices from IOJpeg
     2. Indices from IOUtil
     3. Coords from IOJpeg
     4. Coords from IOUtil
 */
#define GET_SLICE(IOJpeg_param, IOUtil_param, index, coord)                   \
        {                                                                     \
          index = IOJpeg_param##i >= 0 ? IOJpeg_param##i : IOUtil_param##i;   \
          coord = IOJpeg_param != -424242 ? IOJpeg_param : IOUtil_param;      \
        }


/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/

/*@@
   @routine   IOJpeg_ChooseOutput
   @date      Thu 18 April 2002
   @author    Thomas Radke
   @desc
              Use parameters to choose the 2D slices through the output data.
   @enddesc

   @calls     IOUtil_2DPlanes

   @var       GH
   @vdesc     Pointer to CCTK grid hierarchy
   @vtype     const cGH *
   @vio       in
   @endvar
 @@*/
void IOJpeg_ChooseOutput (CCTK_ARGUMENTS)
{
  int i, maxdim;
  ioJpegGH *myGH;
  int origin_index[3];
  CCTK_REAL origin_phys[3];
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;


  GET_SLICE (out2D_xyplane_z, out_xyplane_z, origin_index[0], origin_phys[0]);
  GET_SLICE (out2D_xzplane_y, out_xzplane_y, origin_index[1], origin_phys[1]);
  GET_SLICE (out2D_yzplane_x, out_yzplane_x, origin_index[2], origin_phys[2]);

  maxdim = CCTK_MaxDim ();
  myGH = (ioJpegGH *) CCTK_GHExtension (cctkGH, "IOJpeg");
  myGH->sp2xyz = (int **) malloc (3 * sizeof (int *));

  for (i = 0; i < maxdim; i++)
  {
    myGH->sp2xyz[i] = (int *) calloc (i + 1, sizeof (int));

    if (i > 0 && i < 3)
    {
      IOUtil_2DPlanes (cctkGH, i + 1, origin_index, origin_phys, myGH->sp2xyz[i]);
    }
  }
}
