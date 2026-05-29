/*@@
   @file      Write.c
   @date      Thu 18 April 2002
   @author    Thomas Radke
   @desc
              Function which implements output method.
   @enddesc
   @version   $Id$
 @@*/


#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "ioSampleIOMethodGH.h"

/* the rcs ID and its dummy function to use it (this is so that you can grep
   for a version number of this source file in a Cactus executable) */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusExamples_SampleIO_Write_c)


/*@@
   @routine   SampleIO_Write
   @date      Thu 18 April 2002
   @author    Thomas Radke
   @desc
              Gets the value of a grid function/array at a chosen location
              as a hyperslab and prints it to screen
   @enddesc
   @calls     Hyperslab_GlobalMappingByIndex
              Hyperslab_FreeMapping
              Hyperslab_Get

   @var       GH
   @vdesc     pointer to CCTK GH
   @vtype     const cGH *
   @vio       in
   @endvar
   @var       vindex
   @vdesc     index of variable to output
   @vtype     int
   @vio       in
   @endvar

   @returntype int
   @returndesc
                0 for success, or<BR>
               -1 if variable has no storage assigned<BR>
               -2 if hyperslab mapping couldn't be defined<BR>
               -3 if hyperslab point couldn't be extracted
   @endreturndesc
@@*/
int SampleIO_Write (const cGH *GH, int vindex)
{
  int retval;
  char *fullname;
  CCTK_INT mapping;
  CCTK_INT origin[3], extent[1], direction[1*3];
  CCTK_REAL value;
  DECLARE_CCTK_PARAMETERS;


  retval = 0;

  /* get the full variable name (used for error messages) */
  fullname = CCTK_FullName (vindex);

  /* check if variable has storage assigned */
  if (! CCTK_QueryGroupStorageI (GH, CCTK_GroupIndexFromVarI (vindex)))
  {
    CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                "No SampleIO output for '%s' (no storage)", fullname);
    free (fullname);
    return (-1);
  }

  /*
   * set the origin from SampleIO's parameters,
   * the extent being 1 (for a single point),
   * and the direction as a diagonal (spanning in each direction of the 3D grid)
   */
  origin[0] = point_x; origin[1] = point_y; origin[2] = point_z;
  extent[0] = 1;
  direction[0] = 1; direction[1] = direction[2] = 0;

  /* define the point to output at as a 1D hyperslab */
  mapping = Hyperslab_GlobalMappingByIndex (GH, vindex, 1, direction,
                                            origin, extent, NULL, -1,
                                            NULL, NULL);
  if (mapping >= 0)
  {
    /*
     * get the hyperslab data as a REAL value on all processors
     * and print it to screen
     */
    if (! Hyperslab_Get (GH, mapping, -1, vindex, 0, CCTK_VARIABLE_REAL,
                         (CCTK_POINTER) &value))
    {
      CCTK_VInfo (CCTK_THORNSTRING, "'%s' at [x][y][z] = [%d][%d][%d]: %f",
                  fullname,
                  (int) point_x, (int) point_y, (int) point_z, (double) value);
    }
    else
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Failed to extract hyperslab for variable '%s'", fullname);
      retval = -3;
    }

    /* release the mapping structure */
    Hyperslab_FreeMapping (mapping);
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Failed to define hyperslab mapping for variable '%s'",
                fullname);
    retval = -2;
  }

  free (fullname);

  return (retval);
}
