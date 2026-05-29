/* startup routine for Exact thorn */
/* $Header$ */

#include "cctk.h"
#include "Slicing.h"

/*
 * prototypes for scheduled routines
 */
int Exact__RegisterSlicing(void);

/******************************************************************************/

int Exact__RegisterSlicing(void) 
{
  int handle;
  handle=Einstein_RegisterSlicing("exact");
  if (handle<0) CCTK_WARN(1,"Cannot register exact slicing");
  return 0;
}
