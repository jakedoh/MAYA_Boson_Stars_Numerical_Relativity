#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include <stdio.h>
#include <stdlib.h>
void InitCoord(CCTK_ARGUMENTS);
voide CCTK_FCALL CCTK_FNAME(getsetl)(int *, int *);

void InitCoord1(CCTK_ARGUMENTS)
{
	CCTK_INT i, j;
	CCTK_REAL dx, dy;
	int *px, *py;
	CCTK_INT lnx, lny;

	DECLARE_CCTK_ARGUMENTS
	DECLARE_CCTK_PARAMETERS

	if ((px = CCTK_Array(cctk_GH, 1, "my2dc"))==NULL)
		call CCTK_Warn("Could not get size");
	if ((py = CCTK_Array(cctk_GH, 2, "my2dc"))==NULL)
		call CCTK_Warn("Could not get size");

	getsetl(px, py);	
}
