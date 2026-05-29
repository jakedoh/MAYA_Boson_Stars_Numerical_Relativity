#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <stdio.h>

void NullGrid_RegisterRange(CCTK_ARGUMENTS)
{
	CCTK_REAL lowerq, lowerp, upperq, upperp;
	CCTK_INT loweri, upperi;
	char message[1000];
	int gsh[2];
	DECLARE_CCTK_ARGUMENTS;
	DECLARE_CCTK_PARAMETERS;
      
	if ( CCTK_GroupgshGN(cctkGH, 2, gsh, "NullGrid::StCrd") !=0)
		CCTK_WARN(0,"Error obtaining gsh");

	lowerq = - *qsize;
	upperq = + *qsize;
	lowerp = - *qsize;
	upperp = + *qsize;

	/*  q registration */
	if (CCTK_CoordRegisterRange(cctkGH,lowerq,upperq,-1,"q","stereo") !=0)
		CCTK_WARN(0,"Error registering q range");

	loweri = 0;   /* physical range is the entire array */
	upperi = gsh[0] - 1;
	if (CCTK_CoordRegisterRangePhysIndex(cctkGH,loweri,upperi,-1,"q","stereo")!=0)
		CCTK_WARN(0,"Error registering q phys index");


	/*  p registration */
	if (CCTK_CoordRegisterRange(cctkGH,lowerp,upperp,-1,"p","stereo") !=0)
		CCTK_WARN(0,"Error registering p range");

	loweri = 0;
	upperi = gsh[1] - 1;
	if (CCTK_CoordRegisterRangePhysIndex(cctkGH,loweri,upperi,-1,"p","stereo")!=0)
		CCTK_WARN(0,"Error registering p phys index");

	CCTK_INFO("Stereographic grid steps:");
	sprintf(message, "dq => %lf , dp => %lf", null_delta[0], null_delta[1]);
	CCTK_INFO(message);

	CCTK_INFO("Stereographic grid ranges:");

	sprintf(message, "q => [%lf,%lf], p => [%lf,%lf]", lowerq, upperq, lowerp, upperp);
	CCTK_INFO(message);

	CCTK_INFO("Characteristic radial grid step:");
	sprintf(message, "dx => %lf", *null_dx);
	CCTK_INFO(message);

	CCTK_INFO("Characteristic radial grid ranges:");
	sprintf(message, "xb => [%lf,%lf], rb => [%lf,+infty]", null_xb[0], null_xb[N_radial_pts-1], null_rb[0]);
	CCTK_INFO(message);

}
