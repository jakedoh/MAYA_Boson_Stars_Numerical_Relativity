/* $Id: Init.c,v 1.4 2005/01/13 17:21:13 herrmann Exp $ */
/* init stuff */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void Init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (verbose>2)
    CCTK_INFO("init not_happy_mass");

  *not_happy_mass=1;

  CCTK_INFO("setting up initial range for binary search");
  /* FIXME : arbitrary */
  *par_m_min=m_range_min_init;
  *par_m_max=m_range_max_init;
}
