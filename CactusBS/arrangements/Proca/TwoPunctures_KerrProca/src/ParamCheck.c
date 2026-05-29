/* TwoPunctures_KerrProca:  File  "TwoPunctures_KerrProca.c"*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "TP_utilities.h"
#include "TwoPunctures.h"

/* -------------------------------------------------------------------*/
void
TPKP_TwoPunctures_ParamCheck (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (use_sources)
  {
    if (! CCTK_IsFunctionAliased ("Set_Rho_ADM"))
      CCTK_WARN (0, "Matter sources have been enabled, "
                 "but there is no aliased function for matter sources.");
  }


  if( par_S_plus[1] != 0 || par_S_plus[2] != 0 || par_S_minus[0] != 0 || par_S_minus[1] != 0 || par_S_minus[2] != 0 )
    { 
        CCTK_WARN( 0, "Initialized wrong spin components. TwoPunctures_KerrProca expects spin in x direction, i.e., only par_S_plus[0] /= 0 ");
    }


}
