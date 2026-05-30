/*
  Set the symmetries for the GiRaFFE variables
*/

#include "cctk.h"
#include <cstdio>
#include <cstdlib>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"
#include "GiRaFFE_headers.h"

extern "C" void GiRaFFE_InitSymBound(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // A-fields CAN EVOLVE, yet the magnetic fields may not change! Thus setting EM_BC==frozen may result in time-changing B-fields at the boundaries (usually glitchy behavior).
  // if( CCTK_EQUALS(Velocity_BC,"frozen") && !CCTK_EQUALS(EM_BC,"frozen") ||
  //     !CCTK_EQUALS(Velocity_BC,"frozen") && CCTK_EQUALS(EM_BC,"frozen") )
  //   CCTK_VError(VERR_DEF_PARAMS,"If Velocity_BC or EM_BC is set to FROZEN, BOTH must be set to frozen!");
  if(CCTK_EQUALS(EM_BC,"frozen")) CCTK_VInfo(CCTK_THORNSTRING,"Warning: EM_BC=frozen WILL NOT NECESSARILY result in frozen B-fields at the boundaries, since A-fields may evolve even with B fixed in time!\n");

  if ((cctk_nghostzones[0]<3 || cctk_nghostzones[1]<3 || cctk_nghostzones[2]<3))
    CCTK_VError(VERR_DEF_PARAMS,"ERROR: The version of PPM in this thorn requires 3 ghostzones. You only have (%d,%d,%d) ghostzones!",cctk_nghostzones[0],cctk_nghostzones[1],cctk_nghostzones[2]);

  if(cctk_iteration==0) {
    CCTK_VInfo(CCTK_THORNSTRING,"Setting Symmetry = %s... at iteration = %d",Symmetry,cctk_iteration);

    int sym[3];

    if(CCTK_EQUALS(Symmetry,"none")) {
      /* FIRST SET NO SYMMETRY OPTION */
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      SetCartSymGN(cctkGH,sym,"GiRaFFE::grmhd_conservatives");
      SetCartSymGN(cctkGH,sym,"GiRaFFE::em_Ax");
      SetCartSymGN(cctkGH,sym,"GiRaFFE::em_Ay");
      SetCartSymGN(cctkGH,sym,"GiRaFFE::em_Az");
      SetCartSymGN(cctkGH,sym,"GiRaFFE::em_psi6phi");
      SetCartSymGN(cctkGH,sym,"GiRaFFE::grmhd_primitives_allbutBi");
    } else if(CCTK_EQUALS(Symmetry,"equatorial")) {
      /* THEN SET EQUATORIAL SYMMETRY OPTION */
      // Set default to no symmetry, which is correct for scalars and most vectors:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      SetCartSymGN(cctkGH,sym,"GiRaFFE::grmhd_conservatives");
      // Don't worry about the wrong sym values since A_{\mu} is staggered
      // and we're going to impose the symmetry separately
      SetCartSymGN(cctkGH,sym,"GiRaFFE::em_Ax");
      SetCartSymGN(cctkGH,sym,"GiRaFFE::em_Ay");
      SetCartSymGN(cctkGH,sym,"GiRaFFE::em_Az");
      SetCartSymGN(cctkGH,sym,"GiRaFFE::em_psi6phi");

      SetCartSymGN(cctkGH,sym,"GiRaFFE::grmhd_primitives_allbutBi");

      // Then set unstaggered B field variables
      sym[2] = -Sym_Bz;
      SetCartSymVN(cctkGH, sym,"GiRaFFE::Bx");
      SetCartSymVN(cctkGH, sym,"GiRaFFE::By");
      sym[2] = Sym_Bz;
      SetCartSymVN(cctkGH, sym,"GiRaFFE::Bz");

      sym[2] = -1;
      SetCartSymVN(cctkGH, sym,"GiRaFFE::mhd_st_z");
      SetCartSymVN(cctkGH, sym,"GiRaFFE::vz");
    } else {
      CCTK_VError(VERR_DEF_PARAMS,"GiRaFFE_initsymbound: Should not be here; picked an impossible symmetry.");
    }
  }
}

