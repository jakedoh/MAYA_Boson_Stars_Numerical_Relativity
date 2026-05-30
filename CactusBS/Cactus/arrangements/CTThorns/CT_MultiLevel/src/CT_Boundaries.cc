#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "cctk_Loop.h"

#include "carpet.hh"
#include "Carpet/Carpet/src/modes.hh"

#include "CT_MultiLevel.hh"

extern "C" void CT_Boundaries(CCTK_ARGUMENTS, const char *varname)
{
  BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
    BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      DECLARE_CCTK_PARAMETERS;

      if (CT_ProcessOwnsData())
      {
        if (CCTK_Equals(verbose, "yes")) CCTK_Info(CCTK_THORNSTRING, "Apply BCs...");

        int ivari = CCTK_FirstVarIndex(varname);
        int ivarf = CCTK_FirstVarIndex(varname) + CCTK_NumVarsInGroup(varname);

        CCTK_INT bndsize[6], is_ghostbnd[6], is_symbnd[6], is_physbnd[6];
        GetBoundarySizesAndTypes(cctkGH, 6, bndsize, is_ghostbnd, is_symbnd, is_physbnd);
        CCTK_INT ivertex = CCTK_GFINDEX3D(cctkGH,cctk_lsh[0]-bndsize[1]-1,cctk_lsh[1]-bndsize[3]-1,cctk_lsh[2]-bndsize[5]-1);

        for (int ivar=ivari; ivar<ivarf; ivar++)
        {
          CCTK_REAL *pvar = (CCTK_REAL *) CCTK_VarDataPtrI(cctkGH, 0, ivar);
          CCTK_REAL coeff = (pvar[ivertex] - 1) * r[ivertex];

          if (CCTK_Equals(boundary_conditions, "Robin"))
          {
            CCTK_LOOP3_BND(robin, cctkGH,
                           i, j, k, ni, nj, nk)
            {
              int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
 
              pvar[index] = 1 + coeff / r[index];
 
            } CCTK_ENDLOOP3_BND (robin);
          } else if (CCTK_Equals(boundary_conditions, "TwoPunctures")) {
            //Insert a check somewhere that TwoPunctures is active and puncture_u has storage
            CCTK_REAL *tpu = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "TwoPunctures::puncture_u");

            CCTK_LOOP3_BND(tp, cctkGH,
                           i, j, k, ni, nj, nk)
            {
              int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

              pvar[index] = 1.0 + tpu[index];

            } CCTK_ENDLOOP3_BND (tp);
          }
        }
      }

    } END_COMPONENT_LOOP;
  } END_MAP_LOOP;

  return;
}
