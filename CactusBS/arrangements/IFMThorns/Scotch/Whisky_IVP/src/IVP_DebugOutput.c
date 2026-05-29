 /*@@
  *    @file      BAM_Calls.c
  *    @date      Sun Feb 04 11:55:46 2003
  *    @author    Frank Loeffler
  *    @desc
  *               Debugoutput
  *    @enddesc 
  *    @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void Call_IVP_Debug_Output(cGH *cctkGH, CCTK_INT kind);
void CCTK_FCALL CCTK_FNAME(Call_IVP_Debug_Output)
                          (cGH **cctkGH, CCTK_INT *kind);

void Call_IVP_Debug_Output(cGH *cctkGH, CCTK_INT kind)
{
  if (kind==0)
  {
    CCTK_OutputVarAs(cctkGH, "Whisky_IVP::IVP_Psi", "0IVP_Psi");
    CCTK_OutputVarAs(cctkGH, "Whisky_IVP::IVP_Mlinear", "0IVP_Mlinear");
    CCTK_OutputVarAs(cctkGH, "Whisky_IVP::IVP_Nsource", "0IVP_NSource");
    CCTK_OutputVarAs(cctkGH, "Whisky_IVP::IVP_Psi5Source", "0IVP_Psi5Source");
    CCTK_OutputVarAs(cctkGH, "Whisky_IVP::IVP_Psim7Source", "0IVP_Psim7Source");
    CCTK_OutputVarAs(cctkGH, "Whisky_IVP::IVP_Psim3Source", "0IVP_Psim3Source");

    CCTK_OutputVarAs(cctkGH, "Whisky_IVP::IVP_uxx", "0IVP_uxx");
    CCTK_OutputVarAs(cctkGH, "Whisky_IVP::IVP_uyy", "0IVP_uyy");
    CCTK_OutputVarAs(cctkGH, "Whisky_IVP::IVP_uzz", "0IVP_uzz");
    CCTK_OutputVarAs(cctkGH, "Whisky_IVP::IVP_uxy", "0IVP_uxy");
    CCTK_OutputVarAs(cctkGH, "Whisky_IVP::IVP_uyz", "0IVP_uyz");
    CCTK_OutputVarAs(cctkGH, "Whisky_IVP::IVP_uxz", "0IVP_uxz");

    CCTK_OutputVarAs(cctkGH, "ADMBase::gxx", "0Base_gxx");
    CCTK_OutputVarAs(cctkGH, "ADMBase::gyy", "0Base_gyy");
    CCTK_OutputVarAs(cctkGH, "ADMBase::gzz", "0Base_gzz");
    CCTK_OutputVarAs(cctkGH, "ADMBase::gxy", "0Base_gxy");
    CCTK_OutputVarAs(cctkGH, "ADMBase::gyz", "0Base_gyz");
    CCTK_OutputVarAs(cctkGH, "ADMBase::gxz", "0Base_gxz");
/*
    CCTK_OutputVarAs(cctkGH, "StaticConformal::Psi", "0SCPsi");
*/
  }
  else
  {
    CCTK_OutputVarAs(cctkGH, "Whisky_IVP::IVP_Psi", "1IVP_Psi");
    CCTK_OutputVarAs(cctkGH, "Whisky_IVP::IVP_Mlinear", "1IVP_Mlinear");
    CCTK_OutputVarAs(cctkGH, "Whisky_IVP::IVP_Nsource", "1IVP_NSource");
    CCTK_OutputVarAs(cctkGH, "Whisky_IVP::IVP_Psi5Source", "1IVP_Psi5Source");
    CCTK_OutputVarAs(cctkGH, "Whisky_IVP::IVP_Psim7Source", "1IVP_Psim7Source");
    CCTK_OutputVarAs(cctkGH, "Whisky_IVP::IVP_Psim3Source", "1IVP_Psim3Source");

    CCTK_OutputVarAs(cctkGH, "Whisky_IVP::IVP_uxx", "1IVP_uxx");
    CCTK_OutputVarAs(cctkGH, "Whisky_IVP::IVP_uyy", "1IVP_uyy");
    CCTK_OutputVarAs(cctkGH, "Whisky_IVP::IVP_uzz", "1IVP_uzz");
    CCTK_OutputVarAs(cctkGH, "Whisky_IVP::IVP_uxy", "1IVP_uxy");
    CCTK_OutputVarAs(cctkGH, "Whisky_IVP::IVP_uyz", "1IVP_uyz");
    CCTK_OutputVarAs(cctkGH, "Whisky_IVP::IVP_uxz", "1IVP_uxz");

    CCTK_OutputVarAs(cctkGH, "ADMBase::gxx", "1Base_gxx");
    CCTK_OutputVarAs(cctkGH, "ADMBase::gyy", "1Base_gyy");
    CCTK_OutputVarAs(cctkGH, "ADMBase::gzz", "1Base_gzz");
    CCTK_OutputVarAs(cctkGH, "ADMBase::gxy", "1Base_gxy");
    CCTK_OutputVarAs(cctkGH, "ADMBase::gyz", "1Base_gyz");
    CCTK_OutputVarAs(cctkGH, "ADMBase::gxz", "1Base_gxz");
/*
    CCTK_OutputVarAs(cctkGH, "StaticConformal::Psi", "1SCPsi");
*/
  }
} 

void CCTK_FCALL CCTK_FNAME(Call_IVP_Debug_Output)
                            (cGH **cctkGH, CCTK_INT *kind)
{
    Call_IVP_Debug_Output(*cctkGH, *kind);
    return;
}
