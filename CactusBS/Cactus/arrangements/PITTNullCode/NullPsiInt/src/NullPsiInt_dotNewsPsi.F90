! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine NullPsiInt_dotNewsPsi(CCTK_ARGUMENTS)
 use cctk
 !use NullInterp
 use NullGrid_Vars
 implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: patch_ID
  CCTK_COMPLEX :: ii
  ii = dcmplx(0.,1.)

!INTEGRATION Int_Psi4 = News
  patch_ID = ip_n-1

  do patch_ID = 1, 2

!Apply mask to the complex Characteristic Fields
    if (CCTK_EQUALS(mask_type, "none")) then

       NewsB_mask(:,:,patch_ID) = NewsB(:,:,patch_ID) 
       Psi4_mask(:,:,patch_ID) = Psi4(:,:,patch_ID) 

    else if (CCTK_EQUALS(mask_type, "EG_mask")) then

       NewsB_mask(:,:,patch_ID) = EG_mask*NewsB(:,:,patch_ID)
       Psi4_mask(:,:,patch_ID) = EG_mask*Psi4(:,:,patch_ID) 

    else if (CCTK_EQUALS(mask_type, "EQ_mask")) then

       NewsB_mask(:,:,patch_ID) = EQ_mask*NewsB(:,:,patch_ID) 
       Psi4_mask(:,:,patch_ID) = EQ_mask*Psi4(:,:,patch_ID)

    else if (CCTK_EQUALS(mask_type, "EV_mask")) then

       NewsB_mask(:,:,patch_ID) = EV_mask*NewsB(:,:,patch_ID) 
       Psi4_mask(:,:,patch_ID) = EV_mask*Psi4(:,:,patch_ID)

    else
      call CCTK_WARN(0, "unsupported circular mask type")
    end if


! the evolution equation

     re_dotNewsB(:,:,patch_ID) = dble(Psi4_mask(:,:,patch_ID))
     im_dotNewsB(:,:,patch_ID) = dimag(Psi4_mask(:,:,patch_ID))

     NewsB_Psi(:,:,patch_ID) = re_PsiInt(:,:,patch_ID) + ii*im_PsiInt(:,:,patch_ID)

  end do

end subroutine NullPsiInt_dotNewsPsi 
