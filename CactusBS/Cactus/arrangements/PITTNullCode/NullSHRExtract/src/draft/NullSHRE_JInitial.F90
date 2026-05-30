! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#define full_initial_data .false.

  ! set J initially on the entire characteristic grid

 subroutine NullSHRE_JInitial(CCTK_ARGUMENTS)

    use cctk
    use NullGrid_Vars
    implicit none

    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_FUNCTIONS

    CCTK_INT :: k 
    

    do k = 1, N_radial_pts

       jcn(:,:,k) = j_wt(:,:,1)*(xb(k) - 1.d0)/(x_wt(:,:,1) - 1.d0)
       jcs(:,:,k) = j_wt(:,:,2)*(xb(k) - 1.d0)/(x_wt(:,:,2) - 1.d0)
       
       jcn_p(:,:,k) = 0.d0
       jcs_p(:,:,k) = 0.d0
       

       if (full_initial_data) then

           ucn(:,:,k) = 0.d0
           ucs(:,:,k) = 0.d0
           bcn(:,:,k) = 0.d0
           bcs(:,:,k) = 0.d0
           wcn(:,:,k) = 0.d0
           wcs(:,:,k) = 0.d0

           ucn_p(:,:,k) = 0.d0
           bcn_p(:,:,k) = 0.d0
           wcn_p(:,:,k) = 0.d0
           ucs_p(:,:,k) = 0.d0
           bcs_p(:,:,k) = 0.d0
           wcs_p(:,:,k) = 0.d0
   
           if (first_order_scheme.ne.0) then
  
              nucn(:,:,k) = 0.d0
              ckcn(:,:,k) = 0.d0
              cbcn(:,:,k) = 0.d0 
              nucs(:,:,k) = 0.d0
              ckcs(:,:,k) = 0.d0
              cbcs(:,:,k) = 0.d0

              ucn_p(:,:,k) = 0.d0  
              ckcn_p(:,:,k) = 0.d0
              cbcn_p(:,:,k) = 0.d0
              nucs_p(:,:,k) = 0.d0
              ckcs_p(:,:,k) = 0.d0
              cbcs_p(:,:,k) = 0.d0
   
           end if
      
       end if

    end do

 end subroutine NullSHRE_JInitial
